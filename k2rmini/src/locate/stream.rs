//! One query pass per genome: SIMD lookup and exact alignment evidence are
//! computed together. Reduction never reads genome bases or compressed caches.
use super::*;
use std::{
    collections::BTreeMap,
    fs::File,
    io::{BufWriter, Write},
    ops::Range,
    sync::atomic::{AtomicUsize, Ordering},
    thread,
    time::Instant,
};

struct ContigMeta {
    name: String,
    start: u64,
    len: usize,
}

struct GenomeMeta {
    input: PathBuf,
    contigs: Vec<ContigMeta>,
}

struct PackedContig {
    forward: PackedSeqVec,
    reverse: PackedSeqVec,
    runs: Vec<Range<usize>>,
}

impl PackedContig {
    fn new(sequence: &[u8], k: usize) -> Self {
        let forward = PackedSeqVec::from_ascii(sequence);
        let reverse = forward.as_slice().to_revcomp();
        let mut runs = Vec::new();
        let mut pos = 0;
        while pos < sequence.len() {
            while pos < sequence.len() && !is_dna(sequence[pos]) {
                pos += 1;
            }
            let start = pos;
            while pos < sequence.len() && is_dna(sequence[pos]) {
                pos += 1;
            }
            if pos - start >= k {
                runs.push(start..pos);
            }
        }
        Self {
            forward,
            reverse,
            runs,
        }
    }

    fn len(&self) -> usize {
        self.forward.len()
    }

    fn bounds(&self, start: usize, len: usize) -> io::Result<Range<usize>> {
        let next = self.runs.partition_point(|run| run.start <= start);
        next.checked_sub(1)
            .and_then(|i| self.runs.get(i))
            .filter(|r| start + len <= r.end)
            .cloned()
            .ok_or_else(|| invalid("hit outside an unambiguous DNA run"))
    }
}

// A proto-MUM has many independently usable positions. Each observed exact
// range retains its proto coordinates and its physical genome occurrence.
#[derive(Clone, Copy, Debug)]
struct Fragment {
    begin: usize,
    end: usize,
    hit: Hit,
}

struct Reference {
    contigs: Vec<PackedContig>,
    meta: GenomeMeta,
    fragments: Vec<Vec<Fragment>>,
}

struct Anchor {
    hit: Hit,
    len: usize,
}

#[derive(Default, Clone, Copy)]
struct ScanStats {
    genomes: usize,
    contexts: usize,
    alignments: usize,
    cache_hits: usize,
    extended_bases: u64,
}

impl ScanStats {
    fn add(&mut self, other: Self) {
        self.genomes += other.genomes;
        self.contexts += other.contexts;
        self.alignments += other.alignments;
        self.cache_hits += other.cache_hits;
        self.extended_bases += other.extended_bases;
    }
}

struct SeedPosting {
    pattern: u32,
    offset: u32,
    next: u32,
    reverse: bool,
}

struct FragmentIndex {
    options: LocateOptions,
    patterns: Vec<PackedContig>,
    heads: FxHashMap<u32, u32>,
    postings: Vec<SeedPosting>,
    probe: Vec<u64>,
    minis: FxHashSet<u64>,
}

impl FragmentIndex {
    fn new(source: ProtoIndex) -> io::Result<Self> {
        let mut this = Self {
            options: source.options,
            patterns: Vec::new(),
            heads: FxHashMap::default(),
            postings: Vec::new(),
            probe: Vec::new(),
            minis: FxHashSet::default(),
        };
        let k = this.options.k;
        let mut hashes = Vec::new();
        let mut positions = Vec::new();
        for (id, pattern) in source.patterns.iter().enumerate() {
            let dna: Vec<_> = (pattern.start..pattern.start + pattern.len)
                .map(|i| source.sequences.as_slice().get_ascii(i))
                .collect();
            let packed = PackedContig::new(&dna, k);
            for reverse in [false, true] {
                if reverse && this.options.forward_only {
                    continue;
                }
                let seq = if reverse {
                    &packed.reverse
                } else {
                    &packed.forward
                };
                hashes.clear();
                NtHasher::<false>::new(k)
                    .hash_kmers_simd(seq.as_slice(), 1)
                    .collect_into(&mut hashes);
                positions.clear();
                minimizers(this.options.m, k - this.options.m + 1)
                    .run(seq.as_slice(), &mut positions);
                for &pos in &positions {
                    this.minis.insert(
                        seq.slice(pos as usize..pos as usize + this.options.m)
                            .as_u64(),
                    );
                }
                for (pos, &hash) in hashes.iter().enumerate() {
                    let posting = u32::try_from(this.postings.len())
                        .map_err(|_| invalid("too many proto positions"))?;
                    let offset = if reverse { pattern.len - pos - k } else { pos };
                    let offset = u32::try_from(offset)
                        .map_err(|_| invalid("proto-MUM exceeds position range"))?;
                    let next = this.heads.insert(hash, posting).unwrap_or(NO_POSTING);
                    this.postings.push(SeedPosting {
                        pattern: id as u32,
                        offset,
                        next,
                        reverse,
                    });
                }
            }
            this.patterns.push(packed);
        }
        let bits = (this.heads.len().max(64) * 16).next_power_of_two();
        this.probe = vec![0; bits / 64];
        for &hash in this.heads.keys() {
            let bit = hash as usize & (bits - 1);
            this.probe[bit / 64] |= 1 << (bit % 64);
        }
        Ok(this)
    }

    // SIMD hashes visit every seed position. A cached exact range consumes all
    // later seed hits on that alignment without extending or emitting it again.
    // For a fixed diagonal, proto positions are monotone during this contig
    // scan, so only its last range needs to be retained.
    fn scan(
        &self,
        query: &PackedContig,
        contig: u32,
        origin: u64,
        mut accept: impl FnMut(Fragment) -> io::Result<()>,
    ) -> io::Result<()> {
        let k = self.options.k;
        let mut hashes = Vec::new();
        let mut mini_pos = Vec::new();
        let mut covered: FxHashMap<(u32, bool, i64), (usize, usize)> = FxHashMap::default();
        for run in &query.runs {
            for start in (run.start..run.end).step_by(self.options.chunk_bases) {
                let core_end = start.saturating_add(self.options.chunk_bases).min(run.end);
                let end = core_end.saturating_add(k - 1).min(run.end);
                if end - start < k {
                    continue;
                }
                let seq = query.forward.slice(start..end);
                mini_pos.clear();
                minimizers(self.options.m, k - self.options.m + 1).run(seq, &mut mini_pos);
                if !mini_pos.iter().any(|&p| {
                    self.minis
                        .contains(&seq.slice(p as usize..p as usize + self.options.m).as_u64())
                }) {
                    continue;
                }
                hashes.clear();
                NtHasher::<false>::new(k)
                    .hash_kmers_simd(seq, 1)
                    .collect_into(&mut hashes);
                for (local, &hash) in hashes.iter().take(core_end - start).enumerate() {
                    let bit = hash as usize & (self.probe.len() * 64 - 1);
                    if self.probe[bit / 64] & (1 << (bit % 64)) == 0 {
                        continue;
                    }
                    let Some(&head) = self.heads.get(&hash) else {
                        continue;
                    };
                    let physical_q = start + local;
                    let mut next = head;
                    while next != NO_POSTING {
                        let posting = &self.postings[next as usize];
                        next = posting.next;
                        let r = posting.offset as usize;
                        let q = if posting.reverse {
                            query.len() - physical_q - k
                        } else {
                            physical_q
                        };
                        let key = (posting.pattern, posting.reverse, q as i64 - r as i64);
                        if covered
                            .get(&key)
                            .is_some_and(|&(a, b)| a <= r && r + k <= b)
                        {
                            continue;
                        }
                        let proto = &self.patterns[posting.pattern as usize];
                        let (qforward, qbackward, qb) = if posting.reverse {
                            (
                                &query.reverse,
                                &query.forward,
                                query.len() - run.end..query.len() - run.start,
                            )
                        } else {
                            (&query.forward, &query.reverse, run.clone())
                        };
                        // Hashes only select candidates; equality and all range
                        // boundaries are checked against packed DNA exactly.
                        if proto
                            .forward
                            .slice(r..r + k)
                            .cmp_lcp(&qforward.slice(q..q + k))
                            .1
                            != k
                        {
                            continue;
                        }
                        let right = proto
                            .forward
                            .slice(r + k..proto.len())
                            .cmp_lcp(&qforward.slice(q + k..qb.end))
                            .1;
                        let left = proto
                            .reverse
                            .slice(proto.len() - r..proto.len())
                            .cmp_lcp(&qbackward.slice(query.len() - q..query.len() - qb.start))
                            .1;
                        let (begin, end) = (r - left, r + k + right);
                        covered.insert(key, (begin, end));
                        let physical_start = if posting.reverse {
                            physical_q - right
                        } else {
                            physical_q - left
                        };
                        accept(Fragment {
                            begin,
                            end,
                            hit: Hit {
                                pattern: posting.pattern,
                                contig,
                                start: origin + physical_start as u64,
                                reverse: posting.reverse,
                            },
                        })?;
                    }
                }
            }
        }
        Ok(())
    }
}

fn header_name(parser: &FastxParser<CONFIG>) -> String {
    String::from_utf8_lossy(parser.get_header())
        .split_whitespace()
        .next()
        .unwrap_or("")
        .to_string()
}

fn read_reference(index: &FragmentIndex, input: &Path) -> io::Result<Reference> {
    let mut reference = Reference {
        contigs: Vec::new(),
        meta: GenomeMeta {
            input: input.into(),
            contigs: Vec::new(),
        },
        fragments: (0..index.patterns.len()).map(|_| Vec::new()).collect(),
    };
    let data = FileInput::new(input)?;
    if data.first_byte() != 0 {
        let mut parser = FastxParser::<CONFIG>::from_input(data)?;
        let mut origin = 0;
        while parser.next().is_some() {
            let id =
                u32::try_from(reference.contigs.len()).map_err(|_| invalid("too many contigs"))?;
            let sequence = parser.get_dna_string();
            let packed = PackedContig::new(sequence, index.options.k);
            reference.meta.contigs.push(ContigMeta {
                name: header_name(&parser),
                start: origin,
                len: sequence.len(),
            });
            index.scan(&packed, id, origin, |fragment| {
                reference.fragments[fragment.hit.pattern as usize].push(fragment);
                Ok(())
            })?;
            origin += packed.len() as u64;
            reference.contigs.push(packed);
        }
    }
    for fragments in &mut reference.fragments {
        fragments
            .sort_unstable_by_key(|f| (f.begin, f.end, f.hit.contig, f.hit.start, f.hit.reverse));
    }
    Ok(reference)
}

// Maximal exact intervals in reference coordinates. These are position validity
// ranges, not a single accept/reject bit for a complete proto-MUM.
#[derive(Clone, Copy, Debug)]
struct Record {
    begin: u64,
    end: u64,
    hit: Hit,
}
const RECORD_BYTES: usize = 32;

struct ContextSink {
    shards: Vec<Vec<Record>>,
    reference_bases: u64,
    k: u64,
}

fn context_path(scratch: &Path, genome: usize, shard: usize) -> PathBuf {
    scratch.join(format!("contexts-{genome}-{shard}.zst"))
}

impl ContextSink {
    fn new(shards: usize, reference: &Reference, k: usize) -> Self {
        Self {
            shards: (0..shards).map(|_| Vec::new()).collect(),
            reference_bases: reference
                .meta
                .contigs
                .iter()
                .map(|c| c.len as u64)
                .sum::<u64>()
                .max(1),
            k: k as u64,
        }
    }
    fn shard(&self, start: u64) -> usize {
        ((start as u128 * self.shards.len() as u128) / self.reference_bases as u128) as usize
    }
    fn push(&mut self, record: Record) {
        let first = self.shard(record.begin);
        let last = self.shard(record.end - self.k);
        for shard in first..=last {
            self.shards[shard].push(record);
        }
    }
    fn finish(mut self, scratch: &Path, genome: usize) -> io::Result<()> {
        let mut compressor = zstd::bulk::Compressor::new(1)?;
        let mut bytes = Vec::new();
        for (shard, records) in self.shards.iter_mut().enumerate() {
            records.sort_unstable_by_key(|r| {
                (r.begin, r.end, r.hit.contig, r.hit.start, r.hit.reverse)
            });
            records.dedup_by_key(|r| (r.begin, r.end, r.hit.contig, r.hit.start, r.hit.reverse));
            bytes.clear();
            bytes.reserve(records.len() * RECORD_BYTES);
            for record in records.iter() {
                bytes.extend_from_slice(&record.begin.to_le_bytes());
                bytes.extend_from_slice(&record.end.to_le_bytes());
                bytes.extend_from_slice(&record.hit.start.to_le_bytes());
                bytes.extend_from_slice(&record.hit.contig.to_le_bytes());
                bytes.extend_from_slice(&u32::from(record.hit.reverse).to_le_bytes());
            }
            std::fs::write(
                context_path(scratch, genome, shard),
                compressor.compress(&bytes)?,
            )?;
        }
        Ok(())
    }
}

fn seed_at(fragment: Fragment, pos: usize, k: usize) -> Hit {
    Hit {
        start: fragment.hit.start
            + if fragment.hit.reverse {
                (fragment.end - pos - k) as u64
            } else {
                (pos - fragment.begin) as u64
            },
        ..fragment.hit
    }
}

fn consume_fragment(
    reference: &Reference,
    fragment: Fragment,
    query: &PackedContig,
    origin: u64,
    k: usize,
    cache: &mut AlignmentCache,
    stats: &mut ScanStats,
    sink: &mut ContextSink,
) -> io::Result<()> {
    for &source in &reference.fragments[fragment.hit.pattern as usize] {
        let pos = source.begin.max(fragment.begin);
        if pos + k > source.end.min(fragment.end) {
            continue;
        }
        let anchor = Anchor {
            hit: seed_at(source, pos, k),
            len: k,
        };
        let hit = seed_at(fragment, pos, k);
        let before = stats.alignments;
        let (left, right) =
            alignment_context(reference, &anchor, query, hit, origin, cache, stats)?;
        if before == stats.alignments {
            continue;
        }
        let begin = anchor.hit.start - if anchor.hit.reverse { right } else { left };
        let end = anchor.hit.start + k as u64 + if anchor.hit.reverse { left } else { right };
        sink.push(Record {
            begin,
            end,
            hit: Hit {
                start: hit.start - if hit.reverse { right } else { left },
                reverse: hit.reverse != anchor.hit.reverse,
                ..hit
            },
        });
        stats.contexts += 1;
    }
    Ok(())
}

fn reference_contexts(
    reference: &Reference,
    shards: usize,
    k: usize,
    scratch: &Path,
) -> io::Result<ScanStats> {
    let mut sink = ContextSink::new(shards, reference, k);
    let mut stats = ScanStats {
        genomes: 1,
        ..Default::default()
    };
    // Group fragment occurrences by contig to reuse the same exact alignment
    // cache across all proto-MUMs, including repeat and reverse occurrences.
    for (id, query) in reference.contigs.iter().enumerate() {
        let mut cache = AlignmentCache::default();
        for fragments in &reference.fragments {
            for &fragment in fragments {
                if fragment.hit.contig as usize == id {
                    consume_fragment(
                        reference,
                        fragment,
                        query,
                        reference.meta.contigs[id].start,
                        k,
                        &mut cache,
                        &mut stats,
                        &mut sink,
                    )?;
                }
            }
        }
    }
    sink.finish(scratch, 0)?;
    Ok(stats)
}

fn scan_query(
    index: &FragmentIndex,
    reference: &Reference,
    input: &Path,
    id: usize,
    shards: usize,
    scratch: &Path,
) -> io::Result<(GenomeMeta, ScanStats)> {
    let mut sink = ContextSink::new(shards, reference, index.options.k);
    let mut stats = ScanStats {
        genomes: 1,
        ..Default::default()
    };
    let mut meta = GenomeMeta {
        input: input.into(),
        contigs: Vec::new(),
    };
    let data = FileInput::new(input)?;
    if data.first_byte() != 0 {
        let mut parser = FastxParser::<CONFIG>::from_input(data)?;
        let mut origin = 0;
        while parser.next().is_some() {
            let contig =
                u32::try_from(meta.contigs.len()).map_err(|_| invalid("too many contigs"))?;
            let sequence = parser.get_dna_string();
            let packed = PackedContig::new(sequence, index.options.k);
            meta.contigs.push(ContigMeta {
                name: header_name(&parser),
                start: origin,
                len: sequence.len(),
            });
            let mut cache = AlignmentCache::default();
            index.scan(&packed, contig, origin, |fragment| {
                consume_fragment(
                    reference,
                    fragment,
                    &packed,
                    origin,
                    index.options.k,
                    &mut cache,
                    &mut stats,
                    &mut sink,
                )
            })?;
            origin += packed.len() as u64;
        }
    }
    sink.finish(scratch, id)?;
    eprintln!(
        "[genome_scan]\t{}\t{} validity ranges\t{} alignments\t{} cached",
        id + 1,
        stats.contexts,
        stats.alignments,
        stats.cache_hits
    );
    Ok((meta, stats))
}

#[derive(Clone, Copy)]
struct Context {
    genome: usize,
    record: Record,
}
struct VerifiedMatch {
    len: u64,
    hits: Vec<Hit>,
}

fn interval_hit(record: Record, begin: u64, end: u64) -> Hit {
    Hit {
        start: record.hit.start
            + if record.hit.reverse {
                record.end - end
            } else {
                begin - record.begin
            },
        ..record.hit
    }
}

// As starts advance, each genome contributes its furthest surviving end.
// Their minimum is the next common range. A second physical occurrence reaching
// that end invalidates uniqueness, while shorter occurrences remain harmless.
fn reduce_ranges(
    contexts: &mut [Context],
    k: u64,
    genomes: usize,
    shard: usize,
    shards: usize,
    bases: u64,
) -> Vec<VerifiedMatch> {
    contexts.sort_unstable_by_key(|c| c.record.begin);
    let mut best = vec![None::<usize>; genomes];
    let mut second = vec![None::<u64>; genomes];
    let mut present = 0;
    let mut previous_end = 0;
    let mut cursor = 0;
    let mut output = Vec::new();
    while cursor < contexts.len() {
        let begin = contexts[cursor].record.begin;
        while cursor < contexts.len() && contexts[cursor].record.begin == begin {
            let c = &contexts[cursor];
            let g = c.genome;
            if let Some(old) = best[g] {
                if c.record.end > contexts[old].record.end {
                    second[g] = Some(second[g].unwrap_or(0).max(contexts[old].record.end));
                    best[g] = Some(cursor);
                } else {
                    second[g] = Some(second[g].unwrap_or(0).max(c.record.end));
                }
            } else {
                best[g] = Some(cursor);
                present += 1;
            }
            cursor += 1;
        }
        if present != genomes {
            continue;
        }
        let end = best
            .iter()
            .map(|i| contexts[i.unwrap()].record.end)
            .min()
            .unwrap();
        if end <= previous_end {
            continue;
        }
        previous_end = end;
        if end < begin + k || ((begin as u128 * shards as u128) / bases as u128) as usize != shard {
            continue;
        }
        let mut hits: Vec<_> = best
            .iter()
            .map(|i| interval_hit(contexts[i.unwrap()].record, begin, end))
            .collect();
        if second.iter().any(|e| e.is_some_and(|e| e >= end)) {
            let mut ambiguous = false;
            for c in &contexts[..cursor] {
                if c.record.end < end {
                    continue;
                }
                let candidate = interval_hit(c.record, begin, end);
                let chosen = &mut hits[c.genome];
                if candidate.contig != chosen.contig || candidate.start != chosen.start {
                    ambiguous = true;
                    break;
                }
                // A reverse-complement palindrome is one physical occurrence.
                if !candidate.reverse {
                    *chosen = candidate;
                }
            }
            if ambiguous {
                continue;
            }
        }
        output.push(VerifiedMatch {
            len: end - begin,
            hits,
        });
    }
    output
}

fn reduce_shard(
    reference: &Reference,
    scratch: &Path,
    shard: usize,
    shards: usize,
    genomes: usize,
    k: usize,
) -> io::Result<Vec<VerifiedMatch>> {
    let mut contexts = Vec::new();
    for genome in 0..genomes {
        let bytes = zstd::stream::decode_all(File::open(context_path(scratch, genome, shard))?)?;
        if bytes.len() % RECORD_BYTES != 0 {
            return Err(invalid("invalid validity range record size"));
        }
        for bytes in bytes.chunks_exact(RECORD_BYTES) {
            contexts.push(Context {
                genome,
                record: Record {
                    begin: u64::from_le_bytes(bytes[..8].try_into().unwrap()),
                    end: u64::from_le_bytes(bytes[8..16].try_into().unwrap()),
                    hit: Hit {
                        pattern: 0,
                        start: u64::from_le_bytes(bytes[16..24].try_into().unwrap()),
                        contig: u32::from_le_bytes(bytes[24..28].try_into().unwrap()),
                        reverse: u32::from_le_bytes(bytes[28..32].try_into().unwrap()) != 0,
                    },
                },
            });
        }
    }
    let bases = reference
        .meta
        .contigs
        .iter()
        .map(|c| c.len as u64)
        .sum::<u64>()
        .max(1);
    let output = reduce_ranges(&mut contexts, k as u64, genomes, shard, shards, bases);
    Ok(output)
}
// For a fixed reference contig, orientation and diagonal, a maximal exact
// alignment interval can serve every seed inside it. No probabilistic LCEs.
type AlignmentCache = FxHashMap<(u32, bool, i64), BTreeMap<usize, usize>>;

fn alignment_context(
    reference: &Reference,
    anchor: &Anchor,
    query: &PackedContig,
    hit: Hit,
    origin: u64,
    cache: &mut AlignmentCache,
    stats: &mut ScanStats,
) -> io::Result<(u64, u64)> {
    let rcontig = anchor.hit.contig as usize;
    let rseq = &reference.contigs[rcontig];
    let r = (anchor.hit.start - reference.meta.contigs[rcontig].start) as usize;
    let physical_q = (hit.start - origin) as usize;
    let reverse = anchor.hit.reverse != hit.reverse;
    let q = if reverse {
        query.len() - physical_q - anchor.len
    } else {
        physical_q
    };
    let key = (anchor.hit.contig, reverse, q as i64 - r as i64);
    let intervals = cache.entry(key).or_default();
    if let Some((&begin, &end)) = intervals.range(..=r).next_back() {
        if end >= r + anchor.len {
            stats.cache_hits += 1;
            let (left, right) = ((r - begin) as u64, (end - r - anchor.len) as u64);
            return Ok(if anchor.hit.reverse {
                (right, left)
            } else {
                (left, right)
            });
        }
    }
    let rb = rseq.bounds(r, anchor.len)?;
    let qb = query.bounds(physical_q, anchor.len)?;
    let qb = if reverse {
        query.len() - qb.end..query.len() - qb.start
    } else {
        qb
    };
    let (qforward, qbackward) = if reverse {
        (&query.reverse, &query.forward)
    } else {
        (&query.forward, &query.reverse)
    };
    // packed-seq compares 29 bases per machine word and stops exactly at the
    // first mismatch. Explicit run bounds prevent extending across ambiguity.
    let right = rseq
        .forward
        .slice(r + anchor.len..rb.end)
        .cmp_lcp(&qforward.slice(q + anchor.len..qb.end))
        .1;
    let left = rseq
        .reverse
        .slice(rseq.len() - r..rseq.len() - rb.start)
        .cmp_lcp(&qbackward.slice(query.len() - q..query.len() - qb.start))
        .1;
    intervals.insert(r - left, r + anchor.len + right);
    stats.alignments += 1;
    stats.extended_bases += (left + right) as u64;
    Ok(if anchor.hit.reverse {
        (right as u64, left as u64)
    } else {
        (left as u64, right as u64)
    })
}

fn write_output_shard(
    mut matches: Vec<VerifiedMatch>,
    genomes: &[GenomeMeta],
    scratch: &Path,
    shard: usize,
) -> io::Result<(usize, usize)> {
    for item in &mut matches {
        if item.hits[0].reverse {
            for hit in &mut item.hits {
                hit.reverse = !hit.reverse;
            }
        }
    }
    matches.sort_unstable_by(|a, b| {
        (a.hits[0].start, a.len)
            .cmp(&(b.hits[0].start, b.len))
            .then_with(|| {
                a.hits
                    .iter()
                    .map(|h| (h.contig, h.start, h.reverse))
                    .cmp(b.hits.iter().map(|h| (h.contig, h.start, h.reverse)))
            })
    });
    let total = matches.len();
    matches.dedup_by_key(|m| (m.hits[0].start, m.len));
    let count = matches.len();
    let file = File::create(scratch.join(format!("output-{shard}.tsv.zst")))?;
    let mut writer = zstd::stream::write::Encoder::new(BufWriter::new(file), 1)?;
    if shard == 0 {
        writeln!(
            writer,
            "# Coordinates: zero-based forward-strand starts within each contig; input genome order"
        )?;
        writeln!(writer, "length\tstarts\tstrands\tcontigs")?;
    }
    for item in matches {
        write!(writer, "{}\t", item.len)?;
        for (i, hit) in item.hits.iter().enumerate() {
            if i != 0 {
                write!(writer, ",")?;
            }
            write!(
                writer,
                "{}",
                hit.start - genomes[i].contigs[hit.contig as usize].start
            )?;
        }
        write!(writer, "\t")?;
        for (i, hit) in item.hits.iter().enumerate() {
            if i != 0 {
                write!(writer, ",")?;
            }
            write!(writer, "{}", if hit.reverse { '-' } else { '+' })?;
        }
        write!(writer, "\t")?;
        for (i, hit) in item.hits.iter().enumerate() {
            if i != 0 {
                write!(writer, ",")?;
            }
            write!(writer, "{}", hit.contig)?;
        }
        writeln!(writer)?;
    }
    writer.finish()?.flush()?;
    Ok((count, total - count))
}

#[derive(Debug, Default)]
pub struct LocateTimings {
    pub index: f64,
    pub reference_scan: f64,
    pub scan_extend: f64,
    pub reduce: f64,
    pub output: f64,
    pub finalize: f64,
}

impl LocateTimings {
    pub fn phases(&self) -> [(&'static str, f64); 6] {
        [
            ("proto_index", self.index),
            ("reference_scan", self.reference_scan),
            ("parallel_genome_scan_extend", self.scan_extend),
            ("parallel_context_reduce", self.reduce),
            ("parallel_dedup_output", self.output),
            ("locator_finalize", self.finalize),
        ]
    }
}

#[derive(Debug, Default)]
pub struct LocateStats {
    pub proto_mums: usize,
    pub verified_mums: usize,
    pub merged_duplicates: usize,
    pub source_genomes_read: usize,
    pub alignment_computations: usize,
    pub alignment_cache_hits: usize,
    pub timings: LocateTimings,
}

/// One source pass per genome after GGCAT seed construction. Genome workers
/// produce exact position validity ranges; independent reducer and output
/// shards use all requested workers without reopening any genome sequence.
pub fn locate_mums(
    patterns: &Path,
    inputs: &[PathBuf],
    output: &Path,
    scratch_parent: Option<&Path>,
    options: LocateOptions,
) -> io::Result<LocateStats> {
    if inputs.is_empty() {
        return Err(invalid("at least one genome is required"));
    }
    let started = Instant::now();
    let index = FragmentIndex::new(ProtoIndex::from_fasta(patterns, options)?)?;
    let mut timings = LocateTimings {
        index: started.elapsed().as_secs_f64(),
        ..Default::default()
    };
    let workers = index.options.threads;
    let shards = workers.min(index.patterns.len().max(1));
    eprintln!(
        "K2Rmini: indexed {} proto-MUMs ({} seed hashes, {} minimizers); {} workers; one query pass per genome",
        index.patterns.len(),
        index.heads.len(),
        index.minis.len(),
        workers
    );
    let scratch = if let Some(parent) = scratch_parent {
        std::fs::create_dir_all(parent)?;
        tempfile::Builder::new()
            .prefix("k2rlocate-")
            .tempdir_in(parent)?
    } else {
        tempfile::Builder::new().prefix("k2rlocate-").tempdir()?
    };
    let started = Instant::now();
    let mut reference = read_reference(&index, &inputs[0])?;
    let mut scan_stats = reference_contexts(&reference, shards, index.options.k, scratch.path())?;
    timings.reference_scan = started.elapsed().as_secs_f64();
    let started = Instant::now();
    let next = AtomicUsize::new(1);
    let parts = thread::scope(|scope| -> io::Result<Vec<(usize, GenomeMeta, ScanStats)>> {
        let mut handles = Vec::new();
        for _ in 0..workers.min(inputs.len().saturating_sub(1)) {
            let next = &next;
            let index = &index;
            let reference = &reference;
            let scratch = scratch.path();
            handles.push(scope.spawn(
                move || -> io::Result<Vec<(usize, GenomeMeta, ScanStats)>> {
                    let mut result = Vec::new();
                    loop {
                        let id = next.fetch_add(1, Ordering::Relaxed);
                        if id >= inputs.len() {
                            break;
                        }
                        let (meta, stats) =
                            scan_query(index, reference, &inputs[id], id, shards, scratch)?;
                        result.push((id, meta, stats));
                    }
                    Ok(result)
                },
            ));
        }
        let mut result = Vec::new();
        for handle in handles {
            result.extend(
                handle
                    .join()
                    .map_err(|_| io::Error::other("genome worker panicked"))??,
            );
        }
        Ok(result)
    })?;
    timings.scan_extend = started.elapsed().as_secs_f64();
    let mut metas: Vec<Option<GenomeMeta>> = (0..inputs.len()).map(|_| None).collect();
    for (id, meta, stats) in parts {
        metas[id] = Some(meta);
        scan_stats.add(stats);
    }
    reference.contigs.clear();
    let started = Instant::now();
    let k = index.options.k;
    let reductions =
        thread::scope(|scope| -> io::Result<Vec<Vec<VerifiedMatch>>> {
            let mut handles = Vec::new();
            for shard in 0..shards {
                let reference = &reference;
                let scratch = scratch.path();
                handles.push(scope.spawn(move || {
                    reduce_shard(reference, scratch, shard, shards, inputs.len(), k)
                }));
            }
            handles
                .into_iter()
                .map(|h| {
                    h.join()
                        .map_err(|_| io::Error::other("context reducer panicked"))?
                })
                .collect()
        })?;
    timings.reduce = started.elapsed().as_secs_f64();
    let reference_bases = reference
        .meta
        .contigs
        .iter()
        .map(|c| c.len as u64)
        .sum::<u64>()
        .max(1);
    metas[0] = Some(reference.meta);
    let genomes: Vec<_> = metas.into_iter().map(Option::unwrap).collect();
    let mut stats = LocateStats {
        proto_mums: index.patterns.len(),
        source_genomes_read: scan_stats.genomes,
        alignment_computations: scan_stats.alignments,
        alignment_cache_hits: scan_stats.cache_hits,
        ..Default::default()
    };
    let started = Instant::now();
    let mut output_shards: Vec<Vec<VerifiedMatch>> = (0..shards).map(|_| Vec::new()).collect();
    for matches in reductions {
        for item in matches {
            let shard =
                ((item.hits[0].start as u128 * shards as u128) / reference_bases as u128) as usize;
            output_shards[shard.min(shards - 1)].push(item);
        }
    }
    let output_stats = thread::scope(|scope| -> io::Result<Vec<(usize, usize)>> {
        let mut handles = Vec::new();
        for (shard, matches) in output_shards.into_iter().enumerate() {
            let genomes = &genomes;
            let scratch = scratch.path();
            handles.push(scope.spawn(move || write_output_shard(matches, genomes, scratch, shard)));
        }
        handles
            .into_iter()
            .map(|h| {
                h.join()
                    .map_err(|_| io::Error::other("output worker panicked"))?
            })
            .collect()
    })?;
    for (count, duplicates) in output_stats {
        stats.verified_mums += count;
        stats.merged_duplicates += duplicates;
    }
    timings.output = started.elapsed().as_secs_f64();
    let started = Instant::now();
    let parent = output
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let mut result = tempfile::NamedTempFile::new_in(parent)?;
    for shard in 0..shards {
        io::copy(
            &mut File::open(scratch.path().join(format!("output-{shard}.tsv.zst")))?,
            result.as_file_mut(),
        )?;
    }
    result.as_file_mut().flush()?;
    let mut mapping_name = output.as_os_str().to_os_string();
    mapping_name.push(".contigs.tsv.zst");
    let mut mapping = tempfile::NamedTempFile::new_in(parent)?;
    {
        let mut writer =
            zstd::stream::write::Encoder::new(BufWriter::new(mapping.as_file_mut()), 1)?;
        writeln!(writer, "genome\tinput\tcontig\tname\tlength")?;
        for (id, genome) in genomes.iter().enumerate() {
            for (contig_id, contig) in genome.contigs.iter().enumerate() {
                writeln!(
                    writer,
                    "{id}\t{}\t{contig_id}\t{}\t{}",
                    genome.input.display(),
                    contig.name,
                    contig.len
                )?;
            }
        }
        writer.finish()?.flush()?;
    }
    mapping.persist(PathBuf::from(mapping_name))?;
    result.persist(output)?;
    timings.finalize = started.elapsed().as_secs_f64();
    stats.timings = timings;
    eprintln!(
        "[scan_summary]\t{} source genome reads\t{} contexts\t{} new exact alignments\t{} alignment cache hits\t{} extended bases",
        scan_stats.genomes,
        scan_stats.contexts,
        scan_stats.alignments,
        scan_stats.cache_hits,
        scan_stats.extended_bases
    );
    Ok(stats)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fragment_offsets_work_across_packed_word_sizes() -> io::Result<()> {
        let mut state = 5613u64;
        let core: Vec<u8> = (0..320)
            .map(|_| {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                b"ACGT"[(state >> 62) as usize]
            })
            .collect();
        let temp = tempfile::tempdir()?;
        let patterns = temp.path().join("proto.fa.zst");
        let reference = temp.path().join("reference.fa.zst");
        let text = format!(">core\n{}\n", String::from_utf8_lossy(&core));
        let compressed = zstd::stream::encode_all(text.as_bytes(), 1)?;
        std::fs::write(&patterns, &compressed)?;
        std::fs::write(&reference, &compressed)?;
        for k in [15, 31, 32, 33, 63, 64, 67] {
            let ranges = [0..120, 121 - k..240, 241 - k..320];
            let mut query_text = String::new();
            for (id, range) in ranges.iter().enumerate() {
                let dna: Vec<_> = core[range.clone()]
                    .iter()
                    .rev()
                    .copied()
                    .map(complement)
                    .collect();
                query_text.push_str(&format!(">c{id}\n{}\n", String::from_utf8_lossy(&dna)));
            }
            let query = temp.path().join("query.fa.zst");
            std::fs::write(&query, zstd::stream::encode_all(query_text.as_bytes(), 1)?)?;
            let output = temp.path().join("mums.tsv.zst");
            let stats = locate_mums(
                &patterns,
                &[reference.clone(), query],
                &output,
                Some(temp.path()),
                LocateOptions {
                    k,
                    m: k.min(21),
                    threads: 4,
                    chunk_bases: 7,
                    forward_only: false,
                },
            )?;
            let text = String::from_utf8(zstd::stream::decode_all(File::open(output)?)?).unwrap();
            assert_eq!(stats.verified_mums, 3, "k={k}: {text}");
            for (id, range) in ranges.iter().enumerate() {
                assert!(
                    text.contains(&format!(
                        "{}\t{},0\t+,-\t0,{id}\n",
                        range.len(),
                        range.start
                    )),
                    "k={k}: {text}"
                );
            }
        }
        Ok(())
    }

    #[test]
    fn colliding_seed_hash_is_checked_against_exact_dna() -> io::Result<()> {
        let mut source = ProtoIndex {
            options: LocateOptions {
                k: 3,
                m: 2,
                threads: 2,
                chunk_bases: 2,
                forward_only: false,
            },
            patterns: Vec::new(),
            sequences: PackedSeqVec::default(),
        };
        source.add_pattern(b"ACGATC")?;
        source.add_pattern(b"TTTAAA")?;
        let mut index = FragmentIndex::new(source)?;
        let sequence = PackedContig::new(b"ACGATC", 3);
        let hash = NtHasher::<false>::new(3).hash_seq(sequence.forward.slice(0..3));
        let next = index.heads[&hash];
        let posting = index.postings.len() as u32;
        index.postings.push(SeedPosting {
            pattern: 1,
            offset: 0,
            next,
            reverse: false,
        });
        index.heads.insert(hash, posting);
        let mut hits = Vec::new();
        index.scan(&sequence, 0, 0, |fragment| {
            hits.push(fragment);
            Ok(())
        })?;
        assert!(
            hits.iter()
                .any(|f| f.hit.pattern == 0 && f.begin == 0 && f.end == 6)
        );
        assert!(hits.iter().all(|f| f.hit.pattern != 1));
        Ok(())
    }

    fn reverse(sequence: &[u8]) -> Vec<u8> {
        sequence.iter().rev().copied().map(complement).collect()
    }

    #[test]
    fn exact_cached_contexts_match_basewise_oracle() -> io::Result<()> {
        let mut state = 91u64;
        let mut sequence: Vec<u8> = (0..1400)
            .map(|_| {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                b"ACGT"[(state >> 62) as usize]
            })
            .collect();
        for pos in [0, 97, 431, 997, 1399] {
            sequence[pos] = b'N';
        }
        for reverse_query in [false, true] {
            let mut query = b"AACCA".to_vec();
            query.extend_from_slice(&sequence);
            query.extend_from_slice(b"TGT");
            for pos in [34 + 5, 451 + 5, 1103 + 5] {
                query[pos] = if query[pos] == b'A' { b'C' } else { b'A' };
            }
            let physical = if reverse_query {
                reverse(&query)
            } else {
                query.clone()
            };
            let packed = PackedContig::new(&physical, 7);
            let reference = Reference {
                contigs: vec![PackedContig::new(&sequence, 7)],
                meta: GenomeMeta {
                    input: "reference".into(),
                    contigs: vec![ContigMeta {
                        name: "chr".into(),
                        start: 37,
                        len: sequence.len(),
                    }],
                },
                fragments: Vec::new(),
            };
            let mut cache = AlignmentCache::default();
            let mut stats = ScanStats::default();
            let mut checked = 0;
            for r in (0..sequence.len() - 7).step_by(3) {
                if !sequence[r..r + 7].iter().copied().all(is_dna)
                    || sequence[r..r + 7] != query[r + 5..r + 12]
                {
                    continue;
                }
                let mut left = 0usize;
                while r > left
                    && is_dna(sequence[r - left - 1])
                    && sequence[r - left - 1] == query[r + 5 - left - 1]
                {
                    left += 1;
                }
                let mut right = 0usize;
                while r + 7 + right < sequence.len()
                    && is_dna(sequence[r + 7 + right])
                    && sequence[r + 7 + right] == query[r + 12 + right]
                {
                    right += 1;
                }
                for reverse_anchor in [false, true] {
                    let anchor = Anchor {
                        hit: Hit {
                            pattern: 0,
                            contig: 0,
                            start: 37 + r as u64,
                            reverse: reverse_anchor,
                        },
                        len: 7,
                    };
                    let q = if reverse_query {
                        physical.len() - (r + 5) - 7
                    } else {
                        r + 5
                    };
                    let hit = Hit {
                        pattern: 0,
                        contig: 0,
                        start: 71 + q as u64,
                        reverse: reverse_query != reverse_anchor,
                    };
                    let expected = if reverse_anchor {
                        (right as u64, left as u64)
                    } else {
                        (left as u64, right as u64)
                    };
                    assert_eq!(
                        alignment_context(
                            &reference, &anchor, &packed, hit, 71, &mut cache, &mut stats
                        )?,
                        expected
                    );
                    checked += 1;
                }
            }
            assert!(checked > 500);
            assert!(
                stats.alignments < 20,
                "{} alignments for {} seeds",
                stats.alignments,
                checked
            );
            assert!(stats.cache_hits > 500);
        }
        Ok(())
    }

    #[test]
    fn long_identical_region_is_extended_once_for_many_seeds() -> io::Result<()> {
        let sequence: Vec<_> = b"ACGATCTG".iter().copied().cycle().take(500_003).collect();
        let query = PackedContig::new(&sequence, 31);
        let reference = Reference {
            contigs: vec![PackedContig::new(&sequence, 31)],
            meta: GenomeMeta {
                input: "reference".into(),
                contigs: vec![ContigMeta {
                    name: "chr".into(),
                    start: 0,
                    len: sequence.len(),
                }],
            },
            fragments: Vec::new(),
        };
        let mut stats = ScanStats::default();
        let mut cache = AlignmentCache::default();
        for pos in (13..sequence.len() - 31).step_by(1237) {
            let hit = Hit {
                pattern: 0,
                contig: 0,
                start: pos as u64,
                reverse: false,
            };
            let anchor = Anchor { hit, len: 31 };
            assert_eq!(
                alignment_context(&reference, &anchor, &query, hit, 0, &mut cache, &mut stats)?,
                (pos as u64, (sequence.len() - pos - 31) as u64)
            );
        }
        assert_eq!(stats.alignments, 1);
        assert!(stats.cache_hits > 400);
        assert_eq!(stats.extended_bases, (sequence.len() - 31) as u64);
        Ok(())
    }
}
