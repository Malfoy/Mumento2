//! Variable-support MUMs and reusable, exact position-validity evidence.
//!
//! A match may be absent from genome zero. The first G-N+1 genomes form a
//! hitting set for every support set of size N; each match is emitted only by
//! its first participating genome. Persisted states retain the union of seeds,
//! including below-threshold regions, so future additions can create matches.
use super::*;
use serde::{Deserialize, Serialize};
use std::{fs, io::Read};

#[cfg(test)]
mod corpus_tests;
#[cfg(test)]
mod tests;

const STATE_VERSION: u32 = 1;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct StateInfo {
    pub version: u32,
    pub k: usize,
    pub m: usize,
    pub forward_only: bool,
    /// None means all genomes, including subsequently added genomes.
    pub minimum_genomes: Option<usize>,
    pub genomes: usize,
    pub patterns: usize,
}

fn generation(state: &Path) -> io::Result<PathBuf> {
    let name = fs::read_to_string(state.join("CURRENT"))?;
    let name = name.trim();
    if !name.starts_with("generation-")
        || Path::new(name).components().count() != 1
        || !matches!(
            Path::new(name).components().next(),
            Some(std::path::Component::Normal(_))
        )
    {
        return Err(invalid("invalid state generation name"));
    }
    Ok(state.join(name))
}

pub fn read_state_info(state: &Path) -> io::Result<StateInfo> {
    read_generation_info(&generation(state)?)
}

fn read_generation_info(generation: &Path) -> io::Result<StateInfo> {
    let info: StateInfo = serde_json::from_reader(File::open(generation.join("state.json"))?)?;
    if info.version != STATE_VERSION || info.genomes == 0 {
        return Err(invalid("unsupported or invalid Mumento2 state"));
    }
    Ok(info)
}

struct StateLock {
    _file: File,
}
impl StateLock {
    fn acquire(state: &Path) -> io::Result<Self> {
        fs::create_dir_all(state)?;
        let path = state.join(".lock");
        let file = fs::OpenOptions::new()
            .write(true)
            .read(true)
            .create(true)
            .truncate(false)
            .open(&path)
            .map_err(|e| {
                io::Error::new(e.kind(), format!("lock state {}: {e}", state.display()))
            })?;
        file.try_lock()
            .map_err(|e| io::Error::other(format!("state {} is locked: {e}", state.display())))?;
        // The OS releases this lock even after an interrupted/crashed update.
        Ok(Self { _file: file })
    }
}

fn write_u64(writer: &mut impl Write, n: usize) -> io::Result<()> {
    writer.write_all(&(n as u64).to_le_bytes())
}
fn read_u64(reader: &mut impl Read) -> io::Result<u64> {
    let mut bytes = [0; 8];
    reader.read_exact(&mut bytes)?;
    Ok(u64::from_le_bytes(bytes))
}
fn read_size(reader: &mut impl Read) -> io::Result<usize> {
    usize::try_from(read_u64(reader)?).map_err(|_| invalid("state length exceeds address space"))
}

// Store only forward packed DNA plus ambiguity-run boundaries. Old source files
// can be moved/deleted after a successful checkpoint; they are never reopened.
fn save_reference(reference: &Reference, path: &Path) -> io::Result<()> {
    let mut out = zstd::stream::write::Encoder::new(BufWriter::new(File::create(path)?), 1)?;
    out.write_all(b"M2REF\0\x01\0")?;
    let meta = serde_json::to_vec(&reference.meta)?;
    write_u64(&mut out, meta.len())?;
    out.write_all(&meta)?;
    for contig in &reference.contigs {
        let mut bytes = vec![0; contig.len().div_ceil(4)];
        let seq = contig.forward.as_slice();
        for i in 0..contig.len() {
            bytes[i / 4] |= seq.get(i) << (2 * (i % 4));
        }
        out.write_all(&bytes)?;
        write_u64(&mut out, contig.runs.len())?;
        for run in &contig.runs {
            write_u64(&mut out, run.start)?;
            write_u64(&mut out, run.end)?;
        }
    }
    write_u64(&mut out, reference.fragments.len())?;
    for fragments in &reference.fragments {
        write_u64(&mut out, fragments.len())?;
        for f in fragments {
            write_u64(&mut out, f.begin)?;
            write_u64(&mut out, f.end)?;
            out.write_all(&f.hit.start.to_le_bytes())?;
            out.write_all(&f.hit.contig.to_le_bytes())?;
            out.write_all(&[u8::from(f.hit.reverse)])?;
        }
    }
    out.finish()?.flush()
}

fn load_reference(path: &Path, patterns: usize) -> io::Result<Reference> {
    let mut input = zstd::stream::read::Decoder::new(File::open(path)?)?;
    let mut magic = [0; 8];
    input.read_exact(&mut magic)?;
    if &magic != b"M2REF\0\x01\0" {
        return Err(invalid("invalid packed genome state"));
    }
    let size = read_size(&mut input)?;
    let mut bytes = Vec::new();
    (&mut input).take(size as u64).read_to_end(&mut bytes)?;
    if bytes.len() != size {
        return Err(invalid("truncated genome metadata"));
    }
    let meta: GenomeMeta = serde_json::from_slice(&bytes)?;
    let mut contigs = Vec::new();
    for c in &meta.contigs {
        let mut packed = vec![0; c.len.div_ceil(4)];
        input.read_exact(&mut packed)?;
        let forward = PackedSeqVec::from_raw_parts(packed, c.len);
        let reverse = forward.as_slice().to_revcomp();
        let count = read_size(&mut input)?;
        let mut runs: Vec<Range<usize>> = Vec::new();
        for _ in 0..count {
            let start = read_size(&mut input)?;
            let end = read_size(&mut input)?;
            if start >= end || end > c.len || runs.last().is_some_and(|r| r.end > start) {
                return Err(invalid("invalid saved ambiguity-run boundaries"));
            }
            runs.push(start..end);
        }
        contigs.push(PackedContig {
            forward,
            reverse,
            runs,
        });
    }
    let count = read_size(&mut input)?;
    if count > patterns {
        return Err(invalid("state proto index is inconsistent"));
    }
    let mut fragments: Vec<Vec<Fragment>> = (0..patterns).map(|_| Vec::new()).collect();
    for (pattern, list) in fragments.iter_mut().enumerate().take(count) {
        let count = read_size(&mut input)?;
        for _ in 0..count {
            let begin = read_size(&mut input)?;
            let end = read_size(&mut input)?;
            let start = read_u64(&mut input)?;
            let mut bytes = [0; 5];
            input.read_exact(&mut bytes)?;
            let contig = u32::from_le_bytes(bytes[..4].try_into().unwrap());
            if begin >= end || contig as usize >= contigs.len() || bytes[4] > 1 {
                return Err(invalid("invalid saved proto validity range"));
            }
            list.push(Fragment {
                begin,
                end,
                hit: Hit {
                    pattern: pattern as u32,
                    contig,
                    start,
                    reverse: bytes[4] != 0,
                },
            });
        }
    }
    let mut extra = [0];
    if input.read(&mut extra)? != 0 {
        return Err(invalid("trailing packed genome state data"));
    }
    Ok(Reference {
        contigs,
        meta,
        fragments,
    })
}

fn link_tree(source: &Path, target: &Path) -> io::Result<()> {
    fs::create_dir_all(target)?;
    for entry in fs::read_dir(source)? {
        let entry = entry?;
        let dest = target.join(entry.file_name());
        if entry.file_type()?.is_dir() {
            link_tree(&entry.path(), &dest)?;
        } else if entry.file_type()?.is_file() {
            if fs::hard_link(entry.path(), &dest).is_err() {
                fs::copy(entry.path(), dest)?;
            }
        } else {
            return Err(invalid("unexpected file type in saved state"));
        }
    }
    Ok(())
}

fn pair_contexts(
    reference: &Reference,
    query: &Reference,
    id: usize,
    scratch: &Path,
    k: usize,
) -> io::Result<ScanStats> {
    let mut sink = ContextSink::new(1, reference, k);
    let mut stats = ScanStats::default();
    let mut grouped: Vec<Vec<Fragment>> = (0..query.contigs.len()).map(|_| Vec::new()).collect();
    for &f in query.fragments.iter().flatten() {
        grouped[f.hit.contig as usize].push(f);
    }
    for (contig, fragments) in grouped.into_iter().enumerate() {
        let mut cache = AlignmentCache::default();
        for fragment in fragments {
            consume_fragment(
                reference,
                fragment,
                &query.contigs[contig],
                query.meta.contigs[contig].start,
                k,
                &mut cache,
                &mut stats,
                &mut sink,
            )?;
        }
    }
    sink.finish(scratch, id)?;
    Ok(stats)
}

struct SharedMatch {
    len: u64,
    hits: Vec<Option<Hit>>,
}

fn read_contexts(scratch: &Path, genomes: usize) -> io::Result<Vec<Context>> {
    let mut contexts = Vec::new();
    for genome in 0..genomes {
        let bytes = zstd::stream::decode_all(File::open(context_path(scratch, genome, 0))?)?;
        if bytes.len() % RECORD_BYTES != 0 {
            return Err(invalid("invalid context record length"));
        }
        for b in bytes.chunks_exact(RECORD_BYTES) {
            contexts.push(Context {
                genome,
                record: Record {
                    begin: u64::from_le_bytes(b[..8].try_into().unwrap()),
                    end: u64::from_le_bytes(b[8..16].try_into().unwrap()),
                    hit: Hit {
                        pattern: 0,
                        start: u64::from_le_bytes(b[16..24].try_into().unwrap()),
                        contig: u32::from_le_bytes(b[24..28].try_into().unwrap()),
                        reverse: u32::from_le_bytes(b[28..32].try_into().unwrap()) != 0,
                    },
                },
            });
        }
    }
    Ok(contexts)
}

// Each distinct end changes the support set. A shorter interval shared by more
// genomes may itself be maximal even when a subset shares a longer interval.
// Uniqueness is checked in EVERY genome where the interval occurs, matching
// strict multi-MUM semantics. Missing genomes do not count as repeats.
fn reduce_shared(
    contexts: &mut [Context],
    k: u64,
    genomes: usize,
    minimum: usize,
    anchor: usize,
) -> Vec<SharedMatch> {
    if minimum == genomes {
        return reduce_ranges(contexts, k, genomes, 0, 1, u64::MAX)
            .into_iter()
            .map(|m| SharedMatch {
                len: m.len,
                hits: m.hits.into_iter().map(Some).collect(),
            })
            .collect();
    }
    contexts.sort_unstable_by_key(|c| c.record.begin);
    let mut best = vec![None::<usize>; genomes];
    let mut active: Vec<BTreeMap<u64, Vec<usize>>> =
        (0..genomes).map(|_| BTreeMap::new()).collect();
    let mut cursor = 0;
    let mut output = Vec::new();
    while cursor < contexts.len() {
        let begin = contexts[cursor].record.begin;
        while cursor < contexts.len() && contexts[cursor].record.begin == begin {
            let c = &contexts[cursor];
            if best[c.genome].is_none_or(|i| contexts[i].record.end < c.record.end) {
                best[c.genome] = Some(cursor);
            }
            active[c.genome]
                .entry(c.record.end)
                .or_default()
                .push(cursor);
            cursor += 1;
        }
        for ranges in &mut active {
            while ranges
                .first_key_value()
                .is_some_and(|(&end, _)| end < begin + k)
            {
                ranges.pop_first();
            }
        }
        let mut ends: Vec<_> = best
            .iter()
            .flatten()
            .map(|&i| contexts[i].record.end)
            .filter(|&end| end >= begin + k)
            .collect();
        if ends.len() < minimum {
            continue;
        }
        ends.sort_unstable_by(|a, b| b.cmp(a));
        // At least N genomes reach this end; all smaller distinct ends also
        // qualify, with larger support sets. Longer ends cannot qualify.
        ends.drain(..minimum - 1);
        ends.dedup();
        for end in ends {
            let mut hits: Vec<_> = best
                .iter()
                .map(|entry| {
                    entry.and_then(|i| {
                        let r = contexts[i].record;
                        (r.end >= end).then(|| interval_hit(r, begin, end))
                    })
                })
                .collect();
            if hits.iter().position(Option::is_some) != Some(anchor) {
                continue;
            }
            let mut leftmost = vec![u64::MAX; genomes];
            let mut ambiguous = false;
            for (g, chosen) in hits.iter_mut().enumerate() {
                let Some(chosen) = chosen else {
                    continue;
                };
                for &i in active[g].range(end..).flat_map(|(_, entries)| entries) {
                    let c = &contexts[i];
                    let candidate = interval_hit(c.record, begin, end);
                    if candidate.contig != chosen.contig || candidate.start != chosen.start {
                        ambiguous = true;
                        break;
                    }
                    leftmost[g] = leftmost[g].min(c.record.begin);
                    // Palindromes count once at a physical position.
                    if !candidate.reverse {
                        *chosen = candidate;
                    }
                }
                if ambiguous {
                    break;
                }
            }
            if ambiguous {
                continue;
            }
            // All occurrences extending left together makes this a suffix of
            // another match with the same support; a lost genome is a boundary.
            if hits
                .iter()
                .enumerate()
                .filter(|(_, h)| h.is_some())
                .all(|(g, _)| leftmost[g] < begin)
            {
                continue;
            }
            if hits[anchor].unwrap().reverse {
                for hit in hits.iter_mut().flatten() {
                    hit.reverse = !hit.reverse;
                }
            }
            output.push(SharedMatch {
                len: end - begin,
                hits,
            });
        }
    }
    output
}

fn write_shared_output(
    mut matches: Vec<SharedMatch>,
    references: &[Reference],
    output: &Path,
) -> io::Result<(usize, usize)> {
    matches.sort_unstable_by(|a, b| {
        a.len.cmp(&b.len).then_with(|| {
            a.hits
                .iter()
                .map(|h| h.map(|h| (h.contig, h.start, h.reverse)))
                .cmp(
                    b.hits
                        .iter()
                        .map(|h| h.map(|h| (h.contig, h.start, h.reverse))),
                )
        })
    });
    let total = matches.len();
    matches.dedup_by(|a, b| a.len == b.len && a.hits == b.hits);
    // Stable coordinate order, with absent initial genomes allowed.
    matches.sort_unstable_by_key(|m| {
        let anchor = m.hits.iter().position(Option::is_some).unwrap();
        (anchor, m.hits[anchor].unwrap().start, m.len)
    });
    let parent = output
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let mut result = tempfile::NamedTempFile::new_in(parent)?;
    {
        let mut out = zstd::stream::write::Encoder::new(BufWriter::new(result.as_file_mut()), 1)?;
        writeln!(
            out,
            "# Coordinates: zero-based forward-strand starts within each contig; input genome order; . = absent"
        )?;
        writeln!(out, "length\tstarts\tstrands\tcontigs")?;
        for m in &matches {
            write!(out, "{}", m.len)?;
            for column in 0..3 {
                write!(out, "\t")?;
                for (g, hit) in m.hits.iter().enumerate() {
                    if g > 0 {
                        write!(out, ",")?;
                    }
                    if let Some(hit) = hit {
                        match column {
                            0 => write!(
                                out,
                                "{}",
                                hit.start - references[g].meta.contigs[hit.contig as usize].start
                            )?,
                            1 => write!(out, "{}", if hit.reverse { '-' } else { '+' })?,
                            _ => write!(out, "{}", hit.contig)?,
                        }
                    } else {
                        write!(out, ".")?;
                    }
                }
            }
            writeln!(out)?;
        }
        out.finish()?.flush()?;
    }
    let mut mapping_name = output.as_os_str().to_os_string();
    mapping_name.push(".contigs.tsv.zst");
    let mut mapping = tempfile::NamedTempFile::new_in(parent)?;
    {
        let mut out = zstd::stream::write::Encoder::new(BufWriter::new(mapping.as_file_mut()), 1)?;
        writeln!(out, "genome\tinput\tcontig\tname\tlength")?;
        for (g, reference) in references.iter().enumerate() {
            for (id, c) in reference.meta.contigs.iter().enumerate() {
                writeln!(
                    out,
                    "{g}\t{}\t{id}\t{}\t{}",
                    reference.meta.input.display(),
                    c.name,
                    c.len
                )?;
            }
        }
        out.finish()?.flush()?;
    }
    mapping.persist(PathBuf::from(mapping_name))?;
    result.persist(output)?;
    Ok((matches.len(), total - matches.len()))
}

/// `patterns` must cover every k-mer retained for discovery. Saving a state
/// requires UNION seeds (presence >= 1), including when N currently means all.
/// Updates build union seeds from NEW documents only and reuse old evidence.
#[allow(clippy::too_many_arguments)]
pub fn locate_shared_mums(
    patterns: &Path,
    inputs: &[PathBuf],
    output: &Path,
    scratch_parent: Option<&Path>,
    options: LocateOptions,
    minimum_genomes: Option<usize>,
    state: Option<&Path>,
    update: Option<&Path>,
) -> io::Result<LocateStats> {
    options.validate()?;
    if inputs.is_empty() {
        return Err(invalid("at least one new genome is required"));
    }
    let destination = state.or(update);
    let _lock = destination.map(StateLock::acquire).transpose()?;
    if update.is_none() && destination.is_some_and(|s| s.join("CURRENT").exists()) {
        return Err(invalid("state already exists; use --update to add genomes"));
    }
    if let (Some(source), Some(dest)) = (update, destination) {
        if source.canonicalize()? != dest.canonicalize()? && dest.join("CURRENT").exists() {
            return Err(invalid("destination state already exists"));
        }
    }
    let previous = update.map(generation).transpose()?;
    let old_info = previous.as_deref().map(read_generation_info).transpose()?;
    if let Some(info) = &old_info {
        if (info.k, info.m, info.forward_only) != (options.k, options.m, options.forward_only) {
            return Err(invalid(
                "updates must retain the state's k, minimizer length and strand mode",
            ));
        }
    }
    let old_count = old_info.as_ref().map_or(0, |s| s.genomes);
    let genomes = old_count + inputs.len();
    let minimum = minimum_genomes.unwrap_or(genomes);
    if minimum == 0 || minimum > genomes {
        return Err(invalid("require 1 <= N <= total genome count"));
    }
    let scratch = if let Some(parent) = destination.or(scratch_parent) {
        fs::create_dir_all(parent)?;
        tempfile::Builder::new()
            .prefix(".mumento2-")
            .tempdir_in(parent)?
    } else {
        tempfile::Builder::new()
            .prefix("mumento2-shared-")
            .tempdir()?
    };
    fs::create_dir_all(scratch.path().join("genomes"))?;
    fs::create_dir_all(scratch.path().join("contexts"))?;
    let started = Instant::now();
    let combined = scratch.path().join("protos.fa.zst");
    let empty_source = || ProtoIndex {
        options: options.clone(),
        patterns: Vec::new(),
        sequences: PackedSeqVec::default(),
    };
    let mut index = FragmentIndex::new(if let Some(previous) = &previous {
        ProtoIndex::from_fasta(&previous.join("protos.fa.zst"), options.clone())?
    } else {
        empty_source()
    })?;
    if old_info
        .as_ref()
        .is_some_and(|s| s.patterns != index.patterns.len())
    {
        return Err(invalid("saved proto index is inconsistent"));
    }
    let mut additions = empty_source();
    {
        let mut out = File::create(&combined)?;
        // Keep the archive in one frame. Short initial frames can look like
        // EOF to parsers which treat a partial first input block as final.
        let mut frame = zstd::stream::write::Encoder::new(&mut out, 1)?;
        if let Some(previous) = &previous {
            let mut old =
                zstd::stream::read::Decoder::new(File::open(previous.join("protos.fa.zst"))?)?;
            io::copy(&mut old, &mut frame)?;
            link_tree(&previous.join("genomes"), &scratch.path().join("genomes"))?;
            link_tree(&previous.join("contexts"), &scratch.path().join("contexts"))?;
        }
        // Preserve old pattern order/IDs and normalize incoming FASTA.
        let input = FileInput::new(patterns)?;
        if input.first_byte() != 0 {
            let mut parser = FastxParser::<CONFIG>::from_input(input)?;
            while parser.next().is_some() {
                let dna = parser.get_dna_string();
                if dna.len() < options.k || !dna.iter().copied().all(is_dna) {
                    return Err(invalid(
                        "proto-MUMs must contain only A/C/G/T and have length >= k",
                    ));
                }
                // Existing seeds already discover any rearrangement made of
                // their k-mers: exact extension is not bounded by proto ends.
                // Append only genuinely new k-mer runs, avoiding duplicate
                // seed postings on every addition of a similar genome.
                let mut covered = vec![false; dna.len() - options.k + 1];
                if previous.is_some() {
                    let packed = PackedContig::new(dna, options.k);
                    index.scan(&packed, 0, 0, |f| {
                        let start = f.hit.start as usize;
                        let end = start + f.end - f.begin - options.k + 1;
                        covered[start..end].fill(true);
                        Ok(())
                    })?;
                }
                let mut pos = 0;
                while pos < covered.len() {
                    if covered[pos] {
                        pos += 1;
                        continue;
                    }
                    let begin = pos;
                    while pos < covered.len() && !covered[pos] {
                        pos += 1;
                    }
                    let novel = &dna[begin..pos + options.k - 1];
                    additions.add_pattern(novel)?;
                    frame.write_all(b">proto\n")?;
                    frame.write_all(novel)?;
                    frame.write_all(b"\n")?;
                }
            }
        }
        frame.finish()?;
    }
    index.extend(additions)?;
    let mut timings = LocateTimings {
        index: started.elapsed().as_secs_f64(),
        ..Default::default()
    };
    eprintln!(
        "K2Rmini: {} union/candidate proto-MUMs; at least {minimum}/{genomes} genomes; {old_count} cached genomes",
        index.patterns.len()
    );
    let started = Instant::now();
    let next = AtomicUsize::new(0);
    let loaded = thread::scope(|scope| -> io::Result<Vec<(usize, Reference)>> {
        let mut handles = Vec::new();
        for _ in 0..options.threads.min(genomes) {
            let next = &next;
            let index = &index;
            let scratch = scratch.path();
            handles.push(scope.spawn(move || -> io::Result<Vec<_>> {
                let mut parts = Vec::new();
                loop {
                    let id = next.fetch_add(1, Ordering::Relaxed);
                    if id >= genomes {
                        break;
                    }
                    let path = scratch.join("genomes").join(format!("{id}.zst"));
                    let reference = if id < old_count {
                        load_reference(&path, index.patterns.len())?
                    } else {
                        let mut r = read_reference(index, &inputs[id - old_count])?;
                        r.meta.input = fs::canonicalize(&r.meta.input)?;
                        if destination.is_some() {
                            save_reference(&r, &path)?;
                        }
                        r
                    };
                    parts.push((id, reference));
                }
                Ok(parts)
            }));
        }
        let mut parts = Vec::new();
        for h in handles {
            parts.extend(
                h.join()
                    .map_err(|_| io::Error::other("genome index worker panicked"))??,
            );
        }
        Ok(parts)
    })?;
    let mut ordered: Vec<Option<Reference>> = (0..genomes).map(|_| None).collect();
    for (id, reference) in loaded {
        ordered[id] = Some(reference);
    }
    let references: Vec<_> = ordered.into_iter().map(Option::unwrap).collect();
    let mut seen = FxHashSet::default();
    for reference in &references {
        if !seen.insert(&reference.meta.input) {
            return Err(invalid("the same genome file was supplied more than once"));
        }
    }
    timings.reference_scan = started.elapsed().as_secs_f64();
    let started = Instant::now();
    let anchors = genomes - minimum + 1;
    let tasks: Vec<_> = (0..anchors)
        .flat_map(|anchor| (0..genomes).map(move |query| (anchor, query)))
        .filter(|&(a, q)| {
            !context_path(&scratch.path().join("contexts").join(a.to_string()), q, 0).exists()
        })
        .collect();
    for a in 0..anchors {
        fs::create_dir_all(scratch.path().join("contexts").join(a.to_string()))?;
    }
    let next = AtomicUsize::new(0);
    let scan_stats = thread::scope(|scope| -> io::Result<ScanStats> {
        let mut handles = Vec::new();
        for _ in 0..options.threads.min(tasks.len()) {
            let next = &next;
            let tasks = &tasks;
            let references = &references;
            let scratch = scratch.path();
            let k = options.k;
            handles.push(scope.spawn(move || -> io::Result<ScanStats> {
                let mut stats = ScanStats::default();
                loop {
                    let id = next.fetch_add(1, Ordering::Relaxed);
                    let Some(&(anchor, query)) = tasks.get(id) else {
                        break;
                    };
                    stats.add(pair_contexts(
                        &references[anchor],
                        &references[query],
                        query,
                        &scratch.join("contexts").join(anchor.to_string()),
                        k,
                    )?);
                }
                Ok(stats)
            }));
        }
        let mut stats = ScanStats::default();
        for h in handles {
            stats.add(
                h.join()
                    .map_err(|_| io::Error::other("validity worker panicked"))??,
            );
        }
        Ok(stats)
    })?;
    timings.scan_extend = started.elapsed().as_secs_f64();
    eprintln!(
        "[update_evidence]\t{} new source reads\t{} new reference/genome pairs\t{} cached pairs reused",
        inputs.len(),
        tasks.len(),
        anchors * genomes - tasks.len()
    );
    let started = Instant::now();
    let next = AtomicUsize::new(0);
    let matches = thread::scope(|scope| -> io::Result<Vec<SharedMatch>> {
        let mut handles = Vec::new();
        for _ in 0..options.threads.min(anchors) {
            let next = &next;
            let scratch = scratch.path();
            let k = options.k as u64;
            handles.push(scope.spawn(move || -> io::Result<Vec<SharedMatch>> {
                let mut result = Vec::new();
                loop {
                    let anchor = next.fetch_add(1, Ordering::Relaxed);
                    if anchor >= anchors {
                        break;
                    }
                    let mut contexts =
                        read_contexts(&scratch.join("contexts").join(anchor.to_string()), genomes)?;
                    result.extend(reduce_shared(&mut contexts, k, genomes, minimum, anchor));
                }
                Ok(result)
            }));
        }
        let mut result = Vec::new();
        for h in handles {
            result.extend(
                h.join()
                    .map_err(|_| io::Error::other("support reducer panicked"))??,
            );
        }
        Ok(result)
    })?;
    timings.reduce = started.elapsed().as_secs_f64();
    let started = Instant::now();
    let (verified_mums, merged_duplicates) = write_shared_output(matches, &references, output)?;
    timings.output = started.elapsed().as_secs_f64();
    let started = Instant::now();
    if let Some(destination) = destination {
        let info = StateInfo {
            version: STATE_VERSION,
            k: options.k,
            m: options.m,
            forward_only: options.forward_only,
            minimum_genomes,
            genomes,
            patterns: index.patterns.len(),
        };
        serde_json::to_writer_pretty(File::create(scratch.path().join("state.json"))?, &info)?;
        let name = format!(
            "generation-{}",
            scratch.path().file_name().unwrap().to_string_lossy()
        );
        let work = scratch.keep();
        fs::rename(work, destination.join(&name))?;
        let mut pointer = tempfile::NamedTempFile::new_in(destination)?;
        writeln!(pointer, "{name}")?;
        pointer.as_file_mut().sync_all()?;
        pointer.persist(destination.join("CURRENT"))?;
        eprintln!(
            "State saved to {} ({genomes} genomes)",
            destination.display()
        );
    }
    timings.finalize = started.elapsed().as_secs_f64();
    Ok(LocateStats {
        proto_mums: index.patterns.len(),
        verified_mums,
        merged_duplicates,
        source_genomes_read: inputs.len(),
        alignment_computations: scan_stats.alignments,
        alignment_cache_hits: scan_stats.cache_hits,
        timings,
    })
}
