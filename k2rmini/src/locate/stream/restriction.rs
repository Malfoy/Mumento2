//! Repeat evidence belongs to individual proto positions. SSHash remains the
//! only k-mer dictionary. Repeated occurrences are partitioned by exact right
//! context until at most two oriented copies remain; no occurrence cross product
//! or suffix array is constructed. Two copies retain physical palindromes.
use super::*;
use sshash_lib::{Kmer, KmerBits};

const NONE: u32 = u32::MAX;
const NO_LOCATION: u64 = u64::MAX;
const LEAF: usize = 1usize << (usize::BITS - 1);

struct ContextGroup {
    depth: usize,
    sample: u64,
    // Internal nodes hold four child IDs. Leaves reuse the first two words
    // for the complete second u64 location; child IDs are never read there.
    edges: [u32; 4],
}
impl ContextGroup {
    fn leaf(depth: usize, sample: u64, second: u64) -> Self {
        Self {
            depth: depth | LEAF,
            sample,
            edges: [second as u32, (second >> 32) as u32, NONE, NONE],
        }
    }
    fn second_location(&self) -> u64 {
        debug_assert!(self.depth & LEAF != 0);
        self.edges[0] as u64 | ((self.edges[1] as u64) << 32)
    }
}
struct RepeatedPosition {
    first: u64,
    second: u64,
    roots: [u32; 2],
    pool: usize,
}
pub(super) struct RestrictionIndex {
    positions: Vec<u64>,
    repeated: Vec<u64>,
    boundaries: Vec<u64>,
    // Direction-specific repeat entry points. A verified alignment can cross
    // proto boundaries through unique seeds and unchanged exact two-copy
    // mapping runs; every new occurrence-set entry remains a stopping point.
    reference_repeats: Vec<u64>,
    reference_repeats_reverse: Vec<u64>,
    groups: Vec<RepeatedPosition>,
    pools: Vec<Vec<ContextGroup>>,
}

// Exact aligned reference intervals for the current query contig. The diagonal
// uses global reference coordinates; interval containment also identifies the
// reference contig and excludes ambiguity boundaries. Collisions only replace
// an entry: the full key and interval are always checked before a hit.
#[derive(Clone, Copy, Default)]
struct RecentAlignment {
    diagonal: i128,
    begin: u64,
    end: u64,
    reverse: bool,
}
struct RecentAlignments {
    entries: [RecentAlignment; 64],
}
impl Default for RecentAlignments {
    fn default() -> Self {
        Self {
            entries: [RecentAlignment::default(); 64],
        }
    }
}
impl RecentAlignments {
    fn key(physical: u64, oriented_query: usize, reverse: bool) -> (usize, i128) {
        let diagonal = oriented_query as i128 - physical as i128;
        let hash = (diagonal as u64 ^ if reverse { 0x517c_c1b7_2722_0a95 } else { 0 })
            .wrapping_mul(0x9e37_79b9_7f4a_7c15);
        ((hash >> 58) as usize, diagonal)
    }
    fn covers(&self, slot: usize, diagonal: i128, reverse: bool, physical: u64, k: usize) -> bool {
        let entry = &self.entries[slot];
        entry.diagonal == diagonal
            && entry.reverse == reverse
            && entry.begin <= physical
            && physical + k as u64 <= entry.end
    }
    fn remember(&mut self, slot: usize, diagonal: i128, reverse: bool, begin: u64, end: u64) {
        self.entries[slot] = RecentAlignment {
            diagonal,
            begin,
            end,
            reverse,
        };
    }
}

struct View<'a> {
    seq: &'a PackedSeqVec,
    start: usize,
    end: usize,
    contig: usize,
    physical: u64,
    reverse: bool,
}
struct SeedView<'a> {
    seq: &'a PackedSeqVec,
    start: usize,
    contig: usize,
    physical: u64,
    reverse: bool,
}
// The caller supplies an indexed, unambiguous seed. Its orientation and
// packed start do not require another search of the reference DNA-run table.
fn seed_view(reference: &Reference, location: u64, direction: bool, k: usize) -> SeedView<'_> {
    let physical = location >> 1;
    let id = if reference.contigs.len() == 1 {
        0
    } else {
        reference
            .meta
            .contigs
            .partition_point(|c| c.start <= physical)
            - 1
    };
    let contig = &reference.contigs[id];
    let local = (physical - reference.meta.contigs[id].start) as usize;
    let reverse = (location & 1 != 0) != direction;
    if reverse {
        SeedView {
            seq: &contig.reverse,
            start: contig.len() - local - k,
            contig: id,
            physical,
            reverse,
        }
    } else {
        SeedView {
            seq: &contig.forward,
            start: local,
            contig: id,
            physical,
            reverse,
        }
    }
}
// Construction and final uniqueness validation require the full DNA-run
// boundary. Lookup paths below already have a proven seed and context depth.
fn view(reference: &Reference, location: u64, direction: bool, k: usize) -> View<'_> {
    let seed = seed_view(reference, location, direction, k);
    let contig = &reference.contigs[seed.contig];
    let local = (seed.physical - reference.meta.contigs[seed.contig].start) as usize;
    let bounds = contig
        .bounds(local, k)
        .expect("indexed seed lies in a DNA run");
    View {
        seq: seed.seq,
        start: seed.start,
        end: if seed.reverse {
            contig.len() - bounds.start
        } else {
            bounds.end
        },
        contig: seed.contig,
        physical: seed.physical,
        reverse: seed.reverse,
    }
}
fn child_base(view: &View<'_>, depth: usize) -> usize {
    view.seq.as_slice().get(view.start + depth) as usize
}

fn refine(
    reference: &Reference,
    locations: Vec<u64>,
    direction: bool,
    known: usize,
    k: usize,
    nodes: &mut Vec<ContextGroup>,
) -> io::Result<u32> {
    let id = u32::try_from(nodes.len()).map_err(|_| invalid("too many proto context groups"))?;
    let sample = locations[0];
    if locations.len() <= 2 {
        nodes.push(ContextGroup::leaf(
            known,
            sample,
            locations.get(1).copied().unwrap_or(NO_LOCATION),
        ));
        return Ok(id);
    }
    // Termination alone cannot resolve a group until only its longest two
    // contexts survive. Ignore shorter ends while finding the first genuine
    // base divergence. This compresses a long tandem run in one group rather
    // than successively removing one overlapping occurrence per level.
    let mut longest = [(0usize, NO_LOCATION); 3];
    for &location in &locations {
        let v = view(reference, location, direction, k);
        let candidate = (v.end - v.start, location);
        if candidate.0 >= longest[2].0 {
            longest[2] = candidate;
            longest.sort_unstable_by(|a, b| b.0.cmp(&a.0));
        }
    }
    let sample = longest[0].1;
    let first = view(reference, sample, direction, k);
    let threshold = longest[2].0 + 1;
    let mut depth = threshold.min(first.end - first.start);
    for &location in &locations {
        if location == sample {
            continue;
        }
        if depth == known {
            break;
        }
        let other = view(reference, location, direction, k);
        let end = (other.end - other.start).min(depth);
        let common = first
            .seq
            .slice(first.start + known..first.start + end)
            .cmp_lcp(&other.seq.slice(other.start + known..other.start + end))
            .1;
        if common < end - known {
            depth = known + common;
        }
    }
    if depth >= threshold {
        nodes.push(ContextGroup::leaf(
            threshold,
            sample,
            if longest[1].0 >= threshold {
                longest[1].1
            } else {
                NO_LOCATION
            },
        ));
        return Ok(id);
    }
    nodes.push(ContextGroup {
        depth,
        sample,
        edges: [NONE; 4],
    });
    let mut parts: [Vec<u64>; 4] = std::array::from_fn(|_| Vec::new());
    for location in locations {
        let v = view(reference, location, direction, k);
        if v.end - v.start > depth {
            parts[child_base(&v, depth)].push(location);
        }
    }
    for (base, part) in parts.into_iter().enumerate() {
        if !part.is_empty() {
            nodes[id as usize].edges[base] =
                refine(reference, part, direction, depth + 1, k, nodes)?;
        }
    }
    Ok(id)
}

impl RestrictionIndex {
    // A normalized oriented edge is used at most once by each copy of
    // either distinct endpoint group. Equal edge/group counts therefore prove
    // an exact bijection of the complete reference occurrence sets. This is
    // full-key equality, not a probabilistic occurrence-set fingerprint.
    fn coalesce_multi_copy_runs(
        copies: &[(usize, u64)],
        positions: &[u64],
        repeated: &[u64],
        groups: usize,
        forward: &mut [u64],
        reverse: &mut [u64],
    ) {
        fn edge(a: usize, b: usize) -> Option<(usize, usize)> {
            (a >> 1 != b >> 1).then(|| (a, b).min((b ^ 1, a ^ 1)))
        }
        let mut counts = vec![0usize; groups];
        let mut blocked = vec![false; groups];
        let mut entries = Vec::new();
        for &(bit, location) in copies {
            if repeated[bit / 64] & (1 << (bit % 64)) == 0 {
                continue;
            }
            let group = positions[bit] as usize;
            counts[group] += 1;
            entries.push((
                (location >> 1) as usize,
                group * 2 + (location & 1) as usize,
            ));
        }
        entries.sort_unstable();
        for adjacent in entries.windows(2) {
            if adjacent[0].0 == adjacent[1].0 {
                blocked[adjacent[0].1 >> 1] = true;
                blocked[adjacent[1].1 >> 1] = true;
            }
        }
        // v21 already handles exact pairs. Omit those and every ambiguous
        // physical entry from the transient multi-copy adjacency graph.
        entries.retain(|&(_, encoded)| counts[encoded >> 1] > 2 && !blocked[encoded >> 1]);
        let mut edges = FxHashMap::<(usize, usize), usize>::default();
        for adjacent in entries.windows(2) {
            let (a, b) = (adjacent[0], adjacent[1]);
            if a.0 + 1 == b.0 {
                if let Some(key) = edge(a.1, b.1) {
                    *edges.entry(key).or_default() += 1;
                }
            }
        }
        let mut coalesced = 0usize;
        for adjacent in entries.windows(2) {
            let (a, b) = (adjacent[0], adjacent[1]);
            if a.0 + 1 != b.0 {
                continue;
            }
            let Some(key) = edge(a.1, b.1) else {
                continue;
            };
            let count = edges[&key];
            if count == counts[a.1 >> 1] && count == counts[b.1 >> 1] {
                let mask = 1 << (b.0 % 64);
                coalesced += usize::from(forward[b.0 / 64] & mask != 0);
                forward[b.0 / 64] &= !mask;
                reverse[a.0 / 64] &= !(1 << (a.0 % 64));
            }
        }
        eprintln!(
            "[multi_copy_reference_runs]\t{} reference copies\t{} oriented edges\t{} unchanged transitions",
            entries.len(),
            edges.len(),
            coalesced
        );
    }

    // Visit the first seed of every exact two-copy mapping run in the scan
    // direction. Once both maximal alignments are recorded there, interior
    // seeds cannot introduce another occurrence. All larger repeat sets and
    // all changes of partner, relative strand or diagonal remain stopping points.
    fn build_reference_repeat_stops(
        reference: &Reference,
        copies: &[(usize, u64)],
        repeated: &[u64],
    ) -> (Vec<u64>, Vec<u64>) {
        struct PairMapping {
            diagonal: i128,
            physical: usize,
            reverse: bool,
        }
        let bases = reference
            .meta
            .contigs
            .last()
            .map_or(0, |c| c.start as usize + c.len);
        let mut forward = vec![0u64; bases.div_ceil(64)];
        let mut blocked = forward.clone();
        let mut pairs = Vec::new();
        let mut start = 0;
        while start < copies.len() {
            let bit = copies[start].0;
            let mut end = start + 1;
            while end < copies.len() && copies[end].0 == bit {
                end += 1;
            }
            if repeated[bit / 64] & (1 << (bit % 64)) != 0 {
                for &(_, location) in &copies[start..end] {
                    let p = (location >> 1) as usize;
                    let mask = 1 << (p % 64);
                    if forward[p / 64] & mask != 0 {
                        blocked[p / 64] |= mask;
                    }
                    forward[p / 64] |= mask;
                }
                if end - start == 2 {
                    let a = copies[start].1;
                    let b = copies[start + 1].1;
                    if a >> 1 != b >> 1 {
                        let reverse = (a ^ b) & 1 != 0;
                        for (own, other) in [(a, b), (b, a)] {
                            let p = (own >> 1) as usize;
                            let q = (other >> 1) as usize;
                            pairs.push(PairMapping {
                                physical: p,
                                reverse,
                                diagonal: if reverse {
                                    p as i128 + q as i128
                                } else {
                                    q as i128 - p as i128
                                },
                            });
                        }
                    }
                }
            }
            start = end;
        }
        let mut reverse = forward.clone();
        pairs.sort_unstable_by_key(|p| p.physical);
        let mut coalesced = 0;
        for adjacent in pairs.windows(2) {
            let (a, b) = (&adjacent[0], &adjacent[1]);
            if a.physical + 1 == b.physical
                && a.diagonal == b.diagonal
                && a.reverse == b.reverse
                && blocked[a.physical / 64] & (1 << (a.physical % 64)) == 0
                && blocked[b.physical / 64] & (1 << (b.physical % 64)) == 0
            {
                forward[b.physical / 64] &= !(1 << (b.physical % 64));
                reverse[a.physical / 64] &= !(1 << (a.physical % 64));
                coalesced += 1;
            }
        }
        eprintln!(
            "[stable_reference_runs]\t{} pair occurrences\t{} unchanged two-copy transitions",
            pairs.len(),
            coalesced
        );
        (forward, reverse)
    }

    fn unique_alignment_span(&self, physical: u64, reverse: bool, extension: u64) -> usize {
        let position = physical as usize;
        let maximum = extension as usize + 1;
        let stops = if reverse {
            &self.reference_repeats_reverse
        } else {
            &self.reference_repeats
        };
        if reverse {
            let lower = position + 1 - maximum;
            let mut upper = position;
            while upper > lower {
                let last = upper - 1;
                let word_index = last / 64;
                let word = stops[word_index] & (u64::MAX >> (63 - last % 64));
                if word != 0 {
                    let boundary = word_index * 64 + 63 - word.leading_zeros() as usize;
                    return if boundary >= lower {
                        position - boundary
                    } else {
                        maximum
                    };
                }
                upper = word_index * 64;
            }
        } else {
            let limit = position + maximum;
            let mut cursor = position + 1;
            while cursor < limit {
                let word_index = cursor / 64;
                let word = stops[word_index] & (u64::MAX << (cursor % 64));
                if word != 0 {
                    let boundary = word_index * 64 + word.trailing_zeros() as usize;
                    return boundary.min(limit) - position;
                }
                cursor = (word_index + 1) * 64;
            }
        }
        maximum
    }

    // A reference fragment maps a proto interval to one physical occurrence.
    // Between fragment boundaries every occurrence continues along the same
    // proto. Once the query matches that proto interval, no interior seed can
    // discover a new left-maximal alignment or a new reference occurrence.
    // This applies to repeats too; every occurrence-set change is visited.
    fn constant_span(&self, bit: usize, backwards: bool, maximum: usize) -> usize {
        if backwards {
            let lower = bit + 2 - maximum;
            let mut upper = bit + 1;
            while upper > lower {
                let last = upper - 1;
                let word_index = last / 64;
                let word = self.boundaries[word_index] & (u64::MAX >> (63 - last % 64));
                if word != 0 {
                    let boundary = word_index * 64 + 63 - word.leading_zeros() as usize;
                    return if boundary >= lower {
                        bit - boundary + 1
                    } else {
                        maximum
                    };
                }
                upper = word_index * 64;
            }
        } else {
            let limit = bit + maximum;
            let mut cursor = bit + 1;
            while cursor < limit {
                let word_index = cursor / 64;
                let word = self.boundaries[word_index] & (u64::MAX << (cursor % 64));
                if word != 0 {
                    let boundary = word_index * 64 + word.trailing_zeros() as usize;
                    return boundary.min(limit) - bit;
                }
                cursor = (word_index + 1) * 64;
            }
        }
        maximum
    }
    pub(super) fn build(index: &FragmentIndex, reference: &Reference) -> io::Result<Self> {
        let started = Instant::now();
        let count = *index.position_offsets.last().unwrap();
        let mut positions = vec![NO_LOCATION; count];
        let mut repeated = vec![0u64; count.div_ceil(64)];
        let mut boundaries = vec![0u64; count.div_ceil(64)];
        let mut copies = Vec::<(usize, u64)>::new();
        for (id, fragments) in reference.fragments.iter().enumerate() {
            for &fragment in fragments {
                for offset in [fragment.begin, fragment.end - index.options.k + 1] {
                    let bit = index.position_offsets[id] + offset;
                    if bit < count {
                        boundaries[bit / 64] |= 1 << (bit % 64);
                    }
                }
                for offset in fragment.begin..=fragment.end - index.options.k {
                    let bit = index.position_offsets[id] + offset;
                    let hit = seed_at(fragment, offset, index.options.k);
                    let location = hit
                        .start
                        .checked_mul(2)
                        .filter(|&v| v != NO_LOCATION)
                        .ok_or_else(|| invalid("reference position exceeds context encoding"))?
                        | u64::from(hit.reverse);
                    if positions[bit] == NO_LOCATION {
                        positions[bit] = location;
                    } else {
                        let mask = 1 << (bit % 64);
                        if repeated[bit / 64] & mask == 0 {
                            copies.push((bit, positions[bit]));
                            repeated[bit / 64] |= mask;
                        }
                        copies.push((bit, location));
                    }
                }
            }
        }
        copies.sort_unstable();
        copies.dedup();
        let mut ranges = Vec::new();
        let mut groups = Vec::new();
        let mut start = 0;
        while start < copies.len() {
            let bit = copies[start].0;
            let mut end = start + 1;
            while end < copies.len() && copies[end].0 == bit {
                end += 1;
            }
            if end - start == 1 {
                repeated[bit / 64] &= !(1 << (bit % 64));
            } else {
                positions[bit] = groups.len() as u64;
                groups.push(RepeatedPosition {
                    first: copies[start].1,
                    second: copies[start + 1].1,
                    roots: [NONE; 2],
                    pool: 0,
                });
                if end - start > 2 {
                    ranges.push((groups.len() - 1, start, end));
                }
            }
            start = end;
        }
        let (mut reference_repeats, mut reference_repeats_reverse) =
            Self::build_reference_repeat_stops(reference, &copies, &repeated);
        Self::coalesce_multi_copy_runs(
            &copies,
            &positions,
            &repeated,
            groups.len(),
            &mut reference_repeats,
            &mut reference_repeats_reverse,
        );
        let next = AtomicUsize::new(0);
        let workers = index.options.threads.min(ranges.len().max(1));
        let results = thread::scope(|scope| -> io::Result<Vec<_>> {
            let mut handles = Vec::new();
            for _ in 0..workers {
                let ranges = &ranges;
                let copies = &copies;
                let next = &next;
                handles.push(scope.spawn(move || -> io::Result<_> {
                    let mut nodes = Vec::new();
                    let mut roots = Vec::new();
                    loop {
                        let job = next.fetch_add(1, Ordering::Relaxed);
                        if job >= ranges.len() {
                            break;
                        }
                        let (group, start, end) = ranges[job];
                        let locations: Vec<_> = copies[start..end].iter().map(|p| p.1).collect();
                        let forward = refine(
                            reference,
                            locations.clone(),
                            false,
                            index.options.k,
                            index.options.k,
                            &mut nodes,
                        )?;
                        let reverse = refine(
                            reference,
                            locations,
                            true,
                            index.options.k,
                            index.options.k,
                            &mut nodes,
                        )?;
                        roots.push((group, [forward, reverse]));
                    }
                    Ok((nodes, roots))
                }));
            }
            handles
                .into_iter()
                .map(|h| {
                    h.join()
                        .map_err(|_| io::Error::other("proto restriction worker panicked"))?
                })
                .collect()
        })?;
        let mut pools = Vec::new();
        for (pool, (nodes, roots)) in results.into_iter().enumerate() {
            for (group, roots) in roots {
                groups[group].roots = roots;
                groups[group].pool = pool;
            }
            pools.push(nodes);
        }
        eprintln!(
            "[proto_restriction]\t{:.6}s\t{} repeated positions\t{} refined positions\t{} context groups",
            started.elapsed().as_secs_f64(),
            groups.len(),
            ranges.len(),
            pools.iter().map(Vec::len).sum::<usize>()
        );
        Ok(Self {
            positions,
            repeated,
            boundaries,
            reference_repeats,
            reference_repeats_reverse,
            groups,
            pools,
        })
    }

    fn lookup(
        &self,
        reference: &Reference,
        bit: usize,
        direction: bool,
        query: &PackedSeqVec,
        start: usize,
        end: usize,
        k: usize,
    ) -> Option<[u64; 2]> {
        let value = *self.positions.get(bit)?;
        if value == NO_LOCATION {
            return None;
        }
        if self.repeated[bit / 64] & (1 << (bit % 64)) == 0 {
            return Some([value, NO_LOCATION]);
        }
        let group = &self.groups[value as usize];
        let mut node = group.roots[usize::from(direction)];
        if node == NONE {
            return Some([group.first, group.second]);
        }
        let mut known = k;
        loop {
            let item = &self.pools[group.pool][node as usize];
            let depth = item.depth & !LEAF;
            if end - start < depth {
                return None;
            }
            if depth > known {
                let v = seed_view(reference, item.sample, direction, k);
                if query
                    .slice(start + known..start + depth)
                    .cmp_lcp(&v.seq.slice(v.start + known..v.start + depth))
                    .1
                    != depth - known
                {
                    return None;
                }
            }
            if item.depth & LEAF != 0 {
                return Some([item.sample, item.second_location()]);
            }
            if start + depth >= end {
                return None;
            }
            node = item.edges[query.as_slice().get(start + depth) as usize];
            if node == NONE {
                return None;
            }
            known = depth + 1;
        }
    }

    pub(super) fn identity(
        reference: &Reference,
        shards: usize,
        k: usize,
        writer: &mut ContextWriter,
    ) -> io::Result<ScanStats> {
        let mut sink = ContextSink::new(shards, reference, k);
        let mut stats = ScanStats {
            genomes: 1,
            ..Default::default()
        };
        for (id, contig) in reference.contigs.iter().enumerate() {
            for run in &contig.runs {
                let begin = reference.meta.contigs[id].start + run.start as u64;
                sink.push(Record {
                    begin,
                    end: reference.meta.contigs[id].start + run.end as u64,
                    hit: Hit {
                        pattern: 0,
                        contig: id as u32,
                        start: begin,
                        reverse: false,
                    },
                });
                stats.contexts += 1;
            }
        }
        sink.finish_to(writer, 0)?;
        Ok(stats)
    }

    pub(super) fn scan(
        &self,
        index: &FragmentIndex,
        reference: &Reference,
        query: &PackedContig,
        contig: u32,
        origin: u64,
        stats: &mut ScanStats,
        sink: &mut ContextSink,
        present: Option<&mut [u64]>,
    ) -> io::Result<()> {
        sshash_lib::dispatch_on_k!(index.options.k,K => self.scan_k::<K>(index,reference,query,contig,origin,stats,sink,present))
    }

    fn scan_k<const K: usize>(
        &self,
        index: &FragmentIndex,
        reference: &Reference,
        query: &PackedContig,
        contig: u32,
        origin: u64,
        stats: &mut ScanStats,
        sink: &mut ContextSink,
        mut present: Option<&mut [u64]>,
    ) -> io::Result<()>
    where
        Kmer<K>: KmerBits,
    {
        let compact = index.compact.as_ref().unwrap();
        let mut stream = compact.query::<K>(&query.forward);
        let mut cache = MonotonicAlignmentCache::default();
        let mut recent = RecentAlignments::default();
        for run in &query.runs {
            stream.reset();
            let mut pos = run.start;
            while pos + K <= run.end {
                let found = stream.lookup(pos, run.end);
                if !found.is_found() {
                    let next = stream.next_after_miss(pos, &found);
                    stats.absent_windows_skipped += next - pos - 1;
                    pos = next;
                    continue;
                }
                let id = found.string_id as usize;
                let offset = found.kmer_id_in_string as usize;
                let bit = index.position_offsets[id] + offset;
                if let Some(bits) = present.as_deref_mut() {
                    bits[bit / 64] |= 1 << (bit % 64);
                }
                let direction = found.kmer_orientation < 0;
                let remaining = if direction {
                    offset + 1
                } else {
                    index.position_offsets[id + 1] - bit
                };
                let mut advance = if remaining == 1 {
                    1
                } else {
                    let proto = &index.patterns[id];
                    let (sequence, seed) = if direction {
                        (&proto.reverse, proto.len() - offset - K)
                    } else {
                        (&proto.forward, offset)
                    };
                    let right = sequence
                        .slice(seed + K..proto.len())
                        .cmp_lcp(&query.forward.slice(pos + K..run.end))
                        .1;
                    self.constant_span(bit, direction, 1 + right)
                };
                let proto_advance = advance;
                let locations =
                    self.lookup(reference, bit, direction, &query.forward, pos, run.end, K);
                // A failed right-context refinement has no anchor to emit.
                // Its full occurrence set still permits a jump through an
                // unchanged reference class, exactly as a constant proto run.
                // The first differing DNA window and every class change are
                // visited; no arbitrary shared right prefix is skipped.
                if locations.is_none() && present.is_none() && pos + K < run.end {
                    let group = self.positions[bit];
                    if group != NO_LOCATION && self.repeated[bit / 64] & (1 << (bit % 64)) != 0 {
                        let sample = self.groups[group as usize].first;
                        let physical = sample >> 1;
                        let reverse = (sample & 1 != 0) != direction;
                        if (!reverse || physical > 0)
                            && self.unique_alignment_span(physical, reverse, 1) > 1
                        {
                            let v = view(reference, sample, direction, K);
                            let limit = (v.end - v.start - K).min(run.end - pos - K);
                            let span = self.unique_alignment_span(physical, reverse, limit as u64);
                            // Compare only up to this class boundary, avoiding
                            // repeated long LCPs at nearby class changes.
                            let right = v
                                .seq
                                .slice(v.start + K..v.start + K + span - 1)
                                .cmp_lcp(&query.forward.slice(pos + K..pos + K + span - 1))
                                .1;
                            advance = advance.max(right + 1);
                        }
                    }
                }
                for location in locations
                    .into_iter()
                    .flatten()
                    .filter(|&l| l != NO_LOCATION)
                {
                    let physical = location >> 1;
                    let reverse = (location & 1 != 0) != direction;
                    let oriented_query = if reverse { query.len() - pos - K } else { pos };
                    let (slot, diagonal) = RecentAlignments::key(physical, oriented_query, reverse);
                    if recent.covers(slot, diagonal, reverse, physical, K) {
                        stats.cache_hits += 1;
                        stats.hot_cache_hits += 1;
                        if present.is_none() {
                            let entry = &recent.entries[slot];
                            let extension = if reverse {
                                physical - entry.begin
                            } else {
                                entry.end - physical - K as u64
                            };
                            advance = advance
                                .max(self.unique_alignment_span(physical, reverse, extension));
                        }
                        continue;
                    }
                    let v = seed_view(reference, location, direction, K);
                    let anchor = Anchor {
                        hit: Hit {
                            pattern: id as u32,
                            contig: v.contig as u32,
                            start: v.physical,
                            reverse: false,
                        },
                        len: K,
                    };
                    let hit = Hit {
                        pattern: id as u32,
                        contig,
                        start: origin + pos as u64,
                        reverse: v.reverse,
                    };
                    let before = stats.alignments;
                    let (left, right) = alignment_context_monotonic(
                        reference, &anchor, query, hit, origin, &mut cache, stats,
                    )?;
                    recent.remember(
                        slot,
                        diagonal,
                        reverse,
                        physical - left,
                        physical + K as u64 + right,
                    );
                    if present.is_none() {
                        advance = advance.max(self.unique_alignment_span(
                            physical,
                            reverse,
                            if reverse { left } else { right },
                        ));
                    }
                    if stats.alignments == before {
                        continue;
                    }
                    sink.push(Record {
                        begin: v.physical - left,
                        end: v.physical + K as u64 + right,
                        hit: Hit {
                            start: hit.start - if hit.reverse { right } else { left },
                            ..hit
                        },
                    });
                    stats.contexts += 1;
                }
                stats.unique_alignment_windows_skipped += advance - proto_advance;
                if advance > 1 {
                    if let Some(bits) = present.as_deref_mut() {
                        let begin = if direction { bit + 1 - advance } else { bit };
                        mark_position_range(bits, begin, begin + advance);
                    }
                    stream.reset();
                }
                pos += advance;
            }
        }
        Ok(())
    }

    pub(super) fn unique(
        &self,
        index: &FragmentIndex,
        reference: &Reference,
        item: &VerifiedMatch,
    ) -> bool {
        sshash_lib::dispatch_on_k!(index.options.k,K => self.unique_k::<K>(index,reference,item))
    }

    fn unique_k<const K: usize>(
        &self,
        index: &FragmentIndex,
        reference: &Reference,
        item: &VerifiedMatch,
    ) -> bool
    where
        Kmer<K>: KmerBits,
    {
        let hit = item.hits[0];
        let contig = &reference.contigs[hit.contig as usize];
        let start = (hit.start - reference.meta.contigs[hit.contig as usize].start) as usize;
        let len = item.len as usize;
        let seed: Vec<_> = (start..start + K)
            .map(|i| contig.forward.as_slice().get_ascii(i))
            .collect();
        let key = Kmer::<K>::from_ascii_unchecked(&seed);
        let found = index.compact.as_ref().unwrap().dictionary.query::<K>(&key);
        if !found.is_found() {
            return false;
        }
        let bit =
            index.position_offsets[found.string_id as usize] + found.kmer_id_in_string as usize;
        let direction = found.kmer_orientation < 0;
        let Some(locations) = self.lookup(
            reference,
            bit,
            direction,
            &contig.forward,
            start,
            start + len,
            K,
        ) else {
            return false;
        };
        let mut physical = None;
        for location in locations.into_iter().filter(|&l| l != NO_LOCATION) {
            let v = view(reference, location, direction, K);
            if v.end - v.start < len
                || contig
                    .forward
                    .slice(start..start + len)
                    .cmp_lcp(&v.seq.slice(v.start..v.start + len))
                    .1
                    != len
            {
                continue;
            }
            let pos = if v.reverse {
                v.physical + K as u64 - item.len
            } else {
                v.physical
            };
            let candidate = (v.contig, pos);
            if physical.is_some_and(|old| old != candidate) {
                return false;
            }
            physical = Some(candidate);
        }
        physical == Some((hit.contig as usize, hit.start))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn multi_copy_stop_maps_match_full_occurrence_set_oracle() {
        let mut rng = 87654321u64;
        let core: Vec<_> = (0..1200)
            .map(|_| {
                rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
                b"ACGT"[(rng >> 62) as usize]
            })
            .collect();
        let mut broken = core[90..950].to_vec();
        broken[333] = b'N';
        broken[555] = if broken[555] == b'A' { b'C' } else { b'A' };
        let dna = [
            core.clone(),
            core[100..800].to_vec(),
            reverse(&core[100..800]),
            core[150..800].to_vec(),
            core[250..900].to_vec(),
            broken,
        ];
        for k in [3, 7, 31] {
            let reference = reference(&dna, k);
            let bases: usize = dna.iter().map(Vec::len).sum();
            let mut canonical = std::collections::BTreeMap::<Vec<u8>, Vec<u64>>::new();
            for (id, sequence) in dna.iter().enumerate() {
                for (p, window) in sequence.windows(k).enumerate() {
                    if !window.iter().copied().all(is_dna) {
                        continue;
                    }
                    let rc = reverse(window);
                    let strand = window > rc.as_slice();
                    let key = if strand { rc } else { window.to_vec() };
                    let physical = reference.meta.contigs[id].start + p as u64;
                    canonical
                        .entry(key)
                        .or_default()
                        .push(physical * 2 + u64::from(strand));
                }
            }
            let mut positions = vec![NO_LOCATION; canonical.len()];
            let mut repeated = vec![0u64; canonical.len().div_ceil(64)];
            let mut copies = Vec::new();
            let mut groups = Vec::<Vec<u64>>::new();
            let mut physical = std::collections::BTreeMap::new();
            let mut forward = vec![0u64; bases.div_ceil(64)];
            for (bit, locations) in canonical.into_values().enumerate() {
                if locations.len() < 2 {
                    continue;
                }
                let group = groups.len();
                positions[bit] = group as u64;
                repeated[bit / 64] |= 1 << (bit % 64);
                for &location in &locations {
                    copies.push((bit, location));
                    let p = (location >> 1) as usize;
                    assert!(physical.insert(p, (group, location & 1 != 0)).is_none());
                    forward[p / 64] |= 1 << (p % 64);
                }
                groups.push(locations);
            }
            let mut reverse = forward.clone();
            RestrictionIndex::coalesce_multi_copy_runs(
                &copies,
                &positions,
                &repeated,
                groups.len(),
                &mut forward,
                &mut reverse,
            );
            // This oracle enumerates complete affine mapping sets. It does not
            // count graph edges or reuse the construction's normalization.
            let signature = |p: usize| {
                let &(group, own_reverse) = &physical[&p];
                let mut mappings: Vec<_> = groups[group]
                    .iter()
                    .map(|&location| {
                        let q = (location >> 1) as i128;
                        let reverse = (location & 1 != 0) != own_reverse;
                        (
                            reverse,
                            if reverse {
                                p as i128 + q
                            } else {
                                q - p as i128
                            },
                        )
                    })
                    .collect();
                mappings.sort_unstable();
                mappings
            };
            let same_class = |a: usize, b: usize| {
                let Some(&(ga, _)) = physical.get(&a) else {
                    return false;
                };
                let Some(&(gb, _)) = physical.get(&b) else {
                    return false;
                };
                ga != gb
                    && groups[ga].len() > 2
                    && groups[gb].len() > 2
                    && signature(a) == signature(b)
            };
            let mut merged = 0;
            for &p in physical.keys() {
                let expected_forward = p == 0 || !same_class(p - 1, p);
                let expected_reverse = !same_class(p, p + 1);
                assert_eq!(
                    forward[p / 64] & (1 << (p % 64)) != 0,
                    expected_forward,
                    "k={k}, forward p={p}"
                );
                assert_eq!(
                    reverse[p / 64] & (1 << (p % 64)) != 0,
                    expected_reverse,
                    "k={k}, reverse p={p}"
                );
                merged += usize::from(!expected_forward);
            }
            if k >= 7 {
                assert!(merged > 100, "exercise full occurrence-class coalescing");
            }
        }
    }

    #[test]
    fn multi_copy_stop_maps_block_conflicting_physical_groups() {
        let copies = [
            (0, 20),
            (0, 200),
            (0, 400),
            (1, 22),
            (1, 202),
            (1, 402),
            (2, 22),
            (2, 300),
            (2, 500),
        ];
        let mut forward = vec![u64::MAX; 5];
        let mut reverse = forward.clone();
        RestrictionIndex::coalesce_multi_copy_runs(
            &copies,
            &[0, 1, 2],
            &[7],
            3,
            &mut forward,
            &mut reverse,
        );
        assert!(forward.iter().all(|&w| w == u64::MAX));
        assert!(reverse.iter().all(|&w| w == u64::MAX));
    }

    #[test]
    fn multi_copy_run_skips_preserve_all_contexts_and_unresolved_repeats() -> io::Result<()> {
        let temp = tempfile::tempdir()?;
        let mut rng = 99887766u64;
        let core: Vec<_> = (0..1800)
            .map(|_| {
                rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
                b"ACGT"[(rng >> 62) as usize]
            })
            .collect();
        let mut with_tail = core[200..900].to_vec();
        with_tail.extend(&core[1200..1350]);
        let mut nested1 = core[250..900].to_vec();
        nested1.extend(&core[1200..1350]);
        let mut nested2 = core[300..900].to_vec();
        nested2.extend(&core[1200..1350]);
        let mut nested3 = core[350..900].to_vec();
        nested3.extend(&core[1200..1350]);
        let mut broken_reference = core[500..1100].to_vec();
        broken_reference[200] = b'N';
        let references = [
            core.clone(),
            core[200..900].to_vec(),
            reverse(&core[200..900]),
            with_tail,
            nested1,
            nested2,
            nested3,
            broken_reference,
        ];
        let mut crossing = core[150..900].to_vec();
        crossing.extend(&core[1200..1350]);
        let unresolved = core[400..650].to_vec();
        let mut broken_query = core.clone();
        broken_query[500] = b'N';
        broken_query[800] = if broken_query[800] == b'A' {
            b'C'
        } else {
            b'A'
        };
        for k in [7usize, 31, 63] {
            let mut source = ProtoIndex {
                options: LocateOptions {
                    k,
                    m: k.min(19),
                    threads: 2,
                    chunk_bases: 127,
                    forward_only: false,
                },
                patterns: Vec::new(),
                sequences: PackedSeqVec::default(),
            };
            let keys: std::collections::BTreeSet<_> = core
                .windows(k)
                .map(|w| w.to_vec().min(reverse(w)))
                .collect();
            for seed in keys {
                source.add_pattern(&seed)?;
            }
            let index = FragmentIndex::from_graph(source, Some(temp.path()))?;
            let path = temp.path().join("multi.fa");
            let mut fasta = Vec::new();
            for (id, sequence) in references.iter().enumerate() {
                fasta.extend(format!(">c{id}\n").bytes());
                fasta.extend(sequence);
                fasta.push(b'\n');
            }
            std::fs::write(&path, fasta)?;
            let reference = read_reference(&index, &path)?;
            let mut restriction = RestrictionIndex::build(&index, &reference)?;
            let forward = restriction.reference_repeats.clone();
            let reverse_stops = restriction.reference_repeats_reverse.clone();
            let queries = [
                core.clone(),
                reverse(&core),
                crossing.clone(),
                reverse(&crossing),
                unresolved.clone(),
                reverse(&unresolved),
                broken_query.clone(),
            ];
            for (query_id, sequence) in queries.iter().enumerate() {
                let query = PackedContig::new(sequence, k);
                let mut outcomes = Vec::new();
                for exhaustive in [false, true] {
                    restriction.reference_repeats.clone_from(&forward);
                    restriction
                        .reference_repeats_reverse
                        .clone_from(&reverse_stops);
                    if exhaustive {
                        restriction.reference_repeats.fill(u64::MAX);
                        restriction.reference_repeats_reverse.fill(u64::MAX);
                    }
                    let mut sink = ContextSink::new(3, &reference, k);
                    let mut stats = ScanStats::default();
                    restriction.scan(
                        &index, &reference, &query, 4, 91, &mut stats, &mut sink, None,
                    )?;
                    let mut rows: Vec<_> = sink
                        .shards
                        .into_iter()
                        .flatten()
                        .map(|r| (r.begin, r.end, r.hit.contig, r.hit.start, r.hit.reverse))
                        .collect();
                    rows.sort_unstable();
                    rows.dedup();
                    if (4..=5).contains(&query_id) {
                        assert!(rows.is_empty(), "unresolved copies cannot become anchors");
                        if !exhaustive && k >= 31 {
                            assert!(
                                stats.unique_alignment_windows_skipped > 100,
                                "exercise unresolved class jump"
                            );
                        }
                    }
                    outcomes.push(rows);
                }
                assert_eq!(outcomes[0], outcomes[1], "k={k}, query={query_id}");
            }
        }
        Ok(())
    }

    #[test]
    fn indexed_seed_views_preserve_physical_coordinates_strands_and_dna() {
        let mut rng = 1234567u64;
        let long: Vec<_> = (0..140)
            .map(|_| {
                rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
                b"ACGT"[(rng >> 62) as usize]
            })
            .collect();
        let dna = [
            b"NNACGTactgGATTACANNTTGCAACGTTNN".to_vec(),
            vec![b'N'],
            vec![],
            long,
        ];
        for k in [3, 7, 31] {
            let reference = reference(&dna, k);
            for (id, sequence) in dna.iter().enumerate() {
                for (local, window) in sequence.windows(k).enumerate() {
                    if !window.iter().copied().all(is_dna) {
                        continue;
                    }
                    let physical = reference.meta.contigs[id].start + local as u64;
                    for canonical_reverse in [false, true] {
                        for direction in [false, true] {
                            let location = physical * 2 + u64::from(canonical_reverse);
                            let seed = seed_view(&reference, location, direction, k);
                            assert_eq!(seed.contig, id);
                            assert_eq!(seed.physical, physical);
                            assert_eq!(seed.reverse, canonical_reverse != direction);
                            let upper: Vec<_> = window.iter().map(u8::to_ascii_uppercase).collect();
                            let expected = if seed.reverse { reverse(&upper) } else { upper };
                            let actual: Vec<_> = (seed.start..seed.start + k)
                                .map(|p| seed.seq.as_slice().get_ascii(p))
                                .collect();
                            assert_eq!(actual, expected);
                            let full = view(&reference, location, direction, k);
                            assert_eq!(full.start, seed.start);
                            assert!(full.end - full.start >= k);
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn compact_context_leaves_preserve_wide_locations_and_absence() {
        for sample in [0, 1, u32::MAX as u64, (1u64 << 32) + 1, u64::MAX - 1] {
            for second in [0, 1, u32::MAX as u64, 1u64 << 32, u64::MAX - 1, NO_LOCATION] {
                for depth in [0, 63, (usize::MAX >> 2), LEAF - 1] {
                    let leaf = ContextGroup::leaf(depth, sample, second);
                    assert_eq!(leaf.sample, sample);
                    assert_eq!(leaf.second_location(), second);
                    assert_eq!(leaf.depth & !LEAF, depth);
                    assert_ne!(leaf.depth & LEAF, 0);
                }
            }
        }
    }

    fn reverse(sequence: &[u8]) -> Vec<u8> {
        sequence
            .iter()
            .rev()
            .map(|base| match base {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                _ => b'N',
            })
            .collect()
    }

    #[test]
    fn directional_repeat_stops_preserve_mapping_entries_and_changes() {
        let reference = reference(&[vec![b'A'; 400]], 3);
        let locations = [
            vec![(10u64, false), (100, false)],
            vec![(11, false), (101, false)],
            vec![(12, true), (102, true)],
            vec![(20, false), (180, true)],
            vec![(21, true), (179, false)],
            vec![(22, false), (178, true), (250, false)],
            vec![(30, false), (130, false)],
            vec![(31, false), (131, true)],
            vec![(40, false), (140, false)],
            vec![(41, false), (142, false)],
            vec![(50, false), (50, true)],
        ];
        let copies: Vec<_> = locations
            .iter()
            .enumerate()
            .flat_map(|(bit, entries)| {
                entries
                    .iter()
                    .map(move |&(p, reverse)| (bit, p * 2 + u64::from(reverse)))
            })
            .collect();
        let (forward, reverse) =
            RestrictionIndex::build_reference_repeat_stops(&reference, &copies, &[u64::MAX]);
        let is_set = |words: &[u64], p: usize| words[p / 64] & (1 << (p % 64)) != 0;
        for p in [
            10, 100, 20, 179, 22, 178, 250, 30, 31, 130, 131, 40, 41, 140, 142, 50,
        ] {
            assert!(is_set(&forward, p), "forward entry/change {p}");
        }
        for p in [11, 12, 101, 102, 21, 180] {
            assert!(!is_set(&forward, p), "forward stable interior {p}");
        }
        for p in [
            12, 102, 21, 180, 22, 178, 250, 30, 31, 130, 131, 40, 41, 140, 142, 50,
        ] {
            assert!(is_set(&reverse, p), "reverse entry/change {p}");
        }
        for p in [10, 11, 100, 101, 20, 179] {
            assert!(!is_set(&reverse, p), "reverse stable interior {p}");
        }
        // Multiple proto groups at a physical position remain conservative.
        let mut overlapping = copies;
        overlapping.extend([(locations.len(), 11 * 2), (locations.len(), 300 * 2)]);
        let (forward, reverse) =
            RestrictionIndex::build_reference_repeat_stops(&reference, &overlapping, &[u64::MAX]);
        for p in [10, 11, 12] {
            assert!(is_set(&forward, p));
            assert!(is_set(&reverse, p));
        }
    }

    #[test]
    fn physical_unique_spans_stop_before_every_repeat_on_both_strands() {
        let mut rng = 0x1234_abcd_9876_5ef0u64;
        for count in [1usize, 2, 63, 64, 65, 127, 128, 129, 513] {
            let mut flags = vec![false; count];
            let mut words = vec![0u64; count.div_ceil(64)];
            for (position, flag) in flags.iter_mut().enumerate() {
                rng ^= rng << 13;
                rng ^= rng >> 7;
                rng ^= rng << 17;
                *flag = rng % 11 == 0;
                if *flag {
                    words[position / 64] |= 1 << (position % 64);
                }
            }
            let reverse_flags: Vec<_> = flags.iter().rev().copied().collect();
            let mut reverse_words = vec![0u64; words.len()];
            for (p, &flag) in reverse_flags.iter().enumerate() {
                if flag {
                    reverse_words[p / 64] |= 1 << (p % 64);
                }
            }
            let index = RestrictionIndex {
                positions: Vec::new(),
                repeated: Vec::new(),
                boundaries: Vec::new(),
                reference_repeats_reverse: reverse_words,
                reference_repeats: words,
                groups: Vec::new(),
                pools: Vec::new(),
            };
            for position in 0..count {
                for reverse in [false, true] {
                    let flags = if reverse { &reverse_flags } else { &flags };
                    let limit = if reverse {
                        position + 1
                    } else {
                        count - position
                    };
                    for maximum in 1..=limit {
                        let expected = (1..maximum)
                            .find(|&step| {
                                flags[if reverse {
                                    position - step
                                } else {
                                    position + step
                                }]
                            })
                            .unwrap_or(maximum);
                        assert_eq!(
                            index.unique_alignment_span(
                                position as u64,
                                reverse,
                                (maximum - 1) as u64
                            ),
                            expected
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn cross_proto_alignment_skips_preserve_all_contexts_with_repeats_and_ambiguity()
    -> io::Result<()> {
        let temp = tempfile::tempdir()?;
        let mut rng = 1234567u64;
        let core: Vec<_> = (0..3000)
            .map(|_| {
                rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
                b"ACGT"[(rng >> 62) as usize]
            })
            .collect();
        let mut broken = core.clone();
        broken[1300] = b'N';
        broken[2300] = if broken[2300] == b'A' { b'C' } else { b'A' };
        let mut rearranged = core[500..1700].to_vec();
        rearranged.extend(reverse(&core[1800..2700]));
        let inside_repeat = core[350..900].to_vec();
        for k in [7usize, 31, 33, 63] {
            let mut source = ProtoIndex {
                options: LocateOptions {
                    k,
                    m: k.min(19),
                    threads: 2,
                    chunk_bases: 127,
                    forward_only: false,
                },
                patterns: Vec::new(),
                sequences: PackedSeqVec::default(),
            };
            // Individual canonical k-mers force every skipped step to cross a
            // proto boundary. No longer proto can justify the new jumps.
            let mut keys = std::collections::BTreeSet::new();
            for seed in core.windows(k) {
                keys.insert(seed.to_vec().min(reverse(seed)));
            }
            for seed in keys {
                source.add_pattern(&seed)?;
            }
            let index = FragmentIndex::from_graph(source, Some(temp.path()))?;
            let path = temp.path().join("reference.fa");
            let mut fasta = b">core\n".to_vec();
            fasta.extend(&core);
            fasta.extend(b"\n>repeat\n");
            fasta.extend(&core[200..700]);
            fasta.extend(b"\n>reverse-repeat\n");
            fasta.extend(reverse(&core[800..1300]));
            fasta.push(b'\n');
            std::fs::write(&path, fasta)?;
            let reference = read_reference(&index, &path)?;
            let mut restriction = RestrictionIndex::build(&index, &reference)?;
            if k >= 31 {
                let is_set = |words: &[u64], p: usize| words[p / 64] & (1 << (p % 64)) != 0;
                assert!(is_set(&restriction.reference_repeats, 200));
                assert!(!is_set(&restriction.reference_repeats, 201));
                assert!(is_set(&restriction.reference_repeats_reverse, 700 - k));
                assert!(!is_set(&restriction.reference_repeats_reverse, 699 - k));
            }
            let mut total_skipped = 0;
            for sequence in [
                &core,
                &reverse(&core),
                &broken,
                &rearranged,
                &inside_repeat,
                &reverse(&inside_repeat),
            ] {
                let query = PackedContig::new(sequence, k);
                let mut rows = Vec::new();
                let repeats = restriction.reference_repeats.clone();
                let reverse_repeats = restriction.reference_repeats_reverse.clone();
                for exhaustive in [false, true] {
                    if exhaustive {
                        restriction.reference_repeats.fill(u64::MAX);
                        restriction.reference_repeats_reverse.fill(u64::MAX);
                    }
                    let mut sink = ContextSink::new(3, &reference, k);
                    let mut stats = ScanStats::default();
                    restriction.scan(
                        &index, &reference, &query, 4, 91, &mut stats, &mut sink, None,
                    )?;
                    if !exhaustive {
                        total_skipped += stats.unique_alignment_windows_skipped;
                    }
                    let mut keys: Vec<_> = sink
                        .shards
                        .into_iter()
                        .flatten()
                        .map(|r| (r.begin, r.end, r.hit.contig, r.hit.start, r.hit.reverse))
                        .collect();
                    keys.sort_unstable();
                    keys.dedup();
                    rows.push(keys);
                }
                restriction.reference_repeats = repeats;
                restriction.reference_repeats_reverse = reverse_repeats;
                assert_eq!(rows[0], rows[1], "k={k}, query length={}", sequence.len());
            }
            assert!(
                total_skipped > 1000,
                "new cross-proto skip must be exercised for k={k}"
            );
        }
        Ok(())
    }

    #[test]
    fn first_level_alignment_cache_checks_exact_mapping_bounds_strand_and_collisions() {
        for physical in [1000, u64::MAX - 200] {
            for query in [128, 123456, usize::MAX - 256] {
                for reverse in [false, true] {
                    let mut cache = RecentAlignments::default();
                    let (slot, diagonal) = RecentAlignments::key(physical, query, reverse);
                    assert!(!cache.covers(slot, diagonal, reverse, physical, 63));
                    cache.remember(slot, diagonal, reverse, physical - 5, physical + 80);
                    for delta in -7i128..=85 {
                        let r = (physical as i128 + delta) as u64;
                        let q = (query as i128 + delta) as usize;
                        let (target, key) = RecentAlignments::key(r, q, reverse);
                        assert_eq!(
                            cache.covers(target, key, reverse, r, 63),
                            delta >= -5 && delta + 63 <= 80
                        );
                        assert!(!cache.covers(target, key, !reverse, r, 63));
                    }
                    // Force a different key into the same slot. Its overlapping
                    // interval must not validate the displaced mapping.
                    cache.remember(slot, diagonal + 1, reverse, physical - 5, physical + 80);
                    assert!(!cache.covers(slot, diagonal, reverse, physical, 63));
                    cache.remember(slot, diagonal, !reverse, physical - 5, physical + 80);
                    assert!(!cache.covers(slot, diagonal, reverse, physical, 63));
                }
            }
        }
    }

    #[test]
    fn constant_occurrence_spans_stop_at_every_boundary_in_both_directions() {
        let mut state = 0x1234_abcd_9876_5ef0u64;
        for count in [1usize, 2, 63, 64, 65, 127, 128, 129, 257] {
            let mut flags = vec![false; count];
            let mut boundaries = vec![0u64; count.div_ceil(64)];
            for bit in 0..count {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                flags[bit] = state % 7 == 0;
                if flags[bit] {
                    boundaries[bit / 64] |= 1 << (bit % 64);
                }
            }
            let index = RestrictionIndex {
                positions: vec![NO_LOCATION; count],
                repeated: vec![0; count.div_ceil(64)],
                boundaries,
                reference_repeats: Vec::new(),
                reference_repeats_reverse: Vec::new(),
                groups: Vec::new(),
                pools: Vec::new(),
            };
            for bit in 0..count {
                for backwards in [false, true] {
                    let limit = if backwards { bit + 1 } else { count - bit };
                    for maximum in 1..=limit.min(65) {
                        let expected = (1..maximum)
                            .find(|&step| {
                                flags[if backwards {
                                    bit + 1 - step
                                } else {
                                    bit + step
                                }]
                            })
                            .unwrap_or(maximum);
                        assert_eq!(index.constant_span(bit, backwards, maximum), expected);
                    }
                }
            }
        }
    }

    fn reference(contigs: &[Vec<u8>], k: usize) -> Reference {
        let mut origin = 0;
        let meta = contigs
            .iter()
            .enumerate()
            .map(|(i, seq)| {
                let record = ContigMeta {
                    name: format!("c{i}"),
                    start: origin,
                    len: seq.len(),
                };
                origin += seq.len() as u64;
                record
            })
            .collect();
        Reference {
            contigs: contigs
                .iter()
                .map(|seq| PackedContig::new(seq, k))
                .collect(),
            meta: GenomeMeta {
                input: "reference".into(),
                contigs: meta,
            },
            fragments: Vec::new(),
        }
    }

    #[test]
    fn nested_overlapping_repeat_contexts_compress_without_losing_short_boundaries()
    -> io::Result<()> {
        let reference = reference(&[vec![b'A'; 20_000]], 3);
        let locations: Vec<_> = (0..=19_997).map(|i| (i as u64) * 2).collect();
        for direction in [false, true] {
            let mut nodes = Vec::new();
            let root = refine(&reference, locations.clone(), direction, 3, 3, &mut nodes)?;
            assert_eq!(nodes.len(), 1);
            assert_ne!(nodes[0].depth & LEAF, 0);
            assert_eq!(nodes[0].depth & !LEAF, 19_999);
            let index = RestrictionIndex {
                positions: vec![0],
                repeated: vec![1],
                boundaries: vec![1],
                reference_repeats: Vec::new(),
                reference_repeats_reverse: Vec::new(),
                groups: vec![RepeatedPosition {
                    first: locations[0],
                    second: locations[1],
                    roots: [root; 2],
                    pool: 0,
                }],
                pools: vec![nodes],
            };
            let query =
                PackedSeqVec::from_ascii(&vec![if direction { b'T' } else { b'A' }; 20_000]);
            for len in [3, 17, 19_997, 19_998, 19_999, 20_000] {
                assert_eq!(
                    index
                        .lookup(&reference, 0, direction, &query, 0, len, 3)
                        .is_some(),
                    len >= 19_999
                );
            }
        }
        Ok(())
    }

    #[test]
    fn identical_terminated_copies_remain_ambiguous_without_a_chain() -> io::Result<()> {
        let sequence = b"ACGATCTGACGCTAGCTCGATCGATCTAGC".to_vec();
        let reference = reference(&[sequence.clone(), sequence.clone(), sequence.clone()], 3);
        let locations: Vec<_> = reference.meta.contigs.iter().map(|c| c.start * 2).collect();
        let mut nodes = Vec::new();
        let root = refine(&reference, locations.clone(), false, 3, 3, &mut nodes)?;
        assert_eq!(nodes.len(), 1);
        assert_eq!(nodes[0].depth & LEAF, 0);
        assert_eq!(nodes[0].edges, [NONE; 4]);
        let index = RestrictionIndex {
            positions: vec![0],
            repeated: vec![1],
            boundaries: vec![1],
            reference_repeats: Vec::new(),
            reference_repeats_reverse: Vec::new(),
            groups: vec![RepeatedPosition {
                first: locations[0],
                second: locations[1],
                roots: [root; 2],
                pool: 0,
            }],
            pools: vec![nodes],
        };
        let query = PackedSeqVec::from_ascii(&sequence);
        for len in 3..=sequence.len() {
            assert!(
                index
                    .lookup(&reference, 0, false, &query, 0, len, 3)
                    .is_none()
            );
        }
        Ok(())
    }
}
