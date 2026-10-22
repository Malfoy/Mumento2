//! Repeat evidence belongs to individual proto positions. SSHash remains the
//! only k-mer dictionary. Repeated occurrences are partitioned by exact right
//! context until at most two oriented copies remain; no occurrence cross product
//! or suffix array is constructed. Two copies retain physical palindromes.
use super::*;
use sshash_lib::{Kmer, KmerBits};

const NONE: u32 = u32::MAX;
const NO_LOCATION: u64 = u64::MAX;
const SINGLETON_EXTRA: usize = 64;
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
    leaves: Range<usize>,
}

#[derive(Clone, Copy, Debug)]
struct LeafOccurrence {
    location: u64,
    leaves: [u32; 2],
}

// A bounded rank directory replaces searches of the per-seed leaf lists.
// Entries still identify an original full-key record in its original group.
// Sparse/wide coordinate spaces and conflicting keys retain binary search.
struct LeafRank {
    words: Vec<u64>,
    ranks: Vec<u32>,
    records: Vec<u32>,
}
impl LeafRank {
    fn build(records: &[LeafOccurrence]) -> Option<Self> {
        if records.is_empty() || records.len() > u32::MAX as usize {
            return None;
        }
        let maximum = usize::try_from(records.iter().map(|r| r.location).max()?).ok()?;
        let count = (maximum / 64).checked_add(1)?;
        // At most 64 MiB for words/prefix counts, and bounded density. This
        // optional accelerator must not allocate across large coordinate gaps.
        if count > (64 * 1024 * 1024) / 12 || count > records.len().saturating_mul(8).max(1024) {
            return None;
        }
        let mut ordered: Vec<u32> = (0..records.len() as u32).collect();
        ordered.sort_unstable_by_key(|&i| records[i as usize].location);
        let mut words = vec![0u64; count];
        for &i in &ordered {
            let location = records[i as usize].location as usize;
            let mask = 1u64 << (location % 64);
            let word = &mut words[location / 64];
            if *word & mask != 0 {
                return None;
            }
            *word |= mask;
        }
        let mut total = 0u32;
        let ranks = words
            .iter()
            .map(|&word| {
                let prefix = total;
                total += word.count_ones();
                prefix
            })
            .collect();
        Some(Self {
            words,
            ranks,
            records: ordered,
        })
    }

    #[inline]
    fn lookup(&self, location: u64) -> Option<usize> {
        let location = usize::try_from(location).ok()?;
        let word_index = location / 64;
        let mask = 1u64 << (location % 64);
        let word = *self.words.get(word_index)?;
        if word & mask == 0 {
            return None;
        }
        let rank = self.ranks[word_index] + (word & (mask - 1)).count_ones();
        Some(self.records[rank as usize] as usize)
    }
}
#[derive(Clone, Copy, Debug)]
struct ReferenceSeed {
    pattern: u32,
    offset: u32,
    reverse: bool,
}

// Query-only reference coordinate directory. Every stored seed comes from an
// original exact reference fragment; SSHash remains the fallback dictionary.
// Dense Vec under the 512 MiB cap; otherwise a sparse map covers only fragment
// coordinates so whole-genome references still get prefetch.
enum SeedStore {
    Dense(Vec<u32>),
    Sparse(std::collections::HashMap<u64, u32>),
}

struct ReferenceSeeds {
    store: SeedStore,
    offset_bits: u32,
    offset_mask: u32,
}
impl ReferenceSeeds {
    const EMPTY: u32 = u32::MAX;
    const UNKNOWN: u32 = u32::MAX - 1;
    /// Soft cap for sparse maps (~512 MiB at ~16 B/entry including HashMap overhead).
    const SPARSE_MAX_ENTRIES: usize = (512 * 1024 * 1024) / 16;

    fn pack_value(id: usize, offset: usize, reverse: bool, offset_bits: u32, offset_mask: u32) -> u32 {
        if offset > offset_mask as usize {
            Self::UNKNOWN
        } else {
            let packed =
                (((id as u32) << offset_bits) | offset as u32) << 1 | u32::from(reverse);
            if packed >= Self::UNKNOWN {
                Self::UNKNOWN
            } else {
                packed
            }
        }
    }

    fn insert_dense(entries: &mut [u32], physical: u64, value: u32) -> Option<()> {
        let entry = entries.get_mut(usize::try_from(physical).ok()?)?;
        if *entry == Self::EMPTY {
            *entry = value;
        } else if *entry != value {
            *entry = Self::UNKNOWN;
        }
        Some(())
    }

    fn insert_sparse(map: &mut std::collections::HashMap<u64, u32>, physical: u64, value: u32) {
        match map.get(&physical).copied() {
            None => {
                map.insert(physical, value);
            }
            Some(old) if old != value => {
                map.insert(physical, Self::UNKNOWN);
            }
            Some(_) => {}
        }
    }

    fn build(index: &FragmentIndex, reference: &Reference) -> Option<Self> {
        let count = index.patterns.len();
        if count == 0 || reference.fragments.iter().all(Vec::is_empty) {
            return None;
        }
        let id_bits = (usize::BITS - (count - 1).leading_zeros()).max(1);
        if id_bits >= 31 {
            return None;
        }
        let offset_bits = 31 - id_bits;
        let offset_mask = (1u32 << offset_bits) - 1;
        let last = reference.meta.contigs.last()?;
        let dense_ok = last
            .start
            .checked_add(last.len as u64)
            .and_then(|end| usize::try_from(end).ok())
            .is_some_and(|bases| bases <= (512 * 1024 * 1024) / 4);

        if dense_ok {
            let bases = usize::try_from(last.start + last.len as u64).ok()?;
            let mut entries = Vec::new();
            entries.try_reserve_exact(bases).ok()?;
            entries.resize(bases, Self::EMPTY);
            for (id, fragments) in reference.fragments.iter().enumerate() {
                for &fragment in fragments {
                    if fragment.end < fragment.begin
                        || fragment.end - fragment.begin < index.options.k
                    {
                        continue;
                    }
                    for offset in fragment.begin..=fragment.end - index.options.k {
                        let hit = seed_at(fragment, offset, index.options.k);
                        let value =
                            Self::pack_value(id, offset, hit.reverse, offset_bits, offset_mask);
                        Self::insert_dense(&mut entries, hit.start, value)?;
                    }
                }
            }
            return Some(Self {
                store: SeedStore::Dense(entries),
                offset_bits,
                offset_mask,
            });
        }

        let mut map = std::collections::HashMap::new();
        for (id, fragments) in reference.fragments.iter().enumerate() {
            for &fragment in fragments {
                if fragment.end < fragment.begin || fragment.end - fragment.begin < index.options.k
                {
                    continue;
                }
                let seeds = fragment.end - fragment.begin - index.options.k + 1;
                if map.len().saturating_add(seeds) > Self::SPARSE_MAX_ENTRIES {
                    return None;
                }
                for offset in fragment.begin..=fragment.end - index.options.k {
                    let hit = seed_at(fragment, offset, index.options.k);
                    let value =
                        Self::pack_value(id, offset, hit.reverse, offset_bits, offset_mask);
                    Self::insert_sparse(&mut map, hit.start, value);
                    if map.len() > Self::SPARSE_MAX_ENTRIES {
                        return None;
                    }
                }
            }
        }
        if map.is_empty() {
            return None;
        }
        Some(Self {
            store: SeedStore::Sparse(map),
            offset_bits,
            offset_mask,
        })
    }

    #[inline]
    fn lookup(&self, index: &FragmentIndex, physical: u64, reverse: bool) -> Option<ReferenceSeed> {
        let value = match &self.store {
            SeedStore::Dense(entries) => *entries.get(usize::try_from(physical).ok()?)?,
            SeedStore::Sparse(map) => *map.get(&physical)?,
        };
        if value >= Self::UNKNOWN {
            return None;
        }
        let payload = value >> 1;
        let id = (payload >> self.offset_bits) as usize;
        let offset = (payload & self.offset_mask) as usize;
        let proto = index.patterns.get(id)?;
        if offset + index.options.k > proto.len() {
            return None;
        }
        Some(ReferenceSeed {
            pattern: id as u32,
            offset: offset as u32,
            reverse: (value & 1 != 0) != reverse,
        })
    }

    #[cfg(test)]
    fn dense_entry(&self, physical: usize) -> Option<u32> {
        match &self.store {
            SeedStore::Dense(entries) => entries.get(physical).copied(),
            SeedStore::Sparse(_) => None,
        }
    }
}

pub(super) struct RestrictionIndex {
    positions: Vec<u64>,
    repeated: Vec<u64>,
    boundaries: Vec<u64>,
    // Direction-specific repeat entry points. A verified alignment can cross
    // proto boundaries through unique seeds and complete oriented occurrence
    // coverage. Copies may disappear; each new mapping remains a stopping point.
    reference_repeats: Vec<u64>,
    reference_repeats_reverse: Vec<u64>,
    context_singletons_forward: Vec<u64>,
    context_singletons_reverse: Vec<u64>,
    groups: Vec<RepeatedPosition>,
    pools: Vec<Vec<ContextGroup>>,
    leaf_occurrences: Vec<LeafOccurrence>,
    leaf_rank: Option<LeafRank>,
    reference_seeds: Option<ReferenceSeeds>,
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
fn collect_leaf_occurrences(
    nodes: &[ContextGroup],
    ranges: [Range<usize>; 2],
) -> Vec<LeafOccurrence> {
    let mut entries = Vec::new();
    for (direction, range) in ranges.into_iter().enumerate() {
        for node in range {
            let item = &nodes[node];
            if item.depth & LEAF == 0 {
                continue;
            }
            for location in [item.sample, item.second_location()] {
                if location != NO_LOCATION {
                    entries.push((location, direction, node as u32));
                }
            }
        }
    }
    entries.sort_unstable();
    let mut result = Vec::new();
    let mut begin = 0;
    while begin < entries.len() {
        let location = entries[begin].0;
        let mut end = begin;
        let mut leaves = [NONE; 2];
        let mut conflict = [false; 2];
        while end < entries.len() && entries[end].0 == location {
            let (_, direction, node) = entries[end];
            conflict[direction] |= leaves[direction] != NONE && leaves[direction] != node;
            leaves[direction] = node;
            end += 1;
        }
        for direction in 0..2 {
            if conflict[direction] {
                leaves[direction] = NONE;
            }
        }
        if leaves != [NONE; 2] {
            result.push(LeafOccurrence { location, leaves });
        }
        begin = end;
    }
    result
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
    // either distinct endpoint group. Matching the target group's complete
    // count proves every target copy continues from a source copy. A scan may
    // cross losses, but must stop before any new oriented mapping enters.
    // This uses full-key equality, not an occurrence-set fingerprint.
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
        // Include pairs too: transitions between repeat multiplicities can
        // preserve coverage in one direction. Conflicting physical entries
        // cannot establish an injective correspondence and stay blocked.
        entries.retain(|&(_, encoded)| !blocked[encoded >> 1]);
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
            if count == counts[b.1 >> 1] {
                let mask = 1 << (b.0 % 64);
                coalesced += usize::from(forward[b.0 / 64] & mask != 0);
                forward[b.0 / 64] &= !mask;
            }
            if count == counts[a.1 >> 1] {
                reverse[a.0 / 64] &= !(1 << (a.0 % 64));
            }
        }
        eprintln!(
            "[multi_copy_reference_runs]\t{} reference copies\t{} oriented edges\t{} covered forward transitions",
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
        let stops = if reverse {
            &self.reference_repeats_reverse
        } else {
            &self.reference_repeats
        };
        Self::stop_span(stops, physical, reverse, extension)
    }

    // Every removed stop has an original singleton leaf requiring at most
    // K+64 exact right-context bases. Only use it inside an already verified
    // reference alignment, retaining 64 bases beyond the skipped seed.
    fn known_alignment_span(&self, physical: u64, reverse: bool, extension: u64) -> usize {
        if extension <= SINGLETON_EXTRA as u64 {
            return self.unique_alignment_span(physical, reverse, extension);
        }
        let stops = if reverse {
            &self.context_singletons_reverse
        } else {
            &self.context_singletons_forward
        };
        if stops.is_empty() {
            return self.unique_alignment_span(physical, reverse, extension);
        }
        let proven = extension - SINGLETON_EXTRA as u64;
        let span = Self::stop_span(stops, physical, reverse, proven);
        if span as u64 <= proven {
            return span;
        }
        // The relaxed scan checked every seed through this endpoint. Beyond
        // it the 64-base context reserve is exhausted, so use original stops.
        let boundary = if reverse {
            physical - proven
        } else {
            physical + proven
        };
        proven as usize + self.unique_alignment_span(boundary, reverse, SINGLETON_EXTRA as u64)
    }

    fn stop_span(stops: &[u64], physical: u64, reverse: bool, extension: u64) -> usize {
        let position = physical as usize;
        let maximum = extension as usize + 1;
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

    fn build_context_singleton_stops(
        positions: &[u64],
        repeated: &[u64],
        copies: &[(usize, u64)],
        groups: &[RepeatedPosition],
        pools: &[Vec<ContextGroup>],
        records: &[LeafOccurrence],
        forward: &[u64],
        reverse: &[u64],
        k: usize,
    ) -> (Vec<u64>, Vec<u64>) {
        let mut seen = vec![0u64; forward.len()];
        let mut blocked = seen.clone();
        let mut observe = |location: u64| {
            let physical = (location >> 1) as usize;
            let word = physical / 64;
            if word >= seen.len() {
                return;
            }
            let mask = 1u64 << (physical % 64);
            blocked[word] |= seen[word] & mask;
            seen[word] |= mask;
        };
        // Multiple proto/strand entries at one physical seed cannot justify
        // this shortcut. Include unique and repeated indexed positions.
        for (bit, &location) in positions.iter().enumerate() {
            if location != NO_LOCATION && repeated[bit / 64] & (1 << (bit % 64)) == 0 {
                observe(location);
            }
        }
        for &(_, location) in copies {
            observe(location);
        }
        drop(observe);
        let mut result = [forward.to_vec(), reverse.to_vec()];
        let mut removed = [0usize; 2];
        for group in groups {
            for record in &records[group.leaves.clone()] {
                let physical = (record.location >> 1) as usize;
                let word = physical / 64;
                let mask = 1u64 << (physical % 64);
                if blocked.get(word).is_none_or(|&bits| bits & mask != 0) {
                    continue;
                }
                for reverse in [false, true] {
                    let direction = (record.location & 1 != 0) != reverse;
                    let node = record.leaves[usize::from(direction)];
                    if node == NONE {
                        continue;
                    }
                    let leaf = &pools[group.pool][node as usize];
                    if leaf.depth & LEAF != 0
                        && (leaf.depth & !LEAF) <= k + SINGLETON_EXTRA
                        && leaf.sample == record.location
                        && leaf.second_location() == NO_LOCATION
                    {
                        let flags = &mut result[usize::from(reverse)];
                        removed[usize::from(reverse)] += usize::from(flags[word] & mask != 0);
                        flags[word] &= !mask;
                    }
                }
            }
        }
        eprintln!(
            "[context_singleton_stops]\t{} forward\t{} reverse stops covered by exact <=k+64 singleton leaves",
            removed[0], removed[1]
        );
        let [forward, reverse] = result;
        (forward, reverse)
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
                    leaves: 0..0,
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
                        let forward_begin = nodes.len();
                        let forward = refine(
                            reference,
                            locations.clone(),
                            false,
                            index.options.k,
                            index.options.k,
                            &mut nodes,
                        )?;
                        let reverse_begin = nodes.len();
                        let reverse = refine(
                            reference,
                            locations,
                            true,
                            index.options.k,
                            index.options.k,
                            &mut nodes,
                        )?;
                        let leaves = collect_leaf_occurrences(
                            &nodes,
                            [forward_begin..reverse_begin, reverse_begin..nodes.len()],
                        );
                        roots.push((group, [forward, reverse], leaves));
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
        let mut leaf_occurrences = Vec::new();
        for (pool, (nodes, roots)) in results.into_iter().enumerate() {
            for (group, roots, leaves) in roots {
                let begin = leaf_occurrences.len();
                leaf_occurrences.extend(leaves);
                groups[group].leaves = begin..leaf_occurrences.len();
                groups[group].roots = roots;
                groups[group].pool = pool;
            }
            pools.push(nodes);
        }
        let (context_singletons_forward, context_singletons_reverse) =
            Self::build_context_singleton_stops(
                &positions,
                &repeated,
                &copies,
                &groups,
                &pools,
                &leaf_occurrences,
                &reference_repeats,
                &reference_repeats_reverse,
                index.options.k,
            );
        let reference_seeds = ReferenceSeeds::build(index, reference);
        eprintln!(
            "[reference_seed_directory]\t{} bytes",
            reference_seeds.as_ref().map_or(0, |map| match &map.store {
                SeedStore::Dense(entries) => entries.len() * 4,
                SeedStore::Sparse(h) => h.capacity() * 16,
            })
        );
        let leaf_rank = LeafRank::build(&leaf_occurrences);
        eprintln!(
            "[context_leaf_rank]\t{} leaf occurrences\t{} rank bytes",
            leaf_occurrences.len(),
            leaf_rank.as_ref().map_or(0, |rank| rank.words.len() * 8
                + rank.ranks.len() * 4
                + rank.records.len() * 4)
        );
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
            context_singletons_forward,
            context_singletons_reverse,
            groups,
            pools,
            leaf_occurrences,
            leaf_rank,
            reference_seeds,
        })
    }

    // A previous exact alignment proves the entire right context to a leaf.
    // The full encoded occurrence key identifies that same leaf; otherwise the
    // original DNA-checked trie lookup remains authoritative.
    fn witnessed_leaf(
        &self,
        group: &RepeatedPosition,
        direction: bool,
        start: usize,
        end: usize,
        query_len: usize,
        k: usize,
        witness: RecentAlignment,
    ) -> Option<[u64; 2]> {
        let oriented = if witness.reverse {
            query_len.checked_sub(start.checked_add(k)?)?
        } else {
            start
        };
        let physical = u64::try_from((oriented as i128).checked_sub(witness.diagonal)?).ok()?;
        let seed_end = physical.checked_add(k as u64)?;
        if physical < witness.begin || seed_end > witness.end {
            return None;
        }
        let location = physical
            .checked_mul(2)?
            .checked_add(u64::from(witness.reverse != direction))?;
        let record = if let Some(rank) = &self.leaf_rank {
            let id = rank.lookup(location)?;
            if !group.leaves.contains(&id) {
                return None;
            }
            let record = self.leaf_occurrences.get(id)?;
            if record.location != location {
                return None;
            }
            record
        } else {
            let records = self.leaf_occurrences.get(group.leaves.clone())?;
            &records[records
                .binary_search_by_key(&location, |r| r.location)
                .ok()?]
        };
        let node = record.leaves[usize::from(direction)];
        if node == NONE {
            return None;
        }
        let leaf = &self.pools[group.pool][node as usize];
        debug_assert_ne!(leaf.depth & LEAF, 0);
        let depth = leaf.depth & !LEAF;
        let proven = if witness.reverse {
            seed_end - witness.begin
        } else {
            witness.end - physical
        };
        if end.checked_sub(start)? < depth || proven < depth as u64 {
            return None;
        }
        Some([leaf.sample, leaf.second_location()])
    }

    fn lookup_with_witness(
        &self,
        reference: &Reference,
        bit: usize,
        direction: bool,
        query: &PackedSeqVec,
        start: usize,
        end: usize,
        k: usize,
        witnesses: &[Option<RecentAlignment>; 2],
        witness_hits: &mut usize,
    ) -> Option<[u64; 2]> {
        let value = *self.positions.get(bit)?;
        if value != NO_LOCATION && self.repeated[bit / 64] & (1 << (bit % 64)) != 0 {
            let group = &self.groups[value as usize];
            if !group.leaves.is_empty() {
                for witness in witnesses.iter().flatten() {
                    if let Some(locations) =
                        self.witnessed_leaf(group, direction, start, end, query.len(), k, *witness)
                    {
                        *witness_hits += 1;
                        return Some(locations);
                    }
                }
            }
        }
        self.lookup(reference, bit, direction, query, start, end, k)
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
        let mut recent = RecentAlignments::default();
        let mut witnesses = [None; 2];
        for run in &query.runs {
            stream.reset();
            let mut pos = run.start;
            let mut prefetched: Option<ReferenceSeed> = None;
            while pos + K <= run.end {
                let (id, offset, direction) = if let Some(seed) = prefetched.take() {
                    // A skipped SSHash call leaves its stream stale. Reset before
                    // any later fallback, including an adjacent seed.
                    stream.reset();
                    stats.reference_seed_hits += 1;
                    (seed.pattern as usize, seed.offset as usize, seed.reverse)
                } else {
                    let found = stream.lookup(pos, run.end);
                    if !found.is_found() {
                        let next = stream.next_after_miss(pos, &found);
                        stats.absent_windows_skipped += next - pos - 1;
                        pos = next;
                        continue;
                    }
                    (
                        found.string_id as usize,
                        found.kmer_id_in_string as usize,
                        found.kmer_orientation < 0,
                    )
                };
                let bit = index.position_offsets[id] + offset;
                if let Some(bits) = present.as_deref_mut() {
                    bits[bit / 64] |= 1 << (bit % 64);
                }
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
                let locations = self.lookup_with_witness(
                    reference,
                    bit,
                    direction,
                    &query.forward,
                    pos,
                    run.end,
                    K,
                    &witnesses,
                    &mut stats.context_witness_hits,
                );
                // A failed right-context refinement has no anchor to emit.
                // Complete occurrence coverage still permits a jump: every
                // later mapping already exists at this seed, even if other
                // copies disappear. The first differing DNA window and each
                // new mapping are visited; no arbitrary prefix is skipped.
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
                let mut seed_alignment: Option<(u64, bool, u64)> = None;
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
                        witnesses[usize::from(reverse)] = Some(recent.entries[slot]);
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
                                .max(self.known_alignment_span(physical, reverse, extension));
                            if seed_alignment.is_none_or(|(_, _, old)| extension > old) {
                                seed_alignment = Some((physical, reverse, extension));
                            }
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
                    let (left, right) =
                        alignment_context_uncached(reference, &anchor, query, hit, origin, stats)?;
                    recent.remember(
                        slot,
                        diagonal,
                        reverse,
                        physical - left,
                        physical + K as u64 + right,
                    );
                    witnesses[usize::from(reverse)] = Some(recent.entries[slot]);
                    if present.is_none() {
                        let extension = if reverse { left } else { right };
                        advance =
                            advance.max(self.known_alignment_span(physical, reverse, extension));
                        if seed_alignment.is_none_or(|(_, _, old)| extension > old) {
                            seed_alignment = Some((physical, reverse, extension));
                        }
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
                if pos + advance + K <= run.end {
                    if let (Some(map), Some((physical, reverse, extension))) =
                        (&self.reference_seeds, seed_alignment)
                    {
                        if advance as u64 <= extension {
                            let next = if reverse {
                                physical - advance as u64
                            } else {
                                physical + advance as u64
                            };
                            prefetched = map.lookup(index, next, reverse);
                        }
                    }
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
    #[test]
    fn reference_seed_directory_matches_independent_dictionary_on_both_strands() -> io::Result<()> {
        fn check<const K: usize>() -> io::Result<()>
        where
            Kmer<K>: KmerBits,
        {
            let temp = tempfile::tempdir()?;
            let mut rng = 555661133u64;
            let core: Vec<_> = (0..2100)
                .map(|_| {
                    rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
                    b"ACGT"[(rng >> 62) as usize]
                })
                .collect();
            let mut source = ProtoIndex {
                options: LocateOptions {
                    k: K,
                    min_length: K,

                    m: K.min(19),
                    sshash_m: None,
                    threads: 2,
                    chunk_bases: 127,
                    forward_only: false,
                    mumemto_format: false,
                    binary: false,
                },
                patterns: Vec::new(),
                sequences: PackedSeqVec::default(),
            };
            for range in [0..550, 600..1200, 1250..2100] {
                source.add_pattern(&core[range])?;
            }
            let index = FragmentIndex::from_graph(source, Some(temp.path()))?;
            let path = temp.path().join("directory.fa");
            let mut fasta = b">original\n".to_vec();
            fasta.extend(&core);
            fasta.extend(b"\n>reverse\n");
            fasta.extend(reverse(&core));
            fasta.extend(b"\n>ambiguous\n");
            fasta.extend(&core[..800]);
            fasta.extend(b"NNNNN");
            fasta.extend(&core[800..]);
            fasta.extend(b"\n>repeat\n");
            fasta.extend(&core[100..900]);
            fasta.push(b'\n');
            std::fs::write(&path, fasta)?;
            let mut reference = read_reference(&index, &path)?;
            let map = ReferenceSeeds::build(&index, &reference).expect("small directory");
            let mut checked = [0usize; 2];
            for (c, contig) in reference.contigs.iter().enumerate() {
                for pos in 0..contig.len() {
                    let physical = reference.meta.contigs[c].start + pos as u64;
                    for reverse_query in [false, true] {
                        let Some(seed) = map.lookup(&index, physical, reverse_query) else {
                            continue;
                        };
                        assert!(
                            contig.bounds(pos, K).is_ok(),
                            "must not cross an N or contig boundary"
                        );
                        let window: Vec<_> = (pos..pos + K)
                            .map(|p| contig.forward.as_slice().get_ascii(p))
                            .collect();
                        let window = if reverse_query {
                            reverse(&window)
                        } else {
                            window
                        };
                        let key = Kmer::<K>::from_ascii_unchecked(&window);
                        let found = index.compact.as_ref().unwrap().dictionary.query::<K>(&key);
                        assert!(found.is_found());
                        assert_eq!(
                            (seed.pattern as u64, seed.offset as u64, seed.reverse),
                            (
                                found.string_id,
                                found.kmer_id_in_string,
                                found.kmer_orientation < 0
                            )
                        );
                        checked[usize::from(reverse_query)] += 1;
                    }
                }
            }
            assert!(checked.iter().all(|&n| n > 5000));
            assert!(map.lookup(&index, u64::MAX, false).is_none());
            // Wide contig starts disable the dense array; sparse still indexes
            // the covered fragment coordinates.
            reference.meta.contigs.last_mut().unwrap().start += 1u64 << 40;
            let wide = ReferenceSeeds::build(&index, &reference).expect("sparse directory");
            assert!(matches!(wide.store, SeedStore::Sparse(_)));
            Ok(())
        }
        check::<31>()?;
        check::<33>()?;
        check::<63>()
    }

    #[test]
    fn reference_seed_directory_conflicts_are_permanent_and_empty_is_not_evidence() -> io::Result<()>
    {
        let temp = tempfile::tempdir()?;
        let mut source = ProtoIndex {
            options: LocateOptions {
                k: 7,
                min_length: 7,

                m: 5,
                sshash_m: None,
                threads: 1,
                chunk_bases: 127,
                forward_only: false,
                mumemto_format: false,
                binary: false,
            },
            patterns: Vec::new(),
            sequences: PackedSeqVec::default(),
        };
        // Canonical 7-mers here are distinct; synthetic fragment conflicts are
        // introduced only after the independent dictionary has been built.
        source.add_pattern(b"ACGATCC")?;
        source.add_pattern(b"GGATGCA")?;
        let index = FragmentIndex::from_graph(source, Some(temp.path()))?;
        let mut reference = reference(&[b"ACGATCCGGATGCA".to_vec()], 7);
        let fragment = |pattern, start, reverse| Fragment {
            begin: 0,
            end: 7,
            hit: Hit {
                pattern,
                contig: 0,
                start,
                reverse,
            },
        };
        reference.fragments = vec![
            vec![fragment(0, 0, false)],
            vec![fragment(1, 0, false), fragment(1, 0, false)],
        ];
        let map = ReferenceSeeds::build(&index, &reference).unwrap();
        assert_eq!(map.dense_entry(0), Some(ReferenceSeeds::UNKNOWN));
        assert!(map.lookup(&index, 0, false).is_none());
        assert_eq!(map.dense_entry(1), Some(ReferenceSeeds::EMPTY));
        assert!(map.lookup(&index, 1, false).is_none());
        // Repeated identical entries are usable; a different strand is ambiguous.
        reference.fragments = vec![vec![fragment(0, 0, false), fragment(0, 0, false)], vec![]];
        let map = ReferenceSeeds::build(&index, &reference).unwrap();
        assert_eq!(map.lookup(&index, 0, false).unwrap().pattern, 0);
        assert!(!map.lookup(&index, 0, false).unwrap().reverse);
        assert!(map.lookup(&index, 0, true).unwrap().reverse);
        reference.fragments[0].push(fragment(0, 0, true));
        reference.fragments[0].push(fragment(0, 0, false));
        let map = ReferenceSeeds::build(&index, &reference).unwrap();
        assert_eq!(map.dense_entry(0), Some(ReferenceSeeds::UNKNOWN));
        assert!(map.lookup(&index, 0, false).is_none());
        Ok(())
    }

    use super::*;

    #[test]
    fn multi_copy_stop_maps_match_full_occurrence_coverage_oracle() {
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
            let covered = |a: usize, b: usize| {
                let Some(&(ga, _)) = physical.get(&a) else {
                    return false;
                };
                let Some(&(gb, _)) = physical.get(&b) else {
                    return false;
                };
                ga != gb && {
                    let source = signature(a);
                    signature(b)
                        .into_iter()
                        .all(|mapping| source.binary_search(&mapping).is_ok())
                }
            };
            let mut merged = 0;
            for &p in physical.keys() {
                let expected_forward = p == 0 || !covered(p - 1, p);
                let expected_reverse = !covered(p + 1, p);
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
    fn reference_run_coverage_skips_copy_losses_but_visits_additions() {
        // Three A copies, two B copies, then three C copies. The second
        // chain is reverse-oriented: C- at 98, B- at 99, A- at 100.
        let copies = [
            (0, 20),
            (0, 201),
            (0, 400),
            (1, 22),
            (1, 199),
            (2, 24),
            (2, 197),
            (2, 600),
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
        let set = |map: &[u64], p: usize| map[p / 64] & (1 << (p % 64)) != 0;
        assert!(!set(&forward, 11)); // A -> B loses a copy.
        assert!(set(&reverse, 10)); // B -> A introduces a copy.
        assert!(set(&forward, 12)); // B -> C introduces a copy.
        assert!(!set(&reverse, 11)); // C -> B loses a copy.
        assert!(!set(&forward, 99));
        assert!(set(&reverse, 98));
        assert!(set(&forward, 100));
        assert!(!set(&reverse, 99));
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
                    min_length: k,
                    m: k.min(19),
                    sshash_m: None,
                    threads: 2,
                    chunk_bases: 127,
                    forward_only: false,
                    mumemto_format: false,
                    binary: false,
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
            let singleton_forward = restriction.context_singletons_forward.clone();
            let singleton_reverse = restriction.context_singletons_reverse.clone();
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
                let saved_leaves = restriction
                    .groups
                    .iter()
                    .map(|g| g.leaves.clone())
                    .collect::<Vec<_>>();
                let mut saved_seeds = None;
                for exhaustive in [false, true] {
                    if exhaustive {
                        saved_seeds = restriction.reference_seeds.take();
                    }
                    if exhaustive {
                        for group in &mut restriction.groups {
                            group.leaves = 0..0;
                        }
                    }
                    restriction.reference_repeats.clone_from(&forward);
                    restriction
                        .context_singletons_forward
                        .clone_from(&singleton_forward);
                    restriction
                        .context_singletons_reverse
                        .clone_from(&singleton_reverse);
                    restriction
                        .reference_repeats_reverse
                        .clone_from(&reverse_stops);
                    if exhaustive {
                        restriction.reference_repeats.fill(u64::MAX);
                        restriction.reference_repeats_reverse.fill(u64::MAX);
                        restriction.context_singletons_forward.fill(u64::MAX);
                        restriction.context_singletons_reverse.fill(u64::MAX);
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
                restriction.reference_seeds = saved_seeds;
                for (group, leaves) in restriction.groups.iter_mut().zip(saved_leaves) {
                    group.leaves = leaves;
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
                context_singletons_forward: Vec::new(),
                context_singletons_reverse: Vec::new(),
                reference_repeats: words,
                groups: Vec::new(),
                pools: Vec::new(),
                leaf_occurrences: Vec::new(),
                leaf_rank: None,
                reference_seeds: None,
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
    fn singleton_context_reserve_is_preserved_across_word_boundaries() {
        let count = 257usize;
        let ordinary = vec![u64::MAX; count.div_ceil(64)];
        let mut relaxed = ordinary.clone();
        for p in 0..count {
            if p % 31 != 0 {
                relaxed[p / 64] &= !(1u64 << (p % 64));
            }
        }
        let mut index = RestrictionIndex {
            positions: Vec::new(),
            repeated: Vec::new(),
            boundaries: Vec::new(),
            reference_repeats: ordinary.clone(),
            reference_repeats_reverse: ordinary,
            context_singletons_forward: relaxed.clone(),
            context_singletons_reverse: relaxed.clone(),
            groups: Vec::new(),
            pools: Vec::new(),
            leaf_occurrences: Vec::new(),
            leaf_rank: None,
            reference_seeds: None,
        };
        for position in 0..count {
            for reverse in [false, true] {
                let limit = if reverse {
                    position + 1
                } else {
                    count - position
                };
                for maximum in 1..=limit {
                    let extension = maximum - 1;
                    let original = 1;
                    let expected = if extension <= SINGLETON_EXTRA {
                        original
                    } else {
                        let safe = extension - SINGLETON_EXTRA + 1;
                        (1..safe)
                            .find(|&step| {
                                let p = if reverse {
                                    position - step
                                } else {
                                    position + step
                                };
                                relaxed[p / 64] & (1u64 << (p % 64)) != 0
                            })
                            .unwrap_or(safe)
                            .max(original)
                    };
                    assert_eq!(
                        index.known_alignment_span(position as u64, reverse, extension as u64),
                        expected,
                        "position={position} reverse={reverse} extension={extension}"
                    );
                    if expected > 1 {
                        assert!(expected - 1 + SINGLETON_EXTRA <= extension);
                    }
                }
            }
        }
        // With no remaining ordinary stops, the suffix fallback must cross
        // the reserve boundary and preserve the complete verified alignment.
        index.reference_repeats.fill(0);
        index.reference_repeats_reverse.fill(0);
        index.context_singletons_forward.fill(0);
        index.context_singletons_reverse.fill(0);
        for position in [0usize, 63, 64, 65, 127, 128, 256] {
            for reverse in [false, true] {
                let limit = if reverse {
                    position + 1
                } else {
                    count - position
                };
                for maximum in 1..=limit {
                    assert_eq!(
                        index.known_alignment_span(position as u64, reverse, (maximum - 1) as u64),
                        maximum
                    );
                }
            }
        }
    }

    #[test]
    fn singleton_context_stops_require_full_keys_bounded_depth_and_no_conflicting_mapping() {
        let k = 7;
        let locations = [20u64, 40, 60, 80, 101, 120, 140, 160];
        let mut nodes = Vec::new();
        let mut records = Vec::new();
        for (id, &location) in locations.iter().enumerate() {
            let first = nodes.len() as u32;
            let depth = if id == 1 {
                k + SINGLETON_EXTRA + 1
            } else {
                k + SINGLETON_EXTRA
            };
            let sample = if id == 2 { location + 2 } else { location };
            let second = if id == 3 { location + 200 } else { NO_LOCATION };
            nodes.push(ContextGroup::leaf(depth, sample, second));
            nodes.push(ContextGroup::leaf(k, location, NO_LOCATION));
            records.push(LeafOccurrence {
                location,
                leaves: if id == 6 {
                    [NONE, NONE]
                } else {
                    [first, first + 1]
                },
            });
        }
        // Position60 has a second repeated mapping. Position80 has a unique
        // mapping in addition to its repeated occurrence. Both remain blocked.
        let mut copies: Vec<_> = locations.iter().enumerate().map(|(i, &l)| (i, l)).collect();
        copies.push((100, 121));
        let mut positions = vec![0; locations.len() + 1];
        positions[locations.len()] = 160;
        let repeated = vec![(1u64 << locations.len()) - 1];
        let groups = vec![RepeatedPosition {
            first: 20,
            second: NO_LOCATION,
            roots: [NONE; 2],
            pool: 0,
            leaves: 0..records.len(),
        }];
        let original = vec![u64::MAX; 2];
        let (forward, reverse) = RestrictionIndex::build_context_singleton_stops(
            &positions,
            &repeated,
            &copies,
            &groups,
            &[nodes],
            &records,
            &original,
            &original,
            k,
        );
        let set =
            |words: &[u64], physical: usize| words[physical / 64] & (1u64 << (physical % 64)) != 0;
        for (id, &location) in locations.iter().enumerate() {
            for backwards in [false, true] {
                let direction = (location & 1 != 0) != backwards;
                let expected_blocked =
                    id == 5 || id == 6 || id == 7 || (!direction && [1, 2, 3].contains(&id));
                assert_eq!(
                    set(
                        if backwards { &reverse } else { &forward },
                        (location >> 1) as usize
                    ),
                    expected_blocked,
                    "id={id} backwards={backwards}"
                );
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
                    min_length: k,
                    m: k.min(19),
                    sshash_m: None,
                    threads: 2,
                    chunk_bases: 127,
                    forward_only: false,
                    mumemto_format: false,
                    binary: false,
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
            for variant in 0..2 {
                let mut copy = core[200..1400].to_vec();
                for p in (61 + variant * 17..copy.len()).step_by(149) {
                    copy[p] = if copy[p] == b'A' { b'C' } else { b'A' };
                }
                fasta.extend(format!("\n>context-copy-{variant}\n").bytes());
                fasta.extend(copy);
            }
            fasta.push(b'\n');
            std::fs::write(&path, fasta)?;
            let reference = read_reference(&index, &path)?;
            let mut restriction = RestrictionIndex::build(&index, &reference)?;
            assert!(
                restriction
                    .reference_repeats
                    .iter()
                    .zip(&restriction.context_singletons_forward)
                    .any(|(&old, &new)| old & !new != 0),
                "exercise context singleton skips k={k}"
            );
            if k >= 31 {
                let is_set = |words: &[u64], p: usize| words[p / 64] & (1 << (p % 64)) != 0;
                assert!(is_set(&restriction.reference_repeats, 200));
                assert!(!is_set(&restriction.reference_repeats, 201));
                assert!(is_set(&restriction.reference_repeats_reverse, 700 - k));
                assert!(!is_set(&restriction.reference_repeats_reverse, 699 - k));
            }
            let mut total_skipped = 0;
            let mut total_seed_hits = 0;
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
                let singleton_forward = restriction.context_singletons_forward.clone();
                let singleton_reverse = restriction.context_singletons_reverse.clone();
                let saved_leaves = restriction
                    .groups
                    .iter()
                    .map(|g| g.leaves.clone())
                    .collect::<Vec<_>>();
                let mut saved_seeds = None;
                for exhaustive in [false, true] {
                    if exhaustive {
                        saved_seeds = restriction.reference_seeds.take();
                    }
                    if exhaustive {
                        for group in &mut restriction.groups {
                            group.leaves = 0..0;
                        }
                    }
                    if exhaustive {
                        restriction.reference_repeats.fill(u64::MAX);
                        restriction.reference_repeats_reverse.fill(u64::MAX);
                        restriction.context_singletons_forward.fill(u64::MAX);
                        restriction.context_singletons_reverse.fill(u64::MAX);
                    }
                    let mut sink = ContextSink::new(3, &reference, k);
                    let mut stats = ScanStats::default();
                    restriction.scan(
                        &index, &reference, &query, 4, 91, &mut stats, &mut sink, None,
                    )?;
                    if !exhaustive {
                        total_skipped += stats.unique_alignment_windows_skipped;
                        total_seed_hits += stats.reference_seed_hits;
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
                restriction.reference_seeds = saved_seeds;
                for (group, leaves) in restriction.groups.iter_mut().zip(saved_leaves) {
                    group.leaves = leaves;
                }
                restriction.reference_repeats = repeats;
                restriction.reference_repeats_reverse = reverse_repeats;
                restriction.context_singletons_forward = singleton_forward;
                restriction.context_singletons_reverse = singleton_reverse;
                assert_eq!(rows[0], rows[1], "k={k}, query length={}", sequence.len());
            }
            assert!(total_seed_hits > 0, "exercise cached reference seeds k={k}");
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
                context_singletons_forward: Vec::new(),
                context_singletons_reverse: Vec::new(),
                groups: Vec::new(),
                pools: Vec::new(),
                leaf_occurrences: Vec::new(),
                leaf_rank: None,
                reference_seeds: None,
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
                context_singletons_forward: Vec::new(),
                context_singletons_reverse: Vec::new(),
                groups: vec![RepeatedPosition {
                    first: locations[0],
                    second: locations[1],
                    roots: [root; 2],
                    pool: 0,
                    leaves: 0..0,
                }],
                pools: vec![nodes],
                leaf_occurrences: Vec::new(),
                leaf_rank: None,
                reference_seeds: None,
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
            context_singletons_forward: Vec::new(),
            context_singletons_reverse: Vec::new(),
            groups: vec![RepeatedPosition {
                first: locations[0],
                second: locations[1],
                roots: [root; 2],
                pool: 0,
                leaves: 0..0,
            }],
            pools: vec![nodes],
            leaf_occurrences: Vec::new(),
            leaf_rank: None,
            reference_seeds: None,
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
    #[test]
    fn leaf_rank_preserves_full_keys_and_original_record_ids() {
        let locations = [8191u64, 0, 63, 64, 65, 127, 128, 4095, 4096, 3];
        let records: Vec<_> = locations
            .iter()
            .map(|&location| LeafOccurrence {
                location,
                leaves: [NONE; 2],
            })
            .collect();
        let rank = LeafRank::build(&records).unwrap();
        for location in 0..9000 {
            assert_eq!(
                rank.lookup(location),
                locations.iter().position(|&key| key == location),
                "location {location}"
            );
        }
        assert_eq!(rank.lookup(u64::MAX), None);
        let mut state = 0x8a33_2198_c49a_bde1u64;
        let mut locations = std::collections::BTreeSet::new();
        for _ in 0..10000 {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            locations.insert(state % 100000);
        }
        let records: Vec<_> = locations
            .iter()
            .rev()
            .map(|&location| LeafOccurrence {
                location,
                leaves: [NONE; 2],
            })
            .collect();
        let rank = LeafRank::build(&records).unwrap();
        for location in 0..100100 {
            let actual = rank.lookup(location);
            assert_eq!(actual.is_some(), locations.contains(&location));
            if let Some(id) = actual {
                assert_eq!(records[id].location, location);
            }
        }
    }

    #[test]
    fn leaf_rank_falls_back_for_conflicts_and_wide_sparse_coordinates() {
        assert!(LeafRank::build(&[]).is_none());
        for locations in [vec![1u64, 1], vec![0, 1 << 44], vec![u64::MAX]] {
            let records: Vec<_> = locations
                .iter()
                .map(|&location| LeafOccurrence {
                    location,
                    leaves: [NONE; 2],
                })
                .collect();
            assert!(LeafRank::build(&records).is_none());
        }
    }

    #[test]
    fn alignment_witness_resolves_only_the_original_exact_leaf() -> io::Result<()> {
        let mut state = 0x69ac_19d2_7f35_138bu64;
        let core: Vec<_> = (0..220)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[(state & 3) as usize]
            })
            .collect();
        let mut changed = core.clone();
        changed[170] = if changed[170] == b'A' { b'C' } else { b'A' };
        let mut ambiguous = core.clone();
        ambiguous[140..145].fill(b'N');
        let sequences = vec![
            core.clone(),
            changed,
            core[..180].to_vec(),
            reverse(&core),
            core[20..190].to_vec(),
            ambiguous,
            [b"ACG".repeat(45), b"TTGAC".to_vec()].concat(),
            [b"ACG".repeat(42), b"CCATT".to_vec()].concat(),
            b"ACG".repeat(43),
        ];
        for k in [3usize, 7, 31, 33, 63] {
            for origin in [0u64, 1 << 44] {
                let mut reference = reference(&sequences, k);
                for contig in &mut reference.meta.contigs {
                    contig.start += origin;
                }
                // Build occurrence groups independently from the graph scanner.
                let mut canonical = std::collections::BTreeMap::<Vec<u8>, Vec<u64>>::new();
                for (id, seq) in sequences.iter().enumerate() {
                    for run in &reference.contigs[id].runs {
                        for pos in run.start..=run.end - k {
                            let window = &seq[pos..pos + k];
                            let rc = reverse(window);
                            let strand = rc.as_slice() < window;
                            canonical
                                .entry(if strand { rc } else { window.to_vec() })
                                .or_default()
                                .push(
                                    ((reference.meta.contigs[id].start + pos as u64) * 2)
                                        | u64::from(strand),
                                );
                        }
                    }
                }
                let mut groups = Vec::new();
                let mut nodes = Vec::new();
                let mut records = Vec::new();
                for locations in canonical.into_values().filter(|v| v.len() > 2) {
                    let begin = nodes.len();
                    let forward = refine(&reference, locations.clone(), false, k, k, &mut nodes)?;
                    let middle = nodes.len();
                    let reverse = refine(&reference, locations.clone(), true, k, k, &mut nodes)?;
                    let leaf_begin = records.len();
                    records.extend(collect_leaf_occurrences(
                        &nodes,
                        [begin..middle, middle..nodes.len()],
                    ));
                    groups.push(RepeatedPosition {
                        first: locations[0],
                        second: locations[1],
                        roots: [forward, reverse],
                        pool: 0,
                        leaves: leaf_begin..records.len(),
                    });
                }
                let index = RestrictionIndex {
                    positions: (0..groups.len() as u64).collect(),
                    repeated: vec![u64::MAX; groups.len().div_ceil(64)],
                    boundaries: Vec::new(),
                    reference_repeats: Vec::new(),
                    reference_repeats_reverse: Vec::new(),
                    context_singletons_forward: Vec::new(),
                    context_singletons_reverse: Vec::new(),
                    groups,
                    pools: vec![nodes],
                    leaf_rank: LeafRank::build(&records),
                    reference_seeds: None,
                    leaf_occurrences: records,
                };
                assert_eq!(index.leaf_rank.is_some(), origin == 0);
                let mut exercised = [0usize; 2];
                let mut shortened = 0;
                for (bit, group) in index.groups.iter().enumerate() {
                    for record in &index.leaf_occurrences[group.leaves.clone()] {
                        for direction in [false, true] {
                            let node = record.leaves[usize::from(direction)];
                            if node == NONE {
                                continue;
                            }
                            let depth = index.pools[0][node as usize].depth & !LEAF;
                            let v = view(&reference, record.location, direction, k);
                            let local =
                                (v.physical - reference.meta.contigs[v.contig].start) as usize;
                            let bounds = reference.contigs[v.contig].bounds(local, k).unwrap();
                            let oriented = if v.reverse {
                                v.seq.len() - v.start - k
                            } else {
                                v.start
                            };
                            let witness = RecentAlignment {
                                diagonal: oriented as i128 - v.physical as i128,
                                begin: reference.meta.contigs[v.contig].start + bounds.start as u64,
                                end: reference.meta.contigs[v.contig].start + bounds.end as u64,
                                reverse: v.reverse,
                            };
                            for remaining in [
                                k,
                                k + 1,
                                depth.saturating_sub(1),
                                depth,
                                depth + 1,
                                v.end - v.start,
                            ] {
                                if remaining < k || remaining > v.end - v.start {
                                    continue;
                                }
                                let actual = index.witnessed_leaf(
                                    group,
                                    direction,
                                    v.start,
                                    v.start + remaining,
                                    v.seq.len(),
                                    k,
                                    witness,
                                );
                                let expected = index.lookup(
                                    &reference,
                                    bit,
                                    direction,
                                    v.seq,
                                    v.start,
                                    v.start + remaining,
                                    k,
                                );
                                if remaining >= depth {
                                    assert_eq!(
                                        actual, expected,
                                        "k={k}, origin={origin}, direction={direction}"
                                    );
                                    assert!(actual.is_some());
                                    exercised[usize::from(direction)] += 1;
                                } else {
                                    assert!(actual.is_none());
                                    shortened += 1;
                                }
                            }
                            if depth > k {
                                let mut short = witness;
                                if short.reverse {
                                    short.begin = v.physical;
                                } else {
                                    short.end = v.physical + k as u64;
                                }
                                assert!(
                                    index
                                        .witnessed_leaf(
                                            group,
                                            direction,
                                            v.start,
                                            v.end,
                                            v.seq.len(),
                                            k,
                                            short
                                        )
                                        .is_none()
                                );
                            }
                            for diagonal in [i128::MAX, i128::MIN] {
                                assert!(
                                    index
                                        .witnessed_leaf(
                                            group,
                                            direction,
                                            v.start,
                                            v.end,
                                            v.seq.len(),
                                            k,
                                            RecentAlignment {
                                                diagonal,
                                                ..witness
                                            }
                                        )
                                        .is_none()
                                );
                            }
                            assert!(
                                index
                                    .witnessed_leaf(
                                        group,
                                        direction,
                                        v.start,
                                        v.end,
                                        v.seq.len(),
                                        k,
                                        RecentAlignment {
                                            begin: witness.end,
                                            ..witness
                                        }
                                    )
                                    .is_none()
                            );
                        }
                    }
                }
                assert!(exercised.into_iter().all(|n| n > 100));
                assert!(shortened > 100);
            }
        }
        Ok(())
    }

    #[test]
    fn alignment_witness_leaf_metadata_rejects_conflicting_assignments() {
        let nodes = vec![
            ContextGroup::leaf(7, 4, NO_LOCATION),
            ContextGroup::leaf(8, 4, 1 << 45),
            ContextGroup::leaf(9, 4, 1 << 45),
        ];
        let records = collect_leaf_occurrences(&nodes, [0..2, 2..3]);
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].location, 4);
        assert_eq!(records[0].leaves, [NONE, 2]);
        assert_eq!(records[1].location, 1 << 45);
        assert_eq!(records[1].leaves, [1, 2]);
    }
}
