//! GBWT construction.
//!
//! Like the C++ implementation, the construction algorithm uses 32-bit integers internally to save space.
//! This limits the maximum node identifier and the number of visits to a node to approximately [`u32::MAX`].
// FIXME: document

use crate::ENDMARKER;
use crate::{Pos, FullPathName};
use crate::bwt::SmallPos;
use crate::support::SmallRun;

use std::collections::BTreeMap;
use std::iter::FusedIterator;

//-----------------------------------------------------------------------------

// FIXME: builder itself
// Construction has three parameters: bidirectional, metadata, buffer size
// There is a background thread that contains a DynamicGBWT and a MetadataBuilder
// When the buffer gets full, the main thread sends the paths (and the metadata) to the background thread
// When the construction is finished, the main thread sends the remaining buffer, followed by an empty buffer to signal the end of construction
// The background thread inserts each batch into the DynamicGBWT and the MetadataBuilder
// When it receives the end signal, it converts the DynamicGBWT and MetadataBuilder into GBWT and sends it back
// Communication between threads uses channels

//-----------------------------------------------------------------------------

/// A sorted list of edges stored as [`SmallPos`].
///
/// If there are at most [`EdgeList::SMALL_CAPACITY`] edges, they are stored inline in a sorted array.
/// Larger edge sets are stored as a map from node to offset.
///
/// # Examples
///
/// ```
/// use gbz::gbwt::builder::EdgeList;
/// use gbz::Pos;
///
/// let mut edges = EdgeList::new();
/// edges.increment(3, 10);
/// edges.increment(1, 5);
/// edges.increment(2, 7);
/// edges.increment(1, 6);
/// edges.increment(4, 12); // Triggers conversion to large representation.
///
/// assert_eq!(edges.len(), 4);
/// assert_eq!(edges.get(1), Some(11));
/// assert_eq!(edges.get(5), None);
///
/// let truth = vec![
///     Pos::new(1, 11),
///     Pos::new(2, 7),
///     Pos::new(3, 10),
///     Pos::new(4, 12)
/// ];
/// let edges: Vec<Pos> = edges.iter().collect();
/// assert_eq!(edges, truth);
/// ```
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum EdgeList {
    // Up to three edges as (len, edges), with the edges in sorted order.
    Small(u8, [SmallPos; Self::SMALL_CAPACITY]),
    // More than three edges as a map from node to offset.
    Large(BTreeMap<u32, u32>),
}

impl EdgeList {
    /// Capacity of the small representation.
    pub const SMALL_CAPACITY: usize = 3;

    /// Creates an empty edge list.
    pub fn new() -> Self {
        EdgeList::Small(0, [SmallPos::default(); Self::SMALL_CAPACITY])
    }

    /// Returns the number of edges in the list.
    #[inline]
    pub fn len(&self) -> usize {
        match self {
            EdgeList::Small(len, _) => *len as usize,
            EdgeList::Large(map) => map.len(),
        }
    }

    /// Returns `true` if the list is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Increments the offset in the given edge by the given amount.
    ///
    /// Returns the offset before the increment, or `0` if the edge did not exist.
    /// Inserts a new edge if it does not already exist.
    pub fn increment(&mut self, node: usize, amount: usize) -> u32 {
        match self {
            EdgeList::Small(len, edges) => {
                for edge in edges.iter_mut().take(*len as usize) {
                    if edge.node == node as u32 {
                        let val = edge.offset;
                        edge.offset += amount as u32;
                        return val;
                    }
                }
                if *len < Self::SMALL_CAPACITY as u8 {
                    edges[*len as usize] = SmallPos::new(node, amount);
                    // Bubble up the new edge to maintain sorted order.
                    for i in (1..=*len as usize).rev() {
                        if edges[i] < edges[i - 1] {
                            edges.swap(i, i - 1);
                        } else {
                            break;
                        }
                    }
                    *len += 1;
                } else {
                    // Convert to large representation.
                    let mut map = BTreeMap::new();
                    for edge in edges.iter().take(*len as usize) {
                        map.insert(edge.node, edge.offset);
                    }
                    map.insert(node as u32, amount as u32);
                    *self = EdgeList::Large(map);
                }
                return 0;
            }
            EdgeList::Large(map) => {
                let entry = map.entry(node as u32).or_insert(0);
                let val = *entry;
                *entry += amount as u32;
                return val;
            }
        }
    }

    /// Returns the offset in the given edge, or [`None`] if it does not exist.
    pub fn get(&self, node: usize) -> Option<u32> {
        match self {
            EdgeList::Small(len, edges) => {
                for edge in edges.iter().take(*len as usize) {
                    if edge.node == node as u32 {
                        return Some(edge.offset);
                    }
                }
                None
            }
            EdgeList::Large(map) => map.get(&(node as u32)).copied(),
        }
    }

    /// Returns a mutable reference to the offset in the given edge, or [`None`] if it does not exist.
    pub fn get_mut(&mut self, node: usize) -> Option<&mut u32> {
        match self {
            EdgeList::Small(len, edges) => {
                for edge in edges.iter_mut().take(*len as usize) {
                    if edge.node == node as u32 {
                        return Some(&mut edge.offset);
                    }
                }
                None
            }
            EdgeList::Large(map) => map.get_mut(&(node as u32)),
        }
    }

    /// Returns an iterator over the edges in the list, in sorted order.
    pub fn iter(&self) -> EdgeListIter<'_> {
        match self {
            EdgeList::Small(_, _) => EdgeListIter {
                parent: self,
                index: 0,
                iter: std::collections::btree_map::Iter::<'_, u32, u32>::default(),
            },
            EdgeList::Large(map) => EdgeListIter {
                parent: self,
                index: 0,
                iter: map.iter(),
            },
        }
    }
}

impl Default for EdgeList {
    fn default() -> Self {
        Self::new()
    }
}

/// An iterator over the edges in an [`EdgeList`].
///
/// The value of `Item` is [`Pos`].
pub struct EdgeListIter<'a> {
    parent: &'a EdgeList,
    index: usize,
    iter: std::collections::btree_map::Iter<'a, u32, u32>,
}

impl<'a> Iterator for EdgeListIter<'a> {
    type Item = Pos;

    fn next(&mut self) -> Option<Self::Item> {
        match self.parent {
            EdgeList::Small(len, edges) => {
                if self.index < *len as usize {
                    let pos = edges[self.index];
                    self.index += 1;
                    Some(Pos::from(pos))
                } else {
                    None
                }
            }
            EdgeList::Large(map) => {
                if self.index < map.len() {
                    self.index += 1;
                }
                self.iter.next().map(|(node, offset)| Pos::new(*node as usize, *offset as usize))
            }
        }
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        let len = self.parent.len() - self.index;
        (len, Some(len))
    }
}

impl ExactSizeIterator for EdgeListIter<'_> {}

impl FusedIterator for EdgeListIter<'_> {}

//-----------------------------------------------------------------------------

// FIXME: examples, tests
/// A mutable node record corresponding to [`crate::bwt::Record`] used for GBWT construction.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MutableRecord {
    // Number of sequence visits to the node (BWT fragment length).
    len: usize,
    // Incoming edges as (source node id, number of visits from that node).
    incoming: EdgeList,
    // Outgoing edges as (target node id, BWT offset of the first visit to that node).
    outgoing: EdgeList,
    // Run-length encoded BWT fragment as a vector of (node id, run length) pairs.
    bwt: Vec<SmallRun>,
}

impl MutableRecord {
    /// Creates an empty record.
    pub fn new() -> Self {
        Self::default()
    }

    /// Returns the number of visits to the node (BWT fragment length).
    pub fn len(&self) -> usize {
        self.len
    }

    /// Returns `true` if the record is empty (this should not happen).
    pub fn empty(&self) -> bool {
        self.len == 0
    }

    /// Sets the number of sequence visits from the given node.
    pub fn set_visits(&mut self, from: usize, count: usize) {
        if let Some(edge_count) = self.incoming.get_mut(from) {
            *edge_count = count as u32;
        } else {
            self.incoming.increment(from, count);
        }
    }

    /// Increases the number of sequence visits from the given node.
    pub fn add_visits(&mut self, from: usize, count: usize) {
        self.incoming.increment(from, count);
    }

    /// Adds an outgoing edge to the given node with an unknown offset.
    ///
    /// Does nothing if the edge already exists.
    pub fn add_edge(&mut self, to: usize) {
        self.outgoing.increment(to, 0);
    }

    /// Sets the offset of the outgoing edge to the given node.
    ///
    /// Does nothing if the edge does not exist.
    pub fn set_edge_offset(&mut self, to: usize, offset: usize) {
        if let Some(edge_offset) = self.outgoing.get_mut(to) {
            *edge_offset = offset as u32;
        }
    }
}

//-----------------------------------------------------------------------------

// A sequence iterator used for GBWT construction.
// Sort order is `(next.node, pos)`, which is consistent with LF-mapping with `next.node`.
#[derive(Clone, Debug, PartialEq, Eq)]
struct Sequence<'a> {
    // Sequence as a slice of node identifiers, excluding the endmarker.
    // We do not store the sequence offset, as it is shared between all sequences.
    nodes: &'a [u32],
    // GBWT position for the current sequence offset.
    curr: SmallPos,
    // GBWT position for the next sequence offset.
    // Offset in the position may be undetermined.
    next: SmallPos,
}

impl<'a> PartialOrd for Sequence<'a> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some((self.next.node, self.curr).cmp(&(other.next.node, other.curr)))
    }
}

impl<'a> Ord for Sequence<'a> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        (self.next.node, self.curr).cmp(&(other.next.node, other.curr))
    }
}

// A support structure for building a sequence of maximal runs incrementally.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
struct RunBuilder {
    source_len: usize,
    target_len: usize,
    runs: Vec<SmallRun>,
    ranks: EdgeList,
}

impl RunBuilder {
    fn new() -> Self {
        Self::default()
    }

    // Inserts a source run until the given offset, and returns the remaining run.
    fn insert_source_run(&mut self, source: SmallRun, offset: usize) -> SmallRun {
        let (to_insert, to_return) = if offset >= self.source_len + source.len as usize {
            (source, SmallRun::default())
        } else {
            let to_insert = SmallRun::new(source.value as usize, offset - self.source_len);
            let to_return = SmallRun::new(source.value as usize, self.source_len + source.len as usize - offset);
            (to_insert, to_return)
        };
        if let Some(last) = self.runs.last_mut() && last.value == to_insert.value {
            last.len += to_insert.len;
        } else {
            self.runs.push(to_insert);
        }
        self.source_len += to_insert.len as usize;
        self.target_len += to_insert.len as usize;
        self.ranks.increment(to_insert.value as usize, to_insert.len as usize);
        to_return
    }

    // Inserts a target value.
    fn insert_target_value(&mut self, value: usize) {
        if let Some(last) = self.runs.last_mut() && last.value == value as u32 {
            last.len += 1;
        } else {
            self.runs.push(SmallRun::new(value, 1));
        }
        self.target_len += 1;
        self.ranks.increment(value, 1);
    }
}

//-----------------------------------------------------------------------------

// FIXME: implement, document, tests, examples
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MutableGBWT {
    // Total length of the sequences, including endmarkers.
    len: usize,
    // Is the GBWT bidirectional?
    bidirectional: bool,
    // Endmarker record as an array of starting positions for each sequence.
    endmarker: Vec<SmallPos>,
    // Outgoing edges from the endmarker as (node id, number of sequences starting with that node).
    endmarker_edges: EdgeList,
    // Other records as a map from node id to mutable record.
    records: BTreeMap<usize, MutableRecord>,
    // Optional metadata as a vector of path names.
    metadata: Option<Vec<FullPathName>>,
}

/// Queries.
impl MutableGBWT {
    /// Returns the total length of the sequences, including endmarkers.
    pub fn len(&self) -> usize {
        self.len
    }

    /// Returns `true` if the GBWT is empty (this should not happen).
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Returns `true` if the GBWT is bidirectional.
    pub fn is_bidirectional(&self) -> bool {
        self.bidirectional
    }

    /// Returns `true` if the GBWT has path metadata.
    pub fn has_metadata(&self) -> bool {
        self.metadata.is_some()
    }

    /// Returns the total number of nodes, excluding the endmarker.
    pub fn node_count(&self) -> usize {
        self.records.len()
    }

    /// Returns the minimum node identifier, or [`None`] if there are no nodes.
    pub fn min_node(&self) -> Option<usize> {
        self.records.keys().next().copied()
    }

    /// Returns the maximum node identifier, or [`None`] if there are no nodes.
    pub fn max_node(&self) -> Option<usize> {
        self.records.keys().next_back().copied()
    }

    /// Returns the total number of sequences.
    pub fn sequence_count(&self) -> usize {
        self.endmarker.len()
    }
}

/// Construction.
impl MutableGBWT {
    /// Creates an empty GBWT.
    pub fn new(bidirectional: bool, with_metadata: bool) -> Self {
        Self {
            len: 0,
            bidirectional,
            endmarker: Vec::new(),
            endmarker_edges: EdgeList::new(),
            records: BTreeMap::new(),
            metadata: if with_metadata { Some(Vec::new()) } else { None },
        }
    }

    // Returns the sequences in the buffer without determining the GBWT positions.
    fn sequences_in_buffer<'a>(&'_ self, buffer: &'a [u32]) -> Result<Vec<Sequence<'a>>, String> {
        let mut sequences: Vec<Sequence<'a>> = Vec::new();
        let mut start = 0;
        for (i, &node) in buffer.iter().enumerate() {
            if node == ENDMARKER as u32 {
                sequences.push(Sequence {
                    nodes: &buffer[start..i],
                    curr: SmallPos::default(),
                    next: SmallPos::default(),
                });
                start = i + 1;
            }
        }
        if start < buffer.len() {
            return Err(String::from("MutableGBWT: Trailing nodes in the buffer"));
        }
        if self.bidirectional && sequences.len() % 2 != 0 {
            return Err(String::from("MutableGBWT: Odd number of sequences in the buffer for a bidirectional GBWT"));
        }
        Ok(sequences)
    }

    // Appends the given path names to the metadata.
    fn append_metadata(&mut self, seq_count: usize, names: Option<&[FullPathName]>) -> Result<(), String> {
        match (&mut self.metadata, names) {
            (Some(metadata), Some(names)) => {
                let expected = if self.bidirectional { seq_count / 2 } else { seq_count };
                if names.len() != expected {
                    return Err(format!("MutableGBWT: Expected {} path names, got {}", expected, names.len()));
                }
                metadata.extend_from_slice(names);
                Ok(())
            }
            (None, Some(_)) => Err(String::from("MutableGBWT: Path names provided for a GBWT without metadata")),
            (Some(_), None) => Err(String::from("MutableGBWT: No path names provided for a GBWT with metadata")),
            (None, None) => Ok(()),
        }
    }

    // Appends sequence starts to the endmarker record and updates the incoming edges of the first nodes.
    fn append_sequence_starts(&mut self, sequences: &mut [Sequence]) {
        for sequence in sequences.iter_mut() {
            let start = sequence.nodes.first().copied().unwrap_or(ENDMARKER as u32) as usize;
            let offset = self.endmarker_edges.increment(start as usize, 1) as usize;
            let pos = SmallPos::new(start, offset);
            sequence.curr = pos;
            self.endmarker.push(pos);
        }
        for edge in self.endmarker_edges.iter().filter(|edge| edge.node != ENDMARKER) {
            self.records.entry(edge.node).or_default().set_visits(ENDMARKER, edge.offset);
        }
    }

    // Inserts `next.node` to record offset `curr.offset` for node `curr.node` in the GBWT.
    // Sets `next.offset` to the rank of `next.node` at `curr.offset`.
    // Ensures that the outgoing edge and the successor record exist.
    // Increments the number of visits from `curr.node` in the successor record.
    fn update_records(&mut self, sequences: &mut [Sequence]) {
        // Update the BWT fragments.
        let mut i = 0;
        let mut edges: Vec<(u32, u32)> = Vec::with_capacity(sequences.len());
        while i < sequences.len() {
            let curr = sequences[i].curr.node;
            let current = self.records.get_mut(&(curr as usize)).unwrap();
            let mut run_iter = current.bwt.iter().copied();
            let mut remaining = SmallRun::default();
            let mut new_bwt = RunBuilder::new();
            while i < sequences.len() && sequences[i].curr.node == curr {
                let seq = &mut sequences[i];
                edges.push((seq.curr.node, seq.next.node));
                while new_bwt.target_len < seq.curr.offset as usize {
                    if remaining.len == 0 {
                        remaining = run_iter.next().unwrap(); // FIXME: error handling?
                    }
                    remaining = new_bwt.insert_source_run(remaining, seq.curr.offset as usize);
                }
                // We need to get the rank before inserting the new value.
                seq.next.offset = new_bwt.ranks.get(seq.next.node as usize).unwrap_or(0);
                new_bwt.insert_target_value(seq.next.node as usize);
                i += 1;
            }
            if remaining.len > 0 {
                new_bwt.insert_source_run(remaining, usize::MAX);
            }
            while let Some(run) = run_iter.next() {
                new_bwt.insert_source_run(run, usize::MAX);
            }
            current.bwt = new_bwt.runs;
        }

        // Ensure that the outgoing edges exist.
        edges.sort_unstable();
        for (i, &(from, to)) in edges.iter().enumerate() {
            if i > 0 && edges[i - 1].0 == from {
                continue;
            }
            let predecessor = self.records.get_mut(&(from as usize)).unwrap();
            predecessor.add_edge(to as usize);
        }

        // Ensure that the successor records exist and increment the visits.
        edges.sort_unstable_by_key(|&(_, to)| to);
        let mut start = 0;
        for (i, &(from, to)) in edges.iter().enumerate() {
            if i + 1 >= edges.len() || edges[i + 1].1 != to {
                if to != ENDMARKER as u32 {
                    let successor = self.records.entry(to as usize).or_default();
                    successor.add_visits(from as usize, i - start + 1);
                }
                start = i + 1;
            }
        }
    }

    // Sorts the active sequences by `(next.node, curr)` and returns a new slice of active sequences.
    // Active sequences are those that have not reached the endmarker at `next`.
    fn sort_sequences<'a>(&mut self, sequences: &'a mut [Sequence<'a>]) -> &'a mut [Sequence<'a>] {
        sequences.sort_unstable();
        let mut i = 0;
        while i < sequences.len() && sequences[i].next.node == ENDMARKER as u32 {
            i += 1;
        }
        &mut sequences[i..]
    }

    // Updates the offsets of the outgoing edges from the current nodes of the active sequences.
    // Then determines the record offsets in the GBWT positions for the next nodes.
    fn determine_offsets(&mut self, sequences: &mut [Sequence]) {
        let mut next = 0; // We never use the offsets in outgoing edges to the endmarker.
        for sequence in sequences.iter() {
            if sequence.next.node == next {
                continue;
            }
            next = sequence.next.node;
            let mut offset = 0;
            let successor = self.records.get(&(next as usize)).unwrap();
            let incoming: Vec<Pos> = successor.incoming.iter().collect();
            for edge in incoming {
                if edge.node == ENDMARKER {
                    // The offset is always 0 in an outgoing edge from the endmarker.
                    offset += self.endmarker_edges.get(edge.node as usize).unwrap_or(0) as usize;
                } else {
                    let predecessor = self.records.get_mut(&(edge.node)).unwrap();
                    predecessor.set_edge_offset(next as usize, offset);
                    offset += edge.offset;
                }
            }
        }

        for sequence in sequences.iter_mut() {
            let current = self.records.get(&(sequence.curr.node as usize)).unwrap();
            let start_offset = current.outgoing.get(sequence.next.node as usize).unwrap_or(0);
            // The offset was the rank of `next.node` at offset `curr.offset`.
            sequence.next.offset += start_offset;
        }
    }

    // Advances the sequences to the next offset.
    // The sequences are now sorted by `curr`.
    // Determines the node but not the offset in the new successor GBWT position.
    fn advance_sequences(&mut self, sequences: &mut [Sequence], next_offset: usize) {
        for sequence in sequences.iter_mut() {
            let next_node = sequence.nodes.get(next_offset).copied().unwrap_or(ENDMARKER as u32) as usize;
            sequence.curr = sequence.next;
            sequence.next = SmallPos::new(next_node, 0);
        }
    }

    /// Inserts a batch of sequences into the GBWT.
    ///
    /// The sequences are concatenated in the buffer, with each terminated by [`ENDMARKER`].
    /// If the GBWT is bidirectional, each sequence is assumed to be followed by its reverse complement.
    ///
    /// If the GBWT contains metadata, path names must be provided.
    /// Each name corresponds to a sequence (in unidirectional GBWT) or a pair of sequences (in bidirectional GBWT).
    ///
    /// # Errors
    ///
    /// Returns an error, if:
    ///
    /// * There are trailing nodes after the last endmarker.
    /// * The number of sequences is odd in a bidirectional GBWT.
    /// * The number of path names is not as expected.
    ///
    /// Does not modify the GBWT if an error is returned.
    pub fn insert(&mut self, buffer: &[u32], names: Option<&[FullPathName]>) -> Result<(), String> {
        let mut sequences = self.sequences_in_buffer(buffer)?;
        self.append_metadata(sequences.len(), names)?;
        self.append_sequence_starts(&mut sequences);
        let mut active_sequences = self.sort_sequences(&mut sequences);
        self.len += buffer.len();

        // Invariants:
        // * GBWT position `curr` corresponds to `nodes[offset - 1]`.
        // * The sequences are sorted by `curr`.
        // * `nodes[0..offset]` have been inserted.
        // * `next.node` (`nodes[offset]`) should be inserted to GBWT offset `curr.offset` in the record for `curr.node`.
        // * If `nodes[offset]` exists, incoming edges for `next.node` already include the visit from `curr.node`.
        let mut offset = 0;
        while !active_sequences.is_empty() {
            self.update_records(active_sequences);
            offset += 1;
            active_sequences = self.sort_sequences(active_sequences);
            self.determine_offsets(active_sequences);
            self.advance_sequences(active_sequences, offset);
        }

        Ok(())
    }
}

//-----------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[ignore]
    fn assumptions() {
        assert_eq!(std::mem::size_of::<EdgeList>(), 32, "EdgeList should be 32 bytes");
        assert_eq!(std::mem::size_of::<MutableRecord>(), 96, "MutableRecord should be 96 bytes");
        assert_eq!(std::mem::size_of::<Sequence>(), 32, "Sequence should be 32 bytes");
    }

    #[test]
    fn edge_list_empty() {
        let edges = EdgeList::new();
        assert!(edges.is_empty(), "Edge list should be empty");
        assert_eq!(edges.len(), 0, "Edge list should have length 0");
        assert_eq!(edges.get(1), None, "Edge list should not contain any edges");
        assert_eq!(edges.iter().count(), 0, "Edge list iterator should be empty");
    }

    fn create_edge_list_and_sort(truth: &mut [Pos]) -> EdgeList {
        let mut edges = EdgeList::new();
        for edge in truth.iter() {
            edges.increment(edge.node, edge.offset);
        }
        truth.sort();
        edges
    }

    // Assumes that the truth is sorted and non-empty.
    fn check_edge_list(edges: &EdgeList, truth: &[Pos]) {
        assert!(!edges.is_empty(), "Edge list should not be empty");
        assert_eq!(edges.len(), truth.len(), "Edge list should have length {}", truth.len());

        let mut max_node = 0;
        for edge in truth.iter() {
            let offset = edge.offset as u32;
            assert_eq!(edges.get(edge.node), Some(offset), "Edge list should contain edge ({}, {})", edge.node, edge.offset);
            if edge.node > max_node {
                max_node = edge.node;
            }
        }
        assert!(edges.get(max_node + 1).is_none(), "Edge list should not contain edge ({}, _)", max_node + 1);

        // Increment the offsets in a copy and check them again.
        let mut copy = edges.clone();
        for edge in truth.iter() {
            let offset = copy.get_mut(edge.node);
            assert!(offset.is_some(), "Edge list should contain edge ({}, _)", edge.node);
            let offset = offset.unwrap();
            *offset += 1;
        }
        for edge in truth.iter() {
            let offset = edge.offset as u32 + 1;
            assert_eq!(copy.get(edge.node), Some(offset), "Incremented edge list should contain edge ({}, {})", edge.node, offset);
        }
    }

    #[test]
    fn edge_list_small() {
        let mut truth = vec![
            Pos::new(3, 10),
            Pos::new(1, 5),
            Pos::new(2, 7),
        ];
        let edges = create_edge_list_and_sort(&mut truth);
        check_edge_list(&edges, &truth);
    }

    #[test]
    fn edge_list_large() {
        let mut truth = vec![
            Pos::new(3, 10),
            Pos::new(1, 5),
            Pos::new(2, 7),
            Pos::new(4, 12),
            Pos::new(10, 18),
        ];
        let edges = create_edge_list_and_sort(&mut truth);
        check_edge_list(&edges, &truth);
    }

    #[test]
    fn edge_list_increment() {
        let mut truth = vec![
            Pos::new(3, 10),
            Pos::new(1, 6),
            Pos::new(2, 8),
            Pos::new(4, 12),
            Pos::new(10, 18),
        ];
        let mut edges = EdgeList::new();
        for edge in truth.iter() {
            edges.increment(edge.node, edge.offset / 2);
            edges.increment(edge.node, edge.offset / 2);
        }
        truth.sort();
        check_edge_list(&edges, &truth);
    }
}

//-----------------------------------------------------------------------------
