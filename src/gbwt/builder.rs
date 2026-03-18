//! GBWT construction.
// FIXME: document; uses 32-bit integers internally

use crate::Pos;
use crate::support::Run;

use std::collections::BTreeMap;
use std::iter::FusedIterator;

//-----------------------------------------------------------------------------

// FIXME: builder itself

//-----------------------------------------------------------------------------

/// A [`Run`] encoded using 32-bit integers to save space.
///
/// This is intended for situations, where we store a large number of runs explicitly.
/// The C++ implementation also uses 32-bit integers in similar situations.
#[derive(Copy, Clone, Debug, Default, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct SmallRun {
    // Value in the run.
    pub value: u32,
    // Length of the run.
    pub len: u32,
}

impl SmallRun {
    /// Creates a new run.
    pub fn new(value: usize, len: usize) -> Self {
        SmallRun {
            value: value as u32,
            len: len as u32,

        }
    }
}

impl From<(usize, usize)> for SmallRun {
    fn from(run: (usize, usize)) -> Self {
        Self::new(run.0, run.1)
    }
}

impl From<SmallRun> for Run {
    fn from(run: SmallRun) -> Self {
        Run::new(run.value as usize, run.len as usize)
    }
}

impl From<Run> for SmallRun {
    fn from(run: Run) -> Self {
        Self::new(run.value, run.len)
    }
}

/// A [`Pos`] encoded using 32-bit integers to save space.
///
/// This is intended for situations, where we store a large number of edges explicitly.
/// The C++ implementation also uses 32-bit integers in similar situations.
#[derive(Copy, Clone, Debug, Default, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct SmallPos {
    /// GBWT node identifier.
    pub node: u32,
    /// BWT offset within the node.
    pub offset: u32,
}

impl SmallPos {
    /// Creates a new position.
    pub fn new(node: usize, offset: usize) -> Self {
        SmallPos {
            node: node as u32,
            offset: offset as u32,
        }
    }
}

impl From<(usize, usize)> for SmallPos {
    fn from(pos: (usize, usize)) -> Self {
        Self::new(pos.0, pos.1)
    }
}

impl From<SmallPos> for Pos {
    fn from(pos: SmallPos) -> Self {
        Pos::new(pos.node as usize, pos.offset as usize)
    }
}

impl From<Pos> for SmallPos {
    fn from(pos: Pos) -> Self {
        Self::new(pos.node, pos.offset)
    }
}

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
/// edges.insert(3, 10);
/// edges.insert(1, 5);
/// edges.insert(2, 7);
/// edges.insert(1, 6); // Duplicate edge will be ignored.
/// edges.insert(4, 12); // Triggers conversion to large representation.
///
/// assert_eq!(edges.len(), 4);
/// assert_eq!(edges.get(1), Some(5)); // Original offset is preserved.
/// assert_eq!(edges.get(5), None);
///
/// let truth = vec![
///     Pos::new(1, 5),
///     Pos::new(2, 7),
///     Pos::new(3, 10),
///     Pos::new(4, 12)
/// ];
/// let edges: Vec<Pos> = edges.iter().collect();
/// assert_eq!(edges, truth);
/// ```
#[derive(Clone, Debug)]
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

    // Returns `true` if the given list contains the given node.
    fn contains(list: &[SmallPos], len: usize, node: usize) -> bool {
        for edge in list.iter().take(len) {
            if edge.node == node as u32 {
                return true;
            }
        }
        false
    }

    /// Inserts a new edge into the list.
    ///
    /// Does nothing if an edge to the same node already exists.
    pub fn insert(&mut self, node: usize, offset: usize) {
        let pos = SmallPos::new(node, offset);
        match self {
            EdgeList::Small(len, edges) => {
                if Self::contains(edges, *len as usize, node) {
                    return;
                }
                if *len < Self::SMALL_CAPACITY as u8 {
                    edges[*len as usize] = pos;
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
                    map.insert(pos.node, pos.offset);
                    *self = EdgeList::Large(map);
                }
            }
            EdgeList::Large(map) => {
                map.entry(pos.node).or_insert(pos.offset);
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

#[cfg(test)]
mod tests {
    use super::*;

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
            edges.insert(edge.node, edge.offset);
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
    fn edge_list_duplicates() {
        let mut truth = vec![
            Pos::new(3, 10),
            Pos::new(1, 5),
            Pos::new(2, 7),
            Pos::new(4, 12),
            Pos::new(10, 18),
        ];
        let mut edges = EdgeList::new();
        for edge in truth.iter() {
            edges.insert(edge.node, edge.offset);
            edges.insert(edge.node, edge.offset + 1); // Duplicate edge with different offset will be ignored.
        }
        truth.sort();
        check_edge_list(&edges, &truth);
    }
}

//-----------------------------------------------------------------------------
