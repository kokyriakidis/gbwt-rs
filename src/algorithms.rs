//! Algorithms using GBWT and GBZ.

use crate::{GBZ, Orientation, support};
use crate::support::Chains;

use std::collections::{BTreeMap, HashMap, HashSet, VecDeque};
use std::fmt::{Display, Formatter};
use std::{cmp, fmt};

#[cfg(test)]
mod tests;

//-----------------------------------------------------------------------------

/// Returns the longest common subsequence of integer sequences `a` and `b`, weighted by the given function.
///
/// The subsequence is returned as pairs of positions, and the second return value is the total weight of the LCS.
/// Weights are applied to the elements of the sequences, and they are assumed to be non-zero.
/// This version of the algorithm uses naive dynamic programming.
/// It is not suitable for long sequences.
/// See [`fast_weighted_lcs`] for an implementation based on Myers' O(nd) algorithm.
pub fn naive_weighted_lcs<F: Fn(usize) -> usize>(a: &[usize], b: &[usize], weight: F) -> (Vec<(usize, usize)>, usize) {
    let mut dp = vec![vec![0; b.len() + 1]; a.len() + 1];
    for (i, a_val) in a.iter().enumerate() {
        for (j, b_val) in b.iter().enumerate() {
            dp[i + 1][j + 1] = cmp::max(dp[i + 1][j], dp[i][j + 1]);
            if a_val == b_val {
                dp[i + 1][j + 1] = cmp::max(dp[i + 1][j + 1], dp[i][j] + weight(a[i]));
            }
        }
    }

    let mut result = Vec::new();
    let mut i = a.len();
    let mut j = b.len();
    while i > 0 && j > 0 {
        if a[i - 1] == b[j - 1] {
            result.push((i - 1, j - 1));
            i -= 1;
            j -= 1;
        } else if dp[i][j] == dp[i - 1][j] {
            i -= 1;
        } else {
            j -= 1;
        }
    }
    result.reverse();
    (result, dp[a.len()][b.len()])
}

//-----------------------------------------------------------------------------

// TODO: We do not need to store the weight. Due to the invariant, we can derive it from
// edits and the prefix sums.
// TODO: We could make this more space-efficient by using 32-bit integers.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct DPPoint {
    // Twice the total weight of the LCS up to this point.
    // This guarantees an invariant: `weight + edits = a_prefix_sum + b_prefix_sum`.
    weight: usize,

    // Non-weighted offset in the first sequence.
    a_offset: usize,

    // Non-weighted offset in the second sequence.
    b_offset: usize,

    // Length of the non-weighted run of matches.
    matches: usize,
}

impl DPPoint {
    fn new(weight: usize, a_offset: usize, b_offset: usize) -> DPPoint {
        DPPoint {
            weight, a_offset, b_offset, matches: 0,
        }
    }
}

impl Display for DPPoint {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "(weight {}, a {}, b {}, matches {})", self.weight, self.a_offset, self.b_offset, self.matches)
    }
}

struct DPMatrix<'a> {
    // First sequence.
    a: &'a [usize],

    // Second sequence.
    b: &'a [usize],

    // Prefix sum of weights on the first sequence.
    a_weights: Vec<usize>,

    // Prefix sum of weights on the second sequence.
    b_weights: Vec<usize>,

    // Furthest point for the given number of weighted edits on the weighted diagonal.
    // A diagonal is a_prefix_sum - b_prefix_sum.
    // The point refers to the first unprocessed pair of offsets.
    points: BTreeMap<(usize, isize), DPPoint>,
}

impl<'a> DPMatrix<'a> {
    fn prefix_sum<F: Fn(usize) -> usize>(sequence: &[usize], weight: F) -> Vec<usize> {
        let mut result = Vec::with_capacity(sequence.len() + 1);
        result.push(0);
        for i in 0..sequence.len() {
            result.push(result[i] + weight(sequence[i]));
        }
        result
    }

    fn new<F: Fn(usize) -> usize>(a: &'a [usize], b: &'a [usize], weight: F) -> Self {
        let a_weights = Self::prefix_sum(a, &weight);
        let b_weights = Self::prefix_sum(b, &weight);
        let mut points = BTreeMap::new();
        points.insert((0, 0), DPPoint::new(0, 0, 0));
        DPMatrix {
            a, b, a_weights, b_weights, points
        }
    }

    fn diagonals_for(&self, edits: usize) -> Vec<isize> {
        let mut result = Vec::new();
        for ((_, diagonal), _) in self.points.range((edits, isize::MIN)..(edits + 1, isize::MIN)) {
            result.push(*diagonal);
        }
        result
    }

    // Inserts the point if it is better than the existing one.
    fn try_insert(&mut self, edits: usize, diagonal: isize, point: DPPoint) {
        if let Some(existing) = self.points.get_mut(&(edits, diagonal)) {
            if point.weight > existing.weight {
                *existing = point;
            }
        } else {
            self.points.insert((edits, diagonal), point);
        }
    }

    fn a_weight(&self, offset: usize) -> usize {
        self.a_weights[offset + 1] - self.a_weights[offset]
    }

    fn b_weight(&self, offset: usize) -> usize {
        self.b_weights[offset + 1] - self.b_weights[offset]
    }

    // Extends matches over all diagonals for the given number of weighted edits.
    // Also adds successor points reachable with a single edit from each extension.
    // Returns the point if we reached the end.
    fn extend(&mut self, edits: usize) -> Option<DPPoint> {
        let diagonals = self.diagonals_for(edits);
        for diagonal in diagonals.into_iter() {
            let point = self.points.get(&(edits, diagonal)).copied();
            if point.is_none() {
                continue;
            }
            let mut point = point.unwrap();
            while point.a_offset < self.a.len() && point.b_offset < self.b.len() && self.a[point.a_offset] == self.b[point.b_offset] {
                point.weight += 2 * self.a_weight(point.a_offset);
                point.a_offset += 1;
                point.b_offset += 1;
                point.matches += 1;
            }
            if point.matches > 0 {
                self.points.insert((edits, diagonal), point);
            }
            if point.a_offset == self.a.len() && point.b_offset == self.b.len() {
                return Some(point);
            }
            if point.a_offset < self.a.len() {
                let weight = self.a_weight(point.a_offset);
                let new_point = DPPoint::new(point.weight, point.a_offset + 1, point.b_offset);
                self.try_insert(edits + weight, diagonal + (weight as isize), new_point);
            }
            if point.b_offset < self.b.len() {
                let weight = self.b_weight(point.b_offset);
                let new_point = DPPoint::new(point.weight, point.a_offset, point.b_offset + 1);
                self.try_insert(edits + weight, diagonal - (weight as isize), new_point);
            }
        }
        None
    }

    // Returns the next possible number of edits after the given number.
    // This assumes that `extend` has been called for the given number of edits.
    fn next_edits(&self, edits: usize) -> Option<usize> {
        if let Some(((value, _), _)) = self.points.range((edits + 1, isize::MIN)..).next() {
            Some(*value)
        } else {
            None
        }
    }

    // Returns the predecessor point and number of edits for the given point and number of edits.
    fn predecessor(&self, a_offset: usize, b_offset: usize, edits: usize) -> Option<(DPPoint, usize)> {
        let diagonal = (self.a_weights[a_offset] as isize) - (self.b_weights[b_offset] as isize);
        let prev = if a_offset > 0 && self.a_weight(a_offset - 1) <= edits {
            let weight = self.a_weight(a_offset - 1);
            self.points.get(&(edits - weight, diagonal - (weight as isize))).copied()
        } else {
            None
        };
        let next = if b_offset > 0 && self.b_weight(b_offset - 1) <= edits {
            let weight = self.b_weight(b_offset - 1);
            self.points.get(&(edits - weight, diagonal + (weight as isize))).copied()
        } else {
            None
        };
        match (prev, next) {
            (Some(p), Some(n)) => {
                if p.weight > n.weight {
                    Some((p, edits - self.a_weight(a_offset - 1)))
                } else {
                    Some((n, edits - self.b_weight(b_offset - 1)))
                }
            },
            (Some(p), None) => Some((p, edits - self.a_weight(a_offset - 1))),
            (None, Some(n)) => Some((n, edits - self.b_weight(b_offset - 1))),
            (None, None) => None,
        }
    }
}

/// Returns the longest common subsequence of integer sequences `a` and `b`, weighted by the given function.
///
/// The subsequence is returned as pairs of positions, and the second return value is the total weight of the LCS.
/// Weights are applied to the elements of the sequences, and they are assumed to be non-zero.
/// This version of the algorithm is based on Myers' O(nd) algorithm.
/// It is efficient with long sequences, as long as they are similar.
/// See [`naive_weighted_lcs`] for an implementation using naive dynamic programming.
///
/// The following specializations are available:
///
/// * [`lcs`] for unweighted LCS.
/// * [`path_lcs`] for paths in a graph, using sequence lengths as weights.
pub fn fast_weighted_lcs<F: Fn(usize) -> usize>(a: &[usize], b: &[usize], weight: F) -> (Vec<(usize, usize)>, usize) {
    if a.is_empty() || b.is_empty() {
        return (Vec::new(), 0);
    }

    // Find the furthest point on each diagonal with each possible number of edits, until we reach the end.
    let mut matrix = DPMatrix::new(a, b, &weight);
    let mut edits = 0;
    let mut point = DPPoint::new(0, 0, 0);
    loop {
        if let Some(next_point) = matrix.extend(edits) {
            point = next_point;
            break;
        }
        if let Some(next_edits) = matrix.next_edits(edits) {
            edits = next_edits;
        } else {
            // This should not happen.
            break;
        }
    }

    // Trace back the LCS.
    let mut result = Vec::new();
    let final_weight = point.weight / 2;
    loop {
        for _ in 0..point.matches {
            point.a_offset -= 1;
            point.b_offset -= 1;
            result.push((point.a_offset, point.b_offset));
        }
        if let Some((p, e)) = matrix.predecessor(point.a_offset, point.b_offset, edits) {
            point = p;
            edits = e;
        } else {
            break;
        }
    }
    result.reverse();
    (result, final_weight)
}

//-----------------------------------------------------------------------------

/// Returns the longest common subsequence of `a` and `b`.
///
/// The return value consists of pairs of positions in the input vectors.
///
/// # Examples
///
/// ```
/// use gbz::algorithms::lcs;
///
/// let a = vec![1, 2, 3, 4, 5];
/// let b = vec![2, 4, 6, 8, 10];
/// let truth = vec![(1, 0), (3, 1)];
/// assert_eq!(lcs(&a, &b), truth);
/// ``````
pub fn lcs(a: &[usize], b: &[usize]) -> Vec<(usize, usize)> {
    fast_weighted_lcs(a, b, |_| 1).0
}

/// Returns the longest common subsequence of paths `a` and `b` in the graph.
///
/// The LCS is weighted by the length of the node and returned as a vector of pairs of positions.
/// The second return value is the total weight of the LCS.
/// If a node is not found in the graph, its length is assumed to be zero.
/// In such cases, the LCS may not be meaningful.
///
/// # Examples
///
/// ```
/// use gbz::{GBZ, Orientation, algorithms};
/// use gbz::support;
/// use simple_sds::serialize;
///
/// let filename = support::get_test_data("translation.gbz");
/// let gbz: GBZ = serialize::load_from(&filename).unwrap();
///
/// fn get_path(graph: &GBZ, path_id: usize) -> Vec<usize> {
///     graph.path(path_id, Orientation::Forward)
///         .unwrap()
///         .map(|(id, o)| support::encode_node(id, o))
///         .collect()
/// }
///
/// // (1: GA), (2: T), (3: T), (5: CA), (6: G), (9: A), (11: TA)
/// let a = get_path(&gbz, 1);
///
/// // (1: GA), (2: T), (4: A), (5: CA), (6: G), (10: T), (11: TA)
/// let b = get_path(&gbz, 2);
///
/// let truth = vec![
///     (0, 0), (1, 1), (3, 3), (4, 4), (6, 6)
/// ];
/// let len = 8;
///
/// assert_eq!(algorithms::path_lcs(&a, &b, &gbz), (truth, len));
/// ```
pub fn path_lcs(a: &[usize], b: &[usize], graph: &GBZ) -> (Vec<(usize, usize)>, usize) {
    let weight = |handle| graph.sequence_len(support::node_id(handle)).unwrap_or(0);
    fast_weighted_lcs(a, b, weight)
}

//-----------------------------------------------------------------------------

// TODO: chains.rs with the `Chains` struct and related code.
// TODO: We could have a version that returns the actual chains (including the non-boundary nodes?).
// We could also compute the minimum distance (in bp) from the start of the chain to each node.
// That would be good enough to replace the distance index in haplotype sampling.
// TODO: Chain as a vector of (handle, is_snarl_start, is_snarl_end, distance_from_start, distance_to_end) tuples?

/// Finds top-level chains in the graph.
///
/// A top-level chain is a sequence of boundary nodes bordering snarls in a weakly connected component.
/// The result contains links between consecutive boundary nodes in each chain.
///
/// Each component must have exactly two tip nodes.
/// All articulation points in underlying undirected graph of the component must be between the tips.
/// If these assumptions do not hold, the component is skipped and not included in the result.
/// In more complex cases, the chains can be extracted from a snarl decomposition or a distance index using `vg chains`.
///
/// # Examples
///
/// ```
/// use gbz::{GBZ, Orientation};
/// use gbz::{algorithms, support};
/// use simple_sds::serialize;
///
/// let filename = support::get_test_data("translation.gbz");
/// let gbz: GBZ = serialize::load_from(&filename).unwrap();
///
/// let chains = algorithms::find_chains(&gbz);
/// assert_eq!(chains.len(), 1);
/// assert_eq!(chains.links(), 3);
/// assert_eq!(
///     chains.next(support::encode_node(2, Orientation::Forward)),
///     Some(support::encode_node(5, Orientation::Forward))
/// );
/// assert_eq!(
///     chains.next(support::encode_node(5, Orientation::Forward)),
///     Some(support::encode_node(6, Orientation::Forward))
/// );
/// assert_eq!(
///    chains.next(support::encode_node(6, Orientation::Forward)),
///    Some(support::encode_node(11, Orientation::Forward))
/// );
/// ```
pub fn find_chains(graph: &GBZ) -> Chains {
    let mut chains = Chains::new();
    let components = graph.weakly_connected_components();

    for component in components.iter() {
        let tips = find_tips(graph, component);
        if tips.len() != 2 {
            continue;
        }

        // Find articulation points (bridges) in the underlying undirected graph of the component.
        // Then add the tips to get the set of chain nodes.
        let mut chain_node_ids = find_articulation_points(graph, component);
        for tip in tips.iter() {
            chain_node_ids.insert(support::node_id(*tip));
        }

        // This will return empty vectors if some chain nodes were not visited.
        let (chain_handles, is_snarl) = walk_chain(graph, tips[0], &chain_node_ids);
        if chain_handles.is_empty() {
            continue;
        }

        // Extract boundary nodes: chain nodes where at least one adjacent edge is a snarl.
        // is_snarl[i] indicates whether the link between chain_handles[i] and chain_handles[i + 1] is a snarl.
        let mut is_boundary = is_snarl;
        for i in (0..(is_boundary.len() - 1)).rev() {
            if is_boundary[i] {
                is_boundary[i + 1] = true;
            }
        }

        // Add links between consecutive boundary nodes.
        let boundary: Vec<usize> = (0..chain_handles.len())
            .filter(|&i| is_boundary[i])
            .map(|i| chain_handles[i])
            .collect();
        for pair in boundary.windows(2) {
            chains.add_link(pair[0], pair[1]);
        }
    }

    chains.count_chains();
    chains
}

// Returns inward-oriented handles for tip nodes in the component.
// A tip node has successors on one side but not the other.
fn find_tips(graph: &GBZ, component: &[usize]) -> Vec<usize> {
    let mut tips = Vec::new();
    for &node_id in component {
        let fwd_succ = graph.successors(node_id, Orientation::Forward)
            .map_or(0, |iter| iter.count());
        let fwd_pred = graph.predecessors(node_id, Orientation::Forward)
            .map_or(0, |iter| iter.count());

        if fwd_succ == 0 && fwd_pred == 0 {
            continue;
        }
        if fwd_succ > 0 && fwd_pred > 0 {
            continue;
        }
        // Tip: has edges on one side only. Inward handle is the orientation with successors.
        if fwd_succ > 0 {
            tips.push(support::encode_node(node_id, Orientation::Forward));
        } else {
            tips.push(support::encode_node(node_id, Orientation::Reverse));
        }
    }
    tips
}

// Finds articulation points in the underlying undirected graph of a component using Tarjan's algorithm.
// Returns node ids (not handles).
fn find_articulation_points(graph: &GBZ, component: &[usize]) -> HashSet<usize> {
    let component_set: HashSet<usize> = component.iter().copied().collect();
    let mut index_map: HashMap<usize, usize> = HashMap::new(); // Node id to index in component vector.
    for (i, &node_id) in component.iter().enumerate() {
        index_map.insert(node_id, i);
    }

    let n = component.len();
    let mut discovered = vec![0usize; n]; // Discovery time (rank) of each node in DFS.
    let mut low = vec![0usize; n]; // Lowest discovery time reachable from the node via DFS tree edges and back edges.
    let mut visited = vec![false; n];
    let mut parent = vec![usize::MAX; n];
    let mut result = HashSet::new();
    let mut timer = 0usize;

    // Get undirected neighbors of a node (all nodes reachable via successors or predecessors in either orientation).
    let get_neighbors = |node_id: usize| -> Vec<usize> {
        let mut nbrs = HashSet::new();
        for orientation in [Orientation::Forward, Orientation::Reverse] {
            if let Some(iter) = graph.successors(node_id, orientation) {
                for (next_id, _) in iter {
                    if component_set.contains(&next_id) {
                        nbrs.insert(next_id);
                    }
                }
            }
        }
        nbrs.into_iter().collect()
    };

    // Iterative Tarjan's algorithm for articulation points.
    // Stack entries: (node index, neighbor list, next neighbor).
    let start_index = 0;
    let mut stack: Vec<(usize, Vec<usize>, usize)> = Vec::new();
    visited[start_index] = true;
    discovered[start_index] = timer;
    low[start_index] = timer;
    timer += 1;
    let nbrs = get_neighbors(component[start_index]);
    stack.push((start_index, nbrs, 0));

    while let Some((from_index, neighbors, next_neighbor)) = stack.last_mut() {
        if *next_neighbor < neighbors.len() {
            let to_id = neighbors[*next_neighbor];
            *next_neighbor += 1;
            let to_index = index_map[&to_id];
            if !visited[to_index] {
                visited[to_index] = true;
                discovered[to_index] = timer;
                low[to_index] = timer;
                timer += 1;
                parent[to_index] = *from_index;
                let to_nbrs = get_neighbors(component[to_index]);
                stack.push((to_index, to_nbrs, 0));
            } else if to_index != parent[*from_index] {
                low[*from_index] = cmp::min(low[*from_index], discovered[to_index]);
            }
        } else {
            // Done processing node from, pop and update parent.
            let from_index = *from_index;
            stack.pop();
            if let Some((parent_index, _, _)) = stack.last() {
                low[*parent_index] = cmp::min(low[*parent_index], low[from_index]);
                // Check if parent_index is an articulation point via child from_index.
                if parent[*parent_index] == usize::MAX {
                    // parent_index is the root; it's an articulation point (or a tip) if it has 2+ children in DFS tree.
                    let child_count = parent.iter().filter(|&&p| p == *parent_index).count();
                    if child_count > 1 {
                        result.insert(component[*parent_index]);
                    }
                } else if low[from_index] >= discovered[*parent_index] {
                    result.insert(component[*parent_index]);
                }
            }
        }
    }

    result
}

// Walks the chain from the start tip to the other tip, classifying links between nodes as simple or snarl.
// Returns (chain_handles, is_snarl) where chain_handles are the oriented node visits.
// is_snarl[i] indicates whether the link between chain_handles[i] and chain_handles[i+1] is a snarl.
fn walk_chain(
    graph: &GBZ,
    start_handle: usize,
    chain_node_ids: &HashSet<usize>
) -> (Vec<usize>, Vec<bool>) {
    let mut chain_handles = Vec::with_capacity(chain_node_ids.len());
    let mut is_snarl = Vec::with_capacity(chain_node_ids.len());

    let mut current = start_handle;
    chain_handles.push(current);

    loop {
        let cur_id = support::node_id(current);
        let cur_orient = support::node_orientation(current);

        let succs: Vec<(usize, Orientation)> = graph.successors(cur_id, cur_orient)
            .map_or(Vec::new(), |iter| iter.collect());
        if succs.is_empty() {
            break;
        }

        // Check if this is a simple bridge edge: single successor with single predecessor.
        if succs.len() == 1 {
            let (next_id, next_o) = succs[0];
            let next_handle = support::encode_node(next_id, next_o);
            let pred_count = graph.predecessors(next_id, next_o)
                .map_or(0, |iter: crate::gbz::EdgeIter<'_>| iter.len());
            if pred_count == 1 {
                is_snarl.push(false);
                chain_handles.push(next_handle);
                current = next_handle;
                continue;
            }
        }

        // This is a snarl. BFS to find the next chain node where all paths converge.
        let next_chain_handle = find_snarl_end(graph, current, chain_node_ids);
        is_snarl.push(true);
        chain_handles.push(next_chain_handle);
        current = next_chain_handle;
    }
    is_snarl.push(false); // We will reuse this vector to mark boundary nodes, so add a dummy value for the last node.

    // If we did not visit every tip and bridge node, the assumptions do not hold.
    // Return empty vectors to skip this component.
    if chain_handles.len() == chain_node_ids.len() {
        (chain_handles, is_snarl)
    } else {
        (Vec::new(), Vec::new())
    }
}

// BFS from the current handle's successors to find the next chain/tip node.
fn find_snarl_end(
    graph: &GBZ,
    current_handle: usize,
    chain_node_ids: &HashSet<usize>
) -> usize {
    let cur_id = support::node_id(current_handle);
    let cur_orient = support::node_orientation(current_handle);

    let mut queue: VecDeque<usize> = VecDeque::new();
    let mut visited: HashSet<usize> = HashSet::new();
    visited.insert(current_handle);

    // Seed BFS with successors of the current handle.
    if let Some(iter) = graph.successors(cur_id, cur_orient) {
        for (next_id, next_o) in iter {
            let next_handle = support::encode_node(next_id, next_o);
            if !visited.contains(&next_handle) {
                visited.insert(next_handle);
                queue.push_back(next_handle);
            }
        }
    }

    while let Some(handle) = queue.pop_front() {
        let h_id = support::node_id(handle);
        let h_orient = support::node_orientation(handle);

        // If this is a chain node (bridge or tip) different from current, it's our target.
        if chain_node_ids.contains(&h_id) && h_id != cur_id {
            return handle;
        }

        if let Some(iter) = graph.successors(h_id, h_orient) {
            for (next_id, next_o) in iter {
                let next_handle = support::encode_node(next_id, next_o);
                if !visited.contains(&next_handle) {
                    visited.insert(next_handle);
                    queue.push_back(next_handle);
                }
            }
        }
    }

    // Should not reach here in a valid graph with exactly 2 tips.
    current_handle
}

//-----------------------------------------------------------------------------
