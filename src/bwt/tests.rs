use super::*;

use simple_sds::serialize;

use std::collections::HashSet;

//-----------------------------------------------------------------------------

fn get_edge_lists(edges: Vec<Vec<Pos>>) -> Vec<EdgeList> {
    edges.iter().map(|record_edges| {
        let mut edge_list = EdgeList::new();
        for edge in record_edges {
            edge_list.increment(edge.node, edge.offset);
        }
        edge_list
    }).collect()
}

// GBWT example from the paper: (endmarker edges, endmarker, edges, runs, invalid_node)
fn get_edges_runs() -> (EdgeList, Vec<Pos>, Vec<EdgeList>, Vec<Vec<Run>>, usize) {
    let mut endmarker_edges = EdgeList::new();
    endmarker_edges.increment(1, 0);
    let endmarker = vec![Pos::new(1, 0), Pos::new(1, 1), Pos::new(1, 2)];

    let edges = vec![
        vec![Pos::new(2, 0), Pos::new(3, 0)],
        vec![Pos::new(4, 0), Pos::new(5, 0)],
        vec![Pos::new(4, 1)],
        vec![Pos::new(5, 1), Pos::new(6, 0)],
        vec![Pos::new(7, 0)],
        vec![Pos::new(7, 2)],
        vec![Pos::new(0, 0)],
    ];
    let edges = get_edge_lists(edges);

    let runs = vec![
        vec![Run::new(2, 2), Run::new(3, 1)],
        vec![Run::new(4, 1), Run::new(5, 1)],
        vec![Run::new(4, 1)],
        vec![Run::new(5, 1), Run::new(6, 1)],
        vec![Run::new(7, 2)],
        vec![Run::new(7, 1)],
        vec![Run::new(0, 3)],
    ];

    (endmarker_edges, endmarker, edges, runs, 8)
}

// Bidirectional version of the example.
fn get_bidirectional() -> (EdgeList, Vec<Pos>, Vec<EdgeList>, Vec<Vec<Run>>, usize) {
    let mut endmarker_edges = EdgeList::new();
    endmarker_edges.increment(2, 0);
    endmarker_edges.increment(15, 0);
    let endmarker = vec![
        Pos::new(2, 0), Pos::new(2, 1), Pos::new(2, 2),
        Pos::new(15, 0), Pos::new(15, 1), Pos::new(15, 2),
    ];

    let edges = vec![
        // 1
        vec![Pos::new(4, 0), Pos::new(6, 0)],
        vec![Pos::new(0, 0)],
        // 2
        vec![Pos::new(8, 0), Pos::new(10, 0)],
        vec![Pos::new(3, 0)],
        // 3
        vec![Pos::new(8, 1)],
        vec![Pos::new(3, 2)],
        // 4
        vec![Pos::new(10, 1), Pos::new(12, 0)],
        vec![Pos::new(5, 0), Pos::new(7, 0)],
        // 5
        vec![Pos::new(14, 0)],
        vec![Pos::new(5, 1), Pos::new(9, 0)],
        // 6
        vec![Pos::new(14, 2)],
        vec![Pos::new(9, 1)],
        // 7
        vec![Pos::new(0, 0)],
        vec![Pos::new(11, 0), Pos::new(13, 0)],
    ];
    let edges = get_edge_lists(edges);

    let runs = vec![
        // 1
        vec![Run::new(4, 2), Run::new(6, 1)],
        vec![Run::new(0, 3)],
        // 2
        vec![Run::new(8, 1), Run::new(10, 1)],
        vec![Run::new(3, 2)],
        // 3
        vec![Run::new(8, 1)],
        vec![Run::new(3, 1)],
        // 4
        vec![Run::new(10, 1), Run::new(12, 1)],
        vec![Run::new(5, 1), Run::new(7, 1)],
        // 5
        vec![Run::new(14, 2)],
        vec![Run::new(5, 1), Run::new(9, 1)],
        // 6
        vec![Run::new(14, 1)],
        vec![Run::new(9, 1)],
        // 7
        vec![Run::new(0, 3)],
        vec![Run::new(11, 1), Run::new(13, 2)],
    ];

    (endmarker_edges, endmarker, edges, runs, 16)
}

fn create_bwt(endmarker_edges: &EdgeList, endmarker: &[Pos], edges: &[EdgeList], runs: &[Vec<Run>]) -> BWT {
    let mut expected_records = if endmarker.is_empty() { 0 } else { 1 };
    let mut builder = BWTBuilder::new(endmarker_edges, endmarker);
    assert_eq!(builder.len(), expected_records, "Wrong number of records in a new builder");
    assert_eq!(builder.is_empty(), endmarker.is_empty(), "Builder should be empty iff the endmarker is empty");

    for i in 0..edges.len() {
        builder.append(&edges[i], runs[i].iter().copied());
        expected_records += 1;
    }
    assert_eq!(builder.len(), expected_records, "Wrong number of records in the builder");
    assert_eq!(builder.is_empty(), endmarker.is_empty(), "Builder should be empty after construction iff the endmarker is empty");

    BWT::from(builder)
}

// Check records in the BWT, using the provided edges as the source of truth.
// Also checks that `id()` works correctly.
fn check_records(bwt: &BWT, endmarker_edges: &EdgeList, edges: &[EdgeList]) {
    let expected_records = edges.len() + if endmarker_edges.is_empty() { 0 } else { 1 };
    assert_eq!(bwt.len(), expected_records, "Invalid number of records in the BWT");
    assert_eq!(bwt.is_empty(), expected_records == 0, "Invalid BWT emptiness");

    // Edges.
    for i in 0..bwt.len() {
        let record = bwt.record(i);
        let curr_edges = if i == 0 { endmarker_edges } else { &edges[i - 1] };
        assert_eq!(record.is_none(), curr_edges.is_empty(), "Invalid record {} existence", i);
        if let Some(record) = record {
            assert_eq!(record.id(), i, "Invalid id for record {}", i);
            assert_eq!(record.outdegree(), curr_edges.len(), "Invalid outdegree in record {}", i);
            for (j, edge) in curr_edges.iter().enumerate() {
                assert_eq!(record.successor(j), edge.node, "Invalid successor {} in record {}", j, i);
                assert_eq!(record.offset(j), edge.offset, "Invalid offset {} in record {}", j, i);
            }
        }

        // Compressed record.
        let compressed = bwt.compressed_record(i);
        assert_eq!(compressed.is_none(), curr_edges.is_empty(), "Invalid compressed record {} existence", i);
        if let Some((edge_bytes, bwt_bytes)) = compressed {
            let decompressed = Record::decompress_edges(edge_bytes);
            assert!(decompressed.is_some(), "Could not decompress edges for record {}", i);
            let (edges, offset) = decompressed.unwrap();
            assert_eq!(offset, edge_bytes.len(), "Invalid offset after edge list for record {}", i);
            assert!(edges.iter().copied().eq(curr_edges.iter()), "Invalid edges in compressed record {}", i);
            let record = bwt.record(i).unwrap();
            assert_eq!(bwt_bytes, record.bwt, "Invalid BWT in compressed record {}", i);
        }
    }
}

// Check that the iterator finds the correct records and that the id iterator finds the same ids.
fn check_iter(bwt: &BWT) {
    let mut iter = bwt.iter();
    let mut id_iter = bwt.id_iter();
    for i in 0..bwt.len() {
        if let Some(truth) = bwt.record(i) {
            if let Some(record) = iter.next() {
                assert_eq!(record.id(), truth.id(), "Invalid record id from the iterator");
            } else {
                panic!("Iterator did not find record {}", i);
            }
            assert_eq!(id_iter.next(), Some(truth.id()), "Invalid id from id iterator");
        }
    }
    assert!(iter.next().is_none(), "Iterator found a record past the end");
    assert!(id_iter.next().is_none(), "Id iterator found a record past the end");
}

// Check all `lf()` results in the BWT, using the provided edges and runs as the source of truth.
// Then check that decompressing the record works correctly.
// Also checks that `offset_to()` works in positive cases and that `len()` is correct.
fn check_lf(bwt: &BWT, endmarker_edges: &EdgeList, endmarker: &[Pos], edges: &[EdgeList], runs: &[Vec<Run>]) {
    let endmarker_runs = BWTBuilder::get_runs(endmarker);

    // `lf()` at each offset of each record.
    for i in 0..bwt.len() {
        if let Some(record) = bwt.record(i) {
            let mut offset = 0;
            let (mut curr_edges, curr_runs) = if i == 0 {
                (endmarker_edges.clone(), &endmarker_runs)
            } else {
                (edges[i - 1].clone(), &runs[i - 1])
            };
            let decompressed = record.decompress();
            assert_eq!(decompressed.len(), record.len(), "Invalid decompressed record {} length", i);
            for run in curr_runs {
                for _ in 0..run.len {
                    let edge = Pos::new(run.value, curr_edges.get(run.value).unwrap() as usize);
                    let expected = if edge.node == ENDMARKER { None } else { Some(edge) };
                    assert_eq!(record.lf(offset), expected, "Invalid lf({}) in record {}", offset, i);
                    assert_eq!(decompressed[offset], edge, "Invalid decompressed lf({}) in record {}", offset, i);
                    let expected = if edge.node == ENDMARKER { None } else { Some(offset) };
                    assert_eq!(record.offset_to(edge), expected, "Invalid offset_to(({}, {})) in record {}", edge.node, edge.offset, i);
                    offset += 1;
                    curr_edges.increment(edge.node, 1);
                }
            }
            assert_eq!(record.len(), offset, "Invalid record {} length", i);
            assert_eq!(record.lf(offset), None, "Got an lf() result past the end in record {}", i);
        }
    }
}

// Check all `follow()` results in the BWT, using `lf()` as the source of truth.
// Also checks that `bd_follow()` returns the same ranges.
// The tests for bidirectional search in `GBWT` make sure that the second return values are correct.
fn check_follow(bwt: &BWT, invalid_node: usize) {
    for record in bwt.iter() {
        let i = record.id();
        // Check all ranges, including empty and past-the-end ones.
        let len = record.len();
        for start in 0..len + 1 {
            for limit in start..len + 1 {
                // With an endmarker.
                assert_eq!(record.follow(start..limit, ENDMARKER), None, "Got a follow({}..{}, endmarker) result in record {}", start, limit, i);
                assert_eq!(record.bd_follow(start..limit, ENDMARKER), None, "Got a bd_follow({}..{}, endmarker) result in record {}", start, limit, i);

                // With each successor node.
                for rank in 0..record.outdegree() {
                    let successor = record.successor(rank);
                    if successor == ENDMARKER {
                        continue;
                    }
                    if let Some(result) = record.follow(start..limit, successor) {
                        let mut found = result.start..result.start;
                        for j in start..limit {
                            if let Some(pos) = record.lf(j) {
                                if pos.node == successor && pos.offset == found.end {
                                    found.end += 1;
                                }
                            }
                        }
                        assert_eq!(result, found, "follow({}..{}, {}) did not find the correct range in record {}", start, limit, successor, i);
                        if let Some((bd_result, _)) =  record.bd_follow(start..limit, successor) {
                            assert_eq!(bd_result, result, "bd_follow({}..{}, {}) did not find the same range as follow() in record {}", start, limit, successor, i);
                        } else {
                            panic!("bd_follow({}..{}, {}) did not find a result in record {}", start, limit, successor, i);
                        }
                    } else {
                        for j in start..limit {
                            if let Some(pos) = record.lf(j) {
                                assert_ne!(pos.node, successor, "follow({}..{}, {}) did not follow offset {} in record {}", start, limit, successor, j, i);
                            }
                            assert_eq!(record.bd_follow(start..limit, successor), None, "Got a bd_follow({}..{}, {}) result in record {}", start, limit, successor, i);
                        }
                    }
                }

                // With an invalid node.
                assert_eq!(record.follow(start..limit, invalid_node), None, "Got a follow({}..{}, invalid) result in record {}", start, limit, i);
                assert_eq!(record.bd_follow(start..limit, invalid_node), None, "Got a bd_follow({}..{}, invalid) result in record {}", start, limit, i);
            }
        }
    }
}

// Check negative cases for `offset_to()`.
fn negative_offset_to(bwt: &BWT, invalid_node: usize) {
    for record in bwt.iter() {
        assert_eq!(record.offset_to(Pos::new(ENDMARKER, 0)), None, "Got an offset to the endmarker from record {}", record.id());
        assert_eq!(record.offset_to(Pos::new(invalid_node, 0)), None, "Got an offset to an invalid node from record {}", record.id());
        for rank in 0..record.outdegree() {
            let successor = record.successor(rank);
            if successor == ENDMARKER {
                continue;
            }
            let offset = record.offset(rank);
            if offset > 0 {
                assert_eq!(record.offset_to(Pos::new(successor, offset - 1)), None, "Got an offset from record {} to a too small position in {}", record.id(), successor);
            }
            let count = record.follow(0..record.len(), successor).unwrap().len();
            assert_eq!(record.offset_to(Pos::new(successor, offset + count)), None, "Got an offset from record {} to a too large position in {}", record.id(), successor);
        }
    }
}

// Check that we can find predecessors for all positions except starting positions.
// The tests for `GBWT::backward()` will make sure that the predecessors are correct.
fn check_predecessor_at(bwt: &BWT) {
    let mut starting_positions = HashSet::<Pos>::new();
    let endmarker = bwt.record(ENDMARKER).unwrap();
    for i in 0..endmarker.len() {
        starting_positions.insert(endmarker.lf(i).unwrap());
    }

    for record in bwt.iter() {
        if record.id() == ENDMARKER {
            continue;
        }
        let reverse_id = ((record.id() + 1) ^ 1) - 1; // Record to node, flip, node to record.
        let reverse_record = bwt.record(reverse_id).unwrap();
        for i in 0..record.len() {
            if starting_positions.contains(&Pos::new(record.id() + 1, i)) {
                assert!(reverse_record.predecessor_at(i).is_none(), "Found a predecessor for a starting position ({}, {})", record.id() + 1, i);
            } else {
                assert!(reverse_record.predecessor_at(i).is_some(), "Did not find a predecessor for position ({}, {})", record.id() + 1, i);
            }
        }
        assert!(reverse_record.predecessor_at(record.len()).is_none(), "Found a predecessor for an invalid offset at node {}", record.id() + 1);
    }
}

//-----------------------------------------------------------------------------

#[test]
fn empty_bwt() {
    let endmarker_edges = EdgeList::new();
    let endmarker = Vec::new();
    let edges = Vec::new();
    let runs = Vec::new();
    let invalid_node = 0;
    let bwt = create_bwt(&endmarker_edges, &endmarker, &edges, &runs);
    check_records(&bwt, &endmarker_edges, &edges);
    check_iter(&bwt);
    check_lf(&bwt, &endmarker_edges, &endmarker, &edges, &runs);
    check_follow(&bwt, invalid_node);
    negative_offset_to(&bwt, invalid_node);
    serialize::test(&bwt, "empty-bwt", None, true);
}

#[test]
fn non_empty_bwt() {
    let (endmarker_edges, endmarker, edges, runs, invalid_node) = get_edges_runs();
    let bwt = create_bwt(&endmarker_edges, &endmarker, &edges, &runs);
    check_records(&bwt, &endmarker_edges, &edges);
    check_iter(&bwt);
    check_lf(&bwt, &endmarker_edges, &endmarker, &edges, &runs);
    check_follow(&bwt, invalid_node);
    negative_offset_to(&bwt, invalid_node);
    serialize::test(&bwt, "non-empty-bwt", None, true);
}

#[test]
fn empty_records() {
    let (endmarker_edges, endmarker, mut edges, mut runs, invalid_node) = get_edges_runs();
    edges[2] = EdgeList::new();
    edges[6] = EdgeList::new();
    runs[2] = Vec::new();
    runs[6] = Vec::new();
 
    let bwt = create_bwt(&endmarker_edges, &endmarker, &edges, &runs);
    check_records(&bwt, &endmarker_edges, &edges);
    check_iter(&bwt);
    check_lf(&bwt, &endmarker_edges, &endmarker, &edges, &runs);
    check_follow(&bwt, invalid_node);
    negative_offset_to(&bwt, invalid_node);
    serialize::test(&bwt, "bwt-with-empty", None, true);
}

#[test]
fn bidirectional_bwt() {
    let (endmarker_edges, endmarker, edges, runs, invalid_node) = get_bidirectional();
    let bwt = create_bwt(&endmarker_edges, &endmarker, &edges, &runs);
    check_records(&bwt, &endmarker_edges, &edges);
    check_iter(&bwt);
    check_lf(&bwt, &endmarker_edges, &endmarker, &edges, &runs);
    check_follow(&bwt, invalid_node);
    negative_offset_to(&bwt, invalid_node);
    check_predecessor_at(&bwt);
    serialize::test(&bwt, "bidirectional-bwt", None, true);
}

//-----------------------------------------------------------------------------
