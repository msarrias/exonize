import pytest
from pathlib import Path
import networkx as nx
from exonize.exonize import Exonize
import portion as P


@pytest.fixture(scope="class")
def exonize_obj():
    # Initialize and return exonize_obj
    return Exonize(
        gff_file_path=Path('mock_gff.gff3'),
        genome_file_path=Path('mock_genome.fa'),
        gene_annot_feature='gene',
        cds_annot_feature='CDS',
        transcript_annot_feature='mRNA',
        sequence_base=1,
        frame_base=0,
        evalue_threshold=0.01,
        self_hit_threshold=0.5,
        query_coverage_threshold=0.8,
        min_exon_length=20,
        exon_clustering_overlap_threshold=0.8,
        targets_clustering_overlap_threshold=0.9,
        output_prefix="mock_specie",
        csv=False,
        enable_debug=False,
        soft_force=False,
        hard_force=False,
        sleep_max_seconds=0,
        cpus_number=1,
        timeout_database=60,
        output_directory_path=Path("."),
    )


@pytest.fixture(scope="class")
def gene_graph(exonize_obj):
    # Setup necessary sets and generate the gene_graph
    gene_cds_set = {
        P.open(0, 100),
        P.open(200, 300),
        P.open(300, 310),
        P.open(320, 380),
        P.open(500, 600),
        P.open(800, 900),
        P.open(1800, 2100),
        P.open(2200, 2350)
    }
    cds_candidates_coordinates = {
        P.open(0, 100),
        P.open(200, 300),
        P.open(320, 380),
        P.open(500, 600),
        P.open(800, 900),
        P.open(1800, 2100),
        P.open(2200, 2350)
    }
    local_records_set = {
        (1, 'gene1', 0, 100, 205, 300, 1e-3),
        (2, 'gene1', 200, 300, 4, 99, 1e-3),
        (3, 'gene1', 800, 900, 6, 101, 1e-3),
        (4, 'gene1', 500, 600, 700, 800, 1e-3),
        (5, 'gene1', 800, 900, 700, 800, 1e-4),
        (6, 'gene1', 1800, 2100, 2210, 2300, 1e-4)
    }
    global_records_set = {
        (1, 0, 100, 200, 300),
        (3, 0, 100, 800, 900),
        (3, 200, 300, 800, 900),
        (3, 0, 100, 500, 600),
        (4, 1800, 2100, 2200, 2350)
    }
    query_coords, target_coords = exonize_obj.event_reconciler.get_hits_query_and_target_coordinates(
        local_records_set=local_records_set
    )
    clusters = exonize_obj.data_container.get_overlapping_clusters(
        target_coordinates_set=target_coords,
        threshold=exonize_obj.environment.targets_clustering_overlap_threshold
    )
    targets_reference_coordinates_dictionary = exonize_obj.event_reconciler.get_matches_reference_mode_dictionary(
        clusters_list=clusters,
        cds_candidates_set=cds_candidates_coordinates,
        gene_cds_set=gene_cds_set
    )
    gene_graph = exonize_obj.event_reconciler.create_events_multigraph(
        targets_reference_coordinates_dictionary=targets_reference_coordinates_dictionary,
        query_local_coordinates_set=query_coords,
        local_records_set=local_records_set,
        global_records_set=global_records_set
    )
    return gene_graph, targets_reference_coordinates_dictionary


def test_build_reference_dictionary(
        exonize_obj
):
    query_coordinates = {
        P.open(0, 100),
        P.open(200, 250)
    }
    target_coordinates = {
        (P.open(0, 50), 0.9),
        (P.open(0, 48), 0.93),
        (P.open(40, 90), 0.8),
        (P.open(200, 250), 0.7),
        (P.open(210, 250), 0.7),
        (P.open(220, 250), 0.7),
        (P.open(220, 270), 0.6),
        (P.open(215, 270), 0.7),
        (P.open(219, 270), 0.8),
        (P.open(400, 450), 0.4),
        (P.open(402, 450), 0.5),
        (P.open(420, 450), 0.5)
    }
    cds_candidates_dictionary = {
        'candidates_cds_coordinates': query_coordinates
    }
    overlapping_targets = exonize_obj.data_container.get_overlapping_clusters(
        target_coordinates_set=target_coordinates,
        threshold=exonize_obj.environment.targets_clustering_overlap_threshold
    )
    expected_output = {
        P.open(200, 250): {
            'reference': P.open(200, 250),
            'mode': exonize_obj.environment.full
        },
        P.open(210, 250): {
            'reference': P.open(200, 250),
            'mode': exonize_obj.environment.full
        },
        P.open(220, 250): {
            'reference': P.open(220, 250),
            'mode': exonize_obj.environment.partial_insertion
        },
        P.open(0, 50): {
            'reference': P.open(0, 50),
            'mode': exonize_obj.environment.partial_insertion
        },
        P.open(0, 48): {
            'reference': P.open(0, 50),
            'mode': exonize_obj.environment.partial_insertion
        },
        P.open(40, 90): {
            'reference': P.open(40, 90),
            'mode': exonize_obj.environment.partial_insertion
        },
        P.open(220, 270): {
            'reference': P.open(220, 270),
            'mode': exonize_obj.environment.inter_boundary
        },
        P.open(215, 270): {
            'reference': P.open(220, 270),
            'mode': exonize_obj.environment.inter_boundary
        },
        P.open(219, 270): {
            'reference': P.open(220, 270),
            'mode': exonize_obj.environment.inter_boundary
        },
        P.open(400, 450): {
            'reference': P.open(400, 450),
            'mode': exonize_obj.environment.intronic
        },
        P.open(402, 450): {
            'reference': P.open(400, 450),
            'mode': exonize_obj.environment.intronic
        },
        P.open(420, 450): {
            'reference': P.open(420, 450),
            'mode': exonize_obj.environment.intronic
        }
        # Assuming no CDS overlap
    }
    assert exonize_obj.event_reconciler.get_matches_reference_mode_dictionary(
        cds_candidates_set=set(cds_candidates_dictionary['candidates_cds_coordinates']),
        clusters_list=overlapping_targets,
        gene_cds_set=set()
    ) == expected_output


def test_get_matches_reference_mode_dictionary(
        exonize_obj
):
    gene_cds_set = {
        P.open(0, 100),
        P.open(200, 300),
        P.open(300, 310),
        P.open(320, 380),
        P.open(500, 600),
        P.open(800, 900)
    }
    cds_candidates_coordinates = {
        P.open(0, 100),
        P.open(200, 300),
        P.open(320, 380),
        P.open(500, 600),
        P.open(800, 900),
        P.open(1800, 2100),

    }
    target_coordinates = {
        (P.open(205, 300), 1e-3),
        (P.open(4, 99), 1e-3),
        (P.open(6, 101), 1e-3),
        (P.open(700, 800), 1e-4),
        (P.open(705, 796), 1e-3),
        (P.open(1900, 2000), 1e-3),
        (P.open(1910, 1997), 1e-3)
    }
    clusters = exonize_obj.data_container.get_overlapping_clusters(
        target_coordinates_set=target_coordinates,
        threshold=exonize_obj.environment.targets_clustering_overlap_threshold
    )
    targets_reference_coordinates_dictionary = exonize_obj.event_reconciler.get_matches_reference_mode_dictionary(
        clusters_list=clusters,
        cds_candidates_set=cds_candidates_coordinates,
        gene_cds_set=gene_cds_set
    )
    expected_res = {
        P.open(205, 300): dict(reference=P.open(200, 300), mode=exonize_obj.environment.full),
        P.open(4, 99): dict(reference=P.open(0, 100), mode=exonize_obj.environment.full),
        P.open(6, 101): dict(reference=P.open(0, 100), mode=exonize_obj.environment.full),
        P.open(700, 800): dict(reference=P.open(700, 800), mode=exonize_obj.environment.intronic),
        P.open(705, 796): dict(reference=P.open(700, 800), mode=exonize_obj.environment.intronic),
        P.open(1900, 2000): dict(reference=P.open(1900, 2000), mode=exonize_obj.environment.partial_insertion),
        P.open(1910, 1997): dict(reference=P.open(1900, 2000), mode=exonize_obj.environment.partial_insertion)
    }
    assert targets_reference_coordinates_dictionary == expected_res


def test_get_hits_query_and_target_coordinates(
        exonize_obj
):
    local_records_set = {
        (1, 'gene1', 0, 100, 205, 300, 1e-3),
        (2, 'gene1', 200, 300, 4, 99, 1e-3),
        (3, 'gene1', 800, 900, 6, 101, 1e-3),
        (4, 'gene1', 500, 600, 700, 800, 1e-3),
        (5, 'gene1', 800, 900, 700, 800, 1e-4)
    }
    query_coordinates = {
        P.open(0, 100),
        P.open(200, 300),
        P.open(800, 900),
        P.open(500, 600)
    }
    target_coordinates = {
        (P.open(205, 300), 1e-3),
        (P.open(4, 99), 1e-3),
        (P.open(6, 101), 1e-3),
        (P.open(700, 800), 1e-4)
    }
    res_query_coords, res_target_coords = exonize_obj.event_reconciler.get_hits_query_and_target_coordinates(
        local_records_set=local_records_set
    )
    assert res_query_coords == query_coordinates
    assert target_coordinates == res_target_coords


def test_create_events_multigraph(
        gene_graph,
):
    gene_graph, _ = gene_graph
    graph_nodes = {
        (0, 100),
        (200, 300),
        (700, 800),
        (800, 900),
        (500, 600),
        (1800, 2100),
        (2200, 2350)
    }
    assert set(gene_graph.nodes) == graph_nodes


def test_find_query_cds_global_matches(
        exonize_obj,
):
    cds_coordinate = P.open(1800, 2100)
    global_records_set = {
        (P.open(0, 100), P.open(200, 300)),
        (P.open(0, 100), P.open(800, 900)),
        (P.open(200, 300), P.open(800, 900)),
        (P.open(0, 100), P.open(500, 600)),
        (P.open(1800, 2100), P.open(2200, 2300)),
        (P.open(2500, 2600), P.open(1800, 2100))
    }
    res = [(P.open(1800, 2100), P.open(2200, 2300)),
           (P.open(2500, 2600), P.open(1800, 2100))]
    assert exonize_obj.event_reconciler.find_query_cds_global_matches(
        global_records_set_pairs_set=global_records_set,
        cds_coordinate=cds_coordinate,
    ) == res


def test_find_local_match_in_global_matches(
        exonize_obj,
):
    cds_coordinate = P.open(1800, 2100)
    global_coordinates = [
        (P.open(1800, 2100), P.open(2200, 2300)),
        (P.open(2500, 2600), P.open(1800, 2100))
    ]
    assert P.open(2200, 2300) == exonize_obj.event_reconciler.find_local_match_in_global_matches(
        global_candidates=global_coordinates,
        cds_coordinate=cds_coordinate,
        reference_coordinate=P.open(2200, 2300)
    )
    assert P.open(2200, 2300) == exonize_obj.event_reconciler.find_local_match_in_global_matches(
        global_candidates=global_coordinates,
        cds_coordinate=cds_coordinate,
        reference_coordinate=P.open(2210, 2305)
    )
    assert P.open(2200, 2300) == exonize_obj.event_reconciler.find_local_match_in_global_matches(
        global_candidates=global_coordinates,
        cds_coordinate=cds_coordinate,
        reference_coordinate=P.open(2190, 2305)
    )


def test_build_events_list(
        exonize_obj,
        gene_graph
):
    gene_graph, targets_reference_coordinates_dictionary = gene_graph
    mode_dictionary = exonize_obj.event_reconciler.build_mode_dictionary(
        targets_reference_coordinates_dictionary=targets_reference_coordinates_dictionary
    )
    event1_component, event2_component = sorted(list(nx.connected_components(gene_graph)), key=lambda x: len(x))
    event_coordinates_dictionary1 = exonize_obj.event_reconciler.build_event_coordinates_dictionary(
        component=event1_component,
        mode_dict=mode_dictionary,
        gene_graph=gene_graph,
    )
    exonize_obj.event_reconciler.assign_cluster_ids_to_event_coordinates(
        event_coordinates_dictionary=event_coordinates_dictionary1
    )
    event_coordinates_dictionary2 = exonize_obj.event_reconciler.build_event_coordinates_dictionary(
        component=event2_component,
        mode_dict=mode_dictionary,
        gene_graph=gene_graph,
    )
    exonize_obj.event_reconciler.assign_cluster_ids_to_event_coordinates(
        event_coordinates_dictionary=event_coordinates_dictionary2
    )
    res_events_list1 = {
        ('gene1', 'FULL', 1800, 2100, 1, None, 0),
        ('gene1', 'FULL', 2200, 2350, 1, None, 0)
    }
    res_events_list2 = {
        ('gene1', 'FULL', 0, 100, 3, None, 1),
        ('gene1', 'FULL', 200, 300, 2, None, 1),
        ('gene1', 'FULL', 500, 600, 2, None, 1),
        ('gene1', 'INTRONIC', 700, 800, 2, None, 1),
        ('gene1', 'FULL', 800, 900, 3, None, 1)
    }
    assert res_events_list1 == set(
        exonize_obj.event_reconciler.build_events_list(
            gene_id='gene1',
            event_coordinates_dictionary=event_coordinates_dictionary1,
            expansion_id_counter=0,
            gene_start=0
        )
    )
    assert res_events_list2 == set(
        exonize_obj.event_reconciler.build_events_list(
            gene_id='gene1',
            event_coordinates_dictionary=event_coordinates_dictionary2,
            expansion_id_counter=1,
            gene_start=0
        )
    )


def test_build_event_coordinates_dictionary(
        exonize_obj,
        gene_graph
):
    gene_graph, targets_reference_coordinates_dictionary = gene_graph
    events1 = {
        P.open(1800, 2100): ['FULL', 1, None],
        P.open(2200, 2350): ['FULL', 1, None]
    }
    events2 = {
        P.open(0, 100): ['FULL', 3, None],
        P.open(200, 300): ['FULL', 2, None],
        P.open(500, 600): ['FULL', 2, None],
        P.open(700, 800): ['INTRONIC', 2, None],
        P.open(800, 900): ['FULL', 3, None]
    }
    event1_component, event2_component = sorted(
        list(nx.connected_components(gene_graph)),
        key=lambda x: len(x)
    )
    mode_dictionary = exonize_obj.event_reconciler.build_mode_dictionary(
        targets_reference_coordinates_dictionary=targets_reference_coordinates_dictionary
    )
    assert events1 == exonize_obj.event_reconciler.build_event_coordinates_dictionary(
        component=event1_component,
        mode_dict=mode_dictionary,
        gene_graph=gene_graph
    )
    assert events2 == exonize_obj.event_reconciler.build_event_coordinates_dictionary(
        component=event2_component,
        mode_dict=mode_dictionary,
        gene_graph=gene_graph
    )


def test_gene_non_reciprocal_fragments(
        exonize_obj,
        gene_graph,
):
    #  local_records_set = {
    #     (1, 'gene1', 0, 100, 205, 300, 1e-3),  # full - there's a global match so we ignore
    #     (2, 'gene1', 200, 300, 4, 99, 1e-3),  # reciprocal full - there's a global match so we ignore
    #     (3, 'gene1', 800, 900, 6, 101, 1e-3),  # full - there's a global match so we ignore
    #     (4, 'gene1', 500, 600, 700, 800, 1e-3),  # intronic
    #     (5, 'gene1', 800, 900, 700, 800, 1e-4),  # intronic
    #     (6, 'gene1', 1800, 2100, 2210, 2300, 1e-4)  # full - there's a global match so we ignore
    # }
    # global_records_set = {
    #     (1, 0, 100, 200, 300),
    #     (3, 0, 100, 800, 900),
    #     (3, 200, 300, 800, 900),
    #     (3, 0, 100, 500, 600),
    #     (4, 1800, 2100, 2200, 2350)
    # }
    gene_graph, targets_reference_coordinates_dictionary = gene_graph
    mode_dictionary = exonize_obj.event_reconciler.build_mode_dictionary(
        targets_reference_coordinates_dictionary=targets_reference_coordinates_dictionary
    )
    event1_component, event2_component = sorted(list(nx.connected_components(gene_graph)), key=lambda x: len(x))
    event_coordinates_dictionary1 = exonize_obj.event_reconciler.build_event_coordinates_dictionary(
        component=event1_component,
        mode_dict=mode_dictionary,
        gene_graph=gene_graph,
    )
    exonize_obj.event_reconciler.assign_cluster_ids_to_event_coordinates(
        event_coordinates_dictionary=event_coordinates_dictionary1
    )
    event_coordinates_dictionary2 = exonize_obj.event_reconciler.build_event_coordinates_dictionary(
        component=event2_component,
        mode_dict=mode_dictionary,
        gene_graph=gene_graph,
    )
    exonize_obj.event_reconciler.assign_cluster_ids_to_event_coordinates(
        event_coordinates_dictionary=event_coordinates_dictionary2
    )
    res_events_list1 = set(
        exonize_obj.event_reconciler.build_events_list(
            gene_id='gene1',
            event_coordinates_dictionary=event_coordinates_dictionary1,
            expansion_id_counter=0,
            gene_start=0
        )
    )
    res_events_list2 = set(
        exonize_obj.event_reconciler.build_events_list(
            gene_id='gene1',
            event_coordinates_dictionary=event_coordinates_dictionary2,
            expansion_id_counter=1,
            gene_start=0
        )
    )
    res = [('INTRONIC', 5), ('INTRONIC', 4)]
    assert res == exonize_obj.event_reconciler.gene_non_reciprocal_fragments(
        gene_graph=gene_graph,
        events_list=[*res_events_list1, *res_events_list2],
        gene_start=0
    )


def test_get_full_event_components(
        exonize_obj,
        gene_graph
):
    gene_graph, targets_reference_coordinates_dictionary = gene_graph
    full_expansion_1 = {('gene1', 'FULL', 0, 100, 3, 1), ('gene1', 'FULL', 200, 300, 2, 1),
                        ('gene1', 'FULL', 500, 600, 1, 1), ('gene1', 'FULL', 800, 900, 2, 1)}
    full_expansion_0 = {('gene1', 'FULL', 1800, 2100, 1, 0), ('gene1', 'FULL', 2200, 2350, 1, 0)}
    event1_component, event2_component = sorted(list(nx.connected_components(gene_graph)), key=lambda x: len(x))
    full_expansions0 = exonize_obj.event_reconciler.get_full_event_component(
        component=event1_component,
        gene_graph=gene_graph,
        gene_start=0,
        gene_id='gene1',
        expansion_id_counter=0
    )
    assert set(full_expansions0) == full_expansion_0
    full_expansions1 = exonize_obj.event_reconciler.get_full_event_component(
        component=event2_component,
        gene_graph=gene_graph,
        gene_start=0,
        gene_id='gene1',
        expansion_id_counter=1
    )
    assert set(full_expansions1) == full_expansion_1

# def test_map_edges_to_records():
#     pass
#
# def test_is_tandem_pair():
#     pass
#
# def test_get_gene_full_events_tandemness_tuples():
#     pass
