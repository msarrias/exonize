from exon_duplication_search import *

eds_obj = ExonDupSearch('','','','')
### Unit testing: ExonDupSearch Class
def test_get_small_large_interv():
    test_pairs = [(P.open(2,20), P.open(2,20)), (P.open(5,30), P.open(5,6))]
    test_pairs_sol = [(P.open(2,20), P.open(2,20)), (P.open(5,6), P.open(5,30))]
    for i, j in zip(test_pairs,test_pairs_sol):
        assert j[0], j[1] == eds_obj.get_small_large_interv(i[0], i[1])
    
    
def test_get_coding_exons_across_transcripts():
    test1 = {i:{} for i in [P.open(1,300), P.open(2,50), P.open(280,400), P.open(500,550)]}
    test1_sol = {i:{} for i in [P.open(1,300), P.open(280,400), P.open(500,550)]}
    test2 = {i:{} for i in [P.open(1,250), P.open(76,260), P.open(2,50),  P.open(500,550)]}
    test2_sol = {i:{} for i in [P.open(1,250), P.open(500,550)]}
    test3 = {i:{} for i in [P.open(1,250), P.open(76,260),P.open(77,250), P.open(2,50),  P.open(500,550)]}
    test3_sol = {i:{} for i in [P.open(1,250), P.open(500,550)]}
    test4 = {i:{} for i in [P.open(260, 400),P.open(499, 550),P.open(450, 500), P.open(1,250)]}
    test4_sol = copy.deepcopy(test4)
    for test, solution in zip([test1,test2,test3,test4], [test1_sol,test2_sol,test3_sol,test4_sol]):
        assert solution == eds_obj.get_coding_exons_across_transcripts(eds_obj.sort_interval_dict(test))
        

def test_check_global_range():
    test1 = [P.open(1,5), P.open(6,20), P.open(21,40)]
    test2 = [P.open(1,6), P.open(6,20), P.open(21,40)]
    test3 = [P.open(1,6), P.open(6,20), P.open(21,40), (41,41)]
    eds_obj.check_global_range(test1) == None
    status = 0
    for l in [test2, test3]:
        try:
            eds_obj.check_global_range(l)
        except ValueError:
            status +=1
    assert status == 2


def test_trans_exons_and_introns_coords():
    test1 = {P.open(5,10):{'id': 'exon1', 'type': 'coding_exon'}}
    res_test1 = {P.open(0,4):{'id': 'intron_0', 'type': 'intron'}, 
                 P.open(5,10):{'id': 'exon1', 'type': 'coding_exon'},
                 P.open(11,30):{'id': 'intron_1', 'type': 'intron'}, }
    test1_gene_coord = P.open(0,30)
    my_res1 = eds_obj.trans_exons_and_introns_coords(test1, test1_gene_coord)
    assert res_test1 == my_res1
    
    test2 = {P.open(5,10):{'id': 'exon1', 'type': 'coding_exon'},
             P.open(7,30):{'id': 'exon2', 'type': 'coding_exon'}}
    res_test2 = {P.open(5,10):{'id': 'exon1', 'type': 'coding_exon'},
                 P.open(7,30):{'id': 'exon2', 'type': 'coding_exon'}}
    test2_gene_coord = P.open(5,30)
    my_res2 = eds_obj.trans_exons_and_introns_coords(test2, test2_gene_coord)
    assert res_test2 == my_res2
    
    test3 = {P.open(5,10):{'id': 'exon1', 'type': 'coding_exon'},
             P.open(15,30):{'id': 'exon2', 'type': 'coding_exon'}}
    res_test3 = {P.open(1,4):{'id': 'intron_0', 'type': 'intron'},
                 P.open(5,10):{'id': 'exon1', 'type': 'coding_exon'},
                 P.open(11,14):{'id': 'intron_1', 'type': 'intron'},
                 P.open(15,30):{'id': 'exon2', 'type': 'coding_exon'}}
    test3_gene_coord = P.open(1,30)
    my_res3 = eds_obj.trans_exons_and_introns_coords(test3, test3_gene_coord)
    assert res_test3 == my_res3

    test4 = {P.open(5,10):{'id': 'exon1', 'type': 'coding_exon'},
             P.open(7,30):{'id': 'exon2', 'type': 'coding_exon'}}
    res_test4 = {P.open(5,10):{'id': 'exon1', 'type': 'coding_exon'},
                 P.open(7,30):{'id': 'exon2', 'type': 'coding_exon'}}
    test4_gene_coord = P.open(5,31)
    my_res4 = eds_obj.trans_exons_and_introns_coords(test4, test4_gene_coord)
    assert res_test4 == my_res4
    
    
def test_get_coords_with_coding_exons()
    ## Case 1 Exon overlaps CDS
    # i == k
    for ce in ['exon', 'coding_exon']:
        test_1_1 = {P.open(1,200): {'type': ce}, P.open(1,300): {'type': 'CDS'}}
        test_1_1_res = {P.open(1,200): {'type': 'coding_exon'}, P.open(201,300): {'type': 'CDS'}}
        my_res_1_1 = get_coords_with_coding_exons_new(test_1_1)
        assert test_1_1_res == my_res_1_1
        test_1_2 = {P.open(1,200): {'type': 'CDS'}, P.open(1,300): {'type': ce}}
        test_1_2_res = {P.open(1,200): {'type': 'coding_exon'}, P.open(201,300): {'type': 'exon'}}
        my_res_1_2 = get_coords_with_coding_exons_new(test_1_2)
        assert test_1_2_res == my_res_1_2 
        test_3_1 = {P.open(1,200): {'type': ce}, P.open(2,300): {'type': 'CDS'}}
        test_3_1_res = {P.open(1,200): {'type': 'coding_exon'}, P.open(201,300): {'type': 'CDS'}}
        my_res_3_1 = get_coords_with_coding_exons_new(test_3_1)
        assert test_3_1_res == my_res_3_1
        test_3_2 = {P.open(1,200): {'type': 'CDS'}, P.open(2,300): {'type': ce}}
        test_3_2_res = {P.open(1,200): {'type': 'coding_exon'}, P.open(201,300): {'type': 'exon'}}
        my_res_3_2 = get_coords_with_coding_exons_new(test_3_2)
        assert test_3_2_res == my_res_3_2
        # j == l
        test_1_3 = {P.open(1,300): {'type': ce}, P.open(200,300): {'type': 'CDS'}}
        test_1_3_res = {P.open(1,199): {'type': 'exon'}, P.open(200,300): {'type': 'coding_exon'}}
        my_res_1_3 = get_coords_with_coding_exons_new(test_1_3)
        assert test_1_3_res == my_res_1_3
        test_1_4 = {P.open(1,300): {'type': 'CDS'}, P.open(200,300): {'type': ce}}
        test_1_4_res = {P.open(1,199): {'type': 'CDS'}, P.open(200,300): {'type': 'coding_exon'}}
        my_res_1_4 = get_coords_with_coding_exons_new(test_1_4)
        assert test_1_4_res == my_res_1_4
        test_3_4 = {P.open(1,299): {'type': ce}, P.open(200,300): {'type': 'CDS'}}
        test_3_4_res = {P.open(1,199): {'type': 'exon'}, P.open(200,300): {'type': 'coding_exon'}}
        my_res_3_4 = get_coords_with_coding_exons_new(test_3_4)
        assert test_3_4_res == my_res_3_4
        test_3_5 = {P.open(1,300): {'type': ce}, P.open(200,299): {'type': 'CDS'}}
        test_3_5_res = {P.open(1,199): {'type': 'exon'}, P.open(200,300): {'type': 'coding_exon'}}
        my_res_3_5 = get_coords_with_coding_exons_new(test_3_5)
        assert test_3_5_res == my_res_3_5
        #k > i and l < j
        test_1_5 = {P.open(1,300): {'type': ce}, P.open(200,400): {'type': 'CDS'}}
        test_1_5_res = {P.open(1,199): {'type': 'exon'},
                        P.open(200,300): {'type': 'coding_exon'}, P.open(301,400): {'type': 'CDS'}}
        my_res_1_5 = get_coords_with_coding_exons_new(test_1_5)
        assert test_1_5_res == my_res_1_5
        test_1_6 = {P.open(1,300): {'type': 'CDS'}, P.open(200,400): {'type': ce}}
        test_1_6_res = {P.open(1,199): {'type': 'CDS'},
                        P.open(200,300): {'type': 'coding_exon'}, P.open(301,400): {'type': 'exon'}}
        my_res_1_6 = get_coords_with_coding_exons_new(test_1_6)
        assert test_1_6_res == my_res_1_6
        test_3_6 = {P.open(1,201): {'type': ce}, P.open(200,400): {'type': 'CDS'}}
        test_3_6_res = {P.open(1,201): {'type': 'exon'}, P.open(202,400): {'type': 'CDS'}}
        my_res_3_6 = get_coords_with_coding_exons_new(test_3_6)
        assert test_3_6_res == my_res_3_6
    ## Case 2 coding region overlaps non coding region and vice versa
    for coding_region in ['exon', 'coding_exon', 'CDS']:
        for utr_intron in ['intron'] + new_Pn_genome.UTR_features:
            # i == k
            test_2_1 = {P.open(1,200): {'type': coding_region}, P.open(1,300): {'type':utr_intron}}
            test_2_1_res = {P.open(1,300): {'type': utr_intron}}
            my_res_2_1 = get_coords_with_coding_exons_new(test_2_1)
            assert test_2_1_res == my_res_2_1
            test_2_2 = {P.open(1,200): {'type': utr_intron}, P.open(1,300): {'type':coding_region}}
            test_2_2_res = {P.open(1,200): {'type': utr_intron}, P.open(201,300) : {'type':coding_region}}
            my_res_2_2 = get_coords_with_coding_exons_new(test_2_2)
            assert test_2_2_res == my_res_2_2
            test_3_8 = {P.open(1,200): {'type': coding_region}, P.open(2,300): {'type':utr_intron}}
            test_3_8_res = {P.open(1,300): {'type': utr_intron}}
            my_res_3_8 = get_coords_with_coding_exons_new(test_3_8)
            assert test_3_8_res == my_res_3_8
            test_3_9 = {P.open(1,200): {'type': utr_intron}, P.open(2,300): {'type':coding_region}}
            test_3_9_res = {P.open(1,200): {'type': utr_intron}, P.open(201,300) : {'type':coding_region}}
            my_res_3_9 = get_coords_with_coding_exons_new(test_3_9)
            assert test_3_9_res == my_res_3_9
            # j == l
            test_2_3 = {P.open(1,300): {'type': coding_region}, P.open(200,300): {'type': utr_intron}}
            test_2_3_res = {P.open(1,199): {'type': coding_region}, P.open(200,300): {'type': utr_intron}}
            my_res_2_3 = get_coords_with_coding_exons_new(test_2_3)
            assert test_2_3_res == my_res_2_3
            test_2_4 = {P.open(1,300): {'type': utr_intron}, P.open(200,300): {'type': coding_region}}
            test_2_4_res = {P.open(1,300): {'type': utr_intron}}
            my_res_2_4 = get_coords_with_coding_exons_new(test_2_4)
            assert test_2_4_res == my_res_2_4
            test_3_10 = {P.open(1,299): {'type': coding_region}, P.open(200,300): {'type': utr_intron}}
            test_3_10_res = {P.open(1,199): {'type': coding_region}, P.open(200,300): {'type': utr_intron}}
            my_res_3_10 = get_coords_with_coding_exons_new(test_3_10)
            assert test_3_10_res == my_res_3_10
            test_3_11 = {P.open(1,299): {'type': utr_intron}, P.open(200,300): {'type': coding_region}}
            test_3_11_res = {P.open(1,300): {'type': utr_intron}}
            my_res_3_11 = get_coords_with_coding_exons_new(test_3_11)
            assert test_3_11_res == my_res_3_11
            # exon overlaps CDS : k > i and l < j
            test_2_5 = {P.open(1,300): {'type': coding_region}, P.open(200,400): {'type': utr_intron}}
            test_2_5_res = {P.open(1,199): {'type': coding_region}, P.open(200,400): {'type': utr_intron}}
            my_res_2_5 = get_coords_with_coding_exons_new(test_2_5)
            assert test_2_5_res == my_res_2_5
            test_2_6 = {P.open(1,300): {'type': utr_intron}, P.open(200,400): {'type': coding_region}}
            test_2_6_res = {P.open(1,300): {'type': utr_intron}, P.open(301,400): {'type': coding_region}}
            my_res_2_6 = get_coords_with_coding_exons_new(test_2_6)
            assert test_2_6_res == my_res_2_6
            test_3_12 = {P.open(1,201): {'type': coding_region}, P.open(200,400): {'type': utr_intron}}
            test_3_12_res = {P.open(1,201): {'type': coding_region}, P.open(202,400): {'type': utr_intron}}
            my_res_3_12 = get_coords_with_coding_exons_new(test_3_12)
            assert test_3_12_res == my_res_3_12
            test_3_13 = {P.open(1,201): {'type': utr_intron}, P.open(200,400): {'type': coding_region}}
            test_3_13_res = {P.open(1,201): {'type': utr_intron}, P.open(202,400): {'type': coding_region}}
            my_res_3_13 = get_coords_with_coding_exons_new(test_3_13)
            assert test_3_13_res == my_res_3_13

