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
        
