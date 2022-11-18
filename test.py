from exon_duplication_search import *

def test_get_coding_exons_across_transcripts():
    eds_obj = ExonDupSearch(None)
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
