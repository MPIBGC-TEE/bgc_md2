#!/usr/bin/env python3
import matplotlib.pyplot as plt
from bgc_md2.resolve.graph_helpers import ( 
    direct_prerequisites,
    sparse_powerset_graph
)
from bgc_md2.resolve.non_graph_helpers import  computable_mvars,directly_computable_mvars,input_mvars,output_mvar,all_mvars,applicable_computers

from unittest import TestCase,skip



class A:
    """A variable we assume to be given"""


class B:
    pass


class C:
    pass


class D:
    pass


class E:
    pass


class F:
    pass


class G:
    pass


class H:
    pass


class I:
    pass



def f_from_b(b: B) -> F:
    """Computes f from b"""
    return 4

def d_from_a_c(a:A,c:C)->D:
    """computes d from a  and c """
    returun (a+3)+c 

def d_from_b_c(b:B,c:C)->D:
    """computes d from a  and c """
    return b+c

def c_from_b(b:B)->C:
    """computes c from b"""
    return 2*b

def b_form_a(a:A)->B:
   """computes b from a"""
   return  a+3

def f_from_e(e:E)->F:
   """computes f from e"""
   return  e**2

class TestResolveWithoutGraphs(TestCase):

    def setUp(self):
        # we produce a small set of Mvars with a loop (b could be something
        # like a CompartmentalSystem that can be computed  in different ways
        # and can also be queried about its constituents.  
        # 
        self.mvars = {
            A,
            B,
            C,
            D,
            E,
            F,
            G,
            H,
            I,
        }
        self.computers = frozenset({
                f_from_b,
                d_from_a_c,
                d_from_b_c,
                c_from_b,
                b_form_a,
                f_from_e
        })
    def test_signature(self):
        self.assertEqual(
            input_mvars(f_from_b),
            frozenset({B})
        )
        self.assertEqual(
            input_mvars(d_from_a_c),
            frozenset({A,C})
        )

    def test_all_mvars(self):
        self.assertEqual(
            all_mvars(self.computers),
            frozenset({A, B, C, D, E, F })
         )

    def test_applicable_computers(self):
        self.assertEqual(
            applicable_computers(self.computers,frozenset({B,C})),
            frozenset({f_from_b,d_from_b_c,c_from_b})
         )

    def test_direct_computability(self):
        self.assertEqual(
            directly_computable_mvars(
                self.computers,
                frozenset({B,C})
            ),
            frozenset({F,D,C})
         )

    def test_computability(self):
        res=computable_mvars(
            allComputers=self.computers,
            available_mvars=frozenset([A,C]) 
        )
        #pe('mvars',locals())
        # e and f are not computable
        self.assertEqual(res,frozenset({
            A,
            B,
            C,
            D,
            F
         }))
