
from unittest import TestCase
from bgc_md2.resolve.graph_helpers import ( 
    #,
    powerlist
)
class TestSets(TestCase):
    
    def test_powerlist(self):
        self.assertEqual(
                 powerlist([])
                ,[[]]
        )
        self.assertEqual(
                 powerlist(['a'])      
                ,[ [] ,['a'] ]
        )
        self.assertEqual(
                 powerlist(['a','b'])      
                ,[[], ['a'], ['b'], ['a', 'b']]
        )
        self.assertEqual(
                 powerlist(['a','b','c'])      
                ,[[], ['a'], ['b'], ['a', 'b'], ['c'], ['a', 'c'], ['b', 'c'], ['a', 'b', 'c']]
        )
