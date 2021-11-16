from unittest import TestCase, skip
from testinfrastructure.InDirTest import InDirTest
from general_helpers import (
        day_2_month_index, 
        month_2_day_index,
        months_by_day_arr,
        TimeStepIterator2,
        respiration_from_compartmental_matrix
)
class Test_general_helpers(InDirTest):
    def test_month_2_day_index(self):
        self.assertEqual(
                month_2_day_index([0]),
                [0]
        ) 
        self.assertEqual(
                month_2_day_index([1]),
                [31]
        ) 
        self.assertEqual(
                month_2_day_index([2]),
                [59]
        ) 
        self.assertEqual(
                month_2_day_index([3]),
                [90]
        ) 
        self.assertEqual(
                month_2_day_index([1,3]),
                [31,90]
        ) 
    
    def test_day_2_month_index(self):
        # note that days are counted from zero so day 30 is January 31.
        self.assertEqual(day_2_month_index( 0), 0) 
        self.assertEqual(day_2_month_index(30), 0) 
        self.assertEqual(day_2_month_index(31), 1) 
        self.assertEqual(day_2_month_index(60), 2) 
