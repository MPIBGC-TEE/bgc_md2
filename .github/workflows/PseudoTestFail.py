from unittest import TestCase

class PseudoTestFail(TestCase):
    def test_fail(self):
        self.assertTrue(False)
