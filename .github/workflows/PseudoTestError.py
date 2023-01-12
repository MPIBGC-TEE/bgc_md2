from unittest import TestCase

class PseudoTestFail(TestCase):
    def test_fail(self):
        raise Exception("borne to be wild")
