import shortesRelativePath as srp
from  unittest import TestCase,skip
from pathlib import Path

class PathTest(TestCase):

    def test_same(self):
        s = Path("/home/mark.mueller/bgc_md2")
        t = Path("/home/mark.mueller/bgc_md2")
        self.assertEqual(srp.rp(s,t),Path("."))

    def test_split_rec(self):
        sp = Path("/home/mark.mueller/bgc_md2").parts
        tp = Path("/home/mark.mueller/").parts
        self.assertEqual(
            srp.split_3_rec((),sp, tp),
            (
                ('/', 'home', 'mark.mueller'),
                ('bgc_md2',),
                ()
            )
        )

    def test_split_3(self):
        sp = Path("/home/mark.mueller/bgc_md2").parts
        tp = Path("/home/mark.mueller/").parts
        self.assertEqual(
            srp.split_3(sp, tp),
            (
                ('/', 'home', 'mark.mueller'),
                ('bgc_md2',),
                ()
            )
        )
    
    def test_rp_t_up(self):
        s = Path("/home/mark.mueller/bgc_md2")
        t = Path("/home/mark.mueller/")
        self.assertEqual(srp.rp(s,t),Path(".."))
    
    def test_rp_t_down(self):
        s = Path("/home/mark.mueller/bgc_md2/")
        t = Path("/home/mark.mueller/bgc_md2/scripts")
        self.assertEqual(srp.rp(s,t),Path("scripts"))
    
    def test_rp_t_up_down(self):
        s = Path("/home/mark.mueller/bgc_md2/scripts")
        t = Path("/home/mark.mueller/bgc_md2/tests")
        self.assertEqual(srp.rp(s,t),Path("../tests"))

