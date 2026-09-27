import copy,json,sys,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'benchmarks'),str(ROOT/'src')]
from direct_nm_model import evaluator,subject_profile
from direct_missing_marker_ties import apply_missingness
from direct_marker_ties import predict,crossings
class SubjectModelTests(unittest.TestCase):
    def setUp(self):
        root=ROOT/'paper/calculator_missing_ties_20260916/H100'
        self.c,self.s,self.g=[json.loads((root/f).read_text()) for f in ('components.json','supplement.json','missing_components.json')]
    def test_anchor_preserved(self):
        c,_=apply_missingness(self.c,.001,self.g['gram_records'],self.g['decode_seconds_per_variant'])
        for fraction in (0.,.25):
            f,_=evaluator(self.c,self.s,self.g,8192,fixed_fraction=fraction)
            for m in (1024,65536,500000):self.assertEqual(f(m),predict(c,self.s,m))
    def test_work_and_immutability(self):
        original=copy.deepcopy((self.c,self.s,self.g))
        c,s,g=subject_profile(self.c,self.s,self.g,16384)
        self.assertEqual((self.c,self.s,self.g),original)
        self.assertEqual(s['binary_writer'],self.s['binary_writer'])
        self.assertEqual(s['pgen_bytes'],2*self.s['pgen_bytes'])
        self.assertEqual(g['gram_records'][0]['seconds_per_variant']['solve'],self.g['gram_records'][0]['seconds_per_variant']['solve'])
        self.assertEqual(g['gram_records'][0]['seconds_per_variant']['syrk'],2*self.g['gram_records'][0]['seconds_per_variant']['syrk'])
    def test_actual_sign_changes(self):
        products=[]
        for n in (2048,8192,32768):
            f,b=evaluator(self.c,self.s,self.g,n)
            self.assertAlmostEqual(b['complete'],.999**n)
            for tool in ('PLINK2','fastGWA-joint'):
                roots=crossings(f,tool,tolerance=1)
                self.assertTrue(roots)
                for root in roots:
                    lo,hi=root['lower_markers'],root['upper_markers']
                    self.assertLessEqual((f(lo)['torchGWAS']-f(lo)[tool])*(f(hi)['torchGWAS']-f(hi)[tool]),0)
                if tool=='PLINK2':products.append(n*roots[0]['upper_markers'])
        self.assertGreater(max(products)/min(products),1.05)
if __name__=='__main__':unittest.main()
