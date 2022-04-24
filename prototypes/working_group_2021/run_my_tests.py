import sys
sys.path.insert(0, "..")  # necessary to import 
from TestSymbolic import TestSymbolic
import os
from pathlib import Path
from unittest import TestSuite, TextTestRunner

folder_name=Path(os.getcwd()).name
print(folder_name)

args=sys.argv

single_test_name=str(args[1]) if len(args)>1 else None 


class A(TestSymbolic):
    @property
    def model_folders(self):
        return [folder_name]

method_list = [method for method in dir(A) if method.startswith("test")]
s = TestSuite()
if single_test_name is None:
    for name in method_list:
        s.addTest(A(name))
else:
    if single_test_name in method_list:
        s.addTest(A(single_test_name))
    else:
        print("no test with name {}".format(single_test_name))


runner = TextTestRunner()
os.chdir(Path('..'))
runner.run(s)
os.chdir(folder_name)

