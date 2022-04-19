import sys
sys.path.insert(0, "..")  # necessary to import 
from TestSymbolic import TestSymbolic
import os
from pathlib import Path
from unittest import TestSuite, TextTestRunner
folder_name=Path(os.getcwd()).name
print(folder_name)

class A(TestSymbolic):
    model_folders = [folder_name]

method_list = [method for method in dir(A) if method.startswith("test")]
s = TestSuite()
for name in method_list:
    s.addTest(A(name))
runner = TextTestRunner()
os.chdir(Path('..'))
runner.run(s)
os.chdir(folder_name)

