from pathlib import Path
from testinfrastructure.test_notebooks import test_notebooks

file_names = [
    "model_comparison_presentation.py",
    "test.py",
    "test2.py",
    "test3.py",
]
test_notebooks([Path(fn) for fn in file_names])
