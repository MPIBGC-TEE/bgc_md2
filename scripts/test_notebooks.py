# notebook_runner.py
import re
import sys
import os

from glob import glob
from pathlib import Path
from functools import reduce
from nbconvert.preprocessors import ExecutePreprocessor
from jupytext import jupytext


def run_nbnode(nb, run_path):
    # from
    # https://www.blog.pythonlibrary.org/2018/10/16/testing-jupyter-notebooks/
    proc = ExecutePreprocessor(timeout=600, kernel_name='python3')
    proc.allow_errors = True

    proc.preprocess(nb, {'metadata': {'path': run_path}})

    errors = []
    for cell in nb.cells:
        if 'outputs' in cell:
            for output in cell['outputs']:
                if output.output_type == 'error':
                    errors.append(output)

    return nb, errors


def test_py_notebook(py_p):
    with py_p.open("r") as f:
        nb = jupytext.read(f)

    print("#########################################################")
    print(str(py_p))
    nb, errors = run_nbnode(nb, run_path=py_p.parent)
    return str(py_p), errors


def is_notebook(py_p):
    with py_p.open("r") as f:
        text = f.read()

    pattern = "formats: ipynb,py:light"
    return (re.search(pattern, text) is not None)


if __name__ == '__main__':
    target_dir=sys.argv[1]

    notebook_paths = [Path(s) for s in glob(f"{target_dir}**/*.py", recursive=True) if is_notebook(Path(s))]
    print([str(p) for p in notebook_paths])

    def add_errors(acc, p):
        tup = test_py_notebook(p)
        s, errors = tup
        return acc + [tup] if len(errors) > 0 else acc

    all_errors = reduce(add_errors, notebook_paths, [])

    le = len(all_errors)
    if le > 0:
        print(f"{le} errors")
        print(all_errors)
        sys.exit(1)
