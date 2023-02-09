from subprocess import run
from pathlib import Path
from collections import namedtuple
from sys import exit


Task = namedtuple(
    "Task",
    field_names=["args", "cwd"],
    defaults=[None]  # counting from right so referring to cwd
)


tasks = [
    Task(
        args=[
            "python",
            Path("scripts").joinpath("test_notebooks.py"),
            Path("tests").joinpath("notebooks"),
        ]
    ),
    Task(
        args=[
            "python",
            str(Path("scripts").joinpath("test_notebooks.py")),
            str(Path("binder_notebooks/")),
        ]
    ),
    Task(
        args=[
            "python", "-m", "unittest", "discover", "-t", ".", "-p", "Test*",
        ],
        cwd=Path("tests")
    ),
    Task(
        args=[
            "python", "-m", "unittest", "Test_general_helpers.py" 
        ],
        cwd=Path("prototypes").joinpath("working_group_2021")
    )
]
results = [run(**t._asdict()) for t in tasks]
print("summary")
for cp in results:
    print(cp)

errors = [cp for cp in results if cp.returncode != 0]
if len(errors)>0:
    for cp in errors:
        print("################")
        print(cp)
        print(cp.stderr)
    exit(1)

