#!/usr/bin/env python3
from string import Template
from pathlib import Path
import sys
import os
dir_path=Path(__file__).parent
file_name=Path(os.path.basename(__file__))
sys.path.insert(0,dir_path)
import shortesRelativePath as srp
from difflib import Differ

out_path=Path(".")
srp=srp.rp(s=out_path,t=dir_path)

t=Template(#"""
"""# automatically  procuded by ${fn} , do not change manually
${command} install -y -c conda-forge --file requirements.test --file requirements.doc --file requirements.non_src pip
pip install -r requirements.src
""")
for suffix in ["sh","bat"]:
    for command in ["conda","mamba","micromamba"]:
        txt=t.substitute(
            command=command if suffix =="sh" else f"call {command}",
            fn=srp.joinpath(file_name)
        )
        script_file_name=f"install_developer_{command}.{suffix}"
        with Path(script_file_name).open("w") as f:
            f.write(txt)
            #from IPython import embed; embed()
