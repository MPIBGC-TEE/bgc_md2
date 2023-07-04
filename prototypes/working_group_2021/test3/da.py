# %load_ext autoreload
# %autoreload 2
from pathlib import Path

from trendy9helpers import general_helpers as gh

p=Path(".")
mf="kv_visit2"
da_scheme="da_1"
par_dir="par_2"
Cs,Js,epa_opt,cp=gh.gm_da_from_folder_names(p,mf,da_scheme,par_dir)

from IPython import embed; embed()
cp.write(Path(f"/home/mm/bgc_md2/src/bgc_md2/models/{mf}/parameterization_from_test_args/")) #uncomment if the parameterization should be available without the original driver data (which would be duplicated)

