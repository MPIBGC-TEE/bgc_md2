# If you want to run the script from this location (where the file lives)
# without installing the package you have to commit a certain institutianalized
# crime by adding the parent directory to the front of the python path.
# import sys
import shutil
from unittest.case import TestCase, skip
from pathlib import Path
from importlib import import_module
from importlib.resources import files as mod_files
from sympy import (
    Symbol,
    Function,
    sympify,
    simplify,
    lambdify,
    Function,
    Symbol,
    diff,
    exp,
    diag,
)
from plotly.offline import plot
import matplotlib.pyplot as plt
import numpy as np

from CompartmentalSystems.smooth_model_run import SmoothModelRun
from CompartmentalSystems.start_distributions import (
    start_age_moments_from_empty_spinup,
    start_age_moments_from_steady_state,
    start_age_moments_from_zero_initial_content,
    compute_fixedpoint_numerically,
    start_age_distributions_from_steady_state,
    start_age_distributions_from_empty_spinup,
    start_age_distributions_from_zero_initial_content,
)
from CompartmentalSystems.ArrayDictResult import ArrayDictResult
import CompartmentalSystems.helpers_reservoir as hr

from ComputabilityGraphs.CMTVS import CMTVS
from testinfrastructure.InDirTest import InDirTest

from bgc_md2.resolve.mvars import (
    NumericParameterization,
    NumericSimulationTimes,
    NumericStartValueArray,
    InFluxesBySymbol,
    OutFluxesBySymbol,
    InternalFluxesBySymbol,
    TimeSymbol,
    StateVariableTuple,
)
import bgc_md2.helper as h

from trendy9helpers import general_helpers as gh

model_names = {
    "ab_classic": "CLASSIC",
    "clm5": "CLM5.0",
    "kv_ft_dlem": "DLEM",
    "bian_ibis2": "IBIS",
    "cj_isam": "ISAM",
    "isba-ctrip": "ISBA-CTRIP",
    "jsbach": "JSBACH",
    "yz_jules": "JULES-ES-1p0",
    "lpj-guess": "LPJ-GUESS",
    "lpjwsl": "LPJ",
    "lpx-bern": "LPX-Bern",
    "ORCHIDEE-V2": "OCN",
    "ORCHIDEE": "ORCHIDEE",
    "ORCHIDEE-CNP": "ORCHIDEE-CNP",
    "ORCHIDEEv3": "ORCHIDEEv3",
    "Aneesh_SDGVM": "SDGVM",
    "kv_visit2": "VISIT",
    "jon_yib": "YIBs",
}
experiment_names = {k: v + "_S2_" for k, v in model_names.items()}


class TestsWithData(InDirTest):

    @property
    def model_folders(self):
        return [
            # first tier (best shape)
            "kv_visit2",
            "jon_yib",
            "yz_jules",
           # ##
            "Aneesh_SDGVM",  # second tier (not quite ready)
            # "kv_ft_dlem",
            ##
            ##third tier
            ##"cj_isam", # has problems with ODE solution probably due to wrong parameters
            ## msh.numericX0 also yields a negative pool value for the last pool
            # "bian_ibis2",#
            ##"cable-pop",
        ]

    def test_make_model_coord_transforms(self):
        for mf in set(self.model_folders):
            with self.subTest(mf=mf):
                th = gh.th(mf)
                lats,lons=th.lats_lons()

                msh = gh.msh(mf)
                # check that we predict the
                # right latitude values for
                # a given array index
                ctr = msh.make_model_coord_transforms()
                # print(lats)
                print([ctr.lat2LAT(lat) for lat in lats])
                print([ctr.lon2LON(lon) for lon in lons])

    def test_make_model_index_transforms(self):
        for mf in set(self.model_folders):
            with self.subTest(mf=mf):
                th = gh.th(mf)
                lats,lons=th.lats_lons()
                msh = gh.msh(mf)
                # check that we predict the
                # right latitude values for
                # a given array index
                tr = msh.make_model_index_transforms()
                n_lats = len(lats)
                print("n_lats={}".format(n_lats))
                for i in range(n_lats):
                    self.assertEqual(lats[i], tr.i2lat(i))

                n_lons = len(lons)
                print("n_lons={}".format(n_lons))
                for i in range(n_lons):
                    self.assertEqual(lons[i], tr.i2lon(i))

                ## inverse operation
                ## check that we can predict the index from a given
                ## latitude
                # for i in range(n_lats):
                #    self.assertEqual(tr.lat2i(lats[i]), i)
                ## or longitude
                # for i in range(n_lons):
                #    self.assertEqual(tr.lon2i(lons[i]), i)

                ## check the interpretation of the pixel boundaries
                # for i in range(n_lats - 1):
                #    # minimum belongs to pixel
                #    self.assertEqual(tr.lat2i(tr.i2lat_min_max(i)[0]), i)
                #    # maximum belongs to next pixel (if there is one)
                #    self.assertEqual(tr.lat2i(tr.i2lat_min_max(i)[1]), i + 1)

                # # the last lat pixel contains also the pole
                # last = n_lats - 1
                # lat_min, lat_max = tr.i2lat_min_max(last)
                # print(lat_max)
                # self.assertEqual(tr.lat2i(lat_max), last)


    def test_global_mean_cache(self):
        # make sure that the global means in the repo are
        # identical to the ones computed from the data (up to date)
        def equal(v1,v2):
            results=[
                np.allclose(
                    v1.__getattribute__(f),
                    v2.__getattribute__(f)
                )
                for f in v1._fields
            ] 
            return np.all(results)
            

        for mf in set(self.model_folders):
            with self.subTest(mf=mf):
                msh = gh.msh(mf)
                mdp = mod_files(f"trendy9helpers.{mf}")
                # where to look for global_means
                ref_cache_path = mdp.joinpath("global_means")
                test_cache_path = Path(mf).joinpath("global_means")
                conf_dict=gh.confDict(mf)
                svs, dvs = msh.get_global_mean_vars(
                    dataPath=Path(conf_dict['dataPath']), 
                    targetPath=test_cache_path, 
                    flash_cache=False
                )
                svs_ref, dvs_ref = msh.get_global_mean_vars(
                    dataPath=None,
                    targetPath=ref_cache_path, 
                    flash_cache=False
                )
                 
                for tup in [("svs", "svs_ref"),("dvs", "dvs_ref")]:
                    with self.subTest(tup=tup):
                        # from IPython import embed; embed()
                        t1,t2=tup
                        v1,v2=map(eval,tup)
                        res=equal(v1,v2)
                        if not res:
                            fig=plt.figure()
                            fig.suptitle(str(type(v1)).split(".")[-1])
                            n=len(v1._fields)
                            axs=fig.subplots(n,1)
                            for i,f in enumerate(v1._fields):
                                ax=axs[i]
                                ax.plot(
                                    v1.__getattribute__(f),
                                    label=t1
                                )
                                ax.plot(
                                    v2.__getattribute__(f),
                                    label=t2
                                )
                                ax.legend()
                                ax.set_title(f)
                            
                            fig.savefig(Path(mf).joinpath(f"compare_{tup[0]}.pdf"))
                        self.assertTrue(res)   

