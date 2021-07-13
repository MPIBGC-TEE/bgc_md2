import dask.array
from pathlib import Path
from unittest import skip
from testinfrastructure.InDirTest import InDirTest

from bgc_md2.models.cable_all.cableHelpers import (
        cacheWrapper,
        val_or_default,
        cable_ds
)
import bgc_md2.models.cable_all.cableCache as cC


# helpers
def time_dict(zarr_cache_path):
    subPaths = zarr_cache_path.iterdir()
    return {
        p.name: (p.stat()).st_mtime
        for p in subPaths
    }


class TestCableCache(InDirTest):
    def setUp(self):
        cable_out_path = Path('/home/data/cable-data/example_runs/parallel_1901_2004_with_spinup/output/new4')

        zarr_cache_path = Path("cache")
        self.cable_out_path = cable_out_path
        self.zarr_cache_path = zarr_cache_path
        self.time_slice = slice(0,100)
        self.landpoint_slice = slice(0,4)
        self.cable_data_set = cable_ds(cable_out_path)
        
        self.kwargs={
            'cable_data_set': self.cable_data_set,
            'zarr_cache_path': self.zarr_cache_path,
            'landpoint_slice': self.landpoint_slice,
            'time_slice': self.time_slice,
            'batch_size': 32
        }
        def check_presence(self,expected_names):
            zarr_dirs = [f.name for f in self.zarr_cache_path.iterdir()]
            self.assertEqual(
                frozenset(zarr_dirs),
                frozenset(expected_names)
            )
        def time_dict(self):
            return time_dict(self.zarr_cache_path)


    #@skip('takes long')
    def test_root(self):

        # compute variables directly from the original data
        funcs = [
            cC.time,
            cC.iveg,
            cC.Clitter,
            cC.Cplant,
            cC.Csoil,
            cC.NPP,
            cC.fracCalloc,
            cC.fromCWDtoS,
            cC.fromLeaftoL,
            cC.fromMettoS,
            cC.fromRoottoL,
            cC.fromSOMtoSOM,
            cC.fromStrtoS,
            cC.fromWoodtoL,
            cC.kplant,
            cC.xkNlimiting,
            cC.xktemp,
            cC.xkwater,
        ]
        results = map(lambda func: cacheWrapper( func, **self.kwargs), funcs)

        [r for r in results]
    
    def test_root0(self):
        funcs = [
            cC.Clitter0,
            cC.Cplant0,
            cC.Csoil0,
        ]
        results = map(lambda func: cacheWrapper( func, **self.kwargs), funcs)

        [r for r in results]
        
    def test_indirectly_dependent_B_org(self):

        res_B = cacheWrapper(
            cC.B_org,
            **self.kwargs
        )

        self.check_presence([
            "kplant",
            "fromLeaftoL",
            "fromRoottoL",
            "fromWoodtoL",
            "fromMettoS",
            "fromStrtoS",
            "fromCWDtoS",
            "fromSOMtoSOM",
            "xktemp",
            "xkwater",
            "xkNlimiting"
        ])


    def test_indirectly_dependent_B_val(self):

        func = cC.B_val
        res_B = cacheWrapper(
            func,
            **self.kwargs
        )
        self.check_presence([
            "iveg",  # for finding valid patch landpoint combies
            "Cplant",# for finding points with vegetation
            "kplant",
            "fromLeaftoL",
            "fromRoottoL",
            "fromWoodtoL",
            "fromMettoS",
            "fromStrtoS",
            "fromCWDtoS",
            "fromSOMtoSOM",
            "xktemp",
            "xkwater",
            "xkNlimiting"
        ])


    def test_indirectly_dependent_x_org(self):
        func = cC.x_org
        res_x = cacheWrapper(
            func,
            **self.kwargs
        )

        self.check_presence([
            'Clitter',
            'Cplant',
            'Csoil',
            'x_org'
        ])


    def test_indirectly_dependent_x0_org(self):
        func = cC.x0_org
        res_x = cacheWrapper(
            func,
            **self.kwargs
        )

        zarr_dirs = [f.name for f in zarr_cache_path.iterdir()]
        self.check_presence([
            'Clitter0',
            'Cplant0',
            'Csoil0'
            'x0_org'
        ])

    def test_indirectly_dependent_x0_val(self):
        func = cC.x0_val
        res_x = cacheWrapper(
            func,
            **self.kwargs
        )

        self.check_presence([
            'Clitter0',
            'Cplant0',
            'Csoil0',
            'x0_org',
            'nz',
            'x0_val',
        ])

    def test_indirectly_dependent_u_org(self):
        func = cC.u_org

        res_u = cacheWrapper(
            func,
            **self.kwargs
        )
        #self.assertTrue((res_D == 20 * raw_data).all().compute())

        self.check_presence([
            'iveg',
            'NPP',
            'fracCalloc',
            'u_org',
        ])

    def test_indirectly_dependent_u_val(self):
        func = cC.u_val

        res_u = cacheWrapper(
            func,
            **self.kwargs
        )
        #self.assertTrue((res_D == 20 * raw_data).all().compute())

        self.check_presence([
            'iveg',
            'NPP',
            'fracCalloc',
            'u_org',
            'nz',
            'u_val'
        ])

    def test_indirectly_dependent_sol_val(self):
        func = cC.sol_val

        res_sol_val = cacheWrapper(
            func,
            **self.kwargs
        )
        #self.assertTrue((res_D == 20 * raw_data).all().compute())

        self.check_presence([
            'NPP',
            'fracCalloc',
            'u_org',
            'u_val',
            'Clitter0',
            'Cplant0',
            'Csoil0',
            'x0_org',
            'nz',
            'x0_val',
            "iveg",  # for finding valid patch landpoint combies
            "Cplant",# for finding points with vegetation
            "kplant",
            "fromLeaftoL",
            "fromRoottoL",
            "fromWoodtoL",
            "fromMettoS",
            "fromStrtoS",
            "fromCWDtoS",
            "fromSOMtoSOM",
            "xktemp",
            "xkwater",
            "xkNlimiting"
        ])
