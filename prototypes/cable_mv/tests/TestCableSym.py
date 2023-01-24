# import unittest
import matplotlib.pyplot as plt
import numpy as np
from string import ascii_lowercase, ascii_uppercase
from pathlib import Path
from sympy import Symbol, Function, latex
from frozendict import frozendict

from CompartmentalSystems.smooth_reservoir_model import SmoothReservoirModel
from CompartmentalSystems.smooth_model_run import SmoothModelRun

from testinfrastructure.InDirTest import InDirTest

from bgc_md2.resolve.helpers import (
    bgc_md2_computers,
    bgc_md2_computer_aliases,
    bgc_md2_mvar_aliases,
)
import bgc_md2.models.cable_all.cableCache as cC
import bgc_md2.models.cable_all.cablePaths as cP
import bgc_md2.models.cable_all.cableHelpers as cH

#from bgc_md2.models.helpers import (
#    provided_mvars,
#    computable_mvars,
#    path_dict_to_single_mvar,
#    get_single_mvar_value,
#    bgc_md2_computers,
#    bgc_md2_computer_aliases,
#    bgc_md2_mvar_aliases,
#)
from bgc_md2.resolve.graph_helpers import sparse_powerset_graph
from bgc_md2.resolve.non_graph_helpers import (
    all_mvars
)
from bgc_md2.resolve.mvars import (
    TimeSymbol,
    StateVariableTuple,
    CompartmentalMatrix,
    VegetationCarbonInputPartitioningTuple,
    VegetationCarbonInputTuple,
    VegetationCarbonInputScalar,
    VegetationCarbonStateVariableTuple,
    InputTuple,
    NumericParameterization,
    NumericStartValueDict,
    NumericSimulationTimes,
    NumericParameterizedSmoothReservoirModel,
)
from bgc_md2.models.BibInfo import BibInfo 

from bgc_md2.resolve.computers import numeric_parameterized_smooth_reservoir_model_1
from bgc_md2.resolve.graph_plotting import (
    AGraphComputerSetMultiDiGraph,
    AGraphComputerMultiDiGraph,
    draw_update_sequence,
    draw_ComputerSetMultiDiGraph_matplotlib,
    # ,draw_Graph_with_computers_svg
)
from testinfrastructure.helpers import pp
from bgc_md2.resolve.MVarSet import MVarSet


class TestCableSym(InDirTest):
    def setUp(self):
        self.mn = "cable_all"
        self.mvs = MVarSet.from_model_name(self.mn)

#    @unittest.skip
    def test_computable_mvars(self):
        #spsg=sparse_powerset_graph(bgc_md2_computers())
        mvs=self.mvs

        mvars = mvs.computable_mvar_types()
        list_str = "\n".join(["<li> " + str(var.__name__) + " </li>" for var in mvars])
        print(list_str)
        #for var in mvars:
        #    print("########################################")
        #    print(str(var.__name__))
        #    print(mvs._get_single_mvar_value(var))

    def test_numeric_input_tuple(self):
        # setup the paths to the testdata
        cable_out_path = Path('/home/data/cable-data/example_runs/parallel_1901_2004_with_spinup/output/new4')
        # cable_data_set = cH.cable_ds(cable_out_path)
        time_slice = slice(0, None)
        landpoint_slice = slice(1590, 1637) # cheated
        # landpoint_slice = slice(None,None)
        # time_slice=slice(None,None,None)

        zarr_cache_path = cP.slice_dir_path(
            cable_out_path,
            sub_dir_trunk="zarr_mm11",
            time_slice=time_slice,
            landpoint_slice=landpoint_slice,
        )
        if "cable_data_set" not in dir():
            cable_data_set = cH.cable_ds(cable_out_path)
        args = {
            "cable_data_set": cable_data_set,
            "zarr_cache_path": zarr_cache_path,
            "landpoint_slice": landpoint_slice,
            "time_slice": time_slice,
            #'batch_size': 128,
            "batch_size": 12,
            #'rm': True
        }

        x_org_iveg = cH.cacheWrapper(cC.x_org_iveg, **args)
        time  = cH.cacheWrapper(cC.time, **args)

        patches, landpoints = cC.all_pools_vary_cond_nz(**args)
        pcs = patches.compute()
        lpcs = landpoints.compute()
        pcs,lpcs
        p = Path('plots')
        p.mkdir(exist_ok=True)
        ind = 0 # first pair that has a nonconstant solution for 
        lp = lpcs[ind]
        patch = lpcs[ind]

        # get the symbolic represantation from the database
        mvs = self.mvs
        # define some stuff to extend it with

        def default(t):
            return 1

        leaf = Symbol('leaf')
        fine_root = Symbol('fine_root')
        Npp = Function("Npp")
        bvec_leaf = Function("bvec_leaf")
        bvec_fine_root = Function("bvec_fine_root")
        xk_leaf_cold = Function("xk_leaf_cold")
        xk_leaf_dry = Function("xk_leaf_dry")
        kleaf = Function("kleaf")
        kfroot = Function("kfroot")
        # bvec_wood = Function("bvec_wood")

        np1 = NumericParameterization(
            par_dict={},
            func_dict=frozendict(
                {
                    Npp: default,
                    bvec_fine_root: default,
                    bvec_leaf: default,
                    xk_leaf_cold: default,
                    kleaf: default,
                    kfroot: default,
                    xk_leaf_dry: default,
                }
            ),
        )
        nsv1 = NumericStartValueDict({leaf: 0.3, fine_root: 3.96})
        ntimes1 = NumericSimulationTimes(np.linspace(0, 1, 11))
        
        # extend the symbolice version with the new stuff
        pvs = mvs.provided_mvar_values
        #from IPython import embed; embed()
        pvs1= pvs.union(frozenset({np1, nsv1, ntimes1}))
        mvs1 = MVarSet(pvs1)

        x = mvs1.get_StateVariableTuple()
        Input = mvs1.get_InputTuple()
        #B = mvs1.get_CompartmentalMatrix()
        sym_times = mvs1.get_NumericSimulationTimes()
        sol_smooth = mvs1.get_NumericSolutionArray()

        comp_slice=slice(0,100)
        n = sol_smooth.shape[1]
        fig = plt.figure()
        for pool in range(n):
            ax = fig.add_subplot(n+1, 1, 2+pool)
            title ="\$" + latex(x[pool]) + "\$"
            #ax.plot(
            #    sym_times[comp_slice],
            #    sol_smooth[comp_slice, pool],
            #    color='r'
            #)
            ax.plot(
                time[comp_slice],
                x_org_iveg[comp_slice, pool, patch,lp],
                color='b'
            )
            fontsize = 10
            ax.set_title(title, fontsize=fontsize)
                 
        fig.savefig('solution.pdf')
