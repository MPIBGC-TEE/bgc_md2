import numpy as np

from sympy import Matrix, symbols

from CompartmentalSystems.smooth_reservoir_model import SmoothReservoirModel

class ModelStructureException(Exception):
    def __str__(self): return "ModelStructureException"


class ModelStructure(object):
    def __init__(self, **keywords):
        self.pool_structure            = keywords['pool_structure']
        self.external_input_structure  = keywords['external_input_structure']
        self.horizontal_structure      = keywords['horizontal_structure']
        self.vertical_structure        = keywords.get('vertical_structure', dict())
        self.external_output_structure = keywords['external_output_structure']


        ## check that pools do not take mass from below/above and
        ## additionally get mass send from above/below

        if self.vertical_structure is not None:
            for pool_name, flux_structure in self.vertical_structure.items():
                to_below   = flux_structure.get('to_below'  , None)
                from_below = flux_structure.get('from_below', None)
                to_above   = flux_structure.get('to_above'  , None)
                from_above = flux_structure.get('from_above', None)
    
                if (to_below and to_above) or (from_below and from_above):
                    raise(ModelStructureException("Vertical fluxes overlap"))

        name2nr = dict()
        nr2name = dict()
        name2nr_lys = dict()
        pool_nr = 0
        for entry in self.pool_structure:
            pn = entry['pool_name']
            nr_lys = entry.get('nr_layers', 1)
            name2nr_lys[pn] = nr_lys
            for ly in range(nr_lys):
                name2nr[(pn, ly)] = pool_nr
                nr2name[pool_nr] = {'pool_name': pn, 'layer_nr': ly}
                pool_nr += 1

        self._name2nr = name2nr
        self._nr2name = nr2name
        self._name2nr_lys = name2nr_lys
        self.nr_pools = pool_nr


    def get_nr_layers(self, pool_name):
        return self._name2nr_lys[pool_name]


    def get_pool_nr(self, pool_name, layer_nr):
        return self._name2nr[(pool_name, layer_nr)]


    def get_pool_name_and_layer_nr(self, pool_nr):
        return self._nr2name[pool_nr]


    def get_pool_nrs(self, pool_name):
        pool_nrs = []
        nr_layers = self.get_nr_layers(pool_name)
        for ly in range(nr_layers):
            pool_nrs.append(self.get_pool_nr(pool_name, ly))

        return np.array(pool_nrs)


    def get_pool_nrs_set(self, pool_names, layers):
        pool_nrs = []
        for pool_name in pool_names:
            for ly in layers:
                pool_nrs.append(
                    self.get_pool_nr(
                        pool_name,
                        ly
                    )
                )

        return pool_nrs


    def get_nr_pools(self):
        nr_pools = 0
        for item in self.pool_structure:
            nr_layers = item.get('nr_layers', 1)
            nr_pools += nr_layers

        return nr_pools


    def get_flux_var_names(self):
        flux_list = []
        for d in [
            self.external_input_structure,
            self.horizontal_structure,
            self.vertical_structure,
            self.external_output_structure
        ]:
            for key, flux_sublist in d.items():
                flux_list.extend(flux_sublist)

        return flux_list

