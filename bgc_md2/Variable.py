import re
import numpy as np
from cf_units import Unit


def FixDumbUnits(unit):
    r"""Try to fix the dumb units people insist on using.
    Parameters
    ----------
    unit : str
        the trial unit
    Returns
    -------
    unit : str
        the fixed unit
    """
    # Various synonyms for 1
    if unit.lower().strip() in ["unitless",
                                "n/a",
                                "none"]: unit = "1"
    # Remove the C which so often is used to mean carbon but actually means coulomb
    # also remove C14 which is even more ridiculous
    tokens = re.findall(r"[\w']+", unit)
    for i,token in enumerate(tokens):
        if token.endswith("C"):
            try:
                if Unit(token[:-1]).is_convertible(Unit("g")):
                    unit = unit.replace(token,token[:-1])
                elif (i > 0) and Unit(tokens[i-1]).is_convertible(Unit("g")):
                    unit = unit.replace(" C","")
            except:
                pass

        if token.endswith("C14"):
            try:
                if Unit(token[:-3]).is_convertible(Unit("g")):
                    unit = unit.replace(token,token[:-3])
                elif (i > 0) and Unit(tokens[i-1]).is_convertible(Unit("g")):
                    unit = unit.replace(" C14","")
            except:
                pass
    return unit


class Variable(object):

    def __init__(self, **keywords):
        self.name = keywords.get('name', None)
        self.data = keywords['data']
        self.unit = FixDumbUnits(keywords.get('unit', '1'))

        if not isinstance(self.data, np.ma.MaskedArray):
            self.data = np.ma.masked_array(
                data = self.data,
                mask = np.zeros_like(self.data)
            )

        if isinstance(self.data.mask, np.bool_):
            self.data.mask = np.ones_like(self.data.data) * self.data.mask 


    def __str__(self):
        s = str(type(self))
        s += '-' * len(str(type(self))) + '\n'
        s +='data: {}\n'.format(self.data)
        s += 'unit: {}\n'.format(self.unit)

        return s


    def data_mult(self, dim_var, given_axes):
        given_axes = np.array(given_axes).ravel().astype(int)

        dim_array = np.ones((1, self.data.ndim), int).ravel()
        dim_array[given_axes] = [self.data.shape[ax] for ax in given_axes]
        dim_var_reshaped = dim_var.data.reshape(dim_array)
        data_res = self.data*dim_var_reshaped
    
        unit_self = Unit(self.unit)
        unit_dim_var = Unit(dim_var.unit)
        unit0 = unit_self*unit_dim_var
        unit_res = Unit(unit0.format().split()[-1])
        unit0.convert(data_res, unit_res, inplace=True)

        data_res.set_fill_value(self.data.get_fill_value())
        return self.__class__(
            data = data_res,
            unit = '%s' % unit_res
        )


    def aggregateInTime(self, nstep):
        data = self.data.data[::nstep,...]
        mask = self.data.mask[::nstep,...]
        if (self.data.shape[0]-1) % nstep != 0:
            data = np.append(data, [self.data[-1,...]],axis=0)
            mask = np.append(mask, [self.data.mask[-1,...]],axis=0)

        return self.__class__(
            name = self.name,
            data = np.ma.masked_array(
                data       = data,
                mask       = mask,
                fill_value = self.data.get_fill_value()
            ),
            unit = self.unit
        )

    def convert(self, tar_unit):
        unit0 = Unit(self.unit)
        unit1 = Unit(tar_unit)
        #assert(unit0.is_convertible(unit1))
        #assert(unit0.is_convertible(unit1), "Units not compatible: %s and %s" % (unit0, unit1))
        data = self.data
        self.data = unit0.convert(data, unit1)
        self.unit = tar_unit

        return self


    def __add__(self, other):
        assert(self.unit == other.unit)

        data_res = self.data + other.data
        unit_res = self.unit

        return self.__class__(
            data = data_res,
            unit = unit_res
        )


    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return self.__add__(other)

    
    def absolute_error(self, v_right):
        assert(self.data.shape==v_right.data.shape)
        abs_err = self.__class__(
            data = np.abs(self.convert(v_right.unit).data-v_right.data),
            unit = v_right.unit
        )

        return abs_err


    def relative_error(self, v_right):
        assert(self.data.shape==v_right.data.shape)
        rel_err = self.__class__(
            data = 100 * self.absolute_error(v_right).data\
                     / np.abs(v_right.data),
            unit = '%'
        )

        return rel_err


class StockVariable(Variable):
    pass


class FluxVariable(Variable):

    def aggregateInTime(self, nstep):
        data_ind = np.arange(0, self.data.shape[0], nstep)
        data = np.add.reduceat(self.data.data, data_ind)
        mask = np.maximum.reduceat(self.data.mask, data_ind)

        masked_data = np.ma.masked_array(
            data       = data,
            mask       = mask,
            fill_value = self.data.get_fill_value()
        )
        return self.__class__(
            name = self.name,
            data = masked_data,
            unit = self.unit
        )


