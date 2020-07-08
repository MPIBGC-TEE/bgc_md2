# The Idea is to create an object representing the given information 
# that has methods depending on the computable properties 
# this will hopefully integrate with ipython and jupyter so that typing a 
# <instanceName><.><tab> shows the available methods

from bgc_md2.resolve.mvars import (
    InFluxesBySymbol
    ,OutFluxesBySymbol
    ,InternalFluxesBySymbol
    ,TimeSymbol
    ,StateVariableTuple
    ,CompartmentalMatrix
)


from bgc_md2.models.helpers import (
    computable_mvars
    ,get_single_mvar_value
)


# Idea is to create  the class directly by a call to type since type accepts 
# a dictionary for the class attributes (including) methods.

#ModelDescriptor=type(
#        'ModelDescriptor',
#        (),
#        {
#            'modelName' : 'testVectorFree'
#            ,
#            'getTimeSymbol'    : lambda obj :get_single_mvar_value(TimeSymbol,obj.modelName)
#        }
#)
#md=ModelDescriptor()
#
#md.getTimeSymbol()

# now we create the dictionary by a comprehension
def funcMaker(var,modelName):
    def meth(obj):
        return  get_single_mvar_value(var,modelName) 

    return meth

modelName = 'testVectorFree'
meth_dict={"get_"+var.__name__: funcMaker(var,modelName) for var in computable_mvars(modelName)}
ModelDescriptor=type( 'ModelDescriptor', (),meth_dict)
md=ModelDescriptor()

md.get_TimeSymbol()
