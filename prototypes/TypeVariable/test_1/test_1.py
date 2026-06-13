# Extension of test_0.py to an potentially unlimited number
# of subtypes without the necessity to mention the subtypes
# explicityly.
# The EPT type variable
# does not have to be redefined if a new Subtype is added.
# The subtype relationship is declared by the new type: EP_1
# and EP_2 both inherit from EP, EPT is the type that have
# to be a subtype of EP.  problem: Mixes class (inheritance)
# and subtyping relationship (the  class EP only exist to
# declare EPT as the Type of subtypes...)

from typing import TypeVar, Callable
import trendy9helpers.model_a.mod1 
import trendy9helpers.model_a.mod2 
import trendy9helpers.model_a.Obs#, model_b
from trendy9helpers.general_helpers import mcmc
#import sys
#sys.path.insert(0, "..")  # necessary to import general_helpers from subdirectories


res1 = mcmc(
    trendy9helpers.model_a.mod1.EstimatedParameters(a=5,b=5),
    trendy9helpers.model_a.Obs.Observations(x=1,y=2),
    trendy9helpers.model_a.mod1.param2res
)
print(res1)

res2 = mcmc(
    trendy9helpers.model_a.mod2.EstimatedParameters(a=3,b=2,c=5,d=3),
    trendy9helpers.model_a.Obs.Observations(x=1,y=2),
    trendy9helpers.model_a.mod2.param2res
)
print(res2)

# now we mix and EstimatedParameter instance with an incompatibel param2res
# but we do not get a runtime warning (since by malicious desing
# both EstimatedParameter classes 
# have "a" and "b" as fields , which could easily happen in practice
# If you run mypy on this file it will detect that mcmc has be called with 
# the wrong param2res function 
# 
res3 = mcmc(
    trendy9helpers.model_a.mod2.EstimatedParameters(a=3,b=2,c=5,d=3),
    trendy9helpers.model_a.Obs.Observations(x=1,y=2),
    trendy9helpers.model_a.mod1.param2res
)
print("This result is wrong. Install and run mypy on this file to see why! ", res3)


# lets check if we can use data installed with the package

trendy9helpers.model_a.mod1.show_file_content()
