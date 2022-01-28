from IPython.display import Math
from sympy import latex
def mass_balance_equation(mvs):
    try:
        ip = mvs.get_InputTuple()
        cm = mvs.get_CompartmentalMatrix()
        sv = mvs.get_StateVariableTuple()
    
        eq = Math(
            r'\frac{dx}{dt}='+
            rf'{latex(ip)}'+
            r'+'+
            rf'{latex(cm)}'+
            rf'{latex(sv)}'
        )
        return eq
    except Exception():
        print(
            """The provided argument mvs does not allow all necessary computations to show the mass 
            balance equation:
             mvs.get_InputTuple()
             mvs.get_CompartmentalMatrix()
             mvs.get_StateVariableTuple()
            """
        )



