
from IPython.display import display_pretty, display
from sympy import Matrix,Function,Symbol,diff, Rational
from sympy.printing import latex
from IPython.display import Latex,Math
x_1,x_2=map(Symbol,("x_1","x_2"))
M=Matrix([[1,-1/2],[-1/2,1/2]])
M=Matrix([[5,0],[-Rational(9,2),1]])
I=Matrix([1,0])
sv=Matrix([x_1,x_2])

X_c=M.inverse_LU()*I

# show the sign flip  in one of the components
X_0=Matrix([1,1])
X_p_0=X_c-X_0
X_dot_0=I-M*X_0
X_c

display(X_c, X_dot_0,X_p_0, sum(X_dot_0), sum(X_p_0))

eq =  r'\frac{d}{dt}'+ rf'{latex(sv)}'+ "="+ rf'{latex(I)}'+ r'+'+ rf'{latex(-M)}'+ rf'{latex(sv)}'
print(eq)
print("##############")

text = rf"""
$\tau={latex(M.inverse_LU())}$,
$\X_c={latex(X_c)}$,
$\X_0={latex(X_0)}$,
$\X_p={latex(X_p_0)}$,
$\X^\prime={latex(X_dot_0)}$
"""
print(text)

# now show the sign flip in the sum (just for Yiqi)

I=Matrix([1,1])
X_c=M.inverse_LU()*I
X_0=Matrix([Rational(1,3),Rational(9,5)])
X_p_0=X_c-X_0
X_dot_0=I-M*X_0
display(X_c, X_dot_0,X_p_0, sum(X_dot_0), sum(X_p_0))
text = rf"""
$I={latex(I)}$,
$\X_c={latex(X_c)}$,
$\X_0={latex(X_0)}$,
$\X_p={latex(X_p_0)}$,
$\X^\prime={latex(X_dot_0)}$
"""
print(text)

eq =  r'\operatorname{sign}(\sum_i \mathbf{X_p})= \operatorname{sign}('+\
rf'{latex(sum(X_p_0))})'+\
r' \ne \operatorname{sign}('+rf'{latex(sum(X_dot_0))})='+\
r'\operatorname{sign}(\sum_i (\mathbf{X}^{\prime})_i)'
print(eq)
