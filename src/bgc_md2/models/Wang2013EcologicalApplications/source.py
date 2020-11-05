from sympy import var, ImmutableMatrix
from frozendict import frozendict
from bgc_md2.resolve.mvars import (
    CompartmentalMatrix,
    InputTuple,
    TimeSymbol,
    StateVariableTuple,
)
from ..BibInfo import BibInfo 
#from bgc_md2.resolve.MVarSet import MVarSet
from bgc_md2.helper import MVarSet

sym_dict = {
        'P': 'particulate organic carbon pool' # "mgC g^{-1}\\text{ soil}"
        ,'M': 'mineral associated organic carbon pool' # "mgC g^{-1}\\text{ soil}"
        ,'Q': '"active layer of $M$ interacting with dissolved organic carbon through adsorption and desorption"' # "mgC g^{-1}\\text{ soil}"
        ,'B': 'microbial biomass pool' # "mgC g^{-1}\\text{ soil}"
        ,'D': 'dissolved organic carbon pool' # "mgC g^{-1}\\text{ soil}"
        ,'EP': '"enzymes for the decomposition of $P$"' # "mgC g^{-1}\\text{ soil}"
        ,'EM': '"enzymes for the decompsotion of $M$"' # "mgC g^{-1}\\text{ soil}"
        ,'V_P': '"maximum specific decomposition rate for $P$ by $EP$"' # "h^{-1}"
        ,'K_P': '"half-saturation constant for decomposition of $P$"' # "mgC g^{-1}\\text{ soil}"
        ,'V_M': '"maximum specific decomposition rate for $M$ by $EM$"' # "h^{-1}"
        ,'K_M': '"half-saturation constant for decomposition of $M$"' # "mgC g^{-1}\\text{ soil}"
        ,'V_D': '"maximum specific uptake rate of $D$ for growth of $B$"' # "h^{-1}"
        ,'K_D': '"half-saturation constant for uptake of $D$ for growth of $B$"' # "mgC g^{-1}\\text{ soil}"
        ,'m_R': 'specific maintenance factor or rate' # "h^{-1}"
        ,'E_C': 'carbon use efficiency' #
        ,'Q_max': 'maximum dissolved organic carbon sorption capacity' # "mgC g^{-1}\\text{ soil}"
        ,'K_ads': 'specific adsorption rate' # "h^{-1}"
        ,'K_des': 'desorption rate' # "h^{-1}"
        ,'r_EP': '"turnover rate of $EP$"' # "h^{-1}"
        ,'r_EM': '"turnover rate of $EM$"' # "h^{-1}"
        ,'g_D': '"fraction of dead $B$ allocated to $D$"' # 
        ,'p_EP': '"fraction of $m_R$ for production of $EP$"' #
        ,'p_EM': '"fraction of $m_R$ for production of $EM$"' #
        ,'f_D': '"fraction of decomposed $P$ allocated to $D$"' # 
        ,'F_E': 'enzyme production rate'
        ,'F_R': 'microbial respiration rate'
        ,'F_U': 'dissolved organic matter uptakte rate'
        ,'F_A': 'adsorption rate of dissolved organic matter'
        ,'I_P': 'soil organic carbon input rate' # "mgC g^{-1}\\text{ soil } h^{-1}"
        ,'I_D': 'dissolved organic carbon input rate' # "mgC g^{-1}\\text{ soil } h^{-1}"
}

for name in sym_dict.keys():
    var(name)
F_E = m_R # "h^{-1}"   
F_R = (1/E_C-1)*(V_D+m_R)*D/(K_D+D) # "h^{-1}"
F_U = 1/E_C*(V_D+m_R)*B/(K_D+D) # "h^{-1}"
F_A = K_ads*(1-Q/Q_max) # "h^{-1}"
t = TimeSymbol("t") # unit: "hour"
x = StateVariableTuple((P, M, Q, B, D, EP, EM))
u = InputTuple((I_P, 0, 0, 0, I_D, 0, 0))
T = ImmutableMatrix([[   -1,  0,  0, (1-g_D)*(1-p_EP-p_EM)*F_E/(F_E+F_R),             0,  0,  0],
                     [1-f_D, -1,  0,                                   0,             0,  0,  0],
                     [    0,  0, -1,                                   0, F_A/(F_U+F_A),  0,  0],
                     [    0,  0,  0,                                  -1, F_U/(F_U+F_A),  0,  0],
                     [  f_D,  1,  1,     g_D*(1-p_EP-p_EM)*F_E/(F_E+F_R),            -1,  1,  1],
                     [    0,  0,  0,                  p_EP*F_E/(F_E+F_R),             0, -1,  0],
                     [    0,  0,  0,                  p_EM*F_E/(F_E+F_R),             0,  0, -1]])
N = ImmutableMatrix([[V_P*EP/(K_P+P),              0,           0,        0,       0,    0,    0],
                     [             0, V_M*EM/(K_M+M),           0,        0,       0,    0,    0],
                     [             0,              0, K_des/Q_max,        0,       0,    0,    0],
                     [             0,              0,           0,  F_E+F_R,       0,    0,    0],
                     [             0,              0,           0,        0, F_U+F_A,    0,    0],
                     [             0,              0,           0,        0,       0, r_EP,    0],
                     [             0,              0,           0,        0,       0,    0, r_EM]])
B = CompartmentalMatrix(T*N)

mvs = MVarSet({
    BibInfo(# Bibliographical Information
        name="MEND",
        longName="Microbial-Enzyme-Mediated Decomposition model", 
        version="1",
        entryAuthor="Holger Metzler",
        entryAuthorOrcid="0000-0002-8239-1601",
        entryCreationDate="21/03/2016",
        doi="10.1890/12-0681.1",
        further_references=BibInfo(doi="10.2307/1313568"),#Li2014Biogeochemistry
        sym_dict=sym_dict
    ),
    B,  # the overall compartmental matrix
    u,  # the overall input
    t,  # time for the complete system
    x,  # state vector of the complete system
})
