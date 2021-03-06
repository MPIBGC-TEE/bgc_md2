from sympy import symbols, Matrix


# time symbol 
t = symbols("t")

# 1 g_gluc/g_dw, g_gluc_to_gC, g_dw_to_gC
zeta, zeta_gluc, zeta_dw = symbols("zeta zeta_gluc zeta_dw")

# active carbon pool and GPP flux
E = symbols("E")
GPP = symbols("GPP")

# state variables for leaves, roots, other, and trunk
B_L, C_L = symbols("B_L C_L")
C_R, B_R = symbols("C_R B_R")
C_S = symbols("C_S")
B_OS, B_OH = symbols("B_OS B_OH")
B_TS, B_TH = symbols("B_TS B_TH")

# costs of growth for leaves, roots, wood
# maintenance respration and senescence rates
C_gL, delta_L, R_mL, S_L = symbols("C_gL delta_L R_mL S_L")
C_gR, delta_R, R_mR, S_R = symbols("C_gR delta_R R_mR S_R")
C_gW, delta_W, delta_S, R_mS, S_O = symbols("C_gW delta_W delta_S R_mS S_O")

# dynamic carbon allocation coefficients to leaves, coarse roots,
# fine roots + branches, and trunk
f_L = symbols("f_L")
f_R = symbols("f_R")
f_O = symbols("f_O")
f_T = symbols("f_T")

# CUEs in production of leaves, roots, wood
eta_L = zeta_dw / (zeta_gluc * C_gL)
eta_R = zeta_dw / (zeta_gluc * C_gR)
eta_W = zeta_dw / (zeta_gluc * C_gW)

# maintenance respiration fluxes for leaves, roots, sapwood
ML = zeta_gluc/zeta_dw * R_mL * B_L
MR = zeta_gluc/zeta_dw * R_mR * B_R
B_S_star = symbols("B_S_star")
MS = R_mS * B_S_star

# growth respiration fluxes for leaves, roots, sapwood
GL = f_L * C_gL/(C_gL+delta_L) * (1-eta_L) * E
GR = f_R * C_gR/(C_gR+delta_R) * (1-eta_R) * E
GS = (f_O+f_T) * C_gW/(C_gW+delta_W) * (1-eta_W) * E

# growth respiration fluxes for leaves, roots, sapwood
GS_O = f_O * C_gW/(C_gW+delta_W) * (1-eta_W) * E
GS_T = f_T * C_gW/(C_gW+delta_W) * (1-eta_W) * E
    
# heartwood production rates from sapwood
v_O, v_T = symbols("v_O v_T") 




