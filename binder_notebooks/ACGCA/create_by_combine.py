# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.4
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# # ACGCA (C only)

# We can either construct ACGCA from scratch by stepwise combination of its leaves, fine roots, coarse roots + branches, and trunk parts.

# +
import matplotlib.pyplot as plt
from sympy import Matrix

from CompartmentalSystems.smooth_reservoir_model import SmoothReservoirModel
from CompartmentalSystems import helpers_reservoir as hr
# -

from bgc_md2.models.ACGCA.leaves.source import mvs as leaves_mvs
from bgc_md2.models.ACGCA.roots.source import mvs as roots_mvs

leaves_srm = leaves_mvs.get_SmoothReservoirModel()
roots_srm = roots_mvs.get_SmoothReservoirModel()

# +
fig, axes = plt.subplots(ncols=2, figsize=(18, 12))

axes[0].set_title("Leaves", fontsize=24)
axes[1].set_title("Roots", fontsize=24)
leaves_srm.plot_pools_and_fluxes(axes[0], legend=False)
roots_srm.plot_pools_and_fluxes(axes[1], legend=False)

# +
from bgc_md2.models.ACGCA.__init__ import GPP, E, B_L, C_L, C_R, B_R, ML, GL, MR, GR

t = leaves_mvs.get_TimeSymbol()

leaves = (
    set(leaves_mvs.get_StateVariableTuple()),
    leaves_mvs.get_InFluxesBySymbol(),
    leaves_mvs.get_OutFluxesBySymbol(),
    leaves_mvs.get_InternalFluxesBySymbol()
)

roots = (
    set(roots_mvs.get_StateVariableTuple()),
    roots_mvs.get_InFluxesBySymbol(),
    roots_mvs.get_OutFluxesBySymbol(),
    roots_mvs.get_InternalFluxesBySymbol()
)

intersect = (
    {E: GPP},
    {E: ML + GL + MR + GR}
)
LR = hr.combine(leaves, roots, {}, {}, intersect)

LR_srm = SmoothReservoirModel.from_state_variable_indexed_fluxes(
    Matrix([E, B_L, C_L, C_R, B_R]),
    t,
    LR[1],
    LR[2],
    LR[3]
)

fig, ax = plt.subplots(figsize=(12, 12))
ax.set_title("Leaves and roots", fontsize=24)
LR_srm.plot_pools_and_fluxes(ax, legend=False)
# -

from bgc_md2.models.ACGCA.other.source import mvs as other_mvs
from bgc_md2.models.ACGCA.trunk.source import mvs as trunk_mvs

other_srm = other_mvs.get_SmoothReservoirModel()
trunk_srm = trunk_mvs.get_SmoothReservoirModel()

# +
fig, axes = plt.subplots(ncols=2, figsize=(18, 12))

axes[0].set_title("Other (coarse roots + branches)", fontsize=24)
axes[1].set_title("Trunk", fontsize=24)
other_srm.plot_pools_and_fluxes(axes[0], legend=False)
trunk_srm.plot_pools_and_fluxes(axes[1], legend=False)

# +
# this cell does not work because of a problem with the __init__.py 
# probably only some variable names changed

#from bgc_md2.models.ACGCA.__init__ import C_S, MS, GS_O, GS_T, B_OS, B_OH, B_TH, B_TS
#
#other = (
#    set(other_mvs.get_StateVariableTuple()),
#    other_mvs.get_InFluxesBySymbol(),
#    other_mvs.get_OutFluxesBySymbol(),
#    other_mvs.get_InternalFluxesBySymbol()
#)
#
#trunk = (
#    set(trunk_mvs.get_StateVariableTuple()),
#    trunk_mvs.get_InFluxesBySymbol(),
#    trunk_mvs.get_OutFluxesBySymbol(),
#    trunk_mvs.get_InternalFluxesBySymbol()
#)
#
#intersect = (
#    {E: GPP, C_S: 0},
#    {E: MS + GS_O + GS_T, C_S: 0}
#)
#
#OT = hr.combine(other, trunk, {}, {}, intersect)
#
#OT_srm = SmoothReservoirModel.from_state_variable_indexed_fluxes(
#    Matrix([E, B_OS, B_OH, C_S, B_TH, B_TS]),
#    t,
#    OT[1],
#    OT[2],
#    OT[3]
#)
#fig, ax = plt.subplots(figsize=(12, 12))
#ax.set_title("Other and trunk", fontsize=24)
#OT_srm.plot_pools_and_fluxes(ax, legend=False)

# +
#intersect = (
#    {E: GPP, C_S: 0},
#    {E: ML + GL + MR + GR + MS + GS_O + GS_T, C_S: 0}
#)
#
#LROT = hr.combine(LR, OT, {}, {}, intersect)
#
#LROT_srm = SmoothReservoirModel.from_state_variable_indexed_fluxes(
#    Matrix([E, B_L, C_L, B_OS, B_OH, C_S, B_TH, B_TS, C_R, B_R]),
#    t,
#    LROT[1],
#    LROT[2],
#    LROT[3]
#)
#fig, ax = plt.subplots(figsize=(18, 18))
#ax.set_title("Complete C model (constructed)", fontsize=24)
#LROT_srm.plot_pools_and_fluxes(ax, legend=False)
# -

# Or we load it directly from the database.

# +
from bgc_md2.models.ACGCA.source import mvs as mvs
srm = mvs.get_SmoothReservoirModel()

fig, ax = plt.subplots(figsize=(18, 18))
ax.set_title("Complete C model (from database)", fontsize=24)
srm.plot_pools_and_fluxes(ax, legend=False)
# -


