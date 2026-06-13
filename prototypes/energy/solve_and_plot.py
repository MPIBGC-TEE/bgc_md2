import helpers as h

for mod_name in [
        "P2TherMic_dormancy_2_Monod",
        "P2TherMic_dormancy_2_MTS",
        "P2TherMic_dormancy_3_Monod",
        "P2TherMic_dormancy_3_MTS"
    ]:    
    h.solve_and_plot(mod_name)

