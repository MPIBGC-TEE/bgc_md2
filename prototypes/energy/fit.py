import sys

#sys.path.insert(0,"DA")
sys.path.insert(0,".")
import helpers as h

for mod_name in [
        (
            #"DA_P2TherMic_dormancy_2_Monod"
            "DA2_P2TherMic_dormancy_2_Monod"
        ),
        #"P2TherMic_dormancy_2_MTS",
        #"P2TherMic_dormancy_3_Monod",
        #"P2TherMic_dormancy_3_MTS"
    ]:    
    h.fit(mod_name)
