

from bgc_md2.resolve.MVarSet import MVarSet
# If  
# just uncomment the one model you are working on and comment the others
model_names=[
    # "Pavlick2013Biogeosciences", # Ecosystem/Community model with equations and descriptions of fluxes and symbols. However, it can be ambiguous and some implementations are not fully compatible with the compartmental model
    # "Thomas2014GeosciModelDev", #Also has equations for N, but C fluxes also depend on them
    # "Luo2012TE",
    # "Wang2010Biogeosciences",
    # "Arora2005GCB-1",
    # "Comins1993EA",
    # "Running1988EcolModel",
    # "TECO", #Same as Luo2012TE (comes from the same paper)
    # "Emanuel1981",
    # "Rasmussen2016JMB",
    ##################################################################
    ####### SOIL MODELS: 
    ## #  *: Models with compartmental matrix = T*N 
    # "Fontaine2005Ecologyletters_4_2",
    # "Fontaine2005Ecologyletters_4_1",
    # "Fontaine2005Ecologyletters_3_2",
    # "Fontaine2005Ecologyletters_3_1",
    # "Fontaine2005Ecologyletters_2",
    # "Fontaine2005Ecologyletters_1",
    # "Fontaine2005Ecologyletters",
    # "sixPairsModel", #  *
    # "Wang2014BG3p",
    # "Wang2014BG2p",
    # "Wang2013EcologicalApplications", #  *
    # "Andren1997EA", #  *
    # "Jenkinson1977SoilScience", #  *
    # "Zelenev2000MicrobialEcology", #  *
    # "Allison2010NG", #  *
    # "Parton1987SoilSciSocAmJ", #  *
    # "Henin1945AA",
    ##################################################################
    ####### VEGETATION MODEL WITHOUT VegetationCarbonInputPartitioningTuple & VegetationCarbonInputScalar
    ######################################################################
    ######################################################################
    # "Hilbert1991AnnBot",
    ######################################################################
    ####### VEGETATION MODELS
    ######################################################################
    # "ElMasri2013AgricForMeteorol", #Paper shows results for soil carbon, and fig. 1 has litter and C pools, but the equations on table A2 (although very detailed) don't include them. There are also equations for Phenology, but they were not included in the source.py
    # "Scheiter2009GlobalChangeBiol", #No soil compartments
    # "Turgman2018EcologyLetters", #No soil compartments
    # "Haverd2016Biogeosciences", #No soil compartments
    # "Foley1996GBC", #No equations for litter and soil, but the figure has those compartments
    # "Gu2010EcologicalComplexity", #No equations for litter and soil, but the model description (CEVSA) mentions them
    # "King1993TreePhysiol", #No soil compartments
    # "DeAngelis2012TheorEcol", #No soil compartments (model based on G’Day, but removed litter and soil compartments)
    # "Potter1993GlobalBiogeochemicalCycles", #No equations for litter and soil, but the model description mentions them
    # "testVectorFree",
    # "Williams2005GCB",
    "CARDAMOM",
]
for mn in model_names:
    mvs = MVarSet.from_model_name(mn)
    mvars = mvs.computable_mvar_types()
    list_str = "\n".join(["<li> " + str(var.__name__) + " </li>" for var in mvars])
    print(list_str)
    for var in mvars:
        print("########################################")
        print(str(var.__name__))
        print(mvs._get_single_mvar_value(var))
