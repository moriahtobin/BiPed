i3_project(BiPed)

i3_add_library(BiPed
	private/LikBiped.cxx
	private/ParaBiped.cxx
	USE_TOOLS log4cplus boost python suitesparse
	USE_PROJECTS icetray dataclasses millipede lilliput gulliver 
	    photonics-service phys-services tableio
)
