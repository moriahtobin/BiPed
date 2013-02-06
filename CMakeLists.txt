i3_project(biped)

i3_add_library(biped
	private/biped/BipedLikelihood.cxx
	# private/biped/BipedParametarization.cxx
	USE_TOOLS log4cplus boost python suitesparse
	USE_PROJECTS icetray dataclasses millipede lilliput gulliver 
	    photonics-service phys-services tableio
)