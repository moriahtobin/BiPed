#!/usr/bin/env python
#Hello Laura :)
from I3Tray import *
import sys
from icecube import icetray, dataio, dataclasses, photonics_service, gulliver
from icecube import millipede
load('gulliver-modules')
load('biped')

if len(sys.argv) < 3:
	print 'Usage: %s output.i3 input1.i3 [input2.i3] ...'
	sys.exit(1)

files = sys.argv[2:]

muon_service = photonics_service.I3PhotoSplineService('/home/nwhitehorn/photon-tables/emu_abs.fits', '/home/nwhitehorn/photon-tables/emu_prob.fits', 0)
cascade_service = photonics_service.I3PhotoSplineService('/home/nwhitehorn/photon-tables/ems_spice1_z20_a10.abs.fits', '/home/nwhitehorn/photon-tables/ems_spice1_z20_a10.prob.fits', 0)

tray = I3Tray()
tray.AddModule('I3Reader', 'reader', FilenameList=files)

tray.AddService('MuMillipedeParametrizationFactory', 'millipedeparam',
    MuonSpacing=1, ShowerSpacing=1, TimeStepsize=0, VertexStepSize=2,
    # DirectionStepsize=0.02, FitBiped=True)
    DirectionStepsize=0.02, LinLengthStepSize=1, FitBiped=True)
tray.AddService('BipedLikelihoodFactory', 'bipedllh',
    MuonPhotonicsService=muon_service, CascadePhotonicsService=cascade_service,
    PhotonsPerBin=1, Pulses='OfflinePulses')
tray.AddService('I3GSLRandomServiceFactory','I3RandomService')
tray.AddService('I3GulliverMinuit2Factory', 'minuit',
    MaxIterations=2000, Algorithm='SIMPLEX')
tray.AddService('I3BasicSeedServiceFactory', 'seed', FirstGuess='MPEFit',
    TimeShiftType='TNone')

tray.AddModule('I3SimpleFitter', 'MillipedeFit', SeedService='seed',
    Parametrization='millipedeparam', LogLikelihood='bipedllh',
    Minimizer='minuit')

tray.AddModule('I3Writer', 'writer', filename=sys.argv[1])
tray.AddModule('TrashCan','can')
tray.Execute(6)
tray.Finish()

