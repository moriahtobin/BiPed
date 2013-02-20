#!/usr/bin/env python
#Hello Laura :)
from I3Tray import *
import sys
from icecube import icetray, dataio, dataclasses, photonics_service, gulliver
from icecube import millipede, wavedeform
load('gulliver-modules')
load('BiPed')
load('WaveCalibrator')

if len(sys.argv) < 3:
	print 'Usage: %s output.i3 input1.i3 [input2.i3] ...'
	sys.exit(1)

files = sys.argv[2:]

muon_service = photonics_service.I3PhotoSplineService('/net/user/mntobin/IceRec/Tables/emu_abs150.fits', '/net/user/mntobin/IceRec/Tables/emu_prob150.fits', 0)
cascade_service = photonics_service.I3PhotoSplineService('/net/user/mntobin/IceRec/Tables/ems_spice1_z20_a10_150.abs.fits', '/net/user/mntobin/IceRec/Tables/ems_spice1_z20_a10_150.prob.fits', 0)
cascade_service_mie = photonics_service.I3PhotoSplineService('/net/user/mntobin/IceRec/Tables/ems_mie_z20_a10_150.abs.fits', '/net/user/mntobin/IceRec/Tables/ems_mie_z20_a10_150.prob.fits', 0)

tray = I3Tray()
tray.AddModule('I3Reader', 'reader', FilenameList=files)
tray.AddModule("I3WaveformTimeRangeCalculator","WaveformTimeRange")
tray.AddModule(wavedeform.AddMissingTimeWindow,"PulseTimeRange",Streams=[icetray.I3Frame.DAQ])


tray.AddService('BipedParametrizationFactory', 'bipedparam',
    MuonSpacing=1, ShowerSpacing=1, StepT=0, StepX=2, StepY=2, StepZ=2,    
# DirectionStepsize=0.02, FitBiped=True)
    StepAzimuth=0.02, StepZenith=0.02, FitContainedTrack=True
    #LinLengthStepSize=1
)
tray.AddService('BipedLikelihoodFactory', 'bipedllh',
    MuonPhotonicsService=muon_service, CascadePhotonicsService=cascade_service_mie,
    PhotonsPerBin=1, Pulses='OfflinePulses')
tray.AddService('I3GSLRandomServiceFactory','I3RandomService')
tray.AddService('I3GulliverMinuit2Factory', 'minuit',
    MaxIterations=2000, Algorithm='SIMPLEX')
tray.AddService('I3BasicSeedServiceFactory', 'seed', FirstGuess='SPEFit4_DC',
    TimeShiftType='TNone')
tray.AddModule('I3SimpleFitter', 'MillipedeFit', SeedService='seed',
    Parametrization='bipedparam', LogLikelihood='bipedllh',
    Minimizer='minuit')


global i
i=1
def count(frame):
	global i 
	print "ran on", i, "frames"
	i = i+1
tray.AddModule(count,"mycounter")

tray.AddModule('I3Writer', 'writer', filename=sys.argv[1])
tray.AddModule('TrashCan','can')
print "Got the tray put together, start running..."
tray.Execute()
tray.Finish()

