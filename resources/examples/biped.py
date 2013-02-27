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
#cascade_service = photonics_service.I3PhotoSplineService('/net/user/mntobin/IceRec/Tables/ems_spice1_z20_a10_150.abs.fits', '/net/user/mntobin/IceRec/Tables/ems_spice1_z20_a10_150.prob.fits', 0)
cascade_service_mie = photonics_service.I3PhotoSplineService('/net/user/mntobin/IceRec/Tables/ems_mie_z20_a10_150.abs.fits', '/net/user/mntobin/IceRec/Tables/ems_mie_z20_a10_150.prob.fits', 0)
print "Yum Photonics Tables!"
tray = I3Tray()
tray.AddModule('I3Reader', 'reader', FilenameList=files)
#fixup code from JVS
def flag_borked_slc(frame):
        # Find SLC launches affected by DOMsimulator off-by-one bug,
        # and tell Millipede to ignore them
        bad_keys = set()
        import numpy
        for om, dls in frame['InIceRawData'].iteritems():
                for dl in dls:
                        if dl.lc_bit == False and dl.charge_stamp_highest_sample == 0:
                                print '%s, she is the bad!' % om
                                bad_keys.add(om)
        if len(bad_keys) > 0:
                frame['BorkedSLC'] = dataclasses.I3VectorOMKey(bad_keys)
                frame['OfflinePulses_NoBorkedSLC'] = dataclasses.I3RecoPulseSeriesMapMask(frame, 'OfflinePulses',
                     lambda om, idx, pulse: om not in bad_keys)
        else:
                frame['OfflinePulses_NoBorkedSLC'] = dataclasses.I3RecoPulseSeriesMapMask(frame, 'OfflinePulses')
tray.AddModule(flag_borked_slc, 'slc_filter', Streams=[icetray.I3Frame.DAQ])

tray.AddModule("I3WaveformTimeRangeCalculator","WaveformTimeRange")
tray.AddModule(wavedeform.AddMissingTimeWindow,"PulseTimeRange")
def copytimerange(frame):
	frame['OfflinePulses_NoBorkedSLCTimeRange']=frame['OfflinePulsesTimeRange']
tray.AddModule(copytimerange, "CopyTimeRange")

tray.AddService('BipedParametrizationFactory', 'bipedparam',
    MuonSpacing=1, ShowerSpacing=1, StepT=0, StepX=2, StepY=2, StepZ=2,    
# DirectionStepsize=0.02, FitBiped=True)
    StepAzimuth=0.02, StepZenith=0.02
    #LinLengthStepSize=1
)
tray.AddService('BipedLikelihoodFactory', 'bipedllh',
    MuonPhotonicsService=muon_service, CascadePhotonicsService=cascade_service_mie,
    PhotonsPerBin=1, Pulses='OfflinePulses_NoBorkedSLC')
tray.AddService('I3GSLRandomServiceFactory','I3RandomService')
tray.AddService('I3GulliverMinuit2Factory', 'minuit',
    MaxIterations=2000, Algorithm='SIMPLEX')
tray.AddService('I3BasicSeedServiceFactory', 'seed', FirstGuess='MPEFitEuler_Contained',
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
tray.Execute(11)
tray.Finish()

