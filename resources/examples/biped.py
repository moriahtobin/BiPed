#!/usr/bin/env python
#Hello Laura :)
from I3Tray import *
import sys
from icecube import icetray, dataio, dataclasses, photonics_service, gulliver
from icecube import millipede, wavedeform
load('gulliver-modules')
load('BiPed')
load('WaveCalibrator')
load('libparticleforge')

if len(sys.argv) < 3:
	print 'Usage: %s output.i3 input1.i3 [input2.i3] ...'
	sys.exit(1)

files = sys.argv[2:]

muon_service = photonics_service.I3PhotoSplineService('/net/user/mntobin/IceRec/Tables/ZeroLengthMuons_z20_a10_150.abs.fits', '/net/user/mntobin/IceRec/Tables/ZeroLengthMuons_z20_a10_150.prob.fits', 0)
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


def Hybridforge(frame, seed, lengthseed, output):
    if frame.Has(seed):
        source = frame[seed]
	lseed = frame[lengthseed]
#        energy = 0.
#        if frame.Has(seed + 'Params'):
#            for particle in frame[seed + 'Params'] : energy += particle.energy
#        else:
#            energy = source.energy
        forged_particle = dataclasses.I3Particle()
        forged_particle.energy = source.energy
	forged_particle.length = (4.5*source.energy + lseed.length)/2.0
        forged_particle.dir.set_direction(source.dir.zenith, source.dir.azimuth)
        forged_particle.pos.x = source.pos.x
        forged_particle.pos.y = source.pos.y
        forged_particle.pos.z = source.pos.z
        forged_particle.time = source.time
        forged_particle.speed = 0.29979245799999998
        forged_particle.shape = dataclasses.I3Particle.ParticleShape.ContainedTrack
        forged_particle.fit_status = dataclasses.I3Particle.FitStatus.OK
        frame.Put(output, forged_particle)



tray.AddModule('I3ParticleForgeModule', 'most_energetic_primary',
    Shape=        'MC',
    Type=        'MC',
    Time=        'MC',
    Position=     'MC',
    Direction=    'MC',
    Energy=       'MC',
    MCTree=       'I3MCTree',
    MCMethod=     'MOSTENERGETICPRIMARY', # this is the method you have to choose, can be the following:
        # MOSTENERGETICCASCADE, MOSTENERGETICPRIMARY, MOSTENERGETICTRACK, REFERENCECASCADE_DEPOSITED, REFERENCECASCADE_VISIBLE
    output=       'I3MCPrimary')

tray.AddModule(Hybridforge,seed='I3MCPrimary',lengthseed='MPEFitEuler_Contained',output='lenSeed')

tray.AddService('BipedParametrizationFactory', 'bipedparam',
    StepT=0,
    StepX=5, RelativeBoundsX=[-50.0,50.0],
    StepY=5, RelativeBoundsY=[-50.0,50.0],
    StepZ=5, RelativeBoundsZ=[-50.0,50.0],
    StepAzimuth=0.3, RelativeBoundsAzimuth=[-3,3],
    StepZenith=0.2, RelativeBoundsZenith=[-6,6],
    StepLogE=0.05, BoundsLogE=[0,4],
    StepLogL=0.1, BoundsLogL=[0,3],
#    StepAzimuth=0.2, BoundsAzimuth=[-1,7.5],
#    StepZenith=0.05, BoundsZenith=[-0.11,3.25],
    StartingCascadeStepSize=0.4
)
tray.AddService('BipedLikelihoodFactory', 'bipedllh',
    MuonPhotonicsService=muon_service, CascadePhotonicsService=cascade_service_mie,
    PhotonsPerBin=5, MuonSpacing=3, Pulses='OfflinePulses_NoBorkedSLC')
tray.AddService('I3GSLRandomServiceFactory','I3RandomService')
tray.AddService('I3GulliverMinuit2Factory', 'minuit',
    MaxIterations=3000, Algorithm='MIGRAD', MinuitStrategy=2, 
    WithGradients=True,
    FlatnessCheck=True,
    IgnoreEDM=True,
    CheckGradient=False,
#    WithGradients=True,
#    CheckGradient=True,
    Tolerance=0.1)
tray.AddService('I3BasicSeedServiceFactory', 'seed', 
    FirstGuess='lenSeed',
#    FirstGuess='MPEFitEuler_Contained',
# FirstGuess='I3MCPrimary',
    TimeShiftType='TNone')
tray.AddModule('I3SimpleFitter', 'BiPedHardCodeFit', SeedService='seed',
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

