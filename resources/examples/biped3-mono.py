#!/usr/bin/env python
from I3Tray import *
import sys
from icecube import icetray, dataio, dataclasses, photonics_service, gulliver
from icecube import millipede, wavedeform
load('gulliver-modules')
load('WaveCalibrator')
load('libparticleforge')

if len(sys.argv) < 3:
	print 'Usage: %s output.i3 input1.i3 [input2.i3] ...'
	sys.exit(1)

files = sys.argv[2:]

muon_service = photonics_service.I3PhotoSplineService('/net/user/mntobin/IceRec/Tables/ZeroLengthMuons_z20_a10_150.abs.fits', '/net/user/mntobin/IceRec/Tables/ZeroLengthMuons_z20_a10_150.prob.fits', 0)
#muon_service = photonics_service.I3PhotoSplineService('/net/user/mntobin/IceRec/Tables/emu_abs150.fits', '/net/user/mntobin/IceRec/Tables/emu_prob150.fits', 0)
#cascade_service = photonics_service.I3PhotoSplineService('/net/user/mntobin/IceRec/Tables/ems_spice1_z20_a10_150.abs.fits', '/net/user/mntobin/IceRec/Tables/ems_spice1_z20_a10_150.prob.fits', 0)
cascade_service_mie = photonics_service.I3PhotoSplineService('/net/user/mntobin/IceRec/Tables/ems_mie_z20_a10_150.abs.fits', '/net/user/mntobin/IceRec/Tables/ems_mie_z20_a10_150.prob.fits', 0)
print "Yum Photonics Tables!"

tray = I3Tray()
tray.AddModule('I3Reader', 'reader', FilenameList=files)

def Hybridforge(frame, seed, factor, output):
    if frame.Has(seed):
        source = frame[seed]
#        energy = 0.
#        if frame.Has(seed + 'Params'):
#            for particle in frame[seed + 'Params'] : energy += particle.energy
#        else:
#            energy = source.energy
        forged_particle = dataclasses.I3Particle()
        forged_particle.energy = source.energy
	forged_particle.length = source.energy*factor
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

#tray.AddModule(Hybridforge,seed='SPEFitSingle_DC',lengthseed='MPEFitEuler_Contained',output='lenSeed')
tray.AddModule(Hybridforge,seed='Monopod', factor=4.5, output='LongMuon')
tray.AddModule(Hybridforge,seed='Monopod',factor=2.5, output='EqualMuon')

tray.AddService('MuMillipedeParametrizationFactory', 'MuMillipede',
    StepT=5, 
    StepX=5, RelativeBoundsX=[-50.0,50.0],
    StepY=5, RelativeBoundsY=[-50.0,50.0],
    StepZ=5, RelativeBoundsZ=[-50.0,50.0],
    StepAzimuth=0.3, BoundsAzimuth=[-0.61,7.0],
    StepZenith=0.2, BoundsZenith=[-0.41,3.55],
    StepLogE=0.05, BoundsLogE=[0,4],
    StepLogL=0.05, BoundsLogL=[0,3],
    MuonSpacing=3, ShowerSpacing=100000000, StartingCascadeStepSize=0.4)
tray.AddService('MillipedeLikelihoodFactory', 'Mil-llh',
    MuonPhotonicsService=muon_service, CascadePhotonicsService=cascade_service_mie,
    PhotonsPerBin=5, Pulses='OfflinePulses_NoBorkedSLC')
tray.AddService('I3GSLRandomServiceFactory','I3RandomService')
tray.AddService('I3GulliverMinuit2Factory', 'minuit',
    MaxIterations=1000, 
#    Algorithm="SIMPLEX", MinuitPrintLevel=-1,
    Algorithm='MIGRAD', MinuitStrategy=0, 
    WithGradients=True,
    FlatnessCheck=True,
    IgnoreEDM=True,
    CheckGradient=False,
#    Tolerance=0.0001)
    Tolerance=0.005)

#tray.AddModule('I3ParticleForgeModule', 'LowEn',
 #   Shape=        'ContainedTrack',
 #   Type=        'MC',
#    Time=        'MPEFitEuler_Contained',
#    Position=     'MPEFitEuler_Contained',
#    Direction=    'MPEFitEuler_Contained',
#    Speed= 'MPEFitEuler_Contained',
 #   Energy=       25, #monopod?
 #   MCTree=       'I3MCTree',
 #   MCMethod=     'MOSTENERGETICPRIMARY', # this is the method you have to choose, can be the following:
        # MOSTENERGETICCASCADE, MOSTENERGETICPRIMARY, MOSTENERGETICTRACK, REFERENCECASCADE_DEPOSITED, REFERENCECASCADE_VISIBLE
#    output=       'NuSeed')

tray.AddService('I3BasicSeedServiceFactory', 'seed', 
   FirstGuess='LongMuon',
# FirstGuess='MPEFitEuler_Contained',
    TimeShiftType='TNone')
tray.AddModule('I3SimpleFitter', 'ActualBipedInMillipedeFit', SeedService='seed',
    Parametrization='MuMillipede', LogLikelihood='Mil-llh',
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

