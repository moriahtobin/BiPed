#!/usr/bin/env python
from I3Tray import *
import sys
from icecube import icetray, dataio, dataclasses, photonics_service, gulliver
from icecube import millipede, wavedeform
load('gulliver-modules')
load('gulliver')
load('WaveCalibrator')
load('libparticleforge')

if len(sys.argv) < 3:
	print 'Usage: %s output.i3 input1.i3 [input2.i3] ...'
	sys.exit(1)

files = sys.argv[2:]

#muon_service = photonics_service.I3PhotoSplineService('/net/user/mntobin/IceRec/Tables/ZeroLengthMuons_z20_a10_150.abs.fits', '/net/user/mntobin/IceRec/Tables/ZeroLengthMuons_z20_a10_150.prob.fits', 0)
muon_service = photonics_service.I3PhotoSplineService('/net/user/mntobin/IceRec/Tables/emu_abs150.fits', '/net/user/mntobin/IceRec/Tables/emu_prob150.fits', 0)
#cascade_service = photonics_service.I3PhotoSplineService('/net/user/mntobin/IceRec/Tables/ems_spice1_z20_a10_150.abs.fits', '/net/user/mntobin/IceRec/Tables/ems_spice1_z20_a10_150.prob.fits', 0)
cascade_service_mie = photonics_service.I3PhotoSplineService('/net/user/mntobin/IceRec/Tables/ems_mie_z20_a10_150.abs.fits', '/net/user/mntobin/IceRec/Tables/ems_mie_z20_a10_150.prob.fits', 0)
print "Yum Photonics Tables!"

tray = I3Tray()


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


class BipedGuessSeedService(gulliver.I3SeedService):
        def __init__(self, context):
                gulliver.I3SeedService.__init__(self, context)
                self.AddParameter("FirstSeed", "First seed name", "")
                self.AddParameter("SecondSeed", "Second seed name", "")

        def Configure(self):
                self.name1=GetParameter("FirstSeed")
                self.name2=GetParameter("SecondSeed")

        def SetEvent(frame):
                self.seed1=frame[self.name1]
                self.seed2=frame[self.name2]
                return 2

        def GetSeed(num):
                if num==0:
                        return self.seed1
                elif num==1:
                        return self.seed2
                else:
                        log_fatal("bad seed index")



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

tray.AddModule(Hybridforge,seed='Monopod', factor=4.5, output='LongMuon')
tray.AddModule(Hybridforge,seed='Monopod',factor=2.5, output='EqualMuon')
tray.AddService(BipedGuessSeedService, 'BestSeed',FirstSeed='LongMuon',SecondSeed='EqualMuon')


tray.AddService('MuMillipedeParametrizationFactory', 'MuMillipede',
    StepT=5, 
    StepX=5, RelativeBoundsX=[-50.0,50.0],
    StepY=5, RelativeBoundsY=[-50.0,50.0],
    StepZ=5, RelativeBoundsZ=[-50.0,50.0],
    StepAzimuth=0.3, BoundsAzimuth=[-0.61,7.],
    StepZenith=0.2, BoundsZenith=[-0.41,3.55],
    StepLogE=0.05, BoundsLogE=[0,4],
    StepLogL=0.05, BoundsLogL=[0,3],
    MuonSpacing=15, ShowerSpacing=15, 
#    StartingCascadeStepSize=0.4
)
tray.AddService('MuMillipedeParametrizationFactory', 'Static',
    StepT=0, 
    StepX=0, 
    StepY=0, 
    StepZ=0, 
    StepAzimuth=0.0, 
    StepZenith=0.0, 
    StepLogE=0.0, 
    StepLogL=0.0, 
    MuonSpacing=15, ShowerSpacing=15, 
)
tray.AddService('MillipedeLikelihoodFactory', 'Mil-llh',
    MuonPhotonicsService=muon_service, CascadePhotonicsService=cascade_service_mie,
    PhotonsPerBin=5, Pulses='OfflinePulses_NoBorkedSLC')
tray.AddService('I3GSLRandomServiceFactory','I3RandomService')
tray.AddService('I3GulliverMinuit2Factory', 'minuit',
    MaxIterations=3000, 
#    Algorithm='SIMPLEX',
#    MinuitPrintLevel=-1, 
    Algorithm='MIGRAD', MinuitStrategy=0, 
    WithGradients=True,
    IgnoreEDM=True,
    FlatnessCheck=True,
 #   CheckGradient=True,
    Tolerance=0.005)

	
tray.AddService('I3BasicSeedServiceFactory', 'seed', 
   FirstGuess='BestSeed',
# FirstGuess='MPEFitEuler_Contained',
    TimeShiftType='TNone')
tray.AddModule('I3SimpleFitter', 'PreFit', SeedService='BestSeed',
    Parametrization='Static', LogLikelihood='Mil-llh',
    Minimizer='minuit')
tray.AddModule('I3SimpleFitter', 'MillipedeFit', SeedService='seed',
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


