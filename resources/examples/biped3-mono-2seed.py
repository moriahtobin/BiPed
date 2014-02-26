#!/usr/bin/env python
from I3Tray import *
import sys
from icecube import icetray, dataio, dataclasses, photonics_service, gulliver
from icecube import millipede, wavedeform
load('gulliver-modules')
load('libparticleforge')

if len(sys.argv) < 3:
	print 'Usage: %s output.i3 input1.i3 [input2.i3] ...'
	sys.exit(1)

files = sys.argv[2:]



#Photonics!
muon_service = photonics_service.I3PhotoSplineService('/net/user/mntobin/IceRec/Tables/ZeroLengthMuons_z20_a10_150.abs.fits', '/net/user/mntobin/IceRec/Tables/ZeroLengthMuons_z20_a10_150.prob.fits', 0)
#muon_service = photonics_service.I3PhotoSplineService('/net/user/mntobin/IceRec/Tables/emu_abs150.fits', '/net/user/mntobin/IceRec/Tables/emu_prob150.fits', 0)
#cascade_service = photonics_service.I3PhotoSplineService('/net/user/mntobin/IceRec/Tables/ems_spice1_z20_a10_150.abs.fits', '/net/user/mntobin/IceRec/Tables/ems_spice1_z20_a10_150.prob.fits', 0)
cascade_service_mie = photonics_service.I3PhotoSplineService('/net/user/mntobin/IceRec/Tables/ems_mie_z20_a10_150.abs.fits', '/net/user/mntobin/IceRec/Tables/ems_mie_z20_a10_150.prob.fits', 0)
print "Yum Photonics Tables!"





#Begin:
tray = I3Tray()
tray.AddModule('I3Reader', 'reader', FilenameList=files)






#Seed Creation Module:
def Hybridforge(frame, seed, factor, output):
    if frame.Has(seed):
        source = frame[seed]
        forged_particle = dataclasses.I3Particle()
        forged_particle.energy = source.energy
	forged_particle.length = source.energy*factor
        forged_particle.dir = source.dir
        forged_particle.pos.x = source.pos.x
        forged_particle.pos.y = source.pos.y
        forged_particle.pos.z = source.pos.z
        forged_particle.time = source.time
        forged_particle.speed = 0.29979245799999998
        forged_particle.shape = dataclasses.I3Particle.ParticleShape.ContainedTrack
        forged_particle.fit_status = dataclasses.I3Particle.FitStatus.OK
        frame.Put(output, forged_particle)


#MC seed creation module:
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



#Module for determining best seed:
def BiPedChooser(frame, winner="",seedlist=[]):
	bestlogl=1e100
	winseed=""
	for seedname in seedlist:
		if frame.Has(seedname) and frame.Has(seedname+"FitParams"):
			fit=frame[seedname]
			fitparams=frame[seedname+"FitParams"]
			if fit.fit_status==dataclasses.I3Particle.FitStatus.OK and fitparams.logl<bestlogl:
				bestlogl=fitparams.logl
				winseed=seedname+'Muon'
	if winseed:
		print "%s is the winning seed" % winseed
		frame.Put(winner,frame[winseed])
	else:
		for seedname in seedlist:
			if not frame.Has(seedname):
				print "fit %s not found" % seedname
			elif not frame.Has(seedname+"FitParams"):
				print "fitparams of %s not found" %s
			else:
				fit=frame[seedname]
				fitparams=fram[seedname+"FitParams"]
				if fit.fit_status!=dataclasses.I3Particle.FitStatus.OK:
					print "%s has bad fitstatus %d (=> '%s')" % (seedname, fit.fit_status,fit.fit_status_string)
				else:
					print "perchance %s has a weird llh value? logl=%f" % (seedname,fitparams.logl)



#Seeds with most common muon/cascade energy sharing patterns (most energy in muon vs. equal sharing between muon and cascade)
#llh does not change for length minimization with StepL < MuonSpacing, so we want our seed to be as close to the actual value as possible
tray.AddModule(Hybridforge,seed='MonopodFit', factor=4.5, output='LongMuon')
tray.AddModule(Hybridforge,seed='MonopodFit',factor=2.3, output='EqualMuon')
#tray.AddModule(Hybridforge,seed='Monopod', factor=4.5, output='LongMuon')
#tray.AddModule(Hybridforge,seed='Monopod',factor=2.3, output='EqualMuon')



#Parametrizations:
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
tray.AddService('MuMillipedeParametrizationFactory', 'static',
    StepT=0, 
    StepX=0, 
    StepY=0, 
    StepZ=0, 
    StepAzimuth=0.0, 
    StepZenith=0.0, 
    StepLogE=0.0, 
    StepLogL=0.0, 
    MuonSpacing=3, ShowerSpacing=100000000, 
)


#Likelihood:
tray.AddService('MillipedeLikelihoodFactory', 'Mil-llh',
    MuonPhotonicsService=muon_service, CascadePhotonicsService=cascade_service_mie,
    PhotonsPerBin=5, Pulses='OfflinePulses_NoBorkedSLC')


#Minimizers:
tray.AddService('I3GSLRandomServiceFactory','I3RandomService')
tray.AddService('I3GulliverMinuit2Factory', 'minuit',
    MaxIterations=1000, 
    Algorithm='MIGRAD', MinuitStrategy=0, 
    WithGradients=True,
    FlatnessCheck=True,
    IgnoreEDM=True,
    CheckGradient=False,
    Tolerance=0.005
)
tray.AddService('I3GulliverMinuit2Factory', 'NoEDM',
    MaxIterations=1000, 
    Algorithm='MIGRAD', MinuitStrategy=0, 
    WithGradients=True,
    FlatnessCheck=True,
    IgnoreEDM=True,
    CheckGradient=False
)

#Set up seeds for test:
tray.AddService('I3BasicSeedServiceFactory', 'long', 
    FirstGuess='LongMuon',
    TimeShiftType='TNone')
tray.AddService('I3BasicSeedServiceFactory', 'equal', 
    FirstGuess='EqualMuon',
    TimeShiftType='TNone')
tray.AddModule('I3SimpleFitter', 'Long', SeedService='long',
    Parametrization='static',LogLikelihood='Mil-llh',
    Minimizer='NoEDM')
tray.AddModule('I3SimpleFitter', 'Equal', SeedService='equal',
    Parametrization='static',LogLikelihood='Mil-llh',
    Minimizer='NoEDM')



#Find best seed hypothesis:
bestseed='BestSeed'
tray.AddModule(BiPedChooser, 'CoinToss', winner=bestseed, seedlist=['Long', 'Equal'])
tray.AddService('I3BasicSeedServiceFactory', 'seed', 
    FirstGuess=bestseed,
    TimeShiftType='TNone')




#Fit with best seed:
tray.AddModule('I3SimpleFitter', 'ActualBipedInMillipedeFit', SeedService='seed',
    Parametrization='MuMillipede', LogLikelihood='Mil-llh',
    Minimizer='minuit')

#Progress:
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
tray.Execute(10)
tray.Finish()

