#include <millipede/Millipede.h>

#include <lilliput/parametrization/I3HalfSphereParametrization.h>
#include <gulliver/I3ParametrizationBase.h>
#include <gulliver/I3EventHypothesis.h>
#include <icetray/I3ServiceFactory.h>
#include <icetray/I3SingleServiceFactory.h>

static void MuMillipedeHypothesis(I3ParticleConstPtr track,
    std::vector<I3Particle> &hypothesis, double boundary,
    double muonspacing, double showerspacing, bool startingtrack,
    double slantstart, double slantstop);

class MuMillipede : public I3MillipedeConditionalModule
{
	public:
		MuMillipede(const I3Context &);
		void Configure();
		void Hypothesis(I3FramePtr frame,
		    std::vector<I3Particle> &hypothesis);
		void Physics(I3FramePtr frame);
	private:
		SET_LOGGER("MuMillipede");

		std::string seed_;
		double boundary_;
		double muonspacing_;
		double showerspacing_;
		std::string output_name;
		bool fit_starting_track_;
		double slantstart_;
		double slantstop_;
};

class MuMillipedeParametrization : public I3HalfSphereParametrization {
	public:
		MuMillipedeParametrization(const I3Context &);
		void Configure();

		double boundary_;
		double muonspacing_;
		double showerspacing_;
		std::string name_;
		bool fit_starting_track_;
		double slantstart_;
		double slantstop_;

		void UpdatePhysicsVariables();
		void UpdateParameters();
		void ApplyChainRule();
	private:
		SET_LOGGER("MuMillipedeParametrization");
};

typedef
    I3SingleServiceFactory<MuMillipedeParametrization, I3ParametrizationBase>
    MuMillipedeParametrizationFactory;
I3_SERVICE_FACTORY(MuMillipedeParametrizationFactory);
I3_MODULE(MuMillipede);

MuMillipede::MuMillipede(const I3Context &context) :
    I3MillipedeConditionalModule(context)
{
	AddOutBox("OutBox");

	AddParameter("SeedTrack", "Track to fit", "");
	AddParameter("Boundary", "Segment boundary, in meters (fits segments "
	    "in [-x,x] in x,y,z in detector coordinates)", 600);
	AddParameter("MuonSpacing", "Spacing of muon (ionization) sources "
	    "along the track, in meters. MUST match source extension in "
	    "muon tables. Set to 0 to not include ionization sources.", 15);
	AddParameter("ShowerSpacing", "Spacing of shower (radiative) "
	    "sources along the track, in meters. Set to 0 to not include "
	    "radiative losses.", 15);
	AddParameter("Output", "Name of output", "MillipededEdx");
	AddParameter("FitStartingTrack", "Parameterize a track starting at "
	    "the fit time at the vertex instead of a throughgoing track",
	    false);
	AddParameter("StartSlantDepth", "Start of segment boundary in slant "
	    "depth bins (fits segments in [Xstart,Xend] in detector "
	    "coordinates). If negative, boundary used instead.", -1);
	AddParameter("EndSlantDepth",  "End of Segment boundary in slant depth "
	    "bins (fits segments in [Xstart,Xend] in detector coordinates). "
	    "If negative, boundary used instead.", -1);
}

void
MuMillipede::Configure()
{
	I3MillipedeConditionalModule::Configure();

	GetParameter("SeedTrack", seed_);
	GetParameter("Boundary", boundary_);
	GetParameter("MuonSpacing", muonspacing_);
	GetParameter("ShowerSpacing", showerspacing_);
	GetParameter("Output", output_name);
	GetParameter("FitStartingTrack", fit_starting_track_);
	GetParameter("StartSlantDepth", slantstart_);
	GetParameter("EndSlantDepth", slantstop_);
	if (timeBinPhotons_ >= 0 && timeBinPhotons_ < 5)
		log_warn("PhotonsPerBin set to %f, but should not be "
		    "less than 5!", timeBinPhotons_);
}

void
MuMillipede::Physics(I3FramePtr frame)
{
	if (!frame->Has(pulses_name_) || !frame->Has(seed_)) {
		PushFrame(frame);
		return;
	}

	boost::shared_ptr<I3Vector<I3Particle> > sources(new
	    I3Vector<I3Particle>);
	MuMillipedeHypothesis(frame->Get<I3ParticleConstPtr>(seed_), *sources,
	    boundary_, muonspacing_, showerspacing_, fit_starting_track_,
	    slantstart_, slantstop_);
	if (sources->size() == 0) {
		PushFrame(frame);
		return;
	}

	cholmod_sparse *response_matrix;
	DatamapFromFrame(*frame);
	response_matrix = GetResponseMatrix(*sources);
	if (response_matrix == NULL)
		log_fatal("Null basis matrix");
	SolveEnergyLosses(*sources, response_matrix);

	MillipedeFitParamsPtr params(new MillipedeFitParams);
	Millipede::FitStatistics(domCache_, *sources, 0, response_matrix,
	    &*params, &c);
	cholmod_l_free_sparse(&response_matrix, &c);

	frame->Put(output_name, sources);
	frame->Put(output_name + "FitParams", params);

	PushFrame(frame);
}

MuMillipedeParametrization::MuMillipedeParametrization(const I3Context& context)
    : I3HalfSphereParametrization(context)
{
	AddParameter("Boundary", "Segment boundary, in meters (fits segments "
	    "within this number of meters of the vertex)", 600);
	AddParameter("MuonSpacing", "Spacing of muon (ionization) sources "
	    "along the track, in meters. MUST match source extension in "
	    "muon tables. Set to 0 to not include ionization sources.", 15);
	AddParameter("ShowerSpacing", "Spacing of shower (radiative) "
	    "sources along the track, in meters. Set to 0 to not include "
	    "radiative losses.", 15);
	AddParameter("FitStartingTrack", "Parameterize a track starting at "
	    "the fit time at the vertex instead of a throughgoing track",
	    false);
	AddParameter("StartSlantDepth", "Start of segment boundary in slant "
	    "depth bins (fits segments in [Xstart,Xend] in detector "
	    "coordinates). If negative, boundary used instead.", -1);
	AddParameter("EndSlantDepth",  "End of Segment boundary in slant depth "
	    "bins (fits segments in [Xstart,Xend] in detector coordinates). "
	    "If negative, boundary used instead.", -1);
}
	
void
MuMillipedeParametrization::UpdatePhysicsVariables()
{
	// Do the hard work
	I3HalfSphereParametrization::UpdatePhysicsVariables();

	// Add the hypothesis vector
	boost::shared_ptr<I3Vector<I3Particle> > sources(new
	    I3Vector<I3Particle>);
	MuMillipedeHypothesis(hypothesis_->particle, *sources,
	    boundary_, muonspacing_, showerspacing_,
	    fit_starting_track_, slantstart_, slantstop_);

	hypothesis_->nonstd = sources;

	if (gradient_) {
		boost::shared_ptr<I3Vector<I3Particle> > gradsources(new
		    I3Vector<I3Particle>);
		gradsources->resize(sources->size());
		gradient_->nonstd = gradsources;
	}
}

void
MuMillipedeParametrization::UpdateParameters()
{
	I3HalfSphereParametrization::UpdateParameters();
	UpdatePhysicsVariables();
}

void
MuMillipedeParametrization::ApplyChainRule()
{
	I3Particle& gradient = *(gradient_->particle);
	boost::shared_ptr<I3Vector<I3Particle> > gradsources =
	  boost::dynamic_pointer_cast<I3Vector<I3Particle> >(gradient_->nonstd);
	boost::shared_ptr<I3Vector<I3Particle> > sources =
	  boost::dynamic_pointer_cast<I3Vector<I3Particle> >
	  (hypothesis_->nonstd);

	double grad_x = 0, grad_y = 0, grad_z = 0, grad_t = 0;
	double grad_zen = 0, grad_azi = 0;
	for (unsigned i = 0; i < gradsources->size(); i++) {
		I3Particle &gradpart = (*gradsources)[i];
		I3Particle &part = (*sources)[i];

		grad_x += gradpart.GetX();
		grad_y += gradpart.GetY();
		grad_z += gradpart.GetZ();
		grad_t += gradpart.GetTime();
		grad_zen += gradpart.GetZenith();
		grad_azi += gradpart.GetAzimuth();

		// X, Y, Z of sub-cascades depend on theta, phi
		double zenith_xyz = 0, azimuth_xyz = 0;
		double theta = part.GetZenith();
		double phi = part.GetAzimuth();
		zenith_xyz += gradpart.GetX()*cos(theta)*cos(phi);
		azimuth_xyz += -gradpart.GetX()*sin(theta)*sin(phi);
		zenith_xyz += gradpart.GetY()*cos(theta)*sin(phi);
		azimuth_xyz += gradpart.GetY()*cos(theta)*cos(phi);
		zenith_xyz += -gradpart.GetZ()*sin(theta);

		// Shower (X, Y, Z) is track *minus* c_(x y z) t
		grad_zen -= part.GetSpeed()*(part.GetTime() -
		    hypothesis_->particle->GetTime())*zenith_xyz;
		grad_azi -= part.GetSpeed()*(part.GetTime() -
		    hypothesis_->particle->GetTime())*azimuth_xyz;
	}

	gradient.SetPos(grad_x, grad_y, grad_z);
	gradient.SetTime(grad_t);
	gradient.SetDir(grad_zen, grad_azi);

	I3HalfSphereParametrization::ApplyChainRule();
}

void
MuMillipedeParametrization::Configure()
{
	I3HalfSphereParametrization::Configure();

	GetParameter("Boundary", boundary_);
	GetParameter("MuonSpacing", muonspacing_);
	GetParameter("ShowerSpacing", showerspacing_);
	GetParameter("FitStartingTrack", fit_starting_track_);
	GetParameter("StartSlantDepth", slantstart_);
	GetParameter("EndSlantDepth", slantstop_);
}

static void
MuMillipedeHypothesis(I3ParticleConstPtr track,
    std::vector<I3Particle> &hypothesis, double boundary,
    double muonspacing, double showerspacing, bool startingtrack,
    double slantstart, double slantstop)
{
	double cross_times[6]; // x1, x2, y1, y2, z1, z2
	double xspeed, yspeed, zspeed;
	double entry_t, exit_t;

	// Projections of particle speed on each axis
	xspeed = track->GetSpeed()*sin(track->GetZenith() / I3Units::radian)*
	    cos(track->GetAzimuth() / I3Units::radian);
	yspeed = track->GetSpeed()*sin(track->GetZenith() / I3Units::radian)*
	    sin(track->GetAzimuth() / I3Units::radian);
	zspeed = track->GetSpeed()*cos(track->GetZenith() / I3Units::radian);

	// Get the time points where it crosses the "detector" boundaries
	// (defined as a cube of width 2*boundary_, centered at the detector
	// origin)

	cross_times[0] = track->GetTime() + (track->GetX() - boundary)/xspeed;
	cross_times[1] = track->GetTime() + (track->GetX() + boundary)/xspeed;
	cross_times[2] = track->GetTime() + (track->GetY() - boundary)/yspeed;
	cross_times[3] = track->GetTime() + (track->GetY() + boundary)/yspeed;
	cross_times[4] = track->GetTime() + (track->GetZ() - boundary)/zspeed;
	cross_times[5] = track->GetTime() + (track->GetZ() + boundary)/zspeed;

	// Sort the crossing points. The middle two are when it enters and
	// exits the volume.
	std::sort(&cross_times[0], &cross_times[6]);
	entry_t = cross_times[2];
	exit_t = cross_times[3];

	if (startingtrack) {
		entry_t = track->GetTime();
	} else if (slantstart >= 0 && slantstop > slantstart) {
		// Track equation : r = r0 + v*(t-t0)
		// For slantstart, (zIceTop - z) = slantstart* cos(zenith),
		// thus we know entry_t : -zspeed, because this speed was
		// pointing in the opposite direction
		entry_t = track->GetTime() +
		    (I3Constants::zIceTop -
		     slantstart*cos(track->GetZenith()/I3Units::radian) 
		     - track->GetZ())/-zspeed;  
		exit_t  = track->GetTime() + 
		    (I3Constants::zIceTop -
		     slantstop * cos(track->GetZenith()/I3Units::radian)  
		     - track->GetZ())/-zspeed;
	} else {
		// Shift the entry point to be an integer number of the largest 
		// segment spacing from the vertex. This simplifies gradient
		// calculation, since it makes x - x0 fixed under small
		// displacements.
		double seg_duration = ((muonspacing > showerspacing) ?
		    muonspacing : showerspacing)/track->GetSpeed();
		entry_t = round((entry_t - track->GetTime())/seg_duration)*
		    seg_duration + track->GetTime();
	}
	
	if (muonspacing > 0) for (double t = entry_t; t < exit_t;
	    t += muonspacing/track->GetSpeed()) {
		I3Particle muon;
		muon.SetDir(track->GetDir());
		muon.SetLength(muonspacing);
		muon.SetType(track->GetType());
		muon.SetShape(I3Particle::ContainedTrack);
		muon.SetSpeed(track->GetSpeed());
		muon.SetTime(t);
		muon.SetPos(track->GetX() - xspeed*(t - track->GetTime()),
			    track->GetY() - yspeed*(t - track->GetTime()),
			    track->GetZ() - zspeed*(t - track->GetTime()));
		hypothesis.push_back(muon);
	}
	if (showerspacing > 0) for (double t = entry_t; t < exit_t;
	    t += showerspacing/track->GetSpeed()) {
		I3Particle shower;
		shower.SetDir(track->GetDir());
		shower.SetLength(showerspacing);
		shower.SetType(track->GetType());
		shower.SetSpeed(track->GetSpeed());
		shower.SetShape(I3Particle::Cascade);
		shower.SetTime(t);
		shower.SetPos(track->GetX() - xspeed*(t - track->GetTime()),
			      track->GetY() - yspeed*(t - track->GetTime()),
			      track->GetZ() - zspeed*(t - track->GetTime()));
		hypothesis.push_back(shower);
	}
}

