#include <millipede/Millipede.h>

#include <phys-services/I3Calculator.h>
#include <lilliput/parametrization/I3SimpleParametrization.h>
#include <gulliver/I3ParametrizationBase.h>
#include <gulliver/I3EventHypothesis.h>
#include <icetray/I3ServiceFactory.h>
#include <icetray/I3SingleServiceFactory.h>

static void BipedHypothesis(I3ParticleConstPtr track,
    std::vector<I3Particle> &hypothesis, double boundary,
    double muonspacing, double showerspacing,
    double slantstart, double slantstop);



class BipedParametrization : public I3SimpleParametrization {
	public:
		BipedParametrization(const I3Context &);
		void Configure();

		double boundary_;
		double muonspacing_;
		double showerspacing_;
		std::string name_;
		double starting_cascade_dirstep_;
		double slantstart_;
		double slantstop_;

		void UpdatePhysicsVariables();
		void UpdateParameters();
		void ApplyChainRule();
	private:
		SET_LOGGER("BipedParametrization");
};

typedef
    I3SingleServiceFactory<BipedParametrization, I3ParametrizationBase>
    BipedParametrizationFactory;
I3_SERVICE_FACTORY(BipedParametrizationFactory);



BipedParametrization::BipedParametrization(const I3Context& context)
    : I3SimpleParametrization(context)
{
	AddParameter("Boundary", "Segment boundary, in meters (fits segments "
	    "within this number of meters of the vertex)", 600*I3Units::m);
	AddParameter("MuonSpacing", "Spacing of muon (ionization) sources "
	    "along the track, in meters. MUST match source extension in "
	    "muon tables. Set to 0 to not include ionization sources.",
	    15*I3Units::m);
	AddParameter("ShowerSpacing", "Spacing of shower (radiative) "
	    "sources along the track, in meters. Set to 0 to not include "
	    "radiative losses.", 15*I3Units::m);
	AddParameter("StartingCascadeStepSize", "Allow the initial cascade to "
	    "point in a different direction than the track, moving it in steps "
	    "this size (starting direction fixed if 0).", 0);
	AddParameter("StartSlantDepth", "Start of segment boundary in slant "
	    "depth bins (fits segments in [Xstart,Xend] in detector "
	    "coordinates). If negative, boundary used instead.", -1);
	AddParameter("EndSlantDepth",  "End of Segment boundary in slant depth "
	    "bins (fits segments in [Xstart,Xend] in detector coordinates). "
	    "If negative, boundary used instead.", -1);
}
	
void
BipedParametrization::UpdatePhysicsVariables()
{
	// Do the hard work
	I3SimpleParametrization::UpdatePhysicsVariables();

	// Add the hypothesis vector
	boost::shared_ptr<I3Vector<I3Particle> > sources(new
	    I3Vector<I3Particle>);
	BipedHypothesis(hypothesis_->particle, *sources,
	    boundary_, muonspacing_, showerspacing_, slantstart_, slantstop_);

	if (starting_cascade_dirstep_ > 0) {
		double p0 = par_[par_.size() - 2];
		double p1 = par_[par_.size() - 1];
		double seedX_, seedY_, seedZ_;
		double perp1X_, perp1Y_, perp1Z_;
		double perp2X_, perp2Y_, perp2Z_;

		const I3Direction& dir = hypothesis_->particle->GetDir();

		seedX_ = dir.GetX();
		seedY_ = dir.GetY();
		seedZ_ = dir.GetZ();
		std::pair<I3Direction,I3Direction> sideways =
		    I3Calculator::GetTransverseDirections(dir);
		const I3Direction &perp1 = sideways.first;
		const I3Direction &perp2 = sideways.second;
		perp1X_ = perp1.GetX();
		perp1Y_ = perp1.GetY();
		perp1Z_ = perp1.GetZ();
		perp2X_ = perp2.GetX();
		perp2Y_ = perp2.GetY();
		perp2Z_ = perp2.GetZ();
		
		double newdirx = seedX_ + p0 * perp1X_ + p1 * perp2X_;
		double newdiry = seedY_ + p0 * perp1Y_ + p1 * perp2Y_;
		double newdirz = seedZ_ + p0 * perp1Z_ + p1 * perp2Z_;
		(*sources)[0].SetDir(newdirx, newdiry, newdirz);
	}

	hypothesis_->nonstd = sources;

	if (gradient_) {
		boost::shared_ptr<I3Vector<I3Particle> > gradsources(new
		    I3Vector<I3Particle>);
		gradsources->resize(sources->size());
		gradient_->nonstd = gradsources;
	}
}

void
BipedParametrization::UpdateParameters()
{
	I3SimpleParametrization::UpdateParameters();
	if (starting_cascade_dirstep_ > 0) {
		par_[par_.size() - 2] = 0;
		par_[par_.size() - 1] = 0;
	}
	UpdatePhysicsVariables();
}

void
BipedParametrization::ApplyChainRule()
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

	I3SimpleParametrization::ApplyChainRule();

	if (starting_cascade_dirstep_ > 0) {
		// Gradient for initial cascade angular difference to track
		// The following code is mostly cribbed from the half-sphere
		// parametrization
		I3Particle &gradpart = (*gradsources)[0];

		double p0 = par_[par_.size() - 2];
		double p1 = par_[par_.size() - 1];
		double seedX_, seedY_, seedZ_;
		double perp1X_, perp1Y_, perp1Z_;
		double perp2X_, perp2Y_, perp2Z_;

		const I3Direction& dir = hypothesis_->particle->GetDir();
		const I3Direction& graddir = gradpart.GetDir();

		seedX_ = dir.GetX();
		seedY_ = dir.GetY();
		seedZ_ = dir.GetZ();
		std::pair<I3Direction,I3Direction> sideways =
		    I3Calculator::GetTransverseDirections(dir);
		const I3Direction &perp1 = sideways.first;
		const I3Direction &perp2 = sideways.second;
		perp1X_ = perp1.GetX();
		perp1Y_ = perp1.GetY();
		perp1Z_ = perp1.GetZ();
		perp2X_ = perp2.GetX();
		perp2Y_ = perp2.GetY();
		perp2Z_ = perp2.GetZ();
		
		double newdirx = seedX_ + p0 * perp1X_ + p1 * perp2X_;
		double newdiry = seedY_ + p0 * perp1Y_ + p1 * perp2Y_;
		double newdirz = seedZ_ + p0 * perp1Z_ + p1 * perp2Z_;
		double r = sqrt(newdirx*newdirx + newdiry*newdiry +
		    newdirz*newdirz);
		par_gradient_[par_gradient_.size()-2] = graddir.GetZenith()*
		    (-newdirz*(newdirx*perp1X_ + newdiry*perp1Y_) +
		    perp1Z_*(newdirx*newdirx + newdiry*newdiry))/
		    (hypot(newdirx, newdiry)*r*r) +
		    graddir.GetAzimuth()*(perp1Y_*newdirx - perp1X_*newdiry)/
		    (newdirx*newdirx + newdiry*newdiry);
		par_gradient_[par_gradient_.size()-1] = graddir.GetZenith()*
		    (-newdirz*(newdirx*perp2X_ + newdiry*perp2Y_) +
		    perp2Z_*(newdirx*newdirx + newdiry*newdiry))/
		    (hypot(newdirx, newdiry)*r*r) +
		    graddir.GetAzimuth()*(perp2Y_*newdirx - perp2X_*newdiry)/
		    (newdirx*newdirx + newdiry*newdiry);
	}
}

void
BipedParametrization::Configure()
{
	I3SimpleParametrization::Configure();

	GetParameter("Boundary", boundary_);
	GetParameter("MuonSpacing", muonspacing_);
	GetParameter("ShowerSpacing", showerspacing_);
	GetParameter("StartingCascadeStepSize", starting_cascade_dirstep_);
	GetParameter("StartSlantDepth", slantstart_);
	GetParameter("EndSlantDepth", slantstop_);

	if (starting_cascade_dirstep_ > 0) {
		I3FitParameterInitSpecs specs("dir");
		specs.minval_ = 0;
		specs.maxval_ = 0;
		specs.stepsize_ = starting_cascade_dirstep_;

		specs.name_ = "cascadedir1";
		parspecs_.push_back(specs);
		
		specs.name_ = "cascadedir2";
		parspecs_.push_back(specs);

		par_.resize(parspecs_.size());
	}
}

static void
BipedHypothesis(I3ParticleConstPtr track,
    std::vector<I3Particle> &hypothesis, double boundary, double muonspacing,
    double showerspacing, double slantstart, double slantstop)
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

	if (track->GetShape() == I3Particle::StartingTrack) {
		entry_t = track->GetTime();
	} else if (track->GetShape() == I3Particle::StoppingTrack) {
		exit_t = track->GetTime();
	} else if (track->GetShape() == I3Particle::ContainedTrack) {
		entry_t = track->GetTime();
		exit_t = track->GetTime() +
		    track->GetLength()/track->GetSpeed();
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
}

