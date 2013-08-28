#include <millipede/Millipede.h>

#include <phys-services/I3Calculator.h>
#include <lilliput/parametrization/I3SimpleParametrization.h>
#include <gulliver/I3ParametrizationBase.h>
#include <gulliver/I3EventHypothesis.h>
#include <icetray/I3ServiceFactory.h>
#include <icetray/I3SingleServiceFactory.h>

static void BipedHypothesis(I3ParticleConstPtr track,
    std::vector<I3Particle> &hypothesis, double boundary,
    double muonspacing);



class BipedParametrization : public I3SimpleParametrization {
	public:
		BipedParametrization(const I3Context &);
		void Configure();

		double boundary_;
		double muonspacing_;
		std::string name_;
		double starting_cascade_dirstep_;

		void UpdatePhysicsVariables();
		void UpdateParameters();
		bool InitChainRule(bool grad); 
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
	    1*I3Units::m);
	AddParameter("StartingCascadeStepSize", "Allow the initial cascade to "
	    "point in a different direction than the track, moving it in steps "
	    "this size (starting direction fixed if 0).", 0);
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
	    boundary_, muonspacing_);
	//log_info("%f is the muon spacing in Para", muonspacing_); 

	if (starting_cascade_dirstep_ > 0 && sources->size() > 0) {
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
		log_info("We have initialized Gradients in Para");
		boost::shared_ptr<I3Vector<I3Particle> > gradsources(new
		    I3Vector<I3Particle>);
		gradsources->resize(sources->size());
		for (unsigned i = 0; i < gradsources->size(); ++i) { 
	                        I3Particle &gradpart = (*gradsources)[i]; 
	                        gradpart.SetPos(0., 0., 0.); 
	                        gradpart.SetDir(0., 0.); 
	                        gradpart.SetTime(0.); 
	                        gradpart.SetEnergy(0.); 
	                        gradpart.SetLength(0.); 
	                        gradpart.SetSpeed(0.); 
	                }
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

bool 
BipedParametrization::InitChainRule(bool want_grad) 
{ 
	I3SimpleParametrization::InitChainRule(want_grad); 

	if (want_grad && starting_cascade_dirstep_ > 0) 
		par_gradient_.resize(parspecs_.size()); 

	return true; 
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
		azimuth_xyz += gradpart.GetY()*sin(theta)*cos(phi);
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
		//log_info("StartingCascadeDirStepPositive");
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
	GetParameter("StartingCascadeStepSize", starting_cascade_dirstep_);

	if (starting_cascade_dirstep_ > 0) {
		//log_info("OneMoreTimePositiveCascadeDirStep");
		I3FitParameterInitSpecs specs("dir");
		specs.minval_ = 0;
		specs.maxval_ = 0;
		specs.stepsize_ = starting_cascade_dirstep_;

		specs.name_ = "cascadedir1";
		//log_info("What is CascadeDir1?");
		parspecs_.push_back(specs);
		
		specs.name_ = "cascadedir2";
		parspecs_.push_back(specs);
		//log_info("PushingBackTheSpecsForCascadeDir2");

		par_.resize(parspecs_.size());
		//log_info("Hath Thou Resized the Parspecs?");
	}
		//log_info("%f is the muon spacing AGAIN o.O", muonspacing_); 
}

static void
BipedHypothesis(I3ParticleConstPtr track,
    std::vector<I3Particle> &hypothesis, double boundary, double muonspacing)
{
	I3Particle muon, cascade;
	muon.SetType(I3Particle::MuMinus);
	muon.SetShape(I3Particle::ContainedTrack);
	muon.SetLength(track->GetLength());
	muon.SetDir(track->GetDir());
	muon.SetPos(track->GetX(),track->GetY(),track->GetZ());
	muon.SetTime(track->GetTime());
	muon.SetSpeed(track->GetSpeed());
	muon.SetEnergy(track->GetEnergy()*9/10);
//	muon.SetEnergy(track->GetLength()/5.00);
	cascade.SetType(I3Particle::Hadrons);
	cascade.SetShape(I3Particle::Cascade);
	cascade.SetDir(track->GetDir());
	cascade.SetPos(track->GetX(),track->GetY(),track->GetZ());
	cascade.SetTime(track->GetTime());
	cascade.SetSpeed(track->GetSpeed());
	cascade.SetEnergy(track->GetEnergy()/10);
//	cascade.SetEnergy(track->GetLength()/50.00);
	hypothesis.push_back(cascade);
	hypothesis.push_back(muon);
}

