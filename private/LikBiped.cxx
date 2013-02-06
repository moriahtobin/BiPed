#include <millipede/Millipede.h>

#include <gulliver/I3EventLogLikelihoodBase.h>
#include <gulliver/I3EventHypothesis.h>
#include <icetray/I3SingleServiceFactory.h>

#include <boost/make_shared.hpp>

typedef I3Vector<I3Particle> I3VectorI3Particle;
typedef boost::shared_ptr<I3VectorI3Particle> I3VectorI3ParticlePtr;

class MillipedeLikelihood : public I3MillipedeService,
    public I3EventLogLikelihoodBase {
	public:
		MillipedeLikelihood(const I3Context &context);

		// Gulliver interface
		void SetEvent(const I3Frame &);
                void SetGeometry(const I3Geometry &geo) {};
		bool HasGradient() { return true; }
		double GetLogLikelihood(const I3EventHypothesis &);
		double GetLogLikelihoodWithGradient(const I3EventHypothesis &,
		    I3EventHypothesis &gradient, double weight);
		unsigned int GetMultiplicity();
		
		double GetLogLikelihood(const I3EventHypothesis&,
		    I3EventHypothesis*, bool solve_for_energies, double weight);
		
		I3FrameObjectPtr GetDiagnostics(const I3EventHypothesis &);
				
		virtual const std::string GetName() const {
			return I3ServiceBase::GetName();
		}
			
	private:
		SET_LOGGER("MillipedeLikelihood");
		
		I3VectorI3ParticlePtr ExtractHypothesis(const I3EventHypothesis&);
};

I3_POINTER_TYPEDEFS(MillipedeLikelihood);

MillipedeLikelihood::MillipedeLikelihood(const I3Context &context) : I3MillipedeService(context), I3EventLogLikelihoodBase()
{}

unsigned int
MillipedeLikelihood::GetMultiplicity()
{
	unsigned j = 0;

	for (MillipedeDOMCacheMap::const_iterator i = domCache_.begin();
	    i != domCache_.end(); i++)
		j += i->second.valid.count();

	return j;
}

void
MillipedeLikelihood::SetEvent(const I3Frame &frame)
{
	DatamapFromFrame(frame);
}

double
MillipedeLikelihood::GetLogLikelihood(const I3EventHypothesis &hypo)
{
	return GetLogLikelihood(hypo, NULL, false, 1.);
}

double
MillipedeLikelihood::GetLogLikelihoodWithGradient(const I3EventHypothesis &hypo,
    I3EventHypothesis &gradient, double weight)
{
	return GetLogLikelihood(hypo, &gradient, true, 1.);
}

I3VectorI3ParticlePtr
MillipedeLikelihood::ExtractHypothesis(const I3EventHypothesis &hypo)
{
	I3VectorI3ParticlePtr sources = 
	    boost::dynamic_pointer_cast<I3VectorI3Particle>(hypo.nonstd);
	if (sources == NULL) {
		sources = boost::make_shared<I3VectorI3Particle>();
		sources->push_back(*hypo.particle);
	}
	return sources;
}

I3FrameObjectPtr
MillipedeLikelihood::GetDiagnostics(const I3EventHypothesis &hypo)
{
	cholmod_sparse *response_matrix;
	I3VectorI3ParticlePtr sources = ExtractHypothesis(hypo);
	
	response_matrix = Millipede::GetResponseMatrix(domCache_, *sources,
	    domEfficiency_, muon_p, cascade_p, NULL, &c);
	if (response_matrix == NULL)
		log_fatal("Null basis matrix");
	
	MillipedeFitParamsPtr params =
	    boost::make_shared<MillipedeFitParams>();
	Millipede::FitStatistics(domCache_, *sources, I3Units::MeV,
		    response_matrix, params.get(), &c);
	cholmod_l_free_sparse(&response_matrix, &c);
	
	return params;
}

double
MillipedeLikelihood::GetLogLikelihood(const I3EventHypothesis &hypo,
    I3EventHypothesis *gradient, bool fit_energy, double weight)
{
	cholmod_sparse *response_matrix, *gradients;
	double llh = 0.;

	I3VectorI3ParticlePtr sources = ExtractHypothesis(hypo);

	I3VectorI3ParticlePtr gradsources;
	if (gradient != NULL) {
		gradsources = ExtractHypothesis(*gradient);
		assert(sources->size() == gradsources->size());
	}
	
	response_matrix = Millipede::GetResponseMatrix(domCache_, *sources,
	    domEfficiency_, muon_p, cascade_p,
	    (gradient == NULL) ? NULL : &gradients, &c);
	if (response_matrix == NULL)
		log_fatal("Null basis matrix");
	
	if (fit_energy) {
		SolveEnergyLosses(*sources, response_matrix,
		    (gradient == NULL ? NULL : gradients));
		if (sources->size() == 1)
			hypo.particle->SetEnergy((*sources)[0].GetEnergy());
	}
	
	llh = Millipede::FitStatistics(domCache_, *sources, I3Units::MeV,
	    response_matrix, NULL, &c);

	if (gradient != NULL) {
		Millipede::LLHGradient(domCache_, *sources, *gradsources,
		    I3Units::MeV, weight, response_matrix, gradients, &c);
		cholmod_l_free_sparse(&gradients, &c);
		if (sources->size() == 1) {
			// NB: the requested weight is already applied in
			// the call to LLHGradient()
			I3Particle &p = *gradient->particle;
			p.SetPos(I3Position(
			    p.GetPos().GetX() +
			      (*gradsources)[0].GetPos().GetX(),
			    p.GetPos().GetY() +
			      (*gradsources)[0].GetPos().GetY(),
			    p.GetPos().GetZ() +
			      (*gradsources)[0].GetPos().GetZ()
			));
			p.SetTime(p.GetTime() + (*gradsources)[0].GetTime());
			p.SetDir(I3Direction(
			    p.GetDir().GetZenith() +
			      (*gradsources)[0].GetDir().GetZenith(),
			    p.GetDir().GetAzimuth() +
			      (*gradsources)[0].GetDir().GetAzimuth()
			));			
		}
	}

	cholmod_l_free_sparse(&response_matrix, &c);
	
	return llh;
}

typedef I3SingleServiceFactory<MillipedeLikelihood, I3EventLogLikelihoodBase>
    MillipedeLikelihoodFactory;
I3_SERVICE_FACTORY(MillipedeLikelihoodFactory);

