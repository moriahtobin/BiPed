#include <millipede/Millipede.h>

#include <gulliver/I3EventLogLikelihoodBase.h>
#include <gulliver/I3EventHypothesis.h>
#include <icetray/I3SingleServiceFactory.h>

#include <boost/make_shared.hpp>

typedef I3Vector<I3Particle> I3VectorI3Particle;
typedef boost::shared_ptr<I3VectorI3Particle> I3VectorI3ParticlePtr;

class BipedLikelihood : public I3MillipedeService,
    public I3EventLogLikelihoodBase {
	public:
		BipedLikelihood(const I3Context &context);

		void Configure();
		double muonspacing_;
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
		SET_LOGGER("BipedLikelihood");

		
		I3VectorI3ParticlePtr ExtractHypothesis(const I3EventHypothesis&);
};

I3_POINTER_TYPEDEFS(BipedLikelihood);

BipedLikelihood::BipedLikelihood(const I3Context &context) : I3MillipedeService(context), I3EventLogLikelihoodBase()
{
	AddParameter("MuonSpacing", "Spacing of muon (ionization) sources "
	    "along the track, in meters. MUST match source extension in "
	    "muon tables. Set to 0 to not include ionization sources.",
	    1*I3Units::m);
}

unsigned int
BipedLikelihood::GetMultiplicity()
{
	unsigned j = 0;

	for (MillipedeDOMCacheMap::const_iterator i = domCache_.begin();
	    i != domCache_.end(); i++)
		j += i->second.valid.count();

	return j;
}

void
BipedLikelihood::SetEvent(const I3Frame &frame)
{
	DatamapFromFrame(frame);
}

double
BipedLikelihood::GetLogLikelihood(const I3EventHypothesis &hypo)
{
	return GetLogLikelihood(hypo, NULL, false, 1.);
}

double
BipedLikelihood::GetLogLikelihoodWithGradient(const I3EventHypothesis &hypo,
    I3EventHypothesis &gradient, double weight)
{
	return GetLogLikelihood(hypo, &gradient, true, 1.);
}

I3VectorI3ParticlePtr
BipedLikelihood::ExtractHypothesis(const I3EventHypothesis &hypo)
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
BipedLikelihood::GetDiagnostics(const I3EventHypothesis &hypo)
{
	cholmod_sparse *response_matrix, *many_response_matrix;
	I3VectorI3ParticlePtr sources = ExtractHypothesis(hypo);

	//most of what follows is a compact form of many of the things that happen (with narrative)
	//in the GetLogLikelihood function further down so look in GetLogLikelihood for running commentary
	double trackLength = (*sources)[1].GetLength();
	unsigned int MuSeg = ceil(trackLength / muonspacing_);
	boost::shared_ptr<I3Vector<I3Particle> > microSources(new I3Vector<I3Particle>);
	microSources->push_back((*sources)[0]);
	I3Particle PrevParticle = (*sources)[1];
	I3Particle NextParticle = PrevParticle;
	for (double d = 0; d < MuSeg*muonspacing_; d = d+muonspacing_){
		microSources->push_back(PrevParticle);
		NextParticle.SetPos(PrevParticle.ShiftAlongTrack(muonspacing_));
		NextParticle.SetTime(PrevParticle.GetTime()+muonspacing_/I3Constants::c);
		PrevParticle = NextParticle;
	}
	double lastSeg = trackLength/muonspacing_ - (MuSeg-1);
	double MuEnFact = muonspacing_/trackLength;
	double endscale = MuEnFact*lastSeg;
	cholmod_triplet *pen_trip = cholmod_l_allocate_triplet(
	    (*microSources).size(), (*sources).size(), (*microSources).size(), 0,
	    CHOLMOD_REAL, &c);
	pen_trip->nnz = 0;
	unsigned i = 0; 
		((long *)(pen_trip->i))[pen_trip->nnz] = i;
		((long *)(pen_trip->j))[pen_trip->nnz] = i;
		((double *)(pen_trip->x))[pen_trip->nnz] = 1;
		pen_trip->nnz++;
	for (unsigned i = 1; i < pen_trip->nrow-1; i++) {
		((long *)(pen_trip->i))[pen_trip->nnz] = i;
		((long *)(pen_trip->j))[pen_trip->nnz] = (i != 0);
		((double *)(pen_trip->x))[pen_trip->nnz] = MuEnFact;
		pen_trip->nnz++;
	}
		((long *)(pen_trip->i))[pen_trip->nnz] = i;
		((long *)(pen_trip->j))[pen_trip->nnz] = pen_trip->nrow-1;
		((double *)(pen_trip->x))[pen_trip->nnz] = endscale;
	cholmod_sparse *collapser = cholmod_l_triplet_to_sparse(pen_trip, 0, &c);
	cholmod_l_free_triplet(&pen_trip, &c);
	many_response_matrix = Millipede::GetResponseMatrix(domCache_, *microSources,
	    domEfficiency_, muon_p, cascade_p, NULL, &c);
	if (many_response_matrix == NULL)
		log_fatal("Null basis matrix");
	response_matrix = cholmod_l_ssmult(many_response_matrix, 
	    collapser, 0, 1, 0, &c);
	cholmod_l_free_sparse(&collapser, &c);
	cholmod_l_free_sparse(&many_response_matrix, &c);
	if (response_matrix == NULL)
		log_fatal("Null basis matrix (little response)");
	MillipedeFitParamsPtr params =
	    boost::make_shared<MillipedeFitParams>();
	Millipede::FitStatistics(domCache_, *sources, I3Units::MeV,
		    response_matrix, params.get(), &c);
	log_info("created Fit Statistics");
	cholmod_l_free_sparse(&response_matrix, &c);

	return params;
}

void
BipedLikelihood::Configure()
{
	I3MillipedeService::Configure();

	GetParameter("MuonSpacing", muonspacing_);
	//log_info("%f is the muon spacing in LikConfig", muonspacing_); 
}

double
BipedLikelihood::GetLogLikelihood(const I3EventHypothesis &hypo,
    I3EventHypothesis *gradient, bool fit_energy, double weight)
{
	cholmod_sparse *response_matrix, *many_response_matrix, *gradients, *manygradients;
	double llh = 0.;

	I3VectorI3ParticlePtr sources = ExtractHypothesis(hypo);


	I3VectorI3ParticlePtr gradsources;
	if (gradient != NULL) {
		gradsources = ExtractHypothesis(*gradient);
		assert(sources->size() == gradsources->size());
	}

	// Shortcut to noise-only likelihood if no elements
	if (sources->size() == 0)
		return Millipede::FitStatistics(domCache_, *sources,
		    I3Units::MeV, NULL, NULL, &c);

// Make a vector of particles out of the cascade & muon:
	// cascade is the starting cascade, and muon is replaced
	// by equal energy, closely spaced cascades to simulate
	// min ionizing. 
	
	// Start with cascade as opening particle:
	boost::shared_ptr<I3Vector<I3Particle> > microSources(new I3Vector<I3Particle>);
	microSources->push_back((*sources)[0]);


	// get total track length from hypothesis muon
	//collect diagnostic variables
	double trackLength = (*sources)[1].GetLength();
	double zen_cascade = (*sources)[0].GetZenith();
	double zen_track = (*sources)[1].GetZenith();
	double azi_track = (*sources)[1].GetAzimuth();
	double azi_casc = (*sources)[0].GetAzimuth();
	double x_v = (*sources)[0].GetX();
	double y_v = (*sources)[0].GetY();
	double z_v = (*sources)[0].GetZ();
	double MuSeg = ceil(trackLength / muonspacing_);
	log_info("%f is the number of ceil(L/muonspacing)", MuSeg);

	//create segmented muon by grabbing original muon
	I3Particle PrevParticle = (*sources)[1];
	I3Particle NextParticle = PrevParticle;

	//create muon segments using user specified muonspacing
	for (double d = 0; d < MuSeg*muonspacing_; d = d+muonspacing_){
		microSources->push_back(PrevParticle); 
		NextParticle.SetPos(PrevParticle.ShiftAlongTrack(muonspacing_));
		NextParticle.SetTime(PrevParticle.GetTime()+muonspacing_/I3Constants::c);
		PrevParticle = NextParticle;
	}

	//Create energy scaling factor for muon segments, since the 	
	//energy solver assumes it is dealing with the full particle energy
	//Further, the last muon segment is not a full segment so here we
	//create a energy scaling factor for the last segment (endscale)
	double lastSeg = trackLength/muonspacing_ - (MuSeg - 1.0);		
	double MuEnFact = muonspacing_/trackLength;
	double endscale = MuEnFact*lastSeg;


	// Make a matrix of ones and zeros for collapsing the multi-muon response matrix
	// into a single muon response matrix
	// (a triplet is: matrix position i,j; and value x)
	cholmod_triplet *pen_trip = cholmod_l_allocate_triplet(
	    (*microSources).size(), (*sources).size(), (*microSources).size(), 0,
	    CHOLMOD_REAL, &c);
	// starting number of non-zero entries in this matrix:
	pen_trip->nnz = 0;

	unsigned i = 0; 
		((long *)(pen_trip->i))[pen_trip->nnz] = i;
		((long *)(pen_trip->j))[pen_trip->nnz] = i;
		((double *)(pen_trip->x))[pen_trip->nnz] = 1;
		pen_trip->nnz++;
	//Energy Correction for Composite Muon
	for (unsigned i = 1; i < pen_trip->nrow-1; i++) {
		((long *)(pen_trip->i))[pen_trip->nnz] = i;
		((long *)(pen_trip->j))[pen_trip->nnz] = (i != 0);
		((double *)(pen_trip->x))[pen_trip->nnz] = MuEnFact;
		pen_trip->nnz++;
	}
	((long *)(pen_trip->i))[pen_trip->nnz] = pen_trip->nrow-1;
	((long *)(pen_trip->j))[pen_trip->nnz] = 1;
	((double *)(pen_trip->x))[pen_trip->nnz]=endscale;
	pen_trip->nnz++;
	cholmod_sparse *collapser = cholmod_l_triplet_to_sparse(pen_trip, (*microSources).size(), &c);
	cholmod_l_free_triplet(&pen_trip, &c);
//	If you want to print out the matrix to check out,
//	you have to convert it into a dense matrix as follows:
//	cholmod_dense *one_muon_mult =
//	cholmod_l_sparse_to_dense(collapser, &c);
//	for(unsigned row = 0; row < collapser->nrow; row++){
//		for(unsigned col = 0; col < collapser->ncol; col++){
//		double elem = ((double*)(one_muon_mult->x))
//			[row + col*one_muon_mult->nrow];
//		std::cout<<row<<","<<col<<" : "<<elem<<std::endl;
//		}
//	}

	//Now we creat a matrix which grabs out the gradient for the length parameter from the response matrix
	cholmod_sparse *length_grad_part = cholmod_l_allocate_sparse(
		(*microSources).size(), 14u, 1u, true, true, 0, CHOLMOD_REAL, &c);
	        long *cols = (long*)(length_grad_part->p);
	        long *rows = (long*)(length_grad_part->i);
	        double *x = (double*)(length_grad_part->x);
	        std::fill(&cols[0], &cols[14], 0);
	        // The last column has one entry
	        cols[14] = 1;
	        // in the last row
	        rows[0] = (*microSources).size()-1;
	        // dP[nsources-1,1]/dL
	        x[0] = MuEnFact/muonspacing_;


	many_response_matrix = Millipede::GetResponseMatrix(domCache_, *microSources,
	    domEfficiency_, muon_p, cascade_p,
	    (gradient == NULL) ? NULL : &manygradients, &c);
	if (many_response_matrix == NULL){
		log_fatal("Null basis matrix");
	}

	// Collapse the resulting big response matrix (with all muon segments) back into a 2-particle matrix
	//AKA: We're just taking the light that was generated by all the muon segments and attributing them to the same composite muon
	response_matrix = cholmod_l_ssmult(many_response_matrix, 
	    collapser, 0, 1, 0, &c);
	cholmod_l_free_sparse(&collapser, &c);
	//Get the length gradient
	cholmod_sparse *length_grad = cholmod_l_ssmult(many_response_matrix, 
	    length_grad_part, 0, 1, 0, &c);
	cholmod_l_free_sparse(&length_grad_part, &c);
	cholmod_l_free_sparse(&many_response_matrix, &c);

	// Make matrix for converting multi-muon gradients to single muon little_gradients
	cholmod_triplet *grad_trip = cholmod_l_allocate_triplet(
	    (*microSources).size()*7, 14u, (*microSources).size()*12, 0,
	    CHOLMOD_REAL, &c);
	// starting number of non-zero entries in this matrix:
	grad_trip->nnz = 0;
	//Keep the same gradients for the cascade
	for (unsigned i = 0; i < 6; i++) { 
		((long *)(grad_trip->i))[grad_trip->nnz] = i;
		((long *)(grad_trip->j))[grad_trip->nnz] = i;
		((double *)(grad_trip->x))[grad_trip->nnz] = 1;
		grad_trip->nnz++;
	}
	//Combine muon gradients
	//x, y, z, t
	for (unsigned i = 7; i < grad_trip->nrow-7; i+=7) {
		((long *)(grad_trip->i))[grad_trip->nnz] = i;
		((long *)(grad_trip->j))[grad_trip->nnz] = 7;
		((double *)(grad_trip->x))[grad_trip->nnz] = MuEnFact;
		grad_trip->nnz++;
	}
	((long *)(grad_trip->i))[grad_trip->nnz] = grad_trip->nrow-7;
	((long *)(grad_trip->j))[grad_trip->nnz] = 7;
	((double *)(grad_trip->x))[grad_trip->nnz] = endscale;
	grad_trip->nnz++;
	for (unsigned i = 8; i < grad_trip->nrow-7; i+=7) {
		((long *)(grad_trip->i))[grad_trip->nnz] = i;
		((long *)(grad_trip->j))[grad_trip->nnz] = 8;
		((double *)(grad_trip->x))[grad_trip->nnz] = MuEnFact;
		grad_trip->nnz++;
	}
	((long *)(grad_trip->i))[grad_trip->nnz] = grad_trip->nrow-6;
	((long *)(grad_trip->j))[grad_trip->nnz] = 8;
	((double *)(grad_trip->x))[grad_trip->nnz] = endscale;
	grad_trip->nnz++;
	for (unsigned i = 9; i < grad_trip->nrow-7; i+=7) {
		((long *)(grad_trip->i))[grad_trip->nnz] = i;
		((long *)(grad_trip->j))[grad_trip->nnz] = 9;
		((double *)(grad_trip->x))[grad_trip->nnz] = MuEnFact;
		grad_trip->nnz++;
	}
	((long *)(grad_trip->i))[grad_trip->nnz] = grad_trip->nrow-5;
	((long *)(grad_trip->j))[grad_trip->nnz] = 9;
	((double *)(grad_trip->x))[grad_trip->nnz] = endscale;
	grad_trip->nnz++;
	for (unsigned i = 10; i < grad_trip->nrow-7; i+=7) {
		((long *)(grad_trip->i))[grad_trip->nnz] = i;
		((long *)(grad_trip->j))[grad_trip->nnz] = 10;
		((double *)(grad_trip->x))[grad_trip->nnz] = MuEnFact;
		grad_trip->nnz++;
	}
	((long *)(grad_trip->i))[grad_trip->nnz] = grad_trip->nrow-4;
	((long *)(grad_trip->j))[grad_trip->nnz] = 10;
	((double *)(grad_trip->x))[grad_trip->nnz] = endscale;
	grad_trip->nnz++;
	//zenith	
	//includes lever arm effect
        std::vector<double> muonLengths((*microSources).size(), 0.);
                for (unsigned i = 1; i <(*microSources).size()-1; i++){
                	muonLengths[i] = MuEnFact*i*muonspacing_;
		}
		muonLengths.back()= endscale*trackLength;
	double zen = (*microSources)[1].GetZenith();
	double azi = (*microSources)[1].GetAzimuth();
	unsigned nom = 1;
	for (unsigned i = 11; i < grad_trip->nrow-7; i+=7) {
		((long *)(grad_trip->i))[grad_trip->nnz] = i;
		((long *)(grad_trip->j))[grad_trip->nnz] = 11;
		((double *)(grad_trip->x))[grad_trip->nnz] = MuEnFact;
		grad_trip->nnz++;
		((long *)(grad_trip->i))[grad_trip->nnz] = i-4;
		((long *)(grad_trip->j))[grad_trip->nnz] = 11;
		((double *)(grad_trip->x))[grad_trip->nnz] = -std::cos(zen)*std::cos(azi)*muonLengths[nom];
		grad_trip->nnz++;
		((long *)(grad_trip->i))[grad_trip->nnz] = i-3;
		((long *)(grad_trip->j))[grad_trip->nnz] = 11;
		((double *)(grad_trip->x))[grad_trip->nnz] = -std::cos(zen)*std::sin(azi)*muonLengths[nom];
		grad_trip->nnz++;
		((long *)(grad_trip->i))[grad_trip->nnz] = i-2;
		((long *)(grad_trip->j))[grad_trip->nnz] = 11;
		((double *)(grad_trip->x))[grad_trip->nnz] = std::sin(zen)*muonLengths[nom];
		grad_trip->nnz++;
		nom++;
	}
	((long *)(grad_trip->i))[grad_trip->nnz] = grad_trip->nrow-3;
	((long *)(grad_trip->j))[grad_trip->nnz] = 11;
	((double *)(grad_trip->x))[grad_trip->nnz] = endscale;
	grad_trip->nnz++;
	((long *)(grad_trip->i))[grad_trip->nnz] = grad_trip->nrow-7;
	((long *)(grad_trip->j))[grad_trip->nnz] = 11;
	((double *)(grad_trip->x))[grad_trip->nnz] = -std::cos(zen)*std::cos(azi)*muonLengths.back();
	grad_trip->nnz++;
	((long *)(grad_trip->i))[grad_trip->nnz] = grad_trip->nrow-6;
	((long *)(grad_trip->j))[grad_trip->nnz] = 11;
	((double *)(grad_trip->x))[grad_trip->nnz] = -std::cos(zen)*std::sin(azi)*muonLengths.back();
	grad_trip->nnz++;
	((long *)(grad_trip->i))[grad_trip->nnz] = grad_trip->nrow-5;
	((long *)(grad_trip->j))[grad_trip->nnz] = 11;
	((double *)(grad_trip->x))[grad_trip->nnz] = std::sin(zen)*muonLengths.back();
	grad_trip->nnz++;
	//azimuth
	//includes lever arm effect
	nom = 1;
	for (unsigned i = 12; i < grad_trip->nrow-7; i+=7) {
		((long *)(grad_trip->i))[grad_trip->nnz] = i;
		((long *)(grad_trip->j))[grad_trip->nnz] = 12;
		((double *)(grad_trip->x))[grad_trip->nnz] = MuEnFact;
		grad_trip->nnz++;
		((long *)(grad_trip->i))[grad_trip->nnz] = i-5;
		((long *)(grad_trip->j))[grad_trip->nnz] = 12;
		((double *)(grad_trip->x))[grad_trip->nnz] = std::sin(zen)*std::sin(azi)*muonLengths[nom];
		grad_trip->nnz++;
		((long *)(grad_trip->i))[grad_trip->nnz] = i-4;
		((long *)(grad_trip->j))[grad_trip->nnz] = 12;
		((double *)(grad_trip->x))[grad_trip->nnz] = -std::sin(zen)*std::cos(azi)*muonLengths[nom];
		grad_trip->nnz++;
		nom++;
	}
	((long *)(grad_trip->i))[grad_trip->nnz] = grad_trip->nrow-2;
	((long *)(grad_trip->j))[grad_trip->nnz] = 12;
	((double *)(grad_trip->x))[grad_trip->nnz] = endscale;
	grad_trip->nnz++;
	((long *)(grad_trip->i))[grad_trip->nnz] = grad_trip->nrow-7;
	((long *)(grad_trip->j))[grad_trip->nnz] = 12;
	((double *)(grad_trip->x))[grad_trip->nnz] = std::sin(zen)*std::sin(azi)*muonLengths.back();
	grad_trip->nnz++;
	((long *)(grad_trip->i))[grad_trip->nnz] = grad_trip->nrow-6;
	((long *)(grad_trip->j))[grad_trip->nnz] = 12;
	((double *)(grad_trip->x))[grad_trip->nnz] = -std::sin(zen)*std::cos(azi)*muonLengths.back();
	//Convert to sparse matrix for type matching
	cholmod_sparse *grad_collapser = cholmod_l_triplet_to_sparse(grad_trip, (*microSources).size()*12, &c);
	cholmod_l_free_triplet(&grad_trip, &c);

	//use the gradient collapser matrix to create our 2-particle spatial gradients matrix from the muon segement many gradients matrix
	cholmod_sparse *spatial_grad = cholmod_l_ssmult(manygradients, 
	    grad_collapser, 0, 1, 0, &c);
	
	//Combine spatial and length gradient matrices to produce the full gradients matrix
	double alpha[2] = {1,0}, beta[2] = {1,0};
	gradients = cholmod_l_add(spatial_grad, length_grad, alpha, beta, true, true, &c);
	cholmod_l_free_sparse(&spatial_grad, &c);
	cholmod_l_free_sparse(&length_grad, &c);
	cholmod_l_free_sparse(&grad_collapser, &c);


	//Now fit for the energy using our collapsed gradients and response matrix
	if (fit_energy) {
		SolveEnergyLosses(*sources, response_matrix,
		    (gradient == NULL ? NULL : gradients));
		if (sources->size() == 1)
			hypo.particle->SetEnergy((*sources)[0].GetEnergy());
	}


	llh = Millipede::FitStatistics(domCache_, *sources, I3Units::MeV,
	    response_matrix, NULL, &c);
		//check that the fitting is sensible by printing out the llh value for each set of checked parameters
		log_info("[%f m mu, vertex (%f, %f, %f)] + ", trackLength, x_v, y_v, z_v);
		log_info("[ (%f zen_mu, %f zen_casc), (%f azi_mu, %f azi_casc)] -> (llh=%f)", zen_track, zen_cascade, azi_track, azi_casc, llh);

	if (gradient != NULL) {
		log_info("Begin non-null gradients loop");
//		assert(sources->size() == gradsources->size());
		Millipede::LLHGradient(domCache_, *sources, *gradsources,
		    I3Units::MeV, weight, response_matrix, gradients, &c);
		cholmod_l_free_sparse(&manygradients, &c);
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
		for(int i=0; i<gradsources->size(); i++) log_info_stream((*gradsources)[i]);
	}
	if (gradient == NULL) {
	log_debug("NULL GRADIENT MATRIX");
	}

	cholmod_l_free_sparse(&response_matrix, &c);
	return llh;
}

typedef I3SingleServiceFactory<BipedLikelihood, I3EventLogLikelihoodBase>
    BipedLikelihoodFactory;
I3_SERVICE_FACTORY(BipedLikelihoodFactory);
