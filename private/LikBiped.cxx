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
	log_info("InGetLLHLoop");
	return GetLogLikelihood(hypo, NULL, false, 1.);
}

double
BipedLikelihood::GetLogLikelihoodWithGradient(const I3EventHypothesis &hypo,
    I3EventHypothesis &gradient, double weight)
{
	log_info("InGetLLHwGradsLoop");
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
	log_info("We Have Sources!");
}

I3FrameObjectPtr
BipedLikelihood::GetDiagnostics(const I3EventHypothesis &hypo)
{
	cholmod_sparse *response_matrix, *little_response_matrix;
	I3VectorI3ParticlePtr sources = ExtractHypothesis(hypo);
	
	log_info("Diagnostics Function");
	double trackLength = (*sources)[1].GetLength();
	unsigned int check2 = trackLength / muonspacing_;
	boost::shared_ptr<I3Vector<I3Particle> > microSources(new I3Vector<I3Particle>);
	microSources->push_back((*sources)[0]);
	I3Particle PrevParticle = (*sources)[1];
	I3Particle NextParticle = PrevParticle;
	if (trackLength >= 1.5*muonspacing_){
	microSources->push_back(PrevParticle);
		unsigned int checkend = abs(check2-2);
		double endspace = (trackLength-(checkend)*muonspacing_)/2;
		for (double d = 0; d < checkend; d = d+muonspacing_){
			NextParticle.SetPos(PrevParticle.ShiftAlongTrack(muonspacing_));
			PrevParticle = NextParticle;
			microSources->push_back(PrevParticle);
		}
		NextParticle.SetPos(PrevParticle.ShiftAlongTrack(endspace));
	}
	else {
	microSources->push_back(PrevParticle);
//		(*sources).erase((*sources).end());
	}
	double MuSeg = (*microSources).size() -1.0;
	double MuEnFact = 1.0;
	if(MuSeg > 0){
		MuEnFact = 1.0/MuSeg;
	}
	cholmod_triplet *pen_trip = cholmod_l_allocate_triplet(
	    (*microSources).size(), (*sources).size(), (*microSources).size(), 0,
	    CHOLMOD_REAL, &c);
	pen_trip->nnz = 0;
	unsigned i = 0; 
		((long *)(pen_trip->i))[pen_trip->nnz] = i;
		((long *)(pen_trip->j))[pen_trip->nnz] = (i != 0);
		((double *)(pen_trip->x))[pen_trip->nnz] = 1;
		pen_trip->nnz++;
	for (unsigned i = 1; i < pen_trip->nrow; i++) {
		((long *)(pen_trip->i))[pen_trip->nnz] = i;
		((long *)(pen_trip->j))[pen_trip->nnz] = (i != 0);
		((double *)(pen_trip->x))[pen_trip->nnz] = MuEnFact;
		pen_trip->nnz++;
	}
	cholmod_sparse *collapser = cholmod_l_triplet_to_sparse(pen_trip, 0, &c);
	cholmod_l_free_triplet(&pen_trip, &c);
	response_matrix = Millipede::GetResponseMatrix(domCache_, *microSources,
	    domEfficiency_, muon_p, cascade_p, NULL, &c);
	if (response_matrix == NULL)
		log_fatal("Null basis matrix");
	little_response_matrix = cholmod_l_ssmult(response_matrix, 
	    collapser, 0, 1, 0, &c);
	cholmod_l_free_sparse(&collapser, &c);
	cholmod_l_free_sparse(&response_matrix, &c);
	if (little_response_matrix == NULL)
		log_fatal("Null basis matrix (little response)");
	MillipedeFitParamsPtr params =
	    boost::make_shared<MillipedeFitParams>();
	Millipede::FitStatistics(domCache_, *sources, I3Units::MeV,
		    little_response_matrix, params.get(), &c);
	log_info("created Fit Statistics");
	cholmod_l_free_sparse(&little_response_matrix, &c);

//	response_matrix = Millipede::GetResponseMatrix(domCache_, *sources,
//	    domEfficiency_, muon_p, cascade_p, NULL, &c);
//	if (response_matrix == NULL)
//		log_fatal("Null basis matrix");
	
//	MillipedeFitParamsPtr params =
//	    boost::make_shared<MillipedeFitParams>();
//	Millipede::FitStatistics(domCache_, *sources, I3Units::MeV,
//		    response_matrix, params.get(), &c);
//	cholmod_l_free_sparse(&response_matrix, &c);
	
	log_info("We Have Parameters!-Now With Longer Matching Response Matrix Calculations!");
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
	//log_info("%f is the muon spacing", muonspacing_); 
	log_info("In the GetLLH Function");
	cholmod_sparse *response_matrix, *little_response_matrix, *gradients, *little_gradients;
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
	unsigned int check2 = trackLength / muonspacing_;
	log_info("%d is the number of int(L/muonspacing)", check2);
//	double check = trackLength - 0.5*muonspacing_;
//	double d=0;

	//create muon segment by grabbing original muon
	I3Particle PrevParticle = (*sources)[1];
	log_info("%d is the particle type of PrevParticle for Muon Looping", PrevParticle.GetType());
	I3Particle NextParticle = PrevParticle;


	//Create Muon that consists of muon points with muon spacing of >= Muonspacing for last segments
	//In case of Muon Length <= 2 >= 1.5 muonspacings, create hypthesis out of 2 muon points
	if (trackLength >= 1.5*muonspacing_){
	microSources->push_back(PrevParticle);
		unsigned int checkend = abs(check2-2);
		log_info("%d is the number of int(abs(L/muonspacing)-2))", checkend);
		double endspace = (trackLength-(checkend)*muonspacing_)/2;
		for (double d = 0; d < checkend; d = d+muonspacing_){
//		for (double d=0; d<check; d = d+muonspacing_){
//		for (; d<check; d = d+muonspacing_){
			//log_info("%f is the length of our Muon now", d);
			//log_info("%f is the total length of our muon", trackLength);
			//log_info("%f is the muon spacing", muonspacing_); 
			NextParticle.SetPos(PrevParticle.ShiftAlongTrack(muonspacing_));
			PrevParticle = NextParticle;
			microSources->push_back(PrevParticle);
		}
		NextParticle.SetPos(PrevParticle.ShiftAlongTrack(endspace));
		log_info("We just made the composite muon");
	}
	//Test Cascade Only Hypothesis for very small Muons
	else {
		(*sources).erase((*sources).end());
		//log_info("%d is the particle type of Source Particle", (*sources)[0].GetType());
		log_info("We just made a cascade-only source");
	}
		
	double MuSeg = (*microSources).size() -1.0;
	double MuEnFact = 1.0;
	if(MuSeg > 0){
		MuEnFact = 1.0/MuSeg;
	}



//	if (d > muonspacing_){
//		(*sources)[1].SetLength(d*I3Units::m);
//	}

//	double MuLen = (*sources)[1].GetLength();

	// Make a matrix of ones and zeros for collapsing the multi-muon response matrix
	// into a single muon response matrix (little_response_matrix)
	// (a triplet is: matrix position i,j; and value x)
	cholmod_triplet *pen_trip = cholmod_l_allocate_triplet(
	    (*microSources).size(), (*sources).size(), (*microSources).size(), 0,
	    CHOLMOD_REAL, &c);
	// starting number of non-zero entries in this matrix:
	pen_trip->nnz = 0;
//	for(unsigned row = 0; row < pen_trip->nrow; row++){
//		for(unsigned col = 0; col < pen_trip -> ncol; col ++){
//			if(row==0 && col==1){
//				((double *)(pen_trip->x))[row + col*pen_trip->nrow] = 1;
//			((double *)(pen_trip->x))[row + col*pen_trip->nrow] = 1;
//				pen_trip->nnz++;
//			} else if (row !=0 && col !=0){
//				((double *)(pen_trip->x))[row + col*pen_trip->nrow] = MuEnFact;
//				pen_trip->nnz++;
//			}
//		} else {
//			((double *)(pen_trip->x))[pen_trip->nnz] = MuEnFact;
//			pen_trip->nnz++;
//		}
//	}


	unsigned i = 0; 
		((long *)(pen_trip->i))[pen_trip->nnz] = i;
		((long *)(pen_trip->j))[pen_trip->nnz] = (i != 0);
		((double *)(pen_trip->x))[pen_trip->nnz] = 1;
		pen_trip->nnz++;
	//Energy Correction for Composite Muon
	for (unsigned i = 1; i < pen_trip->nrow; i++) {
		((long *)(pen_trip->i))[pen_trip->nnz] = i;
		((long *)(pen_trip->j))[pen_trip->nnz] = (i != 0);
		((double *)(pen_trip->x))[pen_trip->nnz] = MuEnFact;
		pen_trip->nnz++;
	}
	cholmod_sparse *collapser = cholmod_l_triplet_to_sparse(pen_trip, 0, &c);
	cholmod_l_free_triplet(&pen_trip, &c);
//	cholmod_dense *one_muon_mult =
//	    cholmod_l_sparse_to_dense(collapser, &c);
//	for(unsigned row = 0; row < collapser->nrow; row++){
//		for(unsigned col = 0; col < collapser->ncol; col++){
//		double elem = ((double*)(one_muon_mult->x))
//			[row + col*one_muon_mult->nrow];
//		std::cout<<row<<","<<col<<" : "<<elem<<std::endl;
//		}
//	}

	// Make matrix for converting multi-muon gradients to single muon little_gradients
	cholmod_triplet *grad_trip = cholmod_l_allocate_triplet(
	    (*microSources).size()*6, 12, (*microSources).size()*6, 0,
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
	for (unsigned i = 6; i < grad_trip->nrow; i+=6) {
		((long *)(grad_trip->i))[grad_trip->nnz] = i;
		((long *)(grad_trip->j))[grad_trip->nnz] = 6;
		((double *)(grad_trip->x))[grad_trip->nnz] = 1;
		grad_trip->nnz++;
	}
	for (unsigned i = 7; i < grad_trip->nrow; i+=6) {
		((long *)(grad_trip->i))[grad_trip->nnz] = i;
		((long *)(grad_trip->j))[grad_trip->nnz] = 7;
		((double *)(grad_trip->x))[grad_trip->nnz] = 1;
		grad_trip->nnz++;
	}
	for (unsigned i = 8; i < grad_trip->nrow; i+=6) {
		((long *)(grad_trip->i))[grad_trip->nnz] = i;
		((long *)(grad_trip->j))[grad_trip->nnz] = 8;
		((double *)(grad_trip->x))[grad_trip->nnz] = 1;
		grad_trip->nnz++;
	}
	for (unsigned i = 9; i < grad_trip->nrow; i+=6) {
		((long *)(grad_trip->i))[grad_trip->nnz] = i;
		((long *)(grad_trip->j))[grad_trip->nnz] = 9;
		((double *)(grad_trip->x))[grad_trip->nnz] = 1;
		grad_trip->nnz++;
	}
	for (unsigned i = 10; i < grad_trip->nrow; i+=6) {
		((long *)(grad_trip->i))[grad_trip->nnz] = i;
		((long *)(grad_trip->j))[grad_trip->nnz] = 10;
		((double *)(grad_trip->x))[grad_trip->nnz] = 1;
		grad_trip->nnz++;
	}
	for (unsigned i = 11; i < grad_trip->nrow; i+=6) {
		((long *)(grad_trip->i))[grad_trip->nnz] = i;
		((long *)(grad_trip->j))[grad_trip->nnz] = 11;
		((double *)(grad_trip->x))[grad_trip->nnz] = 1;
		grad_trip->nnz++;
	}
	cholmod_sparse *grad_collapser = cholmod_l_triplet_to_sparse(grad_trip, 0, &c);
	cholmod_l_free_triplet(&grad_trip, &c);
//	cholmod_dense *one_muon_grad_mult =
//	    cholmod_l_sparse_to_dense(grad_collapser, &c);
//	for(unsigned row = 0; row < grad_collapser->nrow; row++){
//		for(unsigned col = 0; col < grad_collapser->ncol; col++){
//		double elem = ((double*)(one_muon_grad_mult->x))
//			[row + col*one_muon_grad_mult->nrow];
//		std::cout<<row<<","<<col<<" : "<<elem<<std::endl;
//		}
//	}

	response_matrix = Millipede::GetResponseMatrix(domCache_, *microSources,
	    domEfficiency_, muon_p, cascade_p,
	    (gradient == NULL) ? NULL : &gradients, &c);
	log_info("Response_matrix has been defined");
	if (response_matrix == NULL)
		log_fatal("Null basis matrix");
//	cholmod_dense *response =
//	    cholmod_l_sparse_to_dense(response_matrix, &c);
//	for(unsigned row = 0; row < response_matrix->nrow; row++){
//		for(unsigned col = 0; col < response_matrix->ncol; col++){
//		double elem = ((double*)(response->x))
//			[row + col*response->nrow];
//		std::cout<<"response_matrix"<<row<<","<<col<<" : "<<elem<<std::endl;
//		}
//	}



	little_gradients = cholmod_l_ssmult(gradients, 
	    grad_collapser, 0, 1, 0, &c);

	log_info("Uno");

//	cholmod_dense *little_gradient_fuckery =
//	    cholmod_l_sparse_to_dense(little_gradients, &c);
//	for(unsigned row = 0; row < little_gradients->nrow; row++){
//		for(unsigned col = 0; col < little_gradients->ncol; col++){
//		double elem = ((double*)(little_gradient_fuckery->x))
//			[row + col*little_gradient_fuckery->nrow];
//		std::cout<<row<<","<<col<<" : "<<elem<<std::endl;
//		}
//	}
	log_info("Dos");
//	cholmod_l_free_dense(&little_gradient_fuckery, &c);

	log_info("Combined Muon Response Matrix Created");

	cholmod_l_free_sparse(&grad_collapser, &c);


	log_info("I have a bad feeling about this"); //this is the magic seg-fault resolving log_info statement

	// Colapse the resulting big matrix back into a 2-particle matrix
	little_response_matrix = cholmod_l_ssmult(response_matrix, 
	    collapser, 0, 1, 0, &c);
	cholmod_l_free_sparse(&collapser, &c);
	cholmod_l_free_sparse(&response_matrix, &c);

	

//	cholmod_dense *tiny_response =
//	    cholmod_l_sparse_to_dense(little_response_matrix, &c);
//	for(unsigned row = 0; row < little_response_matrix->nrow; row++){
//		for(unsigned col = 0; col < little_response_matrix->ncol; col++){
//		double elem = ((double*)(tiny_response->x))
//			[row + col*tiny_response->nrow];
//		std::cout<<"little_response_matrix"<<row<<","<<col<<" : "<<elem<<std::endl;
//		}
//	}
	
	if (fit_energy) {
		log_info("Solve for energy losses");
		SolveEnergyLosses(*sources, little_response_matrix,
	//Needs to either be micro and resp or src and little NOT micro and little
		    (gradient == NULL ? NULL : little_gradients));
		if (sources->size() == 1)
			hypo.particle->SetEnergy((*sources)[0].GetEnergy());
		log_info("We should now have a nice energy");
	}


	llh = Millipede::FitStatistics(domCache_, *sources, I3Units::MeV,
	    little_response_matrix, NULL, &c);
		//Dec 4th: above was micro and response....why did I do that? changed to sources and little, let's see if the llh goes crazy?
		log_info("[%f m mu, vertex (%f, %f, %f)] + ", trackLength, x_v, y_v, z_v);
		log_info("[ (%f zen_mu, %f zen_casc), (%f azi_mu, %f azi_casc)] -> (llh=%f)", zen_track, zen_cascade, azi_track, azi_casc, llh);

	if (gradient != NULL) {
		log_info("Begin non-null gradients loop");
//		assert(sources->size() == gradsources->size());
		Millipede::LLHGradient(domCache_, *sources, *gradsources,
		    I3Units::MeV, weight, little_response_matrix, little_gradients, &c);
	//Jun 20, 2013 replace and resp w/little_resp
	//SOMETHING IN THE ABOVE STATEMENT IS CRYING OUT TO CHOLMOD
		log_info("Just finished LLHGradient calculation");
		cholmod_l_free_sparse(&little_gradients, &c);
		cholmod_l_free_sparse(&gradients, &c);
		if (sources->size() == 1) {
			// NB: the requested weight is already applied in
			// the call to LLHGradient()
			I3Particle &p = *gradient->particle;
			log_info("Only Cascade Grads Calculation");
			p.SetPos(I3Position(
			    p.GetPos().GetX() +
			      (*gradsources)[0].GetPos().GetX(),
			    p.GetPos().GetY() +
			      (*gradsources)[0].GetPos().GetY(),
			    p.GetPos().GetZ() +
			      (*gradsources)[0].GetPos().GetZ()
			));
			log_info("Increased Position");
			p.SetTime(p.GetTime() + (*gradsources)[0].GetTime());
			log_info("Increased Time");
			p.SetDir(I3Direction(
			    p.GetDir().GetZenith() +
			      (*gradsources)[0].GetDir().GetZenith(),
			    p.GetDir().GetAzimuth() +
			      (*gradsources)[0].GetDir().GetAzimuth()
			));
			log_info("Increased Zen and Azi");			
		}
		for(int i=0; i<gradsources->size(); i++) log_warn_stream((*gradsources)[i]); //info->warn for error checking
	}
	if (gradient == NULL) {
	log_info("NULL GRADIENT MATRIX");
	}

	cholmod_l_free_sparse(&little_response_matrix, &c);
	log_info("End of LLH Function");
	return llh;
}

typedef I3SingleServiceFactory<BipedLikelihood, I3EventLogLikelihoodBase>
    BipedLikelihoodFactory;
I3_SERVICE_FACTORY(BipedLikelihoodFactory);
