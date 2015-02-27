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
	cholmod_sparse *response_matrix, *many_response_matrix;
	I3VectorI3ParticlePtr sources = ExtractHypothesis(hypo);
	
	log_info("Diagnostics Function");
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
	double MuEnFact = 1.0/MuSeg;
	double endscale = MuEnFact*(trackLength - (MuSeg-1)*muonspacing_)/muonspacing_;
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
		pen_trip->nnz++;
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
	unsigned int MuSeg = ceil(trackLength / muonspacing_);
	log_warn("%d is the number of ceil(L/muonspacing)", MuSeg);
//	double check = trackLength - 0.5*muonspacing_;
//	double d=0;

	//create muon segment by grabbing original muon
	I3Particle PrevParticle = (*sources)[1];
	log_info("%d is the particle type of PrevParticle for Muon Looping", PrevParticle.GetType());
	I3Particle NextParticle = PrevParticle;
        std::vector<double> muonLengths;
	//Create Muon that consists of muon points with muon spacing of >= Muonspacing for last segments
	//In case of Muon Length <= 2 >= 1.5 muonspacings, create hypthesis out of 2 muon points
//	if (trackLength >= 1.5*muonspacing_){
//		microSources->push_back(PrevParticle);
//		unsigned int checkend = abs(check2-2);
//		log_info("%d is the number of int(abs(L/muonspacing)-2))", checkend);
//		double endspace = (trackLength-(checkend)*muonspacing_)/2;
	for (double d = 0; d < MuSeg*muonspacing_; d = d+muonspacing_){
		microSources->push_back(PrevParticle);
//		for (double d=0; d<check; d = d+muonspacing_){
//		for (; d<check; d = d+muonspacing_){
			//log_info("%f is the length of our Muon now", d);
			//log_info("%f is the total length of our muon", trackLength);
			//log_info("%f is the muon spacing", muonspacing_); 
		NextParticle.SetPos(PrevParticle.ShiftAlongTrack(muonspacing_));
		NextParticle.SetTime(PrevParticle.GetTime()+muonspacing_/I3Constants::c);
		PrevParticle = NextParticle;
                muonLengths.push_back(d);
	}
//	NextParticle.SetTime(PrevParticle.GetTime()+(muonLengths.back() + endspace)/I3Constants::c);
//	muonLengths.push_back(muonLengths.back()+endspace);
//	NextParticle.SetPos(PrevParticle.ShiftAlongTrack(endspace));
//	microSources->push_back(NextParticle);
	log_info("We just made the composite muon");
//	std::vector<double> muonLengths((*microSources).size()-1, muonspacing_);
//	for (unsigned i = 1; i < (*microSources).size()-2; i++)
//		muonLengths[i] = muonLengths[i-1]+muonspacing_;
//	for (unsigned i = (*microSources).size()-2; i <= (*microSources).size(); i++)
//		muonLengths[i] = muonLengths[i-1]+endspace;
//	}
	//Only one muon segment for very short Muons
//	else {
//	microSources->push_back(PrevParticle);
//	muonLengths.push_back(0.0);
//		(*sources).erase((*sources).end());
//		log_info("%d is the particle type of Source Particle", (*sources)[0].GetType());
//		log_info("We just made a cascade-only source");
//	}

	//Create a energy scaling factor		
	double MuEnFact = 1.0/MuSeg;
	double endscale = MuEnFact*(trackLength - (MuSeg-1)*muonspacing_)/muonspacing_;
	log_warn("%f is the endscale", endscale);




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
	cholmod_sparse *collapser = cholmod_l_triplet_to_sparse(pen_trip, 0, &c);
	cholmod_l_free_triplet(&pen_trip, &c);
//	cholmod_dense *one_muon_mult =
//	cholmod_l_sparse_to_dense(collapser, &c);
//	for(unsigned row = 0; row < collapser->nrow; row++){
//		for(unsigned col = 0; col < collapser->ncol; col++){
//		double elem = ((double*)(one_muon_mult->x))
//			[row + col*one_muon_mult->nrow];
//		std::cout<<row<<","<<col<<" : "<<elem<<std::endl;
//		}
//	}
	log_info("created matrix to collapse response matrix");

//	cholmod_triplet *length_grad_part_trip = cholmod_l_allocate_triplet(
//	    (*microSources).size(), 14, 1, 0, CHOLMOD_REAL, &c);
	// starting number of non-zero entries in this matrix:
//	length_grad_part_trip->nnz = 0;
//	((long *)(length_grad_part_trip->i))[length_grad_part_trip->nnz] = length_grad_part_trip->nrow-1;
//	((long *)(length_grad_part_trip->j))[length_grad_part_trip->nnz] = length_grad_part_trip->ncol-1;
//	((double *)(length_grad_part_trip->x))[length_grad_part_trip->nnz] = MuEnFact/muonspacing_;
//	length_grad_part_trip->nnz = 1;
//	cholmod_sparse *length_grad_part = cholmod_l_triplet_to_sparse(length_grad_part_trip, 0, &c);
//	cholmod_l_free_triplet(&length_grad_part_trip, &c);
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
	log_info("created matrix to grab out length gradient from the response matrix");


	many_response_matrix = Millipede::GetResponseMatrix(domCache_, *microSources,
	    domEfficiency_, muon_p, cascade_p,
	    (gradient == NULL) ? NULL : &manygradients, &c);
	log_info("Response_matrix has been defined");
	if (many_response_matrix == NULL){
		log_fatal("Null basis matrix");
	}
//	cholmod_dense *response =
//	    cholmod_l_sparse_to_dense(response_matrix, &c);
//	for(unsigned row = 0; row < response_matrix->nrow; row++){
//		for(unsigned col = 0; col < response_matrix->ncol; col++){
//		double elem = ((double*)(response->x))
//			[row + col*response->nrow];
//		std::cout<<"response_matrix"<<row<<","<<col<<" : "<<elem<<std::endl;
//		}
//	}
	// Colapse the resulting big matrix back into a 2-particle matrix
	response_matrix = cholmod_l_ssmult(many_response_matrix, 
	    collapser, 0, 1, 0, &c);
	cholmod_l_free_sparse(&collapser, &c);
	//Get the length gradient
	cholmod_sparse *length_grad = cholmod_l_ssmult(many_response_matrix, 
	    length_grad_part, 0, 1, 0, &c);
	cholmod_l_free_sparse(&length_grad_part, &c);
	cholmod_l_free_sparse(&many_response_matrix, &c);
	log_info("Made length gradient");



	// Make matrix for converting multi-muon gradients to single muon little_gradients
	cholmod_triplet *grad_trip = cholmod_l_allocate_triplet(
	    (*microSources).size()*7, 14, (*microSources).size()*12, 0,
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
	for (unsigned i = 7; i < grad_trip->nrow-1; i+=7) {
		((long *)(grad_trip->i))[grad_trip->nnz] = i;
		((long *)(grad_trip->j))[grad_trip->nnz] = 7;
		((double *)(grad_trip->x))[grad_trip->nnz] = MuEnFact;
		grad_trip->nnz++;
	}
	((long *)(grad_trip->i))[grad_trip->nnz] = grad_trip->nrow-1;
	((long *)(grad_trip->j))[grad_trip->nnz] = 7;
	((double *)(grad_trip->x))[grad_trip->nnz] = endscale;
	grad_trip->nnz++;
	for (unsigned i = 8; i < grad_trip->nrow-1; i+=7) {
		((long *)(grad_trip->i))[grad_trip->nnz] = i;
		((long *)(grad_trip->j))[grad_trip->nnz] = 8;
		((double *)(grad_trip->x))[grad_trip->nnz] = MuEnFact;
		grad_trip->nnz++;
	}
	((long *)(grad_trip->i))[grad_trip->nnz] = grad_trip->nrow-1;
	((long *)(grad_trip->j))[grad_trip->nnz] = 8;
	((double *)(grad_trip->x))[grad_trip->nnz] = endscale;
	grad_trip->nnz++;
	for (unsigned i = 9; i < grad_trip->nrow-1; i+=7) {
		((long *)(grad_trip->i))[grad_trip->nnz] = i;
		((long *)(grad_trip->j))[grad_trip->nnz] = 9;
		((double *)(grad_trip->x))[grad_trip->nnz] = MuEnFact;
		grad_trip->nnz++;
	}
	((long *)(grad_trip->i))[grad_trip->nnz] = grad_trip->nrow-1;
	((long *)(grad_trip->j))[grad_trip->nnz] = 9;
	((double *)(grad_trip->x))[grad_trip->nnz] = endscale;
	grad_trip->nnz++;
	for (unsigned i = 10; i < grad_trip->nrow-1; i+=7) {
		((long *)(grad_trip->i))[grad_trip->nnz] = i;
		((long *)(grad_trip->j))[grad_trip->nnz] = 10;
		((double *)(grad_trip->x))[grad_trip->nnz] = MuEnFact;
		grad_trip->nnz++;
	}
	((long *)(grad_trip->i))[grad_trip->nnz] = grad_trip->nrow-1;
	((long *)(grad_trip->j))[grad_trip->nnz] = 10;
	((double *)(grad_trip->x))[grad_trip->nnz] = endscale;
	grad_trip->nnz++;

	//zenith	
	//includes lever arm effect
	double zen((*microSources)[1].GetZenith()), azi((*microSources)[1].GetAzimuth());
	for (unsigned i = 11; i < grad_trip->nrow-1; i+=7) {
		unsigned nom=0;
		((long *)(grad_trip->i))[grad_trip->nnz] = i;
		((long *)(grad_trip->j))[grad_trip->nnz] = 11;
		((double *)(grad_trip->x))[grad_trip->nnz] = MuEnFact;
		grad_trip->nnz++;
		((long *)(grad_trip->i))[grad_trip->nnz] = i-4;
		((long *)(grad_trip->j))[grad_trip->nnz] = 11;
		((double *)(grad_trip->x))[grad_trip->nnz] = -std::cos(zen)*std::cos(azi)*MuEnFact*muonLengths[nom];
		grad_trip->nnz++;
		((long *)(grad_trip->i))[grad_trip->nnz] = i-3;
		((long *)(grad_trip->j))[grad_trip->nnz] = 11;
		((double *)(grad_trip->x))[grad_trip->nnz] = -std::cos(zen)*std::sin(azi)*MuEnFact*muonLengths[nom];
		grad_trip->nnz++;
		((long *)(grad_trip->i))[grad_trip->nnz] = i-2;
		((long *)(grad_trip->j))[grad_trip->nnz] = 11;
		((double *)(grad_trip->x))[grad_trip->nnz] = std::sin(zen)*MuEnFact*muonLengths[nom];
		grad_trip->nnz++;
		nom++;
	}
	((long *)(grad_trip->i))[grad_trip->nnz] = grad_trip->nrow-1;
	((long *)(grad_trip->j))[grad_trip->nnz] = 11;
	((double *)(grad_trip->x))[grad_trip->nnz] = endscale;
	grad_trip->nnz++;
	((long *)(grad_trip->i))[grad_trip->nnz] = grad_trip->nrow-5;
	((long *)(grad_trip->j))[grad_trip->nnz] = 11;
	((double *)(grad_trip->x))[grad_trip->nnz] = -std::cos(zen)*std::cos(azi)*endscale*muonLengths.back();
	grad_trip->nnz++;
	((long *)(grad_trip->i))[grad_trip->nnz] = grad_trip->nrow-4;
	((long *)(grad_trip->j))[grad_trip->nnz] = 11;
	((double *)(grad_trip->x))[grad_trip->nnz] = -std::cos(zen)*std::sin(azi)*endscale*muonLengths.back();
	grad_trip->nnz++;
	((long *)(grad_trip->i))[grad_trip->nnz] = grad_trip->nrow-3;
	((long *)(grad_trip->j))[grad_trip->nnz] = 11;
	((double *)(grad_trip->x))[grad_trip->nnz] = std::sin(zen)*endscale*muonLengths.back();
	grad_trip->nnz++;
	//azimuth
	//includes lever arm effect
	for (unsigned i = 12; i < grad_trip->nrow-1; i+=7) {
		unsigned nom = 0;
		((long *)(grad_trip->i))[grad_trip->nnz] = i;
		((long *)(grad_trip->j))[grad_trip->nnz] = 12;
		((double *)(grad_trip->x))[grad_trip->nnz] = MuEnFact;
		grad_trip->nnz++;
		((long *)(grad_trip->i))[grad_trip->nnz] = i-5;
		((long *)(grad_trip->j))[grad_trip->nnz] = 12;
		((double *)(grad_trip->x))[grad_trip->nnz] = std::sin(zen)*std::sin(azi)*MuEnFact*muonLengths[nom];
		grad_trip->nnz++;
		((long *)(grad_trip->i))[grad_trip->nnz] = i-4;
		((long *)(grad_trip->j))[grad_trip->nnz] = 12;
		((double *)(grad_trip->x))[grad_trip->nnz] = -std::sin(zen)*std::cos(azi)*MuEnFact*muonLengths[nom];
		grad_trip->nnz++;
		nom++;
	}
	((long *)(grad_trip->i))[grad_trip->nnz] = grad_trip->nrow-1;
	((long *)(grad_trip->j))[grad_trip->nnz] = 12;
	((double *)(grad_trip->x))[grad_trip->nnz] = endscale;
	grad_trip->nnz++;
	((long *)(grad_trip->i))[grad_trip->nnz] = grad_trip->nrow-6;
	((long *)(grad_trip->j))[grad_trip->nnz] = 12;
	((double *)(grad_trip->x))[grad_trip->nnz] = std::sin(zen)*std::sin(azi)*endscale*muonLengths.back();
	grad_trip->nnz++;
	((long *)(grad_trip->i))[grad_trip->nnz] = grad_trip->nrow-5;
	((long *)(grad_trip->j))[grad_trip->nnz] = 12;
	((double *)(grad_trip->x))[grad_trip->nnz] = -std::sin(zen)*std::cos(azi)*endscale*muonLengths.back();
	grad_trip->nnz++;
	//length
	//length depends upon the placement of and spacing of the muon segments
	//changes in the length only change the last two muon segments positions
	//all other muons remain in their previous location for epsilon changes in length
//	std::vector<double> xLContribution(MuSeg, abs (std::sin(zen)*std::cos(azi)));
//	std::vector<double> yLContribution(MuSeg, abs (std::sin(zen)*std::sin(azi)));
//	std::vector<double> zLContribution(MuSeg, abs (std::cos(zen)));
//	if ((*microSources).size() > 2){
//		for(unsigned i = 7+7*((*microSources).size()-3); i < grad_trip->nrow; i+=7) {
//			unsigned nom = 0;
//			((long *)(grad_trip->i))[grad_trip->nnz] = i;
//			((long *)(grad_trip->j))[grad_trip->nnz] = 13;
//			((double *)(grad_trip->x))[grad_trip->nnz] = MuEnFact*xLContribution[nom];
//			grad_trip->nnz++;
//			((long *)(grad_trip->i))[grad_trip->nnz] = i+1;
//			((long *)(grad_trip->j))[grad_trip->nnz] = 13;
//			((double *)(grad_trip->x))[grad_trip->nnz] = MuEnFact*yLContribution[nom];
//			grad_trip->nnz++;
//			((long *)(grad_trip->i))[grad_trip->nnz] = i+2;
//			((long *)(grad_trip->j))[grad_trip->nnz] = 13;
//			((double *)(grad_trip->x))[grad_trip->nnz] = MuEnFact*zLContribution[nom];
//			grad_trip->nnz++;
//			nom ++;
//		}
//	}
	//for the case of just one muon, the effective length (end point) of the muon segment
	//changes for epsilon scale changes in L even though the vertex remains fixed
//	else{
//		unsigned nom = 0;
//		((long *)(grad_trip->i))[grad_trip->nnz] = 7;
//		((long *)(grad_trip->j))[grad_trip->nnz] = 13;
//		((double *)(grad_trip->x))[grad_trip->nnz] = xLContribution[nom];
//		grad_trip->nnz++;
//		((long *)(grad_trip->i))[grad_trip->nnz] = 8;
//		((long *)(grad_trip->j))[grad_trip->nnz] = 13;
//		((double *)(grad_trip->x))[grad_trip->nnz] = yLContribution[nom];
//		grad_trip->nnz++;
//		((long *)(grad_trip->i))[grad_trip->nnz] = 9;
//		((long *)(grad_trip->j))[grad_trip->nnz] = 13;
//		((double *)(grad_trip->x))[grad_trip->nnz] = zLContribution[nom];
//		grad_trip->nnz++;
//	}
	
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
	log_info("Made the spatial gradient collapser");




	cholmod_sparse *spatial_grad = cholmod_l_ssmult(manygradients, 
	    grad_collapser, 0, 1, 0, &c);
	cholmod_l_free_sparse(&manygradients, &c);
	log_info("Made the spatial gradient");
	double alpha[2] = {1,0}, beta[2] = {1,0};
	gradients = cholmod_l_add(spatial_grad, length_grad, alpha, beta, true, true, &c);
	log_info("Combined the spatial and length gradients");
	cholmod_l_free_sparse(&spatial_grad, &c);
	cholmod_l_free_sparse(&length_grad, &c);


//	log_info("Uno");

//	cholmod_dense *little_gradient_fuckery =
//	    cholmod_l_sparse_to_dense(little_gradients, &c);
//	for(unsigned row = 0; row < little_gradients->nrow; row++){
//		for(unsigned col = 0; col < little_gradients->ncol; col++){
//		double elem = ((double*)(little_gradient_fuckery->x))
//			[row + col*little_gradient_fuckery->nrow];
//		std::cout<<row<<","<<col<<" : "<<elem<<std::endl;
//		}
//	}
//	log_info("Dos");
//	cholmod_l_free_dense(&little_gradient_fuckery, &c);

	log_info("Combined Muon Response Matrix Created");

	cholmod_l_free_sparse(&grad_collapser, &c);


	log_info("I have a bad feeling about this"); 



	

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
		SolveEnergyLosses(*sources, response_matrix,
	//Needs to either be micro and resp or src and little NOT micro and little
		    (gradient == NULL ? NULL : gradients));
		if (sources->size() == 1)
			hypo.particle->SetEnergy((*sources)[0].GetEnergy());
		log_info("We should now have a nice energy");
	}


	llh = Millipede::FitStatistics(domCache_, *sources, I3Units::MeV,
	    response_matrix, NULL, &c);
		//Dec 4th: above was micro and response....why did I do that? changed to sources and little, let's see if the llh goes crazy?
		log_info("[%f m mu, vertex (%f, %f, %f)] + ", trackLength, x_v, y_v, z_v);
		log_info("[ (%f zen_mu, %f zen_casc), (%f azi_mu, %f azi_casc)] -> (llh=%f)", zen_track, zen_cascade, azi_track, azi_casc, llh);

	if (gradient != NULL) {
		log_info("Begin non-null gradients loop");
//		assert(sources->size() == gradsources->size());
		Millipede::LLHGradient(domCache_, *sources, *gradsources,
		    I3Units::MeV, weight, response_matrix, gradients, &c);
	//Jun 20, 2013 replace and resp w/little_resp
	//SOMETHING IN THE ABOVE STATEMENT IS CRYING OUT TO CHOLMOD
		log_info("Just finished LLHGradient calculation");
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

	cholmod_l_free_sparse(&response_matrix, &c);
	log_info("End of LLH Function");
	return llh;
}

typedef I3SingleServiceFactory<BipedLikelihood, I3EventLogLikelihoodBase>
    BipedLikelihoodFactory;
I3_SERVICE_FACTORY(BipedLikelihoodFactory);
