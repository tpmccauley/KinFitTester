#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicVertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"

#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include <iostream>

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

class KinFitTester : public edm::one::EDAnalyzer<edm::one::SharedResources>  
{
public:
  explicit KinFitTester(const edm::ParameterSet&);
  ~KinFitTester();
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  edm::InputTag patMuonInputTag_;
  edm::InputTag prunedGenParticleInputTag_;
  edm::InputTag packedGenParticleInputTag_;

  bool isMC_;

  double ptMin_;
  double etaMax_;

  edm::EDGetTokenT<std::vector<pat::Muon> > patMuonToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> prunedGenParticleToken_;
  edm::EDGetTokenT<pat::PackedGenParticleCollection> packedGenParticleToken_;

  bool isAncestor(const reco::Candidate* ancestor, const reco::Candidate* particle);
  bool inKinematicAcceptanceRegion(const pat::Muon& muon);
  bool isGoodMuon(const pat::Muon& muon);
  bool isKalmanVertexFit(std::vector<reco::TransientTrack>);

  void KinematicParticleVertexFit(std::vector<RefCountedKinematicParticle>);

  void KinematicParticleVertexFit(std::vector<RefCountedKinematicParticle>,
                                  std::vector<RefCountedKinematicParticle>);

  void KinematicConstrainedVertexFit(std::vector<RefCountedKinematicParticle>);

  void printout(const RefCountedKinematicVertex& vertex) const;
  void printout(const RefCountedKinematicParticle& particle) const;
  void printout(const RefCountedKinematicTree& tree) const;
};

KinFitTester::KinFitTester(const edm::ParameterSet& iConfig)
  : patMuonInputTag_(iConfig.getParameter<edm::InputTag>("patMuonTag")),
    prunedGenParticleInputTag_(iConfig.getParameter<edm::InputTag>("prunedGenParticleTag")),
    packedGenParticleInputTag_(iConfig.getParameter<edm::InputTag>("packedGenParticleTag")),
    isMC_(iConfig.getParameter<bool>("isMC")),
    ptMin_(iConfig.getParameter<double>("ptMin")),
    etaMax_(iConfig.getParameter<double>("etaMax"))
{
  patMuonToken_ = consumes<std::vector<pat::Muon> >(patMuonInputTag_);
  prunedGenParticleToken_ = consumes<reco::GenParticleCollection>(prunedGenParticleInputTag_);
  packedGenParticleToken_ = consumes<pat::PackedGenParticleCollection>(packedGenParticleInputTag_);
}

KinFitTester::~KinFitTester()
{}

void KinFitTester::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<std::vector<pat::Muon> > pat_muons;
  iEvent.getByToken(patMuonToken_, pat_muons);
  std::cout<< pat_muons->size() << " pat::Muons in the event" <<std::endl;

  edm::ESHandle<TransientTrackBuilder> ttb;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", ttb);

  std::vector<reco::TransientTrack> ttrks;

  if ( isMC_ ) 
  {
    edm::Handle<reco::GenParticleCollection> pruned;
    iEvent.getByToken(prunedGenParticleToken_, pruned);
    std::cout<< pruned->size()<<" reco::GenParticles in the event"<<std::endl;

    edm::Handle<pat::PackedGenParticleCollection> packed;
    iEvent.getByToken(packedGenParticleToken_, packed);
    std::cout<< packed->size() <<" pat::PackedGenParticles in the event"<<std::endl;

    std::vector<const reco::Candidate*> cands;    

    for ( size_t i = 0; i < pruned->size(); i++ )
    {
      if ( abs((*pruned)[i].pdgId()) == 531 )
      {
        const reco::Candidate* Bs = &(*pruned)[i];
        
        for ( size_t j = 0; j < packed->size(); j++ )
        { 
          const reco::Candidate * motherInPrunedCollection = (*packed)[j].mother(0);
 
          if ( motherInPrunedCollection != nullptr && isAncestor(Bs, motherInPrunedCollection) )
          {
            if ( abs((*packed)[j].pdgId()) == 13 || abs((*packed)[j].pdgId()) == 321 ) 
            {
              std::cout<<"  pdgId, pt, eta, phi: "
                       << (*packed)[j].pdgId() <<",  "
                       << (*packed)[j].pt() <<", "
                       << (*packed)[j].eta() <<", "
                       << (*packed)[j].phi() <<std::endl;
                
              cands.push_back((*packed)[j].clone());
            }            
          }     
        }
      }
    }   

    /*
      Fix this
    if ( cands.size() > 4 )
      return;

    for ( auto cand: cands )
    {
      if ( ! cand )
        continue;
    
      // Dear God C++
      ttrks.push_back((*ttb).build(cand->bestTrack()));
    }
 
    std::cout<< ttrks.size() <<" TransientTracks"<<std::endl;

    KinematicParticleFactoryFromTransientTrack kp_factory;

    float muon_mass  = 0.1056583;
    float muon_sigma = 0.0000001;

    float kaon_mass  = 0.493677;
    float kaon_sigma = 0.000016;
  
    float chi2 = 0.0;
    float ndf  = 0.0;

    std::vector<RefCountedKinematicParticle> particles;
    std::vector<RefCountedKinematicParticle> muons;
    std::vector<RefCountedKinematicParticle> kaons;
  
    particles.push_back(kp_factory.particle(ttrks[0], muon_mass, chi2, ndf, muon_sigma));
    particles.push_back(kp_factory.particle(ttrks[1], muon_mass, chi2, ndf, muon_sigma));
    particles.push_back(kp_factory.particle(ttrks[2], kaon_mass, chi2, ndf, kaon_sigma));
    particles.push_back(kp_factory.particle(ttrks[3], kaon_mass, chi2, ndf, kaon_sigma));

    muons.push_back(kp_factory.particle(ttrks[0], muon_mass, chi2, ndf, muon_sigma));
    muons.push_back(kp_factory.particle(ttrks[1], muon_mass, chi2, ndf, muon_sigma));

    kaons.push_back(kp_factory.particle(ttrks[2], kaon_mass, chi2, ndf, kaon_sigma));
    kaons.push_back(kp_factory.particle(ttrks[3], kaon_mass, chi2, ndf, kaon_sigma));
  
    KinematicConstrainedVertexFit(particles);  

    KinematicParticleVertexFit(muons, kaons);

    KinematicParticleVertexFit(particles);
    */

    /*
    unsigned int n_gen_muons = 0;
    unsigned int n_gen_kaons = 0;

    for ( pat::PackedGenParticleCollection::const_iterator gp = packed->begin();
          gp != packed->end(); ++gp ) 
    {
      if ( abs(gp->pdgId()) == 13 )
      {
        n_gen_muons++;    
      } 

      if ( abs(gp->pdgId()) == 321 )
      {
        n_gen_kaons++;
      }
    }

    std::cout<< n_gen_muons <<" PackedGenParticles are muons"<<std::endl;
    std::cout<< n_gen_kaons <<" PackedGenParticles are kaons"<<std::endl;
    */
  }

  /*
  std::vector<pat::Muon> gms;

  for ( std::vector<pat::Muon>::const_iterator m = pat_muons->begin(), mend = pat_muons->end();
        m != mend; ++m ) 
  {   
    if ( isGoodMuon(*m) && inKinematicAcceptanceRegion(*m) ) 
    {
      trks.push_back(m->innerTrack().get());
      gms.push_back(*m);
      
      std::cout<<"good muon:"
               <<" pt = "<< m->pt()
               <<" charge = "<< m->charge()
               <<std::endl;
    }
    
  }
 
  for ( auto trk: trks )
  {
    if ( ! trk )
      continue;
    
    ttrks.push_back((*ttb).build(trk));
  }
 
  std::cout<< ttrks.size() <<" TransientTracks"<<std::endl;

  if ( ttrks.size() != 2 )
    return;

  if ( ! isKalmanVertexFit(ttrks) ) 
    return;
  
  double M = 2*gms[0].pt()*gms[1].pt();
  M *= cosh(gms[0].eta()-gms[1].eta())-cos(gms[0].phi()-gms[1].phi());
  M = sqrt(M);
  
  if ( ! ( M < 3.467 && M > 2.947 ) ) 
    return;
  */ 
  /*
  KinematicParticleFactoryFromTransientTrack kp_factory;

  float muon_mass  = 0.1056583;
  float muon_sigma = 0.0000001;

  float kaon_mass  = 0.493677;
  float kaon_sigma = 0.000016;
  
  float chi2 = 0.0;
  float ndf  = 0.0;

  std::vector<RefCountedKinematicParticle> particles;
  std::vector<RefCountedKinematicParticle> muons;
  std::vector<RefCountedKinematicParticle> kaons;
  
  particles.push_back(kp_factory.particle(ttrks[0], muon_mass, chi2, ndf, muon_sigma));
  particles.push_back(kp_factory.particle(ttrks[1], muon_mass, chi2, ndf, muon_sigma));
  particles.push_back(kp_factory.particle(ttrks[2], kaon_mass, chi2, ndf, kaon_sigma));
  particles.push_back(kp_factory.particle(ttrks[3], kaon_mass, chi2, ndf, kaon_sigma));

  muons.push_back(kp_factory.particle(ttrks[0], muon_mass, chi2, ndf, muon_sigma));
  muons.push_back(kp_factory.particle(ttrks[1], muon_mass, chi2, ndf, muon_sigma));

  kaons.push_back(kp_factory.particle(ttrks[2], kaon_mass, chi2, ndf, kaon_sigma));
  kaons.push_back(kp_factory.particle(ttrks[3], kaon_mass, chi2, ndf, kaon_sigma));
  
  KinematicConstrainedVertexFit(particles);  

  KinematicParticleVertexFit(muons, kaons);

  KinematicParticleVertexFit(particles);
  */

}

bool KinFitTester::isAncestor(const reco::Candidate* ancestor, const reco::Candidate* particle)
{
  if ( ancestor == particle ) 
    return true;

  for ( size_t i=0; i < particle->numberOfMothers(); i++ )
  {
    if ( isAncestor(ancestor, particle->mother(i))) 
      return true;
  }
  
  return false;
}


bool KinFitTester::inKinematicAcceptanceRegion(const pat::Muon& muon)
{  
  if ( muon.pt() < ptMin_ )
    return false;
  
  if ( fabs(muon.eta()) > etaMax_ )
    return false;
  
  if ( fabs(muon.eta()) < 1.3 )
  {
    if ( muon.pt() > 3.3 ) 
      return true;
    else
      return false;
  }
  
  if ( 1.3 < fabs(muon.eta()) && fabs(muon.eta()) < 2.2 )
  {
    if ( muon.pt() > 2.9 )
      return true;
    else
      return false;
  }
  
  if ( 2.2 < fabs(muon.eta()) && fabs(muon.eta()) < 2.4 ) 
  {
    if ( muon.pt() > 0.8 )
      return true;
    else
      return false;
  }       

  return false;
}

bool KinFitTester::isGoodMuon(const pat::Muon& muon) 
{
  if ( ! muon.isLooseMuon() )
    return false;

  if ( ! muon.isTrackerMuon() )
    return false;

  if ( ! muon.innerTrack()->quality(reco::Track::highPurity) ) 
    return false;

  return true;
}

bool KinFitTester::isKalmanVertexFit(std::vector<reco::TransientTrack> tts)
{  
  KalmanVertexFitter kvf(false);
  TransientVertex tv = kvf.vertex(tts);
  
  if ( ! tv.isValid() )
    std::cout<<"KalmanVertexFitter failed"<<std::endl;
  else 
  {
    std::cout<<"KalmanVertexFit: "<<std::endl;
    std::cout<<"  position: " << reco::Vertex::Point(tv.position()) << std::endl;
    std::cout<<"  chi2pdof: "<< tv.normalisedChiSquared() <<std::endl; 
    
    return true;
  }

  return false;
}

void KinFitTester::KinematicConstrainedVertexFit(std::vector<RefCountedKinematicParticle> particles)
{
  ParticleMass JPsi_mass = 3.09687;

  MultiTrackKinematicConstraint* constraint = new TwoTrackMassKinematicConstraint(JPsi_mass);

  KinematicConstrainedVertexFitter kcv_fitter;

  try 
  {  
    RefCountedKinematicTree tree = kcv_fitter.fit(particles, constraint);
  
    std::cout<<"KinematicTree from KCVF:"<<std::endl;
    printout(tree);
  }
  catch (const std::exception& e) 
  {
    std::cout<<"exception thrown in KCVF "<< e.what() <<std::endl;
    return;
  }  
}

void KinFitTester::KinematicParticleVertexFit(std::vector<RefCountedKinematicParticle> particles)
{
  KinematicParticleVertexFitter kpv_fitter;
  RefCountedKinematicTree tree = kpv_fitter.fit(particles);
  
  std::cout<<"KinematicTree from KPVF w/o constraints"<<std::endl;
  printout(tree);
}

void KinFitTester::KinematicParticleVertexFit(std::vector<RefCountedKinematicParticle> muons,
                                              std::vector<RefCountedKinematicParticle> kaons)
{
  /*
    - Fit two final state muons to a common vertex, reconstructing J/psi parameters
      at this vertex
    - Constrain the invariant mass of two muons to be equal to J/psi mass
    - Fit the J/psi and the K to a common vertex, reconstructing the B parameters
  */

  try 
  {
    // The vertex fitter
    KinematicParticleVertexFitter kpv_fitter;
  
    RefCountedKinematicTree tree = kpv_fitter.fit(muons);

    // The particle fitter
    KinematicParticleFitter kp_fitter;

    ParticleMass JPsi_mass = 3.09687;
    KinematicConstraint* constraint = new MassKinematicConstraint(JPsi_mass, 0.00004);

    tree = kp_fitter.fit(constraint, tree);

    tree->movePointerToTheTop();
    RefCountedKinematicParticle JPsi_particle = tree->currentParticle();
    kaons.push_back(JPsi_particle); // Push this to this?

    RefCountedKinematicTree btree = kpv_fitter.fit(kaons);

    std::cout<<"KinematicTree from KPVF:"<<std::endl;
    printout(btree);
  }
  catch (const std::exception& e) 
  {
    std::cout<<"exception thrown in KPVF "<< e.what() <<std::endl;
    return;
  }
}

void KinFitTester::printout(const RefCountedKinematicVertex& vertex) const
{
  if ( ! vertex->vertexIsValid() )
  {
    std::cout<<"Vertex not valid"<<std::endl;
    return;
  } 

  std::cout<< "   Decay vertex:" 
           << vertex->position() 
           << "chi2/ndf = "
           << vertex->chiSquared() <<"/"<< vertex->degreesOfFreedom() <<std::endl;
}

void KinFitTester::printout(const RefCountedKinematicParticle& particle) const
{
  std::cout<< "   Particle:"<<std::endl;

  /*
    The 7 parameters
    - (x,y,z): reference position in the global frame
    - (px, py, pz): particle momentum at reference position
    - m: particle mass

    typedef ROOT::Math::SVector<double, 7> AlgebraicVector7

    AlgebraicVector7 particle->currentState().kinematicParameters().vector();


    The 7x7 covariance matrix
    
    AlgebraicSymMatrix77 particle->currentState().kinematicParametersError().matrix();
  */

  double x,y,z;
  double px, py, pz;
  double m;

  AlgebraicVector7 pars = particle->currentState().kinematicParameters().vector();
  
  x = pars[0];
  y = pars[1];
  z = pars[2];

  px = pars[3];
  py = pars[4];
  pz = pars[5];

  m = pars[6];

  std::cout<< "     Momentum at vertex:" << particle->currentState().globalMomentum() <<std::endl;
  std::cout<< "     Parameters at vertex: "
           <<"(x,y,z) = ("<< x <<","<< y <<","<< z <<") "
           <<"(px,py,pz) = ("<< px <<","<< py <<","<< pz <<") "
           <<"m = "<< m 
           <<std::endl;
}

void KinFitTester::printout(const RefCountedKinematicTree& tree) const
{
  if ( ! tree->isValid() ) 
  {
    std::cout<<"RefCountedKinematicTree is invalid"<<std::endl;
    return;
  }
  
  tree->movePointerToTheTop();
  // Now on top of the decay tree
  
  RefCountedKinematicParticle particle = tree->currentParticle();
  printout(particle);

  RefCountedKinematicVertex vertex = tree->currentDecayVertex();
  printout(vertex);

  std::vector<RefCountedKinematicParticle> children = tree->finalStateParticles();

  for ( unsigned int i = 0; i < children.size(); ++i ) 
  {
    printout(children[i]);
  }
  
  // Move pointer down the tree
  bool have_child = tree->movePointerToTheFirstChild();

  if ( have_child ) while ( tree->movePointerToTheNextChild() ) {
      RefCountedKinematicParticle child = tree->currentParticle();
      printout(child);
    }
}

void KinFitTester::beginJob()
{}

void KinFitTester::endJob()
{}

DEFINE_FWK_MODULE(KinFitTester);
