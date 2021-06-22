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
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include <iostream>

namespace 
{
  float muon_mass  = 0.10565837;
  float muon_sigma = 3.5e-9;

  float kaon_mass  = 0.493677;
  float kaon_sigma = 1.6e-5;

  float jpsi_mass = 3.0969;
  float jpsi_sigma = 92.9e-6;
}

class KinFitTester : public edm::one::EDAnalyzer<edm::one::SharedResources>  
{
public:
  explicit KinFitTester(const edm::ParameterSet&);
  ~KinFitTester();
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  edm::InputTag recoMuonInputTag_;
  edm::InputTag patMuonInputTag_;
  edm::InputTag candidateInputTag_;
  edm::InputTag prunedGenParticleInputTag_;
  edm::InputTag packedGenParticleInputTag_;

  bool isMC_;
  bool isAOD_;
  bool isMiniAOD_;

  double ptMin_;
  double etaMax_;

  edm::EDGetTokenT<std::vector<reco::Muon> > recoMuonToken_;
  edm::EDGetTokenT<std::vector<pat::Muon> > patMuonToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> candidateToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> prunedGenParticleToken_;
  edm::EDGetTokenT<pat::PackedGenParticleCollection> packedGenParticleToken_;

  bool isAncestor(const reco::Candidate* ancestor, const reco::Candidate* particle);
  bool inKinematicAcceptanceRegion(const pat::Muon&);
  bool isGoodMuon(const pat::Muon&);

  TransientVertex KalmanVertexFit(std::vector<reco::TransientTrack>);

  void KinematicParticleVertexFit(std::vector<RefCountedKinematicParticle>);

  void KinematicParticleVertexFit(std::vector<RefCountedKinematicParticle>,
                                  std::vector<RefCountedKinematicParticle>);

  void KinematicConstrainedVertexFit(std::vector<RefCountedKinematicParticle>);

  void printout(const RefCountedKinematicVertex& vertex) const;
  void printout(const RefCountedKinematicParticle& particle) const;
  void printout(const RefCountedKinematicTree& tree) const;

  double deltaR(const double& eta1, const double& phi1,
                const double& eta2, const double& phi2);

};

KinFitTester::KinFitTester(const edm::ParameterSet& iConfig)
  : recoMuonInputTag_(iConfig.getParameter<edm::InputTag>("recoMuonTag")),      
    patMuonInputTag_(iConfig.getParameter<edm::InputTag>("patMuonTag")),
    candidateInputTag_(iConfig.getParameter<edm::InputTag>("packedCandidateTag")),
    prunedGenParticleInputTag_(iConfig.getParameter<edm::InputTag>("prunedGenParticleTag")),
    packedGenParticleInputTag_(iConfig.getParameter<edm::InputTag>("packedGenParticleTag")),
    isMC_(iConfig.getParameter<bool>("isMC")),
    isAOD_(iConfig.getParameter<bool>("isAOD")),
    isMiniAOD_(iConfig.getParameter<bool>("isMiniAOD")),
    ptMin_(iConfig.getParameter<double>("ptMin")),
    etaMax_(iConfig.getParameter<double>("etaMax"))
{
  recoMuonToken_ = consumes<std::vector<reco::Muon> >(recoMuonInputTag_);
  patMuonToken_ = consumes<std::vector<pat::Muon> >(patMuonInputTag_);
  candidateToken_ = consumes<pat::PackedCandidateCollection>(candidateInputTag_);
  prunedGenParticleToken_ = consumes<reco::GenParticleCollection>(prunedGenParticleInputTag_);
  packedGenParticleToken_ = consumes<pat::PackedGenParticleCollection>(packedGenParticleInputTag_);
}

KinFitTester::~KinFitTester()
{}

double KinFitTester::deltaR(const double& eta1, const double& phi1,
                            const double& eta2, const double& phi2) 
{
  double deltaEta = eta1-eta2;
  double deltaPhi = phi1-phi2;
  
  if ( fabs(deltaPhi) > M_PI ) 
    deltaPhi = deltaPhi < 0 ? 2*M_PI + deltaPhi : deltaPhi - 2*M_PI; 

  return sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
}

void KinFitTester::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  /*
    First let's get the muons.
    pat::Muons for miniAOD and 
    reco::Muons for AOD.
  */
  edm::Handle<std::vector<pat::Muon> > pat_muons;
  edm::Handle<reco::MuonCollection> reco_muons;

  if ( isMiniAOD_ ) 
  {
    iEvent.getByToken(patMuonToken_, pat_muons);

    if ( ! pat_muons.isValid() ) {
      std::string error = "### Error: pat::Muons "
                          + edm::TypeID (typeid (std::vector<pat::Muon>)).friendlyClassName() + ":" 
                          + patMuonInputTag_.label() + ":"
                          + patMuonInputTag_.instance() + ":" 
                          + patMuonInputTag_.process() + " are not found.";

      std::cout<< error <<std::endl;
      return;
    }
  } 
  else if ( isAOD_ ) 
  {
    iEvent.getByToken(recoMuonToken_, reco_muons);
    return; // Don't do anything with AOD yet

    if ( ! reco_muons.isValid() )
    {  
      std::string error = "### Error: Muons "
                          + edm::TypeID (typeid (reco::MuonCollection)).friendlyClassName() + ":" 
                          + recoMuonInputTag_.label() + ":"
                          + recoMuonInputTag_.instance() + ":" 
                          + recoMuonInputTag_.process() + " are not found.";
      
      std::cout<< error <<std::endl;
      return;
    }      
  } 
  else 
  {     
    return;
  }
  
  /* 
     Get the TransientTrackBuilder.
     We will build TransientTracks from 
     reco::Tracks which in turn will be used
     to build ReferenceCountedKinematicParticles
     which will be the input to the fitters.
  */
  edm::ESHandle<TransientTrackBuilder> ttb;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", ttb);

  std::vector<reco::TransientTrack> ttrks;

  edm::Handle<pat::PackedCandidateCollection> pcands;
  iEvent.getByToken(candidateToken_, pcands);

  if ( isMC_ ) 
  {
    edm::Handle<reco::GenParticleCollection> pruned;
    iEvent.getByToken(prunedGenParticleToken_, pruned);
  
    edm::Handle<pat::PackedGenParticleCollection> packed;
    iEvent.getByToken(packedGenParticleToken_, packed);

    /*
      Find the muons and kaons from the Bs
      and fill the vector of Candidates.
    */
    std::vector<const reco::Candidate*> cands;    
    const reco::Candidate* Bs;

    for ( unsigned int i = 0; i < pruned->size(); ++i )
    {
      if ( abs((*pruned)[i].pdgId()) == 531 ) // Bs0
      {
        Bs = &(*pruned)[i];
        break;
      }
    }
    
    if ( ! Bs )
      return;
        
    for ( unsigned int j = 0; j < packed->size(); ++j )
    { 
      const reco::Candidate* motherInPrunedCollection = (*packed)[j].mother(0);
 
      if ( motherInPrunedCollection != nullptr && isAncestor(Bs, motherInPrunedCollection) )
      {
        if ( abs((*packed)[j].pdgId()) == 13 || abs((*packed)[j].pdgId()) == 321 ) // muon or K+/- 
        {
          cands.push_back((*packed)[j].clone());          
        }            
      }     
    }
       
    /*
      Muons are in vector<pat::Muon> and 
      kaons are in vector<pat::PackedCandidate>
      We know the truth from the Candidates so here
      we iterate over them and match Muons and PackedCandidates.
      When there is a match we build a TransientTrack and
      add it to our collection.
    */ 
    for ( auto cand: cands )
    {        
      if ( abs(cand->pdgId()) == 13 ) // muons
      {
        for ( std::vector<pat::Muon>::const_iterator m = pat_muons->begin(), mend = pat_muons->end();
              m != mend; ++m ) 
        {   
          if ( (m->charge() == cand->charge()) &&       
               (deltaR(m->eta(), m->phi(), cand->eta(), cand->phi()) < 0.01) &&
               (fabs((m->pt() - cand->pt()) / cand->pt()) < 0.1 )  
            ) 
          { 
            if ( m->innerTrack().isNonnull() && 
                 isGoodMuon(*m) && 
                 inKinematicAcceptanceRegion(*m) ) 
            {
              ttrks.push_back((*ttb).build(m->innerTrack().get()));
            }
          }  
        }
      } 
          
      if ( abs(cand->pdgId()) == 321 ) // kaons
      {
        for ( std::vector<pat::PackedCandidate>::const_iterator c = pcands->begin(), cend = pcands->end();
              c != cend; ++c )      
        {
          if ( (c->charge() == cand->charge()) &&
               (deltaR(c->eta(), c->phi(), cand->eta(), cand->phi()) < 0.01) &&
               (fabs((c->pt() - cand->pt()) / cand->pt()) < 0.1 )
            )
          {
            if ( c->bestTrack() != nullptr ) 
            {               
              ttrks.push_back((*ttb).build(c->bestTrack()));
            }
          }   
        }
      }
     
    } // Iterated through the candidates
  }  
  else 
  { // Not MC
   
    /*
      If data just fill all the muons
      (not realistic of course but works in the sense that things run through)
    */
    for ( std::vector<pat::Muon>::const_iterator m = pat_muons->begin(), mend = pat_muons->end();
          m != mend; ++m ) 
    {     
      if ( m->innerTrack().isNonnull() &&
           isGoodMuon(*m) && 
           inKinematicAcceptanceRegion(*m))
      { 
        ttrks.push_back((*ttb).build(m->innerTrack().get()));
      }      
    }
  }
    
  if ( ttrks.size() < 4 )
    return;

  KinematicParticleFactoryFromTransientTrack kp_factory;
  
  double chi2 = 0.0;
  double ndf  = 0.0;

  std::vector<RefCountedKinematicParticle> particles;
  std::vector<RefCountedKinematicParticle> kp_muons;
  std::vector<RefCountedKinematicParticle> kp_kaons;

  // For MC this is correct by design. Not so much for data
  
  particles.push_back(kp_factory.particle(ttrks[0], muon_mass, chi2, ndf, muon_sigma));
  particles.push_back(kp_factory.particle(ttrks[1], muon_mass, chi2, ndf, muon_sigma));
  particles.push_back(kp_factory.particle(ttrks[2], kaon_mass, chi2, ndf, kaon_sigma));
  particles.push_back(kp_factory.particle(ttrks[3], kaon_mass, chi2, ndf, kaon_sigma));
 
  kp_muons.push_back(kp_factory.particle(ttrks[0], muon_mass, chi2, ndf, muon_sigma));
  kp_muons.push_back(kp_factory.particle(ttrks[1], muon_mass, chi2, ndf, muon_sigma));
  kp_kaons.push_back(kp_factory.particle(ttrks[2], kaon_mass, chi2, ndf, kaon_sigma));
  kp_kaons.push_back(kp_factory.particle(ttrks[3], kaon_mass, chi2, ndf, kaon_sigma));
  
  /*
    The KinematicConstrainedVertexFit is the "global" strategy
  */
  std::cout<<"kcvf * KinematicConstrainedVertexFit * kcvf"<<std::endl;
  KinematicConstrainedVertexFit(particles);
  std::cout<<"kcvfkcvfkcvfkcvfkcvfkcvfkcvfkcvfkcvfkcvfkcvf"<<std::endl;
  
  /*
    The KinematicParticleVertexFit is the "sequential" strategy
  */
  std::cout<<"kpvf * KinematicParticleVertexFit * kpvf"<<std::endl;
  //KinematicConstrainedVertexFit(kp_muons);
  KinematicParticleVertexFit(kp_muons, kp_kaons); 
  //KinematicParticleVertexFit(particles);
  std::cout<<"kpvfkpvfkpvfkpvfkpvfkpvfkpvfkpvfkpvfkpvf"<<std::endl;

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

  return true;
}

bool KinFitTester::isGoodMuon(const pat::Muon& muon) 
{  
  if ( ! muon.isTrackerMuon() )
    return false;

  if ( ! muon.innerTrack()->quality(reco::Track::highPurity) ) 
    return false;

  return true;
}

TransientVertex KinFitTester::KalmanVertexFit(std::vector<reco::TransientTrack> tts)
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
  }

  return tv;
}

void KinFitTester::KinematicConstrainedVertexFit(std::vector<RefCountedKinematicParticle> particles)
{
  ParticleMass mass = jpsi_mass;
  MultiTrackKinematicConstraint* constraint = new TwoTrackMassKinematicConstraint(mass);
  KinematicConstrainedVertexFitter kcv_fitter;

  try 
  {  
    RefCountedKinematicTree tree = kcv_fitter.fit(particles, constraint);
  
    std::cout<<"  KinematicTree from KCVF:"<<std::endl;
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
  try 
  {
    KinematicParticleVertexFitter kpv_fitter;
    RefCountedKinematicTree tree = kpv_fitter.fit(particles);
  
    std::cout<<"   KinematicTree from KPVF w/o constraints"<<std::endl;
    printout(tree);
  }
  catch (const std::exception& e) 
  {
    std::cout<<"exception thrown in KPVF w/o constraints "<< e.what() <<std::endl;
    return;
  }
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
    //std::cout<<"   KinematicTre from KPVF (vertex fit):"<<std::endl;
    //printout(tree); 

    // The particle fitter
    KinematicParticleFitter kp_fitter;
    KinematicConstraint* constraint = new MassKinematicConstraint(jpsi_mass, jpsi_sigma);

    tree = kp_fitter.fit(constraint, tree);

    tree->movePointerToTheTop();
    RefCountedKinematicParticle JPsi_particle = tree->currentParticle();
    kaons.push_back(JPsi_particle); // Push this to this?

    RefCountedKinematicTree btree = kpv_fitter.fit(kaons);

    std::cout<<"  KinematicTree from KPVF: "<<std::endl;
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
