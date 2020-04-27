// -*- C++ -*-
//
// Package:    KinFitTester/KinFitTester
// Class:      KinFitTester
//
/**\class KinFitTester KinFitTester.cc KinFitTester/KinFitTester/plugins/KinFitTester.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/

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

#include <iostream>

using namespace reco;

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

  double ptMin_;
  double etaMax_;

  edm::EDGetTokenT<std::vector<pat::Muon> > patMuonToken_;

  bool isGoodMuon(const pat::Muon& muon);

  void KalmanVertexFit(std::vector<TransientTrack>);

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
    ptMin_(iConfig.getParameter<double>("ptMin")),
    etaMax_(iConfig.getParameter<double>("etaMax"))
{

  patMuonToken_ = consumes<std::vector<pat::Muon> >(patMuonInputTag_);

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

  std::vector<const reco::Track*> trks;    
  std::vector<reco::TransientTrack> ttrks;

  for ( std::vector<pat::Muon>::const_iterator m = pat_muons->begin(), mend = pat_muons->end();
        m != mend; ++m ) 
  {
    if ( isGoodMuon(*m) )
      trks.push_back(m->innerTrack().get());
  }
  
  for ( auto trk: trks )
  {
    if ( ! trk )
      continue;
    
    ttrks.push_back((*ttb).build(trk));
  }
 
  std::cout<< ttrks.size() <<" TransientTracks"<<std::endl;

  if ( ttrks.size() < 4 )
    return;

  KalmanVertexFit(ttrks);

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
  
  /*
    This scenario is not realistic, but let's proceed and get the chain running all 
    the way through.
  */
  particles.push_back(kp_factory.particle(ttrks[0], muon_mass, chi2, ndf, muon_sigma));
  particles.push_back(kp_factory.particle(ttrks[1], muon_mass, chi2, ndf, muon_sigma));
  particles.push_back(kp_factory.particle(ttrks[2], kaon_mass, chi2, ndf, kaon_sigma));
  particles.push_back(kp_factory.particle(ttrks[3], kaon_mass, chi2, ndf, kaon_sigma));

  muons.push_back(kp_factory.particle(ttrks[0], muon_mass, chi2, ndf, muon_sigma));
  muons.push_back(kp_factory.particle(ttrks[1], muon_mass, chi2, ndf, muon_sigma));

  kaons.push_back(kp_factory.particle(ttrks[2], kaon_mass, chi2, ndf, kaon_sigma));
  kaons.push_back(kp_factory.particle(ttrks[3], kaon_mass, chi2, ndf, kaon_sigma));
  
  /*
    This is the global fit
  */
  KinematicConstrainedVertexFit(particles);  

  /*
    This is the sequential fit
  */
  KinematicParticleVertexFit(muons, kaons);

}

bool KinFitTester::isGoodMuon(const pat::Muon& muon) 
{
  if ( ! muon.isLooseMuon() )
    return false;

  if ( ! muon.isTrackerMuon() )
    return false;

  if ( ! muon.innerTrack()->quality(reco::Track::highPurity) ) 
    return false;

  if ( muon.pt() < ptMin_ )
    return false;

  if ( fabs(muon.eta()) > etaMax_ )
    return false;

  return true;
}

void KinFitTester::KalmanVertexFit(std::vector<TransientTrack> transient_tracks)
{  
  KalmanVertexFitter kvf(false);
  TransientVertex tv = kvf.vertex(transient_tracks);
  
  if ( ! tv.isValid() )
    std::cout<<"KalmanVertexFitter failed"<<std::endl;
  else
    std::cout << "KalmanVertexFitter fit position: " << Vertex::Point(tv.position()) << std::endl;
}

void KinFitTester::KinematicConstrainedVertexFit(std::vector<RefCountedKinematicParticle> particles)
{
  ParticleMass JPsi_mass = 3.09687;

  MultiTrackKinematicConstraint* constraint = new TwoTrackMassKinematicConstraint(JPsi_mass);

  KinematicConstrainedVertexFitter kcv_fitter;
  
  RefCountedKinematicTree tree = kcv_fitter.fit(particles, constraint);
  
  std::cout<<"KinematicTree from KCVF:"<<std::endl;
  printout(tree);
}

void KinFitTester::KinematicParticleVertexFit(std::vector<RefCountedKinematicParticle> particles)
{
  KinematicParticleVertexFitter kpv_fitter;
  RefCountedKinematicTree tree = kpv_fitter.fit(particles);
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

void KinFitTester::printout(const RefCountedKinematicVertex& vertex) const
{
  if ( ! vertex->vertexIsValid() )
  {
    std::cout<<"Vertex not valid"<<std::endl;
    return;
  } 

  std::cout<< "    Decay vertex:" 
           << vertex->position() << vertex->chiSquared() <<" " 
           << vertex->degreesOfFreedom() <<std::endl;
}

void KinFitTester::printout(const RefCountedKinematicParticle& particle) const
{
  std::cout<< "   Particle:"<<std::endl;

  //AlgebraicVector7 par = particle->currentState().kinematicParameters().vector();
  //AlgebraicSymMatrix77 err = particle->currentState().kinematicParametersError().matrix();

  std::cout<< "     Momentum at vertex:" << particle->currentState().globalMomentum() <<std::endl;
  std::cout<< "     Parameters at vertex:" 
           << particle->currentState().kinematicParameters().vector() <<" " 
           << particle->currentState().kinematicParametersError().matrix() <<std::endl;
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
