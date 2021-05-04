// -*- C++ -*-
//
// Package:    Demo/DemoAnalyzer
// Class:      DemoAnalyzer
//
/**\class DemoAnalyzer DemoAnalyzer.cc Demo/DemoAnalyzer/plugins/DemoAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Dener De Souza Lemos
//         Created:  Wed, 12 Feb 2020 13:05:50 GMT
//
//



#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
//v0 candidates
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/V0Candidate/interface/V0Candidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"

#include "DataFormats/TrackReco/interface/DeDxData.h"

#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"

#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"

//
// class declaration
//
#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>

//
#include <iostream>
#include <cmath>
#include <math.h>
#include <fstream>
#include "TVector3.h"
#include <vector>
#include <map>
#include "trackingEfficiency2018PbPb.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::TrackCollection;

class DemoAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit DemoAnalyzer(const edm::ParameterSet&);
      ~DemoAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      Int_t getHiBinFromhiHF(const Double_t hiHF);

      static bool vtxSort( const reco::Vertex &  a, const reco::Vertex & b );
      
	  void initHistos(const edm::Service<TFileService> & fs);  

      bool passesTrackCuts(const reco::Track & track, const reco::Vertex & vertex);
      
      // ----------member data ---------------------------
      edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
      edm::EDGetTokenT<reco::TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
	  edm::EDGetTokenT<reco::GenParticleCollection> genParticleToken_;
	  edm::EDGetTokenT<TrackingParticleCollection> tpEffSrc_;
  	  edm::EDGetTokenT<reco::SimToRecoCollection> associatorMapSTR_;
  	  edm::EDGetTokenT<CaloTowerCollection> _towerSrc;
  	  
  	  
  	  double vtx_x, vtx_y, vtx_z, vtx_rho, vtx_xError, vtx_yError, vtx_zError;  //primary vertex
  	  
TH1D* Cent;

TH2D* pT_reco;
TH2D* pT_reco_corr;
TH2D* pT_gen;
TH2D* pT_tp;
TH2D* pT_reco_onepix;

TH2D* pT_reco_pos;
TH2D* pT_reco_corr_pos;
TH2D* pT_gen_pos;
TH2D* pT_tp_pos;
   
TH2D* pT_reco_neg;
TH2D* pT_reco_corr_neg;
TH2D* pT_gen_neg;
TH2D* pT_tp_neg;

TH2D* eta_reco;
TH2D* eta_reco_corr;
TH2D* eta_gen;
TH2D* eta_tp;

TH2D* eta_reco_neg;
TH2D* eta_reco_corr_neg;
TH2D* eta_gen_neg;
TH2D* eta_tp_neg;

TH2D* eta_reco_pos;
TH2D* eta_reco_corr_pos;
TH2D* eta_gen_pos;
TH2D* eta_tp_pos;

TH2D* phi_reco;
TH2D* phi_reco_corr;
TH2D* phi_gen;
TH2D* phi_tp;

TH2D* phi_reco_neg;
TH2D* phi_reco_corr_neg;
TH2D* phi_gen_neg;
TH2D* phi_tp_neg;

TH2D* phi_reco_pos;
TH2D* phi_reco_corr_pos;
TH2D* phi_gen_pos;
TH2D* phi_tp_pos;

TrkEff2018PbPb trkEff =  TrkEff2018PbPb("general", false, "");

TrkEff2018PbPb trkEff_plus =  TrkEff2018PbPb("generalMB+", false, "");

TrkEff2018PbPb trkEff_minus =  TrkEff2018PbPb("generalMB-", false, "");

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
DemoAnalyzer::DemoAnalyzer(const edm::ParameterSet& iConfig)
 :
  vertexToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollection"))),
  tracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"))),
  genParticleToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticleTag"))),
  tpEffSrc_(consumes<TrackingParticleCollection>(iConfig.getParameter<edm::InputTag>("tpEffSrc"))),
  associatorMapSTR_(consumes<reco::SimToRecoCollection>(iConfig.getParameter<edm::InputTag>("associatorMap"))),
  _towerSrc(consumes<CaloTowerCollection>(iConfig.getParameter<edm::InputTag>("towerSrc")))
{
   //now do what ever initialization is needed
}


DemoAnalyzer::~DemoAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
DemoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace reco;
   using namespace std;

//Vertex selection

  edm::Handle<std::vector<reco::Vertex> > vertexCollection;
  iEvent.getByToken(vertexToken_,vertexCollection);
  std::vector<reco::Vertex> vtx_sorted = *vertexCollection;
//  const reco::Vertex & vtx_0 = (*vertexCollection)[0];
  std::sort( vtx_sorted.begin(), vtx_sorted.end(), DemoAnalyzer::vtxSort );
  if(vtx_sorted.size() == 0) return;
  vtx_x = (double)vtx_sorted.begin()->position().x();// VertexX_->Fill(vtx_x);
  vtx_y = (double)vtx_sorted.begin()->position().y();// VertexY_->Fill(vtx_y);
  vtx_z = (double)vtx_sorted.begin()->position().z(); 
  vtx_rho= (double)vtx_sorted.begin()->position().Rho(); 
  vtx_xError = (double)vtx_sorted.begin()->xError();
  vtx_yError = (double)vtx_sorted.begin()->yError();
  vtx_zError = (double)vtx_sorted.begin()->zError();
  math::XYZPoint vtxx(vtx_x,vtx_y,vtx_z);

   edm::Handle<reco::SimToRecoCollection > simtorecoCollectionH;
   iEvent.getByToken(associatorMapSTR_,simtorecoCollectionH);
   reco::SimToRecoCollection simRecColl;
   simRecColl= *(simtorecoCollectionH.product());

  if(vtx_z > 15.)return;

  double etHFtowerSumPlus=0;
  double etHFtowerSumMinus=0;
  double etHFtowerSum=0;
  Handle<CaloTowerCollection> towers;
  iEvent.getByToken(_towerSrc,towers);
  for( size_t i = 0; i<towers->size(); ++ i){
  const CaloTower & tower = (*towers)[ i ];
  double eta = tower.eta();
  bool isHF = tower.ietaAbs() > 29;
  if(isHF && eta > 0){
  etHFtowerSumPlus += tower.pt();
  }
  if(isHF && eta < 0){
  etHFtowerSumMinus += tower.pt();
  }
  }
  etHFtowerSum=etHFtowerSumPlus + etHFtowerSumMinus; 
  
  Int_t cent = getHiBinFromhiHF(etHFtowerSum);
  Cent->Fill(cent);
  

  
    Handle<TrackCollection> tracks;
    iEvent.getByToken(tracksToken_, tracks);
    for(TrackCollection::const_iterator iter_tk = tracks->begin(); iter_tk != tracks->end();++iter_tk) {
      double aux_tk_dz_vtx = (double)iter_tk->dz(vtxx);
      double aux_tk_dzError_vtx  = (double)sqrt(iter_tk->dzError()*iter_tk->dzError()+vtx_zError*vtx_zError);
      double aux_tk_dxy_vtx = (double)iter_tk->dxy(vtxx);
      double aux_tk_dxyError_vtx  = (double)sqrt(iter_tk->dxyError()*iter_tk->dxyError()+vtx_xError*vtx_yError);
      if(iter_tk->pt()<0.5 || iter_tk->pt()>8.0)continue;
      if(fabs(iter_tk->eta())>2.4)continue;
      if(!iter_tk->quality(reco::TrackBase::highPurity))continue;
      if(fabs(iter_tk->ptError())/iter_tk->pt()>0.1)continue;
      if(fabs(aux_tk_dz_vtx/aux_tk_dzError_vtx)>3)continue;
      if(fabs(aux_tk_dxy_vtx/aux_tk_dxyError_vtx)>3)continue;
      if(iter_tk->numberOfValidHits()<11)continue;
      if((iter_tk->normalizedChi2()/iter_tk->hitPattern().trackerLayersWithMeasurement())>0.18)continue;
      
	  pT_reco->Fill(iter_tk->pt(),cent);
	  eta_reco->Fill(iter_tk->eta(),cent);
	  phi_reco->Fill(iter_tk->phi(),cent);

	  float corr = trkEff.getCorrection(iter_tk->pt(),iter_tk->eta(),cent);
	  pT_reco_corr->Fill(iter_tk->pt(),cent,corr);
	  eta_reco_corr->Fill(iter_tk->eta(),cent,corr);
	  phi_reco_corr->Fill(iter_tk->phi(),cent,corr);	  

	  if(iter_tk->charge()>0){


	  pT_reco_pos->Fill(iter_tk->pt(),cent);
	  eta_reco_pos->Fill(iter_tk->eta(),cent);
	  phi_reco_pos->Fill(iter_tk->phi(),cent);

	  float corr_plus = trkEff_plus.getCorrection(iter_tk->pt(),iter_tk->eta(),cent);
	  
	  pT_reco_corr_pos->Fill(iter_tk->pt(),cent,corr_plus);
	  eta_reco_corr_pos->Fill(iter_tk->eta(),cent,corr_plus);
	  phi_reco_corr_pos->Fill(iter_tk->phi(),cent,corr_plus);	

	  }else if(iter_tk->charge()<0){
	  
	  pT_reco_neg->Fill(iter_tk->pt(),cent);
	  eta_reco_neg->Fill(iter_tk->eta(),cent);
	  phi_reco_neg->Fill(iter_tk->phi(),cent);
	  
	  float corr_minus = trkEff_minus.getCorrection(iter_tk->pt(),iter_tk->eta(),cent);
	  
	  pT_reco_corr_neg->Fill(iter_tk->pt(),cent,corr_minus);
	  eta_reco_corr_neg->Fill(iter_tk->eta(),cent,corr_minus);
	  phi_reco_corr_neg->Fill(iter_tk->phi(),cent,corr_minus);	
	  
	  }

      if(iter_tk->hitPattern().pixelLayersWithMeasurement()<=0)continue;
	  double corrx = trkEff.getCorrection(iter_tk->pt(),iter_tk->eta(),cent);
      pT_reco_onepix->Fill(iter_tk->pt(),cent,corrx);

    }

	edm::Handle<reco::GenParticleCollection> genpars;
	iEvent.getByToken(genParticleToken_, genpars);
	if(genpars->size() >= 1 ) {          
	  for(reco::GenParticleCollection::const_iterator ig = genpars->begin(); ig != genpars->end(); ig++){
      int status = ig->status();
      if(status != 1 || ig->charge()==0) continue; //check is target 
      if(ig->pt()<0.5 || ig->pt()>8.0)continue;   	
      if(fabs(ig->eta())>2.4)continue;  
	  pT_gen->Fill(ig->pt(),cent);
	  eta_gen->Fill(ig->eta(),cent);
	  phi_gen->Fill(ig->phi(),cent);
	  if(ig->charge()>0){
	  pT_gen_pos->Fill(ig->pt(),cent);
	  eta_gen_pos->Fill(ig->eta(),cent);
	  phi_gen_pos->Fill(ig->phi(),cent);
	  }else if(ig->charge()<0){
	  pT_gen_neg->Fill(ig->pt(),cent);
	  eta_gen_neg->Fill(ig->eta(),cent);
	  phi_gen_neg->Fill(ig->phi(),cent);
	  } 	
	}
  }
  
   // obtain collections of simulated particles 
   edm::Handle<TrackingParticleCollection>  TPCollectionHeff;
   iEvent.getByToken(tpEffSrc_,TPCollectionHeff);
if(TPCollectionHeff->size()>=1){
   for(TrackingParticleCollection::size_type i=0; i<TPCollectionHeff->size(); i++)
   {
     TrackingParticleRef tpr(TPCollectionHeff, i);
     TrackingParticle* tp=const_cast<TrackingParticle*>(tpr.get());
     
     if(tp->status() < 0 || tp->charge()==0) continue; //only charged primaries
     if(tp->pt()<0.5 || tp->pt()>8.0)continue;
     if(fabs(tp->eta())>2.4)continue;  

     std::vector<std::pair<edm::RefToBase<reco::Track>, double> > rt;
	 size_t nrec=0;


     if(simRecColl.find(tpr) != simRecColl.end())
     {
       rt = (std::vector<std::pair<edm::RefToBase<reco::Track>, double> >) simRecColl[tpr];
       std::vector<std::pair<edm::RefToBase<reco::Track>, double> >::const_iterator rtit;
       for (rtit = rt.begin(); rtit != rt.end(); ++rtit)
       {
     //    const reco::Track* tmtr = rtit->first.get();
    //     if( ! passesTrackCuts(*tmtr, vtx_sorted[0]) ) continue;
         nrec++;
       }
     }

     // if(nrec <= 0)continue;
     
	  pT_tp->Fill(tp->pt(),cent);
	  eta_tp->Fill(tp->eta(),cent);
	  phi_tp->Fill(tp->phi(),cent);
	  if(tp->charge()>0){
	  pT_tp_pos->Fill(tp->pt(),cent);
	  eta_tp_pos->Fill(tp->eta(),cent);
	  phi_tp_pos->Fill(tp->phi(),cent);
	  }else if(tp->charge()<0){
	  pT_tp_neg->Fill(tp->pt(),cent);
	  eta_tp_neg->Fill(tp->eta(),cent);
	  phi_tp_neg->Fill(tp->phi(),cent);
	  } 
    }
 }
  
}


// ------------ method called once each job just before starting event loop  ------------
void
DemoAnalyzer::beginJob()
{
  std::cout<<"  This is called once for each job: Begin Job " << std::endl;
  edm::Service<TFileService> fs;
  initHistos(fs);
}

// ------------ method called once each job just after ending the event loop  ------------
void
DemoAnalyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DemoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

Int_t DemoAnalyzer::getHiBinFromhiHF(const Double_t hiHF){
  Int_t binPos = -1;
  const Int_t nBins = 200; // table of bin edges
 const Double_t binTable[nBins+1] = {0, 12.2187, 13.0371, 13.7674, 14.5129, 15.2603, 16.0086, 16.7623, 17.5335, 18.3283, 19.1596, 19.9989, 20.8532, 21.7297, 22.6773, 23.6313, 24.6208, 25.6155, 26.6585, 27.7223, 28.8632, 30.041, 31.2865, 32.5431, 33.8655, 35.2539, 36.6912, 38.2064, 39.7876, 41.4818, 43.2416, 45.0605, 46.9652, 48.9918, 51.1, 53.2417, 55.5094, 57.9209, 60.3817, 62.9778, 65.6099, 68.4352, 71.3543, 74.4154, 77.6252, 80.8425, 84.1611, 87.7395, 91.3973, 95.1286, 99.0571, 103.185, 107.482, 111.929, 116.45, 121.178, 126.081, 130.995, 136.171, 141.612, 147.298, 153.139, 159.419, 165.633, 172.114, 178.881, 185.844, 192.845, 200.244, 207.83, 215.529, 223.489, 231.878, 240.254, 249.319, 258.303, 267.508, 277.037, 286.729, 296.845, 307.458, 317.882, 328.787, 340.074, 351.295, 362.979, 375.125, 387.197, 399.604, 412.516, 425.683, 439.001, 452.667, 466.816, 481.007, 495.679, 510.588, 526.138, 541.782, 557.641, 574.141, 591.071, 608.379, 626.068, 643.616, 661.885, 680.288, 699.449, 718.925, 738.968, 758.983, 779.459, 800.376, 821.638, 843.555, 865.771, 888.339, 911.031, 934.979, 958.56, 982.582, 1007.02, 1031.9, 1057.81, 1084.01, 1111.71, 1138.21, 1165.72, 1193.73, 1221.65, 1251.51, 1281.23, 1311.01, 1341.1, 1372.4, 1404.29, 1436.52, 1468.65, 1501.91, 1535.56, 1569.69, 1604.69, 1640.65, 1676.05, 1712.62, 1749.28, 1787.43, 1825.89, 1866.07, 1906.58, 1947.84, 1989.66, 2031.4, 2072.8, 2115.32, 2159.5, 2205.23, 2252.68, 2298.58, 2345.65, 2393.36, 2442.87, 2491.45, 2541.04, 2592.81, 2645.52, 2699.1, 2753.29, 2807.93, 2864.37, 2922.6, 2979.42, 3038.68, 3098.72, 3159.29, 3221.66, 3285.9, 3350.95, 3415.81, 3482.69, 3552.62, 3623.61, 3694.63, 3767.25, 3840.28, 3917.04, 3993.66, 4073.36, 4154.33, 4238.13, 4322.21, 4409.83, 4498.89, 4589.72, 4681.56, 4777.09, 4877.95, 4987.05, 5113.04, 5279.58, 6242.82};
  for(int i = 0; i < nBins; ++i){
    if(hiHF >= binTable[i] && hiHF < binTable[i+1]){
      binPos = i;
      break;
    }
  }

  binPos = nBins - 1 - binPos;

  return (Int_t)(200*((Double_t)binPos)/((Double_t)nBins));
}

bool DemoAnalyzer::vtxSort( const reco::Vertex &  a, const reco::Vertex & b ){
   if( a.tracksSize() != b.tracksSize() )
      return  a.tracksSize() > b.tracksSize() ? true : false ;
   else
      return  a.chi2() < b.chi2() ? true : false ;
}

void DemoAnalyzer::initHistos(const edm::Service<TFileService> & fs){

   TFileDirectory Inf = fs->mkdir( "Info" );

   Cent=Inf.make<TH1D>("Cent","",200,0.0,200.0);

   pT_reco=Inf.make<TH2D>("pT_reco","",100,0.0,10.0,200,0.0,200.0);
   pT_reco_corr=Inf.make<TH2D>("pT_reco_corr","",100,0.0,10.0,200,0.0,200.0);
   pT_gen=Inf.make<TH2D>("pT_gen","",100,0.0,10.0,200,0.0,200.0);
   pT_tp=Inf.make<TH2D>("pT_tp","",100,0.0,10.0,200,0.0,200.0);
   pT_reco_onepix=Inf.make<TH2D>("pT_reco_onepix","",100,0.0,10.0,200,0.0,200.0);

   pT_reco_pos=Inf.make<TH2D>("pT_reco_pos","",100,0.0,10.0,200,0.0,200.0);
   pT_reco_corr_pos=Inf.make<TH2D>("pT_reco_corr_pos","",100,0.0,10.0,200,0.0,200.0);
   pT_gen_pos=Inf.make<TH2D>("pT_gen_pos","",100,0.0,10.0,200,0.0,200.0);
   pT_tp_pos=Inf.make<TH2D>("pT_tp_pos","",100,0.0,10.0,200,0.0,200.0);


   pT_reco_neg=Inf.make<TH2D>("pT_reco_neg","",100,0.0,10.0,200,0.0,200.0);
   pT_reco_corr_neg=Inf.make<TH2D>("pT_reco_corr_neg","",100,0.0,10.0,200,0.0,200.0);
   pT_gen_neg=Inf.make<TH2D>("pT_gen_neg","",100,0.0,10.0,200,0.0,200.0);
   pT_tp_neg=Inf.make<TH2D>("pT_tp_neg","",100,0.0,10.0,200,0.0,200.0);

   double etabins[19]= {-2.4, -2.0, -1.6, -1.4, -1.3, -1.2, -1.0, -0.8, -0.4, 0.0, 0.4, 0.8, 1.0, 1.2, 1.3, 1.4, 1.6, 2.0, 2.4};



   eta_reco=Inf.make<TH2D>("eta_reco","",18,etabins,200,0.0,200.0);
   eta_reco_corr=Inf.make<TH2D>("eta_reco_corr","",18,etabins,200,0.0,200.0);
   eta_gen=Inf.make<TH2D>("eta_gen","",18,etabins,200,0.0,200.0);
   eta_tp=Inf.make<TH2D>("eta_tp","",18,etabins,200,0.0,200.0);

   eta_reco_neg=Inf.make<TH2D>("eta_reco_neg","",18,etabins,200,0.0,200.0);
   eta_reco_corr_neg=Inf.make<TH2D>("eta_reco_corr_neg","",18,etabins,200,0.0,200.0);
   eta_gen_neg=Inf.make<TH2D>("eta_gen_neg","",18,etabins,200,0.0,200.0);
   eta_tp_neg=Inf.make<TH2D>("eta_tp_neg","",18,etabins,200,0.0,200.0);


   eta_reco_pos=Inf.make<TH2D>("eta_reco_pos","",18,etabins,200,0.0,200.0);
   eta_reco_corr_pos=Inf.make<TH2D>("eta_reco_corr_pos","",18,etabins,200,0.0,200.0);
   eta_gen_pos=Inf.make<TH2D>("eta_gen_pos","",18,etabins,200,0.0,200.0);
   eta_tp_pos=Inf.make<TH2D>("eta_tp_pos","",18,etabins,200,0.0,200.0);


   phi_reco=Inf.make<TH2D>("phi_reco","",64,-3.2,3.2,200,0.0,200.0);
   phi_reco_corr=Inf.make<TH2D>("phi_reco_corr","",64,-3.2,3.2,200,0.0,200.0);
   phi_gen=Inf.make<TH2D>("phi_gen","",64,-3.2,3.2,200,0.0,200.0);
   phi_tp=Inf.make<TH2D>("phi_tp","",64,-3.2,3.2,200,0.0,200.0);


   phi_reco_neg=Inf.make<TH2D>("phi_reco_neg","",64,-3.2,3.2,200,0.0,200.0);
   phi_reco_corr_neg=Inf.make<TH2D>("phi_reco_corr_neg","",64,-3.2,3.2,200,0.0,200.0);
   phi_gen_neg=Inf.make<TH2D>("phi_gen_neg","",64,-3.2,3.2,200,0.0,200.0);
   phi_tp_neg=Inf.make<TH2D>("phi_tp_neg","",64,-3.2,3.2,200,0.0,200.0);


   phi_reco_pos=Inf.make<TH2D>("phi_reco_pos","",64,-3.2,3.2,200,0.0,200.0);
   phi_reco_corr_pos=Inf.make<TH2D>("phi_reco_corr_pos","",64,-3.2,3.2,200,0.0,200.0);
   phi_gen_pos=Inf.make<TH2D>("phi_gen_pos","",64,-3.2,3.2,200,0.0,200.0);
   phi_tp_pos=Inf.make<TH2D>("phi_tp_pos","",64,-3.2,3.2,200,0.0,200.0);


  
}

bool DemoAnalyzer::passesTrackCuts(const reco::Track & track, const reco::Vertex & vertex)
{
   math::XYZPoint vtxPoint(0.0,0.0,0.0);
   double vzErr =0.0, vxErr=0.0, vyErr=0.0;
   vtxPoint=vertex.position();
   vzErr=vertex.zError();
   vxErr=vertex.xError();
   vyErr=vertex.yError();

   double dxy=0.0, dz=0.0, dxysigma=0.0, dzsigma=0.0;
   dxy = track.dxy(vtxPoint);
   dz = track.dz(vtxPoint);
   dxysigma = sqrt(track.d0Error()*track.d0Error()+vxErr*vyErr);
   dzsigma = sqrt(track.dzError()*track.dzError()+vzErr*vzErr);

   double chi2n = track.normalizedChi2();
   double nlayers = track.hitPattern().trackerLayersWithMeasurement();
   chi2n = chi2n/nlayers;
   int nhits = track.numberOfValidHits();
//   int algo  = track.algo();
   if(track.pt() < 0.4)return false;   	
   if(fabs(track.eta()) > 2.4)return false;  
   if(!track.quality(reco::TrackBase::highPurity) ) return false;
   if(fabs(dxy/dxysigma) > 3.) return false;
   if(fabs(dz/dzsigma) > 3.) return false;
   if(fabs(track.ptError()) / track.pt() > 0.1) return false;
   if(nhits < 11 ) return false;
   if(chi2n >= 0.18 ) return false;
   
   return true;
}

//define this as a plug-in
DEFINE_FWK_MODULE(DemoAnalyzer);
