// -*- C++ -*-
//
// Package:    PFAnalysis/PFAnalyzers
// Class:      PFTrackHFAnalyzer
//
/**\class PFTrackHFAnalyzer PFTrackHFAnalyzer.cc PFAnalysis/PFAnalyzers/plugins/PFTrackHFAnalyzer.cc

  Description: Analyzer of PFTracks/clusters in HF region
          The input step3 files will need: --outputCommand 'keep recoPFRecHits_particleFlow*_*_*','keep *_*pfTrack*_*_*'

  Implementation:
    [Notes on implementation]
*/
//
// Original Author:  Kenichi Hatakeyama
//         Created:  Wed, 17 Jul 2019 15:24:34 GMT
//
//

// system include files
#include <memory>
#include <map>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticleFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrackFwd.h"
#include "DataFormats/Common/interface/Ref.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElement.h"  
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementCluster.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
//#include "Geometry/HcalTowerAlgo/interface/HcalGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "boost/format.hpp"

#include "TH2.h"
#include "TH1.h"
#include "TProfile2D.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::TrackCollection;

bool mysort(const reco::PFCandidate& elem1, const reco::PFCandidate& elem2){ return elem1.energy() > elem2.energy() ; }

const int ETA_BINS = 82;

float etabins[ETA_BINS+1] =
  {-5.191, -4.889, -4.716, -4.538, -4.363, -4.191, -4.013, -3.839, -3.664, -3.489, -3.314, -3.139, -2.964, -2.853, -2.65,
   -2.5, -2.322, -2.172, -2.043, -1.93, -1.83, -1.74, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218, -1.131, -1.044, -0.957,
   -0.879, -0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087, 0,
   0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479,
   1.566, 1.653, 1.74, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013,
   4.191, 4.363, 4.538, 4.716, 4.889, 5.191};



class PFTrackHFAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:
    explicit PFTrackHFAnalyzer(const edm::ParameterSet&);
    ~PFTrackHFAnalyzer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    // ----------member data ---------------------------
    edm::EDGetTokenT<reco::GenParticleCollection> genparToken_; 
    edm::EDGetTokenT<CaloParticleCollection> caloparToken_; 
    edm::EDGetTokenT<reco::VertexCollection> vertexToken_; 
    edm::EDGetTokenT<std::vector<reco::PFCandidate>> pfcandToken_; 
    edm::EDGetTokenT<std::vector<reco::PFCluster>> pfclusterHFToken_; 
    edm::EDGetTokenT<std::vector<reco::PFRecHit>> pfrechitHFToken_; 
    edm::EDGetTokenT<std::vector<reco::PFRecTrack>> pftrackToken_; 
    edm::EDGetTokenT<reco::TrackCollection> trackToken_;  //used to select what tracks to read from configuration file
    edm::EDGetTokenT<HFRecHitCollection> hfrechitToken_;
    edm::EDGetTokenT< vector<PileupSummaryInfo>> muToken_;

    bool debug_;
    bool debugRecHit_;

    int nev = 0 ;
    std::vector<int> EventsToScan_;

    double ptlow_;
    double pthigh_;
    double etalow_;
    double etahigh_;
  
    TH1I * histo;
    map<TString, TH1*> m_Histos1D;
    map<TString, TH2*> m_Histos2D;

    map<TString, TProfile*> m_Profiles;
    map<TString, TProfile2D*> m_Profiles2D;
    void FillHist1D(const TString& histName, const Double_t& value, const double& weight);
    void FillHist2D(const TString& histName, const Double_t& value1, const Double_t& value2, const double& weight);
    void FillProfile(const TString& histName, const Double_t& value1, const Double_t& value2, const double& weight);
    void FillProfile2D(const TString& histName, const Double_t& value1, const Double_t& value2, const Double_t& value3, const double& weight);
    bool IsInCell(double eta, double phi, const CaloCellGeometry::CornersVec& CV);  
    int idphi(int iphiA, int iphiB);        
    int getEtaIndex(float eta);   
    enum Flavor{chm = 0, chu, nh, ne, hfh, hfe, lep, untrk, numFlavors, X };
    TString ids[8] = {"chm", "chu", "nh", "ne", "hfh", "hfe", "lep", "untrk"};
    Flavor getFlavor(reco::PFCandidate::ParticleType id);
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
PFTrackHFAnalyzer::PFTrackHFAnalyzer(const edm::ParameterSet& iConfig)
  :
  genparToken_(consumes<reco::GenParticleCollection>(iConfig.getUntrackedParameter<edm::InputTag>("source_genpars"))),
  caloparToken_(consumes<CaloParticleCollection>(iConfig.getUntrackedParameter<edm::InputTag>("source_calopars"))),
  vertexToken_(consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("source_vertices"))),
  pfcandToken_(consumes<std::vector<reco::PFCandidate>>(iConfig.getUntrackedParameter<edm::InputTag>("source_pfcands"))),
  pfclusterHFToken_(consumes<std::vector<reco::PFCluster>>(iConfig.getUntrackedParameter<edm::InputTag>("source_pfclustersHF"))),
  pfrechitHFToken_(consumes<std::vector<reco::PFRecHit>>(iConfig.getUntrackedParameter<edm::InputTag>("source_pfrechitsHF"))),
  pftrackToken_(consumes<std::vector<reco::PFRecTrack>>(iConfig.getUntrackedParameter<edm::InputTag>("source_pftracks"))),
  trackToken_(consumes<reco::TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("source_tracks"))),
  hfrechitToken_(consumes<HFRecHitCollection>(iConfig.getUntrackedParameter<edm::InputTag>("source_hfrechits"))),
  muToken_(consumes<std::vector<PileupSummaryInfo>>(iConfig.getUntrackedParameter<edm::InputTag>("source_pileup"))),

  debug_(iConfig.getUntrackedParameter<bool>("debug")),
  debugRecHit_(iConfig.getUntrackedParameter<bool>("debugRecHit")),
  EventsToScan_(iConfig.getParameter<std::vector<int>>("EventsToScan")),
  ptlow_(iConfig.getParameter<double>("ptlow")),
  pthigh_(iConfig.getParameter<double>("pthigh")),
  etalow_(iConfig.getParameter<double>("etalow")),
  etahigh_(iConfig.getParameter<double>("etahigh"))

{
  //now do what ever initialization is needed
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  TString hname;

  hname = "gen_pt";   
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 250 , 0. , 250. );
  hname = "gen_eta";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 120 , 0. , 6. );

  hname = "mu";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , -0.5 , 199.5 );
  hname = "mua";  
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , -0.5 , 199.5 );

  hname = "pv_n";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , -0.5 , 199.5 );
  hname = "pv_z";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 60 , -30. , 30. );

  hname = "pftrackHF_n";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 21 , -0.5 , 20.5 );
  hname = "pftrackHF_n_range";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 300 , 0. , 3000. );
  hname = "pfcandCH1ele_pftrackHF_n";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 21 , -0.5 , 20.5 );


  hname = "pftrack_pt";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 250 , 0. , 250. );
  hname = "pftrack_eta";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 120 , -6. , 6. );
  hname = "pftrack_pull";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 100 , -50. , 50. );
//  hname = "pftrackHF_pt";
//  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 100 , 0. , 20. );
  hname = "pfPUtrack_pt";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 100 , 0. , 20. );
  hname = "pfPUtrackHF_pt";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 100 , 0. , 20. );
  hname = "pfPUtrack_eta";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 120 , -6. , 6. );
  hname = "pfPUtrack_etapt";
  m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 120 , -6. , 6.,  100 , 0. , 20. );


  hname = "DietaDiphiHFEM";
  m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 37, -18.5, 18.5, 37 , -18.5 , 18.5 );
  hname = "DietaDiphiHFHAD"; 
  m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 37, -18.5, 18.5, 37 , -18.5 , 18.5 );

  hname = "pfrechitHFEM_n";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 500 , -0.5 , 999.5 );
  hname = "pfrechitHFHAD_n";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 500 , -0.5 , 999.5 );
  hname = "pfrechitHFEM_occupancy";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 100 , 0. , 1. );
  hname = "pfrechitHFHAD_occupancy";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 100 , 0. , 1. );
  hname = "pfrechitHFEM_occupancy_vs_absieta";
  m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 21, 24.5, 45.5, 100 , 0. , 1. );
  hname = "pfrechitHFHAD_occupancy_vs_absieta";
  m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 21, 24.5, 45.5, 100 , 0. , 1. );


  hname = "pfrechitHFEM_E";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , 0. , 200. );
  hname = "pfrechitHFHAD_E";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , 0. , 200. );

  hname = "pfrechitHFEM_Emax";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 205 , -5. , 200. );
  hname = "pfrechitHFHAD_Emax";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 205 , -5. , 200. );
  hname = "pfrechitHF_Emax";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 205 , -5. , 200. );

  hname = "GProfileHFEMP_E";
  m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 37, -18.5, 18.5, 37 , -18.5 , 18.5 );
  hname = "GProfileHFHADP_E";
  m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 37, -18.5, 18.5, 37 , -18.5 , 18.5 );
  hname = "GProfileHFEMN_E";
  m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 37, -18.5, 18.5, 37 , -18.5 , 18.5 );
  hname = "GProfileHFHADN_E";
  m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 37, -18.5, 18.5, 37 , -18.5 , 18.5 );
  hname = "TProfileHFEMP_E";
  m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 37, -18.5, 18.5, 37 , -18.5 , 18.5 );
  hname = "TProfileHFHADP_E";
  m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 37, -18.5, 18.5, 37 , -18.5 , 18.5 );
  hname = "TProfileHFEMN_E";
  m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 37, -18.5, 18.5, 37 , -18.5 , 18.5 );
  hname = "TProfileHFHADN_E";
  m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 37, -18.5, 18.5, 37 , -18.5 , 18.5 );

  TFileDirectory subDir = fs->mkdir( "ScanEvents" );
  for (unsigned int i=0; i < EventsToScan_.size(); i++){
    hname = Form("GenEMP%i", EventsToScan_[i]) ;
    m_Histos2D[hname] = subDir.make<TH2F>(hname, hname, 15, 27.5, 42.5,  75 , -0.5 , 74.5 );
    hname = Form("GenHADP%i", EventsToScan_[i]) ;
    m_Histos2D[hname] = subDir.make<TH2F>(hname, hname, 15, 27.5, 42.5,  75 , -0.5 , 74.5 );
    hname = Form("ClustersEMP%i", EventsToScan_[i]) ;
    m_Histos2D[hname] = subDir.make<TH2F>(hname, hname, 15, 27.5, 42.5,  75 , -0.5 , 74.5 );
    hname = Form("ClustersHADP%i", EventsToScan_[i]) ;
    m_Histos2D[hname] = subDir.make<TH2F>(hname, hname, 15, 27.5, 42.5,  75 , -0.5 , 74.5 );

    hname = Form("hCandidatesEMP%i", EventsToScan_[i]) ;
    m_Histos2D[hname] = subDir.make<TH2F>(hname, hname, 15, 27.5, 42.5,  75 , -0.5 , 74.5 );
    hname = Form("hCandidatesHADP%i", EventsToScan_[i]) ;
    m_Histos2D[hname] = subDir.make<TH2F>(hname, hname, 15, 27.5, 42.5,  75 , -0.5 , 74.5 );
    hname = Form("ERhCandidatesEMP%i", EventsToScan_[i]) ;
    m_Histos2D[hname] = subDir.make<TH2F>(hname, hname, 15, 27.5, 42.5,  75 , -0.5 , 74.5 );
    hname = Form("ERhCandidatesHADP%i", EventsToScan_[i]) ;
    m_Histos2D[hname] = subDir.make<TH2F>(hname, hname, 15, 27.5, 42.5,  75 , -0.5 , 74.5 );


    hname = Form("EClustersEMP%i", EventsToScan_[i]) ;
    m_Histos2D[hname] = subDir.make<TH2F>(hname, hname, 15, 27.5, 42.5,  75 , -0.5 , 74.5 );
    hname = Form("EClustersHADP%i", EventsToScan_[i]) ;
    m_Histos2D[hname] = subDir.make<TH2F>(hname, hname, 15, 27.5, 42.5,  75 , -0.5 , 74.5 );
    hname = Form("ERClustersEMP%i", EventsToScan_[i]) ;
    m_Histos2D[hname] = subDir.make<TH2F>(hname, hname, 15, 27.5, 42.5,  75 , -0.5 , 74.5 );
    hname = Form("ERClustersHADP%i", EventsToScan_[i]) ;
    m_Histos2D[hname] = subDir.make<TH2F>(hname, hname, 15, 27.5, 42.5,  75 , -0.5 , 74.5 );
  }
  
  hname = "pfclusHFEM_nhits";   
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 12 , -1.5 , 10.5 );
  hname = "pfclusHFEM_nhits_threshold";   
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 12 , -1.5 , 10.5 );
  hname = "pfclusHFHAD_nhits";  
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 12 , -1.5 , 10.5 );

  hname = "pfclusHFEM_n";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 500 , -0.5 , 999.5 );
  hname = "pfclusHFHAD_n";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 500 , -0.5 , 999.5 );

  hname = "pfclusHFEM_E";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , 0. , 200. );
  hname = "pfclusHFHAD_E";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , 0. , 200. );

  hname = "pfclusHFEM_pt";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , 0. , 200. );
  hname = "pfclusHFHAD_pt";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , 0. , 200. );

  hname = "pfclusHFEM_Emax";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 205 , -5. , 200. );
  hname = "pfclusHFHAD_Emax";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 205 , -5. , 200. ); 
  hname = "pfclusHF_Emax";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 205 , -5. , 200. ); 

  hname = "pfclusHFEM_ptmax";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 205 , -5. , 200. );
  hname = "pfclusHFHAD_ptmax";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 205 , -5. , 200. );
  hname = "pfclusHF_ptmax";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 205 , -5. , 200. );

  hname = "pfclusHF_nmatch";
  m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 22 , -1.5 , 20.5, 22 , -1.5 , 20.5);
  hname = "pfclusHF_nmatch5";
  m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 22 , -1.5 , 20.5, 22 , -1.5 , 20.5 );  
  hname = "pfclusHF_nmatch9";
  m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 22 , -1.5 , 20.5, 22 , -1.5 , 20.5 );

  hname = "match_ptfrac";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 100, 0. , 2.0 );   
  hname = "match_ptfrac5";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 100, 0. , 2.0 );
  hname = "match_ptfrac9";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 100, 0. , 2.0 );

  hname = "blocksHFEM_n";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 6 , -0.5 , 5.5 );
  hname = "blocksHFHAD_n";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 6 , -0.5 , 5.5 );
  hname = "blocksHFCH_n";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 6 , -0.5 , 5.5 );

  hname = "pfcandHFEM_n";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 500 , -0.5 , 499.5 ); 
  hname = "pfcandHFHAD_n";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 500 , -0.5 , 499.5 );
  hname = "pfcandHFCH_n";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 1001 , -0.5 , 1000.5 );
  hname = "PFcand_vs_PFtrack";
  m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 1001 , -0.5 , 1000.5 , 1001 , -0.5 , 1000.5  ); 

  hname = "pfcandHFEM_pt";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , 0. , 100. );
  hname = "pfcandHFHAD_pt";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , 0. , 100. );
  hname = "pfcandHFCH_pt";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , 0. , 100. );
  hname = "pfcandNoHFCH_pt";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , 0. , 100. );
  hname = "pfcandNoHFCH_eta";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 120 , -6. , 6. );
  hname = "pfcandNoHFCHHF_pt";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , 0. , 100. );
  hname = "pfcandNoHFCHHF_eta";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 120 , -6. , 6. );


  hname = "pfcandHFEM_E";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 150 , 0. , 1500. );
  hname = "pfcandHFHAD_E";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 150 , 0. , 1500. );
  hname = "pfcandHFCH_E";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 150 , 0. , 1500. );


  hname = "pfcandHFEM_nelements" ;
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 201 , -0.5, 200.5 );
  hname = "pfcandHFHAD_nelements" ;
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 201 , -0.5, 200.5 );

// ....
  hname = "pfcandCH_nelements" ;
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 2001 , -0.5, 2000.5 );

  hname = "pfcandCH_nelementsHF" ;
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 2001 , -0.5, 2000.5 );

  hname = "pfcandCH_nelementsHFEM" ;
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 2001 , -0.5, 2000.5 );

  hname = "pfcandCH_nelementsHFHAD" ;
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 2001 , -0.5, 2000.5 );

  hname = "pfcandCH_brem_nelements" ;
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 2001 , -0.5, 2000.5 );

  hname = "pfcandCH_nobrem_nelements" ;
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 2001 , -0.5, 2000.5 );

  hname = "pfcandCHPU_nelements";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 2001 , -0.5, 2000.5 );

  hname = "pfcandCHPU_nelementsHF";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 2001 , -0.5, 2000.5 );

  hname = "pfcandCHPU_nelementsHFHAD";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 2001 , -0.5, 2000.5 );

  hname = "pfcandCHPU_nelementsHFEM";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 2001 , -0.5, 2000.5 );
  
  hname = "pfcandCHPUinHF_nelements";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 201 , -0.5, 200.5 ); 

  hname = "pfcandCHPUinHF_pt_nelements";
  m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 200 , 0. , 10.,  2001 , -0.5, 2000.5  );

// ....

  hname = "pfcandHFEM_nhits" ;
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 41 , -0.5, 40.5 );
  hname = "pfcandHFHAD_nhits" ;
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 41 , -0.5, 40.5 );
  hname = "pfcandCH_nhits" ;
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 41 , -0.5, 40.5 );
  hname = "pfcandCHPUinHF_nhits" ;
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 41 , -0.5, 40.5 );


  hname = "pfcandHFEM_n9" ;
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 16 , -0.5, 15.5 );
  hname = "pfcandHFHAD_n9" ;
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 16 , -0.5, 15.5 );  
  hname = "pfcandHFCH_n9" ;
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 16 , -0.5, 15.5 );

  hname = "Eratio_PFcand1toGen";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , 0. , 2. );
  hname = "Eratio_PFcand2toGen";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , 0. , 2. );
  hname = "Eratio_PFcand3toGen";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , 0. , 2. );

  hname = "Eratio_PFcandAlltoGen";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , 0. , 2. );
  hname = "Eratio_PFcandNonTtoGen";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , 0. , 2. );

  hname = "PovGenE_brem";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , 0. , 2. );
  hname = "PovGenE_nobrem";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , 0. , 2. );

  hname = "EovP";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 100 , 0. , 5. );
  hname = "EovP_brem";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 100 , 0. , 5. );
  hname = "EovP_nobrem";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 100 , 0. , 5. );

  hname = "EovPPU_brem";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 100 , 0. , 5. );
  hname = "EovPPU_nobrem";
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 100 , 0. , 5. );

  hname = "HFCH_DR_track_HFclusters" ;
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , 0. , 4. );

  hname = "pfcandCH1ele_drmin" ;
  m_Histos1D[hname] = fs->make<TH1F>(hname, hname , 200 , 0. , 4. );

  TFileDirectory subDir_offset = fs->mkdir( "OffsetPlots" );
  for (int i_id=0; i_id<numFlavors; i_id++){
    for (int i_nPU=0; i_nPU<25; i_nPU++){
      hname = Form("edensity_eta_nPU%i_",i_nPU) + ids[i_id];
      m_Profiles[hname] = subDir_offset.make<TProfile>(hname, hname, ETA_BINS, etabins);
      hname = Form("ptdensity_eta_nPU%i_",i_nPU) + ids[i_id];
      m_Profiles[hname] = subDir_offset.make<TProfile>(hname, hname, ETA_BINS, etabins);
    }
  }
}


PFTrackHFAnalyzer::~PFTrackHFAnalyzer()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
PFTrackHFAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  edm::Handle<reco::GenParticleCollection> genpars; iEvent.getByToken(genparToken_, genpars);
  edm::Handle<CaloParticleCollection> calopars; iEvent.getByToken(caloparToken_, calopars);
  edm::Handle<reco::VertexCollection> vertices; iEvent.getByToken(vertexToken_, vertices);
  edm::Handle<std::vector<reco::PFCandidate>> pfcands; iEvent.getByToken(pfcandToken_, pfcands);
  edm::Handle<std::vector<reco::PFCluster>> pfclustersHF; iEvent.getByToken(pfclusterHFToken_, pfclustersHF);
  edm::Handle<std::vector<reco::PFRecHit>> pfrechitsHF; iEvent.getByToken(pfrechitHFToken_, pfrechitsHF);
  edm::Handle<std::vector<reco::PFRecTrack>> pftracks; iEvent.getByToken(pftrackToken_, pftracks);
  edm::Handle<reco::TrackCollection> tracks; iEvent.getByToken(trackToken_, tracks);
  edm::Handle<HFRecHitCollection> hfRecHits;   iEvent.getByToken(hfrechitToken_, hfRecHits);
  edm::Handle<std::vector<PileupSummaryInfo>> pileups;  iEvent.getByToken(muToken_, pileups);

  TString hname;

  nev++ ;
  std::vector<int>::iterator it;
  bool scan = false ;
  it = std::find (EventsToScan_.begin(), EventsToScan_.end(), nev);
  if (it != EventsToScan_.end()){
    scan = true ;
    if (debug_)
    LogPrint("PFTrackHFAnalyzer") << boost::format(" Will scan this event  %i") % EventsToScan_[it-EventsToScan_.begin()] ;
  }

  if (debug_ && scan) 
  LogPrint("PFTrackHFAnalyzer") << "\n\n ======== genpars: ======== "  << genpars->size();
  double genP_pt = 0., genP_E = 0., genP_eta = 0., genP_phi = 0.,
         genN_pt = 0., genN_E = 0., genN_eta = 0., genN_phi = 0. ;
  int    genP_ieta1 = -1, genP_ieta2 = -1, genP_iphi1 = -1, genP_iphi2 = -1,
         genN_ieta1 = -1, genN_ieta2 = -1, genN_iphi1 = -1, genN_iphi2 = -1;
  for(const auto& genpar : *(genpars.product()) ){
    if (debug_ && scan)
    LogPrint("PFTrackHFAnalyzer") << boost::format("genpar (pt,eta,phi,E): (%6.2f, %6.2f, %6.2f, %6.2f)") % genpar.pt() % genpar.eta() % genpar.phi() % genpar.energy();
    if(genpar.eta()>0.){
      genP_pt  = genpar.pt();
      genP_E   = genpar.energy();
      genP_eta = genpar.eta();
      genP_phi = genpar.phi();
    } else {
      genN_pt  = genpar.pt();
      genN_E   = genpar.energy();
      genN_eta = genpar.eta();
      genN_phi = genpar.phi();
    }
  }
  if(genP_pt < ptlow_ || genP_pt > pthigh_ || fabs(genP_eta) < etalow_ || fabs(genP_eta) > etahigh_) return ;
  FillHist1D("gen_pt",  genP_pt,  1. );
  FillHist1D("gen_eta", genP_eta,  1. );


  if (debug_ && scan)
  LogPrint("PFTrackHFAnalyzer") << "\n =========== calopars: =========== "     << calopars->size();
  for(const auto& calopar : *(calopars.product()) ){
    if (debug_ && scan)
    LogPrint("PFTrackHFAnalyzer") << boost::format("calopar (pt,eta,phi): (%6.1f, %6.2f, %6.2f)") % calopar.pt() % calopar.eta() % calopar.phi();
  }

 
   double mu = 0. ;
   std::vector<PileupSummaryInfo>::const_iterator pileupinfo;
   for(pileupinfo = pileups->begin(); pileupinfo != pileups->end(); ++pileupinfo){
     if(pileupinfo->getBunchCrossing() == 0 ){
       mu = pileupinfo->getTrueNumInteractions() ;
       FillHist1D("mua", pileupinfo->getPU_NumInteractions(),  1. );
     }
//     cout << "BX = "<< pileupinfo->getBunchCrossing()<<endl;
//     cout << "mu = "<<pileupinfo->getTrueNumInteractions() <<endl;
//     cout << "muactual = "<< pileupinfo->getPU_NumInteractions()<<endl;
    }
  FillHist1D("mu", mu,  1. );


  if (debug_ && scan)
  LogPrint("PFTrackHFAnalyzer") << "\n =========== vertices: =========== "     << vertices->size();
  for(const auto& pv : *(vertices.product()) ){
    FillHist1D("pv_z", pv.z(),  1. );
  } 
  FillHist1D("pv_n", vertices->size(),  1. );


  if (debug_ && scan)
  LogPrint("PFTrackHFAnalyzer") << "\n =========== pftracks: =========== "     << pftracks->size();
  int pftrack_n = 0,  pftrack_npos = 0, pftrack_nneg = 0, pftrackHF_n = 0;
  double trkP_P = -1., trkP_pt  = -1., trkP_pterror = -1., trkP_eta = 0., trkP_phi = 0., trkP_etaHF = 0., trkP_phiHF = 0.,
         trkN_P = -1., trkN_pt  = -1., trkN_pterror = -1., trkN_eta = 0., trkN_phi = 0., trkN_etaHF = 0., trkN_phiHF = 0.;
  int    trkP_ieta1 = -1, trkP_ieta2 = -1, trkP_iphi1 = -1, trkP_iphi2 = -1,
         trkN_ieta1 = -1, trkN_ieta2 = -1, trkN_iphi1 = -1, trkN_iphi2 = -1;

  reco::TrackRef trackrefP, trackrefN ;
  for(const auto& genpar : *(genpars.product()) ){
    double deltarMin = 100000.; 
    for(const auto& pftrack : *(pftracks.product()) ){  
      const reco::TrackRef trackref = pftrack.trackRef();
      double deltar = deltaR(genpar, (*trackref));
      if(deltar < deltarMin){
        deltarMin = deltar ;
        constexpr reco::PFTrajectoryPoint::LayerType VFcalEntrance = reco::PFTrajectoryPoint::VFcalEntrance;
        const reco::PFTrajectoryPoint& tkAtHF = pftrack.extrapolatedPoint( VFcalEntrance );
        if(genpar.eta() > 0){
          trkP_P   = trackref->p(); 
          trkP_pt  = trackref->pt(); 
          trkP_pterror = trackref->ptError();
          trkP_eta = trackref->eta();
          //trkP_phi = trackref->phi();
          trkP_etaHF = tkAtHF.positionREP().Eta();
          trkP_phiHF = tkAtHF.positionREP().Phi();
          trackrefP = trackref ;
        } else {
          trkN_P   = trackref->p();
          trkN_pt  = trackref->pt();
          trkN_pterror = trackref->ptError();
          trkN_eta = trackref->eta();
          //trkN_phi = trackref->phi();
          trkN_etaHF = tkAtHF.positionREP().Eta();
          trkN_phiHF = tkAtHF.positionREP().Phi();
          trackrefN = trackref ;
        }
      }
    }
  }
  FillHist1D("pftrack_pt",  trkP_pt,  1. );
  FillHist1D("pftrack_eta", trkP_eta, 1. );
  FillHist1D("pftrack_pull",  ((trkP_pt-genP_pt)/trkP_pterror),  1. );
  FillHist1D("pftrack_pt",  trkN_pt,  1. ); 
  FillHist1D("pftrack_eta", trkN_eta,  1. );
  FillHist1D("pftrack_pull",  ((trkN_pt-genN_pt)/trkN_pterror),  1. );

  for(const auto& pftrack : *(pftracks.product()) ){
    pftrack_n++; 
    const reco::TrackRef trackref = pftrack.trackRef();
    if(trackref->eta() > 0.)   pftrack_npos++ ;
    if(trackref->eta() < 0.)   pftrack_nneg++ ;
    if(fabs(trackref->eta()) > 3.)   pftrackHF_n++ ;
    if(trackref != trackrefP && trackref != trackrefN){
      FillHist1D("pfPUtrack_pt",  trackref->pt(),  1. );
      FillHist1D("pfPUtrack_eta", trackref->eta(),  1. );
      FillHist2D("pfPUtrack_etapt", trackref->eta(), trackref->pt(),  1. );
      if(fabs(trackref->eta()) > 3.) FillHist1D("pfPUtrackHF_pt",  trackref->pt(),  1. );
    } 

    if (debug_ && scan){
      // looping over all points on the trajectory
      std::vector<reco::PFTrajectoryPoint> trajectoryPoints = pftrack.trajectoryPoints();
      for (unsigned ip = 0; ip<trajectoryPoints.size(); ++ip){
      //ctmp  if (trajectoryPoints[ip].isValid()) cout << trajectoryPoints[ip] <<endl;
        if (fabs(trackref->eta()) < 1.0 && trackref->pt() < 3.0 ) cout << trajectoryPoints[ip] <<endl;
      }
    }

    constexpr reco::PFTrajectoryPoint::LayerType VFcalEntrance =
    reco::PFTrajectoryPoint::VFcalEntrance;
    
    const reco::PFTrajectoryPoint& tkAtHF =
    pftrack.extrapolatedPoint( VFcalEntrance );
    
    const double tracketa = tkAtHF.positionREP().Eta();
    const double trackphi = tkAtHF.positionREP().Phi();

    if (debug_ && scan)
    LogPrint("PFTrackHFAnalyzer") << boost::format("pftrack (pt,eta,phi)@origin (eta,phi)@HF: (%6.1f +- %4.1f, %6.2f, %6.2f) (%6.2f, %6.2f)")
    % trackref->pt() % trackref->ptError() % trackref->eta() % trackref->phi() % tracketa % trackphi;

  }
  if (debug_)
  if(pftrack_n> 2) LogPrint("PFTrackHFAnalyzer") << boost::format(" Warning!! More than two tracks. Scan this event # %i") % nev  ;
  FillHist1D("pftrackHF_n",  pftrackHF_n,  1. );
  FillHist1D("pftrackHF_n_range",  pftrackHF_n,  1. );


  if (debug_ && scan)
  LogPrint("PFTrackHFAnalyzer") << "\n  =========== HFRecHits: ===========  " <<  hfRecHits->size();;

  edm::ESHandle<CaloGeometry> geoHandle;
  iSetup.get<CaloGeometryRecord>().get(geoHandle);
  const CaloSubdetectorGeometry* HFGeom = geoHandle->getSubdetectorGeometry(DetId::Hcal, 4);

  for(const auto& rechit : *(hfRecHits.product()) ){
    int ieta = rechit.id().ieta() ;
    int iphi = rechit.id().iphi() ;
    int idep = rechit.id().depth() ;
    double E = rechit.energy() ;
    if (debugRecHit_ && nev==1)
    LogPrint("PFTrackHFAnalyzer") << boost::format(" hfrechit (ieta, iphi, depth, E): (%3d, %3d, %2d, %6.2f)") % ieta % iphi % idep % E ;

    auto thisCell = HFGeom->getGeometry(rechit.id().rawId());
    const CaloCellGeometry::CornersVec& cv = thisCell->getCorners();
    if (debugRecHit_ && nev==1)
    LogPrint("PFTrackHFAnalyzer") << boost::format(" corners (eta, eta, eta, eta, phi, phi, phi, phi): (%6.2f, %6.2f, %6.2f, %6.2f, %6.2f, %6.2f, %6.2f, %6.2f)")
        % cv[0].eta() % cv[1].eta() % cv[2].eta() % cv[3].eta() % cv[0].phi() % cv[1].phi() % cv[2].phi() % cv[3].phi()  ;
    
    bool incell ;
    incell = IsInCell(genP_eta, genP_phi, cv) ;  
    if(incell) {
      if(idep == 1){
        genP_ieta1 = ieta ; genP_iphi1 = iphi ;
      } else {
        genP_ieta2 = ieta ; genP_iphi2 = iphi ;
      }
    } 
    incell = IsInCell(genN_eta, genN_phi, cv) ;
    if(incell) {
      if(idep == 1){ 
        genN_ieta1 = ieta ; genN_iphi1 = iphi ;
      } else {
        genN_ieta2 = ieta ; genN_iphi2 = iphi ;
      }
    }
    incell = IsInCell(trkP_etaHF, trkP_phiHF, cv) ;
    if(incell) {
      if(idep == 1){
        trkP_ieta1 = ieta ; trkP_iphi1 = iphi ;
      } else {
        trkP_ieta2 = ieta ; trkP_iphi2 = iphi ;
      }
    }
    incell = IsInCell(trkN_etaHF, trkN_phiHF, cv) ;
    if(incell) {
      if(idep == 1){
        trkN_ieta1 = ieta ; trkN_iphi1 = iphi ;
      } else {
        trkN_ieta2 = ieta ; trkN_iphi2 = iphi ;
      }
    }  
  } // end of loop over HFRecHits ...
  FillHist2D("DietaDiphiHFEM" , trkP_ieta1-genP_ieta1 , trkP_iphi1-genP_iphi1  , 1.) ;
  FillHist2D("DietaDiphiHFEM" , trkP_ieta1-genN_ieta1 , trkP_iphi1-genN_iphi1  , 1.) ;
  FillHist2D("DietaDiphiHFHAD", trkP_ieta2-genP_ieta2 , trkP_iphi2-genP_iphi2  , 1.) ;
  FillHist2D("DietaDiphiHFHAD", trkP_ieta2-genN_ieta2 , trkP_iphi2-genN_iphi2  , 1.) ;

  if (scan){ 
    hname = Form("GenEMP%i", nev);
    FillHist2D(hname, genP_ieta1,genP_iphi1,genP_E) ;
    hname = Form("GenHADP%i", nev);
    FillHist2D(hname, genP_ieta2,genP_iphi2,genP_E) ;
  }

  if (debug_ && scan){
    LogPrint("PFTrackHFAnalyzer") << boost::format("\nParticle in positive eta with (eta phi): (%6.2f, %6.2f)") % genP_eta % genP_phi;
    LogPrint("PFTrackHFAnalyzer") << boost::format("   Depth1 (ieta iphi): (%3d, %3d)") % genP_ieta1 % genP_iphi1  ;
    LogPrint("PFTrackHFAnalyzer") << boost::format("   Depth2 (ieta iphi): (%3d, %3d)") % genP_ieta2 % genP_iphi2  ;
    LogPrint("PFTrackHFAnalyzer") << boost::format(" Particle in negative eta with (eta phi): (%6.2f, %6.2f)") % genN_eta % genN_phi ;
    LogPrint("PFTrackHFAnalyzer") << boost::format("   Depth1 (ieta iphi): (%3d, %3d)") % genN_ieta1 % genN_iphi1  ;
    LogPrint("PFTrackHFAnalyzer") << boost::format("   Depth2 (ieta iphi): (%3d, %3d)") % genN_ieta2 % genN_iphi2  ;

    LogPrint("PFTrackHFAnalyzer") << boost::format("\nTrack in positive eta with (eta phi): (%6.2f, %6.2f)") % trkP_etaHF % trkP_phiHF;
    LogPrint("PFTrackHFAnalyzer") << boost::format("   Depth1 (ieta iphi): (%3d, %3d)") % trkP_ieta1 % trkP_iphi1  ;
    LogPrint("PFTrackHFAnalyzer") << boost::format("   Depth2 (ieta iphi): (%3d, %3d)") % trkP_ieta2 % trkP_iphi2  ;
    LogPrint("PFTrackHFAnalyzer") << boost::format(" Track in negative eta with (eta phi): (%6.2f, %6.2f)") % trkN_etaHF % trkN_phiHF ;
    LogPrint("PFTrackHFAnalyzer") << boost::format("   Depth1 (ieta iphi): (%3d, %3d)") % trkN_ieta1 % trkN_iphi1  ;
    LogPrint("PFTrackHFAnalyzer") << boost::format("   Depth2 (ieta iphi): (%3d, %3d)") % trkN_ieta2 % trkN_iphi2  ;
  }


  if (debug_ && scan)
  LogPrint("PFTrackHFAnalyzer") << "\n =========== pfrechitsHF: =========== "  << pfrechitsHF->size();
  unsigned pfrechitHFEM_n = 0, pfrechitHFHAD_n = 0 ;
  map<unsigned, unsigned> pfrechitHFEM_n_ieta, pfrechitHFHAD_n_ieta ;
  for (unsigned i=28; i<=41 ; ++i){
    pfrechitHFEM_n_ieta.insert(make_pair(i,0));
    pfrechitHFHAD_n_ieta.insert(make_pair(i,0));
  }
  double pfrechitHFEMP_Emax = -1.,  pfrechitHFHADP_Emax = -1.,
         pfrechitHFEMN_Emax = -1.,  pfrechitHFHADN_Emax = -1. ;

  for(const auto& pfrechit : *(pfrechitsHF.product()) ){
    int ieta = HcalDetId(pfrechit.detId()).ieta() ;
    int iphi = HcalDetId(pfrechit.detId()).iphi() ;
    int idep = HcalDetId(pfrechit.detId()).depth() ;
    double E = pfrechit.energy() ;
    if (debug_ && scan)
    LogPrint("PFTrackHFAnalyzer") << boost::format(" pfrechit (ieta, iphi, depth, E): (%3d, %3d, %2d, %6.2f)")  % ieta % iphi % idep % E;

    if(idep == 1) { // HFEM
      pfrechitHFEM_n++ ;
      pfrechitHFEM_n_ieta[abs(ieta)]++;
      FillHist1D("pfrechitHFEM_E", E, 1. );
      if(ieta > 0){
        if(E > pfrechitHFEMP_Emax ) pfrechitHFEMP_Emax = E ; 
        FillHist2D("GProfileHFEMP_E", ieta-genP_ieta1 , idphi(iphi, genP_iphi1), E/genP_E);
        FillHist2D("TProfileHFEMP_E", ieta-trkP_ieta1 , idphi(iphi, trkP_iphi1), E/trkP_P);
      } else {
        if(E > pfrechitHFEMN_Emax ) pfrechitHFEMN_Emax = E ;
        FillHist2D("GProfileHFEMN_E", ieta-genN_ieta1,  idphi(iphi, genN_iphi1), E/genN_E);
        FillHist2D("TProfileHFEMN_E", ieta-trkN_ieta1,  idphi(iphi, trkN_iphi1), E/trkN_P);
      }
    } else { // HFHAD
      pfrechitHFHAD_n++ ;
      pfrechitHFHAD_n_ieta[abs(ieta)]++;
      FillHist1D("pfrechitHFHAD_E", E, 1. );
      if(ieta > 0){
        if(E > pfrechitHFHADP_Emax ) pfrechitHFHADP_Emax = E ;
        FillHist2D("GProfileHFHADP_E", ieta-genP_ieta2,  idphi(iphi, genP_iphi2), E/genP_E);
        FillHist2D("TProfileHFHADP_E", ieta-trkP_ieta2,  idphi(iphi, trkP_iphi2), E/trkP_P);
      } else {
        if(E > pfrechitHFHADN_Emax ) pfrechitHFHADN_Emax = E ;
        FillHist2D("GProfileHFHADN_E", ieta-genN_ieta2,  idphi(iphi, genN_iphi2), E/genN_E);
        FillHist2D("TProfileHFHADN_E", ieta-trkN_ieta2,  idphi(iphi, trkN_iphi2), E/trkN_P);
      }
    }
  } // end of loop over pfrechitsHF
  if (debug_ && scan){
    LogPrint("PFTrackHFAnalyzer") << "  pfrechitsHFEM:  " << pfrechitHFEM_n;
    LogPrint("PFTrackHFAnalyzer") << "  pfrechitsHFHAD: " << pfrechitHFHAD_n;
  }
  FillHist1D("pfrechitHFEM_n",  pfrechitHFEM_n,  1. );
  FillHist1D("pfrechitHFHAD_n", pfrechitHFHAD_n, 1. );
  FillHist1D("pfrechitHFEM_occupancy",  double(pfrechitHFEM_n)/864.,  1. );
  FillHist1D("pfrechitHFHAD_occupancy", double(pfrechitHFHAD_n)/864., 1. );
  for (unsigned i=28; i<=41 ; ++i){
    double n = i < 40 ? 2*36 : 2*18 ; 
    FillHist2D("pfrechitHFEM_occupancy_vs_absieta",   i, double(pfrechitHFEM_n_ieta[i])/n,  1. );
    FillHist2D("pfrechitHFHAD_occupancy_vs_absieta",  i, double(pfrechitHFHAD_n_ieta[i])/n,  1. );
  }

  FillHist1D("pfrechitHFEM_Emax",  pfrechitHFEMP_Emax,  1. );
  FillHist1D("pfrechitHFEM_Emax",  pfrechitHFEMN_Emax,  1. );
  FillHist1D("pfrechitHFHAD_Emax", pfrechitHFHADP_Emax,  1. );
  FillHist1D("pfrechitHFHAD_Emax", pfrechitHFHADN_Emax,  1. );
  FillHist1D("pfrechitHF_Emax",  pfrechitHFEMP_Emax > pfrechitHFHADP_Emax ? pfrechitHFEMP_Emax : pfrechitHFHADP_Emax ,  1. );
  FillHist1D("pfrechitHF_Emax",  pfrechitHFEMN_Emax > pfrechitHFHADN_Emax ? pfrechitHFEMN_Emax : pfrechitHFHADN_Emax ,  1. );

  if (debug_ && scan)
  LogPrint("PFTrackHFAnalyzer") << "\n =========== pfclustersHF: =========== " << pfclustersHF->size();

  int pfclusHFEM_n  = 0, pfclusHFHAD_n  = 0,
      pfclusHFEMP_n = 0, pfclusHFHADP_n = 0 ;
  double pfclusHFEMP_Emax  = -1., pfclusHFEMN_Emax  = -1.,  pfclusHFHADP_Emax  = -1., pfclusHFHADN_Emax  = -1. ;
  double pfclusHFEMP_ptmax = -1., pfclusHFEMN_ptmax = -1.,  pfclusHFHADP_ptmax = -1., pfclusHFHADN_ptmax = -1. ;

  int pfclusHFEMP_nmatch = 0, pfclusHFHADP_nmatch = 0,  pfclusHFEMP_nmatch5 = 0, pfclusHFHADP_nmatch5 = 0,  pfclusHFEMP_nmatch9 = 0, pfclusHFHADP_nmatch9 = 0,
      pfclusHFEMN_nmatch = 0, pfclusHFHADN_nmatch = 0,  pfclusHFEMN_nmatch5 = 0, pfclusHFHADN_nmatch5 = 0,  pfclusHFEMN_nmatch9 = 0, pfclusHFHADN_nmatch9 = 0;

  double pfclusHFP_pttot = 0., pfclusHFP_pttot5 = 0., pfclusHFP_pttot9 = 0., 
         pfclusHFN_pttot = 0., pfclusHFN_pttot5 = 0., pfclusHFN_pttot9 = 0.;

  for(const auto& pfclus : *(pfclustersHF.product()) ){
    double eta = pfclus.eta() ;
    double phi = pfclus.phi() ;
    double pt  = pfclus.pt() ;
    double E   = pfclus.energy() ;
    int layer  = pfclus.layer() ; 
    int idep   = fabs(pfclus.depth() - 1.) < 0.001 ? 1 : 2 ;   
    const std::vector<reco::PFRecHitFraction> &fracs = pfclus.recHitFractions();   
    const std::vector<std::pair<DetId, float>> &hfracs = pfclus.hitsAndFractions();
    unsigned nhits = fracs.size();

    if (debug_ && scan)
    LogPrint("PFTrackHFAnalyzer") << boost::format("pfclus (pt,E,eta,phi,layer,depth,nhits): (%6.2f, %6.2f, %6.2f, %6.2f, %i, %i, %i)")
      % pt % E % eta % phi % layer % idep % nhits;

    bool matchP = false,  match4P = false,  match8P = false,
         matchN = false,  match4N = false,  match8N = false;
    int nhits_threshold = 0;
    for(unsigned i=0; i<nhits; i++) { 
//    const auto& id = hfracs[i].first.rawId();
      const reco::PFRecHitRef& pfRecHits = fracs[i].recHitRef();
      double rawenergy = pfRecHits->energy();
      double frac      = fracs[i].fraction();
      int ieta  = HcalDetId(pfRecHits->detId()).ieta() ;
      int iphi  = HcalDetId(pfRecHits->detId()).iphi() ;
      if (frac > 0.01) nhits_threshold++ ;
    
      if (scan){ 
        if(idep == 1){
          hname = Form("ClustersEMP%i", nev);
          FillHist2D(hname, ieta,iphi,pfclusHFEMP_n+1) ; 
          hname = Form("EClustersEMP%i", nev);
          FillHist2D(hname, ieta,iphi,E) ;
          hname = Form("ERClustersEMP%i", nev);
          FillHist2D(hname, ieta,iphi,rawenergy*frac) ;
        } else {
          hname = Form("ClustersHADP%i", nev);
          FillHist2D(hname, ieta,iphi,pfclusHFHADP_n+1) ;
          hname = Form("EClustersHADP%i", nev);
          FillHist2D(hname, ieta,iphi,E) ;
          hname = Form("ERClustersHADP%i", nev);
          FillHist2D(hname, ieta,iphi,rawenergy*frac) ;
        } 
      } 

      if (debug_ && scan)
      LogPrint("PFTrackHFAnalyzer") << boost::format(" pfrechit (ieta, iphi, depth, E, frac): (%3d, %3d, %2d, %6.2f, %6.2f)")
      % HcalDetId(pfRecHits->detId()).ieta() % HcalDetId(pfRecHits->detId()).iphi()  % HcalDetId(pfRecHits->detId()).depth() % rawenergy % frac ;

      if(idep == 1){
        if( ieta == genP_ieta1 && iphi == genP_iphi1 ) matchP = true ;
        if( ieta == genN_ieta1 && iphi == genN_iphi1 ) matchN = true ;
        if( (abs(ieta - genP_ieta1) == 1 && iphi == genP_iphi1) || (abs(idphi(iphi, genP_iphi1)) == 1 && genP_ieta1 == ieta) ) match4P = true ;
        if( (abs(ieta - genN_ieta1) == 1 && iphi == genN_iphi1) || (abs(idphi(iphi, genN_iphi1)) == 1 && genN_ieta1 == ieta) ) match4N = true ;
        if( abs(ieta - genP_ieta1) == 1 && abs(idphi(iphi, genP_iphi1)) == 1 ) match8P = true ;
        if( abs(ieta - genN_ieta1) == 1 && abs(idphi(iphi, genN_iphi1)) == 1 ) match8N = true ;
      } else {
        if( ieta == genP_ieta2 && iphi == genP_iphi2 ) matchP = true ;
        if( ieta == genN_ieta2 && iphi == genN_iphi2 ) matchN = true ;
        if( (abs(ieta - genP_ieta2) == 1 && iphi == genP_iphi2) || (abs(idphi(iphi, genP_iphi2)) == 1 && genP_ieta2 == ieta) ) match4P = true ;
        if( (abs(ieta - genN_ieta2) == 1 && iphi == genN_iphi2) || (abs(idphi(iphi, genN_iphi2)) == 1 && genN_ieta2 == ieta) ) match4N = true ;
        if( abs(ieta - genP_ieta2) == 1 && abs(idphi(iphi, genP_iphi2)) == 1 ) match8P = true ;
        if( abs(ieta - genN_ieta2) == 1 && abs(idphi(iphi, genN_iphi2)) == 1 ) match8N = true ;
      }

      if (debug_ && scan && (matchP || match4P || match8P) ){
        LogPrint("PFTrackHFAnalyzer") << boost::format(" => Cluster matched to genP particle  with (eta phi E) (ieta1 iphi1) (ieta2 iphi2):  (%6.2f, %6.2f, %6.2f), (%i, %i), (%i, %i)")
        % genP_eta % genP_phi % genP_E % genP_ieta1 % genP_iphi1 % genP_ieta2 % genP_iphi2 ;
      }
      if (debug_ && scan && (matchN || match4N || match8N) ){
        LogPrint("PFTrackHFAnalyzer") << boost::format(" => Cluster matched to genN particle  with (eta phi E) (ieta1 iphi1) (ieta2 iphi2):  (%6.2f, %6.2f, %6.2f), (%i, %i), (%i, %i)")
        % genN_eta % genN_phi % genN_E % genN_ieta1 % genN_iphi1 % genN_ieta2 % genN_iphi2 ;
      }  
    } // end of loop over PFRecHits

    if(matchP) {
      if (idep == 1) pfclusHFEMP_nmatch++ ;
      else           pfclusHFHADP_nmatch++ ;
      pfclusHFP_pttot += pt;
    }
    if(matchN) {
      if (idep == 1) pfclusHFEMN_nmatch++ ;
      else           pfclusHFHADN_nmatch++ ;
      pfclusHFN_pttot += pt;
    }

    if(matchP || match4P ) {
      if (idep == 1) pfclusHFEMP_nmatch5++ ;
      else           pfclusHFHADP_nmatch5++ ;
      pfclusHFP_pttot5 += pt;
    }
      if(matchN || match4N ) {
      if (idep == 1) pfclusHFEMN_nmatch5++ ;
      else           pfclusHFHADN_nmatch5++ ;
      pfclusHFN_pttot5 += pt;
    }

    if(matchP || match4P || match8P ) {
      if (idep == 1) pfclusHFEMP_nmatch9++ ;
      else           pfclusHFHADP_nmatch9++ ;
      pfclusHFP_pttot9 += pt;
    } 
    if(matchN || match4N || match8N ) {
      if (idep == 1) pfclusHFEMN_nmatch9++ ;
      else           pfclusHFHADN_nmatch9++ ;
      pfclusHFN_pttot9 += pt;
    }

    if(idep == 1){
      pfclusHFEM_n++  ;
      FillHist1D("pfclusHFEM_nhits",  nhits, 1.);
      FillHist1D("pfclusHFEM_nhits_threshold",  nhits_threshold, 1.);
      FillHist1D("pfclusHFEM_E", E, 1.);
      FillHist1D("pfclusHFEM_pt",pt, 1.);
      if(eta > 0){
        pfclusHFEMP_n++  ;
        if(E > pfclusHFEMP_Emax ) pfclusHFEMP_Emax  = E ;
        if(pt> pfclusHFEMP_ptmax) pfclusHFEMP_ptmax = pt ;
      } else {
        if(E > pfclusHFEMN_Emax ) pfclusHFEMN_Emax  = E ;
        if(pt> pfclusHFEMN_ptmax) pfclusHFEMN_ptmax = pt ;
      }
    } else {
      pfclusHFHAD_n++ ; 
      FillHist1D("pfclusHFHAD_nhits", nhits, 1.);
      FillHist1D("pfclusHFHAD_E", E, 1.);
      FillHist1D("pfclusHFHAD_pt",pt, 1.);
      if(eta > 0){
        pfclusHFHADP_n++  ;
        if(E > pfclusHFHADP_Emax ) pfclusHFHADP_Emax  = E ;
        if(pt> pfclusHFHADP_ptmax) pfclusHFHADP_ptmax = pt ;
      } else {
        if(E > pfclusHFHADN_Emax ) pfclusHFHADN_Emax  = E ;
        if(pt> pfclusHFHADN_ptmax) pfclusHFHADN_ptmax = pt ;
      }
    } 
  } // end of loop over pfclusHF

  if (debug_ && scan){
    LogPrint("PFTrackHFAnalyzer") << "\n  pfclusHFEM_n:  " << pfclusHFEM_n;
    LogPrint("PFTrackHFAnalyzer") << "  pfclusHFHAD_n: " << pfclusHFHAD_n;

    LogPrint("PFTrackHFAnalyzer") << "\n  pfclusHFEMP_nmatch:  " << pfclusHFEMP_nmatch;
    LogPrint("PFTrackHFAnalyzer") << "  pfclusHFHADP_nmatch: " << pfclusHFHADP_nmatch;
    LogPrint("PFTrackHFAnalyzer") << "  pfclusHFP_pttot:     " << pfclusHFP_pttot;

    LogPrint("PFTrackHFAnalyzer") << "  pfclusHFEMN_nmatch:  " << pfclusHFEMN_nmatch;
    LogPrint("PFTrackHFAnalyzer") << "  pfclusHFHADN_nmatch: " << pfclusHFHADN_nmatch;
    LogPrint("PFTrackHFAnalyzer") << "  pfclusHFN_pttot:     " << pfclusHFN_pttot;

    LogPrint("PFTrackHFAnalyzer") << "\n  pfclusHFEMP_nmatch5:  " << pfclusHFEMP_nmatch5;
    LogPrint("PFTrackHFAnalyzer") << "  pfclusHFHADP_nmatch5: " << pfclusHFHADP_nmatch5;
    LogPrint("PFTrackHFAnalyzer") << "  pfclusHFP_pttot5:     " << pfclusHFP_pttot5;

    LogPrint("PFTrackHFAnalyzer") << "  pfclusHFEMN_nmatch5:  " << pfclusHFEMN_nmatch5;
    LogPrint("PFTrackHFAnalyzer") << "  pfclusHFHADN_nmatch5: " << pfclusHFHADN_nmatch5;
    LogPrint("PFTrackHFAnalyzer") << "  pfclusHFN_pttot5:     " << pfclusHFN_pttot5;

    LogPrint("PFTrackHFAnalyzer") << "\n  pfclusHFEMP_nmatch9:  " << pfclusHFEMP_nmatch9;
    LogPrint("PFTrackHFAnalyzer") << "  pfclusHFHADP_nmatch9: " << pfclusHFHADP_nmatch9;
    LogPrint("PFTrackHFAnalyzer") << "  pfclusHFP_pttot9:     " << pfclusHFP_pttot9;
    
    LogPrint("PFTrackHFAnalyzer") << "  pfclusHFEMN_nmatch9:  " << pfclusHFEMN_nmatch9;
    LogPrint("PFTrackHFAnalyzer") << "  pfclusHFHADN_nmatch9: " << pfclusHFHADN_nmatch9;
    LogPrint("PFTrackHFAnalyzer") << "  pfclusHFN_pttot9:     " << pfclusHFN_pttot9;

    LogPrint("PFTrackHFAnalyzer") << boost::format(" pfclusHFEMP_ptmax  pfclusHFHADP_ptmax: %6.2f, %6.2f ") % pfclusHFEMP_ptmax %  pfclusHFHADP_ptmax ;
    LogPrint("PFTrackHFAnalyzer") << boost::format(" pfclusHFEMN_ptmax  pfclusHFHADN_ptmax: %6.2f, %6.2f ") % pfclusHFEMN_ptmax %  pfclusHFHADN_ptmax ;
  }
  FillHist1D("pfclusHFEM_n",   pfclusHFEM_n, 1.);
  FillHist1D("pfclusHFHAD_n",  pfclusHFHAD_n, 1.);

  FillHist2D("pfclusHF_nmatch",  pfclusHFEMP_nmatch,  pfclusHFHADP_nmatch,  1.);
  FillHist2D("pfclusHF_nmatch",  pfclusHFEMN_nmatch,  pfclusHFHADN_nmatch,  1.);

  FillHist2D("pfclusHF_nmatch5", pfclusHFEMP_nmatch5, pfclusHFHADP_nmatch5, 1.);
  FillHist2D("pfclusHF_nmatch5", pfclusHFEMN_nmatch5, pfclusHFHADN_nmatch5, 1.);

  FillHist2D("pfclusHF_nmatch9", pfclusHFEMP_nmatch9, pfclusHFHADP_nmatch9, 1.);
  FillHist2D("pfclusHF_nmatch9", pfclusHFEMN_nmatch9, pfclusHFHADN_nmatch9, 1.);

  FillHist1D("match_ptfrac", pfclusHFP_pttot/genP_pt, 1.);
  FillHist1D("match_ptfrac", pfclusHFN_pttot/genN_pt, 1.);

  FillHist1D("match_ptfrac5", pfclusHFP_pttot5/genP_pt, 1.);
  FillHist1D("match_ptfrac5", pfclusHFN_pttot5/genN_pt, 1.);

  FillHist1D("match_ptfrac9", pfclusHFP_pttot9/genP_pt, 1.);
  FillHist1D("match_ptfrac9", pfclusHFN_pttot9/genN_pt, 1.);

  FillHist1D("pfclusHFEM_Emax",  pfclusHFEMP_Emax,  1. );
  FillHist1D("pfclusHFEM_Emax",  pfclusHFEMN_Emax,  1. );

  FillHist1D("pfclusHFHAD_Emax",  pfclusHFHADP_Emax,  1. );
  FillHist1D("pfclusHFHAD_Emax",  pfclusHFHADN_Emax,  1. );

  FillHist1D("pfclusHF_Emax",  pfclusHFEMP_Emax > pfclusHFHADP_Emax ? pfclusHFEMP_Emax : pfclusHFHADP_Emax ,  1. );
  FillHist1D("pfclusHF_Emax",  pfclusHFEMN_Emax > pfclusHFHADN_Emax ? pfclusHFEMN_Emax : pfclusHFHADN_Emax ,  1. );

  FillHist1D("pfclusHFEM_ptmax",  pfclusHFEMP_ptmax,  1. );
  FillHist1D("pfclusHFEM_ptmax",  pfclusHFEMN_ptmax,  1. );  
  
  FillHist1D("pfclusHFHAD_ptmax",  pfclusHFHADP_ptmax,  1. );
  FillHist1D("pfclusHFHAD_ptmax",  pfclusHFHADN_ptmax,  1. );
  
  FillHist1D("pfclusHF_ptmax",  pfclusHFEMP_ptmax > pfclusHFHADP_ptmax ? pfclusHFEMP_ptmax : pfclusHFHADP_ptmax ,  1. );
  FillHist1D("pfclusHF_ptmax",  pfclusHFEMN_ptmax > pfclusHFHADN_ptmax ? pfclusHFEMN_ptmax : pfclusHFHADN_ptmax ,  1. );


  if (debug_ && scan)
  LogPrint("PFTrackHFAnalyzer") << "\n =========== pfcands: =========== " << pfcands->size();
  int pfcandHFEM_n = 0, pfcandHFHAD_n = 0, pfcandHFCH_n = 0 ;
  vector<reco::PFCandidate> matchPFcandP, matchPFcandN ;
  int track_poseta_n = 0, track_negeta_n = 0 ;

  for(const auto& pfcand : *(pfcands.product()) ){
//  if (pfcand.trackRef().isNonnull())  cout << "All PFCands: pfcand.trackRef().get() "  << pfcand.trackRef().get()  << endl ;
    double eta = pfcand.eta() ;
    double phi = pfcand.phi() ;
    double pt  = pfcand.pt() ;
    double E   = pfcand.energy() ;
    reco::PFCandidate::ParticleType id = pfcand.particleId();
    const reco::TrackRef trackrefPF = pfcand.trackRef();
    //reco::TrackRef trackrefPF = pfcand.trackRef().get();

    const reco::PFCandidate::ElementsInBlocks& theElements = pfcand.elementsInBlocks();
    int nblocks = theElements.size() ;
    if( nblocks < 1 ){
      if(debug_)
      LogPrint("PFTrackHFAnalyzer") << 
      boost::format(" Warning!! PFcandidate with no blocks:  pfcand (pt,E,eta,phi,id,nblocks): (%6.2f, %6.2f, %6.2f, %6.2f, %i, %i)") %
      pt % E % eta % phi % id  % nblocks  ;
//      continue ;
    }

    const reco::PFBlockRef blockRef = nblocks > 0 ? theElements[0].first : edm::Ref<std::vector<reco::PFBlock> >() ;
    const edm::OwnVector<reco::PFBlockElement>& elements = nblocks > 0 ? blockRef->elements() : 0. ;
    unsigned int nele = elements.size() ;	
    unsigned int neleHFHAD = 0, neleHFEM = 0;

    bool isHF = false, hasBREM = false ;
    vector<reco::PFRecHitRef> PFrechitsHF;
    bool track_poseta = false, track_negeta = false ;

    for (unsigned int el = 0; el < nele; el++) {
     reco::PFBlockElement::Type type = elements[el].type();

      if(type == reco::PFBlockElement::HFEM || type == reco::PFBlockElement::HFHAD){
        if (type == reco::PFBlockElement::HFEM) neleHFEM++;
        else if (type == reco::PFBlockElement::HFHAD) neleHFHAD++;
        const std::vector<reco::PFRecHitFraction> &fracs = elements[el].clusterRef()->recHitFractions();
        isHF = true ; 
        unsigned nhits = fracs.size();
        for(unsigned i=0; i<nhits; i++) {
          const reco::PFRecHitRef& pfRecHits = fracs[i].recHitRef();
          if (std::find (PFrechitsHF.begin(), PFrechitsHF.end(), pfRecHits) == PFrechitsHF.end() )
          PFrechitsHF.push_back(pfRecHits);
        }
      }
      if(type == reco::PFBlockElement::BREM) hasBREM = true ;
    }  // end of loop over elements

    if (id == reco::PFCandidate::h ){
      if(trackrefPF == trackrefP ){
        track_poseta = true ; track_poseta_n++ ;}
      else if(trackrefPF == trackrefN ){
        track_negeta = true ; track_negeta_n++ ;}
      if(track_poseta||track_negeta){
        FillHist1D("pfcandCH_nelements", nele, 1. );
        FillHist1D("pfcandCH_nelementsHF", neleHFHAD+neleHFEM, 1. );
        FillHist1D("pfcandCH_nelementsHFHAD", neleHFHAD, 1. );
        FillHist1D("pfcandCH_nelementsHFEM", neleHFEM, 1. );
        FillHist1D("pfcandCH_nhits", PFrechitsHF.size(), 1. );
        if(hasBREM) FillHist1D("pfcandCH_brem_nelements",   nele, 1. );
        else        FillHist1D("pfcandCH_nobrem_nelements", nele, 1. );

        if(nele == 1){ // study why nele = 1 for "our" track
          if(track_poseta) FillHist1D("pfcandCH1ele_pftrackHF_n", pftrack_npos, 1. );
          if(track_negeta) FillHist1D("pfcandCH1ele_pftrackHF_n", pftrack_nneg, 1. );
          double drmin = 3.5 ;
          for(const auto& pfclus : *(pfclustersHF.product()) ){ // find closest cluster
            double cluseta = pfclus.eta() ;
            double clusphi = pfclus.phi() ;
            double trk_etaHF = track_poseta ? trkP_etaHF : trkN_etaHF ;
            double trk_phiHF = track_poseta ? trkP_phiHF : trkN_phiHF ;
            double dr = deltaR(trk_etaHF, trk_phiHF,  cluseta,  clusphi);
            if(dr < drmin) drmin = dr ;
          }
          FillHist1D("pfcandCH1ele_drmin", drmin, 1);
        }
      } else {
        FillHist1D("pfcandCHPU_nelements", nele, 1. );
        FillHist1D("pfcandCHPU_nelementsHF", neleHFEM+neleHFHAD, 1. );
        FillHist1D("pfcandCHPU_nelementsHFHAD", neleHFHAD, 1. );
        FillHist1D("pfcandCHPU_nelementsHFEM", neleHFEM, 1. );
        if(fabs(eta)>3.){
          FillHist1D("pfcandCHPUinHF_nelements", nele, 1. );
          FillHist1D("pfcandCHPUinHF_nhits", PFrechitsHF.size(), 1. );
          FillHist2D("pfcandCHPUinHF_pt_nelements", pt, nele, 1. );
        }
      }
    }

    if (!isHF){
      if (debug_ && scan) {
        LogPrint("PFTrackHFAnalyzer") << boost::format("[ pfcand (pt,E,eta,phi,id,nblocks, nelements): (%6.2f, %6.2f, %6.2f, %6.2f, %i, %i, %i)") %
        pt % E % eta % phi % id  % nblocks % nele ;
        for (unsigned int el = 0; el < nele; el++) {
          reco::PFBlockElement::Type type = elements[el].type();
          LogPrint("PFTrackHFAnalyzer") << boost::format("  element type :  %i ]") % type ;
        }
      }
    }


    // now only candidates with HF clusters ....
    if (!isHF) continue ;
    if (debug_ && scan)
    LogPrint("PFTrackHFAnalyzer") << boost::format("pfcand (pt,E,eta,phi,id,nblocks, nelements): (%6.2f, %6.2f, %6.2f, %6.2f, %i, %i, %i)") %
    pt % E % eta % phi % id  % nblocks % nele ;

    bool matchCandP = false, matchCandN = false ;
    double ehf = 0., enonhf = 0.  ;
    double tp = -1. ;

    for (unsigned int el = 0; el < nele; el++) {
      reco::PFBlockElement::Type type = elements[el].type();
      if (debug_ && scan)
      LogPrint("PFTrackHFAnalyzer") << boost::format("  element type :  %i ") % type ;

      if(type == reco::PFBlockElement::TRACK){
        if (debug_ && scan)
        LogPrint("PFTrackHFAnalyzer") << boost::format("   track in the block (pt, p, eta, phi, type): (%6.2f, %6.2f, %6.2f, %6.2f, %i)") %
        elements[el].trackRef()->pt() % elements[el].trackRef()->p() % elements[el].trackRef()->eta() % elements[el].trackRef()->phi() % type ;
        tp = elements[el].trackRef()->p() ;
      } else if(type == reco::PFBlockElement::HFEM || type == reco::PFBlockElement::HFHAD){
        if (debug_ && scan)
        LogPrint("PFTrackHFAnalyzer") << boost::format("   cluster in the block (pt, e, eta, phi, type): (%6.2f, %6.2f, %6.2f, %6.2f, %i)") %
        elements[el].clusterRef()->pt() % elements[el].clusterRef()->energy() % elements[el].clusterRef()->eta() % elements[el].clusterRef()->phi() % type ;
        const std::vector<reco::PFRecHitFraction> &fracs = elements[el].clusterRef()->recHitFractions();
        const std::vector<std::pair<DetId, float>> &hfracs = elements[el].clusterRef()->hitsAndFractions();
        ehf += elements[el].clusterRef()->energy() ;
        unsigned nhits = fracs.size();

        if(track_poseta || track_negeta){
          double trk_etaHF = track_poseta ? trkP_etaHF : trkN_etaHF ;
          double trk_phiHF = track_poseta ? trkP_phiHF : trkN_phiHF ;
          double dr = deltaR(trk_etaHF, trk_phiHF,  elements[el].clusterRef()->eta(),  elements[el].clusterRef()->phi());
          FillHist1D("HFCH_DR_track_HFclusters", dr, 1);
        }

        for(unsigned i=0; i<nhits; i++) {
          const reco::PFRecHitRef& pfRecHits = fracs[i].recHitRef();

          double rawenergy = pfRecHits->energy();
          double frac      = fracs[i].fraction();
          int ieta    = HcalDetId(pfRecHits->detId()).ieta() ;
          int iphi    = HcalDetId(pfRecHits->detId()).iphi() ;
          int idep    = HcalDetId(pfRecHits->detId()).depth() ; 

          if (debug_ && scan)
          LogPrint("PFTrackHFAnalyzer") << boost::format("      pfrechit (ieta, iphi, depth, E, frac): (%3d, %3d, %2d, %6.2f, %6.2f)")
          % ieta % iphi % idep % rawenergy % frac ; 

          if( idep == 1 && (abs(ieta - genP_ieta1)<= 1 && abs(idphi(iphi, genP_iphi1)) <=1 ) )   matchCandP = true ;
          if( idep == 2 && (abs(ieta - genP_ieta2)<= 1 && abs(idphi(iphi, genP_iphi2)) <=1 ) )   matchCandP = true ;
          if( idep == 1 && (abs(ieta - genN_ieta1)<= 1 && abs(idphi(iphi, genN_iphi1)) <=1 ) )   matchCandN = true ;
          if( idep == 2 && (abs(ieta - genN_ieta2)<= 1 && abs(idphi(iphi, genN_iphi2)) <=1 ) )   matchCandN = true ;

          if (scan && id == reco::PFCandidate::h ){
            if(idep == 1){
              hname = Form("hCandidatesEMP%i", nev);
              FillHist2D(hname, ieta,iphi,pfcandHFCH_n+1) ;
              hname = Form("ERhCandidatesEMP%i", nev);
              FillHist2D(hname, ieta,iphi,rawenergy*frac) ;
            } else {
              hname = Form("hCandidatesHADP%i", nev);
              FillHist2D(hname, ieta,iphi,pfclusHFHADP_n+1) ;
              hname = Form("ERhCandidatesHADP%i", nev);
              FillHist2D(hname, ieta,iphi,rawenergy*frac) ;
            }
          }


        }
      } else if(type == reco::PFBlockElement::BREM){
//      enonhf += elements[el].clusterRef()->energy() ;
//      const reco::PFBlockElementBrem& BREM = static_cast<const reco::PFBlockElementBrem&>(elemnts[el]);
        elements[el].Dump(std::cout) ;
        enonhf += 0. ;
      }
    } // end of loop over elements

    if (id == reco::PFCandidate::egamma_HF ){
      pfcandHFEM_n++ ;
      FillHist1D("blocksHFEM_n", nblocks, 1);
      FillHist1D("pfcandHFEM_pt", pt, 1. );
      FillHist1D("pfcandHFEM_E",  E,  1. );
      FillHist1D("pfcandHFEM_nelements", nele, 1. );
      FillHist1D("pfcandHFEM_nhits", PFrechitsHF.size(), 1. );
    } else if (id == reco::PFCandidate::h_HF ) {
      pfcandHFHAD_n++ ;
      FillHist1D("blocksHFHAD_n", nblocks, 1);
      FillHist1D("pfcandHFHAD_pt", pt, 1. );
      FillHist1D("pfcandHFHAD_E",  E,  1. );
      FillHist1D("pfcandHFHAD_nelements", nele, 1. );
      FillHist1D("pfcandHFHAD_nhits", PFrechitsHF.size(), 1. );
    } else if (id == reco::PFCandidate::h) {
      pfcandHFCH_n++ ;
      FillHist1D("blocksHFCH_n", nblocks, 1);
      FillHist1D("pfcandHFCH_pt", pt, 1. );
      FillHist1D("pfcandHFCH_E",  E,  1. );
    } else {
      cout << "Special Case: Candidate has HFCluster and is of the id " << id <<endl;
    }

    if(matchCandP)  matchPFcandP.push_back(pfcand);
    if(matchCandN)  matchPFcandN.push_back(pfcand);
    if( id == reco::PFCandidate::h ){
      FillHist1D("EovP",     (ehf+enonhf)/tp, 1.);
      if(hasBREM){
        if(track_poseta || track_negeta){
          FillHist1D("EovP_brem",     (ehf+enonhf)/tp, 1.);
        } else {
          FillHist1D("EovPPU_brem",   (ehf+enonhf)/tp, 1.);
        }
      } else {
        if(track_poseta || track_negeta){ 
          FillHist1D("EovP_nobrem",   ehf/tp, 1.);
        } else {
          FillHist1D("EovPPU_nobrem", ehf/tp, 1.);
        } 
      }

      if(eta > 0. && track_poseta) {
        if(hasBREM) FillHist1D("PovGenE_brem",  tp/genP_E, 1.);
        else        FillHist1D("PovGenE_nobrem",tp/genP_E, 1.); 
      } else if(eta < 0. && track_negeta) {
        if(hasBREM) FillHist1D("PovGenE_brem",  tp/genN_E, 1.);
        else        FillHist1D("PovGenE_nobrem",tp/genN_E, 1.);
      }
    }
  } // end of loop overy PFcandidates ...

  if (debug_ && scan){
    LogPrint("PFTrackHFAnalyzer") << boost::format(" track_poseta_n %i)") % track_poseta_n ;
    LogPrint("PFTrackHFAnalyzer") << boost::format(" track_negeta_n %i)") % track_negeta_n ;
  }

  FillHist2D("PFcand_vs_PFtrack", pftrackHF_n, pfcandHFCH_n, 1.);

  sort(matchPFcandP.begin(),matchPFcandP.end(), mysort);
  sort(matchPFcandN.begin(),matchPFcandN.end(), mysort);
  int nP = matchPFcandP.size() ;
  int nN = matchPFcandN.size() ;

  double etot = 0., etot_nontrk = 0.;
  int pfcandHFHAD_n9 = 0, pfcandHFEM_n9 = 0, pfcandHFCH_n9 = 0 ;
  for (int i=0; i<nP; i++){
    etot        += matchPFcandP[i].energy() ;
    reco::PFCandidate::ParticleType id = matchPFcandP[i].particleId();
    if (id == reco::PFCandidate::h_HF){
      pfcandHFHAD_n9++ ;
      etot_nontrk += matchPFcandP[i].energy() ;
    } else if (id == reco::PFCandidate::egamma_HF ) {
      pfcandHFEM_n9++ ;
      etot_nontrk += matchPFcandP[i].energy() ;
    } else if (id == reco::PFCandidate::h ) {
      pfcandHFCH_n9++ ;
    }
  }
  FillHist1D("pfcandHFEM_n9",  pfcandHFEM_n9, 1.);  
  FillHist1D("pfcandHFHAD_n9", pfcandHFHAD_n9, 1.);
  FillHist1D("pfcandHFCH_n9",  pfcandHFCH_n9, 1.);
  FillHist1D("Eratio_PFcandAlltoGen", etot/genP_E, 1.);
  FillHist1D("Eratio_PFcandNonTtoGen", etot_nontrk/genP_E, 1.);

  etot = 0., etot_nontrk = 0. ;
  pfcandHFHAD_n9 = 0; pfcandHFEM_n9 = 0, pfcandHFCH_n9 = 0 ;
  for (int i=0; i<nN; i++){
    etot += matchPFcandN[i].energy() ;
    reco::PFCandidate::ParticleType id = matchPFcandN[i].particleId();
    if (id == reco::PFCandidate::h_HF){
      pfcandHFHAD_n9++ ;
      etot_nontrk += matchPFcandN[i].energy() ;
    } else if (id == reco::PFCandidate::egamma_HF){
      pfcandHFEM_n9++ ;
      etot_nontrk += matchPFcandN[i].energy() ;
    } else if (id == reco::PFCandidate::h ) {
      pfcandHFCH_n9++ ;
    }
  }  
  FillHist1D("pfcandHFEM_n9",  pfcandHFEM_n9, 1.);
  FillHist1D("pfcandHFHAD_n9", pfcandHFHAD_n9, 1.);
  FillHist1D("pfcandHFCH_n9",  pfcandHFCH_n9, 1.);
  FillHist1D("Eratio_PFcandAlltoGen",  etot/genN_E, 1.);
  FillHist1D("Eratio_PFcandNonTtoGen", etot_nontrk/genN_E, 1.);

  if (nP > 0) {
    double r = matchPFcandP[0].energy()/genP_E ;
    FillHist1D("Eratio_PFcand1toGen", r, 1.);
    if (nP > 1) FillHist1D("Eratio_PFcand2toGen", (matchPFcandP[0].energy()+matchPFcandP[1].energy())/genP_E, 1. );
    if (nP > 2) FillHist1D("Eratio_PFcand3toGen", (matchPFcandP[0].energy()+matchPFcandP[1].energy()+matchPFcandP[2].energy())/genP_E, 1. );
    if(r< 0.05) cout << " Scan this event # " << nev <<  endl ;
    
  } 

  if (nN > 0) {
    double r = matchPFcandN[0].energy()/genN_E ;
    FillHist1D("Eratio_PFcand1toGen", r, 1.);
    if (nN > 1) FillHist1D("Eratio_PFcand2toGen", (matchPFcandN[0].energy()+matchPFcandN[1].energy())/genN_E, 1. );
    if (nN > 2) FillHist1D("Eratio_PFcand3toGen", (matchPFcandN[0].energy()+matchPFcandN[1].energy()+matchPFcandN[2].energy())/genN_E, 1. );
  }

  if (debug_ && scan){
    LogPrint("PFTrackHFAnalyzer") << "  pfcandsHF:    "  << pfcandHFEM_n+pfcandHFHAD_n+pfcandHFCH_n;
    LogPrint("PFTrackHFAnalyzer") << "  pfcandsHFEM:  "  << pfcandHFEM_n;
    LogPrint("PFTrackHFAnalyzer") << "  pfcandsHFHAD: "  << pfcandHFHAD_n;
    LogPrint("PFTrackHFAnalyzer") << "  pfcandsHFCH:  "  << pfcandHFCH_n;
  }

  FillHist1D("pfcandHFEM_n",  pfcandHFEM_n, 1.);
  FillHist1D("pfcandHFHAD_n", pfcandHFHAD_n, 1.);
  FillHist1D("pfcandHFCH_n",  pfcandHFCH_n, 1.);

 
  // second loop over pfcandidates: Offset part 
  float eFlavor[numFlavors][ETA_BINS] = {};
  for(const auto& pfcand : *(pfcands.product()) ){
    // remove single pions: keep only PU-originated pfcandidates ...
    if( pfcand.eta() > 0. ){
      double deltaphiP = fabs(pfcand.phi() - genP_phi) ;
      if( deltaphiP > M_PI ) deltaphiP = 2. * M_PI  - deltaphiP ;
      if( deltaphiP < M_PI/2. )  continue ; 
    }
    if( pfcand.eta() < 0. ){
      double deltaphiN = fabs(pfcand.phi() - genN_phi) ;
      if( deltaphiN > M_PI ) deltaphiN = 2. * M_PI  - deltaphiN ;
      if( deltaphiN < M_PI/2. )  continue ;
    }

    int etaIndex = getEtaIndex(pfcand.eta());
    Flavor flavor = getFlavor(pfcand.particleId());
    if (etaIndex == -1 || flavor == X) continue;

    bool attached1 = false;
    reco::TrackRef trackref = pfcand.trackRef(); 
    if (flavor == chm && !trackref.isNull() ) {
      for(const auto& pv : *(vertices.product()) ){
        if ( !pv.isFake() && pv.ndof() >= 4.0 && fabs(pv.z()) <= 24.0 && fabs(pv.position().rho())<=2.0 ) {
          reco::Vertex::trackRef_iterator i_vtxTrk, endvtxTrk = pv.tracks_end();
          for(i_vtxTrk = pv.tracks_begin(); i_vtxTrk != endvtxTrk && !attached1; ++i_vtxTrk) {
            reco::TrackRef vtxTrk(i_vtxTrk->castTo<reco::TrackRef>());
            if (vtxTrk == trackref) attached1 = true;
          }
        } 
      }
      if (!attached1) flavor = chu; //unmatched charged hadron
    }

    eFlavor[flavor][etaIndex] += pfcand.energy();
  }


//------------ Tracks ------------//

/*
  for(const auto& track : *(tracks.product()) ){

    if ( !track.quality(reco::Track::tight) ) continue;
    bool matched = false;
//    reco::TrackRef vtxTrk(i_trk->castTo<reco::TrackRef>());

    for(const auto& pfcand : *(pfcands.product()) ){ 
      reco::TrackRef trackref = pfcand.trackRef();

      if ( track == pfcand.trackRef().get() )   
      matched = true;      
    }
    if (matched) continue;

    int etaIndex = getEtaIndex( track.eta() );
    if (etaIndex == -1) continue;

    float e = track.p();
    eFlavor[untrk][etaIndex] += e;
  }
*/

  vector<reco::Track>::const_iterator i_trk, endtrk = tracks->end();
  for (i_trk = tracks->begin(); i_trk != endtrk; ++i_trk) {

    if ( !i_trk->quality(reco::Track::tight) ) continue;
    bool matched = false;

    vector<reco::PFCandidate>::const_iterator i_pf, endpf = pfcands->end();
    for (i_pf = pfcands->begin();  i_pf != endpf && !matched; ++i_pf) {
      if ( &(*i_trk) == i_pf->trackRef().get() ) matched = true;      
    }
    if (matched) continue;

    int etaIndex = getEtaIndex( i_trk->eta() );
    if (etaIndex == -1) continue;

    float e = i_trk->p();
    eFlavor[untrk][etaIndex] += e;
  }



  int intmu = mu / 10 ;
  for (int ieta = 0; ieta != ETA_BINS; ++ieta){
    double eta = 0.5*(etabins[ieta] + etabins[ieta+1]);
    for (int ifl = 0; ifl != numFlavors; ++ifl){
      double area = M_PI * (etabins[ieta+1] - etabins[ieta]);  // only half pi
      hname = Form("edensity_eta_nPU%i_", intmu) + ids[ifl];
      FillProfile(hname, eta, eFlavor[ifl][ieta]/area, 1.); 
      hname = Form("ptdensity_eta_nPU%i_", intmu) + ids[ifl];
      FillProfile(hname, eta, eFlavor[ifl][ieta]/area/cosh(eta), 1.); 
      // cout << " ieta, efl, e " << ieta << "  " << ifl << "  " << eFlavor[ifl][ieta] << endl ;
    }
  }

  LogPrint("PFTrackHFAnalyzer") << "\n\n"  ; 

}


// ------------ method called once each job just before starting event loop  ------------
void
PFTrackHFAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
PFTrackHFAnalyzer::endJob()
{
  if (debug_) {
    // LogPrint("PFTrackHFAnalyzer") << "\n";

  }

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PFTrackHFAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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

void PFTrackHFAnalyzer::FillHist1D(const TString& histName, const Double_t& value, const double& weight) 
{
  map<TString, TH1*>::iterator hid=m_Histos1D.find(histName);
  if (hid==m_Histos1D.end())
    cout << "%FillHist1D -- Could not find histogram with name: " << histName << endl;
  else
    hid->second->Fill(value, weight);
}

void PFTrackHFAnalyzer::FillHist2D(const TString& histName, const Double_t& value1, const Double_t& value2, const double& weight) 
{
  map<TString, TH2*>::iterator hid=m_Histos2D.find(histName);
  if (hid==m_Histos2D.end())
    cout << "%FillHist2D -- Could not find histogram with name: " << histName << endl;
  else
    hid->second->Fill(value1, value2, weight);
}

void PFTrackHFAnalyzer::FillProfile(const TString& histName, const Double_t& value1, const Double_t& value2, const double& weight) 
{
  map<TString, TProfile*>::iterator hid=m_Profiles.find(histName);
  if (hid==m_Profiles.end())
    cout << "%FillProfile -- Could not find profile with name: " << histName << endl;
  else
    hid->second->Fill(value1, value2, weight);
}

void PFTrackHFAnalyzer::FillProfile2D(const TString& histName, const Double_t& value1, const Double_t& value2, const Double_t& value3, const double& weight) 
{
  map<TString, TProfile2D*>::iterator hid=m_Profiles2D.find(histName);
  if (hid==m_Profiles2D.end())
    cout << "%FillProfile2D -- Could not find profile with name: " << histName << endl;
  else
    hid->second->Fill(value1, value2, value3, weight);
}


bool PFTrackHFAnalyzer::IsInCell(double eta, double phi, const CaloCellGeometry::CornersVec& CV)
{

  if (eta < CV[0].eta() && eta > CV[2].eta()) {  
    if (phi < CV[0].phi() && phi > CV[2].phi()){
    return true ; 
    } else if (CV[0].phi() < CV[2].phi()) {
    if ( phi < CV[0].phi()) return true ;
    if ( phi > CV[2].phi()) return true ;
    }
  }

  return false ;
}


int PFTrackHFAnalyzer::idphi(int iphiA, int iphiB)
{ 
  int d = (iphiA - iphiB)/2 ;
  if(d >  18) d = 18 - d ;
  if(d < -18) d = 36 + d ;
  return d ; 
}

int PFTrackHFAnalyzer::getEtaIndex(float eta){

  for (int i=0; i != ETA_BINS; ++i){
    if (etabins[i] <= eta && eta < etabins[i+1]) return i;
  }
  if (eta == etabins[ETA_BINS]) return ETA_BINS-1;
  else return -1;
}


PFTrackHFAnalyzer::Flavor PFTrackHFAnalyzer::getFlavor(reco::PFCandidate::ParticleType id)
{
    if (id == reco::PFCandidate::h)
        return chm;     //initially matched charged hadron
    else if (id == reco::PFCandidate::e)
        return lep;
    else if (id == reco::PFCandidate::mu)
        return lep;
    else if (id == reco::PFCandidate::gamma)
        return ne;
    else if (id == reco::PFCandidate::h0)
        return nh;
    else if (id == reco::PFCandidate::h_HF)
        return hfh;
    else if (id == reco::PFCandidate::egamma_HF)
        return hfe;
    else
        return X;
}

//define this as a plug-in
DEFINE_FWK_MODULE(PFTrackHFAnalyzer);
