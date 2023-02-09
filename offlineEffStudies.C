#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TEfficiency.h"
#include "TTree.h"
#include "TBranch.h"
#include "TChain.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <map>
#include "TLorentzVector.h"
#include <math.h>
#include "MuonTrackTree.h"

using namespace std;

double muonmass = 0.10565837;
bool debug = false;

enum Sig { 
  Prompt = 0,
  DiMuon,
  LowPt,
  DisplacedOld,
  DisplacedNew,
};

float deltaR(float, float, float, float); 

std::string getProbeFilter(int);
float getLeadingPtCut(int);
float getTrailingPtCut(int);

void printProgBar(int);

double pt_bins[35]  = {0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.8,0.9,1,1.2,1.4,1.6,1.8,2.0,2.5,3.0,3.5,4.0,5.0,6.0,7.0,8.0,10,12,15,20};
double dz_bins[11]  = {-15, -8, -6, -4, -2, 0, 2, 4, 6, 8, 15};
double eta_bins[16] = {-2.4, -2.1, -1.6, -1.2, -1.04, -0.9, -0.3, -0.2,  0.2, 0.3, 0.9, 1.04, 1.2, 1.6, 2.1, 2.4};
double iso_bins[12] = { 0  , 0.02, 0.04, 0.06, 0.08,  0.1, 0.12, 0.16, 0.2, 0.3, 0.6, 1   };
double offlineIsoCut = 0.15;


void offlineEffStudies(TString inputfilename="nutples.root"){

  TFile* outfile = TFile::Open("offlineEffs.root","RECREATE");
  std::cout << "output file: " << outfile -> GetName() << std::endl;

  //Create histograms  

  TH1F* pt                      = new TH1F("pt"            ,"pt"        ,   200,  0, 20 );
  TH1F* pt1                     = new TH1F("pt1"           ,"pt1"       ,   200,  0, 20 );
  TH1F* pt2                     = new TH1F("pt2"           ,"pt2"       ,   200,  0, 20 );
  TH1F* pt3                     = new TH1F("pt3"           ,"pt3"       ,   200,  0, 20 );

  TH1F* ptGen                      = new TH1F("ptGen"            ,"ptGen"        ,   200,  0, 20 );
  TH1F* ptGen1                     = new TH1F("pt1Gen"           ,"pt1Gen"       ,   200,  0, 20 );
  TH1F* ptGen2                     = new TH1F("pt2Gen"           ,"pt2Gen"       ,   200,  0, 20 );
  TH1F* ptGen3                     = new TH1F("pt3Gen"           ,"pt3Gen"       ,   200,  0, 20 );

  TH1F* ptGenMatched                      = new TH1F("ptGenMatched"            ,"ptGenMatched"        ,   200,  0, 20 );
  TH1F* ptGenMatched1                     = new TH1F("ptGenMatched1"           ,"pt1GenMatched"       ,   200,  0, 20 );
  TH1F* ptGenMatched2                     = new TH1F("ptGen2Matched"           ,"pt2GenMatched"       ,   200,  0, 20 );
  TH1F* ptGenMatched3                     = new TH1F("ptGenMatched3"           ,"pt3GenMatched"       ,   200,  0, 20 );


  TH1F* ptMatched               = new TH1F("ptMatched"      ,"ptMatched"        ,   200,  0, 20 );
  TH1F* pt1Matched              = new TH1F("pt1Matched"     ,"pt1Matched"       ,   200,  0, 20 );
  TH1F* pt2Matched              = new TH1F("pt2Matched"     ,"pt2Matched"       ,   200,  0, 20 );
  TH1F* pt3Matched              = new TH1F("pt3Matched"     ,"pt3Matched"       ,   200,  0, 20 );

  TH1F* eta                      = new TH1F("eta"            ,"eta"        ,   15,  eta_bins );
  TH1F* eta1                     = new TH1F("eta1"           ,"eta1"       ,   15,  eta_bins );
  TH1F* eta2                     = new TH1F("eta2"           ,"eta2"       ,   15,  eta_bins );
  TH1F* eta3                     = new TH1F("eta3"           ,"eta3"       ,   15,  eta_bins );

  TH1F* etaGen                      = new TH1F("etaGen"            ,"etaGen"        ,   15,  eta_bins );
  TH1F* etaGen1                     = new TH1F("etaGen1"           ,"etaGen1"       ,   15,  eta_bins );
  TH1F* etaGen2                     = new TH1F("etaGen2"           ,"etaGen2"       ,   15,  eta_bins );
  TH1F* etaGen3                     = new TH1F("etaGen3"           ,"etaGen3"       ,   15,  eta_bins );

  TH1F* etaGenMatched                      = new TH1F("etaGenMatched"            ,"etaGenMatched"        ,   15,  eta_bins );
  TH1F* etaGenMatched1                     = new TH1F("etaGenMatched1"           ,"etaGenMatched1"       ,   15,  eta_bins );
  TH1F* etaGenMatched2                     = new TH1F("etaGenMatched2"           ,"etaGenMatched2"       ,   15,  eta_bins );
  TH1F* etaGenMatched3                     = new TH1F("etaGenMatched3"           ,"etaGenMatched3"       ,   15,  eta_bins );
 
  TH1F* etaMatched               = new TH1F("etaMatched"      ,"etaMatched"        ,   15,  eta_bins );
  TH1F* eta1Matched              = new TH1F("eta1Matched"     ,"eta1Matched"       ,   15,  eta_bins );
  TH1F* eta2Matched              = new TH1F("eta2Matched"     ,"eta2Matched"       ,   15,  eta_bins );
  TH1F* eta3Matched              = new TH1F("eta3Matched"     ,"eta3Matched"       ,   15,  eta_bins );

  TH1F* nGenMuon                = new TH1F("nGenMuon"        ,"nGenMuon"         ,   501,  -0.5, 500.5 );
  TH1F* nRecoMuon               = new TH1F("nRecoMuon"       ,"nRecoMuon"        ,   501,  -0.5, 500.5 );
  TH1F* nRecoMuonMatched        = new TH1F("nRecoMuonMatched","nRecoMuonMatched" ,   501,  -0.5, 500.5 );
 
  TH1F* m3Mu                    = new TH1F("m3Mu","m3Mu" ,   110,  -0.5, 10.5 );
  TH1F* m3MuSoft                = new TH1F("m3MuSoft","m3MuSoft" ,   110,  -0.5, 10.5 );
  TH1F* m3MuLoose               = new TH1F("m3MuLoose","m3MuLoose" ,   110,  -0.5, 10.5 );
  TH1F* m3MuCombi               = new TH1F("m3MuCombi","m3MuCombi" ,   110,  -0.5, 10.5 );
 
  TEfficiency* muonPt           = new TEfficiency("muonPt"          ,"muonPt"           ,   200,  0, 20 ); 
  TEfficiency* muonEta          = new TEfficiency("muonEta"         ,"muonEta"          ,   15,  eta_bins );
  TEfficiency* muonPhi          = new TEfficiency("muonPhi"         ,"muonPhi"          ,   14, -3.2, 3.2);
  TEfficiency* muonEff          = new TEfficiency("muonEff"         ,"muonEff"          ,    1,   0., 1.0);

  TEfficiency* softMuonPt           = new TEfficiency("softMuonPt"          ,"softMuonPt"           ,   200,  0, 20 ); 
  TEfficiency* softMuonEta          = new TEfficiency("softMuonEta"         ,"softMuonEta"          ,   15,  eta_bins );
  TEfficiency* softMuonPhi          = new TEfficiency("softMuonPhi"         ,"softMuonPhi"          ,   14, -3.2, 3.2);
  TEfficiency* softMuonEff          = new TEfficiency("softMuonEff"         ,"softMuonEff"          ,    1,   0., 1.0);
   
  TEfficiency* looseMuonPt           = new TEfficiency("looseMuonPt"          ,"looseMuonPt"           ,   200,  0, 20 ); 
  TEfficiency* looseMuonEta          = new TEfficiency("looseMuonEta"         ,"looseMuonEta"          ,   15,  eta_bins );
  TEfficiency* looseMuonPhi          = new TEfficiency("looseMuonPhi"         ,"looseMuonPhi"          ,   14, -3.2, 3.2);
  TEfficiency* looseMuonEff          = new TEfficiency("looseMuonEff"         ,"looseMuonEff"          ,    1,   0., 1.0);
   

  TH1F* trackPt                      = new TH1F("trackPt"            ,"trackPt"        ,   200,  0, 20 );
  TH1F* trackPt1                     = new TH1F("trackPt1"           ,"trackPt1"       ,   200,  0, 20 );
  TH1F* trackPt2                     = new TH1F("trackPt2"           ,"trackPt2"       ,   200,  0, 20 );
  TH1F* trackPt3                     = new TH1F("trackPt3"           ,"trackPt3"       ,   200,  0, 20 );
 
  TH1F* trackPtMatched               = new TH1F("trackPtMatched"      ,"trackPtMatched"        ,   200,  0, 20 );
  TH1F* trackPt1Matched              = new TH1F("trackPt1Matched"     ,"trackPt1Matched"       ,   200,  0, 20 );
  TH1F* trackPt2Matched              = new TH1F("trackPt2Matched"     ,"trackPt2Matched"       ,   200,  0, 20 );
  TH1F* trackPt3Matched              = new TH1F("trackPt3Matched"     ,"trackPt3Matched"       ,   200,  0, 20 );

  TH1F* trackEta                      = new TH1F("trackEta"            ,"trackEta"        ,   15,  eta_bins );
  TH1F* trackEta1                     = new TH1F("trackEta1"           ,"trackEta1"       ,   15,  eta_bins );
  TH1F* trackEta2                     = new TH1F("trackEta2"           ,"trackEta2"       ,   15,  eta_bins );
  TH1F* trackEta3                     = new TH1F("trackEta3"           ,"trackEta3"       ,   15,  eta_bins );
 
  TH1F* trackEtaMatched               = new TH1F("trackEtaMatched"      ,"trackEtaMatched"        ,   15,  eta_bins );
  TH1F* trackEta1Matched              = new TH1F("trackEta1Matched"     ,"trackEta1Matched"       ,   15,  eta_bins );
  TH1F* trackEta2Matched              = new TH1F("trackEta2Matched"     ,"trackEta2Matched"       ,   15,  eta_bins );
  TH1F* trackEta3Matched              = new TH1F("trackEta3Matched"     ,"trackEta3Matched"       ,   15,  eta_bins );

  TH1F* nRecoTrack               = new TH1F("nRecoTrack"       ,"nRecoTrack"        ,   50001,  -0.5, 50000.5 );
  TH1F* nRecoTrackMatched        = new TH1F("nRecoTrackMatched","nRecoTrackMatched" ,   50001,  -0.5, 50000.5 );

 
  TH1F* m3Tracks                    = new TH1F("m3Track","m3Track" ,   110,  -0.5, 10.5 );
 
  TEfficiency* trackEffPt             = new TEfficiency("trackEffPt"          ,"trackEffPt"           ,   200,  0, 20 ); 
  TEfficiency* trackEffEta          = new TEfficiency("trackEffEta"         ,"trackEffEta"          ,   15,  eta_bins );
  TEfficiency* trackEffPhi          = new TEfficiency("trackEffPhi"         ,"trackEffPhi"          ,   14, -3.2, 3.2);
  TEfficiency* trackEff          = new TEfficiency("trackEff"         ,"trackEff"          ,    1,   0., 1.0);

  TH1F* m3TracksHP                    = new TH1F("m3TrackHP","m3TrackHP" ,   110,  -0.5, 10.5 );
 
  TEfficiency* trackEffHPPt             = new TEfficiency("trackEffHPPt"          ,"trackEffHPPt"           ,   200,  0, 20 ); 
  TEfficiency* trackEffHPEta          = new TEfficiency("trackEffHPEta"         ,"trackEffHPEta"          ,   15,  eta_bins );
  TEfficiency* trackEffHPPhi          = new TEfficiency("trackEffHPPhi"         ,"trackEffHPPhi"          ,   14, -3.2, 3.2);
  TEfficiency* trackEffHP          = new TEfficiency("trackEffHP"         ,"trackEffHP"          ,    1,   0., 1.0);
 



  //TFile* inputfile = TFile::Open(inputfilename, "READ");
  // std::cout << "input file: " << inputfile -> GetName() << std::endl;

  //TTree *tree = (TTree*) inputfile -> Get("Ntuplizer/MuonTrackTree"); 

  TChain *tree = new TChain("Ntuplizer/MuonTrackTree");
  tree->Add(inputfilename);
  
  if (!tree) {
    std::cout << " *** tree not found *** " << std::endl;
    return;
  }
    
  //TBranch*  evBranch = tree->GetBranch("event"); 
  //evBranch -> SetAddress(&ev);

   UInt_t          n_EMTF_mu;
   vector<float>   *EMTF_mu_pt;
   vector<float>   *EMTF_mu_pt_xml;
   vector<float>   *EMTF_mu_eta;
   vector<float>   *EMTF_mu_theta;
   vector<float>   *EMTF_mu_phi;
   vector<int>     *EMTF_mu_charge;
   vector<int>     *EMTF_mu_mode;
   vector<int>     *EMTF_mu_endcap;
   vector<int>     *EMTF_mu_sector;
   vector<int>     *EMTF_mu_bx;
   vector<int>     *EMTF_mu_hitref1;
   vector<int>     *EMTF_mu_hitref2;
   vector<int>     *EMTF_mu_hitref3;
   vector<int>     *EMTF_mu_hitref4;
   UInt_t          n_barrel_mu;
   vector<float>   *barrel_mu_pt;
   vector<float>   *barrel_mu_eta;
   vector<float>   *barrel_mu_phi;
   vector<int>     *barrel_mu_charge;
   vector<int>     *barrel_mu_qual;
   UInt_t          n_ovrlap_mu;
   vector<float>   *ovrlap_mu_pt;
   vector<float>   *ovrlap_mu_eta;
   vector<float>   *ovrlap_mu_phi;
   vector<int>     *ovrlap_mu_charge;
   vector<int>     *ovrlap_mu_qual;
   UInt_t          n_L1TT_trk;
   vector<float>   *L1TT_trk_pt;
   vector<float>   *L1TT_trk_eta;
   vector<float>   *L1TT_trk_phi;
   vector<int>     *L1TT_trk_charge;
   vector<float>   *L1TT_trk_p;
   vector<float>   *L1TT_trk_z;
   vector<float>   *L1TT_trk_chi2;
   vector<int>     *L1TT_trk_nstubs;
   vector<int>     *L1TT_trk_gen_qual;
   vector<int>     *L1TT_trk_gen_TP_ID;
   vector<float>   *L1TT_trk_gen_TP_pt;
   vector<float>   *L1TT_trk_gen_TP_eta;
   vector<float>   *L1TT_trk_gen_TP_phi;
   vector<float>   *L1TT_trk_gen_TP_m;
   UInt_t          n_L1_TkMu;
   vector<float>   *L1_TkMu_pt;
   vector<float>   *L1_TkMu_eta;
   vector<float>   *L1_TkMu_phi;
   vector<int>     *L1_TkMu_charge;
   vector<float>   *L1_TkMu_p;
   vector<float>   *L1_TkMu_z;
   vector<float>   *L1_TkMu_chi2;
   vector<int>     *L1_TkMu_nstubs;
   vector<int>     *L1_TkMu_mudetID;
   vector<int>     *L1_TkMu_gen_qual;
   vector<int>     *L1_TkMu_gen_TP_ID;
   vector<float>   *L1_TkMu_gen_TP_pt;
   vector<float>   *L1_TkMu_gen_TP_eta;
   vector<float>   *L1_TkMu_gen_TP_phi;
   vector<float>   *L1_TkMu_gen_TP_m;
   UInt_t          n_L1_TkMuStub;
   vector<float>   *L1_TkMuStub_pt;
   vector<float>   *L1_TkMuStub_eta;
   vector<float>   *L1_TkMuStub_phi;
   vector<int>     *L1_TkMuStub_charge;
   vector<float>   *L1_TkMuStub_p;
   vector<float>   *L1_TkMuStub_z;
   vector<float>   *L1_TkMuStub_chi2;
   vector<int>     *L1_TkMuStub_nstubs;
   vector<int>     *L1_TkMuStub_gen_qual;
   vector<int>     *L1_TkMuStub_gen_TP_ID;
   vector<float>   *L1_TkMuStub_gen_TP_pt;
   vector<float>   *L1_TkMuStub_gen_TP_eta;
   vector<float>   *L1_TkMuStub_gen_TP_phi;
   vector<float>   *L1_TkMuStub_gen_TP_m;
   UInt_t          n_gen_mu;
   vector<float>   *gen_mu_pt;
   vector<float>   *gen_mu_eta;
   vector<float>   *gen_mu_phi;
   vector<float>   *gen_mu_e;
   vector<int>     *gen_mu_charge;
   vector<int>     *gen_mu_gentauidx;
   UInt_t          n_mu_hit;
   vector<short>   *mu_hit_endcap;
   vector<short>   *mu_hit_station;
   vector<short>   *mu_hit_ring;
   vector<short>   *mu_hit_sector;
   vector<short>   *mu_hit_subsector;
   vector<short>   *mu_hit_chamber;
   vector<short>   *mu_hit_cscid;
   vector<short>   *mu_hit_bx;
   vector<short>   *mu_hit_type;
   vector<short>   *mu_hit_neighbor;
   vector<short>   *mu_hit_strip;
   vector<short>   *mu_hit_wire;
   vector<short>   *mu_hit_roll;
   vector<short>   *mu_hit_quality;
   vector<short>   *mu_hit_pattern;
   vector<short>   *mu_hit_bend;
   vector<short>   *mu_hit_time;
   vector<short>   *mu_hit_fr;
   vector<int>     *mu_hit_emtf_phi;
   vector<int>     *mu_hit_emtf_theta;
   vector<float>   *mu_hit_sim_phi;
   vector<float>   *mu_hit_sim_theta;
   vector<float>   *mu_hit_sim_eta;
   vector<float>   *mu_hit_sim_r;
   vector<float>   *mu_hit_sim_z;
   UInt_t          n_gen_tau;
   vector<float>   *gen_tau_pt;
   vector<float>   *gen_tau_eta;
   vector<float>   *gen_tau_phi;
   vector<float>   *gen_tau_e;
   vector<int>     *gen_tau_charge;
   UInt_t          n_reco_mu;
   vector<float>   *reco_mu_pt;
   vector<float>   *reco_mu_ptErr;
   vector<float>   *reco_mu_eta;
   vector<float>   *reco_mu_phi;
   vector<float>   *reco_mu_innerPt;
   vector<float>   *reco_mu_innerEta;
   vector<float>   *reco_mu_innerPhi;
   vector<float>   *reco_mu_e;
   vector<int>     *reco_mu_charge;
   vector<int>     *reco_mu_isTight;
   vector<int>     *reco_mu_isMedium;
   vector<int>     *reco_mu_isLoose;
   vector<int>     *reco_mu_isSoft;
   vector<float>   *reco_mu_trackIso;
   vector<float>   *reco_mu_chargedIso;
   vector<float>   *reco_mu_neutralIso;
   vector<float>   *reco_mu_photonIso;
   vector<float>   *reco_mu_PUIso;
   vector<float>   *reco_mu_dB;
   vector<float>   *reco_mu_edB;
   vector<int>     *reco_mu_validMuonHits;
   vector<float>   *reco_mu_normChi2;
   vector<int>     *reco_mu_matchedStations;
   vector<int>     *reco_mu_innerTrackNLayers;
   vector<int>     *reco_mu_innerTrackNPixHit;
   vector<int>     *reco_mu_isGlobal;
   vector<int>     *reco_mu_isPF;
   vector<int>     *reco_mu_isStandalone;
   vector<int>     *reco_mu_isTracker;
   UInt_t          n_reco_track;
   vector<float>   *reco_track_pt;
   vector<float>   *reco_track_ptErr;
   vector<float>   *reco_track_eta;
   vector<float>   *reco_track_phi;
   vector<int>     *reco_track_charge;
   vector<int>     *reco_track_isHighPurity;

   // List of branches
   TBranch        *b_n_EMTF_mu;   //!
   TBranch        *b_EMTF_mu_pt;   //!
   TBranch        *b_EMTF_mu_pt_xml;   //!
   TBranch        *b_EMTF_mu_eta;   //!
   TBranch        *b_EMTF_mu_theta;   //!
   TBranch        *b_EMTF_mu_phi;   //!
   TBranch        *b_EMTF_mu_charge;   //!
   TBranch        *b_EMTF_mu_mode;   //!
   TBranch        *b_EMTF_mu_endcap;   //!
   TBranch        *b_EMTF_mu_sector;   //!
   TBranch        *b_EMTF_mu_bx;   //!
   TBranch        *b_EMTF_mu_hitref1;   //!
   TBranch        *b_EMTF_mu_hitref2;   //!
   TBranch        *b_EMTF_mu_hitref3;   //!
   TBranch        *b_EMTF_mu_hitref4;   //!
   TBranch        *b_n_barrel_mu;   //!
   TBranch        *b_barrel_mu_pt;   //!
   TBranch        *b_barrel_mu_eta;   //!
   TBranch        *b_barrel_mu_phi;   //!
   TBranch        *b_barrel_mu_charge;   //!
   TBranch        *b_barrel_mu_qual;   //!
   TBranch        *b_n_ovrlap_mu;   //!
   TBranch        *b_ovrlap_mu_pt;   //!
   TBranch        *b_ovrlap_mu_eta;   //!
   TBranch        *b_ovrlap_mu_phi;   //!
   TBranch        *b_ovrlap_mu_charge;   //!
   TBranch        *b_ovrlap_mu_qual;   //!
   TBranch        *b_n_L1TT_trk;   //!
   TBranch        *b_L1TT_trk_pt;   //!
   TBranch        *b_L1TT_trk_eta;   //!
   TBranch        *b_L1TT_trk_phi;   //!
   TBranch        *b_L1TT_trk_charge;   //!
   TBranch        *b_L1TT_trk_p;   //!
   TBranch        *b_L1TT_trk_z;   //!
   TBranch        *b_L1TT_trk_chi2;   //!
   TBranch        *b_L1TT_trk_nstubs;   //!
   TBranch        *b_L1TT_trk_gen_qual;   //!
   TBranch        *b_L1TT_trk_gen_TP_ID;   //!
   TBranch        *b_L1TT_trk_gen_TP_pt;   //!
   TBranch        *b_L1TT_trk_gen_TP_eta;   //!
   TBranch        *b_L1TT_trk_gen_TP_phi;   //!
   TBranch        *b_L1TT_trk_gen_TP_m;   //!
   TBranch        *b_n_L1_TkMu;   //!
   TBranch        *b_L1_TkMu_pt;   //!
   TBranch        *b_L1_TkMu_eta;   //!
   TBranch        *b_L1_TkMu_phi;   //!
   TBranch        *b_L1_TkMu_charge;   //!
   TBranch        *b_L1_TkMu_p;   //!
   TBranch        *b_L1_TkMu_z;   //!
   TBranch        *b_L1_TkMu_chi2;   //!
   TBranch        *b_L1_TkMu_nstubs;   //!
   TBranch        *b_L1_TkMu_mudetID;   //!
   TBranch        *b_L1_TkMu_gen_qual;   //!
   TBranch        *b_L1_TkMu_gen_TP_ID;   //!
   TBranch        *b_L1_TkMu_gen_TP_pt;   //!
   TBranch        *b_L1_TkMu_gen_TP_eta;   //!
   TBranch        *b_L1_TkMu_gen_TP_phi;   //!
   TBranch        *b_L1_TkMu_gen_TP_m;   //!
   TBranch        *b_n_L1_TkMuStub;   //!
   TBranch        *b_L1_TkMuStub_pt;   //!
   TBranch        *b_L1_TkMuStub_eta;   //!
   TBranch        *b_L1_TkMuStub_phi;   //!
   TBranch        *b_L1_TkMuStub_charge;   //!
   TBranch        *b_L1_TkMuStub_p;   //!
   TBranch        *b_L1_TkMuStub_z;   //!
   TBranch        *b_L1_TkMuStub_chi2;   //!
   TBranch        *b_L1_TkMuStub_nstubs;   //!
   TBranch        *b_L1_TkMuStub_gen_qual;   //!
   TBranch        *b_L1_TkMuStub_gen_TP_ID;   //!
   TBranch        *b_L1_TkMuStub_gen_TP_pt;   //!
   TBranch        *b_L1_TkMuStub_gen_TP_eta;   //!
   TBranch        *b_L1_TkMuStub_gen_TP_phi;   //!
   TBranch        *b_L1_TkMuStub_gen_TP_m;   //!
   TBranch        *b_n_gen_mu;   //!
   TBranch        *b_gen_mu_pt;   //!
   TBranch        *b_gen_mu_eta;   //!
   TBranch        *b_gen_mu_phi;   //!
   TBranch        *b_gen_mu_e;   //!
   TBranch        *b_gen_mu_charge;   //!
   TBranch        *b_gen_mu_gentauidx;   //!
   TBranch        *b_n_mu_hit;   //!
   TBranch        *b_mu_hit_endcap;   //!
   TBranch        *b_mu_hit_station;   //!
   TBranch        *b_mu_hit_ring;   //!
   TBranch        *b_mu_hit_sector;   //!
   TBranch        *b_mu_hit_subsector;   //!
   TBranch        *b_mu_hit_chamber;   //!
   TBranch        *b_mu_hit_cscid;   //!
   TBranch        *b_mu_hit_bx;   //!
   TBranch        *b_mu_hit_type;   //!
   TBranch        *b_mu_hit_neighbor;   //!
   TBranch        *b_mu_hit_strip;   //!
   TBranch        *b_mu_hit_wire;   //!
   TBranch        *b_mu_hit_roll;   //!
   TBranch        *b_mu_hit_quality;   //!
   TBranch        *b_mu_hit_pattern;   //!
   TBranch        *b_mu_hit_bend;   //!
   TBranch        *b_mu_hit_time;   //!
   TBranch        *b_mu_hit_fr;   //!
   TBranch        *b_mu_hit_emtf_phi;   //!
   TBranch        *b_mu_hit_emtf_theta;   //!
   TBranch        *b_mu_hit_sim_phi;   //!
   TBranch        *b_mu_hit_sim_theta;   //!
   TBranch        *b_mu_hit_sim_eta;   //!
   TBranch        *b_mu_hit_sim_r;   //!
   TBranch        *b_mu_hit_sim_z;   //!
   TBranch        *b_n_gen_tau;   //!
   TBranch        *b_gen_tau_pt;   //!
   TBranch        *b_gen_tau_eta;   //!
   TBranch        *b_gen_tau_phi;   //!
   TBranch        *b_gen_tau_e;   //!
   TBranch        *b_gen_tau_charge;   //!
   TBranch        *b_n_reco_mu;   //!
   TBranch        *b_reco_mu_pt;   //!
   TBranch        *b_reco_mu_ptErr;   //!
   TBranch        *b_reco_mu_eta;   //!
   TBranch        *b_reco_mu_phi;   //!
   TBranch        *b_reco_mu_innerPt;   //!
   TBranch        *b_reco_mu_innerEta;   //!
   TBranch        *b_reco_mu_innerPhi;   //!
   TBranch        *b_reco_mu_e;   //!
   TBranch        *b_reco_mu_charge;   //!
   TBranch        *b_reco_mu_isTight;   //!
   TBranch        *b_reco_mu_isMedium;   //!
   TBranch        *b_reco_mu_isLoose;   //!
   TBranch        *b_reco_mu_isSoft;   //!
   TBranch        *b_reco_mu_trackIso;   //!
   TBranch        *b_reco_mu_chargedIso;   //!
   TBranch        *b_reco_mu_neutralIso;   //!
   TBranch        *b_reco_mu_photonIso;   //!
   TBranch        *b_reco_mu_PUIso;   //!
   TBranch        *b_reco_mu_dB;   //!
   TBranch        *b_reco_mu_edB;   //!
   TBranch        *b_reco_mu_validMuonHits;   //!
   TBranch        *b_reco_mu_normChi2;   //!
   TBranch        *b_reco_mu_matchedStations;   //!
   TBranch        *b_reco_mu_innerTrackNLayers;   //!
   TBranch        *b_reco_mu_innerTrackNPixHit;   //!
   TBranch        *b_reco_mu_isGlobal;   //!
   TBranch        *b_reco_mu_isPF;   //!
   TBranch        *b_reco_mu_isStandalone;   //!
   TBranch        *b_reco_mu_isTracker;   //!
   TBranch        *b_n_reco_track;   //!
   TBranch        *b_reco_track_pt;   //!
   TBranch        *b_reco_track_ptErr;   //!
   TBranch        *b_reco_track_eta;   //!
   TBranch        *b_reco_track_phi;   //!
   TBranch        *b_reco_track_charge;   //!
   TBranch        *b_reco_track_isHighPurity;   //!


   EMTF_mu_pt = 0;
   EMTF_mu_pt_xml = 0;
   EMTF_mu_eta = 0;
   EMTF_mu_theta = 0;
   EMTF_mu_phi = 0;
   EMTF_mu_charge = 0;
   EMTF_mu_mode = 0;
   EMTF_mu_endcap = 0;
   EMTF_mu_sector = 0;
   EMTF_mu_bx = 0;
   EMTF_mu_hitref1 = 0;
   EMTF_mu_hitref2 = 0;
   EMTF_mu_hitref3 = 0;
   EMTF_mu_hitref4 = 0;
   barrel_mu_pt = 0;
   barrel_mu_eta = 0;
   barrel_mu_phi = 0;
   barrel_mu_charge = 0;
   barrel_mu_qual = 0;
   ovrlap_mu_pt = 0;
   ovrlap_mu_eta = 0;
   ovrlap_mu_phi = 0;
   ovrlap_mu_charge = 0;
   ovrlap_mu_qual = 0;
   L1TT_trk_pt = 0;
   L1TT_trk_eta = 0;
   L1TT_trk_phi = 0;
   L1TT_trk_charge = 0;
   L1TT_trk_p = 0;
   L1TT_trk_z = 0;
   L1TT_trk_chi2 = 0;
   L1TT_trk_nstubs = 0;
   L1TT_trk_gen_qual = 0;
   L1TT_trk_gen_TP_ID = 0;
   L1TT_trk_gen_TP_pt = 0;
   L1TT_trk_gen_TP_eta = 0;
   L1TT_trk_gen_TP_phi = 0;
   L1TT_trk_gen_TP_m = 0;
   L1_TkMu_pt = 0;
   L1_TkMu_eta = 0;
   L1_TkMu_phi = 0;
   L1_TkMu_charge = 0;
   L1_TkMu_p = 0;
   L1_TkMu_z = 0;
   L1_TkMu_chi2 = 0;
   L1_TkMu_nstubs = 0;
   L1_TkMu_mudetID = 0;
   L1_TkMu_gen_qual = 0;
   L1_TkMu_gen_TP_ID = 0;
   L1_TkMu_gen_TP_pt = 0;
   L1_TkMu_gen_TP_eta = 0;
   L1_TkMu_gen_TP_phi = 0;
   L1_TkMu_gen_TP_m = 0;
   L1_TkMuStub_pt = 0;
   L1_TkMuStub_eta = 0;
   L1_TkMuStub_phi = 0;
   L1_TkMuStub_charge = 0;
   L1_TkMuStub_p = 0;
   L1_TkMuStub_z = 0;
   L1_TkMuStub_chi2 = 0;
   L1_TkMuStub_nstubs = 0;
   L1_TkMuStub_gen_qual = 0;
   L1_TkMuStub_gen_TP_ID = 0;
   L1_TkMuStub_gen_TP_pt = 0;
   L1_TkMuStub_gen_TP_eta = 0;
   L1_TkMuStub_gen_TP_phi = 0;
   L1_TkMuStub_gen_TP_m = 0;
   gen_mu_pt = 0;
   gen_mu_eta = 0;
   gen_mu_phi = 0;
   gen_mu_e = 0;
   gen_mu_charge = 0;
   gen_mu_gentauidx = 0;
   mu_hit_endcap = 0;
   mu_hit_station = 0;
   mu_hit_ring = 0;
   mu_hit_sector = 0;
   mu_hit_subsector = 0;
   mu_hit_chamber = 0;
   mu_hit_cscid = 0;
   mu_hit_bx = 0;
   mu_hit_type = 0;
   mu_hit_neighbor = 0;
   mu_hit_strip = 0;
   mu_hit_wire = 0;
   mu_hit_roll = 0;
   mu_hit_quality = 0;
   mu_hit_pattern = 0;
   mu_hit_bend = 0;
   mu_hit_time = 0;
   mu_hit_fr = 0;
   mu_hit_emtf_phi = 0;
   mu_hit_emtf_theta = 0;
   mu_hit_sim_phi = 0;
   mu_hit_sim_theta = 0;
   mu_hit_sim_eta = 0;
   mu_hit_sim_r = 0;
   mu_hit_sim_z = 0;
   gen_tau_pt = 0;
   gen_tau_eta = 0;
   gen_tau_phi = 0;
   gen_tau_e = 0;
   gen_tau_charge = 0;
   reco_mu_pt = 0;
   reco_mu_ptErr = 0;
   reco_mu_eta = 0;
   reco_mu_phi = 0;
   reco_mu_innerPt = 0;
   reco_mu_innerEta = 0;
   reco_mu_innerPhi = 0;
   reco_mu_e = 0;
   reco_mu_charge = 0;
   reco_mu_isTight = 0;
   reco_mu_isMedium = 0;
   reco_mu_isLoose = 0;
   reco_mu_isSoft = 0;
   reco_mu_trackIso = 0;
   reco_mu_chargedIso = 0;
   reco_mu_neutralIso = 0;
   reco_mu_photonIso = 0;
   reco_mu_PUIso = 0;
   reco_mu_dB = 0;
   reco_mu_edB = 0;
   reco_mu_validMuonHits = 0;
   reco_mu_normChi2 = 0;
   reco_mu_matchedStations = 0;
   reco_mu_innerTrackNLayers = 0;
   reco_mu_innerTrackNPixHit = 0;
   reco_mu_isGlobal = 0;
   reco_mu_isPF = 0;
   reco_mu_isStandalone = 0;
   reco_mu_isTracker = 0;
   reco_track_pt = 0;
   reco_track_ptErr = 0;
   reco_track_eta = 0;
   reco_track_phi = 0;
   reco_track_charge = 0;
   reco_track_isHighPurity = 0;


   tree->SetBranchAddress("n_EMTF_mu", &n_EMTF_mu, &b_n_EMTF_mu);
   tree->SetBranchAddress("EMTF_mu_pt", &EMTF_mu_pt, &b_EMTF_mu_pt);
   tree->SetBranchAddress("EMTF_mu_pt_xml", &EMTF_mu_pt_xml, &b_EMTF_mu_pt_xml);
   tree->SetBranchAddress("EMTF_mu_eta", &EMTF_mu_eta, &b_EMTF_mu_eta);
   tree->SetBranchAddress("EMTF_mu_theta", &EMTF_mu_theta, &b_EMTF_mu_theta);
   tree->SetBranchAddress("EMTF_mu_phi", &EMTF_mu_phi, &b_EMTF_mu_phi);
   tree->SetBranchAddress("EMTF_mu_charge", &EMTF_mu_charge, &b_EMTF_mu_charge);
   tree->SetBranchAddress("EMTF_mu_mode", &EMTF_mu_mode, &b_EMTF_mu_mode);
   tree->SetBranchAddress("EMTF_mu_endcap", &EMTF_mu_endcap, &b_EMTF_mu_endcap);
   tree->SetBranchAddress("EMTF_mu_sector", &EMTF_mu_sector, &b_EMTF_mu_sector);
   tree->SetBranchAddress("EMTF_mu_bx", &EMTF_mu_bx, &b_EMTF_mu_bx);
   tree->SetBranchAddress("EMTF_mu_hitref1", &EMTF_mu_hitref1, &b_EMTF_mu_hitref1);
   tree->SetBranchAddress("EMTF_mu_hitref2", &EMTF_mu_hitref2, &b_EMTF_mu_hitref2);
   tree->SetBranchAddress("EMTF_mu_hitref3", &EMTF_mu_hitref3, &b_EMTF_mu_hitref3);
   tree->SetBranchAddress("EMTF_mu_hitref4", &EMTF_mu_hitref4, &b_EMTF_mu_hitref4);
   tree->SetBranchAddress("n_barrel_mu", &n_barrel_mu, &b_n_barrel_mu);
   tree->SetBranchAddress("barrel_mu_pt", &barrel_mu_pt, &b_barrel_mu_pt);
   tree->SetBranchAddress("barrel_mu_eta", &barrel_mu_eta, &b_barrel_mu_eta);
   tree->SetBranchAddress("barrel_mu_phi", &barrel_mu_phi, &b_barrel_mu_phi);
   tree->SetBranchAddress("barrel_mu_charge", &barrel_mu_charge, &b_barrel_mu_charge);
   tree->SetBranchAddress("barrel_mu_qual", &barrel_mu_qual, &b_barrel_mu_qual);
   tree->SetBranchAddress("n_ovrlap_mu", &n_ovrlap_mu, &b_n_ovrlap_mu);
   tree->SetBranchAddress("ovrlap_mu_pt", &ovrlap_mu_pt, &b_ovrlap_mu_pt);
   tree->SetBranchAddress("ovrlap_mu_eta", &ovrlap_mu_eta, &b_ovrlap_mu_eta);
   tree->SetBranchAddress("ovrlap_mu_phi", &ovrlap_mu_phi, &b_ovrlap_mu_phi);
   tree->SetBranchAddress("ovrlap_mu_charge", &ovrlap_mu_charge, &b_ovrlap_mu_charge);
   tree->SetBranchAddress("ovrlap_mu_qual", &ovrlap_mu_qual, &b_ovrlap_mu_qual);
   tree->SetBranchAddress("n_L1TT_trk", &n_L1TT_trk, &b_n_L1TT_trk);
   tree->SetBranchAddress("L1TT_trk_pt", &L1TT_trk_pt, &b_L1TT_trk_pt);
   tree->SetBranchAddress("L1TT_trk_eta", &L1TT_trk_eta, &b_L1TT_trk_eta);
   tree->SetBranchAddress("L1TT_trk_phi", &L1TT_trk_phi, &b_L1TT_trk_phi);
   tree->SetBranchAddress("L1TT_trk_charge", &L1TT_trk_charge, &b_L1TT_trk_charge);
   tree->SetBranchAddress("L1TT_trk_p", &L1TT_trk_p, &b_L1TT_trk_p);
   tree->SetBranchAddress("L1TT_trk_z", &L1TT_trk_z, &b_L1TT_trk_z);
   tree->SetBranchAddress("L1TT_trk_chi2", &L1TT_trk_chi2, &b_L1TT_trk_chi2);
   tree->SetBranchAddress("L1TT_trk_nstubs", &L1TT_trk_nstubs, &b_L1TT_trk_nstubs);
   tree->SetBranchAddress("L1TT_trk_gen_qual", &L1TT_trk_gen_qual, &b_L1TT_trk_gen_qual);
   tree->SetBranchAddress("L1TT_trk_gen_TP_ID", &L1TT_trk_gen_TP_ID, &b_L1TT_trk_gen_TP_ID);
   tree->SetBranchAddress("L1TT_trk_gen_TP_pt", &L1TT_trk_gen_TP_pt, &b_L1TT_trk_gen_TP_pt);
   tree->SetBranchAddress("L1TT_trk_gen_TP_eta", &L1TT_trk_gen_TP_eta, &b_L1TT_trk_gen_TP_eta);
   tree->SetBranchAddress("L1TT_trk_gen_TP_phi", &L1TT_trk_gen_TP_phi, &b_L1TT_trk_gen_TP_phi);
   tree->SetBranchAddress("L1TT_trk_gen_TP_m", &L1TT_trk_gen_TP_m, &b_L1TT_trk_gen_TP_m);
   tree->SetBranchAddress("n_L1_TkMu", &n_L1_TkMu, &b_n_L1_TkMu);
   tree->SetBranchAddress("L1_TkMu_pt", &L1_TkMu_pt, &b_L1_TkMu_pt);
   tree->SetBranchAddress("L1_TkMu_eta", &L1_TkMu_eta, &b_L1_TkMu_eta);
   tree->SetBranchAddress("L1_TkMu_phi", &L1_TkMu_phi, &b_L1_TkMu_phi);
   tree->SetBranchAddress("L1_TkMu_charge", &L1_TkMu_charge, &b_L1_TkMu_charge);
   tree->SetBranchAddress("L1_TkMu_p", &L1_TkMu_p, &b_L1_TkMu_p);
   tree->SetBranchAddress("L1_TkMu_z", &L1_TkMu_z, &b_L1_TkMu_z);
   tree->SetBranchAddress("L1_TkMu_chi2", &L1_TkMu_chi2, &b_L1_TkMu_chi2);
   tree->SetBranchAddress("L1_TkMu_nstubs", &L1_TkMu_nstubs, &b_L1_TkMu_nstubs);
   tree->SetBranchAddress("L1_TkMu_mudetID", &L1_TkMu_mudetID, &b_L1_TkMu_mudetID);
   tree->SetBranchAddress("L1_TkMu_gen_qual", &L1_TkMu_gen_qual, &b_L1_TkMu_gen_qual);
   tree->SetBranchAddress("L1_TkMu_gen_TP_ID", &L1_TkMu_gen_TP_ID, &b_L1_TkMu_gen_TP_ID);
   tree->SetBranchAddress("L1_TkMu_gen_TP_pt", &L1_TkMu_gen_TP_pt, &b_L1_TkMu_gen_TP_pt);
   tree->SetBranchAddress("L1_TkMu_gen_TP_eta", &L1_TkMu_gen_TP_eta, &b_L1_TkMu_gen_TP_eta);
   tree->SetBranchAddress("L1_TkMu_gen_TP_phi", &L1_TkMu_gen_TP_phi, &b_L1_TkMu_gen_TP_phi);
   tree->SetBranchAddress("L1_TkMu_gen_TP_m", &L1_TkMu_gen_TP_m, &b_L1_TkMu_gen_TP_m);
   tree->SetBranchAddress("n_L1_TkMuStub", &n_L1_TkMuStub, &b_n_L1_TkMuStub);
   tree->SetBranchAddress("L1_TkMuStub_pt", &L1_TkMuStub_pt, &b_L1_TkMuStub_pt);
   tree->SetBranchAddress("L1_TkMuStub_eta", &L1_TkMuStub_eta, &b_L1_TkMuStub_eta);
   tree->SetBranchAddress("L1_TkMuStub_phi", &L1_TkMuStub_phi, &b_L1_TkMuStub_phi);
   tree->SetBranchAddress("L1_TkMuStub_charge", &L1_TkMuStub_charge, &b_L1_TkMuStub_charge);
   tree->SetBranchAddress("L1_TkMuStub_p", &L1_TkMuStub_p, &b_L1_TkMuStub_p);
   tree->SetBranchAddress("L1_TkMuStub_z", &L1_TkMuStub_z, &b_L1_TkMuStub_z);
   tree->SetBranchAddress("L1_TkMuStub_chi2", &L1_TkMuStub_chi2, &b_L1_TkMuStub_chi2);
   tree->SetBranchAddress("L1_TkMuStub_nstubs", &L1_TkMuStub_nstubs, &b_L1_TkMuStub_nstubs);
   tree->SetBranchAddress("L1_TkMuStub_gen_qual", &L1_TkMuStub_gen_qual, &b_L1_TkMuStub_gen_qual);
   tree->SetBranchAddress("L1_TkMuStub_gen_TP_ID", &L1_TkMuStub_gen_TP_ID, &b_L1_TkMuStub_gen_TP_ID);
   tree->SetBranchAddress("L1_TkMuStub_gen_TP_pt", &L1_TkMuStub_gen_TP_pt, &b_L1_TkMuStub_gen_TP_pt);
   tree->SetBranchAddress("L1_TkMuStub_gen_TP_eta", &L1_TkMuStub_gen_TP_eta, &b_L1_TkMuStub_gen_TP_eta);
   tree->SetBranchAddress("L1_TkMuStub_gen_TP_phi", &L1_TkMuStub_gen_TP_phi, &b_L1_TkMuStub_gen_TP_phi);
   tree->SetBranchAddress("L1_TkMuStub_gen_TP_m", &L1_TkMuStub_gen_TP_m, &b_L1_TkMuStub_gen_TP_m);
   tree->SetBranchAddress("n_gen_mu", &n_gen_mu, &b_n_gen_mu);
   tree->SetBranchAddress("gen_mu_pt", &gen_mu_pt, &b_gen_mu_pt);
   tree->SetBranchAddress("gen_mu_eta", &gen_mu_eta, &b_gen_mu_eta);
   tree->SetBranchAddress("gen_mu_phi", &gen_mu_phi, &b_gen_mu_phi);
   tree->SetBranchAddress("gen_mu_e", &gen_mu_e, &b_gen_mu_e);
   tree->SetBranchAddress("gen_mu_charge", &gen_mu_charge, &b_gen_mu_charge);
   tree->SetBranchAddress("gen_mu_gentauidx", &gen_mu_gentauidx, &b_gen_mu_gentauidx);
   tree->SetBranchAddress("n_mu_hit", &n_mu_hit, &b_n_mu_hit);
   tree->SetBranchAddress("mu_hit_endcap", &mu_hit_endcap, &b_mu_hit_endcap);
   tree->SetBranchAddress("mu_hit_station", &mu_hit_station, &b_mu_hit_station);
   tree->SetBranchAddress("mu_hit_ring", &mu_hit_ring, &b_mu_hit_ring);
   tree->SetBranchAddress("mu_hit_sector", &mu_hit_sector, &b_mu_hit_sector);
   tree->SetBranchAddress("mu_hit_subsector", &mu_hit_subsector, &b_mu_hit_subsector);
   tree->SetBranchAddress("mu_hit_chamber", &mu_hit_chamber, &b_mu_hit_chamber);
   tree->SetBranchAddress("mu_hit_cscid", &mu_hit_cscid, &b_mu_hit_cscid);
   tree->SetBranchAddress("mu_hit_bx", &mu_hit_bx, &b_mu_hit_bx);
   tree->SetBranchAddress("mu_hit_type", &mu_hit_type, &b_mu_hit_type);
   tree->SetBranchAddress("mu_hit_neighbor", &mu_hit_neighbor, &b_mu_hit_neighbor);
   tree->SetBranchAddress("mu_hit_strip", &mu_hit_strip, &b_mu_hit_strip);
   tree->SetBranchAddress("mu_hit_wire", &mu_hit_wire, &b_mu_hit_wire);
   tree->SetBranchAddress("mu_hit_roll", &mu_hit_roll, &b_mu_hit_roll);
   tree->SetBranchAddress("mu_hit_quality", &mu_hit_quality, &b_mu_hit_quality);
   tree->SetBranchAddress("mu_hit_pattern", &mu_hit_pattern, &b_mu_hit_pattern);
   tree->SetBranchAddress("mu_hit_bend", &mu_hit_bend, &b_mu_hit_bend);
   tree->SetBranchAddress("mu_hit_time", &mu_hit_time, &b_mu_hit_time);
   tree->SetBranchAddress("mu_hit_fr", &mu_hit_fr, &b_mu_hit_fr);
   tree->SetBranchAddress("mu_hit_emtf_phi", &mu_hit_emtf_phi, &b_mu_hit_emtf_phi);
   tree->SetBranchAddress("mu_hit_emtf_theta", &mu_hit_emtf_theta, &b_mu_hit_emtf_theta);
   tree->SetBranchAddress("mu_hit_sim_phi", &mu_hit_sim_phi, &b_mu_hit_sim_phi);
   tree->SetBranchAddress("mu_hit_sim_theta", &mu_hit_sim_theta, &b_mu_hit_sim_theta);
   tree->SetBranchAddress("mu_hit_sim_eta", &mu_hit_sim_eta, &b_mu_hit_sim_eta);
   tree->SetBranchAddress("mu_hit_sim_r", &mu_hit_sim_r, &b_mu_hit_sim_r);
   tree->SetBranchAddress("mu_hit_sim_z", &mu_hit_sim_z, &b_mu_hit_sim_z);
   tree->SetBranchAddress("n_gen_tau", &n_gen_tau, &b_n_gen_tau);
   tree->SetBranchAddress("gen_tau_pt", &gen_tau_pt, &b_gen_tau_pt);
   tree->SetBranchAddress("gen_tau_eta", &gen_tau_eta, &b_gen_tau_eta);
   tree->SetBranchAddress("gen_tau_phi", &gen_tau_phi, &b_gen_tau_phi);
   tree->SetBranchAddress("gen_tau_e", &gen_tau_e, &b_gen_tau_e);
   tree->SetBranchAddress("gen_tau_charge", &gen_tau_charge, &b_gen_tau_charge);
   tree->SetBranchAddress("n_reco_mu", &n_reco_mu, &b_n_reco_mu);
   tree->SetBranchAddress("reco_mu_pt", &reco_mu_pt, &b_reco_mu_pt);
   tree->SetBranchAddress("reco_mu_ptErr", &reco_mu_ptErr, &b_reco_mu_ptErr);
   tree->SetBranchAddress("reco_mu_eta", &reco_mu_eta, &b_reco_mu_eta);
   tree->SetBranchAddress("reco_mu_phi", &reco_mu_phi, &b_reco_mu_phi);
   tree->SetBranchAddress("reco_mu_innerPt", &reco_mu_innerPt, &b_reco_mu_innerPt);
   tree->SetBranchAddress("reco_mu_innerEta", &reco_mu_innerEta, &b_reco_mu_innerEta);
   tree->SetBranchAddress("reco_mu_innerPhi", &reco_mu_innerPhi, &b_reco_mu_innerPhi);
   tree->SetBranchAddress("reco_mu_e", &reco_mu_e, &b_reco_mu_e);
   tree->SetBranchAddress("reco_mu_charge", &reco_mu_charge, &b_reco_mu_charge);
   tree->SetBranchAddress("reco_mu_isTight", &reco_mu_isTight, &b_reco_mu_isTight);
   tree->SetBranchAddress("reco_mu_isMedium", &reco_mu_isMedium, &b_reco_mu_isMedium);
   tree->SetBranchAddress("reco_mu_isLoose", &reco_mu_isLoose, &b_reco_mu_isLoose);
   tree->SetBranchAddress("reco_mu_isSoft", &reco_mu_isSoft, &b_reco_mu_isSoft);
   tree->SetBranchAddress("reco_mu_trackIso", &reco_mu_trackIso, &b_reco_mu_trackIso);
   tree->SetBranchAddress("reco_mu_chargedIso", &reco_mu_chargedIso, &b_reco_mu_chargedIso);
   tree->SetBranchAddress("reco_mu_neutralIso", &reco_mu_neutralIso, &b_reco_mu_neutralIso);
   tree->SetBranchAddress("reco_mu_photonIso", &reco_mu_photonIso, &b_reco_mu_photonIso);
   tree->SetBranchAddress("reco_mu_PUIso", &reco_mu_PUIso, &b_reco_mu_PUIso);
   tree->SetBranchAddress("reco_mu_dB", &reco_mu_dB, &b_reco_mu_dB);
   tree->SetBranchAddress("reco_mu_edB", &reco_mu_edB, &b_reco_mu_edB);
   tree->SetBranchAddress("reco_mu_validMuonHits", &reco_mu_validMuonHits, &b_reco_mu_validMuonHits);
   tree->SetBranchAddress("reco_mu_normChi2", &reco_mu_normChi2, &b_reco_mu_normChi2);
   tree->SetBranchAddress("reco_mu_matchedStations", &reco_mu_matchedStations, &b_reco_mu_matchedStations);
   tree->SetBranchAddress("reco_mu_innerTrackNLayers", &reco_mu_innerTrackNLayers, &b_reco_mu_innerTrackNLayers);
   tree->SetBranchAddress("reco_mu_innerTrackNPixHit", &reco_mu_innerTrackNPixHit, &b_reco_mu_innerTrackNPixHit);
   tree->SetBranchAddress("reco_mu_isGlobal", &reco_mu_isGlobal, &b_reco_mu_isGlobal);
   tree->SetBranchAddress("reco_mu_isPF", &reco_mu_isPF, &b_reco_mu_isPF);
   tree->SetBranchAddress("reco_mu_isStandalone", &reco_mu_isStandalone, &b_reco_mu_isStandalone);
   tree->SetBranchAddress("reco_mu_isTracker", &reco_mu_isTracker, &b_reco_mu_isTracker);
   tree->SetBranchAddress("n_reco_track", &n_reco_track, &b_n_reco_track);
   tree->SetBranchAddress("reco_track_pt", &reco_track_pt, &b_reco_track_pt);
   tree->SetBranchAddress("reco_track_ptErr", &reco_track_ptErr, &b_reco_track_ptErr);
   tree->SetBranchAddress("reco_track_eta", &reco_track_eta, &b_reco_track_eta);
   tree->SetBranchAddress("reco_track_phi", &reco_track_phi, &b_reco_track_phi);
   tree->SetBranchAddress("reco_track_charge", &reco_track_charge, &b_reco_track_charge);
   tree->SetBranchAddress("reco_track_isHighPurity", &reco_track_isHighPurity, &b_reco_track_isHighPurity);


//   int n_gen_mu;
//   int n_reco_mu;
  
//   tree-> SetBranchAddress("n_gen_mu",&n_gen_mu);
//   tree-> SetBranchAddress("n_reco_mu",&n_reco_mu);


 
  int nentries = tree->GetEntries();
  std::cout << "Number of entries = " << nentries << std::endl;

  for (Int_t eventNo=0; eventNo < nentries; eventNo++)     {

    tree->GetEntry(eventNo);

    //Int_t IgetEvent   = tree   -> GetEvent(eventNo);
    printProgBar((int)(eventNo*100./nentries));


    if (!(n_gen_tau==1)) continue; 
    //if (!((gen_mu_e->at(0) > 2.5 && gen_mu_e->at(1) > 2.5) || (gen_mu_e->at(1) > 2.5 && gen_mu_e->at(2) > 2.5) || (gen_mu_e->at(0) > 2.5 && gen_mu_e->at(2) > 2.5))) continue;

    nGenMuon  -> Fill(n_gen_mu);
    nRecoMuon -> Fill(n_reco_mu);
    nRecoTrack -> Fill(n_reco_track);

    int nRecoMatched = 0; 
    int nTrackMatched = 0; 

    for (int i = 0; i < n_reco_mu; i++){

        pt -> Fill (reco_mu_pt->at(i));
        if (i == 0) pt1 -> Fill (reco_mu_pt->at(i));
        if (i == 1) pt2 -> Fill (reco_mu_pt->at(i));
        if (i == 2) pt3 -> Fill (reco_mu_pt->at(i));

        eta -> Fill (reco_mu_eta->at(i));
        if (i == 0) eta1 -> Fill (reco_mu_eta->at(i));
        if (i == 1) eta2 -> Fill (reco_mu_eta->at(i));
        if (i == 2) eta3 -> Fill (reco_mu_eta->at(i));

        bool matched = false;
        float dR = 0.05;

        // std::cout << "new muon" << std::endl;
        for (int y = 0; y < n_gen_mu; y++){
            //std::cout << deltaR(reco_mu_eta->at(i), reco_mu_phi->at(i), gen_mu_eta->at(y), gen_mu_phi->at(y)) << std::endl;
            if (deltaR(reco_mu_eta->at(i), reco_mu_phi->at(i), gen_mu_eta->at(y), gen_mu_phi->at(y)) < dR){
                dR = deltaR(reco_mu_eta->at(i), reco_mu_phi->at(i), gen_mu_eta->at(y), gen_mu_phi->at(y));
                matched = true;
                // std::cout << reco_mu_pt->at(i) << " " <<  gen_mu_pt->at(y) << " " << reco_mu_eta->at(i) << " " << gen_mu_eta->at(y) << " " << reco_mu_phi->at(i) << " " << gen_mu_phi->at(y) << std::endl;
            }
            
        }
        if (matched) nRecoMatched++;

        if (matched) ptMatched -> Fill(reco_mu_pt->at(i));
        if (matched && nRecoMatched == 1) pt1Matched -> Fill(reco_mu_pt->at(i));
        if (matched && nRecoMatched == 2) pt2Matched -> Fill(reco_mu_pt->at(i));
        if (matched && nRecoMatched == 3) pt3Matched -> Fill(reco_mu_pt->at(i));

        if (matched) etaMatched -> Fill(reco_mu_eta->at(i));
        if (matched && nRecoMatched == 1) eta1Matched -> Fill(reco_mu_eta->at(i));
        if (matched && nRecoMatched == 2) eta2Matched -> Fill(reco_mu_eta->at(i));
        if (matched && nRecoMatched == 3) eta3Matched -> Fill(reco_mu_eta->at(i));
    }

     for (int i = 0; i < n_reco_track; i++){

        trackPt -> Fill (reco_track_pt->at(i));
        if (i == 0) trackPt1 -> Fill (reco_track_pt->at(i));
        if (i == 1) trackPt2 -> Fill (reco_track_pt->at(i));
        if (i == 2) trackPt3 -> Fill (reco_track_pt->at(i));

        trackEta -> Fill (reco_track_eta->at(i));
        if (i == 0) trackEta1 -> Fill (reco_track_eta->at(i));
        if (i == 1) trackEta2 -> Fill (reco_track_eta->at(i));
        if (i == 2) trackEta3 -> Fill (reco_track_eta->at(i));

        bool matched = false;
        float dR = 0.1;
        // std::cout << "new track" << std::endl;
        for (int y = 0; y < n_gen_mu; y++){
            //std::cout << deltaR(reco_mu_eta->at(i), reco_mu_phi->at(i), gen_mu_eta->at(y), gen_mu_phi->at(y)) << std::endl;
            if (deltaR(reco_track_eta->at(i), reco_track_phi->at(i), gen_mu_eta->at(y), gen_mu_phi->at(y)) < dR){
                dR = deltaR(reco_track_eta->at(i), reco_track_phi->at(i), gen_mu_eta->at(y), gen_mu_phi->at(y));
                
                matched = true;                // std::cout << reco_track_pt->at(i) << " " <<  gen_mu_pt->at(y) << " " << reco_track_eta->at(i) << " " << gen_mu_eta->at(y) << " " << reco_track_phi->at(i) << " " << gen_mu_phi->at(y) << std::endl;
                
            }

        }

        if (matched) nTrackMatched++;
        if (matched) trackPtMatched -> Fill(reco_track_pt->at(i));
        if (matched && nTrackMatched == 1) trackPt1Matched -> Fill(reco_track_pt->at(i));
        if (matched && nTrackMatched == 2) trackPt2Matched -> Fill(reco_track_pt->at(i));
        if (matched && nTrackMatched == 3) trackPt3Matched -> Fill(reco_track_pt->at(i));

        if (matched) trackEtaMatched -> Fill(reco_track_eta->at(i));
        if (matched && nTrackMatched == 1) trackEta1Matched -> Fill(reco_track_eta->at(i));
        if (matched && nTrackMatched == 2) trackEta2Matched -> Fill(reco_track_eta->at(i));
        if (matched && nTrackMatched == 3) trackEta3Matched -> Fill(reco_track_eta->at(i));
    }

    nRecoTrackMatched -> Fill(nTrackMatched);

    std::vector<int> matched_indices;
    std::vector<int> matched_indices_soft;
    std::vector<int> matched_indices_loose;
    std::vector<int> matched_indices_track;
    std::vector<int> matched_indices_trackHP;

    std::vector<float> genMuPts = *gen_mu_pt;
    std::sort(genMuPts.begin(), genMuPts.end(), greater<float>());
    int indices[3];
 
    for (int i =0; i < genMuPts.size(); i++){
        for (int y = 0; y < gen_mu_pt->size();y++){
         if (gen_mu_pt->at(y) == genMuPts[i]) indices[i] = y;   
        }
    }

    ptGen1 -> Fill (gen_mu_pt->at(indices[0]));
    etaGen1 -> Fill (gen_mu_eta->at(indices[0]));
    ptGen2 -> Fill (gen_mu_pt->at(indices[1]));
    etaGen2 -> Fill (gen_mu_eta->at(indices[1]));
    ptGen3 -> Fill (gen_mu_pt->at(indices[2]));
    etaGen3 -> Fill (gen_mu_eta->at(indices[2]));
    
 
    std::vector<float> matchedPts;
    std::vector<float> matchedEtas;

    for (int i = 0; i < n_gen_mu; i++){

        bool matched = false;
        float dR = 0.1;
        int matchedIndex = -999;

        ptGen -> Fill (gen_mu_pt->at(i));
        etaGen -> Fill (gen_mu_eta->at(i));

        for (int y = 0; y < n_reco_mu; y++){
            
            if (deltaR(reco_mu_eta->at(y), reco_mu_phi->at(y), gen_mu_eta->at(i), gen_mu_phi->at(i)) < dR){
                dR = deltaR(reco_mu_eta->at(y), reco_mu_phi->at(y), gen_mu_eta->at(i), gen_mu_phi->at(i));
                matched = true;
                matchedIndex = y;
            }
        }
        if (matched && matchedIndex > 0) matched_indices.push_back(matchedIndex);

        if (matched) ptGenMatched -> Fill (gen_mu_pt->at(i));
        if (matched) etaGenMatched -> Fill (gen_mu_eta->at(i));
        if (matched) matchedPts.push_back(gen_mu_pt->at(i));   
        if (matched) matchedEtas.push_back(gen_mu_eta->at(i));   

        muonPt->Fill(matched, gen_mu_pt->at(i));
        muonEta->Fill(matched, gen_mu_eta->at(i));
        muonPhi->Fill(matched, gen_mu_phi->at(i));
        muonEff->Fill(matched, 0.5);

        matched = false;
        dR = 0.1;
        matchedIndex = -999;
        for (int y = 0; y < n_reco_mu; y++){
            
            if (deltaR(reco_mu_eta->at(y), reco_mu_phi->at(y), gen_mu_eta->at(i), gen_mu_phi->at(i)) < dR){
                dR = deltaR(reco_mu_eta->at(y), reco_mu_phi->at(y), gen_mu_eta->at(i), gen_mu_phi->at(i));
                if (reco_mu_isSoft->at(y)) matched = true;
                if (matched) matchedIndex = y;
            }
        }
        if (matched && matchedIndex > 0) matched_indices_soft.push_back(matchedIndex);
        softMuonPt->Fill(matched, gen_mu_pt->at(i));
        softMuonEta->Fill(matched, gen_mu_eta->at(i));
        softMuonPhi->Fill(matched, gen_mu_phi->at(i));
        softMuonEff->Fill(matched, 0.5);

        ptGenMatched -> Fill(gen_mu_pt->at(i));
        etaGenMatched -> Fill(gen_mu_eta->at(i));

        matched = false;
        dR = 0.1;
        matchedIndex = -999;        
        for (int y = 0; y < n_reco_mu; y++){
            
            if (deltaR(reco_mu_eta->at(y), reco_mu_phi->at(y), gen_mu_eta->at(i), gen_mu_phi->at(i)) < dR){
                dR = deltaR(reco_mu_eta->at(y), reco_mu_phi->at(y), gen_mu_eta->at(i), gen_mu_phi->at(i));
                if (reco_mu_isLoose->at(y)) matched = true;
                if (matched) matchedIndex = y;
            }
        }
        if (matched && matchedIndex > 0) matched_indices_loose.push_back(matchedIndex);
        looseMuonPt->Fill(matched, gen_mu_pt->at(i));
        looseMuonEta->Fill(matched, gen_mu_eta->at(i));
        looseMuonPhi->Fill(matched, gen_mu_phi->at(i));
        looseMuonEff->Fill(matched, 0.5);
        
        matched = false;
        bool matchedHP = false;
        dR = 0.1;
        matchedIndex = -999;
        float dRHP = 0.1;
        int matchedIndexHP = -999;
        for (int y = 0; y < n_reco_track; y++){

            if ((deltaR(reco_track_eta->at(y), reco_track_phi->at(y), gen_mu_eta->at(i), gen_mu_phi->at(i)) < dR) && (std::abs(reco_track_pt->at(y)/gen_mu_pt->at(i)-1) < 0.5)){
                dR = deltaR(reco_track_eta->at(y), reco_track_phi->at(y), gen_mu_eta->at(i), gen_mu_phi->at(i));
                matched = true;
                matchedIndex = y;
            }
            if (!reco_track_isHighPurity->at(y)) continue;
            if ((deltaR(reco_track_eta->at(y), reco_track_phi->at(y), gen_mu_eta->at(i), gen_mu_phi->at(i)) < dRHP) && (std::abs(reco_track_pt->at(y)/gen_mu_pt->at(i)-1) < 0.5)){
                dRHP = deltaR(reco_track_eta->at(y), reco_track_phi->at(y), gen_mu_eta->at(i), gen_mu_phi->at(i));
                matchedHP = true;
                matchedIndexHP = y;
            }

        }
        if (matched && matchedIndex > 0) matched_indices_track.push_back(matchedIndex);
        trackEffPt->Fill(matched, gen_mu_pt->at(i));
        trackEffEta->Fill(matched, gen_mu_eta->at(i));
        trackEffPhi->Fill(matched, gen_mu_phi->at(i));
        trackEff->Fill(matched, 0.5);

        if (matchedHP && matchedIndexHP > 0) matched_indices_trackHP.push_back(matchedIndexHP);
        trackEffHPPt->Fill(matchedHP, gen_mu_pt->at(i));
        trackEffHPEta->Fill(matchedHP, gen_mu_eta->at(i));
        trackEffHPPhi->Fill(matchedHP, gen_mu_phi->at(i));
        trackEffHP->Fill(matchedHP, 0.5);



    }
    

    std::sort(matchedPts.begin(), matchedPts.end(), greater<float>());
    std::vector<int> indicesGen;

    for (int i =0; i < matchedPts.size(); i++){
        for (int y = 0; y < gen_mu_pt->size();y++){
         if (gen_mu_pt->at(y) == matchedPts[i]) indicesGen.push_back(y);   
        }
    }
    if (indicesGen.size() > 0){
        ptGenMatched1 -> Fill (gen_mu_pt->at(indicesGen[0]));
        etaGenMatched1 -> Fill (gen_mu_eta->at(indicesGen[0]));
    }
    if (indicesGen.size() > 1){
        ptGenMatched2 -> Fill (gen_mu_pt->at(indicesGen[1]));
        etaGenMatched2 -> Fill (gen_mu_eta->at(indicesGen[1]));
    }
    if (indicesGen.size() > 2){
        ptGenMatched3 -> Fill (gen_mu_pt->at(indicesGen[2]));
        etaGenMatched3 -> Fill (gen_mu_eta->at(indicesGen[2]));
    }
    // std::cout << matched_indices.size()  << " " << matched_indices_soft.size() << " " << matched_indices_loose.size() << std::endl;

    nRecoMuonMatched -> Fill(matched_indices.size());


    if (matched_indices.size() == 3){

        TLorentzVector mu1;
        TLorentzVector mu2;
        TLorentzVector mu3;
        TLorentzVector tau;

        mu1.SetPtEtaPhiM(reco_mu_pt->at(matched_indices[0]), reco_mu_eta->at(matched_indices[0]), reco_mu_phi->at(matched_indices[0]), 0.105);
        mu2.SetPtEtaPhiM(reco_mu_pt->at(matched_indices[1]), reco_mu_eta->at(matched_indices[1]), reco_mu_phi->at(matched_indices[1]), 0.105);
        mu3.SetPtEtaPhiM(reco_mu_pt->at(matched_indices[2]), reco_mu_eta->at(matched_indices[2]), reco_mu_phi->at(matched_indices[2]), 0.105);

        tau = mu1 + mu2 + mu3;

        m3Mu -> Fill(tau.M()); 

    }    


    if (matched_indices_soft.size() == 3){

        TLorentzVector mu1;
        TLorentzVector mu2;
        TLorentzVector mu3;
        TLorentzVector tau;

        mu1.SetPtEtaPhiM(reco_mu_pt->at(matched_indices_soft[0]), reco_mu_eta->at(matched_indices_soft[0]), reco_mu_phi->at(matched_indices_soft[0]), 0.105);
        mu2.SetPtEtaPhiM(reco_mu_pt->at(matched_indices_soft[1]), reco_mu_eta->at(matched_indices_soft[1]), reco_mu_phi->at(matched_indices_soft[1]), 0.105);
        mu3.SetPtEtaPhiM(reco_mu_pt->at(matched_indices_soft[2]), reco_mu_eta->at(matched_indices_soft[2]), reco_mu_phi->at(matched_indices_soft[2]), 0.105);

        tau = mu1 + mu2 + mu3;

        m3MuSoft -> Fill(tau.M()); 

    }

    if (matched_indices_loose.size() == 3){

        TLorentzVector mu1;
        TLorentzVector mu2;
        TLorentzVector mu3;
        TLorentzVector tau;

        mu1.SetPtEtaPhiM(reco_mu_pt->at(matched_indices_loose[0]), reco_mu_eta->at(matched_indices_loose[0]), reco_mu_phi->at(matched_indices_loose[0]), 0.105);
        mu2.SetPtEtaPhiM(reco_mu_pt->at(matched_indices_loose[1]), reco_mu_eta->at(matched_indices_loose[1]), reco_mu_phi->at(matched_indices_loose[1]), 0.105);
        mu3.SetPtEtaPhiM(reco_mu_pt->at(matched_indices_loose[2]), reco_mu_eta->at(matched_indices_loose[2]), reco_mu_phi->at(matched_indices_loose[2]), 0.105);

        tau = mu1 + mu2 + mu3;

        m3MuLoose -> Fill(tau.M()); 

    }

    if (matched_indices_track.size() == 3){

        TLorentzVector track1;
        TLorentzVector track2;
        TLorentzVector track3;
        TLorentzVector tau;

        track1.SetPtEtaPhiM(reco_track_pt->at(matched_indices_track[0]), reco_track_eta->at(matched_indices_track[0]), reco_track_phi->at(matched_indices_track[0]), 0.105);
        track2.SetPtEtaPhiM(reco_track_pt->at(matched_indices_track[1]), reco_track_eta->at(matched_indices_track[1]), reco_track_phi->at(matched_indices_track[1]), 0.105);
        track3.SetPtEtaPhiM(reco_track_pt->at(matched_indices_track[2]), reco_track_eta->at(matched_indices_track[2]), reco_track_phi->at(matched_indices_track[2]), 0.105);

        tau = track1 + track2 + track3;

        m3Tracks -> Fill(tau.M()); 

    }

    if (matched_indices_trackHP.size() == 3){

        TLorentzVector track1;
        TLorentzVector track2;
        TLorentzVector track3;
        TLorentzVector tau;

        track1.SetPtEtaPhiM(reco_track_pt->at(matched_indices_trackHP[0]), reco_track_eta->at(matched_indices_trackHP[0]), reco_track_phi->at(matched_indices_trackHP[0]), 0.105);
        track2.SetPtEtaPhiM(reco_track_pt->at(matched_indices_trackHP[1]), reco_track_eta->at(matched_indices_trackHP[1]), reco_track_phi->at(matched_indices_trackHP[1]), 0.105);
        track3.SetPtEtaPhiM(reco_track_pt->at(matched_indices_trackHP[2]), reco_track_eta->at(matched_indices_trackHP[2]), reco_track_phi->at(matched_indices_trackHP[2]), 0.105);

        tau = track1 + track2 + track3;

        m3TracksHP -> Fill(tau.M()); 

    }


  }  
  
  //Writing the histograms in a file.
  outfile           -> cd();
  
  muonPt            -> Write();
  muonEta           -> Write();
  muonPhi           -> Write();
  muonEff           -> Write();
  softMuonPt            -> Write();
  softMuonEta           -> Write();
  softMuonPhi           -> Write();
  softMuonEff           -> Write();
  looseMuonPt            -> Write();
  looseMuonEta           -> Write();
  looseMuonPhi           -> Write();
  looseMuonEff           -> Write();
 
  pt -> Write();
  pt1 -> Write();
  pt2 -> Write();
  pt3 -> Write();

  eta -> Write();
  eta1 -> Write();
  eta2 -> Write();
  eta3 -> Write();

  ptGen -> Write();
  ptGen1 -> Write();
  ptGen2 -> Write();
  ptGen3 -> Write();

  etaGen -> Write();
  etaGen1 -> Write();
  etaGen2 -> Write();
  etaGen3 -> Write();

  ptGenMatched -> Write();
  ptGenMatched1 -> Write();
  ptGenMatched2 -> Write();
  ptGenMatched3 -> Write();

  etaGenMatched -> Write();
  etaGenMatched1 -> Write();
  etaGenMatched2 -> Write();
  etaGenMatched3 -> Write();

  ptMatched -> Write();
  pt1Matched -> Write();
  pt2Matched -> Write();
  pt3Matched -> Write();

  etaMatched -> Write();
  eta1Matched -> Write();
  eta2Matched -> Write();
  eta3Matched -> Write();

  nGenMuon -> Write();
  nRecoMuon -> Write();
  nRecoMuonMatched -> Write();

  m3Mu -> Write();
  m3MuSoft -> Write();
  m3MuLoose -> Write();
 

  
  trackEffPt            -> Write();
  trackEffEta           -> Write();
  trackEffPhi           -> Write();
  trackEff           -> Write();
 
  trackEffHPPt            -> Write();
  trackEffHPEta           -> Write();
  trackEffHPPhi           -> Write();
  trackEffHP           -> Write();
 
  trackPt -> Write();
  trackPt1 -> Write();
  trackPt2 -> Write();
  trackPt3 -> Write();

  trackEta -> Write();
  trackEta1 -> Write();
  trackEta2 -> Write();
  trackEta3 -> Write();

  trackPtMatched -> Write();
  trackPt1Matched -> Write();
  trackPt2Matched -> Write();
  trackPt3Matched -> Write();

  trackEtaMatched -> Write();
  trackEta1Matched -> Write();
  trackEta2Matched -> Write();
  trackEta3Matched -> Write();

  nRecoTrack -> Write();
  nRecoTrackMatched -> Write();

  m3Tracks -> Write();
  m3TracksHP -> Write();

  outfile          -> Close();  
  
  return;
}



float deltaR(float eta1, float phi1, float eta2, float phi2) {
    float deta = eta1 - eta2;
    float dphi = std::abs(phi1 - phi2);
    if (dphi > M_PI)
      dphi -= 2 * M_PI;
    return sqrt(deta * deta + dphi * dphi);
}



void printProgBar( int percent ){
  std::string bar;  
  for(int i = 0; i < 50; i++){
    if( i < (percent/2)){
      bar.replace(i,1,"=");
    }else if( i == (percent/2)){
      bar.replace(i,1,">");
    }else{
      bar.replace(i,1," ");
    }
  }

  std::cout<< "\r" "[" << bar << "] ";
  std::cout.width( 3 );
  std::cout<< percent << "%     " << std::flush;
}
