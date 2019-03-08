#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <string>

#include "RawIcrrStationEvent.h"
#include "RawAtriStationEvent.h"
#include "UsefulAraStationEvent.h"
#include "UsefulIcrrStationEvent.h"
#include "UsefulAtriStationEvent.h"
#include "AraGeomTool.h"

#include "TChain.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH3F.h"
#include "TH3D.h"
#include "TGraph2D.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TROOT.h"
#include "TPad.h"
#include "TPolyMarker3D.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLegend.h"

#include "calibrationTools.h"
#include "calibrationToolsVs3.h"
#include "evProcessTools.h"
#include "recoTools.h"
#include "analysisTools.h"

#include "Detector.h"
#include "Trigger.h"
#include "Settings.h"
#include "Report.h"
#include "Event.h"

using namespace std;

float statisticalRadiusReco(int arraySize, const float * const radii, const float * const value, float trueRadius,
                          TH1F * const statRReco_weight, TH1F * const statRReco_max, TH1F * const statRReco_eventStack);

float getMean(const vector<float>& thetas);
float getRMS(const vector<float>& thetas);

int main(int argc, char **argv){

gROOT->ProcessLine("#include <vector>");

//TChain *onionTree=new TChain("onionTree");
TChain *recoSettingsTree=new TChain("recoSettingsTree");
TChain *dataTree=new TChain("dataTree");
TChain *runInfoTree=new TChain("runInfoTree");

for(int i=1; i<argc; i++){

   TFile fp( argv[i] );

   if( fp.IsZombie() ){ cerr<<"File "<<argv[i]<<" is zombie. Skipping..."<<endl; continue; }
   if( fp.TestBit(TFile::kRecovered) ){ cerr<<"File "<<argv[i]<<" is recovered file. Skipping..."<<endl; continue; }

   recoSettingsTree->Add( argv[i] );
   dataTree->Add( argv[i] );
   runInfoTree->Add( argv[i] );

}


//TFile fp(argv[1],"read");

//TTree *onionTree=(TTree*)fp.Get("onionTree");
recoSettings *settings = new recoSettings();;
int nSideExp, nLayer;
//onionTree->SetBranchAddress("nSideExp",&nSideExp);
//onionTree->SetBranchAddress("nLayer",  &nLayer);

recoSettingsTree->SetBranchAddress("settings", &settings);

recoSettingsTree->GetEntry(0);

nSideExp = settings->nSideExp;
nLayer = settings->nLayer;

//onionTree->GetEntry(0);
//nSideExp = 7;
//nLayer= 25;
Healpix_Onion onion(nSideExp, nLayer, settings->layerFirstRadius, settings->layerLastRadius);
int nDir = onion.nDir;
cout<<"nSideExp: "<<onion.nSideExp<<" nLayer: "<<onion.nLayer<<endl;

//TTree *dataTree=(TTree*)fp.Get("dataTree");
recoData *dummyData = new recoData();
//vector<int>   * topMaxPixIdx             = &dummyData->topMaxPixIdx;
//vector<float> * topMaxPixCoherence       = &dummyData->topMaxPixCoherence;
//vector<int>   * maxPixIdxEachLayer       = &dummyData->maxPixIdxEachLayer;
//vector<float> * maxPixCoherenceEachLayer = &dummyData->maxPixCoherenceEachLayer;
dataTree->SetBranchAddress("summary", &dummyData);

int runEventCount, trigEventCount, recoEventCount;
int runStartTime, runEndTime;

runInfoTree->SetBranchAddress("runEventCount", &runEventCount);
runInfoTree->SetBranchAddress("trigEventCount", &trigEventCount);
runInfoTree->SetBranchAddress("recoEventCount", &recoEventCount);
runInfoTree->SetBranchAddress("utime_runStart", &runStartTime);
runInfoTree->SetBranchAddress("utime_runEnd", &runEndTime);

int totalRunEventCount, totalTrigEventCount, totalRecoEventCount;
int totalLiveTime = 0;
totalRunEventCount = totalTrigEventCount = totalRecoEventCount = 0;

for(int run=0; run<runInfoTree->GetEntries(); run++){

   runInfoTree->GetEntry(run);
   totalRunEventCount += runEventCount;
   totalTrigEventCount += trigEventCount;
   totalRecoEventCount += recoEventCount;

   if(runEndTime <= runStartTime){ cerr<<"Run "<<run<<" livetime error. Skipping...\n"; continue; }
   else totalLiveTime += (runEndTime - runStartTime);


}

printf("totalRunEventCount: %d\ttotalTrigEventCount: %d\ttotalRecoEventCount: %d\ntotalLiveTime :%ds\n", totalRunEventCount, totalTrigEventCount, totalRecoEventCount, totalLiveTime);

/*
dataTree->SetBranchAddress("weight",     &dummyData->weight);
dataTree->SetBranchAddress("trueZen",    &dummyData->trueZen);
dataTree->SetBranchAddress("trueAzi",    &dummyData->trueAzi);
dataTree->SetBranchAddress("recoZen",    &dummyData->recoZen);
dataTree->SetBranchAddress("recoAzi",    &dummyData->recoAzi);
dataTree->SetBranchAddress("trueRadius", &dummyData->trueRadius);
dataTree->SetBranchAddress("recoRadius", &dummyData->recoRadius);
dataTree->SetBranchAddress("recoChan",   &dummyData->recoChan);
dataTree->SetBranchAddress("maxPixIdx",  &dummyData->maxPixIdx);
dataTree->SetBranchAddress("coherence",  &dummyData->maxPixCoherence);
dataTree->SetBranchAddress("topN",       &dummyData->topN);

dataTree->SetBranchAddress("topMaxPixIdx",             &topMaxPixIdx);
dataTree->SetBranchAddress("topMaxPixCoherence",       &topMaxPixCoherence);
dataTree->SetBranchAddress("maxPixIdxEachLayer",       &maxPixIdxEachLayer);
dataTree->SetBranchAddress("maxPixCoherenceEachLayer", &maxPixCoherenceEachLayer);
*/
int Nentries = dataTree->GetEntries();
cout<<"Nentries: "<<Nentries<<endl;

TH1F *trueR_all=new TH1F("trueR_all","trueR_all",500,0,5000);
TH1F *trueR_par=new TH1F("trueR_par","trueR_par",500,0,5000);

TH1F *w_all=new TH1F("w_all","w_all",100,0,1);
TH1F *w_par=new TH1F("w_par","w_par",100,0,1);

TH1F *nchnl_all=new TH1F("nchnl_all","nchnl_all",17,-0.5,16.5);
TH1F *nchnl_par=new TH1F("nchnl_par","nchnl_par",17,-0.5,16.5);

TH1F *usedChan_all=new TH1F("usedChan_all","usedChan_all",16,-0.5,15.5);
TH1F *usedChan_par=new TH1F("usedChan_par","usedChan_par",16,-0.5,15.5);

//TH1F *coherence_all=new TH1F("coherence_all","coherence_all",100,0,50);
TH1F *coherence_all=new TH1F("coherence_all","coherence_all",1000,0,1);
TH1F *coherence_par=new TH1F("coherence_par","coherence_par",100,0,50);
TH1F *coherence_rf  =new TH1F("coherence_rf","coherence_rf",1000,0,1);
TH1F *coherence_cal =new TH1F("coherence_cal","coherence_cal",1000,0,1);
TH1F *coherence_soft=new TH1F("coherence_soft","coherence_soft",1000,0,1);

TH1F *coherence25Hist=new TH1F("coherence25Hist","coherence25Hist",1000,0,1);
TH1F *coherence10Hist=new TH1F("coherence10Hist","coherence10Hist",1000,0,1);
TH1F *coherence5Hist=new TH1F("coherence5Hist","coherence5Hist",1000,0,1);
TH1F *coherence2Hist=new TH1F("coherence2Hist","coherence2Hist",1000,0,1);

TH1F *coherence2_rf  =new TH1F("coherence2_rf","coherence2_rf",1000,0,1);
TH1F *coherence2_cal =new TH1F("coherence2_cal","coherence2_cal",1000,0,1);
TH1F *coherence2_soft=new TH1F("coherence2_soft","coherence2_soft",1000,0,1);


TH1F *coherence2Hist_nGoodChanNorm=new TH1F("coherence2Hist_nGoodChanNorm","coherence2Hist_nGoodChanNorm",1000,0,1);

TH1F *dCoherence25Hist=new TH1F("dCoherence25Hist","dCoherence25Hist",1000,-0.5,0.5);
TH1F *dCoherence10Hist=new TH1F("dCoherence10Hist","dCoherence10Hist",1000,-0.5,0.5);
TH1F *dCoherence5Hist=new TH1F("dCoherence5Hist","dCoherence5Hist",1000,-0.5,0.5);
TH1F *dCoherence2Hist=new TH1F("dCoherence2Hist","dCoherence2Hist",1000,-0.5,0.5);

TH1F *zenDist = new TH1F("zenDist","zenDist",375,0,180);
TH1F *aziDist = new TH1F("aziDist","aziDist",750,0,360);
TH2F *zenDist_coherence2 = new TH2F("zenDist_coherence2","zenDist_coherence2",1000,0,1,375,0,180);
TH2F *aziDist_coherence2 = new TH2F("aziDist_coherence2","aziDist_coherence2",1000,0,1,750,0,360);


TH1F *dZen = new TH1F("dZen","dZen",750,-180,180);
TH1F *dAzi = new TH1F("dAzi","dAzi",1500,-360,360);
TH1F *dR   = new TH1F("dR",  "dR",  1000, -5000, 5000);
TH1F *dZen_weight = new TH1F("dZen_weight","dZen_weight",750,-180,180);
TH1F *dAzi_weight = new TH1F("dAzi_weight","dAzi_weight",1500,-360,360);
TH1F *dR_weight   = new TH1F("dR_weight","dR_weight", 1000, -5000, 5000);

TH1F *dZen_weight_nLayer25 = new TH1F("dZen_weight_nLayer25","dZen_weight_nLayer25",750,-180,180);
TH1F *dZen_weight_nLayer10 = new TH1F("dZen_weight_nLayer10","dZen_weight_nLayer10",750,-180,180);
TH1F *dZen_weight_nLayer5 = new TH1F("dZen_weight_nLayer5","dZen_weight_nLayer5",750,-180,180);
TH1F *dZen_weight_nLayer2 = new TH1F("dZen_weight_nLayer2","dZen_weight_nLayer2",750,-180,180);

TH1F *dAzi_weight_nLayer25 = new TH1F("dAzi_weight_nLayer25","dAzi_weight_nLayer25",1500,-360,360);
TH1F *dAzi_weight_nLayer10 = new TH1F("dAzi_weight_nLayer10","dAzi_weight_nLayer10",1500,-360,360);
TH1F *dAzi_weight_nLayer5 = new TH1F("dAzi_weight_nLayer5","dAzi_weight_nLayer5",1500,-360,360);
TH1F *dAzi_weight_nLayer2 = new TH1F("dAzi_weight_nLayer2","dAzi_weight_nLayer2",1500,-360,360);

TH1F *dR_weight_nLayer25 = new TH1F("dR_weight_nLayer25","dR_weight_nLayer25",1000,-5000,5000);
TH1F *dR_weight_nLayer10 = new TH1F("dR_weight_nLayer10","dR_weight_nLayer10",1000,-5000,5000);
TH1F *dR_weight_nLayer5 = new TH1F("dR_weight_nLayer5","dR_weight_nLayer5",1000,-5000,5000);
TH1F *dR_weight_nLayer2 = new TH1F("dR_weight_nLayer2","dR_weight_nLayer2",1000,-5000,5000);

TH1F *dRecoZen_weight_nLayer25 = new TH1F("dRecoZen_weight_nLayer25","dRecoZen_weight_nLayer25",750,-180,180);
TH1F *dRecoZen_weight_nLayer10 = new TH1F("dRecoZen_weight_nLayer10","dRecoZen_weight_nLayer10",750,-180,180);
TH1F *dRecoZen_weight_nLayer5 = new TH1F("dRecoZen_weight_nLayer5","dRecoZen_weight_nLayer5",750,-180,180);
TH1F *dRecoZen_weight_nLayer2 = new TH1F("dRecoZen_weight_nLayer2","dRecoZen_weight_nLayer2",750,-180,180);

TH1F *dRecoAzi_weight_nLayer25 = new TH1F("dRecoAzi_weight_nLayer25","dRecoAzi_weight_nLayer25",1500,-360,360);
TH1F *dRecoAzi_weight_nLayer10 = new TH1F("dRecoAzi_weight_nLayer10","dRecoAzi_weight_nLayer10",1500,-360,360);
TH1F *dRecoAzi_weight_nLayer5 = new TH1F("dRecoAzi_weight_nLayer5","dRecoAzi_weight_nLayer5",1500,-360,360);
TH1F *dRecoAzi_weight_nLayer2 = new TH1F("dRecoAzi_weight_nLayer2","dRecoAzi_weight_nLayer2",1500,-360,360);

TH1F *dRecoR_weight_nLayer25 = new TH1F("dRecoR_weight_nLayer25","dRecoR_weight_nLayer25",1000,-5000,5000);
TH1F *dRecoR_weight_nLayer10 = new TH1F("dRecoR_weight_nLayer10","dRecoR_weight_nLayer10",1000,-5000,5000);
TH1F *dRecoR_weight_nLayer5 = new TH1F("dRecoR_weight_nLayer5","dRecoR_weight_nLayer5",1000,-5000,5000);
TH1F *dRecoR_weight_nLayer2 = new TH1F("dRecoR_weight_nLayer2","dRecoR_weight_nLayer2",1000,-5000,5000);

TH2F *dRecoZen_dCoherence_weight_nLayer25 = new TH2F("dRecoZen_dCoherence_weight_nLayer25","dRecoZen_dCoherence_weight_nLayer25",750,-180,180,1000,-0.5,0.5);
TH2F *dRecoZen_dCoherence_weight_nLayer10 = new TH2F("dRecoZen_dCoherence_weight_nLayer10","dRecoZen_dCoherence_weight_nLayer10",750,-180,180,1000,-0.5,0.5);
TH2F *dRecoZen_dCoherence_weight_nLayer5 = new TH2F("dRecoZen_dCoherence_weight_nLayer5","dRecoZen_dCoherence_weight_nLayer5",750,-180,180,1000,-0.5,0.5);
TH2F *dRecoZen_dCoherence_weight_nLayer2 = new TH2F("dRecoZen_dCoherence_weight_nLayer2","dRecoZen_dCoherence_weight_nLayer2",750,-180,180,1000,-0.5,0.5);

TH2F *dRecoAzi_dCoherence_weight_nLayer25 = new TH2F("dRecoAzi_dCoherence_weight_nLayer25","dRecoAzi_dCoherence_weight_nLayer25",1500,-360,360,1000,-0.5,0.5);
TH2F *dRecoAzi_dCoherence_weight_nLayer10 = new TH2F("dRecoAzi_dCoherence_weight_nLayer10","dRecoAzi_dCoherence_weight_nLayer10",1500,-360,360,1000,-0.5,0.5);
TH2F *dRecoAzi_dCoherence_weight_nLayer5 = new TH2F("dRecoAzi_dCoherence_weight_nLayer5","dRecoAzi_dCoherence_weight_nLayer5",1500,-360,360,1000,-0.5,0.5);
TH2F *dRecoAzi_dCoherence_weight_nLayer2 = new TH2F("dRecoAzi_dCoherence_weight_nLayer2","dRecoAzi_dCoherence_weight_nLayer2",1500,-360,360,1000,-0.5,0.5);

TH2F *dRecoR_dCoherence_weight_nLayer25 = new TH2F("dRecoR_dCoherence_weight_nLayer25","dRecoR_dCoherence_weight_nLayer25",1000,-5000,5000,1000,-0.5,0.5);
TH2F *dRecoR_dCoherence_weight_nLayer10 = new TH2F("dReco_dCoherenceR_weight_nLayer10","dRecoR_dCoherence_weight_nLayer10",1000,-5000,5000,1000,-0.5,0.5);
TH2F *dRecoR_dCoherence_weight_nLayer5 = new TH2F("dRecoR_dCoherence_weight_nLayer5","dRecoR_dCoherence_weight_nLayer5",1000,-5000,5000,1000,-0.5,0.5);
TH2F *dRecoR_dCoherence_weight_nLayer2 = new TH2F("dRecoR_dCoherence_weight_nLayer2","dRecoR_dCoherence_weight_nLayer2",1000,-5000,5000,1000,-0.5,0.5);


TH1F *dZen_cut = new TH1F("dZen_cut","dZen_cut",750,-180,180);
TH1F *dAzi_cut = new TH1F("dAzi_cut","dAzi_cut",1500,-360,360);
TH1F *dR_cut   = new TH1F("dR_cut",  "dR_cut",  1000, -5000, 5000);
TH1F *dZen_weight_cut = new TH1F("dZen_weight_cut","dZen_weight_cut",750,-180,180);
TH1F *dAzi_weight_cut = new TH1F("dAzi_weight_cut","dAzi_weight_cut",1500,-360,360);
TH1F *dR_weight_cut   = new TH1F("dR_weight_cut","dR_weight_cut", 1000, -5000, 5000);

TH1F *zen_all  = new TH1F("zen_all","Reco zenith (All)",450,-90,90);
TH1F *zen_rf   = new TH1F("zen_rf","Reco zenith (RF)",450,-90,90);
TH1F *zen_cal  = new TH1F("zen_cal","Reco zenith (Cal)",450,-90,90);
TH1F *zen_soft = new TH1F("zen_soft","Reco zenith (Soft)",450,-90,90);

TH1F *azi_all  = new TH1F("azi_all","Reco azimuth (All)",900,0,360);
TH1F *azi_rf   = new TH1F("azi_rf","Reco azimuth (RF)",900,0,360);
TH1F *azi_cal  = new TH1F("azi_cal","Reco azimuth (Cal)",900,0,360);
TH1F *azi_soft = new TH1F("azi_soft","Reco azimuth (Soft)",900,0,360);

TH2F *zen_azi_rf   =new TH2F("zen_azi_rf","Reco direction (RF)",900,0,360,450,-90,90);
TH2F *zen_azi_cal  =new TH2F("zen_azi_cal","Reco direction (Cal)",900,0,360,450,-90,90);
TH2F *zen_azi_soft =new TH2F("zen_azi_soft","Reco direction (Soft)",900,0,360,450,-90,90);

TH2F *zen_coherence_rf   =new TH2F("zen_coherence_rf","Reco zenith vs Coherence (RF)",1000,0,1,360,-90,90);
TH2F *zen_coherence_cal  =new TH2F("zen_coherence_cal","Reco zenith vs Coherence (Cal)",1000,0,1,360,-90,90);
TH2F *zen_coherence_soft =new TH2F("zen_coherence_soft","Reco zenith vs Coherence (Soft)",1000,0,1,360,-90,90);

TH2F *azi_coherence_rf   =new TH2F("azi_coherence_rf","Reco azimuth vs Coherence (RF)",1000,0,1,750,0,360);
TH2F *azi_coherence_cal  =new TH2F("azi_coherence_cal","Reco azimuth vs Coherence (Cal)",1000,0,1,750,0,360);
TH2F *azi_coherence_soft =new TH2F("azi_coherence_soft","Reco azimuth vs Coherence (Soft)",1000,0,1,750,0,360);

TH2F *zen_snr_rf   =new TH2F("zen_snr_rf","Reco zenith vs SNR (RF)",400,0,40,360,-90,90);
TH2F *zen_snr_cal  =new TH2F("zen_snr_cal","Reco zenith vs SNR (Cal)",400,0,40,360,-90,90);
TH2F *zen_snr_soft =new TH2F("zen_snr_soft","Reco zenith vs SNR (Soft)",400,0,40,360,-90,90);

TH2F *azi_snr_rf   =new TH2F("azi_snr_rf","Reco azimuth vs SNR (RF)",400,0,40,750,0,360);
TH2F *azi_snr_cal  =new TH2F("azi_snr_cal","Reco azimuth vs SNR (Cal)",400,0,40,750,0,360);
TH2F *azi_snr_soft =new TH2F("azi_snr_soft","Reco azimuth vs SNR (Soft)",400,0,40,750,0,360);

TH2F *coherence_snr_all  = new TH2F("coherence_snr_all","coherence_snr_all",400,0,40,1000,0,1);
TH2F *coherence_snr_rf   = new TH2F("coherence_snr_rf","coherence_snr_rf",400,0,40,1000,0,1);
TH2F *coherence_snr_cal  = new TH2F("coherence_snr_cal","coherence_snr_cal",400,0,40,1000,0,1);
TH2F *coherence_snr_soft = new TH2F("coherence_snr_soft","coherence_snr_soft",400,0,40,1000,0,1);

TH2F *coherence_snr_cw  = new TH2F("coherence_snr_cw","coherence_snr_cw",400,0,40,1000,0,1);

TH2F *coherence2_snr_all  = new TH2F("coherence2_snr_all","coherence2_snr_all",400,0,40,1000,0,1);
TH2F *coherence2_snr_rf   = new TH2F("coherence2_snr_rf","coherence2_snr_rf",400,0,40,1000,0,1);
TH2F *coherence2_snr_cal  = new TH2F("coherence2_snr_cal","coherence2_snr_cal",400,0,40,1000,0,1);
TH2F *coherence2_snr_soft = new TH2F("coherence2_snr_soft","coherence2_snr_soft",400,0,40,1000,0,1);

TH1F *impulsivityHist = new TH1F("impulsivityHist","impulsivityHist", 1000, 0, 1);

//TH2F *coherence_
TH2F *bipolarRatio_dT = new TH2F("bipolarRatio_dT","bipolarRatio_dT",800,0,800,1000,0,1);
TH2F *impulsivity_dT = new TH2F("impulsivity_dT","impulsivity_dT",800,0,800,1000,0,1);
TH2F *impulsivity_bipolarRatio = new TH2F("impulsivity_bipolarRatio","impulsivity_bipolarRatio",1000,0,1,1000,0,1);

TH1F *dTHist = new TH1F("dTHist","dTHist",800,0,800);
TH1F *impAvgHist = new TH1F("impAvgHist","impAvgHist",1000,0,1);
TH1F *impPassThresHist = new TH1F("impPassThresAvgHist","impPassThresAvgHist",1000,0,1);


TH2F *dZen_coherence=new TH2F("dZen_coherence","dZen_coherence",100,0,1,750,-180,180);
TH2F *dZen_nchnl    =new TH2F("dZen_nchnl","dZen_nchnl",17,-0.5,16.5,750,-180,180);
TH2F *dZen_trueR    =new TH2F("dZen_trueR","dZen_trueR",500,0,5000,750,-180,180);
TH2F *dZen_recoZen  =new TH2F("dZen_recoZen","dZen_recoZen",180,0,180,750,-180,180);
TH2F *dZen_trueZen  =new TH2F("dZen_trueZen","dZen_trueZen",180,0,180,750,-180,180);
TH2F *dZen_recoAzi  =new TH2F("dZen_recoAzi","dZen_recoAzi",360,0,360,750,-180,180);
TH2F *dZen_trueAzi  =new TH2F("dZen_trueAzi","dZen_trueAzi",360,0,360,750,-180,180);

TH2F *dZen_coherence_nLayer25=new TH2F("dZen_coherence_nLayer25","dZen_coherence_nLayer25",100,0,1,750,-180,180);
TH2F *dZen_coherence_nLayer10=new TH2F("dZen_coherence_nLayer10","dZen_coherence_nLayer10",100,0,1,750,-180,180);
TH2F *dZen_coherence_nLayer5=new TH2F("dZen_coherence_nLayer5","dZen_coherence_nLayer5",100,0,1,750,-180,180);
TH2F *dZen_coherence_nLayer2=new TH2F("dZen_coherence_nLayer2","dZen_coherence_nLayer2",100,0,1,750,-180,180);

TH2F *dAzi_coherence=new TH2F("dAzi_coherence","dAzi_coherence",100,0,1,1500,-360,360);
TH2F *dAzi_nchnl    =new TH2F("dAzi_nchnl","dAzi_nchnl",17,-0.5,16.5,1500,-360,360);
TH2F *dAzi_trueR    =new TH2F("dAzi_trueR","dAzi_trueR",500,0,5000,1500,-360,360);
TH2F *dAzi_recoZen  =new TH2F("dAzi_recoZen","dAzi_recoZen",180,0,180,1500,-360,360);
TH2F *dAzi_trueZen  =new TH2F("dAzi_trueZen","dAzi_trueZen",180,0,180,1500,-360,360);
TH2F *dAzi_recoAzi  =new TH2F("dAzi_recoAzi","dAzi_recoAzi",360,0,360,1500,-360,360);
TH2F *dAzi_trueAzi  =new TH2F("dAzi_trueAzi","dAzi_trueAzi",360,0,360,1500,-360,360);

TH2F *dAzi_coherence_nLayer25=new TH2F("dAzi_coherence_nLayer25","dAzi_coherence_nLayer25",100,0,1,1500,-360,360);
TH2F *dAzi_coherence_nLayer10=new TH2F("dAzi_coherence_nLayer10","dAzi_coherence_nLayer10",100,0,1,1500,-360,360);
TH2F *dAzi_coherence_nLayer5=new TH2F("dAzi_coherence_nLayer5","dAzi_coherence_nLayer5",100,0,1,1500,-360,360);
TH2F *dAzi_coherence_nLayer2=new TH2F("dAzi_coherence_nLayer2","dAzi_coherence_nLayer2",100,0,1,1500,-360,360);

TH2F *dZen_snr = new TH2F("dZen_snr","dZen_snr",400, 0, 40, 750, -180, 180);
TH2F *dAzi_snr = new TH2F("dAzi_snr","dAzi_snr", 400, 0, 40, 1500, -360, 360);

TH2F *dAlpha_coherence=new TH2F("dAlpha_coherence","dAlpha_coherence",50,0,25,180,0,180);
TH2F *dAlpha_nchnl    =new TH2F("dAlpha_nchnl","dAlpha_nchnl",17,-0.5,16.5,180,0,180);
TH2F *dAlpha_trueR    =new TH2F("dAlpha_trueR","dAlpha_trueR",500,0,5000,180,0,180);
TH2F *dAlpha_recoZen  =new TH2F("dAlpha_recoZen","dAlpha_recoZen",180,0,180,180,0,180);
TH2F *dAlpha_trueZen  =new TH2F("dAlpha_trueZen","dAlpha_trueZen",180,0,180,180,0,180);
TH2F *dAlpha_recoAzi  =new TH2F("dAlpha_recoAzi","dAlpha_recoAzi",360,0,360,180,0,180);
TH2F *dAlpha_trueAzi  =new TH2F("dAlpha_trueAzi","dAlpha_trueAzi",360,0,360,180,0,180);


TH2F *dZen_dAzi     =new TH2F("dZen_dAzi","dZen_dAzi",1500,-360,360,750,-180,180);
TH2F *dZen_dR       =new TH2F("dZen_dR"  ,"dZen_dR",  1000,-5000,5000,750,-180,180);
TH2F *dAzi_dR       =new TH2F("dAzi_dR"  ,"dAzi_dR",  1000,-5000,5000,1500,-360,360);
TH2F *dAlpha_dR     =new TH2F("dAlpha_dR","dAlpha_dR",1000,-5000,5000,180,0,180);

TH2F *coherence_nchnl = new TH2F("coherence_nchnl","coherence_nchnl",17,-0.5,16.5,1000,0,1);
TH2F *coherence_trueR = new TH2F("coherence_trueR","coherence_trueR",500,0,5000,1000,0,1);
TH2F *nchnl_trueR     = new TH2F("nchnl_trueR","nchnl_trueR",500,0,5000,17,-0.5,16.5);

TH2F *dR_coherence=new TH2F("dR_coherence","dR_coherence",100,0,1,1000,-5000,5000);
TH2F *dR_nchnl    =new TH2F("dR_nchnl","dR_nchnl",17,-0.5,16.5,1000,-5000,5000);
TH2F *dR_trueR    =new TH2F("dR_trueR","dR_trueR",10,0,5000,500,-5000,5000);
TH2F *dRRatio_trueR=new TH2F("dRRatio_trueR","dRRatio_trueR",10,0,5000,50,0,1);
TH2F *dR_recoR    =new TH2F("dR_recoR","dR_recoR",25,100,5100,1000,-5000,5000);
TH2F *dR_recoZen  =new TH2F("dR_recoZen","dR_recoZen",180,0,180,1000,-5000,5000);
TH2F *dR_trueZen  =new TH2F("dR_trueZen","dR_trueZen",180,0,180,1000,-5000,5000);
TH2F *dR_recoAzi  =new TH2F("dR_recoAzi","dR_recoAzi",360,0,360,1000,-5000,5000);
TH2F *dR_trueAzi  =new TH2F("dR_trueAzi","dR_trueAzi",360,0,360,1000,-5000,5000);

//TGraph2D *topMaxPixGraph = new TGraph2D(); //Have not come up with a good way to use dummyData->topN while avoiding memory leaks
//TGraph2D *maxPixEachLayerGraph = new TGraph2D();
TPolyMarker3D *topMaxPixGraph;//       = new TPolyMarker3D(settings->topN/*NULL*/, 6);
TPolyMarker3D *trueVertexGraph;
TPolyMarker3D *recoVertexGraph;
//TH3F          *topMaxPixHist        = new TH3F("topMaxPixHist","topMaxPixHist",1000,-5000,5000,1000,-5000,5000,290,-3000,100);
TH2F          *topMaxPixHist        = new TH2F("topMaxPixHist","topMaxPixHist",100,0,5000,250,0,10);
TLine trueDistanceLine;
TPolyMarker3D *maxPixEachLayerGraph = new TPolyMarker3D(onion.nLayer, 3);
TH3D *hist = new TH3D("hist","hist",10,-5000,5000,10,-5000,5000,5,-3000,100);
//TH3D *hist = new TH3D("hist","hist",10,-1000,1000,10,-1000,1000,5,-3000,100);
//maxPixEachLayerGraph->SetHistogram(hist);

TH1F *statRReco_weight     = new TH1F("statRReco_weight","statRReco_weight",1000,-5000,5000);
TH1F *statRReco_max        = new TH1F("statRReco_max","statRReco_max",1000,-5000,5000);
TH1F *statRReco_eventStack = new TH1F("statRReco_eventStack","statRReco_eventStack",1000,-5000,5000);


//Define alpha as space angle between reco and true directionsw
TProfile *dAlphaRProf = new TProfile("dAlphaRProf","Space Angle Difference Profile", 50, 0, 5000, 0, 180, "");
TProfile *dZenRProf = new TProfile("dZenRProf","Zenith Reconsturction Profile", 50, 0, 5000, -180, 180, "");
TProfile *dAziRProf   = new TProfile("dAziRProf","Azimuth Reconsturction Profile", 50, 0, 5000, -360, 360, "");
TProfile *dRRProf     = new TProfile("dRRProf","Vertex Distance Reconsturction Profile", 50, 0, 5000, -5000, 5000, "");
TH2F *dStatR_coherence=new TH2F("dStatR_coherence","dStatR_coherence",50,0,25,1000,-5000,5000);
TH2F *dStatR_nchnl    =new TH2F("dStatR_nchnl","dStatR_nchnl",17,-0.5,16.5,1000,-5000,5000);
TH2F *dStatR_trueR    =new TH2F("dStatR_trueR","dStatR_trueR",500,0,5000,1000,-5000,5000);
TH1F *dStatR[5];
float rArray[5]={5000,4000,3000,2000,1000};//,800,600,400,200};
char histName[200];
for(int j=0; j<5; j++){
   sprintf(histName, "dStatR_%d", (int)rArray[j]);
   dStatR[j] = new TH1F(histName, histName, 100,-5000,5000);
}

TH1F *maxPixEachLayerZen/*=new TH1F("maxPixEachLayerZen","maxPixEachLayerZen",450,-90,90)*/;
TH1F *maxPixEachLayerAzi/*=new TH1F("maxPixEachLayerAzi","maxPixEachLayerAzi",900,0,360)*/;
TH1F *maxPixEachLayerZenSpread=new TH1F("maxPixEachLayerZenSpread","maxPixEachLayerZenSpread",225,0,90);
TH1F *maxPixEachLayerAziSpread=new TH1F("maxPixEachLayerAziSpread","maxPixEachLayerAziSpread",900,0,360);
TH1F *topMaxPixZen/*=new TH1F("maxPixEachLayerZen","maxPixEachLayerZen",450,-90,90)*/;
TH1F *topMaxPixAzi/*=new TH1F("maxPixEachLayerAzi","maxPixEachLayerAzi",900,0,360)*/;
TH1F *topMaxPixZenSpread=new TH1F("topMaxPixZenSpread","topMaxPixZenSpread",225,0,90);
TH1F *topMaxPixAziSpread=new TH1F("topMaxPixAziSpread","topMaxPixAziSpread",900,0,360);

TH1F *recoRecAngleHist[16], *trueRecAngleHist[16], *dRecAngle[16], *dLauAngle[16];
for(int i=0; i<16; i++){

   sprintf(histName, "recoRecAngleHist_%d", i);
   recoRecAngleHist[i] = new TH1F(histName, histName, 450,0,180);
   sprintf(histName, "trueRecAngleHist_%d", i);
   trueRecAngleHist[i] = new TH1F(histName, histName, 450,0,180);
   sprintf(histName, "dRecAngle_%d", i);
   dRecAngle[i] = new TH1F(histName, histName, 900,-180,180);
   sprintf(histName, "dLauAngle_%d", i);
   dLauAngle[i] = new TH1F(histName, histName, 900,-180,180);


}


int plotCount=0;

TH2F *coherence2_maxPixEachLayerAziSpread = new TH2F("coherence2_maxPixEachLayerAziSpread","coherence2_maxPixEachLayerAziSpread",900,0,360,1000,0,0.5);
TH2F *coherence2_maxPixEachLayerZenSpread = new TH2F("coherence2_maxPixEachLayerZenSpread","coherence2_maxPixEachLayerZenSpread",225,0,90,1000,0,0.5);
TH2F *maxPixEachLayerZenAziSpread = new TH2F("maxPixEachLayerZenAziSpread","maxPixEachLayerZenAziSpread",900,0,360,225,0,90);

TH2F *coherence2_topMaxPixAziSpread = new TH2F("coherence2_topMaxPixAziSpread","coherence2_topMaxPixAziSpread",900,0,360,1000,0,0.5);
TH2F *coherence2_topMaxPixZenSpread = new TH2F("coherence2_topMaxPixZenSpread","coherence2_topMaxPixZenSpread",225,0,90,1000,0,0.5);
TH2F *topMaxPixZenAziSpread = new TH2F("topMaxPixZenAziSpread","topMaxPixZenAziSpread",900/5,0,360,225/5,0,90);

TH2F *coherence2_thetaRMSCutXing = new TH2F("coherence2_thetaRMSCutXing","coherence2_thetaRMSCutXing",500,0,500,1000,0,0.5);
TH2F *coherence2_phiRMSCutXing = new TH2F("coherence2_phiRMSCutXing","coherence2_phiRMSCutXing",500,0,500,1000,0,0.5);
TH2F *topMaxPixZenSpread_thetaRMSCutXing = new TH2F("topMaxPixZenSpread_thetaRMSCutXing","topMaxPixZenSpread_thetaRMSCutXing",500,0,500,225,0,90);
TH2F *topMaxPixAziSpread_phiRMSCutXing = new TH2F("topMaxPixAziSpread_phiRMSCutXing","topMaxPixAziSpread_phiRMSCutXing",500,0,500,900,0,360);
TH1F *thetaRatioHist = new TH1F("thetaRatioHist","thetaRatioHist",500,0,500);
TH1F *phiRatioHist = new TH1F("phiRatioHist","phiRatioHist",500,0,500);
TGraph *grThetaRMS, *grPhiRMS;
TGraph *grGrowingThetaRMS, *grGrowingPhiRMS;
int *thetaIdx, *phiIdx;
float *dRMSThetaVec, *dRMSPhiVec;
vector<float> thetaVec, phiVec;
vector<float> growingThetaVec, growingPhiVec;
bool thetaRatioDone, phiRatioDone;
float rms;
const float rmsCut = 1.f;
float thetaXing, phiXing;

int nchnl;
int parCnt=0;
int thetaCnt, phiCnt, thetaPhiCnt;
thetaCnt = phiCnt = thetaPhiCnt = 0;
float r, theta, phi;
float alpha; //space angle between true and reco direction
float statRecoRadius;

float *value, *radii, *pixCount;
char filename[200];
TCanvas *cvs;//("cvs","cvs",800,600);

int rfEventCount, calEventCount, softEventCount;
rfEventCount = calEventCount = softEventCount = 0;

//int entry = atoi(1argv[1] );
float p = 0;

float coherence25[25], coherence10[10], coherence5[5], coherence2[2];
int idx25[25], idx10[10], idx5[5], idx2[2];
int pixIdx25[25], pixIdx10[10], pixIdx5[5], pixIdx2[2];

float recoRecAngle[16], recoLauAngle[16], trueRecAngle[16], trueLauAngle[16];

TH1F *snrGradient_group3 = new TH1F("snrGradient_group3","snrGradient_group3",100,0,1);
TH1F *snrGradient_group4 = new TH1F("snrGradient_group4","snrGradient_group4",100,0,1);
TH1F *snrGradient_7_0 = new TH1F("snrGradient_7_0","snrGradient_7_0",100,0,1);
TH1F *snrGradient_6_1 = new TH1F("snrGradient_6_1","snrGradient_6_1",100,0,1);
float snrArray[16];
int index[16];

//TH1F *timeSeuqenceHist = new TH1F("timeSequenceHist","timeSequenceHist",3e4,0,1e4);
TH1F *iterCWCountHist = new TH1F("iterCWCountHist","iterCWCountHist",11,-0.1,10.5);

TH1F *avgPowerRatioHist[16];

for(int ch=0; ch<16; ch++){
   sprintf(histName,"avgPowerRatioHist_%d",ch);
   avgPowerRatioHist=new TH1F(histName,histName,1000,0,20000);
}

double fftRes;
/*if(type<=3)*/ fftRes = 1/(379e-9)/1e6;
//else        fftRes = 1/(499e-9)/1e6;
//cout<<"fftRes: "<<fftRes<<endl;

for(int entry=0; entry<Nentries; entry++){

   if(Nentries > 100) {  if(  entry % (Nentries/100) == 0  ){ cout<<"Progess: "<<entry / (Nentries/100) <<"%\n"; } }
   dataTree->GetEntry(entry);
   //cout<<"eventTrigType: "<<dummyData->eventTrigType<<endl;
   if(dummyData->eventTrigType == 0) rfEventCount++;
   else if (dummyData->eventTrigType == 1) calEventCount++;
   else if (dummyData->eventTrigType == 2) softEventCount++;
   else { cerr<<"Event "<<entry<<" eventTrigType undefined! Skipping...\n"; continue; }

   //topMaxPixGraph = new TGraph2D(dummyData->topN);
   topMaxPixGraph  = new TPolyMarker3D(settings->topN/*NULL*/, 6);
   trueVertexGraph = new TPolyMarker3D(1,6);
   recoVertexGraph = new TPolyMarker3D(1,6);

   dRMSThetaVec = (float*)calloc(settings->topN, sizeof(float));
   dRMSPhiVec = (float*)calloc(settings->topN, sizeof(float));
   grThetaRMS = new TGraph();
   grPhiRMS   = new TGraph();
   grGrowingThetaRMS = new TGraph();
   grGrowingPhiRMS   = new TGraph();
   thetaIdx = (int*)calloc(settings->topN, sizeof(int));
   phiIdx   = (int*)calloc(settings->topN, sizeof(int));

   bool isVpolCW, isHpolCW, isXpolCW;
   isVpolCW = isHpolCW = isXpolCW = false;

   /* Now inspect why theta reco fails sometimes */

   trueR_all->Fill(dummyData->trueRadius);
   w_all->Fill(dummyData->weight);
   coherence_all->Fill(dummyData->maxPixCoherence);
   if(dummyData->eventTrigType == 0) coherence_rf->Fill(dummyData->maxPixCoherence);
   else if(dummyData->eventTrigType == 1) coherence_cal->Fill(dummyData->maxPixCoherence);
   else coherence_soft->Fill(dummyData->maxPixCoherence);


//   for(int ch=0; ch<16; ch++) snrArray[ch] = dummyData->channelInWindowSNR[ch];
//   TMath::Sort(16, snrArray, index);
//   float top3Avg     = (snrArray[index[0]] + snrArray[index[1]] + snrArray[index[2]]) / 3.;
//   float bottom3Avg = (snrArray[index[5]] + snrArray[index[6]] + snrArray[index[7]]) / 3.;
//   float top4Avg     = (snrArray[index[0]] + snrArray[index[1]] + snrArray[index[2]] + snrArray[index[3]]) / 4.;
//   float bottom4Avg  = (snrArray[index[4]] + snrArray[index[5]] + snrArray[index[6]] + snrArray[index[7]]) / 4.;
//
//   snrGradient_group3->Fill(bottom3Avg/top3Avg, dummyData->weight);
//   snrGradient_group4->Fill(bottom4Avg/top4Avg, dummyData->weight);
//   snrGradient_7_0->Fill(snrArray[index[7]]/snrArray[index[0]], dummyData->weight);
//   snrGradient_6_1->Fill(snrArray[index[6]]/snrArray[index[1]], dummyData->weight);
//
//
//   nchnl=0;
//   for(int ch=0; ch<8; ch++){
//   nchnl+=dummyData->recoChan[ch];
//   if(dummyData->recoChan[ch]) usedChan_all->Fill(ch);
//   }
//   nchnl_all->Fill(nchnl);
//
//   if(dummyData->maxPixCoherence>0.1){
//   for(int i=0; i<16; i++){
//
//      recoRecAngle[i] = dummyData->recoRecAngle[i] /** TMath::RadToDeg()*/;
//      recoLauAngle[i] = dummyData->recoLauAngle[i] /** TMath::RadToDeg()*/;
//      trueRecAngle[i] = dummyData->trueRecAngle[i] /** TMath::RadToDeg()*/;
//      trueLauAngle[i] = dummyData->trueLauAngle[i] /** TMath::RadToDeg()*/;
//
//      recoRecAngleHist[i]->Fill(recoRecAngle[i]*TMath::RadToDeg(), dummyData->weight);
//      trueRecAngleHist[i]->Fill(trueRecAngle[i]*TMath::RadToDeg(), dummyData->weight);
//      dRecAngle[i]->Fill((recoRecAngle[i]-trueRecAngle[i])*TMath::RadToDeg());
//      dLauAngle[i]->Fill((recoLauAngle[i]-trueLauAngle[i])*TMath::RadToDeg());
//
//   }
//   }
   dZen->Fill(dummyData->recoZen - dummyData->trueZen);
   dAzi->Fill(dummyData->recoAzi - dummyData->trueAzi);
   dR->Fill(onion.getLayerRadius(dummyData->maxPixIdx) - dummyData->trueRadius);
   dZen_snr->Fill( dummyData->inWindowSNR_V, dummyData->recoZen - dummyData->trueZen);
   dAzi_snr->Fill( dummyData->inWindowSNR_V, dummyData->recoAzi - dummyData->trueAzi);
//   if(dummyData->maxPixCoherence>0.12){
//   dZen_weight->Fill(dummyData->recoZen - dummyData->trueZen, dummyData->weight);
//   dAzi_weight->Fill(dummyData->recoAzi - dummyData->trueAzi, dummyData->weight);
//   dR_weight->Fill(onion.getLayerRadius(dummyData->maxPixIdx) - dummyData->trueRadius, dummyData->weight);
//   }
   zen_all->Fill(90.f-dummyData->recoZen);
   azi_all->Fill(dummyData->recoAzi);
   coherence_snr_all->Fill(dummyData->unmodSNR, dummyData->maxPixCoherence);
   if(dummyData->eventTrigType == 0){ //RF

      zen_rf->Fill(90.f-dummyData->recoZen);
      azi_rf->Fill(dummyData->recoAzi);
      zen_azi_rf->Fill(dummyData->recoAzi, 90.f-dummyData->recoZen);
      zen_coherence_rf->Fill(dummyData->maxPixCoherence, 90.f-dummyData->recoZen);
      azi_coherence_rf->Fill(dummyData->maxPixCoherence, dummyData->recoAzi);
      zen_snr_rf->Fill(dummyData->unmodSNR, 90.f-dummyData->recoZen);
      azi_snr_rf->Fill(dummyData->unmodSNR, dummyData->recoAzi);
      coherence_snr_rf->Fill(dummyData->unmodSNR, dummyData->maxPixCoherence);

      bool isVpolCW, isHpolCW;
      int maxCountFreqBin_H, maxCountFreqBin_V;
      int cwBinThres = 3;
      int nonZeroCount = 0;
      double avgImpulsivity = 0.;

      if (isCW_coincidence(isVpolCW, isHpolCW, maxCountFreqBin_V, maxCountFreqBin_H, dummyData, cwBinThres)){
/*
         if (isVpolCW && isHpolCW){

            for(int ch=0; ch<16; ch++){
               if(dummyData->impulsivity[ch]>0){
                  avgImpulsivity += dummyData->impulsivity[ch];
                  nonZeroCount ++;
               }
            }

            if(nonZeroCount>0)
            avgImpulsivity /= (double)nonZeroCount;

         } else if (isVpolCW){

            for(int ch=0; ch<8; ch++){
               if(dummyData->impulsivity[ch]>0){
                  avgImpulsivity += dummyData->impulsivity[ch];
                  nonZeroCount ++;
               }
            }

            if(nonZeroCount>0)
            avgImpulsivity /= (double)nonZeroCount;

         } else {

               for(int ch=8; ch<16; ch++){
               if(dummyData->impulsivity[ch]>0){
                  avgImpulsivity += dummyData->impulsivity[ch];
                  nonZeroCount ++;
               }
            }

            if(nonZeroCount>0)
            avgImpulsivity /= (double)nonZeroCount;

         }
*/

         for(int ch=0; ch<8; ch++){
            if(dummyData->impulsivity[ch]>0){
               avgImpulsivity += dummyData->impulsivity[ch];
               nonZeroCount ++;
            }
         }

         if(nonZeroCount>0)
         avgImpulsivity /= (double)nonZeroCount;

         impulsivityHist->Fill(avgImpulsivity, dummyData->weight);



      }

       bool isCW = isCW_freqWindow(isVpolCW, isHpolCW, isXpolCW, dummyData, fftRes);

       if(isCW){
          iterCWCountHist->Fill(dummyData->cwIterCount, dummyData->weight);

          coherence_snr_cw->Fill(dummyData->inWindowSNR_V, (dummyData->maxPixCoherence>dummyData->maxPixCoherence2?dummyData->maxPixCoherence:dummyData->maxPixCoherence2), dummyData->weight);

          for(int ch=0; ch<8; ch++){
             if(dummyData->posPowerPeak[ch]>0){
                double avgPowerRatio = (dummyData->posPowerPeak[ch] - dummyData->negPowerPeak[ch]) / 2./*dummyData->posPowerPeak[ch]*/;
                avgPowerRatioHist[ch]->Fill(avgPowerRatio, dummyData->weight);
             }
          }

       }




/*
      double impSum=0.;
      double impPassThresSum=0.;
      int nNonZero=0;
      int nPassThres=0;

      for(int ch=0; ch<8; ch++){

      double bipolarRatio;
      double pos, neg;
      pos = dummyData->posPowerPeak[ch];
      neg = fabs(dummyData->negPowerPeak[ch]);

      if(pos > neg) bipolarRatio = neg/pos;
      else{ if(neg>0) bipolarRatio = pos/neg; else bipolarRatio = 0.; }

      double dT = fabs(dummyData->powerPeaksDeltaT[ch]);
      bipolarRatio_dT->Fill(dT, bipolarRatio, dummyData->weight);
      impulsivity_dT->Fill(dT, dummyData->impulsivity[ch], dummyData->weight);
      impulsivity_bipolarRatio->Fill(bipolarRatio, dummyData->impulsivity[ch], dummyData->weight);

      dTHist->Fill(dT, dummyData->weight);

      if(dummyData->slidingV2SNR[ch] > settings->nchnlThreshold){

         if(fabs(dummyData->impulsivity[ch]-0)>1e-9){
            impPassThresSum += dummyData->impulsivity[ch];
            nPassThres++;
         }

      }

      if(fabs(dummyData->impulsivity[ch]-0)>1e-9){
         impSum += dummyData->impulsivity[ch];
         nNonZero++;
      }

      }

      double impAvg          = impSum / (double)nNonZero;
      double impPassThresAvg = impPassThresSum / (double)nPassThres;

      impAvgHist->Fill(impAvg, dummyData->weight);
      impPassThresHist->Fill(impPassThresAvg, dummyData->weight);


   } else if(dummyData->eventTrigType == 1){ //Cal
      zen_cal->Fill(90.f-dummyData->recoZen);
      azi_cal->Fill(dummyData->recoAzi);
      zen_azi_cal->Fill(dummyData->recoAzi, 90.f-dummyData->recoZen);
      zen_coherence_cal->Fill(dummyData->maxPixCoherence, 90.f-dummyData->recoZen);
      azi_coherence_cal->Fill(dummyData->maxPixCoherence, dummyData->recoAzi);
      zen_snr_cal->Fill(dummyData->unmodSNR, 90.f-dummyData->recoZen);
      azi_snr_cal->Fill(dummyData->unmodSNR, dummyData->recoAzi);
      coherence_snr_cal->Fill(dummyData->unmodSNR, dummyData->maxPixCoherence);
   } else { //Soft
      zen_soft->Fill(90.f-dummyData->recoZen);
      azi_soft->Fill(dummyData->recoAzi);
      zen_azi_soft->Fill(dummyData->recoAzi, 90.f-dummyData->recoZen);
      zen_coherence_soft->Fill(dummyData->maxPixCoherence, 90.f-dummyData->recoZen);
      azi_coherence_soft->Fill(dummyData->maxPixCoherence, dummyData->recoAzi);
      zen_snr_soft->Fill(dummyData->unmodSNR, 90.f-dummyData->recoZen);
      azi_snr_soft->Fill(dummyData->unmodSNR, dummyData->recoAzi);
      coherence_snr_soft->Fill(dummyData->unmodSNR, dummyData->maxPixCoherence);
   }//end of if RF
*/
//
//   if(dummyData->maxPixCoherence > 5.f){
//
//   dZen_cut->Fill(dummyData->recoZen - dummyData->trueZen);
//   dAzi_cut->Fill(dummyData->recoAzi - dummyData->trueAzi);
//   dR_cut->Fill(onion.getLayerRadius(dummyData->maxPixIdx) - dummyData->trueRadius);
//   dZen_weight_cut->Fill(dummyData->recoZen - dummyData->trueZen, dummyData->weight);
//   dAzi_weight_cut->Fill(dummyData->recoAzi - dummyData->trueAzi, dummyData->weight);
//   dR_weight_cut->Fill(onion.getLayerRadius(dummyData->maxPixIdx) - dummyData->trueRadius, dummyData->weight);
//
//   }
//
//   dZen_coherence->Fill(dummyData->maxPixCoherence, dummyData->recoZen - dummyData->trueZen);
//   dZen_nchnl->Fill(nchnl, dummyData->recoZen - dummyData->trueZen);
//   dZen_trueR->Fill(dummyData->trueRadius, dummyData->recoZen - dummyData->trueZen);
//   dZen_recoZen->Fill(dummyData->recoZen, dummyData->recoZen - dummyData->trueZen);
//   dZen_trueZen->Fill(dummyData->trueZen, dummyData->recoZen - dummyData->trueZen);
//   dZen_recoAzi->Fill(dummyData->recoAzi, dummyData->recoZen - dummyData->trueZen);
//   dZen_trueAzi->Fill(dummyData->trueAzi, dummyData->recoZen - dummyData->trueZen);
//
//   dAzi_coherence->Fill(dummyData->maxPixCoherence, dummyData->recoAzi - dummyData->trueAzi);
//   dAzi_nchnl->Fill(nchnl, dummyData->recoAzi - dummyData->trueAzi);
//   dAzi_trueR->Fill(dummyData->trueRadius, dummyData->recoAzi - dummyData->trueAzi);
//   dAzi_recoZen->Fill(dummyData->recoZen, dummyData->recoAzi - dummyData->trueAzi);
//   dAzi_trueZen->Fill(dummyData->trueZen, dummyData->recoAzi - dummyData->trueAzi);
//   dAzi_recoAzi->Fill(dummyData->recoAzi, dummyData->recoAzi - dummyData->trueAzi);
//   dAzi_trueAzi->Fill(dummyData->trueAzi, dummyData->recoAzi - dummyData->trueAzi);
//
//   dZen_dAzi->Fill(dummyData->recoAzi - dummyData->trueAzi, dummyData->recoZen - dummyData->trueZen);
//
//   coherence_nchnl->Fill(nchnl, dummyData->maxPixCoherence);
//   coherence_trueR->Fill(dummyData->trueRadius, dummyData->maxPixCoherence);
//   nchnl_trueR->Fill(dummyData->trueRadius, nchnl);
//
//   dR_coherence->Fill(dummyData->maxPixCoherence, onion.getLayerRadius(dummyData->maxPixIdx) - dummyData->trueRadius);
//   dR_nchnl->Fill(nchnl, onion.getLayerRadius(dummyData->maxPixIdx) - dummyData->trueRadius);
//   if(dummyData->maxPixCoherence > 0.12)
//   dR_trueR->Fill(dummyData->trueRadius, (onion.getLayerRadius(dummyData->maxPixIdx) - dummyData->trueRadius));
//   dRRatio_trueR->Fill(dummyData->trueRadius, fabs(onion.getLayerRadius(dummyData->maxPixIdx) - dummyData->trueRadius)/5000.);
//   dR_recoR->Fill(onion.getLayerRadius(dummyData->maxPixIdx), onion.getLayerRadius(dummyData->maxPixIdx) - dummyData->trueRadius);
//   dR_recoZen->Fill(dummyData->recoZen, onion.getLayerRadius(dummyData->maxPixIdx) - dummyData->trueRadius);
//   dR_trueZen->Fill(dummyData->trueZen, onion.getLayerRadius(dummyData->maxPixIdx) - dummyData->trueRadius);
//   dR_recoAzi->Fill(dummyData->recoAzi, onion.getLayerRadius(dummyData->maxPixIdx) - dummyData->trueRadius);
//   dR_trueAzi->Fill(dummyData->trueAzi, onion.getLayerRadius(dummyData->maxPixIdx) - dummyData->trueRadius);
//
//
//   snprintf(filename,sizeof(char)*200,"topMaxPixZen_%d",entry);
//   topMaxPixZen=new TH1F(filename,filename,450,0,180);
//   snprintf(filename,sizeof(char)*200,"topMaxPixAzi_%d",entry);
//   topMaxPixAzi=new TH1F(filename,filename,900,0,360);
//
//   thetaVec.clear();
//   phiVec.clear();
//
//
//   for(int i=0; i</*dummyData*/settings->topN; i++){
//
//   //cout<<"topMaxPixIdx: "<<dummyData->topMaxPixIdx.at(i)<<endl;
//   r     = onion.getLayerRadius(dummyData->topMaxPixIdx.at(i));
//   theta = onion.getPointing(dummyData->topMaxPixIdx.at(i)).theta;
//   phi   = onion.getPointing(dummyData->topMaxPixIdx.at(i)).phi;
//
//   thetaVec.push_back(theta*TMath::RadToDeg());
//   phiVec.push_back(phi*TMath::RadToDeg());
//
//   topMaxPixZen->Fill(theta*TMath::RadToDeg());
//   topMaxPixAzi->Fill(phi  *TMath::RadToDeg());
//
//   //printf("i:%d r:%f theta:%f phi:%f\n",i,r,theta,phi);
//
//   //topMaxPixGraph->SetPoint(i, r*sin(theta)*cos(phi), r*sin(theta)*sin(phi), r*cos(theta));
//
//   //topMaxPixHist->Fill(r*sin(theta)*cos(phi), r*sin(theta)*sin(phi), r*cos(theta), topMaxPixCoherence->at(i));
//   /*
//   if( entry == atoi(argv[2]) ){
//   topMaxPixHist->Fill(r, getSpaceAngle(dummyData->trueZen*TMath::DegToRad(), dummyData->trueAzi*TMath::DegToRad(), theta, phi ) * TMath::RadToDeg(), topMaxPixCoherence->at(i));
//   trueDistanceLine=TLine(dummyData->trueRadius, 0., dummyData->trueRadius, 10.);
//   }*/
//
//   }//end of topN
//
//   topMaxPixZenSpread->Fill(topMaxPixZen->GetRMS());
//   topMaxPixAziSpread->Fill(topMaxPixAzi->GetRMS());
//   topMaxPixZenAziSpread->Fill(topMaxPixAzi->GetRMS(), topMaxPixZen->GetRMS());
//
//   cvs = new TCanvas("cvs","cvs",800,800);
//
//
//if(dummyData->eventTrigType == 0){
//
//   float fullThetaRMS = getRMS(thetaVec);
//   float fullPhiRMS   = getRMS(phiVec);
//
//   //printf("Zen rms: (%f,%f)\tAzi rms: (%f,%f)\n",topMaxPixZen->GetRMS(),fullThetaRMS,topMaxPixAzi->GetRMS(),fullPhiRMS);
//
//   for(int i=0; i<settings->topN; i++){
//
//   vector<float> newThetaVec(thetaVec);
//   vector<float> newPhiVec(phiVec);
//
//   newThetaVec.erase(newThetaVec.begin()+i);
//   newPhiVec.erase(newPhiVec.begin()+i);
//
//   if((int)newThetaVec.size() != settings->topN-1 || (int)newPhiVec.size() != settings->topN-1) cerr<<"Vector size mismatch!\n";
//
//   //dRMSThetaVec.push_back( fabs(fullThetaRMS - getRMS(newThetaVec)) );
//   //dRMSPhiVec.push_back(   fabs(fullPhiRMS   - getRMS(newPhiVec))   );
//   dRMSThetaVec[i] = fullThetaRMS - getRMS(newThetaVec);
//   dRMSPhiVec[i]   = fullPhiRMS   - getRMS(newPhiVec);
//
//   /* Deacllocating vectors */
//   vector<float>().swap(newThetaVec);
//   vector<float>().swap(newPhiVec);
//
//   }
//
//   TMath::Sort(settings->topN, dRMSThetaVec, thetaIdx, false); //sort in increasing order
//   TMath::Sort(settings->topN, dRMSPhiVec,   phiIdx  , false);
//
//   //vector<float> newThetaVec(thetaVec);
//   //vector<float> newPhiVec(phiVec);
//   growingThetaVec.clear(); growingPhiVec.clear();
//   thetaRatioDone = phiRatioDone = false;
//
//   growingThetaVec.push_back(thetaVec[thetaIdx[0]]);
//   growingThetaVec.push_back(thetaVec[thetaIdx[1]]);
//   growingPhiVec.push_back(phiVec[phiIdx[0]]);
//   growingPhiVec.push_back(phiVec[phiIdx[1]]);
//
//   rms = getRMS(growingThetaVec);
//   grGrowingThetaRMS->SetPoint(0,2,rms);
//   if(rms > rmsCut && thetaRatioDone == false){ thetaXing = 2.f; thetaRatioHist->Fill(2.f/*/(float)recoSettings->topN*/); thetaRatioDone = true; }
//   //if(entry==0) cout<<getRMS(growingPhiVec)<<endl;
//   rms = getRMS(growingPhiVec);
//   grGrowingPhiRMS->SetPoint(0,2,rms);
//   if(rms > rmsCut && phiRatioDone == false){ phiXing = 2.f; phiRatioHist->Fill(2.f/*/(float)recoSettings->topN*/); phiRatioDone = true; }
//
//
//   for(int i=0; i<settings->topN-2; i++){
//
//   vector<float> newThetaVec(thetaVec);
//   vector<float> newPhiVec(phiVec);
//
//   newThetaVec.erase(newThetaVec.begin()+thetaIdx[i]);
//   newPhiVec.erase(newPhiVec.begin()+phiIdx[i]);
//
//   grThetaRMS->SetPoint(i, i+1, getRMS(newThetaVec));
//   grPhiRMS->SetPoint(i, i+1, getRMS(newPhiVec));
//
//   vector<float>().swap(newThetaVec);
//   vector<float>().swap(newPhiVec);
//
//   }
//
//   for(int i=0; i<=settings->topN-2; i++){
//
//      growingThetaVec.push_back(thetaVec[thetaIdx[2+i]]);
//      growingPhiVec.push_back(phiVec[phiIdx[2+i]]);
//
//      rms = getRMS(growingThetaVec);
//      grGrowingThetaRMS->SetPoint(i+1, i+3, rms);
//      if(rms > rmsCut && thetaRatioDone == false){ thetaXing = i+3.f; thetaRatioHist->Fill((float)(i+3.f)/*/(float)recoSettings->topN*/); thetaRatioDone = true; }
//      rms = getRMS(growingPhiVec);
//      grGrowingPhiRMS->SetPoint(i+1, i+3, rms);
//      if(rms > rmsCut && phiRatioDone == false){ phiXing = i+3.f; phiRatioHist->Fill((float)(i+3.f)/*/(float)recoSettings->topN*/); phiRatioDone = true; }
//      //if(entry==0) cout<<getRMS(growingPhiVec)<<endl;
//   }
//
//   if(thetaRatioDone == false){ thetaXing = settings->topN+1; thetaRatioHist->Fill(settings->topN+1); }
//   if(phiRatioDone == false){   phiXing = settings->topN+1; phiRatioHist->Fill(settings->topN+1); }
//
//   //cvs = new TCanvas("cvs","cvs",800,800);
//   cvs->Divide(1,2);
//   if(entry==1){
//
//   cvs->cd(1);
//   grThetaRMS->Draw("AL");
//   cvs->cd(2);
//   grPhiRMS->Draw("AL");
//   cvs->SaveAs("grThetaPhiRMS.C");
//
//   cvs->cd(1);
//   grGrowingThetaRMS->Draw("AL");
//   cvs->cd(2);
//   grGrowingPhiRMS->Draw("AL");
//   cvs->SaveAs("grGrowingThetaPhiRMS.C");
//
//   }
}//end of RF
/*
   trueVertexGraph->SetPoint(0, dummyData->trueRadius*sin(dummyData->trueZen*TMath::DegToRad())*cos(dummyData->trueAzi*TMath::DegToRad()),
                                dummyData->trueRadius*sin(dummyData->trueZen*TMath::DegToRad())*sin(dummyData->trueAzi*TMath::DegToRad()),
                                dummyData->trueRadius*cos(dummyData->trueZen*TMath::DegToRad()));
   recoVertexGraph->SetPoint(0, dummyData->recoRadius*sin(dummyData->recoZen*TMath::DegToRad())*cos(dummyData->recoAzi*TMath::DegToRad()),
                                dummyData->recoRadius*sin(dummyData->recoZen*TMath::DegToRad())*sin(dummyData->recoAzi*TMath::DegToRad()),
                                dummyData->recoRadius*cos(dummyData->recoZen*TMath::DegToRad()));
   trueVertexGraph->SetMarkerColor(kRed);
   recoVertexGraph->SetMarkerColor(kOrange);
   topMaxPixGraph->SetMarkerColor(kGray);
*/
/*
TCanvas *cvs2 = new TCanvas("cvs2","cvs2",800,800);
if(entry==0){
   cvs2->cd();
   hist->Draw();
   hist->SetStats(0);
   topMaxPixGraph->ls();
   topMaxPixGraph->Draw("same");
   trueVertexGraph->Draw("same");
   recoVertexGraph->Draw("same");
   if(dummyData->trueRadius > 1000)
   snprintf(filename, 200*sizeof(char), "topMaxPixGraph_E18_far_%d.C", entry);
   else
   snprintf(filename, 200*sizeof(char), "topMaxPixGraph_E18_near_%d.C", entry);
   cvs2->SaveAs(filename);
}
*/

//   snprintf(filename,sizeof(char)*200,"maxPixEachLayerZen_%d",entry);
//   maxPixEachLayerZen=new TH1F(filename,filename,450,0,180);
//   snprintf(filename,sizeof(char)*200,"maxPixEachLayerAzi_%d",entry);
//   maxPixEachLayerAzi=new TH1F(filename,filename,900,0,360);
//
//   if(entry%1000==0) cout<<"Entry: "<<entry<<endl;
//
////if(dummyData->maxPixCoherence>0.12){
//   for(int i=0; i<onion.nLayer; i++){
//
//   //cout<<"maxPixIdxEachLayer size: "<<maxPixIdxEachLayer->size()<<endl;
//   r     = onion.getLayerRadius(dummyData->maxPixIdxEachLayer.at(i));
//   theta = onion.getPointing(dummyData->maxPixIdxEachLayer.at(i)).theta * TMath::RadToDeg();
//   phi   = onion.getPointing(dummyData->maxPixIdxEachLayer.at(i)).phi   * TMath::RadToDeg();
//
//   maxPixEachLayerZen->Fill(theta);
//   maxPixEachLayerAzi->Fill(phi);
/*
   if( i%2 == 0 ){
      coherence25[ i/2 ] = dummyData->maxPixCoherenceEachLayer.at(i);
      pixIdx25[ i/2 ] = dummyData->maxPixIdxEachLayer.at(i);
   }

  if( i%5 == 0 ){
      coherence10[ i/5 ] = dummyData->maxPixCoherenceEachLayer.at(i);
      pixIdx10[ i/5 ] = dummyData->maxPixIdxEachLayer.at(i);

   }

   if( i%10 == 0 ){
      coherence5[ i/10 ] = dummyData->maxPixCoherenceEachLayer.at(i);
      pixIdx5[ i/10 ] = dummyData->maxPixIdxEachLayer.at(i);

   }
*/
//   if( i == 0 ){
//      coherence2[ 0 ] = dummyData->maxPixCoherenceEachLayer.at(i);
//      pixIdx2[ 0 ] = dummyData->maxPixIdxEachLayer.at(i);
//   }
//
//   if( i == 30 ){
//      coherence2[ 1 ] = dummyData->maxPixCoherenceEachLayer.at(i);
//      pixIdx2[ 1 ] = dummyData->maxPixIdxEachLayer.at(i);
//   }
//
//
//   //printf("r:%f theta:%f phi:%f\n",r,theta,phi);
//
//   //maxPixEachLayerGraph->SetPoint(i, r*sin(theta)*cos(phi), r*sin(theta)*sin(phi), r*cos(theta));
//
//   //printf("%f\t", maxPixCoherenceEachLayer->at(i));
//   }//end of nLayer
//
//   maxPixEachLayerZenSpread->Fill(maxPixEachLayerZen->GetRMS());
//   maxPixEachLayerAziSpread->Fill(maxPixEachLayerAzi->GetRMS());
//
//   maxPixEachLayerZenAziSpread->Fill(maxPixEachLayerAzi->GetRMS(), maxPixEachLayerZen->GetRMS());
//
//   //TMath::Sort(25, coherence25, idx25);
//   //TMath::Sort(10, coherence10, idx10);
//   //TMath::Sort(5,  coherence5,  idx5);
//   TMath::Sort(2,  coherence2,  idx2);
//
//   //coherence25Hist->Fill(coherence25[idx25[0]]);
//   //coherence10Hist->Fill(coherence10[idx10[0]]);
//   //coherence5Hist->Fill(coherence5[idx5[0]]);
//   coherence2Hist->Fill(coherence2[idx2[0]]);
//
//if(dummyData->eventTrigType == 0){
//
//coherence2_rf->Fill(coherence2[idx2[0]]);
//coherence2_snr_rf->Fill(dummyData->inWindowSNR_V, coherence2[idx2[0]]/*, dummyData->weight*/);
///*
//int nGoodChan = 0;
//for(int i=0; i<8; i++) nGoodChan += dummyData->recoChan[i];
//
////cout<<"nGoodChan: "<<nGoodChan<<endl;
//float normCoherence2 = coherence2[idx2[0]] * 28.f / (float)((nGoodChan*(nGoodChan-1))/2);
//coherence2Hist_nGoodChanNorm->Fill(normCoherence2);
//*/
////if(normCoherence2 > 0.24){ zenDist->Fill(dummyData->recoZen, normCoherence2); aziDist->Fill(dummyData->recoAzi, normCoherence2); }
////if(coherence2[idx2[0]] > 0.24){ zenDist->Fill(dummyData->recoZen); aziDist->Fill(dummyData->recoAzi); zenDist_coherence2->Fill(coherence2[idx2[0]],dummyData->recoZen); aziDist_coherence2->Fill(coherence2[idx2[0]],dummyData->recoAzi); }
//
//} else if( dummyData->eventTrigType == 1) { coherence2_cal->Fill(coherence2[idx2[0]]);  coherence2_snr_cal->Fill(dummyData->inWindowSNR_V, coherence2[idx2[0]]); }
//  else if( dummyData->eventTrigType == 2) { coherence2_soft->Fill(coherence2[idx2[0]]); coherence2_snr_soft->Fill(dummyData->inWindowSNR_V, coherence2[idx2[0]]); }
//
//
// //  dCoherence25Hist->Fill(coherence25[idx25[0]]-dummyData->maxPixCoherence);
// //  dCoherence10Hist->Fill(coherence10[idx10[0]]-dummyData->maxPixCoherence);
// //  dCoherence5Hist->Fill(coherence5[idx5[0]]-dummyData->maxPixCoherence);
//   dCoherence2Hist->Fill(coherence2[idx2[0]]-dummyData->maxPixCoherence);
//
//   coherence2_maxPixEachLayerZenSpread->Fill(maxPixEachLayerZen->GetRMS(), coherence2[idx2[0]]);
//   coherence2_maxPixEachLayerAziSpread->Fill(maxPixEachLayerAzi->GetRMS(), coherence2[idx2[0]]);
//   coherence2_topMaxPixZenSpread->Fill(topMaxPixZen->GetRMS(), coherence2[idx2[0]]);
//   coherence2_topMaxPixAziSpread->Fill(topMaxPixAzi->GetRMS(), coherence2[idx2[0]]);
//
//if(dummyData->eventTrigType == 0){
//   coherence2_thetaRMSCutXing->Fill(thetaXing, coherence2[idx2[0]]);
//   coherence2_phiRMSCutXing->Fill(phiXing, coherence2[idx2[0]]);
//   topMaxPixZenSpread_thetaRMSCutXing->Fill(thetaXing, topMaxPixZen->GetRMS());
//   topMaxPixAziSpread_phiRMSCutXing->Fill(phiXing, topMaxPixAzi->GetRMS());
//}
/*
   //cout<<"idx25[0]: "<<idx25[0]<<" pixIdx25[idx25[0]]: "<<pixIdx25[idx25[0]]<<endl;
   r     = onion.getLayerRadius(pixIdx25[idx25[0]]);
   theta = TMath::RadToDeg()*onion.getPointing(pixIdx25[idx25[0]]).theta;
   phi   = TMath::RadToDeg()*onion.getPointing(pixIdx25[idx25[0]]).phi;
   dZen_weight_nLayer25->Fill(theta-dummyData->trueZen,    dummyData->weight);
   dAzi_weight_nLayer25->Fill(phi  -dummyData->trueAzi,    dummyData->weight);
   dR_weight_nLayer25->Fill(  r    -dummyData->trueRadius, dummyData->weight);
   dRecoZen_weight_nLayer25->Fill(theta-dummyData->recoZen,    dummyData->weight);
   dRecoAzi_weight_nLayer25->Fill(phi  -dummyData->recoAzi,    dummyData->weight);
   dRecoR_weight_nLayer25->Fill(  r    -dummyData->recoRadius, dummyData->weight);
   dRecoZen_dCoherence_weight_nLayer25->Fill(theta-dummyData->recoZen, coherence25[idx25[0]]-dummyData->maxPixCoherence, dummyData->weight);
   dRecoAzi_dCoherence_weight_nLayer25->Fill(phi  -dummyData->recoAzi, coherence25[idx25[0]]-dummyData->maxPixCoherence, dummyData->weight);
   dRecoR_dCoherence_weight_nLayer25->Fill(  r    -dummyData->recoRadius, coherence25[idx25[0]]-dummyData->maxPixCoherence, dummyData->weight);


   //cout<<"idx10[0]: "<<idx10[0]<<" pixIdx10[idx10[0]]: "<<pixIdx10[idx10[0]]<<endl;
   r     = onion.getLayerRadius(pixIdx10[idx10[0]]);
   theta = TMath::RadToDeg()*onion.getPointing(pixIdx10[idx10[0]]).theta;
   phi   = TMath::RadToDeg()*onion.getPointing(pixIdx10[idx10[0]]).phi;
   dZen_weight_nLayer10->Fill(theta-dummyData->trueZen,    dummyData->weight);
   dAzi_weight_nLayer10->Fill(phi  -dummyData->trueAzi,    dummyData->weight);
   dR_weight_nLayer10->Fill(  r    -dummyData->trueRadius, dummyData->weight);
   dRecoZen_weight_nLayer10->Fill(theta-dummyData->recoZen,    dummyData->weight);
   dRecoAzi_weight_nLayer10->Fill(phi  -dummyData->recoAzi,    dummyData->weight);
   dRecoR_weight_nLayer10->Fill(  r    -dummyData->recoRadius, dummyData->weight);
    dRecoZen_dCoherence_weight_nLayer10->Fill(theta-dummyData->recoZen, coherence10[idx10[0]]-dummyData->maxPixCoherence, dummyData->weight);
   dRecoAzi_dCoherence_weight_nLayer10->Fill(phi  -dummyData->recoAzi, coherence10[idx10[0]]-dummyData->maxPixCoherence, dummyData->weight);
   dRecoR_dCoherence_weight_nLayer10->Fill(  r    -dummyData->recoRadius, coherence10[idx10[0]]-dummyData->maxPixCoherence, dummyData->weight);


   //cout<<"idx5[0]: "<<idx5[0]<<" pixIdx5[idx5[0]]: "<<pixIdx5[idx5[0]]<<endl;
   r     = onion.getLayerRadius(pixIdx5[idx5[0]]);
   theta = TMath::RadToDeg()*onion.getPointing(pixIdx5[idx5[0]]).theta;
   phi   = TMath::RadToDeg()*onion.getPointing(pixIdx5[idx5[0]]).phi;
   dZen_weight_nLayer5->Fill(theta-dummyData->trueZen,    dummyData->weight);
   dAzi_weight_nLayer5->Fill(phi  -dummyData->trueAzi,    dummyData->weight);
   dR_weight_nLayer5->Fill(  r    -dummyData->trueRadius, dummyData->weight);
   dRecoZen_weight_nLayer5->Fill(theta-dummyData->recoZen,    dummyData->weight);
   dRecoAzi_weight_nLayer5->Fill(phi  -dummyData->recoAzi,    dummyData->weight);
   dRecoR_weight_nLayer5->Fill(  r    -dummyData->recoRadius, dummyData->weight);
    dRecoZen_dCoherence_weight_nLayer5->Fill(theta-dummyData->recoZen, coherence5[idx5[0]]-dummyData->maxPixCoherence, dummyData->weight);
   dRecoAzi_dCoherence_weight_nLayer5->Fill(phi  -dummyData->recoAzi, coherence5[idx5[0]]-dummyData->maxPixCoherence, dummyData->weight);
   dRecoR_dCoherence_weight_nLayer5->Fill(  r    -dummyData->recoRadius, coherence5[idx5[0]]-dummyData->maxPixCoherence, dummyData->weight);
*/

   //cout<<"idx2[0]: "<<idx2[0]<<" pixIdx2[idx2[0]]: "<<pixIdx2[idx2[0]]<<endl;
//   r     = onion.getLayerRadius(pixIdx2[idx2[0]]);
//   theta = TMath::RadToDeg()*onion.getPointing(pixIdx2[idx2[0]]).theta;
//   phi   = TMath::RadToDeg()*onion.getPointing(pixIdx2[idx2[0]]).phi;
//   dZen_weight_nLayer2->Fill(theta-dummyData->trueZen,    dummyData->weight);
//   dAzi_weight_nLayer2->Fill(phi  -dummyData->trueAzi,    dummyData->weight);
//   dR_weight_nLayer2->Fill(  r    -dummyData->trueRadius, dummyData->weight);
//   dRecoZen_weight_nLayer2->Fill(theta-dummyData->recoZen,    dummyData->weight);
//   dRecoAzi_weight_nLayer2->Fill(phi  -dummyData->recoAzi,    dummyData->weight);
//   dRecoR_weight_nLayer2->Fill(  r    -dummyData->recoRadius, dummyData->weight);
//    dRecoZen_dCoherence_weight_nLayer2->Fill(theta-dummyData->recoZen, coherence2[idx2[0]]-dummyData->maxPixCoherence, dummyData->weight);
//   dRecoAzi_dCoherence_weight_nLayer2->Fill(phi  -dummyData->recoAzi, coherence2[idx2[0]]-dummyData->maxPixCoherence, dummyData->weight);
//   dRecoR_dCoherence_weight_nLayer2->Fill(  r    -dummyData->recoRadius, coherence2[idx2[0]]-dummyData->maxPixCoherence, dummyData->weight);
//
//
//   cout<<"maxPixIdx: "<<dummyData->maxPixIdx<<" maxPixCoherence: "<<dummyData->maxPixCoherence<<endl;

//}
   //cout<<endl;
   //maxPixEachLayerGraph->RemovePoint(0);
/*
   cvs.cd();
   hist->Draw();
   hist->SetStats(0);
   maxPixEachLayerGraph->Draw("same");
   snprintf(filename, sizeof(filename), "maxPixEachLayerGraph_1_%d.C", entry);
   cvs.SaveAs(filename);
*/
/*
  alpha = getSpaceAngle(dummyData->trueZen*TMath::DegToRad(), dummyData->trueAzi*TMath::DegToRad(),
                         onion.getPointing(dummyData->maxPixIdx).theta, onion.getPointing(dummyData->maxPixIdx).phi
                        ) * TMath::RadToDeg();

   dAlpha_coherence->Fill(dummyData->maxPixCoherence, alpha);
   dAlpha_nchnl->Fill(nchnl, alpha);
   dAlpha_trueR->Fill(dummyData->trueRadius, alpha);
   dAlpha_recoZen->Fill(dummyData->recoZen, alpha);
   dAlpha_trueZen->Fill(dummyData->trueZen, alpha);
   dAlpha_recoAzi->Fill(dummyData->recoAzi, alpha);
   dAlpha_trueAzi->Fill(dummyData->trueAzi, alpha);


   dAlphaRProf->Fill(dummyData->trueRadius, alpha);
   dZenRProf->Fill(dummyData->trueRadius,  dummyData->recoZen - dummyData->trueZen);
   dAziRProf->Fill(dummyData->trueRadius,  dummyData->recoAzi - dummyData->trueAzi);
   //cout<<"recoRadius: "<<dummyData->recoRadius<<" maxPix layer radius: "<<onion.getLayerRadius(dummyData->maxPixIdx)<<endl;
   dRRProf->Fill(dummyData->trueRadius, onion.getLayerRadius(dummyData->maxPixIdx) - dummyData->trueRadius);

   dZen_dR->Fill(onion.getLayerRadius(dummyData->maxPixIdx) - dummyData->trueRadius, dummyData->recoZen - dummyData->trueZen);
   dAzi_dR->Fill(onion.getLayerRadius(dummyData->maxPixIdx) - dummyData->trueRadius, dummyData->recoAzi - dummyData->trueAzi);
   dAlpha_dR->Fill(onion.getLayerRadius(dummyData->maxPixIdx) - dummyData->trueRadius, alpha);
*/
/*
   printf("Entry: %d\nTrue Radius: %f Reco Radius: %f True Zen: %f Reco Zen: %f True Azi: %f Reco Azi: %f\n",
          entry, dummyData->trueRadius, dummyData->recoRadius,
          dummyData->trueZen, dummyData->recoZen, dummyData->trueAzi, dummyData->recoAzi);
*/
/*
   dZen_coherence->Fill(dummyData->maxPixCoherence, dummyData->recoZen - dummyData->trueZen, dummyData->weight);
   dZen_nchnl->Fill(nchnl, dummyData->recoZen - dummyData->trueZen, dummyData->weight);
   dZen_r->Fill(dummyData->trueRadius, dummyData->recoZen - dummyData->trueZen, dummyData->weight);
   dZen_recoZen->Fill(dummyData->recoZen, dummyData->recoZen - dummyData->trueZen, dummyData->weight);
   dZen_trueZen->Fill(dummyData->trueZen, dummyData->recoZen - dummyData->trueZen, dummyData->weight);
   dZen_recoAzi->Fill(dummyData->recoAzi, dummyData->recoZen - dummyData->trueZen, dummyData->weight);
   dZen_trueAzi->Fill(dummyData->trueAzi, dummyData->recoZen - dummyData->trueZen, dummyData->weight);

   dAzi_coherence->Fill(dummyData->maxPixCoherence, dummyData->recoAzi - dummyData->trueAzi, dummyData->weight);
   dAzi_nchnl->Fill(nchnl, dummyData->recoAzi - dummyData->trueAzi, dummyData->weight);
   dAzi_r->Fill(dummyData->trueRadius, dummyData->recoAzi - dummyData->trueAzi, dummyData->weight);
   dAzi_recoZen->Fill(dummyData->recoZen, dummyData->recoAzi - dummyData->trueAzi, dummyData->weight);
   dAzi_trueZen->Fill(dummyData->trueZen, dummyData->recoAzi - dummyData->trueAzi, dummyData->weight);
   dAzi_recoAzi->Fill(dummyData->recoAzi, dummyData->recoAzi - dummyData->trueAzi, dummyData->weight);
   dAzi_trueAzi->Fill(dummyData->trueAzi, dummyData->recoAzi - dummyData->trueAzi, dummyData->weight);

   dZen_dAzi->Fill(dummyData->recoAzi - dummyData->trueAzi, dummyData->recoZen - dummyData->trueZen, dummyData->weight);
*/
/*
   if(fabs(dummyData->recoZen - dummyData->trueZen) > 3 ){

   thetaCnt++;
   if( fabs(dummyData->recoAzi - dummyData->trueAzi) > 1 ) thetaPhiCnt++;

   parCnt++;
   trueR_par->Fill(dummyData->trueRadius);
   w_par->Fill(dummyData->weight);
   coherence_par->Fill(dummyData->maxPixCoherence);
   nchnl=0;
   for(int ch=0; ch<8; ch++){
   nchnl+=dummyData->recoChan[ch];
   if(dummyData->recoChan[ch]) usedChan_par->Fill(ch);
   }
   nchnl_par->Fill(nchnl);

   }

   if( fabs(dummyData->recoAzi - dummyData->trueAzi) > 1 ) phiCnt++;
*/
/*
 * Trying to extract more precise radius reconstruction using statistical tools
 */
/*float **///value = (float*)calloc(onion.nLayer, sizeof(float));
/*float **///radii = (float*)calloc(onion.nLayer, sizeof(float));
/*float **///pixCount = (float*)calloc(onion.nLayer, sizeof(float));
float M;
/*
for(int pix=0; pix<dummyData->topN; pix++){

   // 1. Sum of M values on each layer
   value[onion.getLayerNumber(dummyData->topMaxPixIdx.at(pix))] += dummyData->topMaxPixCoherence.at(pix);
   //cout<<"topMaxPixCoherence: "<<topMaxPixCoherence->at(pix)<<endl;
   pixCount[onion.getLayerNumber(dummyData->topMaxPixIdx.at(pix))] += 1.f;

}

for(int layer=0; layer<onion.nLayer; layer++){

   radii[layer] = onion.layerRadii[layer];
   // 2. Average M values on each layer
   if(pixCount[layer] != 0 )   value[layer] /= pixCount[layer];
   // 3. Simple max pix on each layer
   //value[layer] = maxPixCoherenceEachLayer->at(layer);

}

statRecoRadius =
statisticalRadiusReco(onion.nLayer, radii, value, dummyData->trueRadius, statRReco_weight, statRReco_max, statRReco_eventStack);
   //topMaxPixGraph = 0;

dStatR_coherence->Fill(dummyData->maxPixCoherence, statRecoRadius - dummyData->trueRadius);
dStatR_nchnl->Fill(nchnl, statRecoRadius - dummyData->trueRadius);
if( dummyData->maxPixCoherence > 5.f )
dStatR_trueR->Fill(dummyData->trueRadius, statRecoRadius - dummyData->trueRadius);
for(int j=0; j<5; j++){

   if(dummyData->trueRadius < rArray[j] ) dStatR[j]->Fill(statRecoRadius - dummyData->trueRadius);

}
*/
/*
cvs = new TCanvas("cvs","cvs",800,800);
if(plotCount<100){
cvs->Divide(2,1);
cvs->cd(1);
maxPixEachLayerZen->Draw();
cvs->cd(2);
maxPixEachLayerAzi->Draw();
snprintf(filename, sizeof(char)*200, "maxPixEachLayerZenAzi_%d.C",entry);
//cvs->SaveAs(filename);
plotCount++;
}
*/

//timeSequenceHist->Fill(dummyData->timeSequenceParameter, dummyData->weight);

//delete cvs;
//delete maxPixEachLayerZen, maxPixEachLayerAzi;
//delete topMaxPixZen, topMaxPixAzi;
//dummyData->clear();
//delete dummyData;

free(dRMSThetaVec);
free(dRMSPhiVec);
free(thetaIdx);
free(phiIdx);
delete grPhiRMS, grThetaRMS, grGrowingPhiRMS, grGrowingThetaRMS;
//free(value);
//free(radii);
//free(pixCount);
//cout<<"end of entry"<<endl;
//topMaxPixGraph->SetPolyMarker(0,&p,6);
delete topMaxPixGraph;
delete trueVertexGraph;
delete recoVertexGraph;
}//end of entry

cout<<"parCnt: "<<parCnt<<endl;
printf("thetaCnt: %d phiCnt: %d thetaPhiCnt: %d\n",thetaCnt,phiCnt,thetaPhiCnt);

//char histName[200];
TH1D *projY[10];

for(int xbin=1; xbin<=10; xbin++){

snprintf(histName,sizeof(histName),"hist_%d",xbin);
//projY = new TH1D(histName, histName, 500, -5000, 5000);
projY[xbin-1] = dR_trueR->ProjectionY(histName, xbin, xbin);
cout<<"xbin: "<<xbin<<" rms: "<<projY[xbin-1]->GetRMS()<<endl;
}

char legText[200];
TCanvas c1("c1","c1",800,600);
TCanvas c2("c2","c2",800,600);
TCanvas c3("c3","c3",800,600);
TCanvas c4("c4","c4",800,600);
TCanvas c5("c5","c5",800,600);

//c1.Divide(3,2);
c1.Divide(2,1);
c1.cd(1);
/*
r_all->Draw();
r_par->SetLineColor(kRed);
r_par->Draw("same");
*/
//dZen_coherence->Draw("colz");
//dZen_recoZen->Draw("colz");
//dAzi_coherence->Draw("colz");
//dAzi_recoZen->Draw("colz");
//dZen->Draw();
//dZen_cut->Draw();
//coherence_nchnl->Draw("colz");
//dZen_trueR->Draw("colz");
zen_all->SetTitle("Reconstructed zenith;[#circ];Entry");
zen_all->Draw();
zen_rf->SetLineColor(kRed);
zen_rf->Draw("same");
zen_cal->SetLineColor(kMagenta);
zen_cal->Draw("same");
zen_soft->SetLineColor(kOrange);
zen_soft->Draw("same");

TLine *zen_A2D5TH = new TLine(90-108.823,0,90-108.823,1e8);
TLine *zen_A2D5BV = new TLine(90-112.76,0,90-112.76,1e8);
TLine *zen_A2D6TH = new TLine(90-75.9206,0,90-75.9206,1e8);
TLine *zen_A2D6BV = new TLine(90-80.4426,0,90-80.44226,1e8);

zen_A2D5TH->SetLineWidth(2);
zen_A2D5BV->SetLineWidth(2);
zen_A2D6TH->SetLineWidth(2);
zen_A2D6BV->SetLineWidth(2);

zen_A2D5TH->SetLineStyle(7);
zen_A2D5BV->SetLineStyle(7);
zen_A2D6TH->SetLineStyle(7);
zen_A2D6BV->SetLineStyle(7);
/*
zen_A2D5TH->Draw("same");
zen_A2D5BV->Draw("same");
zen_A2D6TH->Draw("same");
zen_A2D6BV->Draw("same");
*/

TLine *zen_IC1Shallow = new TLine(-22.7,0,-22.7,500);
TLine *zen_IC22Shallow = new TLine(-23.1,0,-23.1,500);

zen_IC1Shallow->SetLineWidth(2);
zen_IC22Shallow->SetLineWidth(2);

zen_IC1Shallow->SetLineStyle(7);
zen_IC22Shallow->SetLineStyle(7);

zen_IC1Shallow->Draw("same");
zen_IC22Shallow->Draw("same");

c1.cd(2);
/*
w_all->Draw();
w_par->SetLineColor(kRed);
w_par->Draw("same");
*/
//coherence_all->Draw();
/*
coherence_par->SetLineColor(kRed);
coherence_par->Draw("same");
*/
//dZen_nchnl->Draw("colz");
//dZen_trueZen->Draw("colz");
//dAzi_nchnl->Draw("colz");
//dAzi_trueZen->Draw("colz");
//dAzi->Draw();
//dAzi_cut->Draw();
//coherence_r->Draw("colz");
//dAzi_trueR->Draw("colz");
azi_all->SetTitle("Reconstructed azimuth;[#circ];Entry");
azi_all->Draw();
azi_rf->SetLineColor(kRed);
azi_rf->Draw("same");
azi_cal->SetLineColor(kMagenta);
azi_cal->Draw("same");
azi_soft->SetLineColor(kOrange);
azi_soft->Draw("same");

TLine *azi_A2D6 = new TLine(63.3316,0,63.3316,1e7);
//TLine *azi_A2D5BV = new TLine(63.3316,0,1e7);
TLine *azi_A2D5 = new TLine(334.855,0,334.855,1e7);

azi_A2D6->SetLineWidth(2);
azi_A2D5->SetLineWidth(2);

azi_A2D6->SetLineStyle(7);
azi_A2D5->SetLineStyle(7);
/*
azi_A2D6->Draw("same");
azi_A2D5->Draw("same");
*/

TLine *azi_IC1Shallow = new TLine(260.0,0,260.0,500);
TLine *azi_IC22Shallow = new TLine(266.0,0,266.0,500);

azi_IC1Shallow->SetLineWidth(2);
azi_IC22Shallow->SetLineWidth(2);

azi_IC1Shallow->SetLineStyle(7);
azi_IC22Shallow->SetLineStyle(7);

azi_IC1Shallow->Draw("same");
azi_IC22Shallow->Draw("same");

//c1.cd(3);
//dR->Draw();
//dR_cut->Draw();
//dAlpha_trueR->Draw("colz");
//dAzi_coherence->Draw("colz");
//c1.cd(4);
/*
nchnl_all->Draw();
nchnl_par->SetLineColor(kRed);
nchnl_par->Draw("same");
*/
//dZen_r->Draw("colz");
//dZen_recoAzi->Draw("colz");
//dAzi_r->Draw("colz");
//dAzi_recoAzi->Draw("colz");
//dZen_weight->Draw();
//dZen_weight_cut->Draw();
//dAzi_nchnl->Draw("colz");
//nchnl_r->Draw("colz");
//dR_trueR->Draw("colz");
//c1.cd(5);
/*
usedChan_all->Draw();
usedChan_par->SetLineColor(kRed);
usedChan_par->Draw("same");
*/
//dZen_dAzi->Draw("colz");
//dZen_trueAzi->Draw("colz");
//dAzi_dZen->Draw("colz");
//dAzi_trueAzi->Draw("colz");
//dAzi_weight->Draw();
//dAzi_weight_cut->Draw();
//c1.cd(6);
//dR_weight->Draw();
//dR_weight_cut->Draw();

//c1.SaveAs("recoAnalysis_1.C");

c2.cd();
//c2.Divide(1,3);
//c2.cd(1);
//dR_coherence->Draw("colz");
//dStatR_coherence->Draw("colz");
//dStatR_trueR->Draw("colz");
zen_azi_rf->SetTitle("Reconstructed RF events;Azimuth [#circ]; Zenith [#circ]");
zen_azi_rf->Draw("colz");
//zen_azi_soft->Draw("colzsame");

float IC1DeepPulserPosZen[1] = {-22.7};
float IC22DeepPulserPosZen[1] = {-23.1};
float IC1DeepPulserPosAzi[1] = {260.0};
float IC22DeepPulserPosAzi[1] = {266.0};

TGraph *deepPulserPosition_IC1 = new TGraph(1, IC1DeepPulserPosAzi, IC1DeepPulserPosZen);
TGraph *deepPulserPosition_IC22 = new TGraph(1, IC22DeepPulserPosAzi, IC22DeepPulserPosZen);

deepPulserPosition_IC1->SetMarkerStyle(22);
deepPulserPosition_IC1->SetMarkerSize(0.8);
deepPulserPosition_IC22->SetMarkerStyle(23);
deepPulserPosition_IC22->SetMarkerSize(0.8);

deepPulserPosition_IC1->Draw("psame");
deepPulserPosition_IC22->Draw("psame");
/*
c2.cd(2);
//dR_nchnl->Draw("colz");
//dStatR_nchnl->Draw("colz");
zen_azi_cal->SetTitle("Reconstructed calpulser events;Azimuth [#circ]; Zenith [#circ]");
zen_azi_cal->Draw("colz");

float A2D5PulserPosZen[2] = {90-108.823, 90-112.76};
float A2D5PulserPosAzi[2] = {334.855, 334.855};
float A2D6PulserPosZen[2] = {90-75.9206, 90-80.4426};
float A2D6PulserPosAzi[2] = {63.3316, 63.3316};

TGraph *calpulserPosition_A2D5 = new TGraph(2,A2D5PulserPosAzi,A2D5PulserPosZen);
TGraph *calpulserPosition_A2D6 = new TGraph(2,A2D6PulserPosAzi,A2D6PulserPosZen);

calpulserPosition_A2D5->SetMarkerStyle(22);
calpulserPosition_A2D5->SetMarkerSize(0.8);
calpulserPosition_A2D6->SetMarkerStyle(23);
calpulserPosition_A2D6->SetMarkerSize(0.8);
/*
calpulserPosition_A2D6->Draw("psame");
calpulserPosition_A2D5->Draw("psame");
*/
/*
c2.cd(3);
//dStatR_trueR->Draw("colz");
zen_azi_soft->SetTitle("Reconstructed software trigger events;Azimuth [#circ]; Zenith [#circ]");
zen_azi_soft->Draw("colz");
/*
hist->Draw();
hist->SetStats(0);
//gStyle->SetOptStat(0);
topMaxPixGraph->Draw("same");
*/
//topMaxPixGraph->Draw("same ogl");
//c2.cd(4);
//dStatR[0]->SetLineColor(kBlue);
//dStatR[0]->Draw();
//for(int j=1; j<5; j++){

//   dStatR[j]->SetLineColor(kBlue+j);
//   dStatR[j]->Draw("same");
//}


//c2.SaveAs("recoAnalysis_2.C");


c3.cd();
c3.Divide(3,2);
c3.cd(1);
//dAlpha_coherence->Draw("colz");
zen_coherence_rf->SetTitle("Reconstructed zenith vs coherence (RF);Max skymap coherence [arb. unit];[#circ]");
zen_coherence_rf->Draw("colz");

c3.cd(2);
//dAlpha_nchnl->Draw("colz");
zen_coherence_cal->SetTitle("Reconstructed zenith vs coherence (Cal);Max skymap coherence [arb. unit];[#circ]");
zen_coherence_cal->Draw("colz");

c3.cd(3);
zen_coherence_soft->SetTitle("Reconstructed zenith vs coherence (Soft);Max skymap coherence [arb. unit];[#circ]");
zen_coherence_soft->Draw("colz");

c3.cd(4);
azi_coherence_rf->SetTitle("Reconstructed azimuth vs coherence (RF);Max skymap coherence [arb. unit];[#circ]");
azi_coherence_rf->Draw("colz");

c3.cd(5);
azi_coherence_cal->SetTitle("Reconstructed azimuth vs coherence (Cal);Max skymap coherence [arb. unit];[#circ]");
azi_coherence_cal->Draw("colz");

c3.cd(6);
azi_coherence_soft->SetTitle("Reconstructed azimuth vs coherence (Soft);Max skymap coherence [arb. unit];[#circ]");
azi_coherence_soft->Draw("colz");
/*
hist->Draw();
hist->SetStats(0);
//gStyle->SetOptStat(0);
//hist->Draw("A");
maxPixEachLayerGraph->Draw("same");
//gStyle->SetCanvasPreferGL(1);
*/
//topMaxPixHist->Draw("colz");
//trueDistanceLine.Draw("same");
//c3.SaveAs("recoAnalysis_3.C");

c3.cd(1);
zen_snr_rf->SetTitle("Reconstructed zenith vs SNR (RF);Signal-to-noise ratio [arb. unit];[#circ]");
zen_snr_rf->Draw("colz");

c3.cd(2);
zen_snr_cal->SetTitle("Reconstructed zenith vs SNR (Cal);Signal-to-noise ratio [arb. unit];[#circ]");
zen_snr_cal->Draw("colz");

c3.cd(3);
zen_snr_soft->SetTitle("Reconstructed zenith vs SNR (Soft);Signal-to-noise ratio [arb. unit];[#circ]");
zen_snr_soft->Draw("colz");

c3.cd(4);
azi_snr_rf->SetTitle("Reconstructed azimuth vs SNR (RF);Signal-to-noise ratio [arb. unit];[#circ]");
azi_snr_rf->Draw("colz");

c3.cd(5);
azi_snr_cal->SetTitle("Reconstructed azimuth vs SNR (Cal);Signal-to-noise ratio [arb. unit];[#circ]");
azi_snr_cal->Draw("colz");

c3.cd(6);
azi_snr_soft->SetTitle("Reconstructed azimuth vs SNR (Soft);Signal-to-noise ratio [arb. unit];[#circ]");
azi_snr_soft->Draw("colz");

//c3.SaveAs("recoAnalysis_3.2.C");


//c4.Divide(2,2);
//c4.Divide(1,3);
c4.cd();
//gStyle->SetOptStat(0);
//c4.cd(1);
//dAlphaRProf->Draw();
//dZen_dR->Draw("colz");
//statRReco_max->Draw();
//dR->Draw();
coherence_all->SetTitle("Max skymap coherence;Coherence [arb. unit];Entry");
coherence_all->SetLineWidth(2);
coherence_all->Draw();

//c4.cd(2);
//dZenRProf->Draw();
//dAzi_dR->Draw("colz");
//statRReco_weight->Draw();
//dR_trueR->Draw("colz");

coherence_rf->SetLineWidth(2);
coherence_rf->SetLineColor(kRed);
coherence_rf->Draw("same");

//c4.cd(3);
//dAziRProf->Draw();
//dAlpha_dR->Draw("colz");
//statRReco_eventStack->Draw();
//dStatR_trueR->Draw("colz");
//dR_cut->Draw();

coherence_cal->SetLineWidth(2);
coherence_cal->SetLineColor(kMagenta);
coherence_cal->Draw("same");

//c4.cd(4);
//dRRProf->Draw();

coherence_soft->SetLineWidth(2);
coherence_soft->SetLineColor(kOrange);
coherence_soft->Draw("same");

//gStyle->SetOptStat(0);

TLegend *leg4 = new TLegend(0.1,0.7,0.48,0.9);
leg4->SetHeader("Max Skymap Coherence"/*,"C"*/); // option "C" allows to center the header
leg4->AddEntry(coherence_all,"All events","l");
leg4->AddEntry(coherence_rf,"RF events","l");
leg4->AddEntry(coherence_cal,"Cal events","l");
leg4->AddEntry(coherence_soft,"Soft events","l");
/*
snprintf(legText,sizeof(legText),"RF overflow bin entry: %d",coherence_rf->GetBinContent(coherence_rf->GetNbinsX()+1));
leg4->AddEntry((TObject*)0, legText, "");
snprintf(legText,sizeof(legText),"Cal overflow bin entry: %d",coherence_cal->GetBinContent(coherence_cal->GetNbinsX()+1));
leg4->AddEntry((TObject*)0, legText, "");
snprintf(legText,sizeof(legText),"Soft overflow bin entry: %d",coherence_soft->GetBinContent(coherence_soft->GetNbinsX()+1));
leg4->AddEntry((TObject*)0, legText, "");
*/
//leg4->Draw();

//c4.SaveAs("recoAnalysis_4.C");

c5.cd();
//c5.Divide(2,2);

//c5.Divide(2,2);
//c5.cd(1);
coherence_snr_all->SetTitle("Max skymap coherence vs SNR;Signal-to-noise ratio [arb. unit];Coherence [arb. unit]");
coherence_snr_all->Draw("colz");
/*
c5.cd(2);
coherence_snr_rf->SetTitle("Max skymap coherence vs SNR (RF);Signal-to-noise ratio [arb. unit];Coherence [arb. unit]");
coherence_snr_rf->Draw("colz");

c5.cd(3);
coherence_snr_cal->SetTitle("Max skymap coherence vs SNR (Cal);Signal-to-noise ratio [arb. unit];Coherence [arb. unit]");
coherence_snr_cal->Draw("colz");

c5.cd(4);
coherence_snr_soft->SetTitle("Max skymap coherence vs SNR (Soft);Signal-to-noise ratio [arb. unit];Coherence [arb. unit]");
coherence_snr_soft->Draw("colz");
*/
/*
c5.cd(1);
dR_trueR->Draw("colz");

c5.cd(2);
dR_coherence->Draw("colz");

c5.cd(3);
dZen_coherence->Draw("colz");

c5.cd(4);
dAzi_coherence->Draw("colz");
*/
//c5.SaveAs("recoAnalysis_5.C");

//TCanvas c6("c6","c6",800,600);
//c6.Divide(5,2);
//
//for(int i=1; i<=10; i++){
//   c6.cd(i);
//   projY[i-1]->Draw();
//}
//
////dRRatio_trueR->Draw("colz");
//
//c6.SaveAs("recoAnalysis_6.C");
//
//TCanvas c7("c7","c7",800,600);
//c7.Divide(3,1);
//c7.cd(1);
//dZen_weight->SetLineColor(kBlack);
//dZen_weight_nLayer25->SetLineColor(kRed);
//dZen_weight_nLayer10->SetLineColor(kMagenta+4);
//dZen_weight_nLayer5->SetLineColor(kMagenta+2);
//dZen_weight_nLayer2->SetLineColor(kMagenta);
//dZen_weight->Draw();
//dZen_weight_nLayer25->Draw("same");
//dZen_weight_nLayer10->Draw("same");
//dZen_weight_nLayer5->Draw("same");
//dZen_weight_nLayer2->Draw("same");
//
//printf("dZen (mean, rms) 50:(%f, %f)\t25:(%f, %f)\t10:(%f, %f)\t5:(%f, %f)\t2:(%f, %f)\n"
//      ,dZen_weight->GetMean(), dZen_weight->GetRMS()
//      ,dZen_weight_nLayer25->GetMean(), dZen_weight_nLayer25->GetRMS()
//      ,dZen_weight_nLayer10->GetMean(), dZen_weight_nLayer10->GetRMS()
//      ,dZen_weight_nLayer5->GetMean(), dZen_weight_nLayer5->GetRMS()
//      ,dZen_weight_nLayer2->GetMean(), dZen_weight_nLayer2->GetRMS()
//      );
//
//c7.cd(2);
//dAzi_weight->SetLineColor(kBlack);
//dAzi_weight_nLayer25->SetLineColor(kRed);
//dAzi_weight_nLayer10->SetLineColor(kMagenta+4);
//dAzi_weight_nLayer5->SetLineColor(kMagenta+2);
//dAzi_weight_nLayer2->SetLineColor(kMagenta);
//dAzi_weight->Draw();
//dAzi_weight_nLayer25->Draw("same");
//dAzi_weight_nLayer10->Draw("same");
//dAzi_weight_nLayer5->Draw("same");
//dAzi_weight_nLayer2->Draw("same");
//
//printf("dAzi (mean, rms) 50:(%f, %f)\t25:(%f, %f)\t10:(%f, %f)\t5:(%f, %f)\t2:(%f, %f)\n"
//      ,dAzi_weight->GetMean(), dAzi_weight->GetRMS()
//      ,dAzi_weight_nLayer25->GetMean(), dAzi_weight_nLayer25->GetRMS()
//      ,dAzi_weight_nLayer10->GetMean(), dAzi_weight_nLayer10->GetRMS()
//      ,dAzi_weight_nLayer5->GetMean(), dAzi_weight_nLayer5->GetRMS()
//      ,dAzi_weight_nLayer2->GetMean(), dAzi_weight_nLayer2->GetRMS()
//      );
//
//c7.cd(3);
//dR_weight->SetLineColor(kBlack);
//dR_weight_nLayer25->SetLineColor(kRed);
//dR_weight_nLayer10->SetLineColor(kMagenta+4);
//dR_weight_nLayer5->SetLineColor(kMagenta+2);
//dR_weight_nLayer2->SetLineColor(kMagenta);
//dR_weight->Draw();
//dR_weight_nLayer25->Draw("same");
//dR_weight_nLayer10->Draw("same");
//dR_weight_nLayer5->Draw("same");
//dR_weight_nLayer2->Draw("same");
//
//printf("dR (mean, rms) 50:(%f, %f)\t25:(%f, %f)\t10:(%f, %f)\t5:(%f, %f)\t2:(%f, %f)\n"
//      ,dR_weight->GetMean(), dR_weight->GetRMS()
//      ,dR_weight_nLayer25->GetMean(), dR_weight_nLayer25->GetRMS()
//      ,dR_weight_nLayer10->GetMean(), dR_weight_nLayer10->GetRMS()
//      ,dR_weight_nLayer5->GetMean(), dR_weight_nLayer5->GetRMS()
//      ,dR_weight_nLayer2->GetMean(), dR_weight_nLayer2->GetRMS()
//      );
//
//c7.SaveAs("recoAnalysis_7.C");
//
//TCanvas c8("c8","c8",800,600);
//c8.Divide(3,1);
//c8.cd(1);
////dRecoZen_weight->SetLineColor(kBlack);
//dRecoZen_weight_nLayer25->SetLineColor(kRed);
//dRecoZen_weight_nLayer10->SetLineColor(kMagenta+4);
//dRecoZen_weight_nLayer5->SetLineColor(kMagenta+2);
//dRecoZen_weight_nLayer2->SetLineColor(kMagenta);
////dRecoZen_weight->Draw("AL");
//dRecoZen_weight_nLayer25->Draw();
//dRecoZen_weight_nLayer10->Draw("same");
//dRecoZen_weight_nLayer5->Draw("same");
//dRecoZen_weight_nLayer2->Draw("same");
//
//printf("dRecoZen (mean, rms) 25:(%f, %f)\t10:(%f, %f)\t5:(%f, %f)\t2:(%f, %f)\n"
//      //,dRecoZen_weight->GetMean(), dRecoZen_weight->GetRMS();
//      ,dRecoZen_weight_nLayer25->GetMean(), dRecoZen_weight_nLayer25->GetRMS()
//      ,dRecoZen_weight_nLayer10->GetMean(), dRecoZen_weight_nLayer10->GetRMS()
//      ,dRecoZen_weight_nLayer5->GetMean(), dRecoZen_weight_nLayer5->GetRMS()
//      ,dRecoZen_weight_nLayer2->GetMean(), dRecoZen_weight_nLayer2->GetRMS()
//      );
//
//c8.cd(2);
////dRecoAzi_weight->SetLineColor(kBlack);
//dRecoAzi_weight_nLayer25->SetLineColor(kRed);
//dRecoAzi_weight_nLayer10->SetLineColor(kMagenta+4);
//dRecoAzi_weight_nLayer5->SetLineColor(kMagenta+2);
//dRecoAzi_weight_nLayer2->SetLineColor(kMagenta);
////dRecoAzi_weight->Draw("AL");
//dRecoAzi_weight_nLayer25->Draw();
//dRecoAzi_weight_nLayer10->Draw("same");
//dRecoAzi_weight_nLayer5->Draw("same");
//dRecoAzi_weight_nLayer2->Draw("same");
//
//printf("dRecoAzi (mean, rms) 25:(%f, %f)\t10:(%f, %f)\t5:(%f, %f)\t2:(%f, %f)\n"
//      //,dRecoAzi_weight->GetMean(), dRecoAzi_weight->GetRMS();
//      ,dRecoAzi_weight_nLayer25->GetMean(), dRecoAzi_weight_nLayer25->GetRMS()
//      ,dRecoAzi_weight_nLayer10->GetMean(), dRecoAzi_weight_nLayer10->GetRMS()
//      ,dRecoAzi_weight_nLayer5->GetMean(), dRecoAzi_weight_nLayer5->GetRMS()
//      ,dRecoAzi_weight_nLayer2->GetMean(), dRecoAzi_weight_nLayer2->GetRMS()
//      );
//
//c8.cd(3);
////dRecoR_weight->SetLineColor(kBlack);
//dRecoR_weight_nLayer25->SetLineColor(kRed);
//dRecoR_weight_nLayer10->SetLineColor(kMagenta+4);
//dRecoR_weight_nLayer5->SetLineColor(kMagenta+2);
//dRecoR_weight_nLayer2->SetLineColor(kMagenta);
////dRecoR_weight->Draw("AL");
//dRecoR_weight_nLayer25->Draw();
//dRecoR_weight_nLayer10->Draw("same");
//dRecoR_weight_nLayer5->Draw("same");
//dRecoR_weight_nLayer2->Draw("same");
//
//printf("dRecoR (mean, rms) 25:(%f, %f)\t10:(%f, %f)\t5:(%f, %f)\t2:(%f, %f)\n"
//      //,dRecoR_weight->GetMean(), dRecoR_weight->GetRMS();
//      ,dRecoR_weight_nLayer25->GetMean(), dRecoR_weight_nLayer25->GetRMS()
//      ,dRecoR_weight_nLayer10->GetMean(), dRecoR_weight_nLayer10->GetRMS()
//      ,dRecoR_weight_nLayer5->GetMean(), dRecoR_weight_nLayer5->GetRMS()
//      ,dRecoR_weight_nLayer2->GetMean(), dRecoR_weight_nLayer2->GetRMS()
//      );
//
//c8.SaveAs("recoAnalysis_8.C");
/*
TCanvas c9("c9","c9",800,600);

c9.Divide(2,1);

c9.cd(1);
coherence_all->Draw();

coherence25Hist->SetLineWidth(2);
coherence25Hist->SetLineColor(kRed);
coherence10Hist->SetLineWidth(2);
coherence10Hist->SetLineColor(kMagenta+4);
coherence5Hist->SetLineWidth(2);
coherence5Hist->SetLineColor(kMagenta+2);
coherence2Hist->SetLineWidth(2);
coherence2Hist->SetLineColor(kMagenta);

coherence25Hist->Draw("same");
coherence10Hist->Draw("same");
coherence5Hist->Draw("same");
coherence2Hist->Draw("same");

c9.cd(2);

dCoherence25Hist->SetLineColor(kRed);
dCoherence10Hist->SetLineColor(kMagenta+4);
dCoherence5Hist->SetLineColor(kMagenta+2);
dCoherence2Hist->SetLineColor(kMagenta);
dCoherence25Hist->Draw();
dCoherence10Hist->Draw("same");
dCoherence5Hist->Draw("same");
dCoherence2Hist->Draw("same");

c9.SaveAs("recoAnalysis_9.C");
*/
//TCanvas c10("c10","c10",800,600);
//c10.Divide(2,2);
//
//c10.cd(1);
//dRecoZen_dCoherence_weight_nLayer25->Draw("colz");
//c10.cd(2);
//dRecoZen_dCoherence_weight_nLayer10->Draw("colz");
//c10.cd(3);
//dRecoZen_dCoherence_weight_nLayer5->Draw("colz");
//c10.cd(4);
//dRecoZen_dCoherence_weight_nLayer2->Draw("colz");
//
//c10.SaveAs("recoAnalysis_10.C");
//
//TCanvas c11("c11","c11",800,600);
//c11.Divide(2,2);
//
//c11.cd(1);
//dRecoAzi_dCoherence_weight_nLayer25->Draw("colz");
//c11.cd(2);
//dRecoAzi_dCoherence_weight_nLayer10->Draw("colz");
//c11.cd(3);
//dRecoAzi_dCoherence_weight_nLayer5->Draw("colz");
//c11.cd(4);
//dRecoAzi_dCoherence_weight_nLayer2->Draw("colz");
//
//c11.SaveAs("recoAnalysis_11.C");
//
//TCanvas c12("c12","c12",800,800);
//c12.Divide(2,1);
//c12.cd(1);
//maxPixEachLayerZenSpread->Draw();
//c12.cd(2);
//maxPixEachLayerAziSpread->Draw();
//c12.SaveAs("recoAnalysis_12.C");
//
//TCanvas c13("c13","c13",800,800);
////c12.Divide(2,1);
////c12.cd(1);
//coherence2_maxPixEachLayerZenSpread->Draw("colz");
////c12.cd(2);
////maxPixEachLayerAziSpread->Draw();
//c13.SaveAs("recoAnalysis_13.C");
//
//TCanvas c14("c14","c14",800,800);
//coherence2_maxPixEachLayerAziSpread->Draw("colz");
//c14.SaveAs("recoAnalysis_14.C");
//
//TCanvas c15("c15","c15",800,800);
//maxPixEachLayerZenAziSpread->Draw("colz");
//c15.SaveAs("recoAnalysis_15.C");
//
//TCanvas c16("c16","c16",800,800);
//c16.Divide(2,1);
//c16.cd(1);
//topMaxPixZenSpread->Draw();
//c16.cd(2);
//topMaxPixAziSpread->Draw();
//c16.SaveAs("recoAnalysis_16.C");
//
//TCanvas c17("c17","c17",800,800);
////c12.Divide(2,1);
////c12.cd(1);
//coherence2_topMaxPixZenSpread->Draw("colz");
////c12.cd(2);
////maxPixEachLayerAziSpread->Draw();
//c17.SaveAs("recoAnalysis_17.C");
//
//TCanvas c18("c18","c18",800,800);
//coherence2_topMaxPixAziSpread->Draw("colz");
//c18.SaveAs("recoAnalysis_18.C");
//
//TCanvas c19("c19","c19",800,800);
//topMaxPixZenAziSpread->Draw("colz");
//c19.SaveAs("recoAnalysis_19.C");
//
//TCanvas c20("c20","c20",800,800);
//c20.Divide(1,2);
//c20.cd(1);
//thetaRatioHist->Draw();
//c20.cd(2);
//phiRatioHist->Draw();
//c20.SaveAs("recoAnalysis_20.C");
//
//TCanvas c21("c21","c21",800,800);
//c21.Divide(2,2);
//c21.cd(1);
//coherence2_thetaRMSCutXing->Draw("colz");
//c21.cd(2);
//coherence2_phiRMSCutXing->Draw("colz");
//c21.cd(3);
//topMaxPixZenSpread_thetaRMSCutXing->Draw("colz");
//c21.cd(4);
//topMaxPixAziSpread_phiRMSCutXing->Draw("colz");
//c21.SaveAs("recoAnalysis_21.C");
//
//TCanvas c22("c22","c22",800,800);
//c22.Divide(4,4);
//for(int i=0; i<16; i++){
//c22.cd(i+1);
//trueRecAngleHist[i]->Draw();
//}
//c22.SaveAs("recoAnalysis_22.C");
//
//TCanvas c23("c23","c23",800,800);
//c23.Divide(4,4);
//for(int i=0; i<16; i++){
//c23.cd(i+1);
//recoRecAngleHist[i]->Draw();
//}
//c23.SaveAs("recoAnalysis_23.C");
//
//TCanvas c24("c24","c24",800,800);
//c24.Divide(4,4);
//for(int i=0; i<16; i++){
//c24.cd(i+1);
//dRecAngle[i]->Draw();
//}
//c24.SaveAs("recoAnalysis_24.C");
//
//TCanvas c25("c25","c25",800,800);
//c25.Divide(4,4);
//for(int i=0; i<16; i++){
//c25.cd(i+1);
//dLauAngle[i]->Draw();
//}
//c25.SaveAs("recoAnalysis_25.C");
//
//TCanvas c26("c26","c26",800,800);
////coherence2Hist_nGoodChanNorm->Draw();
////coherence2Hist->SetLineColor(kRed);
////coherence2Hist->Draw();
//
//coherence2Hist->SetTitle("Max skymap C2;C2 [unitless];Entry");
//coherence2Hist->SetLineWidth(2);
//coherence2Hist->Draw();
//
//coherence2_rf->SetLineWidth(2);
//coherence2_rf->SetLineColor(kRed);
//coherence2_rf->Draw("same");
//
//coherence2_cal->SetLineWidth(2);
//coherence2_cal->SetLineColor(kMagenta);
//coherence2_cal->Draw("same");
//
//coherence2_soft->SetLineWidth(2);
//coherence2_soft->SetLineColor(kOrange);
//coherence2_soft->Draw("same");
//
//TLegend *leg5 = new TLegend(0.1,0.7,0.48,0.9);
//leg5->SetHeader("Max Skymap Coherence"/*,"C"*/); // option "C" allows to center the header
//leg5->AddEntry(coherence2Hist,"All events","l");
//leg5->AddEntry(coherence2_rf,"RF events","l");
//leg5->AddEntry(coherence2_cal,"Cal events","l");
//leg5->AddEntry(coherence2_soft,"Soft events","l");
//
//c26.SaveAs("recoAnalysis_26.C");
//
//TCanvas c27("c27","c27",800,800);
//c27.Divide(2,2);
//c27.cd(1);
//zenDist->Draw();
//c27.cd(2);
//aziDist->Draw();
//c27.cd(3);
//zenDist_coherence2->Draw("colz");
//c27.cd(4);
//aziDist_coherence2->Draw("colz");
//c27.SaveAs("recoAnalysis_27.C");
//
//TCanvas c28("c28","c28",800,800);
//coherence2_snr_rf->SetTitle(";SNR [#sigma^{-1}_{noise}];C2 [unitless]");
//coherence2_snr_rf->Draw("colz");
////c28.SaveAs("recoAnalysis_28_2016_ARA02_vnchnl3.C");
//
//TCanvas c29("c29","c29",800,800);
////dZen_snr->SetTitle("dZen vs SNR;SNR [#signma^{-1}];dZen [#circ]");
////dZen_snr->Draw("colz");
//dAzi_snr->SetTitle("dAzi vs SNR;SNR [#signma^{-1}];dAzi [#circ]");
//dAzi_snr->Draw("colz");
////c29.SaveAs("recoAnalysis_29.C");
//
//
//TCanvas c30("c30","c30",800,800);
//c30.Divide(2,2);
//c30.cd(1);
//snrGradient_group3->Draw();
////snrGradient_group3_weirdEvents->SetLineColor(kRed);
////snrGradient_group3_weirdEvents->Draw("same");
//c30.cd(2);
//snrGradient_group4->Draw();
////snrGradient_group4_weirdEvents->SetLineColor(kRed);
////snrGradient_group4_weirdEvents->Draw("same");
//c30.cd(3);
//snrGradient_7_0->Draw();
////snrGradient_7_0_weirdEvents->SetLineColor(kRed);
////snrGradient_7_0_weirdEvents->Draw("same");
//c30.cd(4);
//snrGradient_6_1->Draw();
////snrGradient_6_1_weirdEvents->SetLineColor(kRed);
////snrGradient_6_1_weirdEvents->Draw("same");
//c30.SaveAs("recoAnalysis_30.C");
//
/*
Canvas c31("c31","c31",800,800);
c31.Divide(2,1);
c31.cd(1);
timeSequenceHist->Draw();
c31.cd(2);
TH1F *timeSequenceCumuHist = new TH1F("timeSeequenceCumuHist","timeSeqeucneCumuHist",3e4, 0, 1e4);
double integral = timeSequenceHist->Integral();
for(int i=0; i<=timeSequenceCumuHist->GetNbinsX(); i++){
   timeSequenceCumuHist->SetBinContent(i, timeSequenceHist->Integral(i, timeSeqenceHist->GetNbinsX())/integral);
}
timeSequenceCumuHist->Draw();

c31.SaveAs("recoAnalysis_31.C");
*/

/*
TFile fp("ARA02_run6395_calpulser_maxPixCoherence.root","update");
coherence_cal->SetName("noChan7");
coherence_cal->SetTitle("No chan 7;Coherence [unitless];Entry");
coherence_cal->Rebin();
coherence_cal->SetLineColor(kGray);
coherence_cal->Write();
fp.Close();
*/
/*
TCanvas c32("c32","c32",1200,800);
c32.Divide(3,1);
c32.cd(1);
bipolarRatio_dT->Draw("colz");
//biPolarRatio_dT->SetTitle(";

c32.cd(2);
impulsivity_dT->Draw("colz");

c32.cd(3);
impulsivity_bipolarRatio->Draw("colz");

c32.SaveAs("recoAnalysis_32.C");
*/
TCanvas c33("c33","c33",1200,800);
c33.Divide(3,1);
c33.cd(1);
//bipolarRatio_dT->Draw("colz");
//biPolarRatio_dT->SetTitle(";
dTHist->Draw();

c33.cd(2);
//impulsivity_dT->Draw("colz");
impAvgHist->Draw();

c33.cd(3);
//impulsivity_bipolarRatio->Draw("colz");
impPassThresHist->Draw();

//c33.SaveAs("recoAnalysis_33.C");

TCanvas c34("c34","c34",800,800);
impulsivityHist->Draw();
//c34.SaveAs("recoAnalysis_34.C");

TCanvas c35("c35","c35",800,800);
iterCWCountHist->Draw();
c35.SaveAs("recoAnalysis_35.C");

TCanvas c36("c36","c36",800,800);
c36.Divide(4,4);
for(int ch=0; ch<16; ch++){
   c36.cd(ch+1);
   avgPowerRatioHist[ch]->Draw();
}
c36.SaveAs("recoAnalysis_36.C");

TCanvas c37("c37","c37",800,800);
coherence_snr_cw->Draw("colz");
c37.SaveAs("recoAnalysis_37.C");

printf("Nentries: %d\trfEventCount: %d\tcalEventCount: %d\tsoftEventCount: %d\n", Nentries, rfEventCount, calEventCount, softEventCount);



//delete dummyData;
delete runInfoTree;
delete dataTree;
delete recoSettingsTree;

return 0;
}


float statisticalRadiusReco(int arraySize, const float * const radii, const float * const value, float trueRadius,
                          TH1F * const statRReco_weight, TH1F * const statRReco_max, TH1F * const statRReco_eventStack){

//int arraySize = sizeof( value ) / sizeof( value[0] );
//cout<<"arraySize: "<<arraySize<<endl;

// Plot 1 - weighted radius - true radius
// Plot 2 - max value radius - true radius
// Plot 3 - stacked (layer radius - true radius)
float weightedRadius = 0.f;
float totalValue = 0.f;
float max = 0.f;
float maxRadius = 0.f;

for(int i=0; i<arraySize; i++){
   weightedRadius += value[i] * radii[i];
   totalValue += value[i];
   if( value[i] > max ){ max = value[i]; maxRadius = radii[i]; }

}
weightedRadius /= totalValue;
statRReco_weight->Fill(weightedRadius - trueRadius);
statRReco_max->Fill(maxRadius - trueRadius);

for(int i=0; i<arraySize; i++)
   if(totalValue != 0 ) statRReco_eventStack->Fill(radii[i] - trueRadius, value[i]/totalValue);
   //else cerr<<"totalValue = 0!!\n";


//return 0;
//return weightedRadius;
return maxRadius;
}

float getMean(const vector<float>& thetas){

   int size = (int)thetas.size();
   if(size<1){ cerr<<"Vector size: "<<size<<endl; return -1e6; }
   float sum = 0.f;

   for(int i=0; i<size; i++){

   sum += thetas[i];

   }

   return sum / (float)size;


}

float getRMS(const vector<float>& thetas){



   float mean = getMean(thetas);
   float rms = 0.f;
   int size = (int)thetas.size();
   if(size <= 1){ cerr<<"Vector size: "<<size<<endl; return -1;}

   for(int i=0; i<size; i++){

   rms += ((thetas[i] - mean) * (thetas[i] - mean));

   }

   rms = sqrt( rms / (float)(size-1) );

   return rms;
}
