#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <string>
#include <algorithm>

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

//#define CORRUPT_EVENT_START_TIME 1449028179
#define CORRUPT_EVENT_START_TIME 1448485911 //ARA02 run 3 start time
#define CORRUPT_EVENT_END_EVENT_NUMBER 4
#define ZEN_BAND_MAX -41
#define ZEN_BAND_MIN -46
//#define SURFACE_CUT 30.443
//#define SURFACE_CUT 34.80372
//#define SURFACE_CUT 35.64602
//#define SURFACE_CUT_2 36.77852 //Defined from nchan>=5 iterative reco majority vote
//#define NOISY_RUN_LIST ARA02_noSuE19aceCut_noisyRuns.txt
using namespace std;

float statisticalRadiusReco(int arraySize, const float * const radii, const float * const value, float trueRadius,
                          TH1F * const statRReco_weight, TH1F * const statRReco_max, TH1F * const statRReco_eventStack);

float getMean(const vector<float>& thetas);
float getRMS(const vector<float>& thetas);
bool isNearNoisyRun(string station, const vector<int>& noisyRuns, int runNum, int plusMinusRunNum);
float getZenMaj(const vector<float>& iterZenVec, float zenRange);

int main(int argc, char **argv){

gROOT->ProcessLine("#include <vector>");


string STATION = string(argv[1]);
int type = atoi(argv[2]);
cout<<"STATION: "<<STATION<<" type: "<<type<<endl;
//string ENERGY = string(argv[3]);

ofstream outputFile(argv[3],std::ofstream::out|std::ofstream::app);

/*TChain*/TTree *recoSettingsTree/*=new TChain("recoSettingsTree")*/;
/*TChain*/TTree *dataTree/*=new TChain("dataTree")*/;
TChain *runInfoTree=new TChain("runInfoTree");

for(int i=4; i<argc; i++){

   TFile fp( argv[i] );

   if( fp.IsZombie() ){ cerr<<"File "<<argv[i]<<" is zombie. Skipping..."<<endl; continue; }
   if( fp.TestBit(TFile::kRecovered) ){ cerr<<"File "<<argv[i]<<" is recovered file. Skipping..."<<endl; continue; }

   //recoSettingsTree->Add( argv[i] );
   //dataTree->Add( argv[i] );
   runInfoTree->Add( argv[i] );
   fp.Close();
}


recoSettings *settings = new recoSettings();;
int nSideExp, nLayer;
//recoSettingsTree->SetBranchAddress("settings", &settings);
//recoSettingsTree->GetEntry(0);

//nSideExp = settings->nSideExp;
//nLayer = settings->nLayer;
//int nAnt = (string(settings->recoPolType)=="both"?16:8);
//int numIter = nAnt - settings->nchnlCut + 1;

//Healpix_Onion onion(nSideExp, nLayer, settings->layerFirstRadius, settings->layerLastRadius);
//int nDir = onion.nDir;
//cout<<"nSideExp: "<<onion.nSideExp<<" nLayer: "<<onion.nLayer<<endl;

recoData *dummyData = new recoData();
//dataTree->SetBranchAddress("summary", &dummyData);
/*
int runEventCount, trigEventCount, recoEventCount, cutWaveEventCount, nonIncreasingSampleTimeEventCount, cutWaveAndNonIncreasingEventCount
, mistaggedSoftEventCount, offsetBlockEventCount, nchnlFilteredEventCount, corruptFirst3EventCount, cwFilteredEventCount;
int runStartTime, runEndTime;
int runRFEventCount, runCalEventCount, runSoftEventCount;
double weightedOffsetBlockEventCount, weightedImpulsivityFilteredEventCount, weightedTrigEventCount;
*/
int runEventCount, runRFEventCount, runCalEventCount, runSoftEventCount, trigEventCount, recoEventCount, utime_runStart, utime_runEnd, cutWaveEventCount, nonIncreasingSampleTimeEventCount, cutWaveAndNonIncreasingEventCount, mistaggedSoftEventCount, offsetBlockEventCount, cwFilteredEventCount, nchnlFilteredEventCount, impulsivityFilteredEventCount, corruptFirst3EventCount, corruptD1EventCount;

double weightedTrigEventCount, weightedRecoEventCount, weightedOffsetBlockEventCount, weightedNchnlFilteredEventCount, weightedCWFilteredEventCount, weightedImpulsivityFilteredEventCount;

int runStartTime, runEndTime;

/*
runInfoTree->SetBranchAddress("runEventCount", &runEventCount);
runInfoTree->SetBranchAddress("runRFEventCount", &runRFEventCount);
runInfoTree->SetBranchAddress("runCalEventCount", &runCalEventCount);
runInfoTree->SetBranchAddress("runSoftEventCount", &runSoftEventCount);
runInfoTree->SetBranchAddress("trigEventCount", &trigEventCount);
runInfoTree->SetBranchAddress("recoEventCount", &recoEventCount);
runInfoTree->SetBranchAddress("utime_runStart", &runStartTime);
runInfoTree->SetBranchAddress("utime_runEnd", &runEndTime);
runInfoTree->SetBranchAddress("cutWaveEventCount", &cutWaveEventCount);
runInfoTree->SetBranchAddress("nonIncreasingSampleTimeEventCount", &nonIncreasingSampleTimeEventCount);
runInfoTree->SetBranchAddress("cutWaveAndNonIncreasingEventCount", &cutWaveAndNonIncreasingEventCount);
runInfoTree->SetBranchAddress("mistaggedSoftEventCount", &mistaggedSoftEventCount);
runInfoTree->SetBranchAddress("offsetBlockEventCount", &offsetBlockEventCount);
runInfoTree->SetBranchAddress("nchnlFilteredEventCount", &nchnlFilteredEventCount);
runInfoTree->SetBranchAddress("corruptFirst3EventCount", &corruptFirst3EventCount);
runInfoTree->SetBranchAddress("cwFilteredEventCount", &cwFilteredEventCount);
runInfoTree->SetBranchAddress("weightedTrigEventCount", &weightedTrigEventCount);
runInfoTree->SetBranchAddress("weightedOffsetBlockEventCount", &weightedOffsetBlockEventCount);
runInfoTree->SetBranchAddress("weightedImpulsivityFilteredEventCount", &weightedImpulsivityFilteredEventCount);


double totalRunEventCount, totalTrigEventCount, totalRecoEventCount, totalCutWaveEventCount, totalNonIncreasingSampleTimeEventCount, totalCutWaveAndNonIncreasingEventCount, totalRFEventCount, totalCalEventCount, totalSoftEventCount, totalMistaggedSoftEventCount, totalOffsetBlockEventCount, totalNchnlFilteredEventCount, totalCorruptFirst3EventCount, totalCWFilteredEventCount, totalImpulsivityFilteredEventCount;
int totalLiveTime = 0;
totalRunEventCount = totalTrigEventCount = totalRecoEventCount = totalCutWaveEventCount = totalNonIncreasingSampleTimeEventCount = totalCutWaveAndNonIncreasingEventCount = totalRFEventCount = totalCalEventCount = totalSoftEventCount = totalMistaggedSoftEventCount = totalOffsetBlockEventCount = totalNchnlFilteredEventCount = totalCorruptFirst3EventCount = totalCWFilteredEventCount = totalImpulsivityFilteredEventCount = 0;
*/
runInfoTree->SetBranchAddress("runEventCount", &runEventCount);
runInfoTree->SetBranchAddress("runRFEventCount", &runRFEventCount);
runInfoTree->SetBranchAddress("runCalEventCount", &runCalEventCount);
runInfoTree->SetBranchAddress("runSoftEventCount", &runSoftEventCount);
runInfoTree->SetBranchAddress("trigEventCount", &trigEventCount);
runInfoTree->SetBranchAddress("recoEventCount", &recoEventCount);
runInfoTree->SetBranchAddress("utime_runStart", &runStartTime);
runInfoTree->SetBranchAddress("utime_runEnd", &runEndTime);
runInfoTree->SetBranchAddress("cutWaveEventCount", &cutWaveEventCount);
runInfoTree->SetBranchAddress("nonIncreasingSampleTimeEventCount", &nonIncreasingSampleTimeEventCount);
runInfoTree->SetBranchAddress("cutWaveAndNonIncreasingEventCount", &cutWaveAndNonIncreasingEventCount);
runInfoTree->SetBranchAddress("mistaggedSoftEventCount", &mistaggedSoftEventCount);
runInfoTree->SetBranchAddress("offsetBlockEventCount", &offsetBlockEventCount);
runInfoTree->SetBranchAddress("cwFilteredEventCount", &cwFilteredEventCount);
runInfoTree->SetBranchAddress("nchnlFilteredEventCount", &nchnlFilteredEventCount);
runInfoTree->SetBranchAddress("impulsivityFilteredEventCount", &impulsivityFilteredEventCount);
runInfoTree->SetBranchAddress("corruptFirst3EventCount", &corruptFirst3EventCount);
runInfoTree->SetBranchAddress("corruptD1EventCount", &corruptD1EventCount);
runInfoTree->SetBranchAddress("weightedTrigEventCount", &weightedTrigEventCount);
runInfoTree->SetBranchAddress("weightedRecoEventCount", &weightedRecoEventCount); //
runInfoTree->SetBranchAddress("weightedOffsetBlockEventCount", &weightedOffsetBlockEventCount);
runInfoTree->SetBranchAddress("weightedNchnlFilteredEventCount", &weightedNchnlFilteredEventCount);
runInfoTree->SetBranchAddress("weightedCWFilteredEventCount", &weightedCWFilteredEventCount);//
runInfoTree->SetBranchAddress("weightedImpulsivityFilteredEventCount", &weightedImpulsivityFilteredEventCount);

int totalRunEventCount, totalRFEventCount, totalCalEventCount, totalSoftEventCount, totalTrigEventCount, totalRecoEventCount, totalCutWaveEventCount, totalNonIncreasingSampleTimeEventCount, totalCutWaveAndNonIncreasingEventCount, totalMistaggedSoftEventCount, totalOffsetBlockEventCount, totalCWFilteredEventCount, totalNchnlFilteredEventCount, totalImpulsivityFilteredEventCount, totalCorruptFirst3EventCount, totalCorruptD1EventCount;

totalRunEventCount = totalRFEventCount = totalCalEventCount = totalSoftEventCount = totalTrigEventCount = totalRecoEventCount = totalCutWaveEventCount = totalNonIncreasingSampleTimeEventCount = totalCutWaveAndNonIncreasingEventCount = totalMistaggedSoftEventCount = totalOffsetBlockEventCount = totalCWFilteredEventCount = totalNchnlFilteredEventCount = totalImpulsivityFilteredEventCount = totalCorruptFirst3EventCount = totalCorruptD1EventCount = 0;

double totalWeightedTrigEventCount, totalWeightedRecoEventCount, totalWeightedOffsetBlockEventCount, totalWeightedNchnlFilteredEventCount, totalWeightedCWFilteredEventCount, totalWeightedImpulsivityFilteredEventCount;

totalWeightedTrigEventCount = totalWeightedRecoEventCount = totalWeightedOffsetBlockEventCount = totalWeightedNchnlFilteredEventCount = totalWeightedCWFilteredEventCount = totalWeightedImpulsivityFilteredEventCount = 0.;

int totalLiveTime=0;

for(int run=0; run<runInfoTree->GetEntries(); run++){

   runInfoTree->GetEntry(run);

   /*
   totalRunEventCount += runEventCount;
   totalRFEventCount  += runRFEventCount;
   totalCalEventCount += runCalEventCount;
   totalSoftEventCount += runSoftEventCount;
   totalTrigEventCount += weightedTrigEventCount;
   totalRecoEventCount += recoEventCount;
   totalCutWaveEventCount += cutWaveEventCount;
   totalNonIncreasingSampleTimeEventCount += nonIncreasingSampleTimeEventCount;
   totalCutWaveAndNonIncreasingEventCount += cutWaveAndNonIncreasingEventCount;
   totalMistaggedSoftEventCount += mistaggedSoftEventCount;
   totalOffsetBlockEventCount += weightedOffsetBlockEventCount;
   totalImpulsivityFilteredEventCount += weightedImpulsivityFilteredEventCount;
   totalCWFilteredEventCount += cwFilteredEventCount;
   totalNchnlFilteredEventCount += nchnlFilteredEventCount;
   totalCorruptFirst3EventCount += corruptFirst3EventCount;
*/

   totalRunEventCount                         += runEventCount;
   totalRFEventCount                          += runRFEventCount;
   totalCalEventCount                         += runCalEventCount;
   totalSoftEventCount                        += runSoftEventCount;
   totalTrigEventCount                        += trigEventCount;
   totalRecoEventCount                        += recoEventCount;
   totalCutWaveEventCount                     += cutWaveEventCount;
   totalNonIncreasingSampleTimeEventCount     += nonIncreasingSampleTimeEventCount;
   totalCutWaveAndNonIncreasingEventCount     += cutWaveAndNonIncreasingEventCount;
   totalMistaggedSoftEventCount               += mistaggedSoftEventCount;
   totalOffsetBlockEventCount                 += offsetBlockEventCount;
   totalCWFilteredEventCount                  += cwFilteredEventCount;
   totalNchnlFilteredEventCount               += nchnlFilteredEventCount;
   totalImpulsivityFilteredEventCount         += impulsivityFilteredEventCount;
   totalCorruptFirst3EventCount               += corruptFirst3EventCount;
   totalCorruptD1EventCount                   += corruptD1EventCount;
   totalWeightedTrigEventCount                += weightedTrigEventCount;
   totalWeightedRecoEventCount                += weightedRecoEventCount;
   totalWeightedOffsetBlockEventCount         += weightedOffsetBlockEventCount;
   totalWeightedNchnlFilteredEventCount       += weightedNchnlFilteredEventCount;
   totalWeightedCWFilteredEventCount          += weightedCWFilteredEventCount;
   totalWeightedImpulsivityFilteredEventCount += weightedImpulsivityFilteredEventCount;

   if(runEndTime < runStartTime){ cerr<<"Run "<<run<<" livetime error. Skipping...\n"; continue; }
   else totalLiveTime += (runEndTime - runStartTime);


}

printf("totalRunEventCount: %d\ttotalTrigEventCount: %d\ttotalRecoEventCount: %d\ntotalLiveTime :%ds\n", totalRunEventCount, totalTrigEventCount, totalRecoEventCount, totalLiveTime);
printf("totalRFEventCount: %d\ttotalCalEventCount: %d\ttotalSoftEventCount: %d\n", totalRFEventCount, totalCalEventCount, totalSoftEventCount);
printf("totalCutWaveEventCount: %d\ttotalNonIncreasingSampleTimeEventCount: %d\ttotalCutWaveAndNonIncreasingEventCount: %d\ncutWave ratio: %f\tnonIncreasing ratio: %f\tcutWaveAndNonIncreasing ratio: %f\n", totalCutWaveEventCount, totalNonIncreasingSampleTimeEventCount, totalCutWaveAndNonIncreasingEventCount, (float)totalCutWaveEventCount/(float)totalRunEventCount, (float)totalNonIncreasingSampleTimeEventCount/(float)totalRunEventCount, (float)totalCutWaveAndNonIncreasingEventCount/(float)totalRunEventCount);
printf("totalMistaggedSoftEventCount: %d\ttotalOffsetBlockEventCount: %d\ttotalCWFilteredEventCount: %d\ttotalNchnlFilteredEventCount: %d\ttotalCorruptFirst3EventCount: %d\nmistag soft ratio: %f\toffset ratio: %f\tCW ratio: %f\tnchnl ratio: %f\tfirst3 ratio: %f\n", totalMistaggedSoftEventCount, totalOffsetBlockEventCount, totalCWFilteredEventCount, totalNchnlFilteredEventCount, totalCorruptFirst3EventCount, (float)totalMistaggedSoftEventCount/(float)totalRunEventCount, (float)totalOffsetBlockEventCount/(float)totalRunEventCount, (float)totalCWFilteredEventCount/(float)totalRunEventCount, (float)totalNchnlFilteredEventCount/(float)totalRunEventCount, (float)totalCorruptFirst3EventCount/(float)totalRunEventCount);

//int Nentries = dataTree->GetEntries();
//cout<<"Total number of dataTree events: "<<Nentries<<endl;

double rfEventCount, calEventCount, softEventCount;
rfEventCount = calEventCount = softEventCount = 0.;

ifstream list;
//list.open("ARA02_noSuE19aceCut_noisyRuns.txt");
//list.open("ARA02_suE19aceCut_noisyRuns.txt");
//list.open("ARA02_surfaceEvents_noisyRuns.txt");
//list.open("ARA02_impulsivityFilter1_surfaceEvents_noisyRuns.txt");
//list.open("ARA02_impulsivityFilter1_noMaskSat_surfaceEvents_noisyRuns.txt");
//list.open("ARA02_impulsivityFilter1_extendedThermalCut_noMaskSat_surfaceEvents_noisyRuns.txt");
//list.open("ARA02_impulsivityFilter1_thermalImpulsivityCut_noMaskSat_surfaceEvents_noisyRuns.txt");
//list.open("ARA02_vnchnl3NoMasking_noMaskSat_surfaceEvents_noisyRuns.txt");
//list.open("ARA02_vnchnl3NoMasking_noMaskSat_impCut_surfaceEvents_noisyRuns.txt");
//list.open("ARA02_vnchnl3NoMasking_noMaskSat_impCut_surfaceEvents_noisyRuns.txt");
//list.open("ARA02_vnchnl3NoMasking_beforeImpCut_noMaskSat_surfaceEvents_noisyRuns.txt");
//list.open("ARA02_vnchnl3NoMasking_beforeImpCut_noMaskSat_2SurfaceCut_surfaceEvents_noisyRuns.txt");
list.open("ARA02_vnchnl3NoMasking_beforeImpCut_noMaskSat_snrMode1_ch6Fit2Corr_2SurfaceCut_surfaceEvents_noisyRuns.txt");
vector<int> listOfRuns;
//vector<int> listOfEvents;
string line;
int run, event;
char line_char[200];

if (list.is_open() ){

   while (list.good()){
      getline(list, line, '\n');
      if (line == "") break;
      run = stoi(line);
      //getline(list, line, ',');
      //event = stoi(line);
      //getline(list, line, '\n');
      cout<<"run: "<<run/*<<" event: "<<event*/<<endl;
      listOfRuns.push_back(run);
      //listOfEvents.push_back(event);
   }
}

list.close();

int numNoisyRuns = listOfRuns.size();
//vector< vector<int> > noisyRun;


/*************    START LOOPING THROUGH EVENTS *********************************************************************
 * we will have
 * 1. Deep pulser time cut
 * 2. Thermal cut
 * 3. calpulser cut
 * 4.
 * Then a list fo remaining events will be produced for iterative reconstruction.
 *******************************************************************************************************************/

bool passThermalCut, passSurfaceCut, passCalpulserCut, passDeepPulserCut, passCorruptionCut, passNoisyRunCut, passThermalImpulsivityCut, passCWCut, passNumSatChanCut, passSurfaceCut_2, passCalpulserTimeCut;
double nPassCorruption, nPassThermalCut, nPassSurfaceCut, nPassDeepPulserCut, nPassNoisyRunCut, nPassCalpulserCut, nPassImpulsivityCut, nPassHighPassFilter, nPassThermalImpulsivityCut, nPassCWCut, nPassNumSatChanCut, nPassSurfaceCut_2, nPassCalpulserTimeCut;
double nCut0, nCut1, nCut2, nCut3, nCut4, nCut5, nCut6, nCut7, nCut3p5, nCut1p5, nCut0p5, nCut6p5, nCut4p5;
nPassCorruption = nPassThermalCut = nPassSurfaceCut = nPassDeepPulserCut = nPassNoisyRunCut = nPassCalpulserCut = nPassImpulsivityCut = nPassHighPassFilter = nPassThermalImpulsivityCut = nPassCWCut = nPassNumSatChanCut = nPassSurfaceCut_2 = nPassCalpulserTimeCut = 0.;
nCut0 = nCut1 = nCut2 = nCut3 = nCut4 = nCut5 = nCut6 = nCut7 = nCut3p5 = nCut1p5 = nCut0p5 = nCut6p5 = nCut4p5 = 0.;

bool passImpulsivityCut;
double nRecoveredByImp = 0.;

bool passHighPassFilter;
double highPassFreq = 170; //170MHz

bool passSNRCut;
double nPassSNRCut = 0.;

float theta, phi, r;
float zen_bestHypo, azi_bestHypo;
bool inBand;
double inBand_all, inBand_pass, outOfBand_all, outOfBand_pass;
inBand_all = inBand_pass = outOfBand_all = outOfBand_pass = 0.;
//vnchnl3NoMasking + uniformWindow + type 2 & 5 macPaddedSample 3072
/*
double snrCut_inBand[5]={10.536000, 9.672000, 9.744000,10.752000, 10.536000};
double coherenceCut_inBand[5]={0.114769, 0.116136, 0.114192, 0.103131, 0.102636};
double snrCut_outOfBand[5]={12.336000, 9.888000, 9.384000, 12.192000, 12.120000};
double coherenceCut_outOfBand[5]={0.110100, 0.112820, 0.110953, 0.099807, 0.099458};
*/
/*
double snrCut_inBand[5]={9.179000, 9.341000,8.783000,10.332500, 10.785500};
double coherenceCut_inBand[5]={0.109637, 0.111487,0.106270,0.100413,0.098118};
double snrCut_outOfBand[5]={9.816500, 9.222500,8.762000,10.253000, 10.277000};
double coherenceCut_outOfBand[5]={0.109582,0.105366,0.101697,0.097897, 0.094119};
*/
/*
//vnchnl0 + impulsivityFilter1
double snrCut_inBand[5]={7.266250, 8.107750, 8.029000, 7.693750, 8.827750};
double coherenceCut_inBand[5]={0.114491, 0.110969, 0.105926, 0.100878, 0.100461};
double snrCut_outOfBand[5]={7.468750, 6.944500, 7.950250, 7.113250, 7.673500};
double coherenceCut_outOfBand[5]={0.105733, 0.108887, 0.106223, 0.096189, 0.093352};
*/
/*
//vnchnl0 + impulsivityFilter1 + noMaskSat
double snrCut_inBand[5]={7.423750,7.489000,7.840000,6.994000,7.583500};
double coherenceCut_inBand[5]={0.112664,0.111144,0.105927, 0.099630,0.097887};
double snrCut_outOfBand[5]={9.471250,6.991750,7.797250,7.115500,7.669000};
double coherenceCut_outOfBand[5]={0.108141,0.108556,0.106329,0.096141,0.092980};
*/
//vnchnl3 + noMaskSat + flattenSat(sig)
/*
double snrCut_inBand[5]={8.970000,9.116000,8.640000,9.360000,9.480000 };
double coherenceCut_inBand[5]={0.108971,0.111546,0.106706,0.099455,0.098104};
double snrCut_outOfBand[5]={9.000000,9.128000,8.704000,9.524000,10.028000};
double coherenceCut_outOfBand[5]={0.103402,0.105563,0.101873,0.094399,0.092622};
*/

/* Extended each cut by 20% */
/*
double fac = 1.2;
for(int c=0; c<5; c++){

   snrCut_inBand[c] *= fac;
   coherenceCut_inBand[c] *= fac;
   snrCut_outOfBand[c] *= fac;
   coherenceCut_outOfBand[c] *= fac;

}
*/

//double coherenceCut_inBand[5] = {0.108008, 0.108996,  0.108338, 0.099024, 0.097075};
//double coherenceCut_outOfBand[5] = {0.103715, 0.107195, 0.103552, 0.094714, 0.095006};

double coherenceCut_inBand[5] = {0.10800825822688973, 0.10932663315269384, 0.10838951766838757, 0.09902432045234812, 0.09707451485364714};
double coherenceCut_outOfBand[5] = {0.10371485739304964, 0.10727052346393086, 0.10355164968621144, 0.09471396819592758, 0.09505296821439854};

//double snrCut[5] = {7.603017232977032,  8.23245102065904,  7.724885855976288,  8.016287951854366, 9.13553213761174};
//double snrCut[5] = {7.603815004511399, 8.328126752203419, 7.718046701397855, 8.016287951854366, 8.817297373211487};
double snrCut[5] = {8.065051326524221, 8.937513009286487,  8.377493977106527, 8.551617773166683, 9.624490817371045};

TH1F *impulsivityHist_nMinusCW = new TH1F("impulsivityHist_nMinusCW","impulsivityHist_nMinusCW",1000, -2, 2);
TH2F *c_vs_snr_hist_nMinusThermal = new TH2F("c_vs_snr_hist_nMinusThermal","c_vs_snr_hist_nMinusThermal",400,0,40,1000,0,1);
TH1F *impulsivityHist_nMinusImp = new TH1F("impulsivityHist_nMinusImp","impulsivityHist_nMinusImp",1000, -2, 2);
TH2F *zen_azi_nMinusCal = new TH2F("zen_azi_nMinusCal", "zen_azi_nMinusCal", 360/0.4, 0, 360, 180/0.4, -90, 90);
TH1F *zen_nMinusSurface = new TH1F("zen_nMinusSurface", "zen_nMinusSurface", 180/0.4, -90, 90);

TH1F *coherence_nMinusThermal = new TH1F("coherence_nMinusThermal","coherence_nMinusThermal",1000,0,1);
TH1F *snr_nMinusSNR = new TH1F("snr_nMinusSNR","snr_nMinusSNR",400,0,40);

TH2F *c_vs_snr_hist = new TH2F("c_vs_snr_hist","c_vs_snr_hist",400,0,40,1000,0,1);

TH2F *c_vs_imp = new TH2F("c_vs_imp", "c_vs_imp", 1000, 0, 1, 1000, 0, 1);

double snrCutValue, coherenceCutValue;
float coherence, snr;

int startRun=0;
int endRun  =9000;
TH1I *runHist = new TH1I("runHist","runHist", endRun-startRun+1, startRun-0.5, endRun+0.5);
TH1I *surfaceRunHist = new TH1I("surfaceRunHist","surfaceRunHist", endRun-startRun+1, startRun-0.5, endRun+0.5);

int maxFreqBin[16];
double maxFreqPower[16];
double freqBinWidth_V, freqBinWidth_H;
double maxCountFreq_V, maxCountFreq_H;
int minBin, maxBin, len;
int *freqCount_V, *freqCount_H;
bool isCW;
double isCWCount = 0;
//int passCWCut = 0;

TH2F *coherence_snr_cw = new TH2F("coherence_snr_cw","coherence_snr_cw",400,0,40,1000,0,1);
TH1F *snr_cw = new TH1F("snr_cw","snr_cw",400,0,40);

TH1F *_snrHist = new TH1F("_snrHist","_snrHist",400,0,40);
TH1F *_snrCumuHist;

TH1F *impulsivityHist_max = new TH1F("impulsivityHist_max","impulsivityHist_max",1000, -2, 2);

TH1F *impulsivityHist_3rd = new TH1F("impulsivityHist_3rd","impulsivityHist_3rd",1000, -2, 2);
TH1F *impulsivityHist_avg = new TH1F("impulsivityHist_avg","impulsivityHist_avg",1000, -2, 2);
double impulsivity[16] = {0.};
double avgImpulsivity;
int nonZeroCount = 0;
double thermalCWEventCount = 0.;
double thermalCWEventCount_both, thermalCWEventCount_V, thermalCWEventCount_H;
thermalCWEventCount_both = thermalCWEventCount_V = thermalCWEventCount_H = 0;

//double  impulsivityCut[5] = {0.29984903021901077, 0.3187978138449127, 0.295231029377421, 0.3008712367990549, 0.29173179944622674};

TH1F *maxCountFreq_V_hist = new TH1F("maxCountFreq_V_hist","maxCountFreq_V_hist", (int)(1000/0.4), 0, 1000);
TH1F *maxCountFreq_H_hist = new TH1F("maxCountFreq_H_hist","maxCountFreq_H_hist", (int)(1000/0.4), 0, 1000);

TH1F *dFHist_V = new TH1F("dFHist_V","dFHist_V", (int)(2000/0.4), -1000, 1000);
TH1F *dFHist_H = new TH1F("dFHist_H","dFHist_H", (int)(2000/0.625), -1000, 1000);


//double postThermalAvgImpulsivityCut = 0.25563478/*0.34107248*//*0.391071*//*0.38524159465351615*/;


TH1F *constantNZenHist = new TH1F("constantNZen", "constantNZen", 450, -90, 90);
TH1F *iterMajorityZenHist = new TH1F("iterMajorityZen","iterMajorityZen", 450, -90, 90);
vector<float> iterZenVec;

int iterIndex[50];
float iterMaxPixCoherenceEachLayer[50];
int iterMaxPixIdxEachLayer[50];


char histName[200];
TH1F *snrHist[5];
for(int i=0; i<5; i++){
   sprintf(histName, "hist_%d", i);
   snrHist[i] = new TH1F(histName, histName, 1000, 0, 50);
}


ARA02_cutValues *cutValues = new ARA02_cutValues();

cout<<"impCut: "<<cutValues->impCut.val<<endl;
double postThermalAvgImpulsivityCut = cutValues->impCut.val;
double surfaceCut_1 = /*SURFACE_CUT*/cutValues->surfaceCut_constantN.val;
double surfaceCut_2 = /*SURFACE_CUT_2*/cutValues->surfaceCut_iterReco.val;

TH2F *bipolarRatio_dT = new TH2F("bipolarRatio_dT","bipolarRatio_dT",800,0,800,1000,0,1);
TH2F *impulsivity_dT = new TH2F("impulsivity_dT","impulsivity_dT",800,0,800,1000,0,1);
TH2F *impulsivity_bipolarRatio = new TH2F("impulsivity_bipolarRatio","impulsivity_bipolarRatio",1000,0,1,1000,0,1);

double fftRes;
if(type<=3) fftRes = 1/(379e-9)/1e6;
else        fftRes = 1/(499e-9)/1e6;
cout<<"fftRes: "<<fftRes<<endl;

int vResBin, hResBin, xResBin;

TH1F *thetaXingHist = new TH1F("thetaXingHist","thetaXingHist",500+1,0.5,500+0.5+1);
TH1F *phiXingHist = new TH1F("phiXingHist","phiXingHist",500+1,0.5,500+0.5+1);
TH1F *avgThetaXingHist = new TH1F("avgThetaXingHist","avgThetaXingHist",500+1,0.5,500+0.5+1);
TH1F *avgPhiXingHist = new TH1F("avgPhiXingHist","avgPhiXingHist",500+1,0.5,500+0.5+1);

TH1F *inRangeThetaFracHist = new TH1F("inRangeThetaFracHist","inRangeThetaFracHist",100,0,1);
TH1F *inRangePhiFracHist = new TH1F("inRangePhiFracHist","inRangePhiFracHist",100,0,1);
TH2F *inRangeThetaPhiFracHist = new TH2F("inRangeThetaPhiFracHist","inRangeThetaPhiFracHist",100,0,1,100,0,1);

//for(int entry=0; entry<Nentries; entry++){
for(int i=4; i<argc; i++){

   string fin(argv[i]);
   int runNum = atoi(  fin.substr(fin.find("run")+3, fin.find(".root")-fin.find("run")-2).c_str() );


   //Exclude certain runs:
   //Cal sweep
   if(runNum>=3177 && runNum<=3186) continue;
   if(runNum>=3302 && runNum<=3311) continue;
   //Weird channel delays. Recon of calpulser is bad
   if(runNum>=3465 && runNum<=3504) continue;
   //First run after ~1mo downtime, D1 shows large glitches
   if(runNum==7100)                 continue;
   //Cal sweep
   if(runNum>=7625 && runNum<=7686) continue;


   //int event  = atoi(  fin.substr(fin.find("event")+5, fin.find(".txt")-fin.find("event")-4).c_str() );
   TFile fp1( argv[i] );

   if( fp1.IsZombie() ){ cerr<<"File "<<argv[i]<<" is zombie. Skipping..."<<endl; continue; }
   if( fp1.TestBit(TFile::kRecovered) ){ cerr<<"File "<<argv[i]<<" is recovered file. Skipping..."<<endl; continue; }

   recoSettingsTree = (TTree*)fp1.Get("recoSettingsTree");
   dataTree = (TTree*)fp1.Get("dataTree");

   recoSettingsTree->SetBranchAddress("settings", &settings);
   recoSettingsTree->GetEntry(0);

   nSideExp = settings->nSideExp;
   nLayer = settings->nLayer;
   int nAnt = (string(settings->recoPolType)=="both"?16:8);
   int numIter = nAnt - settings->nchnlCut + 1;

   Healpix_Onion onion(nSideExp, nLayer, settings->layerFirstRadius, settings->layerLastRadius);
   int nDir = onion.nDir;
   //cout<<"nSideExp: "<<onion.nSideExp<<" nLayer: "<<onion.nLayer<<endl;

   dataTree->SetBranchAddress("summary", &dummyData);
   int Nentries = dataTree->GetEntries();

   vector<int> maxFreqBinVec_V;
   vector<int> maxFreqBinVec_H;
   vector<int> maxFreqBinVec;
   double maxFreqArray[16];
   int maxFreqArrayPolType[16];
   int fIndex[16];
   double orderedArray[16];
   int orderedArrayPolType[16];



   for(int entry=0; entry<Nentries; entry++){
   //if(Nentries > 100) {  if(  entry % (Nentries/100) == 0  ){ cout<<"Progess: "<<entry / (Nentries/100) <<"%\n"; } }
   dataTree->GetEntry(entry);
   //cout<<"eventTrigType: "<<dummyData->eventTrigType<<endl;
   //if(dummyData->eventNumber != 127378) continue;

   //Exclude the offset-block events and block-gap events
   if(runNum==2889 && dummyData->eventNumber==108253) continue;
   if(runNum==2759 && dummyData->eventNumber==8625) continue;
   if(runNum==4838 && dummyData->eventNumber==15842) continue;
   if(runNum==6842 && dummyData->eventNumber==324) continue;
   //block-gap
   if(runNum==4429 && dummyData->eventNumber==34200) continue;




   if(dummyData->eventTrigType == 0) rfEventCount+=dummyData->weight;
   else if (dummyData->eventTrigType == 1) calEventCount+=dummyData->weight;
   else if (dummyData->eventTrigType == 2) softEventCount+=dummyData->weight;
   else { cerr<<"Event "<<entry<<" eventTrigType undefined! Skipping...\n"; continue; }

   passThermalCut = passSurfaceCut = passCalpulserCut = passDeepPulserCut = passCorruptionCut = false;
   passImpulsivityCut = false;
   passNoisyRunCut = true;
   passHighPassFilter = false;
   passThermalImpulsivityCut = false;
   passCWCut = false;
   passNumSatChanCut = false;
   passSurfaceCut_2 = false;
   passCalpulserTimeCut = false;

   passSNRCut = false;

   isCW = false;
   bool isVpolCW, isHpolCW, isXpolCW;
   isVpolCW = isHpolCW = isXpolCW = false;
   int maxCountFreqBin_V, maxCountFreqBin_H;
   maxCountFreqBin_V = maxCountFreqBin_H = 0;
   int cwBinThres = 3;

   //int cwPol = -1; //-1: not CW, 0: Vpol CW, 1: Hpol CW
   maxFreqBinVec_V.clear();
   maxFreqBinVec_H.clear();
   maxFreqBinVec.clear();

   //isCW = isCW_coincidence(isVpolCW, isHpolCW, maxCountFreqBin_V, maxCountFreqBin_H, dummyData, cwBinThres);
/*
   for(int i=0; i<8; i++){
      //cout<<"ch: "<<i<<" maxFreqBin: "<<dummyData->maxFreqBin[i]<<" maxFreq: "<<dummyData->maxFreqBin[i] * dummyData->freqBinWidth_V<<" maxFreqPower: "<<dummyData->maxFreqPower[i]<<" ";
      maxFreqBinVec_V.push_back(dummyData->maxFreqBin[i]);
   }
   //cout<<endl;
   for(int i=8; i<16; i++){
      //cout<<"ch: "<<i<<" maxFreqBin: "<<dummyData->maxFreqBin[i]<<" maxFreq: "<<dummyData->maxFreqBin[i] * dummyData->freqBinWidth_H<<" maxFreqPower: "<<dummyData->maxFreqPower[i]<<" ";
      maxFreqBinVec_H.push_back(dummyData->maxFreqBin[i]);
   }
   //cout<<endl;
   //cout<<"maxCountFreq_V: "<<dummyData->maxCountFreq_V<<" maxCountFreq_H: "<<dummyData->maxCountFreq_H<<endl;

   minBin  =  *min_element(maxFreqBinVec_V.begin(), maxFreqBinVec_V.end());
   maxBin  =  *max_element(maxFreqBinVec_V.begin(), maxFreqBinVec_V.end());
   len = maxBin - minBin + 1;
   freqCount_V = new int [len];
   std::fill(&freqCount_V[0], &freqCount_V[len], 0);

   for(int i=0; i<8; i++) freqCount_V[dummyData->maxFreqBin[i]-minBin] += 1;

   for(int i=0; i<len; i++){
      if(freqCount_V[i] > 2) { isVpolCW = true; maxCountFreqBin_V = minBin+i; }
   }

   if(!isVpolCW){
   if (freqCount_V[0] > 1){
      if (freqCount_V[1] > 0 ) { isVpolCW = true; maxCountFreqBin_V = minBin+0; }
   } else if ( freqCount_V[len-1] > 1){
      if (freqCount_V[len-2] > 0 ) { isVpolCW = true; maxCountFreqBin_V = minBin+len-1; }
   }
   }
   if (!isVpolCW){
      for(int i=1; i<len-1; i++){
         if(freqCount_V[i]>1){
            if(freqCount_V[i-1] > 0 || freqCount_V[i+1] > 0){
               isVpolCW = true;
               maxCountFreqBin_V = minBin+i;
               break;
            }
         }
      }
   }


   minBin  =  *min_element(maxFreqBinVec_H.begin(), maxFreqBinVec_H.end());
   maxBin  =  *max_element(maxFreqBinVec_H.begin(), maxFreqBinVec_H.end());
   len = maxBin - minBin + 1;
   freqCount_H = new int [len];
   std::fill(&freqCount_H[0], &freqCount_H[len], 0);

  for(int i=8; i<16; i++) freqCount_H[dummyData->maxFreqBin[i]-minBin] += 1;

   for(int i=0; i<len; i++){
      if(freqCount_H[i] > 2) { isHpolCW = true; maxCountFreqBin_H = minBin+i; }
   }

   if(!isHpolCW){
   if (freqCount_H[0] > 1){
      if (freqCount_H[1] > 0 ) { isHpolCW = true; maxCountFreqBin_H = minBin+0; }
   } else if ( freqCount_H[len-1] > 1){
      if (freqCount_H[len-2] > 0 ) { isHpolCW = true; maxCountFreqBin_H = minBin+len-1; }
   }
   }
   if (!isHpolCW){
      for(int i=1; i<len-1; i++){
         if(freqCount_H[i]>1){
            if(freqCount_H[i-1] > 0 || freqCount_H[i+1] > 0){
               isHpolCW = true;
               maxCountFreqBin_H = minBin+i;
               break;
            }
         }
      }
   }


*/

  //isCW = isCW_freqWindow(isVpolCW, isHpolCW, isXpolCW, dummyData, fftRes);

//      for(int ch=0; ch<16; ch++){
//
//         //maxFreqBinVec.push_back(dummyData->maxFreqBin[ch]);
//         maxFreqArray[ch] = dummyData->maxFreqBin[ch] * (ch<8?dummyData->freqBinWidth_V:dummyData->freqBinWidth_H);
//         if(ch<8) maxFreqArrayPolType[ch] = 0;//maxFreqArray_V[ch] =  dummyData->maxFreqBin[ch] * dummyData->freqBinWidth_V;
//         else     maxFreqArrayPolType[ch] = 1;//maxFreqArray_H[ch] =  dummyData->maxFreqBin[ch] * dummyData->freqBinWidth_H;
//
//      }
//
//
//
//      TMath::Sort(16, maxFreqArray, fIndex, kFALSE);
//      //TMath::Sort(8, maxFreqArray_V, fIndex_V, kFALSE);
//      //TMath::Sort(8, maxFreqArray_H, fIndex_H, kFALSE);
//
//
//      for(int ch=0; ch<16; ch++){
//
//         orderedArray[ch] = maxFreqArray[fIndex[ch]];
//         orderedArrayPolType[ch] = maxFreqArrayPolType[fIndex[ch]];
//         //cout<<orderedArray[ch]<<",";
//
//      }
//      //cout<<endl;
//      int cwCount=0;
//      int cwCount_V, cwCount_H, cwCount_X;
//      cwCount_V = cwCount_H = cwCount_X = 0;
//
//      for(int i=0; i<16; i++){
//         for(int j=i+1; j<16; j++){
//
//            //cout<<"i: "<<i<<" j: "<<j<<endl;
//            double fftResGap;
//            if(orderedArrayPolType[i]+orderedArrayPolType[j] == 0){ //2 Vpol
//               //fftRes = 2. * dummyData->freqBinWidth_V;
//               vResBin = int(fftRes / dummyData->freqBinWidth_V)+1;
//               fftResGap = dummyData->freqBinWidth_V * (double)vResBin;
//            }
//            else if(orderedArrayPolType[i]+orderedArrayPolType[j] == 2){ //2H
//               //fftRes = dummyData->freqBinWidth_V + dummyData->freqBinWidth_H;
//               hResBin = int(fftRes / dummyData->freqBinWidth_H)+1;
//               fftResGap = dummyData->freqBinWidth_H * (double)hResBin;
//            }
//            else{ //1V + 1H
//              //fftRes = 2. * dummyData->freqBinWidth_H;
//              //xResBin = int(fftRes / (dummyData->freqBinWidth_V + dummyData->freqBinWidth_H));
//              //fftResGap = (dummyData->freqBinWidth_V + dummyData->freqBinWidth_H) * (double)xResBin + (dummyData->freqBinWidth_V>dummyData->freqBinWidth_H?dummyData->freqBinWidth_V:dummyData->freqBinWidth_H);
//              fftResGap = fftRes + (dummyData->freqBinWidth_V>dummyData->freqBinWidth_H?dummyData->freqBinWidth_V:dummyData->freqBinWidth_H);
//            }
//
//            //cout<<"poltype: "<<orderedArrayPolType[i]+orderedArrayPolType[j]<<" fftResGap: "<<fftResGap<<" orderedArray[i]: "<<orderedArray[i]<<" orderedArray[j]: "<<orderedArray[j]<<" diff: "<<orderedArray[j]-orderedArray[i]<<endl;
//            //printf("fftResGap: %le diff: %le diff-fftResGap: %le\n", fftResGap, orderedArray[j]-orderedArray[i], orderedArray[j]-orderedArray[i]-fftResGap);
//            if(orderedArray[i] > 1e-6 && orderedArray[j] > 1e-6){ //not zeros
//
//            if( (orderedArray[j] - orderedArray[i]) < fftResGap+1e-6/*fftRes*/) {
//               //cout<<"i: "<<i<<" j: "<<j<<" freq_i: "<<orderedArray[i]<<" freq_j: "<<orderedArray[j]<<endl;
//               if(orderedArrayPolType[i]+orderedArrayPolType[j] == 0){ cwCount_V++; }
//               else if (orderedArrayPolType[i]+orderedArrayPolType[j] == 2){ cwCount_H++; }
//               else {cwCount_X++;}
//               cwCount++;
//               i = j;
//            }
//
//            }
//
//         }
//      }
//
//      if(cwCount_V>=2) isVpolCW = true;
//      if(cwCount_H>=2) isHpolCW = true;
//      if(cwCount_X>=2) isXpolCW = true;
//      if(cwCount>=2) isCW = true;//cout<<"CW EVENT!!!!!"<<endl;
//      else isCW=false;//cout<<"NOT CW!!!!!"<<endl;
//




   bool lowFreqDominance = false;
   int lowFreqCountThres = 4;
   int lowFreqCount_V, lowFreqCount_H;
   lowFreqCount_V = lowFreqCount_H = 0;

   lowFreqDominance = isLowFreqDominance(lowFreqCount_V, lowFreqCount_H, dummyData, highPassFreq, lowFreqCountThres);
/*
   for(int i=0; i<8; i++){
      //cout<<"maxFreqBin: "<<dummyData->maxFreqBin[i]<<" maxFreq_V: "<<dummyData->maxFreqBin[i]   * dummyData->freqBinWidth_V<<endl;
      //cout<<"maxFreqBin: "<<dummyData->maxFreqBin[i+8]<<" maxFreq_H: "<<dummyData->maxFreqBin[i+8]   * dummyData->freqBinWidth_H<<endl;
      if( dummyData->maxFreqBin[i]   * dummyData->freqBinWidth_V < highPassFreq ) lowFreqCount_V += 1;
      if( dummyData->maxFreqBin[i+8] * dummyData->freqBinWidth_H < highPassFreq ) lowFreqCount_H += 1;
   }

   if(lowFreqCount_V >= lowFreqCountThres || lowFreqCount_H >= lowFreqCountThres) lowFreqDominance = true;
*/
   //cout<<"lowFreqCount_V: "<<lowFreqCount_V<<" lowFreqCount_H: "<<lowFreqCount_H<<" lowFreqDominance: "<<lowFreqDominance<<endl;

   //isCW = (isVpolCW || isHpolCW);
   //cout<<"isCW: "<<isCW<<endl;
   //passCWCut = (!isCW);

   maxCountFreq_V_hist->Fill(dummyData->maxCountFreq_V, dummyData->weight);
   maxCountFreq_H_hist->Fill(dummyData->maxCountFreq_H, dummyData->weight);

   if(isVpolCW){

      for(int ch=0; ch<8; ch++){
      dFHist_V->Fill( dummyData->freqBinWidth_V * (dummyData->maxFreqBin[ch] - maxCountFreqBin_V) );
      }
   }

   if(isHpolCW){

      for(int ch=8; ch<16; ch++){
      dFHist_H->Fill( dummyData->freqBinWidth_H * (dummyData->maxFreqBin[ch] - maxCountFreqBin_H) );
      }
   }

   /***** 1. Spikey event rejection *********/
/*
   if( !(dummyData->unixTime >= CORRUPT_EVENT_START_TIME && dummyData->eventNumber < CORRUPT_EVENT_END_EVENT_NUMBER) ){
      passCorruptionCut = true;
   }
*/
   /***** 2. Thermal cut ********************/

//passThermalCut = !isThermal_boxCut(inBand, settings, dummyData, onion,  cutValues->snrCut_inBand[type-1].val, cutValues->coherenceCut_inBand[type-1].val, cutValues->snrCut_outOfBand[type-1].val, cutValues->coherenceCut_outOfBand[type-1].val);


/*
   r     = onion.getLayerRadius(dummyData->maxPixIdx2);
   theta = onion.getPointing(dummyData->maxPixIdx2).theta * TMath::RadToDeg();
   phi   = onion.getPointing(dummyData->maxPixIdx2).phi   * TMath::RadToDeg();
*/
   if(dummyData->maxPixCoherence > dummyData->maxPixCoherence2){

      zen_bestHypo = 90.f-dummyData->recoZen;
      azi_bestHypo = dummyData->recoAzi;
      coherence = dummyData->maxPixCoherence;

   } else {

      zen_bestHypo = 90.f-theta;
      azi_bestHypo = phi;
      coherence = dummyData->maxPixCoherence2;

   }

   //if (zen_bestHypo < ZEN_BAND_MAX && zen_bestHypo > ZEN_BAND_MIN){ snrCutValue = cutValues->snrCut_inBand[type-1].val;    coherenceCutValue = cutValues->coherenceCut_inBand[type-1].val; inBand = true; }
   //else                                                           { snrCutValue = cutValues->snrCut_outOfBand[type-1].val; coherenceCutValue = cutValues->coherenceCut_outOfBand[type-1].val; inBand = false;}

   if (zen_bestHypo < ZEN_BAND_MAX && zen_bestHypo > ZEN_BAND_MIN){ coherenceCutValue = cutValues->coherenceCut_inBand[type-1].val; inBand = true; }
   else                                                           { coherenceCutValue = cutValues->coherenceCut_outOfBand[type-1].val; inBand = false;}

   if(string(settings->recoPolType)=="vpol"){ snr = dummyData->inWindowSNR_V; }
   else if(string(settings->recoPolType)=="hpol"){ snr = dummyData->inWindowSNR_H; }

   //if( snr > snrCutValue || coherence > coherenceCutValue){

   //if(inBand){
      //if( 0.003 * snr + 1.916 * coherence - 0.368 > -0.144889) passThermalCut = true;
   //   if(coherence > coherenceCut_inBand[type-1]/*0.108008*/) passThermalCut = true;
   //   else passThermalCut = false;
   //} else{
      //if( 0.007 * snr + 1.039 * coherence - 0.284 > -0.126811) passThermalCut = true;
   //   if(coherence > coherenceCut_outOfBand[type-1]/*0.103715*/) passThermalCut = true;
   //   else passThermalCut = false;

   //}
   if(coherence > coherenceCutValue) passThermalCut = true;
   else passThermalCut = false;
      //passThermalCut = true;
   //}

   snrCutValue = cutValues->snrCut[type-1].val;
   if(snr > snrCutValue){
      passSNRCut = true;
   }


   /***** 3. Surface cut ********************/

   passSurfaceCut = !isSurface(dummyData, surfaceCut_1);

   //if(90.f-dummyData->constantNZen < /*SURFACE_CUT*/surfaceCut_1){
   //   passSurfaceCut = true;
   //}

   float zenRange = 3.;
   double zenMaj;
   passSurfaceCut_2 = !isIterSurface(zenMaj, dummyData, onion, settings, zenRange, surfaceCut_2);
//
//   cout<<"zenMaj in func: "<<zenMaj<<endl;
//   iterZenVec.clear();
//   for(int iter=0; iter<numIter; iter++){
//
//      if(8-5+iter >= 5){ // has >=5 chans in reco
//
//         for(int layer=0; layer<nLayer; layer++){
//
//            iterMaxPixIdxEachLayer[layer] = dummyData->iterMaxPixIdxEachLayer.at(iter*nLayer+layer);
//            iterMaxPixCoherenceEachLayer[layer] = dummyData->iterMaxPixCoherenceEachLayer.at(iter*nLayer+layer);
//
//         }
//
//         TMath::Sort(50, iterMaxPixCoherenceEachLayer, iterIndex);
//         float iterZen = 90.f - onion.getPointing(iterMaxPixIdxEachLayer[iterIndex[0]]).theta * TMath::RadToDeg();
//         cout<<"iterZen: "<<iterZen<<endl;
//         iterZenVec.push_back(iterZen);
//
//      }//end of if
//   }//end of iter
//
//   zenRange = 3.;
//   zenMaj = getZenMaj(iterZenVec, zenRange);
//   cout<<"zenMaj: "<<zenMaj<<endl;
//
//   if(zenMaj <= 90 ){
//      //iterMajorityZenHist->Fill(zenMaj,dummyData->weight);
//
//      if(zenMaj < /*SURFACE_CUT_2*/surfaceCut_2 ){
//         passSurfaceCut_2 = true;
//      }
//
//   } else passSurfaceCut_2 = true; //if no majority zenith can be found through iter reco, use solely the constantN zenith to check whether passed surface cut or not
//   //{  zenMaj = 90.f - onion.getPointing(iterMaxPixIdxEachLayer[3]).theta * TMath::RadToDeg(); }

   /***** 4. Deep pulser time cut ***********/
   passDeepPulserCut = !isDeepPulser(STATION, dummyData, runNum);
   /*
   if( STATION == "ARA02" ){
      if( !(dummyData->unixTime < 1420.5122e6 && dummyData->unixTime > 1420.50905e6) && !(runNum >= 4795 && runNum <= 4800) && !(runNum == 4787 || runNum==4785 ) ){
         passDeepPulserCut = true;
      }
   }
   */

   /***** Calpulser sweep time cut **********/

   passCalpulserTimeCut = !isCalpulserTime(STATION, dummyData);
/*
   if( STATION == "ARA02" ){
      if( !(dummyData->unixTime < 1393923046 && dummyData->unixTime > 1393917793) && !(dummyData->unixTime < 1395649842 && dummyData->unixTime > 1395648365) ){
         passCalpulserTimeCut = true;
      }
   } else if (STATION == "ARA03" ){
      if( !(dummyData->unixTime < 1393923742 && dummyData->unixTime > 1393922266) && !(dummyData->unixTime < 1395650418 && dummyData->unixTime > 1395648942) ){
         passCalpulserTimeCut = true;
      }
   }
*/
   /***** Calpulser angular cut *************/
   /*
   const int nBoxes = 4;
   float zenMin[nBoxes], zenMax[nBoxes], aziMin[nBoxes], aziMax[nBoxes];
   //ARA02 type 1-5 group 1 calpulser cut values - D6BV
   zenMin[0] = 0.5409;
   zenMax[0] = 10.99;
   aziMin[0] = 54.96;
   aziMax[0] = 66.74;

   ////ARA02 type 1-5 group 2 calpulser cut values - D5BV
   zenMin[1] = -31.40;
   zenMax[1] = -14.26;
   aziMin[1] = 324.05;
   aziMax[1] = 341.01;

   ////ARA02 type 3 group 3 calpulser cut values
   zenMin[2] = -27.78;
   zenMax[2] = -14.22;
   aziMin[2] = 50.80;
   aziMax[2] = 58.32;

   ////ARA02 type 3 group 4 calpulser cut values
   zenMin[3] = 24.07;
   zenMax[3] = 31.45;
   aziMin[3] = 67.30;
   aziMax[3] = 75.25;
   */

   float inBoxTheta, inBoxPhi;
   inBoxTheta = inBoxPhi = 0.f;
   passCalpulserCut = !isCalpulser(inBoxTheta, inBoxPhi, STATION, dummyData, onion, settings, type);

//   bool inBox = false;
//
//   bool iterInBox = false;
//   int inBoxCount = 0;
//   float inBoxTheta = 0.f;
//   float inBoxPhi   = 0.f;
//   for(int iter=0; iter<numIter; iter++){
//
//      //int maxPixIdx = dummyData->iterMaxPixIdx.at(iter);
//      //float maxPixCoherence = dummyData->iterMaxPixCoherence.at(iter);
//      int maxPixIdx = dummyData->iterMaxPixIdxEachLayer.at(iter*nLayer+0);
//      float maxPixCoherence = dummyData->iterMaxPixCoherenceEachLayer.at(iter*nLayer+0);
//      //cout<<"maxPixIdx: "<<maxPixIdx<<" layer: "<<maxPixIdx/nDir<<" maxPixCoherence: "<<maxPixCoherence<<endl;
//      float theta = 90.f-TMath::RadToDeg()*onion.getPointing(maxPixIdx).theta;
//      float phi   = TMath::RadToDeg()*onion.getPointing(maxPixIdx).phi;
//
//
//      iterInBox = false;
//      /*** Box 1 ***/
//      ///if( theta > zenMin[0] && theta < zenMax[0] && phi > aziMin[0] && phi < aziMax[0] ) { inBox = true; iterInBox = true;}
//
//      /*** Box 2 ***/
//      //if( type != 3){
//      //if( (theta > zenMin[1] && theta < zenMax[1] && phi > aziMin[1] && phi < aziMax[1]) ){ inBox = true; iterInBox = true;}
//      //}
//
//      for(int box=0; box<cutValues->nBoxes; box++){
//
//      // 0: D5BV, 1: D6BV, 2: D5BV Mirror, only for type 5
//         if(box<2){
//
//            if( theta > cutValues->zenMin[box].val && theta < cutValues->zenMax[box].val && phi > cutValues->aziMin[box].val && phi < cutValues->aziMax[box].val ) { inBox = true; iterInBox = true;}
//
//         } else {
//
//            if( type == 5 ){
//                if( theta > cutValues->zenMin[box].val && theta < cutValues->zenMax[box].val && phi > cutValues->aziMin[box].val && phi < cutValues->aziMax[box].val ) { inBox = true; iterInBox = true;}
//            }
//         }
//      }
//
//      /*** Box 3 & 4 ***/
//      //if ( type == 3 ){
//      //   if( (theta > zenMin[2] && theta < zenMax[2] && phi > aziMin[2] && phi < aziMax[2])
//      //     ||(theta > zenMin[3] && theta < zenMax[3] && phi > aziMin[3] && phi < aziMax[3])){ inBox = true; iterInBox = true;}
//      //}
//
//      if(iterInBox){
//         inBoxCount++;
//         inBoxTheta += theta;
//         inBoxPhi   += phi;
//      }
//
//   }//end of iter
//
//   if (!inBox) passCalpulserCut = true;
//
//   if(inBoxCount > 0 ){
//
//      inBoxTheta /= (float)inBoxCount;
//      inBoxPhi   /= (float)inBoxCount;
//   }
//
   /******* Noisy run cut ****************/

   /*
   for(int run=0; run<numNoisyRuns; run++){

      if( dummyData->unixTime>noisyRun[run][1] && dummyData->unixTime<noisyRun[run][2] ){
         passNoisyRunCut = false;
         break;
      }
   }
   */
   int plusMinusRunNum = 0; //an event is considered in noisy period if its run number is +-1 run of a known (listed) noisy run
   passNoisyRunCut = !(isNearNoisyRun(listOfRuns, runNum, plusMinusRunNum));
   //if(runNum >= 4795 && runNum <= 4800) passNoisyRunCut = false; //DP
   //if(runNum >= 3 && runNum <=60 && runNum != 50) passNoisyRunCut = false; //Corrupted wf
   //if(runNum == 4787 || runNum==4785 ) passNoisyRunCut = false; //DP


   /* Check if CW-tagged event can be recovered by surviving the impulsivity cut */
   /*
   if ( isCW ){

      isCWCount += dummyData->weight;
      double impCut = impulsivityCut[type-1];

      std::fill(&impulsivity[0], &impulsivity[16], 0.);

      if((isVpolCW && isHpolCW) || isXpolCW){

      nonZeroCount = 0;
      double sum = 0.;
      for(int ch=0; ch<16; ch++){
         if(fabs( dummyData->impulsivity[ch] - 0 ) > 1e-9 ){
            nonZeroCount++;
            impulsivity[ch] = dummyData->impulsivity[ch];
            sum += impulsivity[ch];
         }
      }

      //int index[16];
      //TMath::Sort(16, impulsivity, index);
      avgImpulsivity = sum / (double)nonZeroCount;

      //impulsivityHist_max->Fill(impulsivity[index[0]], dummyData->weight);
      //impulsivityHist_3rd->Fill(impulsivity[index[2]], dummyData->weight);
      //impulsivityHist_avg->Fill(avgImpulsivity, dummyData->weight);
      //outputFile<<avgImpulsivity<<",";
      //thermalCWEventCount_both += dummyData->weight;

      if(dummyData->maxCountFreq_V > highPassFreq || dummyData->maxCountFreq_H > highPassFreq) passHighPassFilter = true;

      }

      else if( isVpolCW && !isHpolCW){

      nonZeroCount = 0;
      double sum = 0.;
      for(int ch=0; ch<8; ch++){
         if(fabs( dummyData->impulsivity[ch] - 0 ) > 1e-9 ){
            nonZeroCount++;
            impulsivity[ch] = dummyData->impulsivity[ch];
            sum += impulsivity[ch];
         }
      }

      //int index[16];
      //TMath::Sort(16, impulsivity, index);
      avgImpulsivity = sum / (double)nonZeroCount;

      //impulsivityHist_max->Fill(impulsivity[index[0]], dummyData->weight);
      //impulsivityHist_3rd->Fill(impulsivity[index[2]], dummyData->weight);
      //impulsivityHist_avg->Fill(avgImpulsivity, dummyData->weight);
      //outputFile<<avgImpulsivity<<",";
      //thermalCWEventCount_V += dummyData->weight;

      if(dummyData->maxCountFreq_V > highPassFreq ) passHighPassFilter = true;

      }

      else {

      nonZeroCount = 0;
      double sum = 0.;
      for(int ch=8; ch<16; ch++){
         if(fabs( dummyData->impulsivity[ch] - 0 ) > 1e-9 ){
            nonZeroCount++;
            impulsivity[ch] = dummyData->impulsivity[ch];
            sum += impulsivity[ch];
         }
      }

      //int index[16];
      //TMath::Sort(16, impulsivity, index);
      avgImpulsivity = sum / (double)nonZeroCount;

      //impulsivityHist_max->Fill(impulsivity[index[0]], dummyData->weight);
      //impulsivityHist_3rd->Fill(impulsivity[index[2]], dummyData->weight);
      //impulsivityHist_avg->Fill(avgImpulsivity, dummyData->weight);
      //outputFile<<avgImpulsivity<<",";
      //thermalCWEventCount_H += dummyData->weight;

      if(dummyData->maxCountFreq_H > highPassFreq ) passHighPassFilter = true;

      }
      //cout<<"maxCountFreq_V: "<<dummyData->maxCountFreq_V<<" maxCountFreq_H: "<<dummyData->maxCountFreq_H<<" avgImp: "<<avgImpulsivity<<" impCut: "<<impCut<<endl;
      if(avgImpulsivity > impCut){
         passImpulsivityCut = true;
         nRecoveredByImp += dummyData->weight;
      }




   }
   */
   double impCut = cutValues->cwImpCut[type-1].val; //impulsivityCut[type-1];
   impCut = 0.2579306/*0.2384656*/;
   //passCWCut = ( !isCW || (isCW && passHighPassFilter && passImpulsivityCut )) && !lowFreqDominance;
   //passCWCut = ( !isCW || (isCW && isRecoverableByImp(isVpolCW, isHpolCW, isXpolCW, dummyData, impCut, highPassFreq) )) && !lowFreqDominance;
   passCWCut = !lowFreqDominance;


   /****** Check if pass thermal impulsivity cut ********/
/*
      std::fill(&impulsivity[0], &impulsivity[16], 0.);
      nonZeroCount = 0;
      double sum = 0.;
      for(int ch=0; ch<8; ch++){
         if(fabs( dummyData->impulsivity[ch] - 0 ) > 1e-9 ){
            nonZeroCount++;
            impulsivity[ch] = dummyData->impulsivity[ch];
            sum += impulsivity[ch];
         }
      }

      int index[16];
      TMath::Sort(16, impulsivity, index);
      avgImpulsivity = sum / (double)nonZeroCount;
      if(avgImpulsivity > postThermalAvgImpulsivityCut){

         //outputFile<<runNum<<","<<dummyData->eventNumber<<","<<dummyData->unixTime<<endl;
         passThermalImpulsivityCut = true;


      }
*/
   postThermalAvgImpulsivityCut = 0.261729;
   passThermalImpulsivityCut = !isBelowThermalImpulsivityCut(avgImpulsivity, dummyData, postThermalAvgImpulsivityCut);
   passThermalImpulsivityCut = true;
   /***Check if have enough un-saturated channels to reconstruct ****/
   //if(dummyData->numSatChan  <= 3){ //Need at least 5 channels in reconstruction
     passNumSatChanCut = true;
   //}

   //passThermalImpulsivityCut=1;
   //passNoisyRunCut = 1;
   //nPassNumSatChanCut += passNumSatChanCut * dummyData->weight;
   //nCut0p5            += (passNumSatChanCut) * dummyData->weight;
   //nPassCorruption    += passCorruptionCut * dummyData->weight;
   //nCut1              += (passCorruptionCut) * dummyData->weight;
   //nCut0              += (passImpulsivityCut) * dummyData->weight;
   //nPassImpulsivityCut+= passImpulsivityCut *dummyData->weight;
   //nCut1              += (passHighPassFilter && passImpulsivityCut) * dummyData->weight;
   //nPassHighPassFilter+= passHighPassFilter * dummyData->weight;
   snrHist[0]->Fill(dummyData->inWindowSNR, dummyData->weight);
   nPassCWCut         += passCWCut * dummyData->weight;
   nCut1p5            += (passNumSatChanCut && passCWCut) * dummyData->weight;
   if (passNumSatChanCut && passCWCut) snrHist[1]->Fill(dummyData->inWindowSNR, dummyData->weight);
   nPassDeepPulserCut += passDeepPulserCut * dummyData->weight;
   nCut2              += (passNumSatChanCut &&/*passHighPassFilter && passImpulsivityCut*/passCWCut &&/*passCorruptionCut && passThermalCut && passSuE19aceCut &&*/ passDeepPulserCut) * dummyData->weight;
   nPassThermalCut    += passThermalCut * dummyData->weight;
   nCut3              += (passNumSatChanCut &&/*passHighPassFilter && passImpulsivityCut*/passCWCut &&/*passCorruptionCut &&*/ passDeepPulserCut && passThermalCut) * dummyData->weight;
   if (passNumSatChanCut &&/*passHighPassFilter && passImpulsivityCut*/passCWCut &&/*passCorruptionCut &&*/ passDeepPulserCut && passThermalCut) snrHist[2]->Fill(dummyData->inWindowSNR, dummyData->weight);
   nPassThermalImpulsivityCut += passThermalImpulsivityCut * dummyData->weight;
   nPassSNRCut += passSNRCut * dummyData->weight;
   nCut3p5            += (passNumSatChanCut &&/*passHighPassFilter && passImpulsivityCut*/passCWCut &&/*passCorruptionCut &&*/ passDeepPulserCut && passThermalCut && passThermalImpulsivityCut && passSNRCut) * dummyData->weight;
   //nPassSuE19aceCut    += passSuE19aceCut * dummyData->weight;
   //nCut3              += (passCorruptionCut && passThermalCut && passSuE19aceCut) * dummyData->weight;
   //nPassDeepPulserCut += passDeepPulserCut * dummyData->weight;
   //nCut4              += (passCorruptionCut && passThermalCut /*&& passSuE19aceCut*/ && passDeepPulserCut) * dummyData->weight;
   //if((passCorruptionCut && passThermalCut /*&& passSuE19aceCut*/ && passDeepPulserCut)){ outputFile<<dummyData->unixTime<<","<<dummyData->eventNumber<<endl; }
   nPassCalpulserCut += passCalpulserCut * dummyData->weight;
   nCut4             += (passNumSatChanCut &&/*passHighPassFilter && passImpulsivityCut*/passCWCut && passDeepPulserCut && passThermalCut && passThermalImpulsivityCut && passSNRCut && passCalpulserCut) * dummyData->weight;
   if (passNumSatChanCut &&/*passHighPassFilter && passImpulsivityCut*/passCWCut && passDeepPulserCut && passThermalCut && passThermalImpulsivityCut && passSNRCut && passCalpulserCut)  snrHist[3]->Fill(dummyData->inWindowSNR, dummyData->weight);
   nPassCalpulserTimeCut += passCalpulserTimeCut * dummyData->weight;
   nCut4p5           += (passNumSatChanCut &&/*passHighPassFilter && passImpulsivityCut*/passCWCut && passDeepPulserCut && passThermalCut && passThermalImpulsivityCut && passSNRCut && passCalpulserCut && passCalpulserTimeCut) * dummyData->weight;
   //nPassNoisyRunCut   += passNoisyRunCut * dummyData->weight;
   //nCut5              += (/*passCorruptionCut &&*/ passThermalCut /*&& passSuE19aceCut*/ && passDeepPulserCut && passCalpulserCut && passNoisyRunCut) * dummyData->weight;
   nPassSurfaceCut    += passSurfaceCut * dummyData->weight;
   nCut6              += (passNumSatChanCut &&/*passHighPassFilter && passImpulsivityCut*/passCWCut &&/*passCorruptionCut &&*/ passThermalCut && passThermalImpulsivityCut && passSNRCut && passDeepPulserCut && passCalpulserCut && passCalpulserTimeCut/*&& passNoisyRunCut*/ && passSurfaceCut) * dummyData->weight;
   nPassSurfaceCut_2  += passSurfaceCut_2 * dummyData->weight;
   nCut6p5            += (passNumSatChanCut &&/*passHighPassFilter && passImpulsivityCut*/passCWCut &&/*passCorruptionCut &&*/ passThermalCut && passThermalImpulsivityCut && passSNRCut && passDeepPulserCut && passCalpulserCut && passCalpulserTimeCut/*&& passNoisyRunCut*/ && passSurfaceCut && passSurfaceCut_2) * dummyData->weight;
   if (passNumSatChanCut &&/*passHighPassFilter && passImpulsivityCut*/passCWCut &&/*passCorruptionCut &&*/ passThermalCut && passThermalImpulsivityCut && passSNRCut && passDeepPulserCut && passCalpulserCut && passCalpulserTimeCut/*&& passNoisyRu     nCut*/ && passSurfaceCut && passSurfaceCut_2) snrHist[4]->Fill(dummyData->inWindowSNR, dummyData->weight);
   nPassNoisyRunCut   += passNoisyRunCut * dummyData->weight;
   nCut7              +=  (passNumSatChanCut &&/*passHighPassFilter && passImpulsivityCut*/passCWCut &&/*passCorruptionCut &&*/ passThermalCut && passThermalImpulsivityCut && passSNRCut && passDeepPulserCut && passCalpulserCut && passCalpulserTimeCut && passNoisyRunCut && passSurfaceCut && passSurfaceCut_2) * dummyData->weight;

//   if(/*passNumSatChanCut &&*//*passHighPassFilter && passImpulsivityCut*/passCWCut && passThermalCut && passThermalImpulsivityCut && passDeepPulserCut && passCalpulserCut && passCalpulserTimeCut/*&& passNoisyRunCut*/ && passSurfaceCut && passSurfaceCut_2 && passNoisyRunCut){
//   //if(passCalpulserCut){
//
//      cout<<endl;
//      cout<<"run: "<<runNum<<" event: "<<dummyData->eventNumber<<" unixtime: "<<dummyData->unixTime<<"snr: "<<dummyData->inWindowSNR_V<<" coherence: "<<(dummyData->maxPixCoherence>dummyData->maxPixCoherence2?dummyData->maxPixCoherence:dummyData->maxPixCoherence2)<<endl;
//      cout<<"inBand: "<<inBand<<" Fisher: "<<(inBand?0.003:0.007)*snr+(inBand?1.916:1.039)*coherence+(inBand?-0.368:-0.284)<<endl;
//      cout<<"constantNZen: "<<90.f-dummyData->constantNZen<<" zenMaj: "<<zenMaj<<endl;
//      cout<<"avg imp: "<<avgImpulsivity<<endl;
//      outputFile<<runNum<<","<<dummyData->eventNumber<<","<<dummyData->unixTime<<endl;
//      //outputFile<<avgImpulsivity<<",";
//      //cout<<avgImpulsivity<<endl;
//      cout<<endl;
//      cout<<"maxCountFreq_V: "<<dummyData->maxCountFreq_V<<" maxCountFreq_H: "<<dummyData->maxCountFreq_H<<endl;
//      cout<<"From recoData channel-by-channel maxFreqBin: "<<endl;
//      cout<<"maxCountFreq_V: "<<maxCountFreqBin_V*dummyData->freqBinWidth_V<<" maxCountFreq_H: "<<maxCountFreqBin_H*dummyData->freqBinWidth_H<<endl;
//      cout<<"freqBinWidth_V: "<<dummyData->freqBinWidth_V<<" freqBinWidth_H: "<<dummyData->freqBinWidth_H<<endl;
//
//      cout<<"maxFreqBin: "<<endl;
//      for(int ch=0; ch<16; ch++) cout<<dummyData->maxFreqBin[ch]<<",";
//      cout<<endl<<"maxFreq: "<<endl;
//      for(int ch=0; ch<8; ch++)  cout<<dummyData->maxFreqBin[ch]*dummyData->freqBinWidth_V<<",";
//      cout<<endl;
//      for(int ch=8; ch<16; ch++) cout<<dummyData->maxFreqBin[ch]*dummyData->freqBinWidth_H<<",";
//      cout<<endl;
//
//      cout<<"recoZen: "<<90.f-dummyData->recoZen<<" azi: "<<dummyData->recoAzi<<endl;
      /*
      for(int iter=0; iter<numIter; iter++){
         int maxPixIdx = dummyData->iterMaxPixIdx.at(iter);
         cout<<"iterMaxPixIdx: "<<maxPixIdx<<" layer: "<<maxPixIdx/nDir<<endl;
         cout<<"zen: "<<90.f-TMath::RadToDeg()*onion.getPointing(maxPixIdx).theta<<" azi: "<<TMath::RadToDeg()*onion.getPointing(maxPixIdx).phi<<endl;
      }
      */
      //outputFile<<"run: "<<runNum<<" event: "<<dummyData->eventNumber<<" unixtime: "<<dummyData->unixTime<<endl;
      //cout<<"isVpolCW: "<<isVpolCW<<" isHpolCW: "<<isHpolCW<<endl;
      //cout<<"maxCountFreq_V: "<<dummyData->maxCountFreq_V<<" maxCountFreq_H: "<<dummyData->maxCountFreq_H<<endl;
//      impulsivityHist_avg->Fill(avgImpulsivity, dummyData->weight);
//      //impulsivityHist_3rd->Fill(impulsivity[index[2]], dummyData->weight);
//      c_vs_snr_hist->Fill(dummyData->inWindowSNR_V, coherence, dummyData->weight);
//
//      for(int ch=0; ch<8; ch++){
//
//         double bipolarRatio;
//         double pos, neg;
//         pos = dummyData->posPowerPeak[ch];
//         neg = fabs(dummyData->negPowerPeak[ch]);
//
//         if(pos > neg) bipolarRatio = neg/pos;
//         else{ if(neg>0) bipolarRatio = pos/neg; else bipolarRatio = 0.; }
//
//         double dT = fabs(dummyData->powerPeaksDeltaT[ch]);
//         if(dummyData->slidingV2SNR[ch]>settings->nchnlThreshold){
//         bipolarRatio_dT->Fill(dT, bipolarRatio, dummyData->weight);
//         impulsivity_dT->Fill(dT, dummyData->impulsivity[ch], dummyData->weight);
//         impulsivity_bipolarRatio->Fill(bipolarRatio, dummyData->impulsivity[ch], dummyData->weight);
//         }
//      }

//      for(int ch=0; ch<16; ch++){
//
//         //maxFreqBinVec.push_back(dummyData->maxFreqBin[ch]);
//         maxFreqArray[ch] = dummyData->maxFreqBin[ch] * (ch<8?dummyData->freqBinWidth_V:dummyData->freqBinWidth_H);
//         if(ch<8) maxFreqArrayPolType[ch] = 0;//maxFreqArray_V[ch] =  dummyData->maxFreqBin[ch] * dummyData->freqBinWidth_V;
//         else     maxFreqArrayPolType[ch] = 1;//maxFreqArray_H[ch] =  dummyData->maxFreqBin[ch] * dummyData->freqBinWidth_H;
//
//      }
//
//
//
//      TMath::Sort(16, maxFreqArray, fIndex, kFALSE);
//      //TMath::Sort(8, maxFreqArray_V, fIndex_V, kFALSE);
//      //TMath::Sort(8, maxFreqArray_H, fIndex_H, kFALSE);
//
//
//      for(int ch=0; ch<16; ch++){
//
//         orderedArray[ch] = maxFreqArray[fIndex[ch]];
//         orderedArrayPolType[ch] = maxFreqArrayPolType[fIndex[ch]];
//         cout<<orderedArray[ch]<<",";
//
//      }
//      cout<<endl;
///*
//      for(int ch=0; ch<8; ch++){
//
//         orderedArray[ch] =
//
//      }
//*/
//      int cwCount=0;
//      for(int i=0; i<16; i++){
//         for(int j=i+1; j<16; j++){
//
//            //cout<<"i: "<<i<<" j: "<<j<<endl;
//            double fftResGap;
//            if(orderedArrayPolType[i]+orderedArrayPolType[j] == 0){ //2 Vpol
//               //fftRes = 2. * dummyData->freqBinWidth_V;
//               vResBin = int(fftRes / dummyData->freqBinWidth_V)+1;
//               fftResGap = dummyData->freqBinWidth_V * (double)vResBin;
//            }
//            else if(orderedArrayPolType[i]+orderedArrayPolType[j] == 2){ //2H
//               //fftRes = dummyData->freqBinWidth_V + dummyData->freqBinWidth_H;
//               hResBin = int(fftRes / dummyData->freqBinWidth_H)+1;
//               fftResGap = dummyData->freqBinWidth_H * (double)hResBin;
//            }
//            else{ //1V + 1H
//              //fftRes = 2. * dummyData->freqBinWidth_H;
//              //xResBin = int(fftRes / (dummyData->freqBinWidth_V + dummyData->freqBinWidth_H));
//              //cout<<"xResBin: "<<xResBin<<" larger binWidth: "<<(dummyData->freqBinWidth_V>dummyData->freqBinWidth_H?dummyData->freqBinWidth_V:dummyData->freqBinWidth_H)<<endl;
//              //cout<<"xResBin*(V+H): "<<(dummyData->freqBinWidth_V + dummyData->freqBinWidth_H) * (double)xResBin<<endl;
//              //fftResGap = (dummyData->freqBinWidth_V + dummyData->freqBinWidth_H) * (double)xResBin + (dummyData->freqBinWidth_V>dummyData->freqBinWidth_H?dummyData->freqBinWidth_V:dummyData->freqBinWidth_H);
//              fftResGap = fftRes + (dummyData->freqBinWidth_V>dummyData->freqBinWidth_H?dummyData->freqBinWidth_V:dummyData->freqBinWidth_H);
//            }
//
//            //cout<<"poltype: "<<orderedArrayPolType[i]+orderedArrayPolType[j]<<" fftResGap: "<<fftResGap<<" orderedArray[i]: "<<orderedArray[i]<<" orderedArray[j]: "<<orderedArray[j]<<" diff: "<<orderedArray[j]-orderedArray[i]<<endl;
//            //printf("fftResGap: %le diff: %le diff-fftResGap: %le\n", fftResGap, orderedArray[j]-orderedArray[i], orderedArray[j]-orderedArray[i]-fftResGap);
//            if(orderedArray[i] > 1e-6 && orderedArray[j] > 1e-6){ //not zeros
//
//            if( (orderedArray[j] - orderedArray[i]) < /*fftResGap+1e-6*/fftRes) {
//               cout<<"i: "<<i<<" j: "<<j<<" freq_i: "<<orderedArray[i]<<" freq_j: "<<orderedArray[j]<<endl;
//               cwCount++;
//               i = j;
//            }
//
//            }
//
//         }
//      }
//
//      if(cwCount>=2) cout<<"CW EVENT!!!!!"<<endl;
//      else cout<<"NOT CW!!!!!"<<endl;
//
//   }
//
//   //c_vs_imp->Fill(avgImpulsivity, coherence, dummyData->weight);
//
   //if(!lowFreqDominance &&/*passCWCut &&*/passThermalCut && passThermalImpulsivityCut && passDeepPulserCut && passCalpulserCut && passCalpulserTimeCut && passSurfaceCut && passSurfaceCut_2 && passNoisyRunCut ) impulsivityHist_nMinusCW->Fill( avgImpulsivity, dummyData->weight);

   if(passCWCut && /*passThermalCut && */ passSNRCut && passThermalImpulsivityCut && passDeepPulserCut && passCalpulserCut && passCalpulserTimeCut && passSurfaceCut && passSurfaceCut_2 && passNoisyRunCut ) //c_vs_snr_hist_nMinusThermal->Fill(snr, coherence, dummyData->weight);
      {
         coherence_nMinusThermal->Fill(coherence, dummyData->weight);
         if(inBand) outputFile<<coherence<<",";
      }

   if(passCWCut && passThermalCut &&/* passThermalImpulsivityCut &&*/ passDeepPulserCut && passCalpulserCut && passCalpulserTimeCut && passSurfaceCut && passSurfaceCut_2 && passNoisyRunCut ) //impulsivityHist_nMinusImp->Fill( avgImpulsivity, dummyData->weight);
   {
      snr_nMinusSNR->Fill(snr, dummyData->weight);
      cout<<"snr: "<<snr<<endl;
   }
   if(passCWCut && passThermalCut && passThermalImpulsivityCut && passSNRCut && passDeepPulserCut /*&& passCalpulserCut*/ && passCalpulserTimeCut && passSurfaceCut && passSurfaceCut_2 && passNoisyRunCut ) zen_azi_nMinusCal->Fill(inBoxPhi, inBoxTheta, dummyData->weight);
   if(passCWCut && passThermalCut && passThermalImpulsivityCut && passSNRCut && passDeepPulserCut && passCalpulserCut && passCalpulserTimeCut /*&& passSurfaceCut && passSurfaceCut_2 */&& passNoisyRunCut ) zen_nMinusSurface->Fill((passSurfaceCut?zenMaj:90.f-dummyData->constantNZen));

   //if(passDeepPulserCut && passThermalCut && passCalpulserCut && passSurfaceCut) runHist->Fill(runNum);
   if( /*passNumSatChanCut && */passCWCut /*!lowFreqDominance*/ && passDeepPulserCut && passThermalCut && /*passThermalImpulsivityCut &&*/ passCalpulserCut && passCalpulserTimeCut && (!passSurfaceCut || !passSurfaceCut_2)) surfaceRunHist->Fill(runNum);
   //if( /*passNumSatChanCut &&*/ /*passCWCut &&*/ passDeepPulserCut && passThermalCut && passSurfaceCut && passSurfaceCut_2){ impulsivityHist_avg->Fill(avgImpulsivity, dummyData->weight); outputFile<<avgImpulsivity<<","; }
//   //if( !lowFreqDominance && passDeepPulserCut){
//   if( passCWCut && passDeepPulserCut ){
//
//
//      if(inBand) inBand_all += dummyData->weight;
//      else       outOfBand_all += dummyData->weight;
//
//      if(passThermalCut){
//         //impulsivityHist_avg->Fill(avgImpulsivity, dummyData->weight);
//         //outputFile<<avgImpulsivity<<",";
//         if(inBand) inBand_pass += dummyData->weight;
//         else       outOfBand_pass += dummyData->weight;
//
//         //constantNZenHist->Fill(90.f-dummyData->constantNZen, dummyData->weight);
//         //outputFile<<90.f-dummyData->constantNZen<<",";
//
//      }
//   }
//
//   if( /*passCWCut &&*/ passDeepPulserCut && passThermalCut && passCalpulserCut && passCalpulserTimeCut){
//
//      constantNZenHist->Fill(90.f-dummyData->constantNZen, dummyData->weight);
//      //outputFile<<90.f-dummyData->constantNZen<<",";
//      /*
//      for(int iter=0; iter<numIter; iter++){
//
//         if(8-5+iter >= 5){
//
//            for(int layer=0; layer<nLayer; layer++){
//
//               iterMaxPixIdxEachLayer[layer] = dummyData->iterMaxPixIdxEachLayer.at(iter*nLayer+layer);
//               iterMaxPixCoherenceEachLayer[layer] = dummyData->iterMaxPixCoherenceEachLayer.at(iter*nLayer+layer);
//
//            }
//
//            TMath::Sort(50, iterMaxPixCoherenceEachLayer, iterIndex);
//            float iterZen = 90.f - onion.getPointing(iterMaxPixIdxEachLayer[iterIndex[0]]).theta * TMath::RadToDeg();
//            iterZenVec.push_back(iterZen);
//
//         }//end of if
//      }//end of iter
//
//      float zenRange = 3.;
//      double zenMaj = getZenMaj(iterZenVec, zenRange);
//      */
//      if(zenMaj <= 90 ) iterMajorityZenHist->Fill(zenMaj,dummyData->weight);
//
//   }//end of pass cut
//
//   if(isCW){
//      coherence_snr_cw->Fill(dummyData->inWindowSNR_V, (dummyData->maxPixCoherence>dummyData->maxPixCoherence2?dummyData->maxPixCoherence:dummyData->maxPixCoherence2), dummyData->weight);
//   }
//
//   int thetaXingPix, phiXingPix, avgThetaXingPix, avgPhiXingPix;
//   double  angThres = 1.;
//   if(passThermalCut && passDeepPulserCut && passSurfaceCut && passSurfaceCut_2 && passCalpulserCut && passCalpulserTimeCut && passNoisyRunCut /*&& !lowFreqDominance */&& passCWCut){
//
//      cout<<"RunNum: "<<runNum<<" event: "<<dummyData->eventNumber<<endl;
//
//      //getAngXingPixels(thetaXingPix, phiXingPix, dummyData, settings, onion, angThres);
//      //getAvgAngXingPixels(avgThetaXingPix, avgPhiXingPix, dummyData, settings, onion, angThres);
//
//      //if(avgPhiXingPix>500){
//      //   cout<<"avgPhiXingPix: "<<avgPhiXingPix<<" runNum: "<<runNum<<" event: "<<dummyData->eventNumber<<endl;
//      //}
//
//      //thetaXingHist->Fill(thetaXingPix, dummyData->weight);
//      //phiXingHist->Fill(phiXingPix, dummyData->weight);
//      //avgThetaXingHist->Fill(avgThetaXingPix, dummyData->weight);
//      //avgPhiXingHist->Fill(avgPhiXingPix, dummyData->weight);
//
//      angThres = 0.5;
//      double inRangeThetaFrac = getZenithInRangeFraction(dummyData, settings, onion, angThres);
//      double inRangePhiFrac   = getAzimuthInRangeFraction(dummyData, settings, onion, angThres);
//      inRangeThetaFracHist->Fill(inRangeThetaFrac, dummyData->weight);
//      inRangePhiFracHist->Fill(inRangePhiFrac, dummyData->weight);
//      inRangeThetaPhiFracHist->Fill(inRangePhiFrac, inRangeThetaFrac, dummyData->weight);
//
//   }
//

//   if(/*isCW &&*/ passThermalCut && passDeepPulserCut && passCalpulserCut && passCalpulserTimeCut && passSurfaceCut && passSurfaceCut_2 && passNoisyRunCut && !lowFreqDominance ){

//      cout<<"run: "<<runNum<<" event: "<<dummyData->eventNumber<<" SNR: "<<snr<<endl;
//      _snrHist->Fill(snr, dummyData->weight);
//      outputFile<<snr<<",";
//
//      std::fill(&impulsivity[0], &impulsivity[16], 0.);
//      int nonZeroCount = 0;
//      double sum = 0.;
//      for(int ch=0; ch<8; ch++){
//         if(fabs( dummyData->impulsivity[ch] - 0 ) > 1e-9 ){
//            nonZeroCount++;
//            impulsivity[ch] = dummyData->impulsivity[ch];
//            sum += impulsivity[ch];
//         }
//      }
//
//      //int index[16];
//      //TMath::Sort(16, impulsivity, index);
//      avgImpulsivity = sum / (double)nonZeroCount;
//
//      impulsivityHist_nMinusCW->Fill(avgImpulsivity, dummyData->weight);
//
//      if(avgImpulsivity > 0.4 ){
//               cout<<endl;
//               cout<<"run: "<<runNum<<" event: "<<dummyData->eventNumber<<" unixtime: "<<dummyData->unixTime<<"snr: "<<snr<<" coherence: "<<coherence<<endl;
//               cout<<"inBand: "<<inBand<<" Fisher: "<<(inBand?0.003:0.007)*snr+(inBand?1.916:1.039)*coherence+(inBand?-0.368:-0.284)<<endl;
//               cout<<"constantNZen: "<<90.f-dummyData->constantNZen<<" zenMaj: "<<zenMaj<<endl;
//               cout<<"avg imp: "<<avgImpulsivity<<endl;
//               outputFile<<runNum<<","<<dummyData->eventNumber<<","<<dummyData->unixTime<<endl;
//               //outputFile<<avgImpulsivity<<",";
//               //cout<<avgImpulsivity<<endl;
//               cout<<endl;
//               cout<<"maxCountFreq_V: "<<dummyData->maxCountFreq_V<<" maxCountFreq_H: "<<dummyData->maxCountFreq_H<<endl;
//               cout<<"From recoData channel-by-channel maxFreqBin: "<<endl;
//               cout<<"maxCountFreq_V: "<<maxCountFreqBin_V*dummyData->freqBinWidth_V<<" maxCountFreq_H: "<<maxCountFreqBin_H*dummyData->freqBinWidth_H<<endl;
//               cout<<"freqBinWidth_V: "<<dummyData->freqBinWidth_V<<" freqBinWidth_H: "<<dummyData->freqBinWidth_H<<endl;
//
//      }

      //if( !isCW ) cerr<<"Not CW! run: "<<runNum<<" ev: "<<entry<<endl;

      //coherence_snr_cw->Fill(dummyData->inWindowSNR_V, (dummyData->maxPixCoherence>dummyData->maxPixCoherence2?dummyData->maxPixCoherence:dummyData->maxPixCoherence2), dummyData->weight);
      //snr_cw->Fill(dummyData->inWindowSNR_V, dummyData->weight);

      //for(int ch=0; ch<16; ch++){
      //   impulsivity[ch] = dummyData->impulsivity[ch];
      //}
      /*
      std::fill(&impulsivity[0], &impulsivity[16], 0.);

      if(isVpolCW && isHpolCW){

      nonZeroCount = 0;
      double sum = 0.;
      for(int ch=0; ch<16; ch++){
         if(fabs( dummyData->impulsivity[ch] - 0 ) > 1e-9 ){
            nonZeroCount++;
            impulsivity[ch] = dummyData->impulsivity[ch];
            sum += impulsivity[ch];
         }
      }

      int index[16];
      TMath::Sort(16, impulsivity, index);
      double avgImpulsivity = sum / (double)nonZeroCount;

      impulsivityHist_max->Fill(impulsivity[index[0]], dummyData->weight);
      impulsivityHist_3rd->Fill(impulsivity[index[2]], dummyData->weight);
      impulsivityHist_avg->Fill(avgImpulsivity, dummyData->weight);

      thermalCWEventCount_both += dummyData->weight;

      }

      else if( isVpolCW && !isHpolCW){

      nonZeroCount = 0;
      double sum = 0.;
      for(int ch=0; ch<8; ch++){
         if(fabs( dummyData->impulsivity[ch] - 0 ) > 1e-9 ){
            nonZeroCount++;
            impulsivity[ch] = dummyData->impulsivity[ch];
            sum += impulsivity[ch];
         }
      }

      int index[16];
      TMath::Sort(16, impulsivity, index);
      double avgImpulsivity = sum / (double)nonZeroCount;

      impulsivityHist_max->Fill(impulsivity[index[0]], dummyData->weight);
      impulsivityHist_3rd->Fill(impulsivity[index[2]], dummyData->weight);
      impulsivityHist_avg->Fill(avgImpulsivity, dummyData->weight);

      thermalCWEventCount_V += dummyData->weight;

      }

      else {

      nonZeroCount = 0;
      double sum = 0.;
      for(int ch=8; ch<16; ch++){
         if(fabs( dummyData->impulsivity[ch] - 0 ) > 1e-9 ){
            nonZeroCount++;
            impulsivity[ch] = dummyData->impulsivity[ch];
            sum += impulsivity[ch];
         }
      }

      int index[16];
      TMath::Sort(16, impulsivity, index);
      double avgImpulsivity = sum / (double)nonZeroCount;

      impulsivityHist_max->Fill(impulsivity[index[0]], dummyData->weight);
      impulsivityHist_3rd->Fill(impulsivity[index[2]], dummyData->weight);
      impulsivityHist_avg->Fill(avgImpulsivity, dummyData->weight);

      thermalCWEventCount_H += dummyData->weight;

      }


      thermalCWEventCount += dummyData->weight;
      */
//   }//if pass all other cuts except CW and thermal


   }//end of entry

delete dataTree;
delete recoSettingsTree;
fp1.Close();

}//end of file

outputFile.close();

printf("totalRecoEventCount: %d\trfEventCount: %f\tcalEventCount: %f\tsoftEventCount: %f\n", totalRecoEventCount, rfEventCount, calEventCount, softEventCount);
//printf("nPassCorruption: %f\tratio: %f\tEvents passed this level: %f\tratio: %f\n", nPassCorruption, (float)nPassCorruption/(float)E19EventCount, nCut1, (float)nCut1/(float)E19EventCount);
//printf("nPassThermalCut: %f\tratio: %f\tEvents passed this level: %f\tratio: %f\n", nPassThermalCut, (float)nPassThermalCut/(float)E19EventCount, nCut2, (float)nCut2/(float)E19EventCount);
//printf("nPassSuE19aceCut: %f\tratio: %f\tEvents passed this level: %f\tratio: %f\n", nPassSuE19aceCut, (float)nPassSuE19aceCut/(float)E19EventCount, nCut3, (float)nCut3/(float)E19EventCount);
//printf("nRecoveredByImp: %f\n", nRecoveredByImp);
//printf("nPassImpulsivityCut: %f\tratio: %f\tEvents passed this level: %f\tratio: %f\n", nPassImpulsivityCut, (float)nPassImpulsivityCut/(float)rfEventCount, nCut0, (float)nCut0/(float)rfEventCount);
//printf("nPassHighPassFilter: %f\tratio: %f\tEvents passed this level: %f\tratio: %f\n", nPassHighPassFilter, (float)nPassHighPassFilter/(float)rfEventCount, nCut1, (float)nCut1/(float)rfEventCount);
//printf("nPassNumSatChanCut %f\tratio: %f\tEvents passed this level: %f\tratio: %f\n", nPassNumSatChanCut, (float)nPassNumSatChanCut/(float)rfEventCount, nCut0p5, (float)nCut0p5/(float)rfEventCount);
printf("nPassCWCut %f\tratio: %f\tEvents passed this level: %f\tratio: %f\n", nPassCWCut, (float)nPassCWCut/(float)rfEventCount, nCut1p5, (float)nCut1p5/(float)rfEventCount);
printf("nPassDeepPulserCut: %f\tratio: %f\tEvents passed this level: %f\tratio: %f\n", nPassDeepPulserCut, (float)nPassDeepPulserCut/(float)rfEventCount, nCut2, (float)nCut2/(float)rfEventCount);
printf("nPassThermalCut: %f\tratio: %f\tEvents passed this level: %f\tratio: %f\n", nPassThermalCut, (float)nPassThermalCut/(float)rfEventCount, nCut3, (float)nCut3/(float)rfEventCount);
//printf("nPassThermalImpulsivityCut: %f\tratio: %f\tEvents passed this level: %f\tratio: %f\n", nPassThermalImpulsivityCut, (float)nPassThermalImpulsivityCut/(float)rfEventCount, nCut3p5, (float)nCut3p5/(float)rfEventCount);
printf("nPassSNRCut: %f\tratio: %f\tEvents passed this level: %f\tratio: %f\n", nPassSNRCut, (float)nPassSNRCut/(float)rfEventCount, nCut3p5, (float)nCut3p5/(float)rfEventCount);
printf("nPassCalpulserCut: %f\tratio: %f\tEvents passed this level: %f\tratio: %f\n", nPassCalpulserCut, (float)nPassCalpulserCut/(float)rfEventCount, nCut4, (float)nCut4/(float)rfEventCount);
printf("nPassCalpulserTimeCut: %f\tratio: %f\tEvents passed this level: %f\tratio: %f\n", nPassCalpulserTimeCut, (float)nPassCalpulserTimeCut/(float)rfEventCount, nCut4p5, (float)nCut4p5/(float)rfEventCount);
//printf("nPassNoisyRunCut: %f\tratio: %f\tEvents passed this level: %f\tratio: %f\n", nPassNoisyRunCut, (float)nPassNoisyRunCut/(float)rfEventCount, nCut5, (float)nCut5/(float)rfEventCount);
printf("nPassSurfaceCut: %f\tratio: %f\tEvents passed this level: %f\tratio: %f\n", nPassSurfaceCut, (float)nPassSurfaceCut/(float)rfEventCount, nCut6, (float)nCut6/(float)rfEventCount);
printf("nPassSurfaceCut_2: %f\tratio: %f\tEvents passed this level: %f\tratio: %f\n", nPassSurfaceCut_2, (float)nPassSurfaceCut_2/(float)rfEventCount, nCut6p5, (float)nCut6p5/(float)rfEventCount);
printf("nPassNoisyRunCut: %f\tratio: %f\tEvents passed this level: %f\tratio: %f\n", nPassNoisyRunCut, (float)nPassNoisyRunCut/(float)rfEventCount, nCut7, (float)nCut7/(float)rfEventCount);
cout<<"rfEventCount: "<<rfEventCount<<" isCWCount: "<<isCWCount<<" ratio: "<<isCWCount/rfEventCount<<" 1e-2 background ratio: "<<1e-2/(isCWCount*10)<<endl;
cout<<"inBand_all: "<<inBand_all<<" inBand_pass: "<<inBand_pass<<" inBand_cut: "<<inBand_all-inBand_pass<<endl;
cout<<"outOfBand_all: "<<outOfBand_all<<" outOfBand_pass: "<<outOfBand_pass<<" outOfBand_cut: "<<outOfBand_all-outOfBand_pass<<endl;
cout<<"pre-thermal: "<<inBand_all+outOfBand_all<<" post-thermal: "<<inBand_pass+outOfBand_pass<<" ratio: "<<(inBand_pass+outOfBand_pass)/(inBand_all+outOfBand_all)<<endl;;
//outputFile<<ENERGY<<","<<totalTrigEventCount<<","<<totalTrigEventCount-totalOffsetBlockEventCount<<","<<totalTrigEventCount-totalOffsetBlockEventCount-totalImpulsivityFilteredEventCount<<",";
//outputFile/*<<nCut0p5<<","*/<<nCut1p5<<","<<nCut2<<","<<nCut3/*<<","<<nCut3p5*/<<","<<nCut4<<","<<nCut6<<","<<nCut7<<endl;
//outputFile.close();

//cout<<"passCWCut: "<<passCWCut<<endl;
//cout<<"thermalCWEventCount: "<<thermalCWEventCount<<endl;
//cout<<"thermalCWEventCount_both: "<<thermalCWEventCount_both<<" thermalCWEventCount_V: "<<thermalCWEventCount_V<<" thermalCWEVentCount_H: "<<thermalCWEventCount_H<<endl;
//delete dummyData;
delete runInfoTree;
//delete dataTree;
//delete recoSettingsTree;

//TFile fout(argv[3], "UPDATE");

char filename[200];
/*
TCanvas c1("c1","c1",1200, 800);
//c1.Divide(2,1);
//c1.cd(1);
maxCountFreq_V_hist->SetLineColor(kBlack);
maxCountFreq_V_hist->Draw();
maxCountFreq_H_hist->SetLineColor(kRed);
maxCountFreq_H_hist->Draw("same");
//sprintf(filename, "%s_type%d_maxCountFreq_cwEvents.C", STATION.c_str(), type);
c1.SetGrid();
c1.SetLogy();
maxCountFreq_V_hist->SetTitle("Vpol;[MHz];entry");
maxCountFreq_V_hist->GetYaxis()->SetRangeUser(1e-3, 1e2);
maxCountFreq_H_hist->SetTitle("Hpol");
c1.BuildLegend();
sprintf(filename,"%s Config %d MC CW Event Max Count Frequency", STATION.c_str(), type);
maxCountFreq_V_hist->SetTitle(filename);
sprintf(filename, "%s_type%d_allEnergies_maxCountFreq_cwEvents.C", STATION.c_str(), type);
c1.SaveAs(filename);
//sprintf(filename, "%s_type%d_maxCountFreq_cwEvents.pdf", STATION.c_str(), type);
//c1.SaveAs(filename);
*/
/*
TCanvas c2("c2","c2",800,800);
//runHist->Draw();
suE19aceRunHist->Draw();
//sprintf(filename,"runHist_suE19aceCut_%s_type%d.C", STATION.c_str(), type);
sprintf(filename,"suE19aceRunHist_%s_type%d.C", STATION.c_str(), type);
c2.SaveAs(filename);
//sprintf(filename,"runHistFile_suE19aceCut_%s.root", STATION.c_str());
sprintf(filename,"suE19aceRunHistFile_suE19aceCut_%s.root", STATION.c_str());
TFile fout(filename, "UPDATE");
sprintf(filename,"type_%d",type);
//runHist->SetName(filename);
//runHist->Write();
suE19aceRunHist->SetName(filename);
suE19aceRunHist->Write();
fout.Close();
*/

//sprintf(filename, "%s_type%d_E18_coherence_snr_cw.C", STATION.c_str(), type);
TCanvas c34("c34","c34",1200,800);
c34.Divide(2,1);
c34.cd(1);
coherence_snr_cw->Draw("colz");
c34.cd(2);
snr_cw->Draw();
//c34.SaveAs(filename);

/*
sprintf(filename, "%s_type%d_E18_impulsivity_cw.C", STATION.c_str(), type);
TCanvas c35("c35","c35",800,800);
impulsivityHist_max->SetLineColor(kBlack);
impulsivityHist_3rd->SetLineColor(kRed);
impulsivityHist_avg->SetLineColor(kBlue);
impulsivityHist_max->Draw();
impulsivityHist_3rd->Draw("same");
impulsivityHist_avg->Draw("same");
c35.SaveAs(filename);
*/
/*
sprintf(filename, "max_E19");
impulsivityHist_max->SetName(filename);
sprintf(filename, "third_E19");
impulsivityHist_3rd->SetName(filename);
sprintf(filename, "avg_E19");
impulsivityHist_avg->SetName(filename);
fout.cd();
impulsivityHist_max->Write();
impulsivityHist_3rd->Write();
impulsivityHist_avg->Write();
fout.Close();
*/

/*
TCanvas c4("c4","c4",800,800);
surfaceRunHist->Draw();
sprintf(filename,"surfaceRunHist_vnchnl3NoMasking_beforeImpCut_noMaskSat_snrMode1_ch6Fit2Corr_2SurfaceCut_%s_type%d.C", STATION.c_str(), type);
c4.SaveAs(filename);

sprintf(filename, "%s_vnchnl3NoMasking_beforeImpCut_noMaskSat_snrMode1_ch6Fit2Corr_2SurfaceCut_surfaceEvents_noisyRuns.txt", STATION.c_str());
ofstream noisyRunFile(filename, ofstream::out|ofstream::app);
for(int bin=1; bin<=surfaceRunHist->GetNbinsX(); bin++){
   if(surfaceRunHist->GetBinContent(bin)>1) noisyRunFile<<surfaceRunHist->GetBinCenter(bin)<<endl;
}
noisyRunFile.close();


TCanvas c5("c5","c5",1200,800);
c5.Divide(2,1);
c5.cd(1);
//impulsivityHist_avg->Draw();
//c5.cd(2);
//impulsivityHist_3rd->Draw();
//sprintf(filename, "%s_type%d_passAllCutsImpulsivity_vnchnl3NoMasking_noMaskSat.C", STATION.c_str(), type);
//sprintf(filename, "%s_type%d_passThermalSurfaceDeepPulserCutOnlyImpulsivity_vnchnl3NoMasking_noMaskSat.C", STATION.c_str(), type);
constantNZenHist->Draw();
constantNZenHist->SetTitle("Quasi-planewave Reco;Reco Zenith [#circ];Entry");
c5.cd(2);
iterMajorityZenHist->Draw();
iterMajorityZenHist->SetTitle("Iter. Majority Reco:Reco Zenith [#circ];Entry");
sprintf(filename, "%s_type%d_vnchnl3NoMasking_noMaskSat_snrMode1_ch6Fit2Corr_constantNZen_iterMajorityZenHist.C", STATION.c_str(), type);
c5.SaveAs(filename);
*/
/*
snrHist[0]->SetLineColor(kBlue);
snrHist[1]->SetLineColor(kRed);
snrHist[2]->SetLineColor(kTeal);
snrHist[3]->SetLineColor(kOrange);
snrHist[4]->SetLineColor(kMagenta);
snrHist[0]->Draw();
for(int i=1; i<5; i++) snrHist[i]->Draw("same");
sprintf(filename, "%s_type%d_signalEffiencyVsSNR.C", STATION.c_str(), type);
c5.SaveAs(filename);
sprintf(filename, "%s_type%d_signalEffiencyVsSNR.root", STATION.c_str(), type);
TFile ff(filename, "RECREATE");
for(int i=0; i<5; i++) snrHist[i]->Write();
ff.Close();
*/

TCanvas c6("c6","c6",800,800);
//impulsivityHist_nMinusCW->Draw();
//impulsivityHist_nMinusCW->SetTitle(";Impulsivity;Entry");
coherence_nMinusThermal->SetTitle(";Coherence;Entry");
coherence_nMinusThermal->Draw();
sprintf(filename, "%s_type%d_nMinusThermal_coherence.C", STATION.c_str(), type);
c6.SaveAs(filename);

TCanvas c7("c7","c7",800,800);
//impulsivityHist_nMinusImp->Draw();
//impulsivityHist_nMinusImp->SetTitle(";Impulsivity;Entry");
snr_nMinusSNR->SetTitle(";SNR;Entry");
snr_nMinusSNR->Draw();
sprintf(filename, "%s_type%d_nMinusSNR_snr.C", STATION.c_str(), type);
c7.SaveAs(filename);

TCanvas c8("c8","c8",800,800);
zen_nMinusSurface->Draw();
zen_nMinusSurface->SetTitle(";Receipt Angle [#circ];Entry");
sprintf(filename, "%s_type%d_nMinusSurface_zen.C", STATION.c_str(), type);
c8.SaveAs(filename);
/*
TCanvas c9("c9","c9",800,800);
c_vs_snr_hist_nMinusThermal->Draw("colz");
c_vs_snr_hist_nMinusThermal->SetTitle(";SNR;Coherence");
sprintf(filename, "%s_type%d_nMinusThermal_c_vs_snr.C", STATION.c_str(), type);
c9.SaveAs(filename);
*/
TCanvas c10("c10","c10",800,800);
zen_azi_nMinusCal->Draw("colz");
zen_azi_nMinusCal->SetTitle(";Azimuth [#circ];Zenith [#circ]");
sprintf(filename, "%s_type%d_nMinusCal_zen_azi.C", STATION.c_str(), type);
c10.SaveAs(filename);
/*
TCanvas c11("c11","c11",800,800);
c_vs_imp->Draw("colz");
c_vs_imp->SetTitle(";Impulsivity;Coherence");
sprintf(filename, "%s_type%d_c_vs_imp.C", STATION.c_str(), type);
c11.SaveAs(filename);
*/
/*
TCanvas c12("c12","c12",1200,800);
c12.Divide(2,1);
c12.cd(1);
impulsivityHist_avg->Draw();
TLine impline(postThermalAvgImpulsivityCut, 0, postThermalAvgImpulsivityCut, impulsivityHist_avg->GetMaximum());
impline.Draw("same");

c12.cd(2);
c_vs_snr_hist->Draw("colz");
TLine snrlineinband(cutValues->snrCut_inBand[type-1].val, 0, cutValues->snrCut_inBand[type-1].val, cutValues->coherenceCut_inBand[type-1].val);
TLine coherencelineinband(0, cutValues->coherenceCut_inBand[type-1].val, cutValues->snrCut_inBand[type-1].val, cutValues->coherenceCut_inBand[type-1].val);
TLine snrlineoutband(cutValues->snrCut_outOfBand[type-1].val, 0, cutValues->snrCut_outOfBand[type-1].val, cutValues->coherenceCut_outOfBand[type-1].val);
TLine coherencelineoutband(0, cutValues->coherenceCut_outOfBand[type-1].val, cutValues->snrCut_outOfBand[type-1].val, cutValues->coherenceCut_outOfBand[type-1].val);


snrlineinband.Draw("same");
coherencelineinband.Draw("same");
snrlineoutband.SetLineColor(kRed);
coherencelineoutband.SetLineColor(kRed);
snrlineoutband.Draw("same");
coherencelineoutband.Draw("same");

c12.SaveAs("recoAnalysis_12.C");
*/
/*
TCanvas c13("c13","c13",1200,800);
c13.Divide(3,1);
c13.cd(1);
bipolarRatio_dT->Draw("colz");

c13.cd(2);
impulsivity_dT->Draw("colz");

c13.cd(3);
impulsivity_bipolarRatio->Draw("colz");

c13.SaveAs("recoAnalysis_13.C");
*/

/*
TCanvas c14("c14","c14",1200,800);
c14.Divide(2,1);
c14.cd(1);
dFHist_V->Draw();

c14.cd(2);
dFHist_H->Draw();

c14.SaveAs("recoAnalysis_14.C");
*/
/*
TCanvas c15("c15","c15",800,800);
c15.Divide(2,2);
c15.cd(1);
thetaXingHist->Draw();
c15.cd(2);
phiXingHist->Draw();
c15.cd(3);
avgThetaXingHist->Draw();
c15.cd(4);
avgPhiXingHist->Draw();
c15.SaveAs("recoAnalysis_15.C");


TH1F *thetaXingHist_cumu = getCumulative(thetaXingHist);
TH1F *phiXingHist_cumu   = getCumulative(phiXingHist);
TH1F *avgThetaXingHist_cumu = getCumulative(avgThetaXingHist);
TH1F *avgPhiXingHist_cumu   = getCumulative(avgPhiXingHist);

TCanvas c16("c16","c16",800,800);
c16.Divide(2,2);
c16.cd(1);
thetaXingHist_cumu->Draw();
c16.cd(2);
phiXingHist_cumu->Draw();
c16.cd(3);
avgThetaXingHist_cumu->Draw();
c16.cd(4);
avgPhiXingHist_cumu->Draw();
c16.SaveAs("recoAnalysis_16.C");
*/
/*
sprintf(filename,"%s_type%d_nMinusThermalImp_inRangeThetaPhiFraction_angThres0.5.C",STATION.c_str(),type);
TCanvas c17("c17","c17",1200,800);
c17.Divide(3,1);
c17.cd(1);
inRangeThetaFracHist->Draw();
c17.cd(2);
inRangePhiFracHist->Draw();
c17.cd(3);
inRangeThetaPhiFracHist->Draw("colz");
c17.SaveAs(filename);
*/

//sprintf(filename,"%s_type%d_coherenceThermalCut_snr.C",STATION.c_str(),type);
//TCanvas c18("c18","c18",/*800*/1200,800);
//c18.Divide(2,1);
//c18.cd(1);
//_snrHist->Draw();
//c18.cd(2);
//_snrCumuHist = getCumulative(_snrHist);
//_snrCumuHist->SetTitle("_snrCumuHist");
//_snrCumuHist->SetName("_snrCumuHist");
//_snrCumuHist->Draw();
//c18.SaveAs(filename);


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
/*
bool isNearNoisyRun(string station, const vector<int>& noisyRuns, int runNum, int plusMinusRunNum){


   bool isNoisy = false;

   for(int run=0; run<noisyRuns.size(); run++){
      if(fabs(noisyRuns[run] - runNum) <= plusMinusRunNum){
         isNoisy=true;
         break;
      }
   }

   return isNoisy;

}
*/
/*
float getZenMaj(const vector<float>& iterZenVec, float zenRange){


   float zenMaj = 100.f;

   vector<int> cntVec;
   int maxCnt=0;
   int maxCntIdx = -1;

   for(int i=0; i<(int)iterZenVec.size(); i++){

      int cnt=0;

      for(int j=0; j<(int)iterZenVec.size(); j++){

         if(fabs(iterZenVec[i]-iterZenVec[j]) < zenRange) cnt++;

      }//end of j
      cntVec.push_back(cnt);
      if(cnt>maxCnt){

         maxCnt = cnt;
         maxCntIdx = i;

      } else if (cnt == maxCnt){

         if(iterZenVec[i] > iterZenVec[maxCntIdx]){ //pick the larger value if two angles have the same count
          maxCntIdx = i;
         }

      }

   }//end of i

   zenMaj = (maxCnt>1?iterZenVec[maxCntIdx]:100.f);
   return zenMaj;
}
*/
