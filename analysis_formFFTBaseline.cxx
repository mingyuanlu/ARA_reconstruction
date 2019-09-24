#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <string>
#include <sys/stat.h>
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
#include "TCanvas.h"
#include "TTree.h"
#include "TObject.h"

#include "calibrationTools.h"
#include "calibrationToolsVs3.h"
#include "evProcessTools.h"
#include "Healpix_Onion.h"
#include "recoTools.h"
#include "recoSettings.h"
#include "recoData.h"
#include "trackEngine.h"
#include "analysisTools.h"

#include "Detector.h"
#include "Trigger.h"
#include "Settings.h"
#include "Report.h"
#include "Event.h"
#include "Position.h"
#include "signal.hh"
#include "IceModel.h"

ClassImp(recoSettings);
ClassImp(recoData);

using namespace std;

RawIcrrStationEvent *rawIcrrEvPtr;
RawAtriStationEvent *rawAtriEvPtr;
RawAraStationEvent *rawEvPtr;

double getMean(TGraph *);
double computeWeight(Settings *settings, Detector *detector, Event *event, IceModel *antarctica, double zCenter, double L_int);//Computes the forc-interaction-corrected weight, which is simply Interaction::weight * L_0 / L_int, where L_0 is the generation volume thickness as seen by the neutrino. In other words, distance from generation volume entry to exit
/*
#ifndef XCORRSUMGRAPH
static TGraph *sillygr = new TGraph();
//TGraph *envelopeSum = new TGraph();
#define XCORRSUMGRAPH
#endif
*/

int main( int argc, char **argv){

gROOT->ProcessLine("#include <vector>");
/*
timer *tmr;
resetTimer(tmr);
recordTime(tmr,0);
*/
time_t t_program_start = time(NULL);
clock_t c_program_start = clock();

if(argc < 4){ cerr<<"Insufficient arguments. Usage: 1. Output directory 2.recoSetupFile 3. Run Number 4. Data ROOT file 5. AraSim: more Data ROOT file. Real: Pedestal File\n"; return -1; }

int err;
//gROOT->ProcessLine("#include <vector>");
/*
 * Specify the channels to be used in the analysis
 * 1: use, 0: don't use
 */
/*
const int chanMask[16] = {  1 //chan 0  D1TV
                           ,1 //chan 1  D2TV
                           ,1 //chan 2  D3TV
                           ,1 //chan 3  D4TV
                           ,1 //chan 4  D1BV
                           ,1 //chan 5  D2BV
                           ,1 //chan 6  D3BV
                           ,1 //chan 7  D4BV
                           ,1 //chan 8  D1TH
                           ,1 //chan 9  D2TH
                           ,1 //chan 10 D3TH
                           ,1 //chan 11 D4TH
                           ,1 //chan 12 D1BH
                           ,1 //chan 13 D2BH
                           ,1 //chan 14 D3BH
                           ,1 //chan 15 D4BH
                         };
*/


recoSettings *settings = new recoSettings();

cout<<"Default settings: "<<endl;
cout<<"programFile: "<<settings->programFile<<endl;
cout<<"nSideExp: "<<settings->nSideExp<<endl;
cout<<"nLayer: "<<settings->nLayer<<endl;
cout<<"beamformMethod: "<<settings->beamformMethod<<endl;
cout<<"recoVertexingMode: "<<settings->recoVertexingMode<<endl;
cout<<"skymapSearchMode: "<<settings->skymapSearchMode<<endl;
cout<<"recoPolType: "<<settings->recoPolType<<endl;
cout<<"nchnlFilter: "<<settings->nchnlFilter<<endl;
cout<<"recoEventIndex: "<<settings->recoEventIndex<<endl;

TFile *outputFile;

string outputDir;
struct stat sb;
if( stat(argv[1], &sb) == 0 && S_ISDIR(sb.st_mode) ){
   outputDir = string( argv[1] ) + "/";
} else { cerr<<"outputDir "<<argv[1]<<" directory does not exist! Aborting...\n"; return -1; }

string recoSetupFile_fullPath = string( argv[2] );
string recoSetupFile = string( basename(argv[2]) );
string runNum = string( argv[3] );
string fitsFile_tmp;
string fitsFileStr;
string evStr;
cout<<"recoStupFile: "<<recoSetupFile<<endl;
if( !settings->readRecoSetupFile( recoSetupFile_fullPath )){

   cerr<<"Error reading the recoSetupFile or invalid parameters !! Aborting now...\n";
   return -1;
   //cerr<<"Will use default reco setup file\n";
   //settings->readRecoSetupFile("recoSetupFile_default.txt");
   //outputFile = new TFile(("recoOut.recoSetupFile_default.txt.run"+runNum+".root").c_str(),"RECREATE","recoOut");
   //fitsFile_tmp = "recoSkymap.recoSetupFile_default.txt.run" + runNum/* + ".fits"*/;

} else {
   cout<<"Obtained new reoSetupFile\n";
   //outputFile = new TFile((outputDir+"fftPedestal."+recoSetupFile+".run"+runNum+".root").c_str(),"RECREATE","fftPedestal");
   fitsFile_tmp = outputDir + "recoSkymap." + recoSetupFile + ".run" + runNum/* + ".fits"*/;
}
char fitsFile[200];
//sprintf(fitsFile, fitsFile_tmp.c_str());
cout<<"New settings: "<<endl;
cout<<"programFile: "<<settings->programFile<<endl;
cout<<"nSideExp: "<<settings->nSideExp<<endl;
cout<<"nLayer: "<<settings->nLayer<<endl;
cout<<"beamformMethod: "<<settings->beamformMethod<<endl;
cout<<"recoVertexingMode: "<<settings->recoVertexingMode<<endl;
cout<<"skymapSearchMode: "<<settings->skymapSearchMode<<endl;
cout<<"recoPolType: "<<settings->recoPolType<<endl;
cout<<"nchnlFilter: "<<settings->nchnlFilter<<endl;
cout<<"recoEventIndex: "<<settings->recoEventIndex<<endl;

TTree *recoSettingsTree = new TTree("recoSettingsTree", "recoSettingsTree");
//TTree *onionTree = new TTree("onionTree", "onionTree");
TTree *dataTree  = new TTree("dataTree",  "dataTree");
TTree *runInfoTree = new TTree("runInfoTree", "runInfoTree");

recoSettingsTree->Branch("settings", &settings);
recoSettingsTree->Fill();

TH1F *dZenDist = new TH1F("recoZenDiff", "recoZenDiff", 360, -180, 180);
TH1F *dAziDist = new TH1F("recoAziDiff", "recoAziDiff", 720, -360, 360);
TH2F *recoTrueZenDist = new TH2F("recoTrueZenDist", "recoTrueZenDist", 180, 0, 180, 180, 0, 180);
TH2F *recoTrueAziDist = new TH2F("recoTrueAziDist", "recoTrueAziDist", 360, 0, 360, 360, 0, 360);

/*
 * Specify the channels to be used in the analysis
 * 1: use, 0: don't use
 */

int chanMask[16];
for(int ch=0; ch<16; ch++){
  chanMask[ch] = settings->chanMask[ch] - '0';
  cout<<"chanMask "<<ch<<": "<<chanMask[ch]<<endl;
}


float radius, r_xy;
float dx, dy, dz;
float r_true, zen_true, azi_true;
float recoRecAngles[16], recoLauAngles[16], trueRecAngles[16], trueLauAngles[16];
bool recoSuccess;
int *snrRank;
bool dropARA02D4BH, dropARA03D4;
dropARA02D4BH = dropARA03D4 = false;

/*
 * Variables used in dataType == 0 case
 */

TChain chain("AraTree"), chain2("AraTree2");
Settings *AraSim_settings = 0;
Detector *detector = 0;
Event    *event    = 0;
Report   *report   = 0;
Trigger  *trigger  = 0;
IceModel *icemodel = 0;
double weight;

/*
 * Variables used in dataType == 1 case
 */

TFile *fp;
TTree *eventTree;
int isIcrrEvent, isAtriEvent;
char dir_char[200];
AraEventCalibrator *calib;
vector<vector<double> > pulserLocation;
AraGeomTool *geom;
int cutWaveAlert, nonIncreasingSampleTimeAlert, mistaggedSoftEventAlert;
double previous_times;
double addDelay;
double times, volts;
double time_1, time_2, time_last;
int utime, utime_runStart, utime_runEnd;
utime_runStart = utime_runEnd = 0;

/*
 * Variables used in detecting offset blocks
 */

 int offsetBlockAlert;
 double maxTime[16];
 double meanMax[16];
 double threshold_V = settings->offsetBlock_threshold_V;
 double threshold_H = settings->offsetBlock_threshold_H;
 double timeRangeCut = settings->offsetBlock_timeRangeCut;
 int nChanBelowThres;
 int nChanBelowThres_Thres;
 double timeRange;
 vector<double> maxTimeVec;

/*
 * Variables used in nchnlFilter > 0 case
 */

double threshold;
int nchnlArray[3];
int nchnl_tmp;

/*
 * Variables used in cwFilter > 0 case
 */

int minCWCoincidence = settings->minCWCoincidence;
int maxCountBin, maxCount_V, maxCount_H;
double maxCountFreq_V, maxCountFreq_H;
int freqCountLen_V, freqCountLen_H;
int *freqCount_V, *freqCount_H;
double freqBinWidth_V, freqBinWidth_H;

/* End of conditional variables pre-declaration */

int runEventCount, trigEventCount, recoEventCount;
runEventCount = trigEventCount = recoEventCount = 0;
int runRFEventCount, runCalEventCount, runSoftEventCount;
runRFEventCount = runCalEventCount = runSoftEventCount = 0;
int cutWaveEventCount, nonIncreasingSampleTimeEventCount, cutWaveAndNonIncreasingEventCount;
cutWaveEventCount = nonIncreasingSampleTimeEventCount = cutWaveAndNonIncreasingEventCount = 0;
int mistaggedSoftEventCount, offsetBlockEventCount;
mistaggedSoftEventCount = offsetBlockEventCount = 0;
int nchnlFilteredEventCount, cwFilteredEventCount, impulsivityFilteredEventCount;
nchnlFilteredEventCount = cwFilteredEventCount = impulsivityFilteredEventCount = 0;
int corruptFirst3EventCount = 0;
int corruptD1EventCount = 0;
double weightedTrigEventCount = 0.;
double weightedRecoEventCount = 0.;
double weightedOffsetBlockEventCount = 0.;
double weightedNchnlFilteredEventCount = 0.;
double weightedCWFilteredEventCount = 0.;
double weightedImpulsivityFilteredEventCount = 0.;
runInfoTree->Branch("runEventCount",  &runEventCount);
runInfoTree->Branch("runRFEventCount", &runRFEventCount);
runInfoTree->Branch("runCalEventCount", &runCalEventCount);
runInfoTree->Branch("runSoftEventCount", &runSoftEventCount);
runInfoTree->Branch("trigEventCount", &trigEventCount);
runInfoTree->Branch("recoEventCount", &recoEventCount);
runInfoTree->Branch("utime_runStart", &utime_runStart);
runInfoTree->Branch("utime_runEnd",   &utime_runEnd);
runInfoTree->Branch("cutWaveEventCount", &cutWaveEventCount);
runInfoTree->Branch("nonIncreasingSampleTimeEventCount", &nonIncreasingSampleTimeEventCount);
runInfoTree->Branch("cutWaveAndNonIncreasingEventCount", &cutWaveAndNonIncreasingEventCount);
runInfoTree->Branch("mistaggedSoftEventCount", &mistaggedSoftEventCount);
runInfoTree->Branch("offsetBlockEventCount", &offsetBlockEventCount);
runInfoTree->Branch("cwFilteredEventCount", &cwFilteredEventCount);
runInfoTree->Branch("nchnlFilteredEventCount", &nchnlFilteredEventCount);
runInfoTree->Branch("impulsivityFilteredEventCount", &impulsivityFilteredEventCount);
runInfoTree->Branch("corruptFirst3EventCount", &corruptFirst3EventCount);
runInfoTree->Branch("corruptD1EventCount", &corruptD1EventCount);
runInfoTree->Branch("weightedTrigEventCount", &weightedTrigEventCount);
runInfoTree->Branch("weightedRecoEventCount", &weightedRecoEventCount);
runInfoTree->Branch("weightedOffsetBlockEventCount", &weightedOffsetBlockEventCount);
runInfoTree->Branch("weightedNchnlFilteredEventCount", &weightedNchnlFilteredEventCount);
runInfoTree->Branch("weightedCWFilteredEventCount", &weightedCWFilteredEventCount);
runInfoTree->Branch("weightedImpulsivityFilteredEventCount", &weightedImpulsivityFilteredEventCount);

if(settings->dataType == 1)//real events
{

   fp = TFile::Open( argv[4] );
   if ( !fp ) { cerr<<"can't open file"<<endl; return -1; }

   eventTree = (TTree*) fp->Get("eventTree");
   if ( !eventTree ){ cerr<<"can't find eventTree"<<endl; return -1; }
   cout<<"evenTree Nentries: "<<eventTree->GetEntries()<<endl;

   //Now check the electronics type of the station
   isIcrrEvent=0;
   isAtriEvent=0;

   //Check an event in the run Tree and see if it is station1 or TestBed (stationId<2)
   eventTree->SetBranchAddress("event",&rawEvPtr);
   eventTree->GetEntry(0);

   if((rawEvPtr->stationId)<2){ isIcrrEvent=1; isAtriEvent=0; }
   else                       { isIcrrEvent=0; isAtriEvent=1; }
   eventTree->ResetBranchAddresses();

   //Now set the appropriate branch addresses
   //The Icrr case
   if(isIcrrEvent){
   eventTree->SetBranchAddress("event", &rawIcrrEvPtr);
   cerr<<"Set Branch address to Icrr\n";
   }
   //The Atri case
   else{
   eventTree->SetBranchAddress("event", &rawAtriEvPtr);
   cerr<<"Set Branch address to Atri\n";
   }
   runEventCount=eventTree->GetEntries();
   cerr<<"isAtri "<<isAtriEvent<<" isIcrr "<<isIcrrEvent<<" number of entries is "<<runEventCount<<endl;

   /*
    * START LOADING GOOD PED
    */
   sprintf(dir_char,argv[5]);
   printf("ped file is %s\n",dir_char);
   calib = AraEventCalibrator::Instance();
   calib->setAtriPedFile(dir_char,rawEvPtr->stationId);
   //************END OF LOADING GOOD PED************
}
else if (settings->dataType == 0)//AraSim events
{

   for(int i=4; i<argc; i++){
      chain.Add( argv[i] );
      chain2.Add( argv[i] );
   }

   chain.SetBranchAddress("settings",&AraSim_settings);
   chain.SetBranchAddress("detector",&detector);
   chain.SetBranchAddress("trigger" ,&trigger);
   chain.SetBranchAddress("icemodel", &icemodel);
   chain2.SetBranchAddress("report" ,&report);
   chain2.SetBranchAddress("event"  ,&event);

   chain.GetEntry(0);
   cout<<"EXPONENT: "<<AraSim_settings->EXPONENT<<endl;
   cout<<"NNU: "<<AraSim_settings->NNU<<endl;
   printf("Station center X: %f Y: %f Z: %f\n",detector->stations[0].GetX(),detector->stations[0].GetY(),detector->stations[0].GetZ());

   runEventCount = chain2.GetEntries();

}
else
{ cerr<<"Undefined dataType!\n"; return -1; }

/*
 * START CALIBRATING STATION GEOMETRY AND DELAYS
 */
double delays[4][4] ={{0}};
double pulserCorr[5]={0};

/* Set station center to be  180m below ice surface */
float stationCenterDepth = 180.f;

vector<vector<double> > antLocation;

if(settings->dataType == 1){

   //vector<vector<double> > pulserLocation;
   err = calibrateGeometryAndDelays(rawEvPtr, delays, pulserCorr, stationCenterDepth, antLocation, pulserLocation);
   if( err<0 ){ cerr<<"Error calibrating geometry and delays\n"; return -1; }

} else {

   err = getAraSimStationGeometry(antLocation, detector, AraSim_settings);
   if( err<0 ){ cerr<<"Error loading AraSim geometry\n"; return -1; }

}

if(settings->applyA2Ch6Correction){
   /* Add correction from Fit 2 of D5BV+D6BV calpulsers to D3BV (ch6) Z */
   antLocation[6][2] += (-0.6);
}

/* Initiate a track engine instance, build baseline vectors */
trackEngine *treg = new trackEngine();
treg->initialize();
treg->buildBaselineTracks(antLocation);


/*
 * Start computing reco delays using RadioSpline
 */

/* Initialize a AraGeomTool instance */

int nAnt;
if(settings->dataType == 1){
   geom = AraGeomTool::Instance();
   nAnt = geom->getStationInfo(rawEvPtr->stationId)->getNumAnts() -   //FIXME this works because, coincidently, # of calpulerse =
          geom->getStationInfo(rawEvPtr->stationId)->getNumCalAnts(); // # of surface antennas = 4
} else {
   nAnt = 16;
}
cout<<"nAnt: "<<nAnt<<endl;

/* Set top N max pixels in whole Healpix_Onion */
int topN = settings->topN;

int nSideExp;
int nLayer, nDir;
float *recoDelays, *recoDelays_V, *recoDelays_H;
float *recoRefracDelays, *recoRefracDelays_V, *recoRefracDelays_H;
Healpix_Onion *onion;

//A2 D5BV
//double cal_r = 42;
//double cal_zen = 112.76;
//double cal_azi = 334.855;
//double calLocation[3] = {cal_r, cal_zen, cal_azi};
//float calRecoDelays[16], calRecoDelays_V[8], calRecoDelays_H[8];

/* Print calpulser location delays */
//compute3DRecoDelaysWithRadioSplineForSinglePoint_sphericalCoordInDeg(nAnt, -1.f*stationCenterDepth, antLocation, calRecoDelays, calRecoDelays_V, calRecoDelays_H, calLocation);

//for(int ch=0; ch<16; ch++) cout<<"calRecoDelays "<<ch<<": "<<calRecoDelays[ch]<<endl;

//recordTime(tmr,1);
time_t t_before_recoDelays = time(NULL);
clock_t c_before_recoDelays = clock();

if( settings->skymapSearchMode == 0){ //No zoom search

/* Set n-side for Healpix and compute reco delays */

   nSideExp = settings->nSideExp;
   nLayer = settings->nLayer;
   nDir = 12 * pow(2, nSideExp) * pow(2, nSideExp);
   onion = new Healpix_Onion(nSideExp, nLayer, settings->layerFirstRadius, settings->layerLastRadius);


   if(nDir*nLayer < topN) {
      cerr<<"topN greater than total number of pixles. Setting topN to equal nDir*nLayer...\n";
      topN = nDir*nLayer; //in case topN < total number of pixels
   }
//
//   recoDelays  = (float*)malloc(nLayer*nDir*nAnt*sizeof(float));
//   recoDelays_V= (float*)malloc(nLayer*nDir*(nAnt/2)*sizeof(float));
//   recoDelays_H= (float*)malloc(nLayer*nDir*(nAnt/2)*sizeof(float));
//
//   if(settings->use2ndRayReco == 1){
//      recoRefracDelays  = (float*)malloc(nLayer*nDir*nAnt*sizeof(float));
//      recoRefracDelays_V= (float*)malloc(nLayer*nDir*(nAnt/2)*sizeof(float));
//      recoRefracDelays_H= (float*)malloc(nLayer*nDir*(nAnt/2)*sizeof(float));
//   }
//
//   if(settings->iceModel == 1){
//   err = computeRecoDelaysWithConstantN(nAnt, -1.f*stationCenterDepth, antLocation,
//                                        //radius, nSideExp,
//                                        onion, recoDelays, recoDelays_V, recoDelays_H);
//   if(settings->use2ndRayReco == 1){ cerr<<"Constant N ice model incompatible with use2ndRayReco == 1\n"; return -1; }
//   } else if(settings->iceModel == 0){
//   if(settings->use2ndRayReco == 0)
//   err = compute3DRecoDelaysWithRadioSpline(nAnt, -1.f*stationCenterDepth, antLocation,
//                                            onion, recoDelays, recoDelays_V, recoDelays_H);
//   else
//   err = compute3DRecoBothDelaysWithRadioSpline(nAnt, -1.f*stationCenterDepth, antLocation,
//                                               onion, recoDelays, recoDelays_V, recoDelays_H,
//                                               recoRefracDelays, recoRefracDelays_V, recoRefracDelays_H);
//
//   } else { cerr<<"Undefined iceModel parameter\n"; return -1; }
//   if( err<0 ){ cerr<<"Error computing reco delays\n"; return -1; }

} else {// zoom search mode

   //Dir = 12 * pow(2, settings->nSideExpEnd) * pow(2, settings->nSideExpEnd);
/*
   nDir = pow(4, settings->nSideExpEnd - settings->nSideExpStart);
   nLayer = settings->nLayer;
   if(nDir*nLayer < topN){
      cerr<<"topN greater than total number of final pixles. Setting topN to equal nDir*nLayer...\n";
      topN = nDir*nLayer; //in case topN < total number of pixels
   }
*/
}

float *constantNDelays, *constantNDelays_V, *constantNDelays_H;
Healpix_Onion *onion_temp;

//if( settings->constantNFilter > 0){
//
//   constantNDelays   = (float*)malloc(1*nDir*nAnt*sizeof(float));
//   constantNDelays_V = (float*)malloc(1*nDir*(nAnt/2)*sizeof(float));
//   constantNDelays_H = (float*)malloc(1*nDir*(nAnt/2)*sizeof(float));
//   onion_temp = new Healpix_Onion(nSideExp, 1, 5000, 5000);
//
//   err = computeRecoDelaysWithNoBoundConstantN(nAnt, -1.f*stationCenterDepth, antLocation,
//                                               onion_temp, constantNDelays, constantNDelays_V, constantNDelays_H);
//
//}


//recordTime(tmr,2);
time_t t_after_recoDelays = time(NULL);
clock_t c_after_recoDelays = clock();

/*
 * Set up reco environment. In this case, an OpenCL environment
 */
recoEnvData clEnv;
err = setupCLRecoEnv(settings, &clEnv, settings->programFile/*.c_str()*/);
if( err<0 ){
   cerr<<"Error setting up reco env\n"; return -1;
}


   TGraph *gr_v_temp[16];
   TGraph *gr_v[16];
   TGraph *grInt[16];
   TGraph *grWinPad[16];
   TGraph *grNormWinPad[16];
   //TGraph *gr_fft[16];
   //TGraph *grHilbert[16];
   TGraph *grMean[16];
   TGraph *grFFT[16];
   TGraph *grCDF[16];
   TGraph *grCumuSum[16];
   TGraph *grCumuSumCDF[16];


if(settings->dataType == 1){
/*
   int cutWaveAlert;
   double addDelay;
   double times, volts;
   double time_1, time_2, time_last;

   int utime, utime_runStart, utime_runEnd;
*/
   eventTree->GetEntry(0);
   utime_runStart=utime_runEnd=rawAtriEvPtr->unixTime;

   if(rawAtriEvPtr->isRFTrigger()){
      if(rawAtriEvPtr->isCalpulserEvent()){ runCalEventCount++; }
      else { runRFEventCount++; }
   } else if (rawAtriEvPtr->isSoftwareTrigger()){ runSoftEventCount++; }
   else { cerr<<"Undefined trigger type!!\n"; /*continue;*/ }
/*
 * Loop over events once to determine run start/end time
 */

   for(int ev=1; ev<runEventCount/*numEntries*/; ev++){
      eventTree->GetEntry(ev);
      if(rawAtriEvPtr->unixTime < utime_runStart) utime_runStart=rawAtriEvPtr->unixTime;
      if(rawAtriEvPtr->unixTime > utime_runEnd  ) utime_runEnd  =rawAtriEvPtr->unixTime;
      if(rawAtriEvPtr->isRFTrigger()){
         if(rawAtriEvPtr->isCalpulserEvent()){ runCalEventCount++; }
         else { runRFEventCount++; }
      } else if (rawAtriEvPtr->isSoftwareTrigger()){ runSoftEventCount++; }
      else { cerr<<"Undefined trigger type!!\n"; continue; }
   }
   cout<<"utime_runStart: "<<utime_runStart<<" dropD4Time: "<<dropD4Time<<endl;
   cout<<"Run time span: "<<utime_runEnd-utime_runStart<<endl;

}//end of if dataType = 1
/*
 * Start looping events for analysis
 */
vector<TGraph *> cleanEvent;
//int recoEventCnt = 0;
int recoFlagCnt = 0;
double t, v, beginTime = 0.;

int *maxPix;
if(settings->skymapSearchMode == 0) maxPix = (int*)calloc(nDir*nLayer, sizeof(int));
else                                maxPix = (int*)calloc(12*pow(2,settings->nSideExpEnd)*pow(2,settings->nSideExpEnd)*settings->nLayer, sizeof(int));

int maxPixIdx = 0;
int maxPixIdx2 = 0;
int constantNMaxPixIdx = 0;
float *mapData = (float*)calloc(nDir*nLayer, sizeof(float));
char histName[200];
TH1F **mapDataHist = (TH1F**)malloc(nDir*nLayer*sizeof(TH1F*));
if(settings->recordMapData == 1){
   //char histName[200];
   //TH1F *mapDataHist[nDir*nLayer];
   for(int pix=0; pix<nDir*nLayer; pix++){
      sprintf(histName, "pix_%d", pix);
      mapDataHist[pix] = new TH1F(histName, histName, 5000, 0.f, 1.f);
   }
}
int index[16]={0};
int index_V[8]={0};
int index_H[8]={0};
float snrArray[16], unmodSNRArray[16], snrArray_V[8], snrArray_H[8];
vector<TGraph *> unpaddedEvent;
TH1F *snrDist = new TH1F("snrDist","snrDist",100,0,50);
int goodChan[16];
int satChan[16];
int numSatChan;
if(settings->nchnlFilter > 0){
   threshold = settings->nchnlThreshold;
   //int nchnlArray[3];
   //int nchnl_tmp;
}
//int trigEvCnt = 0;
int triggerCode[3];
for(int i=0; i<3; i++){ triggerCode[i] = settings->triggerCode[i] - '0'; cout<<"triggerCode "<<i<<": "<<triggerCode[i]<<endl; }

cout<<"runEventCount: "<<runEventCount<<endl;

recoData *summary = new recoData();
dataTree->Branch("summary", &summary);

char histname[200];
//TH1F *dtHist[64];
//for(int b=0; b<64; b++){
//   sprintf(histname, "chan%d_%d", b/8, b%8);
//   dtHist[b] = new TH1F(histname, histname, 840*2, 0, 840);
//}


vector< vector<double> > fftValues[16];
double f, p;
int eventCount = 0;
vector<double>* fftValues_V;
vector<double>* fftValues_H;
int fInt_V, fInt_H;
//vector<double> vec_v;
//vector<double> vec_h;

//recordTime(tmr,3);
time_t t_before_event_loop = time(NULL);
clock_t c_before_event_loop = clock();

if(settings->dataType == 1){

trigEventCount = runEventCount;
for (Long64_t ev=0; ev<runEventCount; ev++){

   if(settings->maxNumberOfReco >= 0)
      if(recoEventCount == settings->maxNumberOfReco) break;

   summary->clear();

   if(ev%1000 == 0) cout<<"*******************************Event got********************************: "<<ev<<endl;

   eventTree->GetEntry(ev);

   if(settings->recoEventIndex > -1){ //check if only want to reconstruct a specified event
   if(rawAtriEvPtr->eventNumber != settings->recoEventIndex) continue;
   }


   //cout<<"Code loop ev: "<<ev<<" eventId: "<<rawAtriEvPtr->eventId<<" eventNumber: "<<rawAtriEvPtr->eventNumber<<endl;
   summary->setEventId(rawAtriEvPtr->eventId);
   summary->setEventNumber(rawAtriEvPtr->eventNumber);
   summary->setEventTime(rawAtriEvPtr->unixTime, rawAtriEvPtr->unixTimeUs, rawAtriEvPtr->timeStamp);


   if(rawAtriEvPtr->isRFTrigger()){
      if(rawAtriEvPtr->isCalpulserEvent()){
         summary->setEventTrigType( 1 );
         if(triggerCode[1] != 1) continue;
         } else { // RF trigger
         summary->setEventTrigType( 0 );
         if(triggerCode[0] != 1) continue;
         }
   } else if (rawAtriEvPtr->isSoftwareTrigger()){
      summary->setEventTrigType( 2 );
      if(triggerCode[2] != 1) continue;
   } else { cerr<<"Undefined trigger type!!\n"; continue; }



   UsefulAtriStationEvent *realAtriEvPtr = new UsefulAtriStationEvent( rawAtriEvPtr, AraCalType::kLatestCalib);
//*************APPLYING DELAYS. CODE FROM T. MEURES*****************
      int stationId = realAtriEvPtr->stationId;
      cutWaveAlert = 0;
      nonIncreasingSampleTimeAlert = 0;
      mistaggedSoftEventAlert = 1;
      previous_times = 0.;
      double stdDelay = 0.;

//***********CHECK ARA02 D1 CORRUPTION*******************************

   if( stationId == 2 && realAtriEvPtr->unixTime >= ARA02D1CorruptionStartTime && realAtriEvPtr->unixTime <= ARA02D1CorruptionEndTime){
      corruptD1EventCount+=1;
      delete realAtriEvPtr;
      continue;
   }

//***********CHECK FIRTS 3 EVENTS CORRUPTION*************************

   if (( stationId == 2 && realAtriEvPtr->unixTime >= corruptFirst3EventStartTime_A2 && realAtriEvPtr->eventNumber < corruptEventEndEventNumber)
   ||
   ( stationId == 3 && realAtriEvPtr->unixTime >= corruptFirst3EventStartTime_A3 && realAtriEvPtr->eventNumber < corruptEventEndEventNumber))
   {
      corruptFirst3EventCount+=1;
      delete realAtriEvPtr;
      continue;
   }

   double average[16]={0.};

   for(int a=0;a<16;a++)
   {//Loop the 16 channels
      //cout<<"*** Channel "<<a<<"***"<<endl;
      addDelay = 0.0;
      //*** Now we add the cable delays from the file and the forgotten antenna feedthrough (4,8,12 ns). ***//
	  if(a/4==0){addDelay+=(4.0  + delays[a%4][3]);}
	  if(a/4==1){addDelay+=(12.0 + delays[a%4][3]);}
	  if(a/4==2){addDelay+=(0.0  + delays[a%4][3]);}
	  if(a/4==3){addDelay+=(8.0  + delays[a%4][3]);}

      //stdDelay= geom->getStationInfo(stationId)->getCableDelay(a);
      //addDelay += stdDelay;
      if(settings->applyA2Ch6Correction) if(a==6) addDelay += 6.8; //Fit 2, ch6 delay+=6.87ns

	  //*** We put the waveform into a graph. ***//
	  gr_v_temp[a] = realAtriEvPtr->getGraphFromRFChan(a);
	  gr_v[a] = new TGraph();
	  //*** A few waveforms are inverted in station 3, which needs to be corrected. I wrote a small method to do that (see above). ***//
	  if(stationId==3 && (a==0||a==4||a==8)){/*cout<<"graph inverted"<<endl;*/invertGraph(gr_v_temp[a]);}
	  //*** The following is to avoid reading corrupted waveforms. ***//
	  //*** I encountered only a few of them, so maybe this is not ***//
	  //*** really neccessary anymore. *******************************//
	  if(gr_v_temp[a]->GetN()<5 ){ cerr<< "BAD EVENT: " << ev << " Channel: " << a << ", points: " << gr_v_temp[a]->GetN() << endl;cutWaveAlert=1; /*cutWaveEventCount++;*/ /*continue;*/}
     //cout<<"Nsamp: "<<gr_v_temp[a]->GetN()<<" Soft trig nSamp: "<<IRS2SamplePerBlock*maxSoftTriggerReadoutBlocks<<endl;
     if(gr_v_temp[a]->GetN()>=(IRS2SamplePerBlock*maxSoftTriggerReadoutBlocks)){ mistaggedSoftEventAlert=0; }
	  int pc = 0;
    gr_v_temp[a]->GetPoint(0, times, volts);
    previous_times = times;
	  //*** The first 20 samples can be corrupted. Therefore, we need to exclude them! ***//
      for(int p=0;p<gr_v_temp[a]->GetN();p++){

         gr_v_temp[a]->GetPoint(p, times, volts);

         //cout<<"a: "<<a<<" p: "<<p<<" times: "<<times<<endl;

         if(/*times*/(times - previous_times)>20.0)
         {
         if(stationId==3 && utime_runStart>=dropD4Time && (a%4==3)){
            gr_v[a]->SetPoint(pc, times-addDelay, 0.); //Drop 2014 ARA03 D4
            dropARA03D4 = true;
         }
         else if(stationId==2 && a==15){
            gr_v[a]->SetPoint(pc, times-addDelay, 0.);//Drop ARA02 D4BH
            dropARA02D4BH = true;
         }
         else {
            gr_v[a]->SetPoint(pc, times - addDelay, volts);
            average[a]+=volts;
         }
		 pc++;
         }
      }
      /*** Zero-mean the waveforms on a channel-by-channel basis ***/
      if( gr_v[a]->GetN() != 0){

         average[a]/=(double)gr_v[a]->GetN();

         gr_v[a]->GetPoint(0, times, volts);
         gr_v[a]->SetPoint(0, times, volts-average[a]);
         previous_times = times;

         for(pc=1; pc<gr_v[a]->GetN(); pc++){

            gr_v[a]->GetPoint(pc, times, volts);

            if( (times - previous_times) > 0.)
            gr_v[a]->SetPoint(pc, times, volts-average[a]);
            else {cerr<< "BAD EVENT Non-increasing sample time: " << event << " Channel: " << a << "this sample time: "<< times << "previous sample time: " << previous_times << endl;nonIncreasingSampleTimeAlert=1; /*nonIncreasingSampleTimeEventCount++; if(cutWaveAlert==1){cutWaveAndNonIncreasingEventCount++;}*/}

            previous_times = times;

         }//end of pc

       } else {cerr<< "BAD EVENT type 2: " << event << " Channel: " << a << ", original number of points: " << gr_v_temp[a]->GetN() << endl; /*if(cutWaveAlert!=1){ cutWaveEventCount++;}*/ cutWaveAlert=1; /*continue;*/
       /*
       for(int p=0; p<gr_v_temp[a]->GetN(); p++){
          gr_v_temp[a]->GetPoint(p, times, volts);
          cout<<"p: "<<p<<" times: "<<times<<" ";
       }
       cout<<endl;
       */
       }

/*
      average[a]/=(double)gr_v[a]->GetN();
      for(pc=0; pc<gr_v[a]->GetN(); pc++){
         gr_v[a]->GetPoint(pc, times, volts);
         gr_v[a]->SetPoint(pc, times, volts-average[a]);
      }
      cout<<"gr_v_temp N:"<<gr_v_temp[a]->GetN()<<endl;
*/
      delete gr_v_temp[a];
   }//End looping channels

   double wInt;
   int maxSamp;
   bool shouldSkip = false;
   if (cutWaveAlert == 1 && nonIncreasingSampleTimeAlert == 1){ cutWaveAndNonIncreasingEventCount++; shouldSkip = true; }
   if (cutWaveAlert == 1) { cerr<<"Event "<<realAtriEvPtr->eventNumber<<" discarded due to cutWaveAlert\n"; cutWaveEventCount++; shouldSkip = true; }
   if (nonIncreasingSampleTimeAlert == 1) { cerr<<"Event "<<realAtriEvPtr->eventNumber<<" discarded due to nonIncreasingSampleTimeAlert\n"; nonIncreasingSampleTimeEventCount++; shouldSkip = true; }
   if (mistaggedSoftEventAlert == 1) { cerr<<"Event "<<realAtriEvPtr->eventNumber<<" discarded due to mistaggedSoftEventAlert\n"; mistaggedSoftEventCount++; shouldSkip = true;}
   if (shouldSkip){
      delete realAtriEvPtr;
      for(int ch=0; ch<16; ch++) delete gr_v[ch];
      continue;
   }

   beginTime = 1e10;
   for(int ch=0; ch<16; ch++){

      gr_v[ch]->GetPoint(0,t,v);
      if( t<beginTime ) beginTime = t ;

   }

   nChanBelowThres = 0;
   maxTimeVec.clear();

   for(int ch=0; ch<16; ch++){

   goodChan[ch] = chanMask[ch];

   if(ch<8){wInt=/*0.4*/settings->wInt_V; maxSamp=/*2048*/settings->maxPaddedSample;}
   else{wInt=/*0.625*/settings->wInt_H; maxSamp=/*2048*/settings->maxPaddedSample;}

   /* Interpolate + apply windowing + zero-pad + equalize wf beginning  to maxSamp */
   //cout<<"N: "<<gr_v[ch]->GetN()<<endl;
   grInt[ch]       = FFTtools::getInterpolatedGraph(gr_v[ch], wInt);
   unpaddedEvent.push_back(grInt[ch]);
   /* Window type is passed in settings now */
   grWinPad[ch]     = evProcessTools::getWindowedAndPaddedEqualBeginGraph(settings->windowingType, grInt[ch], maxSamp, beginTime);
   /* The task of normalizing wf should be the responsibility of each reco method */
   cleanEvent.push_back(grWinPad[ch]);

   grMean[ch] = evProcessTools::getRollingMeanGraph(grInt[ch], IRS2SamplePerBlock);
   meanMax[ch] = evProcessTools::getMax(grMean[ch], &maxTime[ch]);
   if(meanMax[ch]<(ch<8?threshold_V:threshold_H)){
      nChanBelowThres += 1;
      maxTimeVec.push_back(maxTime[ch]);
   }

   delete gr_v[ch];
   }//end of ch

   /* Check for offset block */

   nChanBelowThres_Thres = (dropARA02D4BH?15:(dropARA03D4?12:16));
   if( nChanBelowThres >= nChanBelowThres_Thres ){
      timeRange = *max_element(maxTimeVec.begin(), maxTimeVec.end()) - *min_element(maxTimeVec.begin(), maxTimeVec.end());
      //cout<<"max element: "<<*max_element(maxTimeVec.begin(), maxTimeVec.end())<<" min element: "<<*min_element(maxTimeVec.begin(), maxTimeVec.end())<<" timeRange: "<<timeRange<<endl;
      if(timeRange < timeRangeCut){
         offsetBlockAlert = 1;
         offsetBlockEventCount += 1;
         unpaddedEvent.clear();
         cleanEvent.clear();
         delete realAtriEvPtr;
         for(int ch=0; ch<16; ch++){ delete grInt[ch]; delete grWinPad[ch]; delete grMean[ch]; }
         continue;
      }
   }

   numSatChan = getSaturation(settings, unpaddedEvent, satChan);
   summary->setSaturatedChannels(numSatChan, satChan);
   if(settings->maskSaturatedChannels) for(int ch=0; ch<16; ch++) goodChan[ch] = goodChan[ch] && (!satChan[ch]);




   double imp;
   double bipolarness;
   double posPowerPeak, negPowerPeak, powerPeaksDeltaT;
   int maxFracBin;
   double maxFrac;

   for(int ch=0; ch<16; ch++){

      /* Measure impulsivity */
      grCDF[ch] =  impulsivityMeasure(unpaddedEvent[ch], &imp);
      //summary->setImpulsivityByChannel(ch, impulsivityMeasure(unpaddedEvent[ch], NULL, NULL));
      summary->setImpulsivityByChannel(ch, imp);

      /* Measure bipolarness */
      grCumuSum[ch] = bipolarnessMeasure(unpaddedEvent[ch], &bipolarness/*, grCumuSumCDF[ch]*/);
      //cout<<"n: "<<grCumuSumCDF[ch]->GetN()<<endl;
      summary->setBipolarnessByChannel(ch, bipolarness);

      /* Measure +/- power peaks and dT */
      getPosNegPowerPeakAndDeltaT(unpaddedEvent[ch], &posPowerPeak, &negPowerPeak, &powerPeaksDeltaT);
      summary->setPowerPeaksByChannel(ch, posPowerPeak, negPowerPeak, powerPeaksDeltaT);

      /* Get max freq bin */
      grFFT[ch] = FFTtools::makePowerSpectrumMilliVoltsNanoSecondsdB(grWinPad[ch]);

      if(ch==0){
         freqCountLen_V = grFFT[ch]->GetN();
         freqCount_V = new int [freqCountLen_V];
         fill(&freqCount_V[0], &freqCount_V[freqCountLen_V], 0);
         freqBinWidth_V = evProcessTools::getFFTBinWidth(grFFT[ch]);
         if(eventCount == 0) fftValues_V = new vector<double> [freqCountLen_V*8];
         //if(eventCount == 0) fftValues_V = (vector<double>*)malloc(8*freqCountLen_V*sizeof(vector<double>));
         //fInt_V = freqBinWidth_V;
      }
      else if (ch==8){
         freqCountLen_H = grFFT[ch]->GetN();
         freqCount_H = new int [freqCountLen_H];
         fill(&freqCount_H[0], &freqCount_H[freqCountLen_H], 0);
         freqBinWidth_H = evProcessTools::getFFTBinWidth(grFFT[ch]);
         if(eventCount == 0) fftValues_H = new vector<double> [freqCountLen_H*8];
         //if(eventCount == 0) fftValues_H = (vector<double>*)malloc(8*freqCountLen_H*sizeof(vector<double>));
      }
      //cout<<"926\n";
      //cout<<"fftValues_V size: "<<sizeof(fftValues_V)<<endl;
      //cout<<"vector<double> size: "<<sizeof(vector<double>)<<endl;
      //cout<<"ch: "<<ch<<endl;
      //cout<<"freqCountLen_V: "<<freqCountLen_V<<endl;
      if(ch<8){
         for(int bin=0; bin<freqCountLen_V; bin++){
            //cout<<"929 bin:"<<bin<<"\n";
            grFFT[ch]->GetPoint(bin, f, p);
            //cout<<"931\n";
            fftValues_V[ch*freqCountLen_V+bin].push_back(p);

            //cout<<"933\n";

         }
      }
       else {

         for(int bin=0; bin<freqCountLen_H; bin++){
            grFFT[ch]->GetPoint(bin, f, p);
            fftValues_H[(ch-8)*freqCountLen_H+bin].push_back(p);
         }
         //fftValues[ch][bin].push_back(p);
      }


      maxFrac = FFTtools::getPeakVal(grFFT[ch], &maxFracBin);
      //cout<<"ch: "<<ch<<" macFrac: "<<maxFrac<<" maxFracBin: "<<maxFracBin<<endl;
      summary->setMaxFreqBinByChannel(ch, maxFracBin, maxFrac);
      if(ch<8)
      freqCount_V[maxFracBin]++;
      else
      freqCount_H[maxFracBin]++;

      delete grCDF[ch];
      delete grCumuSum[ch];
      //delete grCumuSumCDF[ch];
      delete grFFT[ch];

   }

   eventCount++;
   //out<<"eventCount: "<<eventCount<<endl;

   summary->setFreqBinWidth(freqBinWidth_V, freqBinWidth_H);

   computeSNR(settings, unpaddedEvent, summary);

   std::fill(&snrArray[0], &snrArray[16], 0.);
   if(settings->snrMode==0) memcpy(snrArray, summary->channelInWindowSNR, sizeof(summary->channelInWindowSNR)); //snrArray[ch] = summary->inWindowSNR[ch];
   else if (settings->snrMode==1) memcpy(snrArray, summary->slidingV2SNR, sizeof(summary->slidingV2SNR)); //snrArray[ch] = summary->slidingV2SNR[ch];
   else if (settings->snrMode==2) memcpy(snrArray, summary->totalPowerSNR, sizeof(summary->totalPowerSNR)); //snrArray[ch] = summary->totalPowerSNR[ch];
   else { cerr<<"Invalid snrMode: "<<settings->snrMode<<endl; return -1; }

   //for(int ch=0; ch<16; ch++){
   //   cout<<"snrArray: "<<snrArray[ch]<<" channelInWindowSNR: "<<summary->channelInWindowSNR[ch]<<" slidingV2SNR: "<<summary->slidingV2SNR[ch]<<" totalPowerSNR: "<<summary->totalPowerSNR[ch]<<endl;
   //}

   for(int ch=0; ch<8; ch++){
      snrArray_V[ch] = snrArray[ch];
      snrArray_H[ch] = snrArray[ch+8];
   }

   TMath::Sort(16, snrArray, index);
   TMath::Sort(8, snrArray_V, index_V);
   TMath::Sort(8, snrArray_H, index_H);

   snrRank = (string(settings->recoPolType)=="both" ? index : (string(settings->recoPolType)=="vpol" ? index_V : index_H));

   //****************************************************
   // FILTER SECTION
   //****************************************************

   int nonZeroChan = 0;
   double avgImp = 0.;

   if(settings->impulsivityFilter > 0){

      for(int ch=(string(settings->recoPolType)=="hpol"?8:0); ch<(string(settings->recoPolType)=="vpol"?8:16); ch++){
         if(fabs(summary->impulsivity[ch]-0)>1e-9){
            nonZeroChan += 1;
            avgImp += summary->impulsivity[ch];
         }
      }

      if(nonZeroChan>0) avgImp /= (double)nonZeroChan;
      else avgImp = 0.;

      if(avgImp < settings->impulsivityThreshold){

         impulsivityFilteredEventCount+=1;
         unpaddedEvent.clear();
         cleanEvent.clear();
         delete realAtriEvPtr;
         for(int ch=0; ch<16; ch++){ delete grInt[ch]; delete grWinPad[ch]; delete grMean[ch]; /*delete grCDF[ch];*/}
         continue;
      }

   }

   /* Nchnl filter */

   if(settings->nchnlFilter > 0){

      float selectedSNR;
      if(settings->nchnlFilter == 3)      selectedSNR = snrArray[index[settings->nchnlCut-1]];
      else if(settings->nchnlFilter == 1) selectedSNR = snrArray_V[index_V[settings->nchnlCut-1]];
      else                                selectedSNR = snrArray_H[index_H[settings->nchnlCut-1]];

      //if(nchnl_tmp < settings->nchnlCut){
      if( selectedSNR < threshold )
      {
         nchnlFilteredEventCount+=1;
         //cerr<<"Failed nchnl cut. nchnl_tmp: "<<nchnl_tmp<<endl;
         unpaddedEvent.clear();
         cleanEvent.clear();
         delete realAtriEvPtr;
         for(int ch=0; ch<16; ch++){ delete grInt[ch]; delete grWinPad[ch]; delete grMean[ch]; /*delete grCDF[ch];*/}
         continue;
      }
      else {
         //if the event passes the filter then we consider whether we want to mask sub-threshold channels
         if(settings->maskSubThresholdChannels){
            for(int ch=0; ch<16; ch++){
               if(snrArray[ch] < threshold) goodChan[ch] = 0;
            }
         }

         //if the events passes the filter, we consider whether it passes another pol as well
         if(settings->nchnlFilter==1){
            if(snrArray_H[index_H[settings->nchnlCut-1]]>=settings->nchnlThreshold_anotherPol) summary->setPassAnotherPolNchnl(true);
         }
         else if(settings->nchnlFilter==2){
            if(snrArray_V[index_V[settings->nchnlCut-1]]>=settings->nchnlThreshold_anotherPol) summary->setPassAnotherPolNchnl(true);
         }

      }
   }//end if if nchnlFilter=1


   /* CW filter */

   //if(settings->cwFilter > 0){

   bool isCW = false;
/*
   for(int ch=0; ch<16; ch++){

      grFFT[ch] = FFTtools::makePowerSpectrumMilliVoltsNanoSecondsdB(grWinPad[ch]);

      if(ch==0){
         freqCountLen_V = grFFT[ch]->GetN();
         freqCount_V = new int [freqCountLen_V];
         fill(&freqCount_V[0], &freqCount_V[freqCountLen_V], 0);
         freqBinWidth_V = evProcessTools::getFFTBinWidth(grFFT[ch]);
      }
      else if (ch==8){
         freqCountLen_H = grFFT[ch]->GetN();
         freqCount_H = new int [freqCountLen_H];
         fill(&freqCount_H[0], &freqCount_H[freqCountLen_H], 0);
         freqBinWidth_H = evProcessTools::getFFTBinWidth(grFFT[ch]);
      }

      int maxFracBin;
      double maxFrac = FFTtools::getPeakVal(grFFT[ch], &maxFracBin);
      //cout<<"ch: "<<ch<<" macFrac: "<<maxFrac<<" maxFracBin: "<<maxFracBin<<endl;
      summary->setMaxFreqBinByChannel(ch, maxFracBin, maxFrac);
      if(ch<8)
      freqCount_V[maxFracBin]++;
      else
      freqCount_H[maxFracBin]++;

   }//end of ch

   summary->setFreqBinWidth(freqBinWidth_V, freqBinWidth_H);
*/
   maxCount_V = evProcessTools::getMaxCount(freqCountLen_V, freqCount_V, &maxCountBin, 2); // a coincidence of >= 2 counts as potentially meaningful
   maxCountFreq_V = freqBinWidth_V * maxCountBin;

   maxCount_H = evProcessTools::getMaxCount(freqCountLen_H, freqCount_H, &maxCountBin, 2);
   maxCountFreq_H = freqBinWidth_H * maxCountBin;

   summary->setMaxCountFreq(maxCountFreq_V, maxCountFreq_H);

   /* If coincidence >= minCWCoincidence */
   if(maxCount_V >= minCWCoincidence || maxCount_H >= minCWCoincidence) isCW = true;

   /* If coincidence = minCWCoincidence - 1 and there is at least one neighboring peak bin */

   //Check Vpol
   if (!isCW){
      if( ((freqCount_V[0] == minCWCoincidence-1) && (freqCount_V[1] > 0) ) ||
      ((freqCount_V[freqCountLen_V-1] == minCWCoincidence-1) && (freqCount_V[freqCountLen_V-2] > 0))
   ) isCW = true;
   }
   if (!isCW){
      for(int i=1; i<freqCountLen_V-1; i++){
         if(freqCount_V[i]==minCWCoincidence-1 && (freqCount_V[i-1]>0 || freqCount_V[i+1]>0)){
            isCW = true;
            break;
         }
      }
   }

   //Check Hpol
   if (!isCW){
      if( ((freqCount_H[0] == minCWCoincidence-1) && (freqCount_H[1] > 0) ) ||
      ((freqCount_H[freqCountLen_H-1] == minCWCoincidence-1) && (freqCount_H[freqCountLen_H-2] > 0))
   ) isCW = true;
   }
   if (!isCW){
      for(int i=1; i<freqCountLen_H-1; i++){
         if(freqCount_H[i]==minCWCoincidence-1 && (freqCount_H[i-1]>0 || freqCount_H[i+1]>0)){
            isCW = true;
            break;
         }
      }
   }

   if(settings->cwFilter>0){
      if(isCW){

         cwFilteredEventCount+=1;
         unpaddedEvent.clear();
         cleanEvent.clear();
         delete realAtriEvPtr;
         for(int ch=0; ch<16; ch++){ delete grInt[ch]; delete grWinPad[ch]; delete grMean[ch]; /*delete grFFT[ch];*/ /*delete grCDF[ch];*/ }
         continue;

      }
   }

   /* Constant N reco for surface filter */

//   recoSuccess = false;
//
//   if(settings->constantNFilter > 0){
//      while( !recoSuccess ){
//
//         stringstream ss;
//         ss << /*ev*/rawAtriEvPtr->eventNumber;
//         evStr = ss.str();
//         fitsFileStr = fitsFile_tmp /*+ ".ev" + evStr*/ + ".constantN.fits";
//         sprintf(fitsFile, fitsFileStr.c_str());
//         constantNMaxPixIdx = reconstruct3DXCorrEnvelopeGetMaxPixAndMapData_constantNFilter(settings, cleanEvent, &clEnv, constantNDelays, constantNDelays_V, constantNDelays_H, goodChan, summary, fitsFile///*argv[5]*/, mapData/*, xCorrAroundPeakHist, sillygr*/
//         );
//
//         if( constantNMaxPixIdx < 0){ cerr<<"Error reconstructing - contant N\n"; return -1; }
//         if(summary->constantNMaxPixCoherence != 0.f) recoSuccess = true; //To catch cases where GPU reco returns coherence value zero
//         else { cout<<"constantNMaxPixCoherence returns 0!! Re-running reco...\n"; }
//
//      }
//
//      if(recordConstantNDir(settings, summary) < 0){ cerr<<"Error recording constant N dir\n"; return -1;}
//
//   }
//
//   /* Track engine object to compute all tracks */
//   treg->computeAllTracks(unpaddedEvent);
//   summary->setTreg(treg);

   //****************************************************
   // END OF FILTER SECTION
   //****************************************************


/*
    if(settings->dataType == 0){
    summary->setWeight(event->Nu_Interaction[0].weight);
    summary->setTrueRadius(r_true);
    summary->setTrueDir(zen_true*180./M_PI, azi_true*180./M_PI);
    }
*/

   /* Real data, set as 1 */
   summary->setWeight(1);
   summary->setProbabilities(1,1);
   summary->setSurvivalProbability(1);

    //summary->setOnion(onion);
   summary->setTopN(topN);
   summary->setRecoChan(goodChan);

   cout<<"*********************************** inWindowSNR_V: "<<summary->inWindowSNR_V<<"*************************************"<<endl;
   recoEventCount++;

   /* Reco with radiospline */

//   recoSuccess = false;
//
//   while( !recoSuccess ){
//   if(settings->beamformMethod == 1){
//   if(settings->getSkymapMode == 0){
//       err = reconstructCSW(settings, cleanEvent, &clEnv, recoDelays, recoDelays_V, recoDelays_H, nDir, chanMask, fitsFile/*argv[5]*/);
//   }
//   else{
//       maxPixIdx = reconstructCSWGetMaxPix(settings, cleanEvent, &clEnv, recoDelays, recoDelays_V, recoDelays_H, nDir, chanMask, summary);
//   }
//   } else {
//   if(settings->getSkymapMode == 0){
//
//      //evStr = std::to_string(ev);
//      stringstream ss;
//      ss << /*ev*/rawAtriEvPtr->eventNumber;
//      evStr = ss.str();
//      fitsFileStr = fitsFile_tmp /*+ ".ev" + evStr*/ + ".fits";
//      sprintf(fitsFile, fitsFileStr.c_str());
//
//      if(settings->skymapSearchMode == 0){ //no zoom mode
//      maxPixIdx = reconstruct3DXCorrEnvelopeGetMaxPixAndMapData(settings, cleanEvent, &clEnv, recoDelays, recoDelays_V, recoDelays_H, goodChan, snrRank, summary, fitsFile/*argv[5]*/, mapData/*, xCorrAroundPeakHist, sillygr*//*, calRecoDelays, dtHist*/);
//
//      if(settings->use2ndRayReco){
//      fitsFileStr = fitsFile_tmp /*+ ".ev" + evStr*/ + ".2ndRay.fits";
//      sprintf(fitsFile, fitsFileStr.c_str());
//      maxPixIdx2 = reconstruct3DXCorrEnvelopeGetMaxPixAndMapData_2ndRayReco(settings, cleanEvent, &clEnv, recoRefracDelays, recoRefracDelays_V, recoRefracDelays_H, goodChan, summary, fitsFile/*argv[5]*/, mapData/*, xCorrAroundPeakHist, sillygr*/);
//      }
//
//      if(settings->recordMapData == 1){
//      for(int pix=0; pix<nDir*nLayer; pix++) mapDataHist[pix]->Fill(mapData[pix]);
//      }
//      } else { //zoom search mode
//      maxPixIdx = reconstruct3DXCorrEnvelopeGetMaxPix_ZoomMode(settings, cleanEvent, &clEnv, stationCenterDepth, antLocation, recoDelays, recoDelays_V, recoDelays_H, goodChan, summary, fitsFile);
//      }
//   } else {
//      maxPixIdx = reconstructXCorrEnvelopeGetMaxPix(settings, cleanEvent, &clEnv, recoDelays, recoDelays_V, recoDelays_H, nDir, goodChan/*chanMask*/, summary);
//
//   }
//
//   }
//   if( err<0 || maxPixIdx<0 || maxPixIdx2<0 ){ cerr<<"Error reconstructing\n"; return -1; }
//   if(summary->maxPixCoherence != 0.f) recoSuccess = true; //To catch cases where GPU reco returns coherence value zero
//   if(settings->use2ndRayReco==1){ if(summary->maxPixCoherence2 == 0.f) recoSuccess=false;}
//   if(recoSuccess!=true) { cout<<"maxPixCoherence/2 returns 0!! Re-running reco...\n"; }
//   }//end of while
//
//   //int recoFlag = record3DDiffGetFlag(summary, outputFile);
//   //if( recoFlag ) recoFlagCnt++;
//   if(settings->use2ndRayReco==0){
//   summary->setFlag( (settings->skymapSearchMode)
//                    ? record3DZoomedDiffGetFlag(settings, summary, dZenDist, dAziDist, recoTrueZenDist, recoTrueAziDist)
//                    : record3DDiffGetFlag(settings, summary, dZenDist, dAziDist, recoTrueZenDist, recoTrueAziDist) );
//   } else {
//   if(settings->skymapSearchMode){ cerr<<"skymapSearchMode==1 & use2ndRayReco==1 incompatible!\n"; return -1;}
//   else summary->setFlag(record3DDiffGetFlag_2ndRayReco(settings, summary, dZenDist, dAziDist, recoTrueZenDist, recoTrueAziDist));
//   }
//   if(summary->flag > 0) recoFlagCnt++;
//   maxPix[maxPixIdx]++;
//
//   if(settings->use2ndRayReco==0)
//   compute3DRecoAnglesWithRadioSplineForSinglePixel(nAnt, -1.f*stationCenterDepth, antLocation, onion, recoLauAngles, recoRecAngles, maxPixIdx);
//   else{
//      if(summary->maxPixCoherence >= summary->maxPixCoherence2)
//      compute3DRecoAnglesWithRadioSplineForSinglePixel(nAnt, -1.f*stationCenterDepth, antLocation, onion, recoLauAngles, recoRecAngles, maxPixIdx);
//      else {
//      cerr<<"2nd ray angles table not yet build! Angles will be default (-1)\n"; //return -1;
//      //compute3DRecoAnglesWithRadioSplineForSinglePixel_2ndRayReco(nAnt, -1.f*stationCenterDepth, antLocation, onion, recoLauAngles, recoRecAngles, maxPixIdx2);
//      }
//   }
//   summary->setRecoAngles(recoRecAngles, recoLauAngles);


   dataTree->Fill();

   unpaddedEvent.clear();
   cleanEvent.clear();
   delete realAtriEvPtr;
   //delete summary;

//   treg->clearForNextEvent();

   for(int ch=0; ch<16; ch++){ delete grInt[ch]; delete grWinPad[ch]; delete grMean[ch]; /*delete grCDF[ch];*/ /*if(settings->cwFilter>0)*/ /*delete grFFT[ch];*/ }

   }//end of ev loop

   fp->Close();

} //if dataType == 1

else {

for (Long64_t ev=0; ev<runEventCount/*numEntries*/; ev++){
   //cout<<"Entering event loop\n";

   if(settings->maxNumberOfReco >= 0)
      if(recoEventCount == settings->maxNumberOfReco) break;

   summary->clear();
   //cout<<"Cleared previous summary\n";
   if(ev%1000 == 0) cout<<"*******************************Event got********************************: "<<ev<<endl;

   if(settings->recoEventIndex > -1){ //check if only want to reconstruct a specified event
   if(ev != settings->recoEventIndex) continue;
   }

   summary->setEventNumber(ev);

   int string_i, antenna_i, AraRootChannel;

   chain2.GetEntry(ev);
   if(report->stations[0].Global_Pass > 0){

      summary->setSurvivalProbability(event->Nu_Interaction[0].weight);
      double L_ice = event->Nu_Interaction[0].len_int_kgm2_total/Signal::RHOICE;
      double p_int = 1.-exp(-1.*(event->Nu_Interaction[0].r_enterice.Distance(event->Nu_Interaction[0].nuexitice)/L_ice)); // probability it interacts in ice along its path
      summary->setProbabilities(p_int, p_int*event->Nu_Interaction[0].weight);
      weight = computeWeight(AraSim_settings, detector, event, icemodel, -1.*stationCenterDepth, L_ice);
      if(weight<0){ cerr<<"Weight compuation error!\n"; return -1; }
      summary->setWeight(weight);

      trigEventCount++;
      weightedTrigEventCount+=weight;
      dx = event->Nu_Interaction[0].posnu.GetX()-detector->stations[0].GetX();
      dy = event->Nu_Interaction[0].posnu.GetY()-detector->stations[0].GetY();
      dz = event->Nu_Interaction[0].posnu.GetZ()-detector->stations[0].GetZ() + stationCenterDepth;
      r_xy   = sqrt( dx*dx + dy*dy );
      r_true = sqrt( dx*dx + dy*dy + dz*dz );
      zen_true = atan( r_xy / dz );
      if(dz<0) zen_true += M_PI;
      azi_true = atan( dy / dx );
      if(dx<0) azi_true += M_PI;
      else if (dy<0) azi_true += 2*M_PI;
      //cout<<"r_true: "<<r_true<<" zen_true: "<<zen_true*180./M_PI<<" azi_true: "<<azi_true*180./M_PI<<endl;

   } else continue;

   double average[16]={0.};

   for(int a=0;a<16;a++)
   {
      string_i  = detector->getStringfromArbAntID(0, a);
      antenna_i = detector->getAntennafromArbAntID(0, a);
      AraRootChannel = detector->GetChannelfromStringAntenna (0, string_i, antenna_i, AraSim_settings);
      gr_v[AraRootChannel-1] = new TGraph();
      int nSamp = (int)report->stations[0].strings[string_i].antennas[antenna_i].time_mimic.size();

      if(AraSim_settings->DETECTOR_STATION==2 && AraRootChannel-1==15 && settings->dropARA02D4BH==1)//ARA02: always drop D4BH
      {
         dropARA02D4BH = true;
         for(int s=64; s<nSamp; s++){
            gr_v[AraRootChannel-1]->SetPoint(s-64, report->stations[0].strings[string_i].antennas[antenna_i].time_mimic[s], 0.);
         }
      }
      else if(AraSim_settings->DETECTOR_STATION==3 && AraRootChannel%4==0 && settings->dropARA03D4==1)
      {
         dropARA03D4 = true;
         for(int s=64; s<nSamp; s++){
            gr_v[AraRootChannel-1]->SetPoint(s-64, report->stations[0].strings[string_i].antennas[antenna_i].time_mimic[s], 0.);
         }
      }
      else
      {
         for(int s=64; s<nSamp; s++){
            gr_v[AraRootChannel-1]->SetPoint(s-64
                                          , report->stations[0].strings[string_i].antennas[antenna_i].time_mimic[s]
                                          , report->stations[0].strings[string_i].antennas[antenna_i].V_mimic[s]
                                             );
            average[AraRootChannel-1] += report->stations[0].strings[string_i].antennas[antenna_i].V_mimic[s];
        }
      }

      int pc = 0;
      /*** Zero-mean the waveforms on a channel-by-channel basis ***/
      if(gr_v[AraRootChannel-1]->GetN() != 0) average[AraRootChannel-1]/=(double)gr_v[AraRootChannel-1]->GetN();
      else cout<<"********Zero samples!!! *** ch:"<<AraRootChannel-1<<" n_samp: "<<gr_v[AraRootChannel-1]->GetN()<<endl;

      for(pc=0; pc<gr_v[AraRootChannel-1]->GetN(); pc++){
         gr_v[AraRootChannel-1]->GetPoint(pc, times, volts);
         gr_v[AraRootChannel-1]->SetPoint(pc, times, volts-average[AraRootChannel-1]);
      }

      if(settings->flattenSaturatedAmplitude){
         for(pc=0; pc<gr_v[AraRootChannel-1]->GetN(); pc++){
            gr_v[AraRootChannel-1]->GetPoint(pc, times, volts);
            if( volts > settings->saturationVoltage_mV ) volts = settings->saturationVoltage_mV;
            else if( volts < -1.*settings->saturationVoltage_mV ) volts = -1.*settings->saturationVoltage_mV;
            gr_v[AraRootChannel-1]->SetPoint(pc, times, volts);
         }
      }

   }//End looping channels

   double wInt;
   int maxSamp;

   beginTime = 1e10;
   for(int ch=0; ch<16; ch++){

      gr_v[ch]->GetPoint(0,t,v);
      if( t<beginTime ) beginTime = t;

   }

   nChanBelowThres = 0;
   maxTimeVec.clear();

   for(int ch=0; ch<16; ch++){

   goodChan[ch] = chanMask[ch];

   if(ch<8){wInt=/*0.4*/settings->wInt_V; maxSamp=/*2048*/settings->maxPaddedSample;/* maxSamp = gr_v[0]->GetN()*4;*/ }
   else{    wInt=/*0.625*/settings->wInt_H; maxSamp=/*2048*/settings->maxPaddedSample;/* maxSamp = gr_v[8]->GetN()*4;*/ }
   grInt[ch]       = FFTtools::getInterpolatedGraph(gr_v[ch], wInt);
   unpaddedEvent.push_back(/*gr_v[ch]*/grInt[ch]);
   /* Window type is passed in settings now */
   grWinPad[ch]     = evProcessTools::getWindowedAndPaddedEqualBeginGraph(settings->windowingType, /*gr_v[ch]*/grInt[ch], maxSamp, beginTime);
   /* The task of normalizing wf should be the responsibility of each reco method */
   cleanEvent.push_back(grWinPad[ch]);

   grMean[ch] = evProcessTools::getRollingMeanGraph(grInt[ch], IRS2SamplePerBlock);
   meanMax[ch] = evProcessTools::getMax(grMean[ch], &maxTime[ch]);
   if(meanMax[ch]<(ch<8?threshold_V:threshold_H)){
      nChanBelowThres += 1;
      maxTimeVec.push_back(maxTime[ch]);
   }

   delete gr_v[ch];
   }//end of ch

   /* Check for offset block */

   nChanBelowThres_Thres = (dropARA02D4BH?15:(dropARA03D4?12:16));
   if( nChanBelowThres >= nChanBelowThres_Thres ){
      timeRange = *max_element(maxTimeVec.begin(), maxTimeVec.end()) - *min_element(maxTimeVec.begin(), maxTimeVec.end());
      //cout<<"max element: "<<*max_element(maxTimeVec.begin(), maxTimeVec.end())<<" min element: "<<*min_element(maxTimeVec.begin(), maxTimeVec.end())<<" timeRange: "<<timeRange<<endl;
      if(timeRange < timeRangeCut){
         offsetBlockAlert = 1;
         offsetBlockEventCount += 1;
         weightedOffsetBlockEventCount += weight;
         unpaddedEvent.clear();
         cleanEvent.clear();
         for(int ch=0; ch<16; ch++){ delete grInt[ch]; delete grWinPad[ch]; delete grMean[ch]; }
         continue;
      }
   }

   numSatChan = getSaturation(settings, unpaddedEvent, satChan);
   summary->setSaturatedChannels(numSatChan, satChan);
   if(settings->maskSaturatedChannels) for(int ch=0; ch<16; ch++) goodChan[ch] = goodChan[ch] && (!satChan[ch]);

   double imp;
   double bipolarness;
   double posPowerPeak, negPowerPeak, powerPeaksDeltaT;
   int maxFracBin;
   double maxFrac;

   for(int ch=0; ch<16; ch++){

      /* Measure impulsivity */
      grCDF[ch] =  impulsivityMeasure(unpaddedEvent[ch], &imp);
      //summary->setImpulsivityByChannel(ch, impulsivityMeasure(unpaddedEvent[ch], NULL, NULL));
      summary->setImpulsivityByChannel(ch, imp);

      /* Measure bipolarness */
      grCumuSum[ch] = bipolarnessMeasure(unpaddedEvent[ch], &bipolarness/*, grCumuSumCDF[ch]*/);
      summary->setBipolarnessByChannel(ch, bipolarness);

      /* Measure +/- power peaks and dT */
      getPosNegPowerPeakAndDeltaT(unpaddedEvent[ch], &posPowerPeak, &negPowerPeak, &powerPeaksDeltaT);
      summary->setPowerPeaksByChannel(ch, posPowerPeak, negPowerPeak, powerPeaksDeltaT);

      /* Get max freq bin */
      grFFT[ch] = FFTtools::makePowerSpectrumMilliVoltsNanoSecondsdB(grWinPad[ch]);

      if(ch==0){
         freqCountLen_V = grFFT[ch]->GetN();
         freqCount_V = new int [freqCountLen_V];
         fill(&freqCount_V[0], &freqCount_V[freqCountLen_V], 0);
         freqBinWidth_V = evProcessTools::getFFTBinWidth(grFFT[ch]);
         if(eventCount == 0) fftValues_V = new vector<double> [freqCountLen_V*8];
         //if(eventCount == 0) fftValues_V = (vector<double>*)malloc(8*freqCountLen_V*sizeof(vector<double>));
      }
      else if (ch==8){
         freqCountLen_H = grFFT[ch]->GetN();
         freqCount_H = new int [freqCountLen_H];
         fill(&freqCount_H[0], &freqCount_H[freqCountLen_H], 0);
         freqBinWidth_H = evProcessTools::getFFTBinWidth(grFFT[ch]);
         if(eventCount == 0) fftValues_H = new vector<double> [freqCountLen_H*8];
         //if(eventCount == 0) fftValues_H = (vector<double>*)malloc(8*freqCountLen_H*sizeof(vector<double>));
      }

      if(ch<8){
         for(int bin=0; bin<freqCountLen_V; bin++){
            //cout<<"929 bin:"<<bin<<"\n";
            grFFT[ch]->GetPoint(bin, f, p);
            //cout<<"931\n";
            fftValues_V[ch*freqCountLen_V+bin].push_back(p);

            //cout<<"933\n";

         }
      }
       else {

         for(int bin=0; bin<freqCountLen_H; bin++){
            grFFT[ch]->GetPoint(bin, f, p);
            fftValues_H[(ch-8)*freqCountLen_H+bin].push_back(p);
         }
         //fftValues[ch][bin].push_back(p);
      }

      eventCount++;

      maxFrac = FFTtools::getPeakVal(grFFT[ch], &maxFracBin);
      //cout<<"ch: "<<ch<<" macFrac: "<<maxFrac<<" maxFracBin: "<<maxFracBin<<endl;
      summary->setMaxFreqBinByChannel(ch, maxFracBin, maxFrac);
      if(ch<8)
      freqCount_V[maxFracBin]++;
      else
      freqCount_H[maxFracBin]++;

      delete grCDF[ch];
      delete grCumuSum[ch];
      //delete grCumuSumCDF[ch];
      delete grFFT[ch];

   }

   summary->setFreqBinWidth(freqBinWidth_V, freqBinWidth_H);

   computeSNR(settings, unpaddedEvent, summary);

   std::fill(&snrArray[0], &snrArray[16], 0.);
   if(settings->snrMode==0) memcpy(snrArray, summary->channelInWindowSNR, sizeof(summary->channelInWindowSNR)); //snrArray[ch] = summary->inWindowSNR[ch];
   else if (settings->snrMode==1) memcpy(snrArray, summary->slidingV2SNR, sizeof(summary->slidingV2SNR)); //snrArray[ch] = summary->slidingV2SNR[ch];
   else if (settings->snrMode==2) memcpy(snrArray, summary->totalPowerSNR, sizeof(summary->totalPowerSNR)); //snrArray[ch] = summary->totalPowerSNR[ch];
   else { cerr<<"Invalid snrMode: "<<settings->snrMode<<endl; return -1; }

   for(int ch=0; ch<8; ch++){
      snrArray_V[ch] = snrArray[ch];
      snrArray_H[ch] = snrArray[ch+8];
   }

   TMath::Sort(16, snrArray, index);
   TMath::Sort(8, snrArray_V, index_V);
   TMath::Sort(8, snrArray_H, index_H);

   snrRank = (string(settings->recoPolType)=="both" ? index : (string(settings->recoPolType)=="vpol" ? index_V : index_H));

   //****************************************************
   // FILTER SECTION
   //****************************************************

   int nonZeroChan = 0;
   double avgImp = 0.;

   /* Measure impulsivity */
/*
   for(int ch=0; ch<16; ch++){
      double imp;
      grCDF[ch] = impulsivityMeasure(unpaddedEvent[ch], &imp);
      //summary->setImpulsivityByChannel(ch, impulsivityMeasure(unpaddedEvent[ch], NULL, NULL));
      summary->setImpulsivityByChannel(ch, imp);
   }
*/
   if(settings->impulsivityFilter > 0){

      for(int ch=(string(settings->recoPolType)=="hpol"?8:0); ch<(string(settings->recoPolType)=="vpol"?8:16); ch++){
         if(fabs(summary->impulsivity[ch]-0)>1e-9){
            nonZeroChan += 1;
            avgImp += summary->impulsivity[ch];
         }
      }

      if(nonZeroChan>0) avgImp /= (double)nonZeroChan;
      else avgImp = 0.;

      if(avgImp < settings->impulsivityThreshold){

         impulsivityFilteredEventCount+=1;
         weightedImpulsivityFilteredEventCount += weight;
         unpaddedEvent.clear();
         cleanEvent.clear();
         for(int ch=0; ch<16; ch++){ delete grInt[ch]; delete grWinPad[ch]; delete grMean[ch]; /*delete grCDF[ch];*/}
         continue;
      }

   }

   /* Nchnl filter */
   //numSatChan = 0;

   if(settings->nchnlFilter > 0){

      float selectedSNR;
      if(settings->nchnlFilter == 3)      selectedSNR = snrArray[index[settings->nchnlCut-1]];
      else if(settings->nchnlFilter == 1) selectedSNR = snrArray_V[index_V[settings->nchnlCut-1]];
      else                                selectedSNR = snrArray_H[index_H[settings->nchnlCut-1]];

      //if(nchnl_tmp < settings->nchnlCut){
      if( selectedSNR < threshold )
      {
         nchnlFilteredEventCount+=1;
         weightedNchnlFilteredEventCount += weight;
         //cerr<<"Failed nchnl cut. nchnl_tmp: "<<nchnl_tmp<<endl;
         unpaddedEvent.clear();
         cleanEvent.clear();
         for(int ch=0; ch<16; ch++){ delete grInt[ch]; delete grWinPad[ch]; delete grMean[ch]; /*delete grCDF[ch];*/}
         continue;
      }
      else {
         //if the event passes the filter then we consider whether we want to mask sub-threshold channels
         if(settings->maskSubThresholdChannels){
            for(int ch=0; ch<16; ch++){
               if(snrArray[ch] < threshold) goodChan[ch] = 0;
            }
         }

         //if the events passes the filter, we consider whether it passes another pol as well
         if(settings->nchnlFilter==1){
            if(snrArray_H[index_H[settings->nchnlCut-1]]>=settings->nchnlThreshold_anotherPol) summary->setPassAnotherPolNchnl(true);
         }
         else if(settings->nchnlFilter==2){
            if(snrArray_V[index_V[settings->nchnlCut-1]]>=settings->nchnlThreshold_anotherPol) summary->setPassAnotherPolNchnl(true);
         }

      }
   }//end if if nchnlFilter=1


   /* CW filter */

   //if(settings->cwFilter > 0){

   bool isCW = false;
/*
   for(int ch=0; ch<16; ch++){

      grFFT[ch] = FFTtools::makePowerSpectrumMilliVoltsNanoSecondsdB(grWinPad[ch]);

      if(ch==0){
         freqCountLen_V = grFFT[ch]->GetN();
         freqCount_V = new int [freqCountLen_V];
         fill(&freqCount_V[0], &freqCount_V[freqCountLen_V], 0);
         freqBinWidth_V = evProcessTools::getFFTBinWidth(grFFT[ch]);
      }
      else if (ch==8){
         freqCountLen_H = grFFT[ch]->GetN();
         freqCount_H = new int [freqCountLen_H];
         fill(&freqCount_H[0], &freqCount_H[freqCountLen_H], 0);
         freqBinWidth_H = evProcessTools::getFFTBinWidth(grFFT[ch]);
      }

      summary->setFreqBinWidth(freqBinWidth_V, freqBinWidth_H);

      int maxFracBin;
      double maxFrac = FFTtools::getPeakVal(grFFT[ch], &maxFracBin);
      //cout<<"ch: "<<ch<<" macFrac: "<<maxFrac<<" maxFracBin: "<<maxFracBin<<endl;
      summary->setMaxFreqBinByChannel(ch, maxFracBin, maxFrac);
      if(ch<8)
      freqCount_V[maxFracBin]++;
      else
      freqCount_H[maxFracBin]++;

   }//end of ch
*/
   maxCount_V = evProcessTools::getMaxCount(freqCountLen_V, freqCount_V, &maxCountBin, 2); // a coincidence of >= 2 counts as potentially meaning
   maxCountFreq_V = freqBinWidth_V * maxCountBin;

   maxCount_H = evProcessTools::getMaxCount(freqCountLen_H, freqCount_H, &maxCountBin, 2);
   maxCountFreq_H = freqBinWidth_H * maxCountBin;

   summary->setMaxCountFreq(maxCountFreq_V, maxCountFreq_H);

   /* If coincidence >= minCWCoincidence */
   if(maxCount_V >= minCWCoincidence || maxCount_H >= minCWCoincidence) isCW = true;

   /* If coincidence = minCWCoincidence - 1 and there is at least one neighboring peak bin */

   //Check Vpol
   if (!isCW){
      if( ((freqCount_V[0] == minCWCoincidence-1) && (freqCount_V[1] > 0) ) ||
      ((freqCount_V[freqCountLen_V-1] == minCWCoincidence-1) && (freqCount_V[freqCountLen_V-2] > 0))
      ) isCW = true;
   }
   if (!isCW){
      for(int i=1; i<freqCountLen_V-1; i++){
         if(freqCount_V[i]==minCWCoincidence-1 && (freqCount_V[i-1]>0 || freqCount_V[i+1]>0)){
            isCW = true;
            break;
         }
      }
   }

   //Check Hpol
   if (!isCW){
      if( ((freqCount_H[0] == minCWCoincidence-1) && (freqCount_H[1] > 0) ) ||
      ((freqCount_H[freqCountLen_H-1] == minCWCoincidence-1) && (freqCount_H[freqCountLen_H-2] > 0))
      ) isCW = true;
   }
   if (!isCW){
      for(int i=1; i<freqCountLen_H-1; i++){
         if(freqCount_H[i]==minCWCoincidence-1 && (freqCount_H[i-1]>0 || freqCount_H[i+1]>0)){
            isCW = true;
            break;
         }
      }
   }

   if(settings->cwFilter > 0){

      if(isCW){

         cwFilteredEventCount+=1;
         weightedCWFilteredEventCount += weight;
         unpaddedEvent.clear();
         cleanEvent.clear();
         for(int ch=0; ch<16; ch++){ delete grInt[ch]; delete grWinPad[ch]; delete grMean[ch]; /*delete grFFT[ch];*/ /*delete grCDF[ch]; */}
         continue;

      }
   }

   /* Constant N reco for surface filter */

//   recoSuccess = false;
//
//   if(settings->constantNFilter > 0){
//      while( !recoSuccess ){
//
//         stringstream ss;
//         ss << ev/*rawAtriEvPtr->eventNumber*/;
//         evStr = ss.str();
//         fitsFileStr = fitsFile_tmp + /*".ev" + evStr +*/ ".constantN.fits";
//         sprintf(fitsFile, fitsFileStr.c_str());
//         constantNMaxPixIdx = reconstruct3DXCorrEnvelopeGetMaxPixAndMapData_constantNFilter(settings, cleanEvent, &clEnv, constantNDelays, constantNDelays_V, constantNDelays_H, goodChan, summary, fitsFile///*argv[5]*/, mapData/*, xCorrAroundPeakHist, sillygr*/
//         );
//
//         if( constantNMaxPixIdx < 0){ cerr<<"Error reconstructing - contant N\n"; return -1; }
//         if(summary->constantNMaxPixCoherence != 0.f) recoSuccess = true; //To catch cases where GPU reco returns coherence value zero
//         else { cout<<"constantNMaxPixCoherence returns 0!! Re-running reco...\n"; }
//
//      }
//
//      if(recordConstantNDir(settings, summary) < 0){ cerr<<"Error recording constant N dir\n"; return -1;}
//
//   }
//
//   /* Track engine object to compute all tracks */
//   treg->computeAllTracks(unpaddedEvent);
//   summary->setTreg(treg);

   //****************************************************
   // END OF FILTER SECTION
   //****************************************************


   for(int a=0;a<16;a++)
   {
      string_i  = detector->getStringfromArbAntID(0, a);
      antenna_i = detector->getAntennafromArbAntID(0, a);
      //AraRootChannel = detector->GetChannelfromStringAntenna (0, string_i, antenna_i, AraSim_settings);

      trueRecAngles[a] = report->stations[0].strings[string_i].antennas[antenna_i].rec_ang[0];
      trueLauAngles[a] = report->stations[0].strings[string_i].antennas[antenna_i].launch_ang[0];

   }

   summary->setTrueAngles(trueRecAngles, trueLauAngles);

//   getChannelSNR(unpaddedEvent, snrArray);
//   summary->setChannelInWindowSNR(snrArray);
//
//   getChannelTotalPowerSNR(unpaddedEvent, (int)settings->powerEnvIntDuration/settings->wInt_V, (int)settings->powerEnvIntDuration/settings->wInt_H, snrArray);
//   summary->setChannelTotalPowerSNR(snrArray);
//
//   getChannelSlidingV2SNR(unpaddedEvent, (int)settings->powerEnvIntDuration/settings->wInt_V, (int)settings->powerEnvIntDuration/settings->wInt_H, snrArray);
//   summary->setChannelSlidingV2SNR(snrArray);
//
//   /* Use sliding V^2 SNR as the event SNR */
//
//   TMath::Sort(16,snrArray,index);
//
//   for(int i=0; i<8; i++){
//      snrArray_V[i] = snrArray[i];
//   }
//   TMath::Sort(8,snrArray_V,index_V);
//
//   for(int i=0; i<8; i++){
//      snrArray_H[i] = snrArray[i+8];
//   }
//   TMath::Sort(8,snrArray_H,index_H);
//
//   getChannelUnmodifiedSNR(unpaddedEvent, unmodSNRArray);
//   //TMath::Sort(16,unmodSNRArray,index);
//
//   snrRank = (string(settings->recoPolType)=="both" ? index : (string(settings->recoPolType)=="vpol" ? index_V : index_H));
//

   //recoData *summary = new recoData();
   //if(settings->dataType == 0){
   //summary->setWeight(event->Nu_Interaction[0].weight);
//   summary->setSurvivalProbability(event->Nu_Interaction[0].weight);
//   double L_ice = event->Nu_Interaction[0].len_int_kgm2_total/Signal::RHOICE;
//   double p_int = 1.-exp(-1.*(event->Nu_Interaction[0].r_enterice.Distance(event->Nu_Interaction[0].nuexitice)/L_ice)); // probability it interacts in ice along its path
//   summary->setProbabilities(p_int, p_int*event->Nu_Interaction[0].weight);
//   double weight = computeWeight(AraSim_settings, detector, event, icemodel, -1.*stationCenterDepth, L_ice);
//   if(weight<0){ cerr<<"Weight compuation error!\n"; return -1; }
//   summary->setWeight(weight);
   summary->setTrueRadius(r_true);
   summary->setTrueDir(zen_true*180./M_PI, azi_true*180./M_PI);
   //}
   //summary->setOnion(onion);
   summary->setTopN(topN);
   summary->setRecoChan(goodChan);
   //summary->setInWindowSNR(snrArray[index[2]]);
   //summary->setInWindowSNRBothPol(snrArray_V[index_V[2]], snrArray_H[index_H[2]]);
   //summary->setUnmodSNR(unmodSNRArray[index[2]]);
//   if(settings->nchnlFilter>0){
//      if(settings->nchnlFilter==1){
//         if(snrArray_H[index_H[2]]>=settings->nchnlThreshold_anotherPol) summary->setPassAnotherPolNchnl(true);
//      }
//      else if(settings->nchnlFilter==2){
//         if(snrArray_V[index_V[2]]>=settings->nchnlThreshold_anotherPol) summary->setPassAnotherPolNchnl(true);
//      }
//   }
//
   recoEventCount++;
   weightedRecoEventCount += weight;
   /* Reco with radiospline */

//   recoSuccess = false;
//
//   while( !recoSuccess ){
//   if(settings->beamformMethod == 1){
//   if(settings->getSkymapMode == 0){
//       err = reconstructCSW(settings, cleanEvent, &clEnv, recoDelays, recoDelays_V, recoDelays_H, nDir, chanMask, fitsFile/*argv[5]*/);
//   }
//   else{
//       maxPixIdx = reconstructCSWGetMaxPix(settings, cleanEvent, &clEnv, recoDelays, recoDelays_V, recoDelays_H, nDir, chanMask, summary);
//   }
//   } else {
//   if(settings->getSkymapMode == 0){
//
//      //evStr = std::to_string(ev);
//      stringstream ss;
//      ss << ev;
//      evStr = ss.str();
//      fitsFileStr = fitsFile_tmp + /*".ev" + evStr +*/ ".fits";
//      sprintf(fitsFile, fitsFileStr.c_str());
//
//      if(settings->skymapSearchMode == 0){ //no zoom mode
//      maxPixIdx = reconstruct3DXCorrEnvelopeGetMaxPixAndMapData(settings, cleanEvent, &clEnv, recoDelays, recoDelays_V, recoDelays_H, goodChan, snrRank, summary, fitsFile/*argv[5]*/, mapData/*, xCorrAroundPeakHist, sillygr*/);
//
//      if(settings->use2ndRayReco){
//      fitsFileStr = fitsFile_tmp + /*".ev" + evStr +*/ ".2ndRay.fits";
//      sprintf(fitsFile, fitsFileStr.c_str());
//      maxPixIdx2 = reconstruct3DXCorrEnvelopeGetMaxPixAndMapData_2ndRayReco(settings, cleanEvent, &clEnv, recoRefracDelays, recoRefracDelays_V, recoRefracDelays_H, goodChan, summary, fitsFile/*argv[5]*/, mapData/*, xCorrAroundPeakHist, sillygr*/);
//      }
//
//      if(settings->recordMapData == 1){
//      for(int pix=0; pix<nDir*nLayer; pix++) mapDataHist[pix]->Fill(mapData[pix]);
//      }
//      } else {//zoom search mode
//      maxPixIdx = reconstruct3DXCorrEnvelopeGetMaxPix_ZoomMode(settings, cleanEvent, &clEnv, stationCenterDepth, antLocation, recoDelays, recoDelays_V, recoDelays_H, goodChan, summary, fitsFile);
//      }
//   } else {
//      maxPixIdx = reconstructXCorrEnvelopeGetMaxPix(settings, cleanEvent, &clEnv, recoDelays, recoDelays_V, recoDelays_H, nDir, goodChan/*chanMask*/, summary);
//
//   }
//   }
//
//   if( err<0 || maxPixIdx<0 || maxPixIdx2<0 ){ cerr<<"Error reconstructing\n"; return -1; }
//   if(summary->maxPixCoherence != 0.f) recoSuccess = true; //To catch cases where GPU reco returns coherence value zero
//   if(settings->use2ndRayReco==1){ if(summary->maxPixCoherence2 == 0.f) recoSuccess=false;}
//   if(recoSuccess!=true) { cout<<"maxPixCoherence/2 returns 0!! Re-running reco...\n"; }
//   //else { cout<<"maxPixCoherence returns 0!! Re-running reco...\n"; }
//   }//end of while
//
//   //int recoFlag = record3DDiffGetFlag(summary, outputFile);
//   //if( recoFlag ) recoFlagCnt++;
//   if(settings->use2ndRayReco==0){
//   summary->setFlag( (settings->skymapSearchMode)
//                    ? record3DZoomedDiffGetFlag(settings, summary, dZenDist, dAziDist, recoTrueZenDist, recoTrueAziDist)
//                    : record3DDiffGetFlag(settings, summary, dZenDist, dAziDist, recoTrueZenDist, recoTrueAziDist) );
//   } else {
//   if(settings->skymapSearchMode){ cerr<<"skymapSearchMode==1 & use2ndRayReco==1 incompatible!\n"; return -1;}
//   else summary->setFlag(record3DDiffGetFlag_2ndRayReco(settings, summary, dZenDist, dAziDist, recoTrueZenDist, recoTrueAziDist));
//   }
//   if(summary->flag > 0) recoFlagCnt++;
//   maxPix[maxPixIdx]++;
//
//   if(settings->use2ndRayReco==0)
//   compute3DRecoAnglesWithRadioSplineForSinglePixel(nAnt, -1.f*stationCenterDepth, antLocation, onion, recoLauAngles, recoRecAngles, maxPixIdx);
//   else{
//      if(summary->maxPixCoherence >= summary->maxPixCoherence2)
//      compute3DRecoAnglesWithRadioSplineForSinglePixel(nAnt, -1.f*stationCenterDepth, antLocation, onion, recoLauAngles, recoRecAngles, maxPixIdx);
//      else {
//      cerr<<"2nd ray angles table not yet build! Angles will be default (-1)\n"; //return -1;
//      //compute3DRecoAnglesWithRadioSplineForSinglePixel_2ndRayReco(nAnt, -1.f*stationCenterDepth, antLocation, onion, recoLauAngles, recoRecAngles, maxPixIdx2);
//      }
//   }
//   //compute3DRecoAnglesWithRadioSplineForSinglePixel(nAnt, -1.f*stationCenterDepth, antLocation, onion, recoLauAngles, recoRecAngles, maxPixIdx);
//
//   summary->setRecoAngles(recoRecAngles, recoLauAngles);

   dataTree->Fill();
   unpaddedEvent.clear();
   cleanEvent.clear();
   //delete summary;
   //treg->clearForNextEvent();
   for(int ch=0; ch<16; ch++){ /*delete gr_v[ch];*/ delete grInt[ch]; delete grWinPad[ch]; delete grMean[ch]; /*delete grCDF[ch];*/ /*if(settings->cwFilter>0)*/ /*delete grFFT[ch];*/ }
   }//end of ev loop

}//end of dataType == 0

//recordTime(tmr,4);
time_t t_after_event_loop = time(NULL);
clock_t c_after_event_loop = clock();

cout<<"runEventCount: "<<runEventCount<<" recoEventCount: "<<recoEventCount<<" trigEventCount: "<<trigEventCount<<endl;
cout<<"cutWaveEventCount: "<<cutWaveEventCount<<" nonIncreasingSampleTimeEventCount: "<<nonIncreasingSampleTimeEventCount<<" cutWaveAndNonIncreasingEventCount: "<<cutWaveAndNonIncreasingEventCount<<" corruptD1EventCount: "<<corruptD1EventCount<<" corruptFirst3EventCount: "<<corruptFirst3EventCount<<endl;
cout<<"mistaggedSoftEventCount: "<<mistaggedSoftEventCount<<" offsetBlockEventCount: "<<offsetBlockEventCount<<" nchnlFilteredEventCount: "<<nchnlFilteredEventCount<<" cwFilteredEventCount: "<<cwFilteredEventCount<<" impulsivityFilteredEventCount: "<<impulsivityFilteredEventCount<<endl;

runInfoTree->Fill();

//outputFile->Write();
//outputFile->Close();

clfftTeardown();
err = tearDown(&clEnv);
delete treg;
delete summary;
//if(settings->skymapSearchMode == 0){
//   delete onion;
//   free(recoDelays);
//   free(recoDelays_V);
//   free(recoDelays_H);
//   if(settings->use2ndRayReco==1){
//   free(recoRefracDelays);
//   free(recoRefracDelays_V);
//   free(recoRefracDelays_H);
//   }
//}
//if(settings->constantNFilter > 0){
//   delete onion_temp;
//   free(constantNDelays);
//   free(constantNDelays_V);
//   free(constantNDelays_H);
//}
//delete settings;
free(mapDataHist);
free(mapData);

if( !settings->readRecoSetupFile( recoSetupFile_fullPath )){

   cerr<<"Error reading the recoSetupFile or invalid parameters !! Aborting now...\n";
   return -1;
   //cerr<<"Will use default reco setup file\n";
   //settings->readRecoSetupFile("recoSetupFile_default.txt");
   //outputFile = new TFile(("recoOut.recoSetupFile_default.txt.run"+runNum+".root").c_str(),"RECREATE","recoOut");
   //fitsFile_tmp = "recoSkymap.recoSetupFile_default.txt.run" + runNum/* + ".fits"*/;

} else {
   cout<<"Obtained new reoSetupFile\n";
   outputFile = new TFile((outputDir+"fftBaseline."+recoSetupFile+".run"+runNum+".root").c_str(),"RECREATE","fftBaseline");
   fitsFile_tmp = outputDir + "recoSkymap." + recoSetupFile + ".run" + runNum/* + ".fits"*/;
}

//double median[16] = {-1000.};
//double percentile_75[16] = {-1000.};
//double percentile_25[16] = {-1000.};
vector<double> tempVec;

TGraph *gr_median[16];
TGraph *gr_75[16];
TGraph *gr_25[16];
char grName[200];

for(int ch=0; ch<8; ch++){
   sprintf(grName,"gr_median_%d",ch);
   gr_median[ch] = new TGraph();
   gr_median[ch]->SetName(grName); gr_median[ch]->SetTitle(grName);
   sprintf(grName,"gr_75_%d",ch);
   gr_75[ch]     = new TGraph();
   gr_75[ch]->SetName(grName); gr_75[ch]->SetTitle(grName);
   sprintf(grName,"gr_25_%d",ch);
   gr_25[ch]     = new TGraph();
   gr_25[ch]->SetName(grName); gr_25[ch]->SetTitle(grName);
   for(int bin=0; bin<freqCountLen_V; bin++){
      tempVec.assign(fftValues_V[ch*freqCountLen_V+bin].begin(), fftValues_V[ch*freqCountLen_V+bin].end());
      //median[ch] = getPercentile(tempVec, 0.5);
      //percentile_75[ch] = getPercentile(tempVec, 0.75);
      //percentile_25[ch] = getPercentile(tempVec, 0.25);
      gr_median[ch]->SetPoint(bin, bin*freqBinWidth_V, getPercentile(tempVec, 0.5));
      gr_75[ch]->SetPoint(bin, bin*freqBinWidth_V, getPercentile(tempVec, 0.75));
      gr_25[ch]->SetPoint(bin, bin*freqBinWidth_V, getPercentile(tempVec, 0.25));
      //cout<<"ch: "<<ch<<" 25%: "<<percentile_25[ch]<<" 50%: "<<median[ch]<<" 75%: "<<percentile_75[ch]<<endl;
      tempVec.clear();
   }
   gr_median[ch]->Write();
   gr_75[ch]->Write();
   gr_25[ch]->Write();
}
for(int ch=8; ch<16; ch++){
   sprintf(grName,"gr_median_%d",ch);
   gr_median[ch] = new TGraph();
   gr_median[ch]->SetName(grName); gr_median[ch]->SetTitle(grName);
   sprintf(grName,"gr_75_%d",ch);
   gr_75[ch]     = new TGraph();
   gr_75[ch]->SetName(grName); gr_75[ch]->SetTitle(grName);
   sprintf(grName,"gr_25_%d",ch);
   gr_25[ch]     = new TGraph();
   gr_25[ch]->SetName(grName); gr_25[ch]->SetTitle(grName);
   for(int bin=0; bin<freqCountLen_H; bin++){
      tempVec.assign(fftValues_H[(ch-8)*freqCountLen_H+bin].begin(), fftValues_H[(ch-8)*freqCountLen_H+bin].end());
      //median[ch] = getPercentile(tempVec, 0.5);
      //percentile_75[ch] = getPercentile(tempVec, 0.75);
      //percentile_25[ch] = getPercentile(tempVec, 0.25);
      gr_median[ch]->SetPoint(bin, bin*freqBinWidth_H, getPercentile(tempVec, 0.5));
      gr_75[ch]->SetPoint(bin, bin*freqBinWidth_H, getPercentile(tempVec, 0.75));
      gr_25[ch]->SetPoint(bin, bin*freqBinWidth_H, getPercentile(tempVec, 0.25));
      //cout<<"ch: "<<ch<<" 25%: "<<percentile_25[ch]<<" 50%: "<<median[ch]<<" 75%: "<<percentile_75[ch]<<endl;
      tempVec.clear();
   }
   gr_median[ch]->Write();
   gr_75[ch]->Write();
   gr_25[ch]->Write();
}
//outputFile->Write();
outputFile->Close();

delete settings;
/*
TH2F *fftValuesHist =new TH2F("fftValuesHist","fftValuesHist",(int)(1000/freqBinWidth_V),0,1000,500,-100,100);
for(int bin=0; bin<freqCountLen_V; bin++){
   for(int val=0; val<fftValues_V[bin].size(); val++){
      fftValuesHist->Fill(bin*freqBinWidth_V,fftValues_V[bin].at(val));
   }
}
*/
/*
TCanvas c1("c1","c1",800,800);
//fftValuesHist->Draw("colz");
c1.Divide(4,4);
for(int ch=0; ch<16; ch++){
   c1.cd(ch+1);
   gr_median[ch]->Draw("AL");
   gr_75[ch]->SetLineColor(kRed);
   gr_75[ch]->Draw("Lsame");
   gr_25[ch]->SetLineColor(kBlue);
   gr_25[ch]->Draw("Lsame");
}
c1.SaveAs("fftValuesMedian.C");
//c1.SaveAs("fftValuesHist.C");
*/
//recordTime(tmr,5);
time_t t_program_end = time(NULL);
clock_t c_program_end = clock();

/*
printf("Chrono time (us)\n");
printf("Total program : %f\tTotal overhead : %f\tGenerating recoDelay: %f\tReco loop: %f\t program - (reco loop + overhead): %f\n",
chrono::duration_cast<chrono::microseconds>(tmr->chronoVec[5] - tmr->chronoVec[0]).count(),
chrono::duration_cast<chrono::microseconds>(tmr->chronoVec[3] - tmr->chronoVec[0]).count() +
chrono::duration_cast<chrono::microseconds>(tmr->chronoVec[5] - tmr->chronoVec[4]).count(),
chrono::duration_cast<chrono::microseconds>(tmr->chronoVec[2] - tmr->chronoVec[1]).count(),
chrono::duration_cast<chrono::microseconds>(tmr->chronoVec[4] - tmr->chronoVec[3]).count(),
chrono::duration_cast<chrono::microseconds>(tmr->chronoVec[5] - tmr->chronoVec[0]).count() -
chrono::duration_cast<chrono::microseconds>(tmr->chronoVec[4] - tmr->chronoVec[3]).count() -
(chrono::duration_cast<chrono::microseconds>(tmr->chronoVec[3] - tmr->chronoVec[0]).count() +
chrono::duration_cast<chrono::microseconds>(tmr->chronoVec[5] - tmr->chronoVec[4]).count())
);

printf("Chrono time (ns)\n");
printf("Total program : %f\tTotal overhead : %f\tGenerating recoDelay: %f\tReco loop: %f\t program - (reco loop + overhead): %f\n",
chrono::duration_cast<chrono::nanoseconds>(tmr->chronoVec[5] - tmr->chronoVec[0]).count(),
chrono::duration_cast<chrono::nanoseconds>(tmr->chronoVec[3] - tmr->chronoVec[0]).count() +
chrono::duration_cast<chrono::nanoseconds>(tmr->chronoVec[5] - tmr->chronoVec[4]).count(),
chrono::duration_cast<chrono::nanoseconds>(tmr->chronoVec[2] - tmr->chronoVec[1]).count(),
chrono::duration_cast<chrono::nanoseconds>(tmr->chronoVec[4] - tmr->chronoVec[3]).count(),
chrono::duration_cast<chrono::nanoseconds>(tmr->chronoVec[5] - tmr->chronoVec[0]).count() -
chrono::duration_cast<chrono::nanoseconds>(tmr->chronoVec[4] - tmr->chronoVec[3]).count() -
(chrono::duration_cast<chrono::nanoseconds>(tmr->chronoVec[3] - tmr->chronoVec[0]).count() +
chrono::duration_cast<chrono::nanoseconds>(tmr->chronoVec[5] - tmr->chronoVec[4]).count())
);
*/
printf("Clock time (s)\n");
printf("Total program : %f\tTotal overhead : %f\tGenerating recoDelay: %f\tReco loop: %f\t program - (reco loop + overhead): %f\n",
/*(double)(tmr->clockVec[5]-tmr->clockVec[0])/CLOCKS_PER_SEC,
(double)(tmr->clockVec[3]-tmr->clockVec[0])/CLOCKS_PER_SEC +
(double)(tmr->clockVec[5]-tmr->clockVec[4])/CLOCKS_PER_SEC,
(double)(tmr->clockVec[2]-tmr->clockVec[1])/CLOCKS_PER_SEC,
(double)(tmr->clockVec[4]-tmr->clockVec[3])/CLOCKS_PER_SEC,
(double)(tmr->clockVec[5]-tmr->clockVec[0])/CLOCKS_PER_SEC -
(double)(tmr->clockVec[4]-tmr->clockVec[3])/CLOCKS_PER_SEC -
((double)(tmr->clockVec[3]-tmr->clockVec[0])/CLOCKS_PER_SEC +
(double)(tmr->clockVec[5]-tmr->clockVec[4])/CLOCKS_PER_SEC)
*/
(double)(c_program_end - c_program_start)/CLOCKS_PER_SEC,
(double)(c_before_event_loop - c_program_start + c_program_end - c_after_event_loop)/CLOCKS_PER_SEC,
(double)(c_after_recoDelays - c_before_recoDelays)/CLOCKS_PER_SEC,
(double)(c_after_event_loop - c_before_event_loop)/CLOCKS_PER_SEC,
(double)(c_program_end - c_program_start)/CLOCKS_PER_SEC -
(double)(c_after_event_loop - c_before_event_loop)/CLOCKS_PER_SEC -
(double)(c_before_event_loop - c_program_start + c_program_end - c_after_event_loop)/CLOCKS_PER_SEC
);

printf("Time_t time (s)\n");
printf("Total program : %f\tTotal overhead : %f\tGenerating recoDelay: %f\tReco loop: %f\t program - (reco loop + overhead): %f\n",
/*difftime(tmr->timeVec[5], tmr->timeVec[0]),
difftime(tmr->timeVec[3], tmr->timeVec[0]) +
difftime(tmr->timeVec[5], tmr->timeVec[4]),
difftime(tmr->timeVec[2], tmr->timeVec[1]),
difftime(tmr->timeVec[4], tmr->timeVec[3]),
difftime(tmr->timeVec[5], tmr->timeVec[0]) -
difftime(tmr->timeVec[4], tmr->timeVec[3]) -
(difftime(tmr->timeVec[3], tmr->timeVec[0]) +
difftime(tmr->timeVec[5], tmr->timeVec[4]))
*/
difftime(t_program_end, t_program_start),
difftime(t_before_event_loop, t_program_start) + difftime(t_program_end, t_after_event_loop),
difftime(t_after_recoDelays, t_before_recoDelays),
difftime(t_after_event_loop, t_before_event_loop),
difftime(t_program_end, t_program_start) -
difftime(t_after_event_loop, t_before_event_loop) -
difftime(t_before_event_loop, t_program_start) + difftime(t_program_end, t_after_event_loop)
);

cout<<"Successfully reached end of main()"<<endl;
return 0;
}

double getMean(TGraph *gr){

double v, t, mean;
mean=0;

for(int i=0; i<gr->GetN(); i++){

gr->GetPoint(i,t,v);
mean+=v;
}

mean = mean / double(gr->GetN());
return mean;
}

double computeWeight(Settings *settings, Detector *detector, Event *event, IceModel *antarctica, double zCenter, double L_int){

   /* The correct weight is W = Wprop * Wint = Interaction::weight * L_0 / L_int.
      The only quatity we need to compute here is L_0.
   */

   if(zCenter > 0){ cerr<<"zCenter should be <0 ! zCenter="<<zCenter<<endl; return -1; }
   //double zCenter = -180.;
   double ox, oy, oz; // coor. of the origin
   ox = detector->stations[0].GetX();
   oy = detector->stations[0].GetY();
   oz = detector->stations[0].GetZ() + zCenter;

   double posx, posy, posz; //vertex position
   double nnux, nnuy, nnuz; //neutrino direction
   //double t;                //r = r_0 + t * nnu
   double t[4], z[4]; //0: surface, 1: bedrock, 2: side solution 1, 3: side solution 2
   int index[4] = {0};

   posx = event->Nu_Interaction[0].posnu.GetX() - ox;
   posy = event->Nu_Interaction[0].posnu.GetY() - oy;
   posz = event->Nu_Interaction[0].posnu.GetZ() - oz;
   nnux = event->Nu_Interaction[0].nnu.GetX();
   nnuy = event->Nu_Interaction[0].nnu.GetY();
   nnuz = event->Nu_Interaction[0].nnu.GetZ();

   double bedrock_z = antarctica->Surface(event->Nu_Interaction[0].posnu) - antarctica->IceThickness(event->Nu_Interaction[0].posnu) - oz;

   double PR = settings->POSNU_RADIUS;

   t[0] = (-zCenter - posz) / nnuz;
   z[0] = -zCenter;
   t[1] = (bedrock_z - posz)/ nnuz;
   z[1] = bedrock_z;

   double a = nnux*nnux + nnuy*nnuy;
   double b = 2*(posx*nnux + posy*nnuy);
   double c = posx*posx + posy*posy - PR*PR;
   t[2] = (-b + sqrt(b*b-4.*a*c)) / (2.*a);
   z[2] = posz + nnuz * t[2];
   t[3] = (-b - sqrt(b*b-4.*a*c)) / (2.*a);
   z[3] = posz + nnuz * t[3];

   //for(int i=0; i<4; i++) t[i] = fabs(t[i]); //only care about the absolute distance
   TMath::Sort(4, z, index); //in descending order. index[1] and index[2] represents generation volume entry and exit, order unimportant
   double dt = fabs(t[index[2]] - t[index[1]]);
   if(dt<0){ cerr<<"Error! dt < 0! \n"; return -1; }

   double L0x, L0y, L0z;
   L0x = nnux * dt;
   L0y = nnuy * dt;
   L0z = nnuz * dt;
   double L_0 = sqrt(L0x*L0x + L0y*L0y + L0z*L0z);

   double weight = event->Nu_Interaction[0].weight * L_0 / L_int;

   printf("zCenter: %f ox: %f oy: %f oz: %f\nposx: %f posy: %f posz: %f\nnnux: %f nnuy: %f nnuz: %f\nbedrock_z: %f PR: %f\nz0: %f z1: %f z2: %f z3: %f\nt0: %f t1: %f t2: %f t3: %f dt: %f\nL0x: %f L0y: %f L0z: %f L_0: %f\nweight: %f\n",
   zCenter, ox, oy, oz,
   posx, posy, posz,
   nnux, nnuy, nnuz,
   bedrock_z, PR,
   z[0], z[1], z[2], z[3],
   t[0], t[1], t[2], t[3], dt,
   L0x, L0y, L0z, L_0,
   weight);

   return weight;
}
