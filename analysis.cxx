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

#include "Detector.h"
#include "Trigger.h"
#include "Settings.h"
#include "Report.h"
#include "Event.h"

ClassImp(recoSettings);
ClassImp(recoData);

//class recoSettings;
//class recoData;

//#define XCorr

//#ifndef XCorr
//#define CSW
//#endif

//#define GetMaxPixMode
//#define GetMaxPixAndMapMode

//#define constantN
//#define nchnl

//#define TOP_N 50 //top N max pixels in whole Healpix_Onion
//#define N_LAYER 1 //Number of Healpix_Onion layers

using namespace std;

RawIcrrStationEvent *rawIcrrEvPtr;
RawAtriStationEvent *rawAtriEvPtr;
RawAraStationEvent *rawEvPtr;
//UsefulIcrrStationEvent *realIcrrEvPtr;
//UsefulAtriStationEvent *realAtriEvPtr;

double getMean(TGraph *);
/*
#ifndef XCORRSUMGRAPH
static TGraph *sillygr = new TGraph();
//TGraph *envelopeSum = new TGraph();
#define XCORRSUMGRAPH
#endif
*/

int main( int argc, char **argv){

if(argc < 3){ cerr<<"Insufficient arguments. Usage: 1.recoSetupFile 2. Run Number 3. Data ROOT file 4. AraSim: more Data ROOT file. Real: Pedestal File\n"; return -1; }

int err;
gROOT->ProcessLine("#include <vector>");
/* 
 * Specify the channels to be used in the analysis
 * 1: use, 0: don't use
 */

int chanMask[16] = {  1 //chan 0  D1TV 
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

TFile *outputFile;
string recoSetupFile = string( argv[1] );
string runNum = string( argv[2] );
string fitsFile_tmp;
string fitsFileStr;
string evStr;
cout<<"recoStupFile: "<<recoSetupFile<<endl;
if( !settings->readRecoSetupFile( recoSetupFile )){ 
   cerr<<"Will use default reco setup file\n"; 
   settings->readRecoSetupFile("recoSetupFile_default.txt");
   outputFile = new TFile(("recoOut.recoSetupFile_default.txt.run"+runNum+".root").c_str(),"RECREATE","recoOut");
   fitsFile_tmp = "recoSkymap.recoSetupFile_default.txt.run" + runNum/* + ".fits"*/;
} else {
   cout<<"Obtained new reoSetupFile\n";
   outputFile = new TFile(("recoOut."+recoSetupFile+".run"+runNum+".root").c_str(),"RECREATE","recoOut");
   fitsFile_tmp = "recoSkymap." + recoSetupFile + ".run" + runNum/* + ".fits"*/; 
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

TTree *recoSettingsTree = new TTree("recoSettingsTree", "recoSettingsTree");
//TTree *onionTree = new TTree("onionTree", "onionTree");
TTree *dataTree  = new TTree("dataTree",  "dataTree");

recoSettingsTree->Branch("settings", &settings);
recoSettingsTree->Fill();

TH1F *dZenDist = new TH1F("recoZenDiff", "recoZenDiff", 360, -180, 180);
TH1F *dAziDist = new TH1F("recoAziDiff", "recoAziDiff", 720, -360, 360); 
TH2F *recoTrueZenDist = new TH2F("recoTrueZenDist", "recoTrueZenDist", 180, 0, 180, 180, 0, 180);
TH2F *recoTrueAziDist = new TH2F("recoTrueAziDist", "recoTrueAziDist", 360, 0, 360, 360, 0, 360);

float radius, r_xy;
float dx, dy, dz;
float r_true, zen_true, azi_true;

/* 
 * Variables used in dataType == 0 case
 */
 
TChain chain("AraTree"), chain2("AraTree2");
Settings *AraSim_settings = 0;
Detector *detector = 0;
Event    *event    = 0;
Report   *report   = 0;
Trigger  *trigger  = 0;

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
int cutWaveAlert, nonIncreasingSampleTimeAlert;
double previous_times;
double addDelay;
double times, volts;
double time_1, time_2, time_last;
int utime, utime_runStart, utime_runEnd;

/*
 * Variables used in nchnlFilter > 0 case
 */

double threshold;
int nchnlArray[3];
int nchnl_tmp;

/* End of conditional variables pre-declaration */
 
Long64_t numEntries;

if(settings->dataType == 1)//real events
{

   fp = TFile::Open( argv[3] );
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
   numEntries=eventTree->GetEntries();
   cerr<<"isAtri "<<isAtriEvent<<" isIcrr "<<isIcrrEvent<<" number of entries is "<<numEntries<<endl;

   /*
    * START LOADING GOOD PED
    */
   sprintf(dir_char,argv[4]);
   printf("ped file is %s\n",dir_char);
   calib = AraEventCalibrator::Instance();
   calib->setAtriPedFile(dir_char,rawEvPtr->stationId);
   //************END OF LOADING GOOD PED************
}
else if (settings->dataType == 0)//AraSim events
{
   
   for(int i=3; i<argc; i++){
      chain.Add( argv[i] );
      chain2.Add( argv[i] );
   }

   chain.SetBranchAddress("settings",&AraSim_settings);
   chain.SetBranchAddress("detector",&detector);
   chain.SetBranchAddress("trigger" ,&trigger);
   chain2.SetBranchAddress("report" ,&report);
   chain2.SetBranchAddress("event"  ,&event);

   chain.GetEntry(0);
   cout<<"EXPONENT: "<<AraSim_settings->EXPONENT<<endl;
   cout<<"NNU: "<<AraSim_settings->NNU<<endl;
   printf("Station center X: %f Y: %f Z: %f\n",detector->stations[0].GetX(),detector->stations[0].GetY(),detector->stations[0].GetZ());

   numEntries = chain2.GetEntries();
  
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
Healpix_Onion *onion;

if( settings->skymapSearchMode == 0){ //No zoom search

/* Set n-side for Healpix and compute reco delays */

   nSideExp = settings->nSideExp;
   nLayer = settings->nLayer;
   nDir = 12 * pow(2, nSideExp) * pow(2, nSideExp);
   onion = new Healpix_Onion(nSideExp, nLayer);


   if(nDir*nLayer < topN) {
      cerr<<"topN greater than total number of pixles. Setting topN to equal nDir*nLayer...\n";
      topN = nDir*nLayer; //in case topN < total number of pixels
   }

   recoDelays  = (float*)malloc(nLayer*nDir*nAnt*sizeof(float));
   recoDelays_V= (float*)malloc(nLayer*nDir*(nAnt/2)*sizeof(float));
   recoDelays_H= (float*)malloc(nLayer*nDir*(nAnt/2)*sizeof(float));

   if(settings->iceModel == 1){
   err = computeRecoDelaysWithConstantN(nAnt, -1.f*stationCenterDepth, antLocation,
                                        //radius, nSideExp,
                                        onion, recoDelays, recoDelays_V, recoDelays_H);
   } else if(settings->iceModel == 0){
   err = compute3DRecoDelaysWithRadioSpline(nAnt, -1.f*stationCenterDepth, antLocation,
                                            onion, recoDelays, recoDelays_V, recoDelays_H);
   } else { cerr<<"Undefined iceModel parameter\n"; return -1; }
   if( err<0 ){ cerr<<"Error computing reco delays\n"; return -1; }
	
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
	
/*
 * Set up reco environment. In this case, an OpenCL environment
 */
recoEnvData clEnv;
err = setupCLRecoEnv(&clEnv, settings->programFile.c_str());
if( err<0 ){
   cerr<<"Error setting up reco env\n"; return -1;
}
	

   TGraph *gr_v_temp[16];
   TGraph *gr_v[16];
   TGraph *grInt[16];
   TGraph *grWinPad[16];
   TGraph *grNormWinPad[16];
   TGraph *gr_fft[16];
   TGraph *grHilbert[16];

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

/* 
 * Loop over events once to determine run start/end time 
 */

   for(int ev=1; ev<numEntries; ev++){
      eventTree->GetEntry(ev);
      if(rawAtriEvPtr->unixTime < utime_runStart) utime_runStart=rawAtriEvPtr->unixTime;
      if(rawAtriEvPtr->unixTime > utime_runEnd  ) utime_runEnd  =rawAtriEvPtr->unixTime;
   }
   cout<<"utime_runStart: "<<utime_runStart<<" dropD4Time: "<<dropD4Time<<endl;
   cout<<"Run time span: "<<utime_runEnd-utime_runStart<<endl;

}//end of if dataType = 1
/*
 * Start looping events for analysis
 */
vector<TGraph *> cleanEvent;
int recoEventCnt = 0;
int recoFlagCnt = 0;
double t, v, beginTime = 0.;

int *maxPix; 
if(settings->skymapSearchMode == 0) maxPix = (int*)calloc(nDir*nLayer, sizeof(int));
else                                maxPix = (int*)calloc(12*pow(2,settings->nSideExpEnd)*pow(2,settings->nSideExpEnd)*settings->nLayer, sizeof(int));

int maxPixIdx = 0;   
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
float snrArray[16], unmodSNRArray[16];  
vector<TGraph *> unpaddedEvent;
TH1F *snrDist = new TH1F("snrDist","snrDist",100,0,50);
int goodChan[16];
int numSatChan;
if(settings->nchnlFilter > 0){
   threshold = settings->nchnlThreshold;
   //int nchnlArray[3];
   //int nchnl_tmp;
}
int trigEvCnt = 0;
int triggerCode[3];
for(int i=0; i<3; i++){ triggerCode[i] = settings->triggerCode[i] - '0'; cout<<"triggerCode "<<i<<": "<<triggerCode[i]<<endl; }

cout<<"numEntries: "<<numEntries<<endl;

recoData *summary = new recoData();
dataTree->Branch("summary", &summary);

if(settings->dataType == 1){

trigEvCnt = numEntries;
for (Long64_t ev=0; ev<numEntries; ev++){

   summary->clear();
   
   if(ev%100 == 0) cout<<"*******************************Event got********************************: "<<ev<<endl;

   eventTree->GetEntry(ev);

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
      previous_times = 0.;
      double stdDelay = 0.;
   

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

      stdDelay= geom->getStationInfo(stationId)->getCableDelay(a);
      addDelay += stdDelay;

	  //*** We put the waveform into a graph. ***//
	  gr_v_temp[a] = realAtriEvPtr->getGraphFromRFChan(a);
	  gr_v[a] = new TGraph();
	  //*** A few waveforms are inverted in station 3, which needs to be corrected. I wrote a small method to do that (see above). ***//
	  if(stationId==3 && (a==0||a==4||a==8)){/*cout<<"graph inverted"<<endl;*/invertGraph(gr_v_temp[a]);}
	  //*** The following is to avoid reading corrupted waveforms. ***//
	  //*** I encountered only a few of them, so maybe this is not ***//
	  //*** really neccessary anymore. *******************************//
	  if(gr_v_temp[a]->GetN()<5 ){ cerr<< "BAD EVENT: " << ev << " Channel: " << a << ", points: " << gr_v_temp[a]->GetN() << endl;cutWaveAlert=1;continue;}
	  int pc = 0;

	  //*** The first 20 samples can be corrupted. Therefore, we need to exclude them! ***//
      for(int p=0;p<gr_v_temp[a]->GetN();p++){
	  gr_v_temp[a]->GetPoint(p, times, volts);
                       
         if(times>20.0) 			
         {
         if(stationId==3 && utime_runStart>=dropD4Time && (a%4==3)) 
         gr_v[a]->SetPoint(pc, times-addDelay, 0.); //Drop 2014 ARA03 D4
         else if(stationId==2 && a==15) gr_v[a]->SetPoint(pc, times-addDelay, 0.);//Drop ARA02 D4BH
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
            else {cerr<< "BAD EVENT Non-increasing sample time: " << event << " Channel: " << a << "this sample time: "<< times << "previous sample time: " << previous_times << endl;nonIncreasingSampleTimeAlert=1;}  
    
            previous_times = times;    

         }//end of pc

         } else {cerr<< "BAD EVENT type 2: " << event << " Channel: " << a << ", original number of points: " << gr_v_temp[a]->GetN() << endl;cutWaveAlert=1;/*continue;*/}

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
   if (cutWaveAlert == 1) { cerr<<"Event "<<ev<<" discarded due to cutWaveAlert\n"; continue; }
   if (nonIncreasingSampleTimeAlert == 1) { cerr<<"Event "<<ev<<" discarded due to nonIncreasingSampleTimeAlert\n"; continue; }

   for(int ch=0; ch<16; ch++){
     
      gr_v[ch]->GetPoint(0,t,v);
      if( t<beginTime ) beginTime = t ;

   }

   for(int ch=0; ch<16; ch++){


   if(ch<8){wInt=0.4; maxSamp=2048;}
   else{wInt=0.625; maxSamp=2048;}

   /* Interpolate + apply windowing + zero-pad + equalize wf beginning  to maxSamp */
   cout<<"N: "<<gr_v[ch]->GetN()<<endl; 
   grInt[ch]       = FFTtools::getInterpolatedGraph(gr_v[ch], wInt);
   unpaddedEvent.push_back(grInt[ch]);
   /* Use a modified Hann window for now */
   grWinPad[ch]     = evProcessTools::getWindowedAndPaddedEqualBeginGraph(grInt[ch], maxSamp, beginTime);
   /* The task of normalizing wf should be the responsibility of each reco method */
   cleanEvent.push_back(grWinPad[ch]);
  
   delete gr_v[ch];
   }//end of ch


   numSatChan = 0;   
   if(settings->nchnlFilter > 0){
   
   getNchnlMaskSat(unpaddedEvent, threshold, nchnlArray, chanMask, goodChan, numSatChan);

   if      (settings->recoPolType == "vpol") nchnl_tmp = nchnlArray[1];
   else if (settings->recoPolType == "hpol") nchnl_tmp = nchnlArray[2];
   else                                      nchnl_tmp = nchnlArray[0];

   if(nchnl_tmp < settings->nchnlCut){ 
   
      cerr<<"Failed nchnl cut. nchnl_tmp: "<<nchnl_tmp<<endl; 
      unpaddedEvent.clear();
      cleanEvent.clear();
      delete realAtriEvPtr;
      for(int ch=0; ch<16; ch++){ delete grInt[ch]; delete grWinPad[ch]; }      
      continue; 
   }
   }
   else {    
   for(int i=0; i<16; i++) goodChan[i] = chanMask[i];
   }
    
    getChannelSNR(unpaddedEvent, snrArray);
    TMath::Sort(16,snrArray,index);

    getChannelUnmodifiedSNR(unpaddedEvent, unmodSNRArray);
    TMath::Sort(16,unmodSNRArray,index);

    //recoData *summary = new recoData();
/*    
    if(settings->dataType == 0){
    summary->setWeight(event->Nu_Interaction[0].weight);
    summary->setTrueRadius(r_true);
    summary->setTrueDir(zen_true*180./M_PI, azi_true*180./M_PI);
    }
*/
    //summary->setOnion(onion);
    summary->setTopN(topN);
    summary->setRecoChan(goodChan);
    summary->setInWindowSNR(snrArray[index[2]]);
    summary->setUnmodSNR(unmodSNRArray[index[2]]);

    recoEventCnt++;



   if(settings->beamformMethod == 1){  
   if(settings->getSkymapMode == 0){ 
       err = reconstructCSW(settings, cleanEvent, &clEnv, recoDelays, recoDelays_V, recoDelays_H, nDir, chanMask, fitsFile/*argv[5]*/); 
   }
   else{
       maxPixIdx = reconstructCSWGetMaxPix(settings, cleanEvent, &clEnv, recoDelays, recoDelays_V, recoDelays_H, nDir, chanMask, summary);
   }
   } else {
   if(settings->getSkymapMode == 0){
 
      evStr = to_string(ev); 
      fitsFileStr = fitsFile_tmp + ".ev" + evStr + ".fits";
      sprintf(fitsFile, fitsFileStr.c_str());

      if(settings->skymapSearchMode == 0){ //no zoom mode
      maxPixIdx = reconstruct3DXCorrEnvelopeGetMaxPixAndMapData(settings, cleanEvent, &clEnv, recoDelays, recoDelays_V, recoDelays_H, goodChan, summary, fitsFile/*argv[5]*/, mapData/*, xCorrAroundPeakHist, sillygr*/);
      if(settings->recordMapData == 1){
      for(int pix=0; pix<nDir*nLayer; pix++) mapDataHist[pix]->Fill(mapData[pix]);
      }
      } else { //zoom search mode
      maxPixIdx = reconstruct3DXCorrEnvelopeGetMaxPix_ZoomMode(settings, cleanEvent, &clEnv, stationCenterDepth, antLocation, recoDelays, recoDelays_V, recoDelays_H, goodChan, summary, fitsFile);      
      }
   } else {
      maxPixIdx = reconstructXCorrEnvelopeGetMaxPix(settings, cleanEvent, &clEnv, recoDelays, recoDelays_V, recoDelays_H, nDir, goodChan/*chanMask*/, summary);
    
   }
  
   }
   if( err<0 || maxPixIdx<0){ cerr<<"Error reconstructing\n"; return -1; }

   //int recoFlag = record3DDiffGetFlag(summary, outputFile);
   //if( recoFlag ) recoFlagCnt++;
   summary->setFlag( (settings->skymapSearchMode)
                    ? record3DZoomedDiffGetFlag(settings, summary, dZenDist, dAziDist, recoTrueZenDist, recoTrueAziDist)
                    : record3DDiffGetFlag(settings, summary, dZenDist, dAziDist, recoTrueZenDist, recoTrueAziDist) );
   if(summary->flag > 0) recoFlagCnt++;
   maxPix[maxPixIdx]++;

   dataTree->Fill();

   unpaddedEvent.clear();
   cleanEvent.clear();
   delete realAtriEvPtr;
   //delete summary;
   for(int ch=0; ch<16; ch++){ delete grInt[ch]; delete grWinPad[ch]; }
   }//end of ev loop

   fp->Close();

} //if dataType == 1

else {

for (Long64_t ev=0; ev<numEntries; ev++){
   //cout<<"Entering event loop\n";
   summary->clear();
   //cout<<"Cleared previous summary\n";
   if(ev%1000 == 0) cout<<"*******************************Event got********************************: "<<ev<<endl;

   int string_i, antenna_i, AraRootChannel;

   chain2.GetEntry(ev);
   if(report->stations[0].Global_Pass > 0){ 
   
      trigEvCnt++;
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
      cout<<"r_true: "<<r_true<<" zen_true: "<<zen_true*180./M_PI<<" azi_true: "<<azi_true*180./M_PI<<endl;
          
   } else continue;
   
   double average[16]={0.};
 
   for(int a=0;a<16;a++)
   {
      string_i  = detector->getStringfromArbAntID(0, a);
      antenna_i = detector->getAntennafromArbAntID(0, a);
      AraRootChannel = detector->GetChannelfromStringAntenna (0, string_i, antenna_i, AraSim_settings);
      gr_v[AraRootChannel-1] = new TGraph();
      int nSamp = (int)report->stations[0].strings[string_i].antennas[antenna_i].time_mimic.size();

      for(int s=0; s<nSamp; s++){
         gr_v[AraRootChannel-1]->SetPoint(s
                                        , report->stations[0].strings[string_i].antennas[antenna_i].time_mimic[s]
                                        , report->stations[0].strings[string_i].antennas[antenna_i].V_mimic[s]
                                         );
        average[AraRootChannel-1] += report->stations[0].strings[string_i].antennas[antenna_i].V_mimic[s];
        }

      int pc = 0;
      /*** Zero-mean the waveforms on a channel-by-channel basis ***/
      if(gr_v[AraRootChannel-1]->GetN() != 0) average[AraRootChannel-1]/=(double)gr_v[AraRootChannel-1]->GetN();
      else cout<<"********Zero samples!!! *** ch:"<<AraRootChannel-1<<" n_samp: "<<gr_v[AraRootChannel-1]->GetN()<<endl;
    
      for(pc=0; pc<gr_v[AraRootChannel-1]->GetN(); pc++){
         gr_v[AraRootChannel-1]->GetPoint(pc, times, volts);
         gr_v[AraRootChannel-1]->SetPoint(pc, times, volts-average[AraRootChannel-1]);
      }

   }//End looping channels

   double wInt;
   int maxSamp;
  
   for(int ch=0; ch<16; ch++){
     
      gr_v[ch]->GetPoint(0,t,v);
      if( t<beginTime ) beginTime = t ;

   }

   for(int ch=0; ch<16; ch++){

   if(ch<8){wInt=0.5; /*maxSamp=2048;*/ maxSamp = gr_v[0]->GetN()*4; }
   else{    wInt=0.5; /*maxSamp=2048;*/ maxSamp = gr_v[8]->GetN()*4; }
   unpaddedEvent.push_back(gr_v[ch]);
   /* Use a modified Hann window for now */
   grWinPad[ch]     = evProcessTools::getWindowedAndPaddedEqualBeginGraph(gr_v[ch], maxSamp, beginTime);
   /* The task of normalizing wf should be the responsibility of each reco method */
   cleanEvent.push_back(grWinPad[ch]);

   }//end of ch

   numSatChan = 0;   
   if(settings->nchnlFilter > 0){
   
   getNchnlMaskSat(unpaddedEvent, threshold, nchnlArray, chanMask, goodChan, numSatChan);

   if      (settings->recoPolType == "vpol") nchnl_tmp = nchnlArray[1];
   else if (settings->recoPolType == "hpol") nchnl_tmp = nchnlArray[2];
   else                                      nchnl_tmp = nchnlArray[0];

   if(nchnl_tmp < settings->nchnlCut){ 
      
      cerr<<"Failed nchnl cut. nchnl_tmp: "<<nchnl_tmp<<endl; 
      unpaddedEvent.clear();
      cleanEvent.clear();
      for(int ch=0; ch<16; ch++){ delete gr_v[ch]; delete grWinPad[ch]; }      
      continue; 
 
   }
   }
   else {    
   for(int i=0; i<16; i++) goodChan[i] = chanMask[i];
   }
    
   getChannelSNR(unpaddedEvent, snrArray);
   TMath::Sort(16,snrArray,index);

   getChannelUnmodifiedSNR(unpaddedEvent, unmodSNRArray);
   TMath::Sort(16,unmodSNRArray,index);

   //recoData *summary = new recoData();
   //if(settings->dataType == 0){
   summary->setWeight(event->Nu_Interaction[0].weight);
   summary->setTrueRadius(r_true);
   summary->setTrueDir(zen_true*180./M_PI, azi_true*180./M_PI);
   //}
   //summary->setOnion(onion);
   summary->setTopN(topN);
   summary->setRecoChan(goodChan);
   summary->setInWindowSNR(snrArray[index[2]]);
   summary->setUnmodSNR(unmodSNRArray[index[2]]);

   recoEventCnt++;

   if(settings->beamformMethod == 1){  
   if(settings->getSkymapMode == 0){ 
       err = reconstructCSW(settings, cleanEvent, &clEnv, recoDelays, recoDelays_V, recoDelays_H, nDir, chanMask, fitsFile/*argv[5]*/); 
   }
   else{
       maxPixIdx = reconstructCSWGetMaxPix(settings, cleanEvent, &clEnv, recoDelays, recoDelays_V, recoDelays_H, nDir, chanMask, summary);
   }
   } else {
   if(settings->getSkymapMode == 0){

      evStr = to_string(ev); 
      fitsFileStr = fitsFile_tmp + ".ev" + evStr + ".fits";
      sprintf(fitsFile, fitsFileStr.c_str());

      if(settings->skymapSearchMode == 0){ //no zoom mode
      maxPixIdx = reconstruct3DXCorrEnvelopeGetMaxPixAndMapData(settings, cleanEvent, &clEnv, recoDelays, recoDelays_V, recoDelays_H, goodChan, summary, fitsFile/*argv[5]*/, mapData/*, xCorrAroundPeakHist, sillygr*/);
      if(settings->recordMapData == 1){
      for(int pix=0; pix<nDir*nLayer; pix++) mapDataHist[pix]->Fill(mapData[pix]);
      }
      } else {//zoom search mode
      maxPixIdx = reconstruct3DXCorrEnvelopeGetMaxPix_ZoomMode(settings, cleanEvent, &clEnv, stationCenterDepth, antLocation, recoDelays, recoDelays_V, recoDelays_H, goodChan, summary, fitsFile);  
      }
   } else {
      maxPixIdx = reconstructXCorrEnvelopeGetMaxPix(settings, cleanEvent, &clEnv, recoDelays, recoDelays_V, recoDelays_H, nDir, goodChan/*chanMask*/, summary);
    
   }
  
   }
   if( err<0 || maxPixIdx<0){ cerr<<"Error reconstructing\n"; return -1; }

   //int recoFlag = record3DDiffGetFlag(summary, outputFile);
   //if( recoFlag ) recoFlagCnt++;
   summary->setFlag( (settings->skymapSearchMode)
                    ? record3DZoomedDiffGetFlag(settings, summary, dZenDist, dAziDist, recoTrueZenDist, recoTrueAziDist)
                    : record3DDiffGetFlag(settings, summary, dZenDist, dAziDist, recoTrueZenDist, recoTrueAziDist) );
   
   if(summary->flag > 0) recoFlagCnt++;
   maxPix[maxPixIdx]++;
   dataTree->Fill();
   unpaddedEvent.clear();
   cleanEvent.clear();
   //delete summary;
   for(int ch=0; ch<16; ch++){ delete gr_v[ch]; delete grWinPad[ch]; }
   }//end of ev loop

}//end of dataType == 0

outputFile->Write();
outputFile->Close();

clfftTeardown();
err = tearDown(&clEnv);   
delete summary;
if(settings->skymapSearchMode == 0){
delete onion;
free(recoDelays);
free(recoDelays_V);
free(recoDelays_H);
}
delete settings;
free(mapDataHist);
free(mapData);

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
