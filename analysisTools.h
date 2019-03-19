#ifndef ANALYSISTOOLS_H
#define ANALYSISTOOLS_H

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

#include "Detector.h"
#include "Trigger.h"
#include "Settings.h"
#include "Report.h"
#include "Event.h"
#include "Position.h"
#include "signal.hh"
#include "IceModel.h"

#define CORRUPT_EVENT_START_TIME 1448485911 //ARA02 run 3 start time
#define CORRUPT_EVENT_END_EVENT_NUMBER 4
#define ZEN_BAND_MAX -41
#define ZEN_BAND_MIN -46

using namespace std;

int getRunType(string STATION, int runNum);

class cutParameter /*: public TObject*/
{

public:



   double val;
   double plus;
   double minus;

   cutParameter();
   ~cutParameter();
   //ClassDef(cutParameter, 1);
};

class ARA02_cutValues /*: public TObject*/
{
private:

protected:

public:

   ARA02_cutValues();
   ~ARA02_cutValues();

   //CW impulsivity
   cutParameter cwImpCut[5];

   //Thermal box cut
   cutParameter coherenceCut_inBand[5];
   cutParameter snrCut_inBand[5];
   cutParameter coherenceCut_outOfBand[5];
   cutParameter snrCut_outOfBand[5];
   cutParameter snrCut[5];

   //Impulsivity cut
   cutParameter impCut;

   //Calpulser cut
   static const int nBoxes = 3;
   cutParameter zenMin[nBoxes];
   cutParameter zenMax[nBoxes];
   cutParameter aziMin[nBoxes];
   cutParameter aziMax[nBoxes];

   //Surface cut
   cutParameter surfaceCut_constantN;
   cutParameter surfaceCut_iterReco;

   void initialize();
   void setValue(cutParameter& param, double _val, double _plus, double _minus);

   //ClassDef(ARA02_cutValues, 1);

};

bool isCW_coincidence(bool &isVpolCW, bool &isHpolCW, int &maxCountFreqBin_V, int &maxCountFreqBin_H, recoData *dummyData, int cwBinThres);
bool isCW_freqWindow(bool &isVpolCW, bool &isHpolCW, bool& isXpolCW, recoData *dummyData, double fftRes);
bool isLowFreqDominance(int& lowFreqCount_V, int& lowFreqCount_H, recoData *dummyData, double highPassFreq, int lowFreqCountThres);
bool isThermal_boxCut(bool &inBand, recoSettings *settings, recoData *dummyData, Healpix_Onion onion, double snrCut_inBand, double coherenceCut_inBand, double snrCut_outOfBand, double coherenceCut_outOfBand);
bool isSurface(recoData *dummyData, double surfaceCut_1);
bool isIterSurface(double &zenMaj, recoData *dummyData, Healpix_Onion onion, recoSettings *settings, double zenRange, double surfaceCut_2);
float getZenMaj(const vector<float>& iterZenVec, float zenRange);
bool isNearNoisyRun(const vector<int>& noisyRuns, int runNum, int plusMinusRunNum);
bool isDeepPulser(string STATION, recoData *dummyData, int runNum);
bool isCalpulserTime(string STATION, recoData *dummyData);
bool isCalpulser(float &inBoxTheta, float &inBoxPhi, string STATION, recoData *dummyData, Healpix_Onion onion, recoSettings *settings, int type);
bool isRecoverableByImp(bool isVpolCW, bool isHpolCW, bool isXpolCW, recoData *dummyData, double impCut, double highPassFreq);
bool isBelowThermalImpulsivityCut(double &avgImpulsivity, recoData *dummyData, double postThermalAvgImpulsivityCut);
double getPercentile(vector<double> fftValues, const double percentile);
void getAngXingPixels(int& thetaPixCount, int& phiPixCount, recoData* dummyData, recoSettings* settings, Healpix_Onion onion, const double angThres);
void getAvgAngXingPixels(int& avgThetaPixCount, int& avgPhiPixCount, recoData* dummyData, recoSettings* settings, Healpix_Onion onion, const double angThres);
TH1F *getCumulative(TH1F *hist);
double getZenithInRangeFraction(recoData* dummyData, recoSettings* settings, Healpix_Onion onion, const double angThres);
double getAzimuthInRangeFraction(recoData* dummyData, recoSettings* settings, Healpix_Onion onion, const double angThres);
#endif
