#ifndef RECOTOOLS_H
#define RECOTOOLS_H

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <cmath>
#include <time.h>
//#include <chrono>

/* OpenCL/CLFFT includes */
//#include <CL/opencl.h>
#include <clFFT.h>

/* radiospline includes */
#include "radiospline/IceGeometry.h"
#include "radiospline/RayDelay.h"
#include "radiospline/FirnShadow.h"

//#define AIR_FILE "InAir1.fits"
//#define ICE_FILE "InIce1.fits"
//#define AIR_FILE "InAirFinal_mod.fits"
//#define ICE_FILE "InIceFinal_mod.fits"
//#define SHADOW_FILE "shadow_5km.fits"
#define AIR_FILE "delay_inair.fits"
#define ICE_FILE "delay_inice.fits"
#define SHADOW_FILE "firn_shadow.fits"

/* CFITSIO include */
#include "fitsio.h"

/* Healpix includes */
/*
#include "healpix_base.h"
#include "healpix_data_io.h"
#include "healpix_map_fitsio.h"
#include "healpix_map.h"
#include "fitshandle.h"
#include "pointing.h"

#include "vec3.h"
#include "arr.h"
*/
/* ROOT includes */
//#include "TGraph.h"
//#include "TCanvas.h"
//#include "TGaxis.h"
//#include "TObject.h"

#include "FFTtools.h"
#include "evProcessTools.h"
//#include "recoData.h"
#include "recoSettings.h"
#include "recoData.h"
#include "Healpix_Onion.h"
#define NPIX_NESTED 4 //number of pixels for nSideExp = n+1 within a single pixel of nSideExp = n

#ifndef TRACKENGINE_H
#define C_INV 3.34
#define speedOfLight 0.3 // m/ns
#define nIce 1.76 // ~-180m
#endif

//#define REFERENCE_MAP_FIT_FILE "testFitFuntFile.root"
//#define REFERENCE_MAP_FIT_FILE "testFitFuncFile_2013_A3_pulserSweep_run410-430.root"

//#define REFERENCE_MAP_FIT_FILE "testFitFuncFile_2014_ARA03_run3623.root"
//#define FIT_FUNC "gaus"
/* If bandpass defined, waveforms will be bandpassed before being cross-correlated in the freq domain */
//#define bandpass

//class recoSettings;
//class Healpix_Onion;

using namespace std;
/*
#ifndef XCORRSUMGRAPH
TGraph *sillygr = new TGraph();
//TGraph *envelopeSum = new TGraph();
#define XCORRSUMGRAPH
#endif
*/
struct recoEnvData{

   cl_platform_id   *platforms;
   cl_device_id     *devices;

   cl_context        context;
   cl_command_queue  queue;
   cl_program        program;
/*
   static cl_context        context;
   static cl_command_queue  queue;
   static cl_program        program;
*/

   cl_kernel        *kernels;

   /* CSW kernels */
   cl_kernel         shiftWf;
   cl_kernel         sumWf;
   cl_kernel         wfPwr;

   /* XCorr kernels */
   cl_kernel         xCorrWf;
   cl_kernel         computeXCorrCoef;
   cl_kernel         computeCoherence;
   cl_kernel         computeNormalizedCoherence;
   cl_kernel         computeIterativeNormalizedCoherence;

   /* Bandpass kernel */
   cl_kernel         bandPassFilter;

   /* 3D reco: get info from Healpix_Onion kernel */
   cl_kernel         getMaxPixInfoEachLayer;


/*
    static cl_kernel         shiftwf;
    static cl_kernel         sumwf;
    static cl_kernel         wfpwr;
*/
   clfftPlanHandle  planHandle;
   clfftSetupData   fftSetup;
};
/*
struct timer{

   vector<clock_t> clockVec;
   //vector<std::chrono::monotonic_clock::time_point> chronoVec;
   vector<time_t> timeVec;


};
*/
/*
class Healpix_Onion : public TObject
{
private:

protected:

public:



   int nSideExp;
   int nLayer;
   int nDir;
   Healpix_Base hpBase;
   vector<float> layerRadii;

   Healpix_Onion(){ initialize(); }
   ~Healpix_Onion(){ }

   void initialize(){
   nSideExp = 2;
   nLayer   = 1;
   hpBase = Healpix_Base(pow(2,nSideExp), RING, SET_NSIDE);
   nDir = hpBase.Npix();
   layerRadii.clear();
   layerRadii.push_back(3000);

   }

   //~Healpix_Onion();

   Healpix_Onion(int _nSideExp, int _nLayer){

   nSideExp = _nSideExp;
   nLayer   = _nLayer;
   hpBase = Healpix_Base(pow(2,nSideExp), RING, SET_NSIDE);
   nDir = hpBase.Npix();
   layerRadii.clear();
   float r=0.f;

   //while( r<= 5000.f ){

   //   r += 5000.f / (float)nLayer;
   //   layerRadii.push_back( r );

   //}


   for(int i=1; i<=nLayer; i++){

   //layerRadii.push_back( (float)i * 5000.f / (float)nLayer);
   layerRadii.push_back( (float)i * 3000.f / (float)nLayer); //for testing 3D calpulser reco

   }

   if( (int)layerRadii.size() != nLayer ) cerr<<"Warning!! layerRadii.size(): "<<layerRadii.size()<<" nLayer: "<<nLayer<<endl;

   }

   pointing getPointing(int pixNum){

   return hpBase.pix2ang( pixNum % nDir );

   }

   int getLayerNumber(int pixNum){

   return pixNum / nDir ;

   }

   float getLayerRadius(int pixNum){

   return layerRadii[ getLayerNumber(pixNum) ];

   }

   ClassDef(Healpix_Onion, 1);
};
*/
/*
class recoData : public TObject
{
private:

protected:

public:

   double weight;
   // Ideally in degrees
   float trueZen, trueAzi, recoZen, recoAzi;
   float trueRadius, recoRadius;
   int recoChan[16]; //channels actually used in the reco
   int maxPixIdx;
   float maxPixCoherence; //max pix coherence value
   int topN;                               //size of topMaxPixIdx
   //Healpix_Onion *onion;
   vector<int> topMaxPixIdx;               //top N max pix index of whole onion
   vector<float> topMaxPixCoherence;       //top N max pix coherence of whole onion
   vector<int> maxPixIdxEachLayer;         //max pix index of each layer
   vector<float> maxPixCoherenceEachLayer; //max pix coherence of each layer
   double likelihood;  //likelihood of the whole skymap compared to the reference map, presumably constructed from thermal events
   double pValue;      //p value of the whole skymap compared to the reference map, presumably constructed from thermal events
   float inWindowSNR;
   float unmodSNR;

   int flag; //whether this event is flagged or not, depending on the flagging condition specified in record3DDiffGetFlag()

   recoData();
   ~recoData();

   void initialize();
   void setAllData(
     double w
   , float zen_true, float azi_true, float zen_reco, float azi_reco, float r_true, float r_reco
   , int *usedChan
   , int idx, float xCorrValue
   , Healpix_Onion *_onion
   , int _topN
   , int *_topMaxPixIdx, float *_topMaxPixCoherence
   , int *_maxPixIdxEachLayer, float *_maxPixCoherenceEachLayer
   , double _likelihood, double _pValue
   , float _inWindowSNR, float _unmodSNR
   , int _flag);
   void setWeight(double w);
   void setTrueRadius(float r_true);
   void setTrueDir(float zen_true, float azi_true);
   void setRecoRadius(float r_reco);
   void setRecoDir(float zen_reco, float azi_reco);
   void setRecoChan(int *usedChan);
   void setMaxPixInfo(int idx, float xCorrValue);
   //void setOnionInfo(Healpix_Onion *_onion);
   void setTopN(int _topN);
   void setTopMaxPixInfo(int *idx, float *xCorrValue);
   void setMaxPixInfoEachLayer(int *idx, float *xCorrValue);
   void setLikelihoodAndPValue(double _likelihood, double _pValue);
   void setInWindowSNR(float _inWindowSNR);
   void setUnmodSNR(float _unmodSNR);
   void setFlag(float _flag);
   void duplicate(recoData *old);
   void clear();

   ClassDef(recoData, 1);
};
*/
float getMeanDelay( vector<float>& solvedDelay);
float getMeanDelay_passByValue( vector<float> solvedDelay);
int setupCLRecoEnv(recoSettings *settings, recoEnvData *clEnv, const char *programFile);
int reconstructCSW(unsigned int dataType, vector<TGraph *>& cleanEvent, recoEnvData *clEnv,
                float *recoDelays, float *recoDelays_V, float *recoDelays_H,
                int nDir, string pol, const int *chanMask, char *filename);
int reconstructCSW(recoSettings *settings, vector<TGraph *>& cleanEvent, recoEnvData *clEnv,
                float *recoDelays, float *recoDelays_V, float *recoDelays_H,
                int nDir, const int *chanMask, char *filename);
int reconstructCSW_Serial(unsigned int dataType, vector<TGraph *>& cleanEvent,
                float *recoDelays, float *recoDelays_V, float *recoDelays_H,
                int nDir, string pol, const int *chanMask, char *filename);
int reconstructXCorr(unsigned int dataType, vector<TGraph *>& cleanEvent, recoEnvData *clEnv,
                    float *recoDelays, float *recoDelays_V, float *recoDelays_H,
                    int nDir, string pol, const int *chanMask, char *filename);
int reconstructXCorrEnvelope(unsigned int dataType, vector<TGraph *>& cleanEvent, recoEnvData *clEnv,
                    float *recoDelays, float *recoDelays_V, float *recoDelays_H,
                    int nDir, string pol, const int *chanMask, char *filename);
int reconstructCSWGetMaxPix(unsigned int dataType, vector<TGraph *>& cleanEvent, recoEnvData *clEnv,
                float *recoDelays, float *recoDelays_V, float *recoDelays_H,
                int nDir, string pol, const int *chanMask, recoData *summary);
int reconstructCSWGetMaxPix(recoSettings *settings, vector<TGraph *>& cleanEvent, recoEnvData *clEnv,
                float *recoDelays, float *recoDelays_V, float *recoDelays_H,
                int nDir, const int *chanMask, recoData *summary);
int reconstructXCorrGetMaxPix(unsigned int dataType, vector<TGraph *>& cleanEvent, recoEnvData *clEnv,
                    float *recoDelays, float *recoDelays_V, float *recoDelays_H,
                    int nDir, string pol, const int *chanMask, recoData *summary);
int reconstructXCorrEnvelopeGetMaxPix(unsigned int dataType, vector<TGraph *>& cleanEvent, recoEnvData *clEnv,
                    float *recoDelays, float *recoDelays_V, float *recoDelays_H,
                    int nDir, string pol, const int *chanMask, recoData *summary);
int reconstructXCorrEnvelopeGetMaxPix(recoSettings *settings, vector<TGraph *>& cleanEvent, recoEnvData *clEnv,
                    float *recoDelays, float *recoDelays_V, float *recoDelays_H,
                    int nDir, const int *chanMask, recoData *summary);
int reconstructXCorrGetMaxPixAndMap(unsigned int dataType, vector<TGraph *>& cleanEvent, recoEnvData *clEnv,
                    float *recoDelays, float *recoDelays_V, float *recoDelays_H,
                    int nDir, string pol, const int *chanMask, recoData *summary, char *filename);
int reconstructXCorrEnvelopeGetMaxPixAndMap(unsigned int dataType, vector<TGraph *>& cleanEvent, recoEnvData *clEnv,
                    float *recoDelays, float *recoDelays_V, float *recoDelays_H,
                    int nDir, string pol, const int *chanMask, recoData *summary, char *filename);
/*
int reconstruct3DXCorrEnvelopeGetMaxPixAndMap(unsigned int dataType, vector<TGraph *>& cleanEvent, recoEnvData *clEnv,
                    float *recoDelays, float *recoDelays_V, float *recoDelays_H,
                    string pol, const int *chanMask, recoData *summary, char *filename,
                    TH1F *xCorrAroundPeakHist[]);
int reconstruct3DXCorrEnvelopeGetMaxPixAndMapData(unsigned int dataType, vector<TGraph *>& cleanEvent, recoEnvData *clEnv,
                    float *recoDelays, float *recoDelays_V, float *recoDelays_H,
                    string pol, const int *chanMask, recoData *summary, char *filename, float *mapData,
                    TH1F *xCorrAroundPeakHist[], TGraph *sillygr[]);
*/
int reconstruct3DXCorrEnvelopeGetMaxPixAndMapData(recoSettings *settings, vector<TGraph *>& cleanEvent, recoEnvData *clEnv,
                    float *recoDelays, float *recoDelays_V, float *recoDelays_H, const int *chanMask, int *index,
                    recoData *summary, char *filename, float *mapData);
int reconstruct3DXCorrEnvelopeGetMaxPixAndMapData(recoSettings *settings, vector<TGraph *>& cleanEvent, recoEnvData *clEnv,
                    float *recoDelays, float *recoDelays_V, float *recoDelays_H, const int *chanMask, /*int *index,*/
                    recoData *summary, char *filename, float *mapData);
int reconstruct3DXCorrEnvelopeGetMaxPixAndMapData_2ndRayReco(recoSettings *settings, vector<TGraph *>& cleanEvent, recoEnvData *clEnv,
                    float *recoDelays, float *recoDelays_V, float *recoDelays_H, const int *chanMask,
                    recoData *summary, char *filename, float *mapData);
int reconstruct3DXCorrEnvelopeGetMaxPixAndMapData_constantNFilter(recoSettings *settings, vector<TGraph *>& cleanEvent,
                                                                  recoEnvData *clEnv, float *recoDelays, float *recoDelays_V, float *recoDelays_H, const int *chanMask, recoData *summary, char *filename/*, float *mapData*/);
int reconstruct3DXCorrEnvelopeGetMaxPix_ZoomMode(recoSettings *settings, vector<TGraph *>& cleanEvent, recoEnvData *clEnv,
                    const float stationCenterDepth, const vector<vector<double> >& antLocation,
                    float *recoDelays, float *recoDelays_V, float *recoDelays_H, const int *chanMask,
                    recoData *summary, char *filename);

int tearDown(recoEnvData *clEnv);
int computeRecoDelaysWithConstantN(const int nAnt, const float zCenter, const vector<vector<double> >& antLoc,
                                     //const float radius, const int nSideExp,
                                     Healpix_Onion *onion,
                                     float *recoDelays, float *recoDelays_V, float *recoDelays_H);
int computeRecoDelaysWithNoBoundConstantN(const int nAnt, const float zCenter, const vector<vector<double> >& antLoc,
                                     //const float radius, const int nSideExp,
                                     Healpix_Onion *onion,
                                     float *recoDelays, float *recoDelays_V, float *recoDelays_H);
int computeRecoDelaysWithConstantNForSinglePixel(const int nAnt, const float zCenter, const vector<vector<double> >& antLoc,
                                     //const float radius, const int nSideExp,
                                     Healpix_Onion *onion,
                                     float *recoDelays, float *recoDelays_V, float *recoDelays_H,
                                     const int pix);

int computeZoomedRecoDelaysWithConstantN(const int nAnt, const float zCenter, const vector<vector<double> >& antLoc,
                                     //const float radius, const int nSideExp,
                                     Healpix_Onion *onion,
                                     float *recoDelays, float *recoDelays_V, float *recoDelays_H,
                                     const int last2DMaxPixIdx, const int nPix_nested);
int computeRecoDelaysWithRadioSpline(const int nAnt, const float zCenter, const vector<vector<double> >& antLoc,
                                     const float radius, const int nSideExp,
                                     //Healpix_Onion *onion,
                                     float *recoDelays, float *recoDelays_V, float *recoDelays_H);
//int computeRecoDelaysWithConstantN(const int nAnt, const float zCenter, const vector<vector<double> >& antLoc,
//                                     const float radius, const int nSideExp,
//                                     float *recoDelays, float *recoDelays_V, float *recoDelays_H);
int compute3DRecoDelaysWithRadioSpline(const int nAnt, const float zCenter, const vector<vector<double> >& antLoc,
                                      //const float radius, const int nSideExp,
                                      Healpix_Onion *onion,
                                      float *recoDelays, float *recoDelays_V, float *recoDelays_H);
int compute3DRecoBothDelaysWithRadioSpline(const int nAnt, const float zCenter, const vector<vector<double> >& antLoc,
                                     //const float radius, const int nSideExp,
                                     Healpix_Onion *onion,
                                     float *recoDelays, float *recoDelays_V, float *recoDelays_H,
                                     float *recoRefracDelays, float *recoRefracDelays_V, float *recoRefracDelays_H);
int compute3DRecoDelaysWithRadioSplineForSinglePixel(const int nAnt, const float zCenter, const vector<vector<double> >& antLoc,
                                      //const float radius, const int nSideExp,
                                      Healpix_Onion *onion,
                                      float *recoDelays, float *recoDelays_V, float *recoDelays_H,
                                      const int pix);

int compute3DRecoAnglesWithRadioSplineForSinglePixel(const int nAnt, const float zCenter, const vector<vector<double> >& antLoc, Healpix_Onion *onion, float *recoLauAngles, float *recoRecAngles, const int pix);

int compute3DRecoAnglesWithRadioSplineForSinglePoint(const int nAnt, const float zCenter, const vector<vector<double> >& antLoc, float *recoLauAngles, float *recoRecAngles, const double *srcLoc);

int compute3DZoomedRecoDelaysWithRadioSpline(const int nAnt, const float zCenter, const vector<vector<double> >& antLoc,
                                     //const float radius, const int nSideExp,
                                     Healpix_Onion *onion,
                                     float *recoDelays, float *recoDelays_V, float *recoDelays_H,
                                     const int last2DMaxPixIdx, const int nPix_nested);

int plotMaxPix(const int nDir, int *maxPix, char *filename);
int plotMaxPixZenAzi(const int nSideExp, int *maxPix, char *rootFilename);
int getMaxBin(TGraph *gr);
void setMeanAndSigmaInNoMax(TGraph *gr, double *stats);
void getNchnl(const vector<TGraph *>& cleanEvent, double threshold, int *nchnlArray);
void getNchnl(recoSettings *settings, const vector<TGraph *>& cleanEvent, double threshold, int *nchnlArray, int *goodChan);
void getNchnlMask(const vector<TGraph *>& cleanEvent, double threshold, int *nchnlArray, const int *chanMask, int *goodChan);
void getNchnlMaskSat(const vector<TGraph *>& cleanEvent, double threshold, int *nchnlArray, const int *chanMask, int *goodChan, int& numSatChan);
int getSaturation(recoSettings* settings, const vector<TGraph *>& cleanEvent, int *satChan);
void getChannelSNR(const vector<TGraph *>& cleanEvent, float *snrArray);
void getChannelUnmodifiedSNR(const vector<TGraph *>& cleanEvent, float *snrArray);
//int recordDiff(int nSideExp, int maxPixIdx, float maxPixValue, double weight, float zen_true, float azi_true, float r_true, int *usedChan, char *rootFilename);
//int recordDiffGetFlag(int nSideExp, int maxPixIdx, float maxPixValue, double weight, float zen_true, float azi_true, float r_true, int *usedChan, char *rootFilename);
//int recordDiff(int nSideExp, recoData *summary, char *rootFilename);
int recordDiffGetFlag(int nSideExp, recoData *summary, char *rootFilename);
//int record3DDiffGetFlag(recoData *summary, char *rootFilename);
int record3DDiffGetFlag(recoSettings *settings, recoData *summary, TH1F *dZenDist, TH1F *dAziDist, TH2F *recoTrueZenDist, TH2F *recoTrueAziDist);
int record3DDiffGetFlag_2ndRayReco(recoSettings *settings, recoData *summary, TH1F *dZenDist, TH1F *dAziDist, TH2F *recoTrueZenDist, TH2F *recoTrueAziDist);
int recordConstantNDir(recoSettings *settings, recoData *summary);
int record3DZoomedDiffGetFlag(recoSettings *settings, recoData *summary, TH1F *dZenDist, TH1F *dAziDist, TH2F *recoTrueZenDist, TH2F *recoTrueAziDist);
float getSpaceAngle(float theta1, float phi1, float theta2, float phi2);
void stackXCorrAroundPeak(const TGraph *gr, TH1F *hist, float plusMinusTime);
int doNchnlScan(const double eventWeight, const vector<TGraph *>& cleanEvent, TH2F *mnMap, int *nchnlArray, const int *chanMask, int *goodChan, const int nThresStep = 250, const double minThres = 0., const double maxThres = 10.);
int doNchnlScan(const double eventWeight, const vector<TGraph *>& cleanEvent, TH2F *mnMap, TH2F *mnMap_V, TH2F *mnMap_H, int *nchnlArray, const int *chanMask, int *goodChan, const int nThresStep = 250, const double minThres = 0., const double maxThres = 10.);
int computeMapLikelihoodAndPValue(const int nDir, const int nLayer, const char *fitFunc, const char *fitFuncFile, float *mapData, double& likelihood, double& pValue);
double getPeakSqValRange(TGraph *gr, int *index, int firstBin, int lastBin);
//int recordTime(timer *t, const int count);
//void resetTimer(timer *t);
double impulsivityMeasure(TGraph * wf, TGraph * distance_cdf = 0, int pt = -1/*, bool    hilbert*/);
TGraph *impulsivityMeasure(TGraph * wf, double *impulsivity/*TGraph * distance_cdf*/,int pt = -1/*, bool    hilbert*/);
TGraph *bipolarnessMeasure(TGraph *wf, double *bipolarness/*, TGraph *grCumuSumCDF*/);
void getPosNegPowerPeakAndDeltaT(TGraph *wf, double *posPowerPeak, double *negPowerPeak, double *deltaT, int nIntSamp=50);
void getChannelSlidingV2SNR(const vector<TGraph *>& cleanEvent, int nIntSamp, float *snrArray);
void getChannelTotalPowerSNR(const vector<TGraph *>& cleanEvent, int nIntSamp, float *snrArray);
#endif
