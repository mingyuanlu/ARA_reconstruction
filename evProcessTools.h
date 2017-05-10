#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <dirent.h>

using namespace std;

//AraRoot Includes
#include "RawIcrrStationEvent.h"

#include "RawAtriStationEvent.h"
#include "UsefulAraStationEvent.h"
#include "UsefulIcrrStationEvent.h"
#include "UsefulAtriStationEvent.h"

//Include FFTtools.h if you want to ask the correlation, etc. tools
#include "FFTtools.h"
#include "AraGeomTool.h"

//ROOT Includes
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TSystem.h"
#include "TMath.h"
#include "TF1.h"
#include "TChain.h"
#include "TH2D.h"
#include "TLine.h"
#include "TFitResult.h"
#include "TVirtualFitter.h"

//selfmade includes
#include "calibrationToolsVs3.h"
//For pdestal file search
#include <sys/types.h>
#include <dirent.h>
#include <vector>

#ifndef EVPROCESSTOOLS_H
#define EVPROCESSTOOLS_H


class evProcessTools{
public:

static std::vector<double> getMaximumCorrelationSum(
	std::vector<TGraph*> gr,
	std::vector<std::vector<double> > ant_loc,
	Double_t *corrSum,
	std::vector<std::vector<double> > *delayMaps,
	std::vector<std::vector< double > > *cohSumWave/*, TGraph * corrExample*/,
	int stationId,
	int plot
	);
static Double_t integratePower(TH1D* histo);
static std::vector < double > vecSum(std::vector< double > v1, std::vector< double > v2, int factor);
static std::vector < double > vecCrossProduct3D(std::vector< double > v1, std::vector< double > v2);
static double vecScalarProduct(std::vector< double > v1, std::vector< double > v2);
static double vecLength(std::vector< double > v1);
static std::vector< double > getXYProjection(std::vector< double > v);
static TGraph *getfftfromgraph(TGraph *gr, Double_t intSample, const int maxSample);
static TGraph *getfftdbfromgraph(TGraph *gr, Double_t intSample, const int maxSample);
static TH1D *makeffthisto(TGraph *gr, char* name, Double_t intSample);
static double getTotalPower(std::vector<TGraph*> gr, std::vector<std::vector<double> > ant_loc);
static std::vector< vector< double > > getFFT(std::vector<TGraph*> gr, std::vector<std::vector<double> > ant_loc, vector<vector< double > > * FFTpower, vector<vector<double> > * FFTfreq, vector<vector< double > > * normFFT,  vector<vector< double > > refArray, double *cwValueV, double *cwValueH, vector<vector< double > > * rayleighV, vector<vector<double> > * rayleighF);
static double getCorehentlySummedWaveAmp(std::vector<TGraph *> gr, vector<vector< double > >delayMap, double intSample, std::vector<std::vector< double > > *cohSumWave);
static std::vector<double> getCrudeReconstruction(vector<vector< double > > delayMap, std::vector<std::vector<double> > ant_loc);
static TGraph *cwCleaner(TGraph *grWaveIn, Double_t minFreq, Double_t maxFreq, double reduction);
static TGraph *getBartlettAndPaddedGraph(TGraph *gr, const int maxSample);
static TGraph *getWindowedAndPaddedEqualBeginGraph(TGraph *gr, const int maxSample, const double beginTime);
static TGraph *getNormalizedGraph(TGraph *gr);
static TGraph *universalCleaner(TGraph *grWaveIn, std::vector< double> cwSingleVec );
static double getSimpleFFT(std::vector<TGraph*> gr, vector<vector< double > > * FFTpower, vector<vector<double> > * FFTfreq);
static int peakFinder(TGraph *gr, int numberOfPeaks, int range, double *peakPos, double *peakHight, double threshold, int precision);
static double getMeanBaseLine(TGraph *gr);
static int spikingString(std::vector<TGraph*> gr, int stationId);
static double modifiedHannWindow(int idx, int numSample, int modFrac);
//static double getPeakSqValRange(TGraph *gr, int *index, int firstBin, int lastBin);

};
#endif
