
#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <fstream>

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
#include "TSystem.h"
#include "TMath.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"

#ifndef CALIBRATIONTOOLSVS3_H
#define CALIBRATIONTOOLSVS3_H

void invertGraph(TGraph *gr);
//TGraphErrors * calibrateElecChan(int elchan, UsefulAtriStationEvent * owncalEvPtr, double *TCorr, int stationId);
//double clearOffsetErr(TGraphErrors *gr);
double clearOffset(TGraph *gr, int stationId);

TGraph *calibrateTimeS(int elchan, UsefulAtriStationEvent *owncalEvPtr, double *tempCorr, int stationId);
TGraphErrors *calibrateTime(int elchan, UsefulAtriStationEvent *owncalEvPtr, double *tempCorr, int stationId);

TGraph * sortGraphS(TGraph * gr);
TGraphErrors * sortGraph(TGraphErrors *gr);


void loadVoltCalValues( int stationID );
void loadVoltCalErrors( int stationID );
void loadHighVoltCalValues( int stationID );
void loadBasicVoltCalValues( int stationID, int elchan );
void loadCableDelays(int stationId);

void loadTimeConstants(int stationId);

TGraph *calibrateTimeAltern(int elchan, UsefulAtriStationEvent * owncalEvPtr, double tempCorr, int stationId);

TGraphErrors * halfSamples(TGraphErrors *ingr, int evodd);
TGraph * halfSamplesS(TGraph *ingr, int evodd);
void writeGraphToFile(TGraph *ingr, string fileName);
void writeHistoToFile(TH1D * inh, string fileName);
int getTempCorr(TTree *eventTree, int runNumber, int stationID, int chipNumber, double *meanSqCfreq, AraEventCalibrator *calib, string rootBaseDir );
int selectSamples(TGraphErrors * SgrOwnCal, TGraphErrors * SgrOwnCal2);
int cutOverTimes(TGraphErrors * inGr, TGraphErrors * outGr, int halfFul=0);


TGraph * calibrateVoltsS(UsefulAtriStationEvent * owncalEvPtr, int chip, int elchan, TGraph *calGrin, int stationId);

TGraphErrors * calibrateVolts(UsefulAtriStationEvent *owncalEvPtr, int chip, int elchan, TGraphErrors *calGrin, int stationId);

TGraph * calibrateBasicVoltsS(UsefulAtriStationEvent * owncalEvPtr, int chip, int elchan, TGraph *calGrin, int stationId);


TGraphErrors * getVoltageError(UsefulAtriStationEvent * owncalEvPtr, int chip, int elchan, TGraphErrors *calGrin, int stationId);

#endif
