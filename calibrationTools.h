#ifndef CALIBRATIONTOOLS_H
#define CALIBRATIONTOOLS_H

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

#include "Detector.h"
#include "Settings.h"
#include "Vector.h"

#include "TChain.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TTree.h"

using namespace std;

#define dropD4Time 1388840284 //Unix timestamp of ARA03 run1946 start time
                              //run1946 from 1388840284 - 1388861817
                              //This run covers 2014 Jan. 4 16:00 GMT. After this time ARA03 D4 should be discarded

//Unix timestamp of 2017 Jan 24 23:48- Jan 25 00:00 (UTC) IC deep pulser string 1 shallow pulser operation. This is in ARA2 run8573
#define deepPulserString1StartTime_2017 1485301680
#define deepPulserString1EndTime_2017 1485302400
//Unix time (sec) of 2017 Jan 25 00:01- Jan 25 00:19 (UTC) IC deep pulser string 22 shallow pulser operation. This is in ARA2 run8573
#define deepPulserString22StartTime_2017 1485302460
#define deepPulserString22EndTime_2017 1485303540
//Start/End timestamp of the above pulses
#define deepPulserString1StartTimeStamp_2017 27847600
#define deepPulserString1EndTimeStamp_2017 27854200
#define deepPulserString22StartTimeStamp_2017 39803400
#define deepPulserString22EndTimeStamp_2017 39809200

//void invertGraph(TGraph *gr);
int calibrateGeometryAndDelays(const RawAraStationEvent *rawAtriEvPtr,
                                double (&posDelayArray)[4][4], double *pulserCorr,
                                const float stationCenterDepth,
                                vector<vector<double> >& ant_loc,
                                vector<vector<double> >& pul_loc);
int calibrateGeometryAndDelaysPlusAdHocDepthOffSet(const RawAraStationEvent *rawAtriEvPtr,
                               double (&posDelayArray)[4][4], double *pulserCorr,
                               const float stationCenterDepth,
                               vector<vector<double> >& ant_loc,
                               vector<vector<double> >& pul_loc,
                               float adhocOffset);
int getAraSimStationGeometry(vector<vector<double> >& ant_loc);
int getAraSimStationGeometry(vector<vector<double> >& ant_loc, Detector *detector, Settings *settings);
int getSeckelStationGeometry(vector<vector<double> >& ant_loc);

#endif
