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

#define maxSoftTriggerReadoutBlocks 16 //Software triggers are mostly 8-1 blocks, except it is 16-1
                                       //ARA02 run4879-4934, run5210-5277
                                       //ARA03 run3918-3974, run4009-4072
                                       //Since we are cutting away mistagged software triggers, we will just use the max number of
                                       //blocks - 15, with +1 as some leeway, to distinguish them from RF triggers. RF triggers have
                                       //typically 21-1 or 26-1 blocks.
#define IRS2SamplePerBlock 64 //Number of samples per IRS2 capacitor block

#define corruptFirst3EventStartTime_A2 1448485911 //ARA02 run 3 start time
#define corruptFirst3EventStartTime_A3 1448835608 //ARA02 run6002 start time
#define corruptEventEndEventNumber 4
#define corruptEventEndEventNumber_A3 5

#define ARA02D1CorruptionStartTime 1448485911 //ARA02 run 3, 11, 59, 60
#define ARA02D1CorruptionEndTime 1448954419

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

#define numDDA 4 //Number of DDAs in an ARA station
#define IRS2NumBlocks 512 //Number of blocks in an IRS2 chip

//void invertGraph(TGraph *gr);
int calibrateGeometryAndDelays(const RawAraStationEvent *rawAtriEvPtr,
                                double (&posDelayArray)[4][4], double *pulserCorr,
                                const float stationCenterDepth,
                                vector<vector<double> >& ant_loc,
                                vector<vector<double> >& pul_loc);
int getAraSimStationGeometry(vector<vector<double> >& ant_loc);
int getAraSimStationGeometry(vector<vector<double> >& ant_loc, Detector *detector, Settings *settings);

#endif
