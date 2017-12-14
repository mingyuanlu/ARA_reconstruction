#ifndef TRACKENGINE_H
#define TRACKENGINE_H

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include "TObject.h"
#include "Vector.h"

#include "RawIcrrStationEvent.h"
#include "RawAtriStationEvent.h"
#include "UsefulAraStationEvent.h"
#include "UsefulIcrrStationEvent.h"
#include "UsefulAtriStationEvent.h"
#include "AraGeomTool.h"

#include "recoTools.h"

using namespace std;

class trackEngine: public TObject
{

private:

protected:

public:

   trackEngine();
   ~trackEngine();

   const Vector zeroVector = Vector(0.,0.,0.);
   vector< vector<Vector> > baselineTracks;
   vector< vector<double> > baselineTrackTimes;
   vector< vector<Vector> > eventOrthoTracks; //conjecture tracks, each orthogonal to its baseline track;
   //vector<union> rank;

   double tolerance;
   void setTolerance(double tol);
/*
   int buildBaselineTracks(const Detector *detector);
   int buildBaselineTracks(const RawAraStationEvent *rawEvPtr);
*/
   int getChannelPeakTime(vector<TGraph *> unpaddedEvent, double *peakTArray);

   int buildBaselineTracks(const vector< vector<double> >& antLocation);
   double checkBaselineTracks();
   int computeCosines(vector<TGraph*> unpaddedEvent);
   int buildEventOrthoTracks(Vector finalTrack);

   int goodTrackCount, badTrackCount;
   int *trackRank;
   float *trackSNRArray;
   double *cosine;

   Vector demoFinalTrack;
   Vector hierFinalTrack;
   Vector demoExtrapFinalTrack;
   Vector hierExtrapFinalTrack;
   Vector iterDemoExtrapFinalTrack;
   Vector iterHierExtrapFinalTrack;

   int computeFinalTracks(vector<TGraph *> unpaddedEvent);
   //int computeExtrapFinalTracksWithDemo();
   //int computeExtrapFinalTracksWithHier();
   //int computeExtrapFinalTracks(Vector finalTrack);
   Vector computeDemoExtrapFinalTrack();
   Vector computeHierExtrapFinalTrack();


   double angleThreshold;
   int maxIteration;
   //int computeIterExtrapFinalTracksWithDemo();
   //int computeIterExtrapFinalTracksWithHier();
   int computeIterExtrapFinalTracks(Vector tempDemoExtrap, Vector tempHierExtrap);
   void setAngleThreshold(double ang_thres);
   void setMaxIteration(int max_iter);

   void initialize();
   void clearForNextEvent();
   void clearAll();

   int computeAllTracks(/*const vector< vector<double> >& antLocation,*/ vector<TGraph *> unpaddedEvent);

   void print();

   //ClassDef(trackEngine, 1);

};

#endif
