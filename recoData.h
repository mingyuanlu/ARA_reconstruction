#ifndef RECODATA_H
#define RECODATA_H

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include "TObject.h"
#include "recoSettings.h"
#include "trackEngine.h"

using namespace std;

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
     int _eventId, int _eventNumber
   , int _unixTime, int _unixTimeUs, int _timeStamp
   , double w
   , int _eventTrigType
   , float zen_true, float azi_true, float zen_reco, float azi_reco, float r_true, float r_reco
   , float *_trueRecAngle, float *_trueLauAngle, float *_recoRecAngle, float *_recoLauAngle
   , int *usedChan
   , int idx, float xCorrValue
   //, Healpix_Onion *_onion
   , recoSettings *settings
   , int _topN
   , int *_topMaxPixIdx, float *_topMaxPixCoherence
   , int *_maxPixIdxEachLayer, float *_maxPixCoherenceEachLayer
   , double _likelihood, double _pValue
   , float _inWindowSNR, float _unmodSNR
   , int _flag
   , trackEngine *_treg);
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
   void setMaxPixInfoEachLayer(recoSettings *settings, int *idx, float *xCorrValue);
   void setLikelihoodAndPValue(double _likelihood, double _pValue);
   void setInWindowSNR(float _inWindowSNR);
   void setUnmodSNR(float _unmodSNR);
   void setFlag(int _flag);
   void duplicate(recoSettings *settings, recoData *old);
   void clear();

//ClassDef 2

   int eventTrigType; //Only applicable for real data. 0: RF trigger (excluding calpulser trigger), 1: Calpulser trigger, 2: Software trigger
   void setEventTrigType(int _eventTrigType);

//ClassDef 3

   int eventId;      //Only applicable for real data. AraRoot eventId.
   int eventNumber;  //Only applicable for real data. AraRoot eventNumber
   void setEventId(int _eventId);
   void setEventNumber(int _eventNumber);

//ClassDef 4

   int unixTime; //unixTime copied from AraRoot AraStationEventHeader_t
   int unixTimeUs; //unixTimeUs as above
   int timeStamp; //timeStamp as above
   void setEventTime(int _unixTime, int _unixTimeUs, int _timeStamp);

//ClassDef 5

   float recoRecAngle[16], recoLauAngle[16];
   float trueRecAngle[16], trueLauAngle[16];
   void setRecoAngles(float *recAngle, float *lauAngle);
   void setTrueAngles(float *recAngle, float *lauAngle);

//ClassDef 6
   trackEngine *treg;
   void setTreg(trackEngine *_treg);


   ClassDef(recoData, 6);
};


#endif
