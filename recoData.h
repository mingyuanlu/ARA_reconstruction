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
   //, double _interactionProbability, double _probability
   , int _eventTrigType
   , float zen_true, float azi_true, float zen_reco, float azi_reco, float r_true, float r_reco
   , float *_trueRecAngle, float *_trueLauAngle, float *_recoRecAngle, float *_recoLauAngle
   , int *usedChan
   , int idx, float xCorrValue
   //, Healpix_Onion *_onion
   , recoSettings *_settings
   , int _topN
   , int *_topMaxPixIdx, float *_topMaxPixCoherence
   , int *_maxPixIdxEachLayer, float *_maxPixCoherenceEachLayer
   , double _likelihood, double _pValue
   , float _inWindowSNR, float _inWindowSNR_V, float _inWindowSNR_H, float _unmodSNR
   , bool passAnotherPolNchnl
   , int _flag
   , trackEngine *_treg
   , int _constantNMaxPixIdx, float _constantNMaxPixCoherence
   , float _constantNZen, float _constantNAzi
   /*, int idx2, float xCorrValue2
   , int *_topMaxPixIdx2, float *_topMaxPixCoherence2
   , int *_maxPixIdxEachLayer2, float *_maxPixCoherenceEachLayer2*/);
   void setAllData2(recoSettings *_settings
   , int idx2, float xCorrValue2
   , int *_topMaxPixIdx2, float *_topMaxPixCoherence2
   , int *_maxPixIdxEachLayer2, float *_maxPixCoherenceEachLayer2);
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

//ClassDef 7
   int constantNMaxPixIdx;
   float constantNMaxPixCoherence;
   float constantNZen, constantNAzi;
   void setConstantNMaxPixInfo(int _constantNMaxPixIdx, float _constantNMaxPixCoherence);
   void setConstantNDir(float _constantNZen, float _constantNAzi);

//ClassDef 8
   float inWindowSNR_V;
   float inWindowSNR_H;
   bool passAnotherPolNchnl;
   void setInWindowSNRBothPol(float _inWindowSNR_V, float _inWindowSNR_H);
   void setPassAnotherPolNchnl(bool _passAnotherPolNchnl);

//ClassDef 9
   //Data members defined here are specific to reconstruction with the 2nd ray table.
   int maxPixIdx2;
   float maxPixCoherence2; //max pix coherence value
   vector<int> topMaxPixIdx2;               //top N max pix index of whole onion
   vector<float> topMaxPixCoherence2;       //top N max pix coherence of whole onion
   vector<int> maxPixIdxEachLayer2;         //max pix index of each layer
   vector<float> maxPixCoherenceEachLayer2; //max pix coherence of each layer
   void setTopMaxPix2Info(int *idx, float *xCorrValue);
   void setMaxPix2InfoEachLayer(recoSettings *settings, int *idx, float *xCorrValue);
   void setMaxPix2Info(int idx, float xCorrValue);

//ClassDef 10
   double interactionProbability;
   double probability; // survival probability (weight) * interaction probability
   void setProbabilities(double p_int, double p);

//ClassDef 11
   int numSatChan;
   int satChan[16];
   void setSaturatedChannels(int _numSatChan, int *_satChan);
   float channelInWindowSNR[16];
   void setChannelInWindowSNR(float *_channelInWindowSNR);

//ClassDef 12
   double survivalProbability; //AraSim Interaction::weight, represents the survival probability up to the vertex
   void setSurvivalProbability(double _survivalProbability);

//ClassDef 13
   vector<int> iterMaxPixIdx;
   vector<float> iterMaxPixCoherence;
   void setIterMaxPixInfo(recoSettings *settings, int *_iterMaxPixIdx, float *_iterMaxPixCoherence);

//ClassDef 14
   int maxFreqBin[16];
   double freqBinWidth_V, freqBinWidth_H;
   double maxCountFreq_V, maxCountFreq_H;
   void setMaxFreqBin(int *_maxFreqBin);
   void setMaxFreqBinByChannel(int ch, int _maxFreqBin);
   void setFreqBinWidth(double _freqBinWidth_V, double _freqBinWidth_H);
   void setMaxCountFreq(double _maxCountFreq_V, double _maxCountFreq_H);

   ClassDef(recoData, 14);

};


#endif
