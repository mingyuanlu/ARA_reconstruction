#include "recoSettings.h"
#include "recoTools.h"
#include "TObject.h"

ClassImp(recoSettings)
ClassImp(recoData);

using namespace std;

   recoData::recoData(){ recoData::initialize(); }
   recoData::~recoData(){ /* default destructor */ }

   void recoData::initialize(){

   weight = 0.;
   trueZen = trueAzi = recoZen = recoAzi = 0.f;
   trueRadius = recoRadius = 0.f;
   std::fill(&recoChan[0], &recoChan[16], 0);
   maxPixIdx = 0;
   maxPixCoherence = 0.f;
   //onion = NULL;
   topN = 0;
   topMaxPixIdx.clear();
   topMaxPixCoherence.clear();
   maxPixIdxEachLayer.clear();
   maxPixCoherenceEachLayer.clear();
   likelihood = 0.;
   pValue = 0.;
   inWindowSNR = 0.f;
   unmodSNR = 0.f;
   flag = 0;

   eventTrigType = 0;

   eventId = 0;
   eventNumber = 0;

   unixTime = 0;
   unixTimeUs = 0;
   timeStamp = 0;

   std::fill(&recoRecAngle[0], &recoRecAngle[16], -1.);
   std::fill(&recoLauAngle[0], &recoLauAngle[16], -1.);
   std::fill(&trueRecAngle[0], &trueRecAngle[16], -1.);
   std::fill(&trueLauAngle[0], &trueLauAngle[16], -1.);

   }



   void recoData::setAllData(
     int _eventId, int _eventNumber
   , int _unixTime, int _unixTimeUs, int _timeStamp
   , double w
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
   , float _inWindowSNR, float _unmodSNR
   , int _flag)
   {

   eventId = _eventId;
   eventNumber = _eventNumber;
   unixTime = _unixTime;
   unixTimeUs = _unixTimeUs;
   timeStamp = _timeStamp;
   weight = w;
   trueZen = zen_true;
   trueAzi = azi_true;
   recoZen = zen_reco;
   recoAzi = azi_reco;
   trueRadius = r_true;
   recoRadius = r_reco;

   for(int i=0; i<16; i++){
      recoChan[i] = usedChan[i];
      trueRecAngle[i] = _trueRecAngle[i];
      trueLauAngle[i] = _trueLauAngle[i];
      recoRecAngle[i] = _recoRecAngle[i];
      recoLauAngle[i] = _recoLauAngle[i];
   }

   maxPixIdx = idx;
   maxPixCoherence = xCorrValue;

   //setOnionInfo(_onion);
   //onion = _onion;
   topN = _topN;

   topMaxPixIdx.clear();
   topMaxPixCoherence.clear();
   for(int i=0; i<topN; i++){
      topMaxPixIdx.push_back( _topMaxPixIdx[i] );
      topMaxPixCoherence.push_back( _topMaxPixCoherence[i] );
   }

   maxPixIdxEachLayer.clear();
   maxPixCoherenceEachLayer.clear();
   for(int i=0; i<_settings->nLayer; i++){
      maxPixIdxEachLayer.push_back( _maxPixIdxEachLayer[i] );
      maxPixCoherenceEachLayer.push_back( _maxPixCoherenceEachLayer[i] );
   }

   likelihood = _likelihood;
   pValue = _pValue;

   inWindowSNR = _inWindowSNR;
   unmodSNR = _unmodSNR;

   flag = _flag;

   eventTrigType = _eventTrigType;

   }

   void recoData::setWeight(double w){

   weight = w;

   }

   void recoData::setTrueRadius(float r_true){

   trueRadius = r_true;

   }

   void recoData::setTrueDir(float zen_true, float azi_true){

   trueZen = zen_true;
   trueAzi = azi_true;

   }

   void recoData::setRecoRadius(float r_reco){

   recoRadius = r_reco;

   }

   void recoData::setRecoDir(float zen_reco, float azi_reco){

   recoZen = zen_reco;
   recoAzi = azi_reco;

   }

   void recoData::setRecoChan(int *usedChan){

   for(int i=0; i<16; i++) recoChan[i] = usedChan[i];

   }

   void recoData::setRecoAngles(float *recAngle, float *lauAngle){

   for(int i=0; i<16; i++){
     recoRecAngle[i] = recAngle[i];
     recoLauAngle[i] = lauAngle[i];
   }

   }
   void recoData::setTrueAngles(float *recAngle, float *lauAngle){

   for(int i=0; i<16; i++){
     trueRecAngle[i] = recAngle[i];
     trueLauAngle[i] = lauAngle[i];
   }

   }

   void recoData::setMaxPixInfo(int idx, float xCorrValue){

   maxPixIdx = idx;
   maxPixCoherence = xCorrValue;

   }
   /*
   void recoData::setOnionInfo(Healpix_Onion *_onion){

   onion = _onion;

   }
   */
   void recoData::setTopN(int _topN){

   topN = _topN;

   }

   void recoData::setTopMaxPixInfo(int *idx, float *xCorrValue){ //pass a sorted array of topN max pix indices

   //if( topN != (int)idx_size ) cerr<<"Warning!! topMaxPixIdx size not correct\n";
   topMaxPixIdx.clear();
   topMaxPixCoherence.clear();
   for(int i=0; i<topN; i++){
      topMaxPixIdx.push_back(idx[i]);
      topMaxPixCoherence.push_back(xCorrValue[i]);
   }
   }

   void recoData::setMaxPixInfoEachLayer(recoSettings *settings, int *idx, float *xCorrValue){

   //if( onion->nLayer != (int)idx_size) cerr<<"Warning!! maxPixIdxEachLayer size not correct\n";
   maxPixIdxEachLayer.clear();
   maxPixCoherenceEachLayer.clear();
   for(int i=0; i<settings->nLayer; i++){
      maxPixIdxEachLayer.push_back(idx[i]);
      maxPixCoherenceEachLayer.push_back(xCorrValue[i]);
   }
   }

   void recoData::setLikelihoodAndPValue(double _likelihood, double _pValue){

   likelihood = _likelihood;
   pValue = _pValue;

   }

   void recoData::setInWindowSNR(float _inWindowSNR){

   inWindowSNR = _inWindowSNR;

   }

   void recoData::setUnmodSNR(float _unmodSNR){

   unmodSNR = _unmodSNR;

   }

   void recoData::setFlag(int _flag){

   flag = _flag;

   }

   void recoData::setEventTrigType(int _eventTrigType){

   eventTrigType = _eventTrigType;

   }

   void recoData::setEventId(int _eventId){

   eventId = _eventId;

   }

   void recoData::setEventNumber(int _eventNumber){

   eventNumber = _eventNumber;

   }

   void recoData::setEventTime(int _unixTime, int _unixTimeUs, int _timeStamp){

    unixTime = _unixTime;
    unixTimeUs = _unixTimeUs;
    timeStamp = _timeStamp;

   }

   void recoData::duplicate(recoSettings *settings, recoData *old){

   int *_topMaxPixIdx         = (int*)calloc(old->topN, sizeof(int));
   float *_topMaxPixCoherence = (float*)calloc(old->topN, sizeof(float));
   for(int i=0; i<old->topN; i++){
      _topMaxPixIdx[i]       = old->topMaxPixIdx[i];
      _topMaxPixCoherence[i] = old->topMaxPixCoherence[i];
   }

   int *_maxPixIdxEachLayer       = (int*)calloc(/*old->onion*/settings->nLayer, sizeof(int));
   float *_maxPixCoherenceEachLayer = (float*)calloc(/*old->onion*/settings->nLayer, sizeof(float));
   for(int i=0; i</*old->onion*/settings->nLayer; i++){
      _maxPixIdxEachLayer[i]       = old->maxPixIdxEachLayer[i];
      _maxPixCoherenceEachLayer[i] = old->maxPixCoherenceEachLayer[i];
   }

   setAllData(old->eventId,     old->eventNumber
             , old->unixTime, old->unixTimeUs, old->timeStamp
             , old->weight
             , old->eventTrigType
             , old->trueZen,    old->trueAzi
             , old->recoZen,    old->recoAzi
             , old->trueRadius, old->recoRadius
             , old->trueRecAngle, old->trueLauAngle
             , old->recoRecAngle, old->recoLauAngle
             , old->recoChan
             , old->maxPixIdx,  old->maxPixCoherence
             //, old->onion
             , settings
             , old->topN
             , _topMaxPixIdx, _topMaxPixCoherence
             , _maxPixIdxEachLayer, _maxPixCoherenceEachLayer
             , old->likelihood, old->pValue
             , old->inWindowSNR, old->unmodSNR
             , old->flag);

/*
   weight = old->weight;

   trueZen = old->trueZen;
   trueAzi = old->trueAzi;
   recoZen = old->recoZen;
   recoAzi = old->recoAzi;
   trueRadius = old->trueRadius;
   recoRadius = old->recoRadius;

   for(int i=0; i<16; i++) recoChan[i] = old->recoChan[i];

   maxPixIdx = old->maxPixIdx;
   maxPixCoherence = old->maxPixCoherence;

   onion = old->onion;

   topN = old->topN;

   topMaxPixIdx.clear();
   topMaxPixCoherence.clear();
   for(int i=0; i<topN; i++){
      topMaxPixIdx.push_back( _topMaxPixIdx[i] );
      topMaxPixCoherence.push_back( _topMaxPixCoherence[i] );
   }

   maxPixIdxEachLayer.clear();
   maxPixCoherenceEachLayer.clear();
   for(int i=0; i<onion->nLayer; i++){
      maxPixIdxEachLayer.push_back( _maxPixIdxEachLayer[i] );
      maxPixCoherenceEachLayer.push_back( _maxPixCoherenceEachLayer[i] );
   }
*/
   free(_topMaxPixIdx);
   free(_topMaxPixCoherence);
   free(_maxPixIdxEachLayer);
   free(_maxPixCoherenceEachLayer);

   }

   void recoData::clear(){

   weight = 0.;
   trueZen = trueAzi = recoZen = recoAzi = 0.f;
   trueRadius = recoRadius = 0.f;
   std::fill(&recoRecAngle[0], &recoRecAngle[16], -1.);
   std::fill(&recoLauAngle[0], &recoLauAngle[16], -1.);
   std::fill(&trueRecAngle[0], &trueRecAngle[16], -1.);
   std::fill(&trueLauAngle[0], &trueLauAngle[16], -1.);
   std::fill(&recoChan[0], &recoChan[16], 0);
   maxPixIdx = 0;
   maxPixCoherence = 0.f;
   //onion = NULL;
   topN = 0;
   topMaxPixIdx.clear();
   topMaxPixCoherence.clear();
   maxPixIdxEachLayer.clear();
   maxPixCoherenceEachLayer.clear();
   likelihood = 0.;
   pValue = 0.;
   inWindowSNR = 0.f;
   unmodSNR = 0.f;
   flag = 0;
   eventTrigType = 0;
   eventId = 0;
   eventNumber = 0;
   unixTime = 0;
   unixTimeUs = 0;
   timeStamp = 0;

   }
