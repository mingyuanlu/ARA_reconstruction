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

   treg = NULL;

   constantNMaxPixIdx = 0;
   constantNMaxPixCoherence = 0.f;
   constantNZen = constantNAzi = 0.f;

   inWindowSNR_V = inWindowSNR_H = 0.;
   passAnotherPolNchnl = false;

   maxPixIdx2 = 0;
   maxPixCoherence2 = 0.f;
   topMaxPixIdx2.clear();
   topMaxPixCoherence2.clear();
   maxPixIdxEachLayer2.clear();
   maxPixCoherenceEachLayer2.clear();

   interactionProbability = 0.;
   probability = 0.;

   numSatChan = 0;
   std::fill(&satChan[0], &satChan[16], 0);

   std::fill(&channelInWindowSNR[0], &channelInWindowSNR[16], 0.f);

   survivalProbability = 0.;

   }



   void recoData::setAllData(
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
   , bool _passAnotherPolNchnl
   , int _flag
   , trackEngine *_treg
   , int _constantNMaxPixIdx, float _constantNMaxPixCoherence
   , float _constantNZen, float _constantNAzi
   /*, int idx2, float xCorrValue2
   , int *_topMaxPixIdx2, float *_topMaxPixCoherence2
   , int *_maxPixIdxEachLayer2, float *_maxPixCoherenceEachLayer2*/
   )
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

   treg = _treg;

   constantNMaxPixIdx = _constantNMaxPixIdx;
   constantNMaxPixCoherence = _constantNMaxPixCoherence;
   constantNZen = _constantNZen;
   constantNAzi = _constantNAzi;

   inWindowSNR_V = _inWindowSNR_V;
   inWindowSNR_H = _inWindowSNR_H;
   passAnotherPolNchnl = _passAnotherPolNchnl;
/*
   maxPixIdx2 = idx2;
   maxPixCoherence2 = xCorrValue2;

   topMaxPixIdx2.clear();
   topMaxPixCoherence2.clear();
   for(int i=0; i<topN; i++){
      topMaxPixIdx2.push_back( _topMaxPixIdx2[i] );
      topMaxPixCoherence2.push_back( _topMaxPixCoherence2[i] );
   }

   maxPixIdxEachLayer2.clear();
   maxPixCoherenceEachLayer2.clear();
   for(int i=0; i<_settings->nLayer; i++){
      maxPixIdxEachLayer2.push_back( _maxPixIdxEachLayer2[i] );
      maxPixCoherenceEachLayer2.push_back( _maxPixCoherenceEachLayer2[i] );
   }
*/

   //interactionProbability = _interactionProbability;
   //probability = _probability;

   }

   void recoData::setAllData2(recoSettings *_settings
   , int idx2, float xCorrValue2
   , int *_topMaxPixIdx2, float *_topMaxPixCoherence2
   , int *_maxPixIdxEachLayer2, float *_maxPixCoherenceEachLayer2){

   maxPixIdx2 = idx2;
   maxPixCoherence2 = xCorrValue2;

   topMaxPixIdx2.clear();
   topMaxPixCoherence2.clear();
   for(int i=0; i<topN; i++){
      topMaxPixIdx2.push_back( _topMaxPixIdx2[i] );
      topMaxPixCoherence2.push_back( _topMaxPixCoherence2[i] );
   }

   maxPixIdxEachLayer2.clear();
   maxPixCoherenceEachLayer2.clear();
   for(int i=0; i<_settings->nLayer; i++){
      maxPixIdxEachLayer2.push_back( _maxPixIdxEachLayer2[i] );
      maxPixCoherenceEachLayer2.push_back( _maxPixCoherenceEachLayer2[i] );
   }

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

   void recoData::setMaxPix2Info(int idx, float xCorrValue){

   maxPixIdx2 = idx;
   maxPixCoherence2 = xCorrValue;

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

   void recoData::setTopMaxPix2Info(int *idx, float *xCorrValue){ //pass a sorted array of topN max pix indices

   //if( topN != (int)idx_size ) cerr<<"Warning!! topMaxPixIdx size not correct\n";
   topMaxPixIdx2.clear();
   topMaxPixCoherence2.clear();
   for(int i=0; i<topN; i++){
      topMaxPixIdx2.push_back(idx[i]);
      topMaxPixCoherence2.push_back(xCorrValue[i]);
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

   void recoData::setMaxPix2InfoEachLayer(recoSettings *settings, int *idx, float *xCorrValue){

   //if( onion->nLayer != (int)idx_size) cerr<<"Warning!! maxPixIdxEachLayer size not correct\n";
   maxPixIdxEachLayer2.clear();
   maxPixCoherenceEachLayer2.clear();
   for(int i=0; i<settings->nLayer; i++){
      maxPixIdxEachLayer2.push_back(idx[i]);
      maxPixCoherenceEachLayer2.push_back(xCorrValue[i]);
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

   void recoData::setTreg(trackEngine *_treg){

      treg = _treg;

   }

   void recoData::setConstantNMaxPixInfo(int _constantNMaxPixIdx, float _constantNMaxPixCoherence){

      constantNMaxPixIdx = _constantNMaxPixIdx;
      constantNMaxPixCoherence = _constantNMaxPixCoherence;

   }

   void recoData::setConstantNDir(float _constantNZen, float _constantNAzi){

      constantNZen = _constantNZen;
      constantNAzi = _constantNAzi;

   }

   void recoData::setInWindowSNRBothPol(float _inWindowSNR_V, float _inWindowSNR_H){

      inWindowSNR_V = _inWindowSNR_V;
      inWindowSNR_H = _inWindowSNR_H;

   }

   void  recoData::setPassAnotherPolNchnl(bool _passAnotherPolNchnl){

      passAnotherPolNchnl = _passAnotherPolNchnl;

   }

   void recoData::setProbabilities(double p_int, double p){

      interactionProbability = p_int;
      probability = p;

   }

   void recoData::setSaturatedChannels(int _numSatChan, int *_satChan){

      numSatChan = _numSatChan;
      for(int i=0; i<16; i++) satChan[i] = _satChan[i];

   }

   void recoData::setChannelInWindowSNR(float *_channelInWindowSNR){

      for(int i=0; i<16; i++) channelInWindowSNR[i] = _channelInWindowSNR[i];

   }

   void recoData::setSurvivalProbability(double _survivalProbability){

      survivalProbability = _survivalProbability;

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

   int *_topMaxPixIdx2         = (int*)calloc(old->topN, sizeof(int));
   float *_topMaxPixCoherence2 = (float*)calloc(old->topN, sizeof(float));
   for(int i=0; i<old->topN; i++){
      _topMaxPixIdx2[i]       = old->topMaxPixIdx2[i];
      _topMaxPixCoherence2[i] = old->topMaxPixCoherence2[i];
   }

   int *_maxPixIdxEachLayer2       = (int*)calloc(/*old->onion*/settings->nLayer, sizeof(int));
   float *_maxPixCoherenceEachLayer2 = (float*)calloc(/*old->onion*/settings->nLayer, sizeof(float));
   for(int i=0; i</*old->onion*/settings->nLayer; i++){
      _maxPixIdxEachLayer2[i]       = old->maxPixIdxEachLayer2[i];
      _maxPixCoherenceEachLayer2[i] = old->maxPixCoherenceEachLayer2[i];
   }

   setAllData(old->eventId,     old->eventNumber
             , old->unixTime, old->unixTimeUs, old->timeStamp
             , old->weight
             //, old->interactionProbability, old->probability
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
             , old->inWindowSNR, old->inWindowSNR_V, old->inWindowSNR_H, old->unmodSNR
             , old->passAnotherPolNchnl
             , old->flag
             , old->treg
             , old->constantNMaxPixIdx, old->constantNMaxPixCoherence
             , old->constantNZen, old->constantNAzi
             /*, old->maxPixIdx2,  old->maxPixCoherence2
             , _topMaxPixIdx2, _topMaxPixCoherence2
             , _maxPixIdxEachLayer2, _maxPixCoherenceEachLayer2*/
          );
   setAllData2(settings
   , old->maxPixIdx2,  old->maxPixCoherence2
   , _topMaxPixIdx2, _topMaxPixCoherence2
   , _maxPixIdxEachLayer2, _maxPixCoherenceEachLayer2);

   setProbabilities(old->interactionProbability, old->probability);
   setSaturatedChannels(old->numSatChan, old->satChan);
   setChannelInWindowSNR(old->channelInWindowSNR);
   setSurvivalProbability(old->survivalProbability);
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
   free(_topMaxPixIdx2);
   free(_topMaxPixCoherence2);
   free(_maxPixIdxEachLayer2);
   free(_maxPixCoherenceEachLayer2);
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
   treg = NULL;
   constantNMaxPixIdx = 0;
   constantNMaxPixCoherence = 0.f;
   constantNZen = constantNAzi = 0.f;
   inWindowSNR_V = inWindowSNR_H = 0.f;
   passAnotherPolNchnl = false;
   maxPixIdx2 = 0;
   maxPixCoherence2 = 0.f;
   topMaxPixIdx2.clear();
   topMaxPixCoherence2.clear();
   maxPixIdxEachLayer2.clear();
   maxPixCoherenceEachLayer2.clear();
   interactionProbability = 0.;
   probability = 0;
   numSatChan = 0;
   std::fill(&satChan[0], &satChan[16], 0);
   std::fill(&channelInWindowSNR[0], &channelInWindowSNR[16], 0.f);
   survivalProbability = 0.;

   }

   double computeWeight(Settings *settings, Detector *detector, Event *event, IceModel *antarctica, double zCenter, double L_int){

      /* The correct weight is W = Wprop * Wint = Interaction::weight * L_0 / L_int.
         The only quatity we need to compute here is L_0.
      */

      if(zCenter > 0){ cerr<<"zCenter should be <0 ! zCenter="<<zCenter<<endl; return -1; }
      //double zCenter = -180.;
      double ox, oy, oz; // coor. of the origin
      ox = detector->stations[0].GetX();
      oy = detector->stations[0].GetY();
      oz = detector->stations[0].GetZ() + zCenter;

      double posx, posy, posz; //vertex position
      double nnux, nnuy, nnuz; //neutrino direction
      //double t;                //r = r_0 + t * nnu
      double t[4]; //0: surface, 1: bedrock, 2: side solution 1, 3: side solution 2
      int index[4] = {0};

      posx = event->Nu_Interaction[0].posnu.GetX() - ox;
      posy = event->Nu_Interaction[0].posnu.GetY() - oy;
      posz = event->Nu_Interaction[0].posnu.GetZ() - oz;
      nnux = event->Nu_Interaction[0].nnu.GetX();
      nnuy = event->Nu_Interaction[0].nnu.GetY();
      nnuz = event->Nu_Interaction[0].nnu.GetZ();

      double bedrock_z = antarctica->Surface(event->Nu_Interaction[0].posnu) - antarctica->IceThickness(event->Nu_Interaction[0].posnu) - oz;

      double PR = settings->POSNU_RADIUS;

      t[0] = (-zCenter - posz) / nnuz;
      t[1] = (bedrock_z - posz)/ nnuz;

      double a = nnux*nnux + nnuy*nnuy;
      double b = posx*nnux + posy*nnuy;
      double c = posx*posx + posy*posy - PR*PR;
      t[2] = (-b + sqrt(b*b-4.*a*c)) / (2.*a);
      t[3] = (-b - sqrt(b*b-4.*a*c)) / (2.*a);

      for(int i=0; i<4; i++) t[i] = fabs(t[i]); //only care about the absolute distance
      TMath::Sort(4, t, index); //in descending order
      double dt = t[index[2]] - t[index[3]];
      if(dt<0){ cerr<<"Error! dt < 0! \n"; return -1; }

      double L0x, L0y, L0z
      L0x = nnux * dt;
      L0y = nnuy * dt;
      L0z = nnuz * dt;
      double L_0 = sqrt(L0x*L0x + L0y*L0y + L0z*L0z);

      double weight = event->Nu_Interaction[0].weight * L_0 / L_int;

      printf("zCenter: %f ox: %f oy: %f oz: %f\nposx: %f posy: %f posz: %f\nnnux: %f nnuy: %f nnuz: %f\nbedrock_z: %f PR: %f\n
      t0: %f t1: %f t2: %f t3: %f dt: %f\nL0x: %f L0y: %f L0z: %f L_0: %f\nweight: %f\n",
      zCenter, ox, oy, oz,
      posx, posy, posz,
      nnux, nnuy, nnuz,
      bedrock_z, PR,
      t[0], t[1], t[2], t[3], dt,
      L0x, L0y, L0z, L_0,
      weight);

      return weight;
   }
