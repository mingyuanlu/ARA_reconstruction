#include "analysisTools.h"

ClassImp(recoSettings);
ClassImp(recoData);
ClassImp(cutParameter);
ClassImp(ARA02_cutValues);

using namespace std;

int getRunType(string STATION, int runNum){

   int type;

   if ( STATION == "ARA02" ) {
      if ( runNum >= 2275 && runNum <= 3463) {
         type = 1;
      } else if ( runNum >= 120 && runNum <= 2274 ) {
         type = 2;
      } else if ( runNum >= 3465 && runNum <= 4027 ){
         type = 3;
      } else if ( ( runNum >= 4029 && runNum <= 6481) || ( runNum >= 3 && runNum <= 60) || ( runNum >= 8100 && runNum <= 8246) ){
         type = 4;
      } else if ( (runNum >= 6500 && runNum <= 8097) || runNum == 0 ){
         type = 5;
      } else {
         type = 0;
      }
   } else if ( STATION == "ARA03" ) {
      if ( runNum >= 1449 && runNum <= 3061 ) {
         type = 1;
      } else if ( runNum >= 470 && runNum <= 1448 ) {
         type = 2;
      } else if ( (runNum >= 3063 && runNum <=6004) || ( runNum >= 7658 && runNum <= 7808 ) ){
         type = 3;
      } else if ( runNum >= 6005 && runNum <= 7653 ){
         type = 4;
      } else {
         type = 0;
      }
   } else {
      type = 0;
   }

   return type;

}

cutParameter::cutParameter() {  }
cutParameter::~cutParameter() { /*default destructor*/ }

ARA02_cutValues::ARA02_cutValues(){ initialize(); }

void ARA02_cutValues::initialize(){


   setValue(cwImpCut[0], 0.29936375,0.00071682,-0.00071682);
   setValue(cwImpCut[1], 0.32188931,0.00097344,-0.00097344);
   setValue(cwImpCut[2], 0.2955359,0.00070074,-0.00070074);
   setValue(cwImpCut[3], 0.30005502,0.00111114,-0.00111114);
   setValue(cwImpCut[4], 0.29206096,0.0019572,-0.0019572);

   setValue(snrCut_inBand[0], 8.97, 0.018104, -0.017988);
   setValue(snrCut_inBand[1], 9.116000, 0.017207, -0.017098);
   setValue(snrCut_inBand[2], 8.640000, 0.093084, -0.089897);
   setValue(snrCut_inBand[3], 9.360000, 0.022266, -0.022081);
   setValue(snrCut_inBand[4], 9.480000, 0.033537, -0.033128);

   setValue(coherenceCut_inBand[0], 0.108971, 0.000276, -0.000271);
   setValue(coherenceCut_inBand[1], 0.111546, 0.000198, -0.000195);
   setValue(coherenceCut_inBand[2], 0.106706, 0.000192, -0.000189);
   setValue(coherenceCut_inBand[3], 0.099455, 0.000143, -0.000141);
   setValue(coherenceCut_inBand[4], 0.098104, 0.000146, -0.000144);

   setValue(snrCut_outOfBand[0],9.000000, 0.033720, -0.033304);
   setValue(snrCut_outOfBand[1],9.128000, 0.021496, -0.021326);
   setValue(snrCut_outOfBand[2],8.704000, 0.096490, -0.093141);
   setValue(snrCut_outOfBand[3],9.524000, 0.033422, -0.033023);
   setValue(snrCut_outOfBand[4],10.028000, 0.184248, -0.174582);

   setValue(coherenceCut_outOfBand[0],0.103402, 0.000165, -0.000162);
   setValue(coherenceCut_outOfBand[1],0.105563, 0.000149, -0.000147);
   setValue(coherenceCut_outOfBand[2],0.101873, 0.000134, -0.000132);
   setValue(coherenceCut_outOfBand[3],0.094399, 0.000217, -0.000212);
   setValue(coherenceCut_outOfBand[4],0.092622, 0.000142, -0.000140);

   setValue(impCut, 0.25563478, 0.00124038, -0.00124038);

   setValue(zenMin[0],0.5409,	-0.001924625,	0.001924625); // +: cut region is larger -: cut regin is smaller
   setValue(zenMin[1],-31.4,	-0.030723404,	0.030723404);
   setValue(zenMin[2],-27.78,	-0.0405,	0.0405);
   setValue(zenMin[3],24.07,	-0.005442857,	0.005442857);

   setValue(zenMax[0],10.99,	0.01,	-0.01);
   setValue(zenMax[1],-14.26,	0.01,	-0.01);
   setValue(zenMax[2],-14.22,	0.0285,	-0.0285);
   setValue(zenMax[3],31.45,	0.050333333,	-0.050333333);

   setValue(aziMin[0],54.96,	-0.01,	0.01);
   setValue(aziMin[1],324.05,	-0.003317073,	0.003317073);
   setValue(aziMin[2],50.8,	-0.041156863,	0.041156863);
   setValue(aziMin[3],67.3,	-0.0022315,	0.0022315);

   setValue(aziMax[0],66.74,	0.0001,	-0.0001);
   setValue(aziMax[1],341.01,	0.027090909,	-0.027090909);
   setValue(aziMax[2],58.32,	0.015909091,	-0.015909091);
   setValue(aziMax[3],75.25,	0.010268085,	-0.010268085);

   setValue(surfaceCut_constantN, 35.648,	-0.693,	0.693); // +: cut region is larger, -: cut region is smaller
   setValue(surfaceCut_iterReco, 36.77852,	-0.588241,	0.588241);

}

ARA02_cutValues::~ARA02_cutValues(){ /* default destructor */ }

void ARA02_cutValues::setValue(cutParameter& param, double _val, double _plus, double _minus){

   param.val = _val;
   param.plus = _plus;
   param.minus = _minus;

}
