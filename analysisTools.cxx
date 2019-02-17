#include "analysisTools.h"
/*
ClassImp(recoSettings);
ClassImp(recoData);
ClassImp(cutParameter);
ClassImp(ARA02_cutValues);
*/
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
ARA02_cutValues::~ARA02_cutValues(){ /* default destructor */ }

void ARA02_cutValues::initialize(){


   setValue(cwImpCut[0], 0.29936375,0.00071682,-0.00071682);
   setValue(cwImpCut[1], 0.32188931,0.00097344,-0.00097344);
   setValue(cwImpCut[2], 0.2955359,0.00070074,-0.00070074);
   setValue(cwImpCut[3], 0.30005502,0.00111114,-0.00111114);
   setValue(cwImpCut[4], 0.29206096,0.0019572,-0.0019572);

//Cut values from crest factor n-chan
/*
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
*/

   //Cut values from sliding V^2 SNR n-chan
   setValue(snrCut_inBand[0], 12.935000, 0.552648, -0.509019);
   setValue(snrCut_inBand[1], 13.365000, 0.051088, -0.050650);
   setValue(snrCut_inBand[2], 12.230000, 0.053964, -0.053508);
   setValue(snrCut_inBand[3], 12.535000, 0.039439, -0.039147);
   setValue(snrCut_inBand[4], 12.725000, 0.041504, -0.041191);

   setValue(coherenceCut_inBand[0], 0.108061, 0.000236, -0.000231);
   setValue(coherenceCut_inBand[1], 0.109367, 0.000201, -0.000198);
   setValue(coherenceCut_inBand[2], 0.108527, 0.000255, -0.000250);
   setValue(coherenceCut_inBand[3], 0.099249, 0.000142, -0.000140);
   setValue(coherenceCut_inBand[4], 0.097097, 0.000132, -0.000131);

   setValue(snrCut_outOfBand[0],11.815000, 0.048386, -0.047920);
   setValue(snrCut_outOfBand[1],13.325000, 0.042806, -0.042496);
   setValue(snrCut_outOfBand[2],12.185000, 0.057194, -0.056679);
   setValue(snrCut_outOfBand[3],13.180000, 0.052173, -0.051708);
   setValue(snrCut_outOfBand[4],12.860000, 0.041406, -0.041101);

   setValue(coherenceCut_outOfBand[0],0.103742, 0.000151, -0.000149);
   setValue(coherenceCut_outOfBand[1],0.107307, 0.000144, -0.000142);
   setValue(coherenceCut_outOfBand[2],0.103737, 0.000132, -0.000130);
   setValue(coherenceCut_outOfBand[3],0.094770, 0.000120, -0.000119);
   setValue(coherenceCut_outOfBand[4],0.095079, 0.000324, -0.000315);

   setValue(impCut, 0.25563478, 0.00124038, -0.00124038);
/*
   setValue(zenMin[0],0.5409,	-0.001924625,	0.001924625); // plus: cut region is larger minus: cut regin is smaller
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
*/

   setValue(zenMin[0], -32.445122,  -0.279150, 0.279150); // plus: cut region is larger minus: cut regin is smaller
   setValue(zenMin[1], -0.5000384756620608, -0.049698, 0.049698);
   setValue(zenMin[2], 26.676257582375015,  -0.014895, 0.014895);

   setValue(zenMax[0], -15.09135433852994, 0.644769, -0.644769);
   setValue(zenMax[1], 11.197893462051383, 0.005116, -0.005116);
   setValue(zenMax[2], 36.74904979886475, 0.002851, -0.002851);

   setValue(aziMin[0], 329.1813651101607, -0.011490, 0.011490);
   setValue(aziMin[1], 55.31723520228566, -0.000069, 0.000069);
   setValue(aziMin[2], 330.64429515858507,  -0.009783, 0.009783);

   setValue(aziMax[0], 341.4794531938459, 0.287608, -0.287608);
   setValue(aziMax[1], 67.22731842523841, 0.018698, -0.018698);
   setValue(aziMax[2], 340.74104652018354, 0.739811, -0.739811);

   setValue(surfaceCut_constantN, 35.648,	-0.693,	0.693); // plus: cut region is larger, minus: cut region is smaller
   setValue(surfaceCut_iterReco, 36.77852,	-0.588241,	0.588241);

}

void ARA02_cutValues::setValue(cutParameter& param, double _val, double _plus, double _minus){

   param.val = _val;
   param.plus = _plus;
   param.minus = _minus;

}
