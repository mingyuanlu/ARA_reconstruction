#include "analysisTools.h"

ClassImp(recoSettings);
ClassImp(recoData);

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
