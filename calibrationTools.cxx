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


#include "TChain.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TTree.h"

#include "calibrationTools.h"
#include "AraGeomTool.h"

using namespace std;

//From T. Meures

//********************************************************************************************************//
//*** A small module to vertically invert waveforms. *****************************************************//
//*** Some channels in ARA03 seem to have inverted amplifier outputs and need to be corrected by this. ***//
//********************************************************************************************************//
/*
void invertGraph(TGraph *gr)
{
        double t1, v1;
        for(int i=0;i<gr->GetN();i++)
        {
                gr->GetPoint(i,t1,v1);
                gr->SetPoint(i,t1,-v1);
        }
}
*/
int calibrateGeometryAndDelays(const RawAraStationEvent *rawEvPtr,
                                double (&posDelayArray)[4][4], double *pulserCorr,
                                const float stationCenterDepth,
                                vector<vector<double> >& ant_loc,
                                vector<vector<double> >& pul_loc){

//From T Meures
//**************************************************************************************************//
////*****here the station correction file is read. the rows are the different strings,****************//
////*****the 4 coloms are the three coordinate-corrections and then the delay-correction.*************//
////**************************************************************************************************//
////*** The structure of the file is:
////***   4 X-corr
////***   4 Y-corr
////***   4 Z-corr
////***   4 delay-corr
////***   5 calpulser corrections:
////***      NON-reference(So far this is always D5):
////***         Z-corr
////***         X-corr
////***         Y-corr
////***      REFERENCE-pulser:
////***         Z-corr
////***         distance-corr ( to be applied as: (1 + distance-corr)*x-coord. or y-coord.  )
////***   1 light speed correction
////***   4 slack corrections.
////*** The last five are not valid!
        ifstream ind;
        char posDelayFile[200];
        int stationId = rawEvPtr->stationId;
        sprintf(posDelayFile, "geometryResultsARA%dE.txt", rawEvPtr->stationId );
        //cout<<"stationId: "<<stationId<<" posDelayFile: "<<posDelayFile<<endl;
        ind.open(posDelayFile);
        //double posDelayArray[4][4]={{0}};
        //double pulserCorr[5] = {0};
        double corrections;
        if(ind.good()){
          for(int i=0;i<4;i++){
            for(int j=0; j<4;j++){
                ind >> corrections;
                posDelayArray[j][i] = corrections;
            }
          }
          for(int i=0;i<5;i++){
                ind >> corrections;
                pulserCorr[i] = corrections;
          }
        }
        else {
        cout << "couldn't read position and delay correction file!!" << endl;
        return -1;
        }
        ind.close();

//***Here the station coordinates are read. the geometry correction is included as well as ***//
//***two specific reading errors, which were not corrected in araroot. ***********************//
//***also is the station center set to 180m under the ice. ***********************************//
//
        AraGeomTool *geom = AraGeomTool::Instance();
        std::vector<double> antl;
        //std::vector<std::vector<double> > ant_loc;
        Double_t *antloc=0;

/*
 * Calibrating antenna positions
 */

        for(int a=0;a<16;a++){
                antloc = geom->getStationInfo(rawEvPtr->stationId)->getAntennaInfo(a)->getLocationXYZ();
                cout << "The location is then: " << antloc[0] << "  " << antloc[1] << "   " << antloc[2] << endl;
                if(stationId==2 && a==0)  antloc[2] = antloc[2] + 1.68;
                if(stationId==3 && a==10) antloc[2] = antloc[2] + 2.01;

                antl.push_back(antloc[0] + posDelayArray[a%4][0]);
                antl.push_back(antloc[1] + posDelayArray[a%4][1]);
                //if((a/4)%2==1)antl.push_back(antloc[2]+180.0 + posDelayArray[a%4][2] + slackArray[a%4]);
                //else
                antl.push_back(antloc[2]+ stationCenterDepth + posDelayArray[a%4][2]);

                cout << "Correcting: " << posDelayArray[a%4][0] << "  " << posDelayArray[a%4][1] << "  " << posDelayArray[a%4][2] << //"  " << slackArray[a%4] <<
                endl;
                cout << "antDepth[1]["<<a<<"] = " << antl[2] << endl;
                printf("Corrected Rx %d X Y Z: %f %f %f\n", a, antl[0], antl[1], antl[2]);

                ant_loc.push_back(antl);
                antl.clear();
         }

/*
 * Calibrating pulser positions
 */
         cout<<"Number of of calpulsers for station "<<rawEvPtr->stationId<<": "
             <<geom->getStationInfo(rawEvPtr->stationId)->getNumCalAnts()<<endl;

         for(int c=0; c<geom->getStationInfo(rawEvPtr->stationId)->getNumCalAnts(); c++){

                 antloc = geom->getStationInfo(rawEvPtr->stationId)->getCalAntennaInfo(c)->getLocationXYZ();
                 string locName(&geom->getStationInfo(rawEvPtr->stationId)->getCalAntennaInfo(c)->locationName[0]);
                 cout<<"locName: "<<locName<<endl;
                 if( locName[locName.length()-1]  == '5' ){
                 antl.push_back( antloc[0] + pulserCorr[1] );
                 antl.push_back( antloc[1] + pulserCorr[2] );
                 antl.push_back( antloc[2] + stationCenterDepth + pulserCorr[0] );
                 } else if ( locName[locName.length()-1] == '6'){
                 antl.push_back( (1.+pulserCorr[4]) * antloc[0] );
                 antl.push_back( (1.+pulserCorr[4]) * antloc[1] );
                 antl.push_back( antloc[2] + stationCenterDepth + pulserCorr[3] );
                 } else {
                 cerr<<"Pulser name undefined\n";
                 //return -1;
                 }
                 printf("Corrected pulser %d X: %f Y: %f Z: %f\n", c, antl[0], antl[1], antl[2]);
                 pul_loc.push_back(antl);
                 antl.clear();
         }
return 0;
}

int calibrateGeometryAndDelaysPlusAdHocDepthOffSet(const RawAraStationEvent *rawEvPtr,
                                double (&posDelayArray)[4][4], double *pulserCorr,
                                const float stationCenterDepth,
                                vector<vector<double> >& ant_loc,
                                vector<vector<double> >& pul_loc,
                                float adhocOffset){

//From T Meures
//**************************************************************************************************//
////*****here the station correction file is read. the rows are the different strings,****************//
////*****the 4 coloms are the three coordinate-corrections and then the delay-correction.*************//
////**************************************************************************************************//
////*** The structure of the file is:
////***   4 X-corr
////***   4 Y-corr
////***   4 Z-corr
////***   4 delay-corr
////***   5 calpulser corrections:
////***      NON-reference(So far this is always D5):
////***         Z-corr
////***         X-corr
////***         Y-corr
////***      REFERENCE-pulser:
////***         Z-corr
////***         distance-corr ( to be applied as: (1 + distance-corr)*x-coord. or y-coord.  )
////***   1 light speed correction
////***   4 slack corrections.
////*** The last five are not valid!
        ifstream ind;
        char posDelayFile[200];
        int stationId = rawEvPtr->stationId;
        sprintf(posDelayFile, "geometryResultsARA%dE.txt", rawEvPtr->stationId );
        //cout<<"stationId: "<<stationId<<" posDelayFile: "<<posDelayFile<<endl;
        ind.open(posDelayFile);
        //double posDelayArray[4][4]={{0}};
        //double pulserCorr[5] = {0};
        double corrections;
        if(ind.good()){
          for(int i=0;i<4;i++){
            for(int j=0; j<4;j++){
                ind >> corrections;
                posDelayArray[j][i] = corrections;
            }
          }
          for(int i=0;i<5;i++){
                ind >> corrections;
                pulserCorr[i] = corrections;
          }
        }
        else {
        cout << "couldn't read position and delay correction file!!" << endl;
        return -1;
        }
        ind.close();

//***Here the station coordinates are read. the geometry correction is included as well as ***//
//***two specific reading errors, which were not corrected in araroot. ***********************//
//***also is the station center set to 180m under the ice. ***********************************//
//
        AraGeomTool *geom = AraGeomTool::Instance();
        std::vector<double> antl;
        //std::vector<std::vector<double> > ant_loc;
        Double_t *antloc=0;

/*
 * Calibrating antenna positions
 */

        for(int a=0;a<16;a++){
                antloc = geom->getStationInfo(rawEvPtr->stationId)->getAntennaInfo(a)->getLocationXYZ();
                cout << "The location is then: " << antloc[0] << "  " << antloc[1] << "   " << antloc[2] << endl;
                if(stationId==2 && a==0)  antloc[2] = antloc[2] + 1.68;
                if(stationId==3 && a==10) antloc[2] = antloc[2] + 2.01;

                antl.push_back(antloc[0] + posDelayArray[a%4][0]);
                antl.push_back(antloc[1] + posDelayArray[a%4][1]);
                //if((a/4)%2==1)antl.push_back(antloc[2]+180.0 + posDelayArray[a%4][2] + slackArray[a%4]);
                //else
                antl.push_back(antloc[2]+ stationCenterDepth + posDelayArray[a%4][2] + adhocOffset);

                cout << "Correcting: " << posDelayArray[a%4][0] << "  " << posDelayArray[a%4][1] << "  " << posDelayArray[a%4][2] << //"  " << slackArray[a%4] <<
                endl;
                cout << "antDepth[1]["<<a<<"] = " << antl[2] << endl;
                printf("Corrected Rx %d X Y Z: %f %f %f\n", a, antl[0], antl[1], antl[2]);

                ant_loc.push_back(antl);
                antl.clear();
         }

/*
 * Calibrating pulser positions
 */
         cout<<"Number of of calpulsers for station "<<rawEvPtr->stationId<<": "
             <<geom->getStationInfo(rawEvPtr->stationId)->getNumCalAnts()<<endl;

         for(int c=0; c<geom->getStationInfo(rawEvPtr->stationId)->getNumCalAnts(); c++){

                 antloc = geom->getStationInfo(rawEvPtr->stationId)->getCalAntennaInfo(c)->getLocationXYZ();
                 string locName(&geom->getStationInfo(rawEvPtr->stationId)->getCalAntennaInfo(c)->locationName[0]);
                 cout<<"locName: "<<locName<<endl;
                 if( locName[locName.length()-1]  == '5' ){
                 antl.push_back( antloc[0] + pulserCorr[1] );
                 antl.push_back( antloc[1] + pulserCorr[2] );
                 antl.push_back( antloc[2] + stationCenterDepth + pulserCorr[0] );
                 } else if ( locName[locName.length()-1] == '6'){
                 antl.push_back( (1.+pulserCorr[4]) * antloc[0] );
                 antl.push_back( (1.+pulserCorr[4]) * antloc[1] );
                 antl.push_back( antloc[2] + stationCenterDepth + pulserCorr[3] );
                 } else {
                 cerr<<"Pulser name undefined\n";
                 //return -1;
                 }
                 printf("Corrected pulser %d X: %f Y: %f Z: %f\n", c, antl[0], antl[1], antl[2]);
                 pul_loc.push_back(antl);
                 antl.clear();
         }
return 0;
}

int getAraSimStationGeometry(vector<vector<double> >& ant_loc){

double zCenter = -180.;
vector<double> xyz;

/* TVpols */
xyz.push_back(-7.358517);xyz.push_back(6.783020);xyz.push_back(-182.999704-zCenter);
ant_loc.push_back(xyz);
xyz.clear();
xyz.push_back(6.783009);xyz.push_back(6.783009);xyz.push_back(-183.021940-zCenter);
ant_loc.push_back(xyz);
xyz.clear();
xyz.push_back(-7.358505);xyz.push_back(-7.358505);xyz.push_back(-182.977468-zCenter);
ant_loc.push_back(xyz);
xyz.clear();
xyz.push_back(6.783020);xyz.push_back(-7.358517);xyz.push_back(-182.999704-zCenter);
ant_loc.push_back(xyz);
xyz.clear();
/* BVpols */
xyz.push_back(-7.385228);xyz.push_back(6.756271);xyz.push_back(-199.999662-zCenter);
ant_loc.push_back(xyz);
xyz.clear();
xyz.push_back(6.756259);xyz.push_back(6.756259);xyz.push_back(-200.021898-zCenter);
ant_loc.push_back(xyz);
xyz.clear();
xyz.push_back(-7.385217);xyz.push_back(-7.385217);xyz.push_back(-199.977426-zCenter);
ant_loc.push_back(xyz);
xyz.clear();
xyz.push_back(6.756271);xyz.push_back(-7.385228);xyz.push_back(-199.999662-zCenter);
ant_loc.push_back(xyz);
xyz.clear();
/* BHpols */
xyz.push_back(-7.382086);xyz.push_back(6.759418);xyz.push_back(-197.999667-zCenter);
ant_loc.push_back(xyz);
xyz.clear();
xyz.push_back(6.759406);xyz.push_back(6.759406);xyz.push_back(-198.021903-zCenter);
ant_loc.push_back(xyz);
xyz.clear();
xyz.push_back(-7.382074);xyz.push_back(-7.382074);xyz.push_back(-197.977431-zCenter);
ant_loc.push_back(xyz);
xyz.clear();
xyz.push_back(6.759418);xyz.push_back(-7.382086);xyz.push_back(-197.999667-zCenter);
ant_loc.push_back(xyz);
xyz.clear();
/* THpols */
xyz.push_back(-7.355374);xyz.push_back(6.786167);xyz.push_back(-180.999709-zCenter);
ant_loc.push_back(xyz);
xyz.clear();
xyz.push_back(6.786156);xyz.push_back(6.786156);xyz.push_back(-181.021945-zCenter);
ant_loc.push_back(xyz);
xyz.clear();
xyz.push_back(-7.355362);xyz.push_back(-7.35536);xyz.push_back(-180.977473-zCenter);
ant_loc.push_back(xyz);
xyz.clear();
xyz.push_back(6.786167);xyz.push_back(-7.355374);xyz.push_back(-180.999709-zCenter);
ant_loc.push_back(xyz);
xyz.clear();


return 0;
}

int getAraSimStationGeometry(vector<vector<double> >& ant_loc, Detector *detector, Settings *settings){

double zCenter = -180.;
vector<double> xyz;
double stationX = detector->stations[0].GetX();
double stationY = detector->stations[0].GetY();
double stationZ = detector->stations[0].GetZ();
int string_i, antenna_i, AraRootChannel;

for(int AraRootChan=0; AraRootChan<16; AraRootChan++){
   for(int AraSimChan=0; AraSimChan<16; AraSimChan++){

   string_i = detector->getStringfromArbAntID(0, AraSimChan);
   antenna_i= detector->getAntennafromArbAntID(0,AraSimChan);
   AraRootChannel = detector->GetChannelfromStringAntenna(0, string_i, antenna_i, settings) - 1;

   //cout<<"AraSimChan: "<<AraSimChan<<" string_i: "<<string_i<<" antenna_i: "<<antenna_i<<" AraRootChannel: "<<AraRootChannel<<endl;

   if( AraRootChannel == AraRootChan ){
   xyz.push_back(detector->stations[0].strings[string_i].antennas[antenna_i].GetX() - stationX);
   xyz.push_back(detector->stations[0].strings[string_i].antennas[antenna_i].GetY() - stationY);
   xyz.push_back(detector->stations[0].strings[string_i].antennas[antenna_i].GetZ() - stationZ - zCenter);
   cout<<"x: "<<xyz[0]<<" y: "<<xyz[1]<<" z: "<<xyz[2]<<endl;
   ant_loc.push_back(xyz);
   xyz.clear();
   }

   }
   }

return 0;
}

int getSeckelStationGeometry(vector<vector<double> >& ant_loc){

   vector<double> antl;
   double stationCenter[3] = {4001.59, -2595.01, -184.267};

   double xyz[3];

   ant_loc.clear();

   /* TVpols */

   xyz[0] = 3993.34; xyz[1] = -2589.51; xyz[2] =  -177.654;
   for(int i=0; i<3; i++) xyz[i] -= stationCenter[i];
   antl.assign(xyz,xyz+3);
   ant_loc.push_back(antl);

   antl.clear();
   xyz[0] = 4006.98; xyz[1] =  -2586.51; xyz[2] = -176.17;
   for(int i=0; i<3; i++) xyz[i] -= stationCenter[i];
   antl.assign(xyz,xyz+3);
   ant_loc.push_back(antl);

   antl.clear();
   xyz[0] = 3996.2; xyz[1] = -2603.94; xyz[2] =  -176.12;
   for(int i=0; i<3; i++) xyz[i] -= stationCenter[i];
   antl.assign(xyz,xyz+3);
   ant_loc.push_back(antl);

   antl.clear();
   xyz[0] = 4009.75; xyz[1] = -2600.1; xyz[2] = -176.359;
   for(int i=0; i<3; i++) xyz[i] -= stationCenter[i];
   antl.assign(xyz,xyz+3);
   ant_loc.push_back(antl);

   /* BVpols */

   antl.clear();
   xyz[0] = 3993.4; xyz[1] = -2589.49; xyz[2] = -195.229;
   for(int i=0; i<3; i++) xyz[i] -= stationCenter[i];
   antl.assign(xyz,xyz+3);
   ant_loc.push_back(antl);

   antl.clear();
   xyz[0] = 4007.04; xyz[1] = -2586.49; xyz[2] = -195.223;
   for(int i=0; i<3; i++) xyz[i] -= stationCenter[i];
   antl.assign(xyz,xyz+3);
   ant_loc.push_back(antl);

   antl.clear();
   xyz[0] = 3996.26; xyz[1] = -2603.92; xyz[2] = -195.173;
   for(int i=0; i<3; i++) xyz[i] -= stationCenter[i];
   antl.assign(xyz,xyz+3);
   ant_loc.push_back(antl);

   antl.clear();
   xyz[0] = 4009.81; xyz[1] = -2600.08; xyz[2] = -195.248;
   for(int i=0; i<3; i++) xyz[i] -= stationCenter[i];
   antl.assign(xyz,xyz+3);
   ant_loc.push_back(antl);

   /* THpols */

   antl.clear();
   xyz[0] = 3993.35; xyz[1] = -2589.5; xyz[2] = -173.219;
   for(int i=0; i<3; i++) xyz[i] -= stationCenter[i];
   antl.assign(xyz,xyz+3);
   ant_loc.push_back(antl);

   antl.clear();
   xyz[0] = 4006.97; xyz[1] = -2586.51; xyz[2] = -173.251;
   for(int i=0; i<3; i++) xyz[i] -= stationCenter[i];
   antl.assign(xyz,xyz+3);
   ant_loc.push_back(antl);

   antl.clear();
   xyz[0] = 3996.19; xyz[1] = -2603.94; xyz[2] = -172.999;
   for(int i=0; i<3; i++) xyz[i] -= stationCenter[i];
   antl.assign(xyz,xyz+3);
   ant_loc.push_back(antl);

   antl.clear();
   xyz[0] = 4009.74; xyz[1] = -2600.11; xyz[2] = -173.402;
   for(int i=0; i<3; i++) xyz[i] -= stationCenter[i];
   antl.assign(xyz,xyz+3);
   ant_loc.push_back(antl);

   /* BHpols */

   antl.clear();
   xyz[0] = 3993.39; xyz[1] = -2589.49; xyz[2] = -192.272;
   for(int i=0; i<3; i++) xyz[i] -= stationCenter[i];
   antl.assign(xyz,xyz+3);
   ant_loc.push_back(antl);

   antl.clear();
   xyz[0] = 4007.03; xyz[1] = -2586.49; xyz[2] = -191.938;
   for(int i=0; i<3; i++) xyz[i] -= stationCenter[i];
   antl.assign(xyz,xyz+3);
   ant_loc.push_back(antl);

   antl.clear();
   xyz[0] = 3996.25; xyz[1] = -2603.92; xyz[2] = -192.052;
   for(int i=0; i<3; i++) xyz[i] -= stationCenter[i];
   antl.assign(xyz,xyz+3);
   ant_loc.push_back(antl);

   antl.clear();
   xyz[0] = 4009.8; xyz[1] = -2600.08; xyz[2] = -191.963;
   for(int i=0; i<3; i++) xyz[i] -= stationCenter[i];
   antl.assign(xyz,xyz+3);
   ant_loc.push_back(antl);

return 0;
}
