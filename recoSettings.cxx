#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include "recoSettings.h"

ClassImp(recoSettings);

using namespace std;

recoSettings::recoSettings(){ initialize(); }

recoSettings::~recoSettings(){ /*default destructor*/ }

void recoSettings::initialize(){

  //programFile = "kernel_3D_analysis.c";
  snprintf(programFile, sizeof(programFile), "kernel_3D_analysis.c");
  nSideExp = 2;
  nLayer   = 1;
  dataType = 0;
  snprintf(triggerCode, sizeof(triggerCode), "111");
  layerAllocationMode = 0;
  skymapSearchMode = 0;
  beamformMethod = 0;
  recoVertexingMode = 0;
  getSkymapMode = 0;
  //recoPolType = "vpol";
  snprintf(recoPolType, sizeof(recoPolType), "vpol");

  nSideExpStart = 2;
  nSideExpEnd   = 8;

  nchnlFilter = 0;
  nchnlCut    = 0;
  nchnlThreshold = 0.;
  nchnlThreshold_A1 = 0.;
  nchnlThreshold_A2 = 0.;
  nchnlThreshold_A3 = 0.;

  iceModel = 0;
  topN     = 50;

  layerFirstRadius = 40;
  layerLastRadius = 5000;

  recordMapData = 0;
  //computePValue = 0;
  computeLLHAndPValue = 0;
  //referenceMapFitFile = "";
  //referenceMapFitFunc = "";
  snprintf(referenceMapFitFile, sizeof(referenceMapFitFile), "");
  snprintf(referenceMapFitFunc, sizeof(referenceMapFitFunc), "");

  snprintf(openCLDeviceType, sizeof(openCLDeviceType), "cpu");
  openCLMaxNumberOfCores = 0;

  maxNumberOfReco = -1;

  constantNFilter = 0;
  surfaceCutAngle = 0.f;

  nchnlThreshold_anotherPol = 0.;

  //remark = "Default.";
  snprintf(remark, sizeof(remark), "Default.");

}

bool recoSettings::readRecoSetupFile(string recoSetupFile){

   if(recoSetupFile == ""){ cerr<<"No reco setup file provided!\n"; return false; }
   ifstream sf(recoSetupFile.c_str());
   string line, label;

   if( sf.is_open() ){
      while( sf.good() ){
         getline(sf, line);
         if(line[0] != "/"[0]){
            label = line.substr(0, line.find_first_of("="));

            if(label == "programFile")              snprintf(programFile, sizeof(programFile), line.substr(line.find_first_of("=")+1).c_str() );
                                                    //programFile         = line.substr(line.find_first_of("=")+1);
            else if(label == "nSideExp")            nSideExp            = atoi( line.substr(line.find_first_of("=")+1).c_str() );
            else if(label == "nLayer")              nLayer              = atoi( line.substr(line.find_first_of("=")+1).c_str() );
            else if(label == "dataType")            dataType            = atoi( line.substr(line.find_first_of("=")+1).c_str() );
            else if(label == "triggerCode"){        triggerCode[0]      = line.substr(line.find_first_of("=")+1).at(0);
                                                    triggerCode[1]      = line.substr(line.find_first_of("=")+1).at(1);
                                                    triggerCode[2]      = line.substr(line.find_first_of("=")+1).at(2); }
            else if(label == "layerAllocationMode") layerAllocationMode = atoi( line.substr(line.find_first_of("=")+1).c_str() );
            else if(label == "skymapSearchMode")    skymapSearchMode    = atoi( line.substr(line.find_first_of("=")+1).c_str() );
            else if(label == "beamformMethod")      beamformMethod      = atoi( line.substr(line.find_first_of("=")+1).c_str() );
            else if(label == "recoVertexingMode")   recoVertexingMode   = atoi( line.substr(line.find_first_of("=")+1).c_str() );
            else if(label == "getSkymapMode")       getSkymapMode       = atoi( line.substr(line.find_first_of("=")+1).c_str() );
            else if(label == "recoPolType")         snprintf(recoPolType, sizeof(recoPolType), line.substr(line.find_first_of("=")+1).c_str() );
                                                    //recoPolType         = line.substr(line.find_first_of("=")+1);
            else if(label == "nSideExpStart")       nSideExpStart       = atoi( line.substr(line.find_first_of("=")+1).c_str() );
            else if(label == "nSideExpEnd")         nSideExpEnd         = atoi( line.substr(line.find_first_of("=")+1).c_str() );
            else if(label == "nchnlFilter")         nchnlFilter         = atoi( line.substr(line.find_first_of("=")+1).c_str() );
            else if(label == "nchnlCut")            nchnlCut            = atoi( line.substr(line.find_first_of("=")+1).c_str() );
            else if(label == "nchnlThreshold")      nchnlThreshold      = atof( line.substr(line.find_first_of("=")+1).c_str() );
            else if(label == "nchnlThreshold_A1")   nchnlThreshold_A1   = atof( line.substr(line.find_first_of("=")+1).c_str() );
            else if(label == "nchnlThreshold_A2")   nchnlThreshold_A2   = atof( line.substr(line.find_first_of("=")+1).c_str() );
            else if(label == "nchnlThreshold_A3")   nchnlThreshold_A3   = atof( line.substr(line.find_first_of("=")+1).c_str() );
            else if(label == "iceModel")            iceModel            = atoi( line.substr(line.find_first_of("=")+1).c_str() );
            else if(label == "topN")                topN                = atoi( line.substr(line.find_first_of("=")+1).c_str() );
            else if(label == "layerFirstRadius")    layerFirstRadius    = atof( line.substr(line.find_first_of("=")+1).c_str() );
            else if(label == "layerLastRadius")     layerLastRadius     = atof( line.substr(line.find_first_of("=")+1).c_str() );
            else if(label == "remark")              snprintf(remark, sizeof(remark), line.substr(line.find_first_of("=")+1).c_str() );
                                                    //remark              = line.substr(line.find_first_of("=")+1);
            else if(label == "recordMapData")       recordMapData       = atoi( line.substr(line.find_first_of("=")+1).c_str() );
            //else if(label == "computePValue")       computePValue       = atoi( line.substr(line.find_first_of("=")+1).c_str() );
            else if(label == "computeLLHAndPValue") computeLLHAndPValue = atoi( line.substr(line.find_first_of("=")+1).c_str() );
            else if(label == "referenceMapFitFile") snprintf(referenceMapFitFile, sizeof(referenceMapFitFile), line.substr(line.find_first_of("=")+1).c_str() );
                                                    //referenceMapFitFile = line.substr(line.find_first_of("=")+1);
            else if(label == "referenceMapFitFunc") snprintf(referenceMapFitFunc, sizeof(referenceMapFitFunc), line.substr(line.find_first_of("=")+1).c_str() );
                                                    //referenceMapFitFunc = line.substr(line.find_first_of("=")+1);
            else if(label == "openCLDeviceType")    snprintf(openCLDeviceType, sizeof(openCLDeviceType), line.substr(line.find_first_of("=")+1).c_str());
            else if(label == "openCLMaxNumberOfCores") openCLMaxNumberOfCores = atoi( line.substr(line.find_first_of("=")+1).c_str() );
            else if(label == "maxNumberOfReco") maxNumberOfReco = atoi( line.substr(line.find_first_of("=")+1).c_str() );
            else if(label == "constantNFilter") constantNFilter = atoi( line.substr(line.find_first_of("=")+1).c_str() );
            else if(label == "surfaceCutAngle") surfaceCutAngle = atof( line.substr(line.find_first_of("=")+1).c_str() );
            else if(label == "nchnlThreshold_anotherPol") nchnlThreshold_anotherPol = atof( line.substr(line.find_first_of("=")+1).c_str() );
            else if(label != "") cerr<<"Undefined parameter detected. Label: "<<label<<endl;
         }
      }

       sf.close();

      /* Parameters sanity checks */

      int errCnt = 0;//error count

      ifstream pf(programFile/*.c_str()*/);
      if( pf.is_open() ) pf.close();
      else{                                       cerr<<"Unable to open "<<programFile<<endl; errCnt++; }
      if( nSideExp < 1 ){                         cerr<<"nSideExp: "<<nSideExp<<endl; errCnt++; }
      if( nLayer < 1 ){                           cerr<<"nLayer: "<<nLayer<<endl; errCnt++; }
      if( (unsigned)(dataType-0) > 1){            cerr<<"dataType: "<<dataType<<endl; errCnt++; }
      //if( triggerCode.length() != 3){             cerr<<"triggerCode: "<<triggerCode<<endl; errCnt++; }
      int code; for(int i=0; i<3; i++){ code = triggerCode[i] - '0'; if( (unsigned)(code-0) > 1){ printf("triggerCode: %s\n", triggerCode); errCnt++;}}
      if( (unsigned)(layerAllocationMode-0) > 1){ cerr<<"layerAllocationMode: "<<layerAllocationMode<<endl; errCnt++; }
      if( (unsigned)(skymapSearchMode-0) > 1){    cerr<<"skymapSearchMode: "<<skymapSearchMode<<endl; errCnt++; }
      if( (unsigned)(beamformMethod-0) > 1){      cerr<<"beamformMethod: "<<beamformMethod<<endl; errCnt++; }
      if( recoVertexingMode != 0){                cerr<<"recoVertexingMode: "<<recoVertexingMode<<endl; errCnt++; }
      if( (unsigned)(getSkymapMode-0) > 1){       cerr<<"getSkymapMode: "<<getSkymapMode<<endl; errCnt++; }
      if( string(recoPolType) != "vpol" && string(recoPolType) != "hpol" && string(recoPolType) != "both" )
        { cerr<<"Undefined recoPolType. Must be \"vpol\", \"hpol\", or \"both\"."<<endl; errCnt++; }
      if( nSideExpStart < 1 ){                    cerr<<"nSideExpStart: "<<nSideExpStart<<endl; errCnt++; }
      if( nSideExpEnd < 1 ){                      cerr<<"nSideExpEnd: "<<nSideExpEnd<<endl; errCnt++; }
      if( (unsigned)(nchnlFilter-0) > 3){         cerr<<"nchnlFilter: "<<nchnlFilter<<endl; errCnt++; }
      if( (unsigned)(nchnlCut-0) > 16){           cerr<<"nchnlCut: "<<nchnlCut<<endl; errCnt++; }
      if( nchnlThreshold < 0){                    cerr<<"nchnlThreshold: "<<nchnlThreshold<<endl; errCnt++; }
      if( nchnlThreshold_A1 < 0){                 cerr<<"nchnlThreshold_A1: "<<nchnlThreshold_A1<<endl; errCnt++; }
      if( nchnlThreshold_A2 < 0){                 cerr<<"nchnlThreshold_A2: "<<nchnlThreshold_A2<<endl; errCnt++; }
      if( nchnlThreshold_A2 < 0){                 cerr<<"nchnlThreshold_A3: "<<nchnlThreshold_A3<<endl; errCnt++; }
      if( (unsigned)(iceModel-0) > 1){            cerr<<"iceModel: "<<iceModel<<endl; errCnt++; }
      if( topN < 0 ){                             cerr<<"topN: "<<topN<<endl; errCnt++; }
      if( layerFirstRadius < 0 ){                 cerr<<"layerFirstRadius: "<<layerFirstRadius<<endl; errCnt++; }
      if( layerLastRadius < 0 ){                  cerr<<"layerLastRadius: "<<layerLastRadius<<endl; errCnt++; }
      if( layerFirstRadius > layerLastRadius ){   cerr<<"layerLastRadius - layerFirstRadius: "<<layerLastRadius - layerFirstRadius<<endl; errCnt++; }
      if( (unsigned)(recordMapData-0) > 1){       cerr<<"recordMapData: "<<recordMapData<<endl; errCnt++; }
      //if( (unsigned)(computePValue-0) > 1){       cerr<<"computePValue: "<<computePValue<<endl; errCnt++; }
      if( (unsigned)(computeLLHAndPValue-0) > 1){ cerr<<"computeLLHAndPValue: "<<computeLLHAndPValue<<endl; errCnt++; }
      if( /*computePValue != 0 ||*/ computeLLHAndPValue != 0){
         ifstream rmff(referenceMapFitFile/*.c_str()*/);
         if( rmff.is_open() ) rmff.close();
         else{                                    cerr<<"referenceMapFitFile: "<<referenceMapFitFile<<endl; errCnt++; }
         if( referenceMapFitFunc == "" ){         cerr<<"referenceMapFitFunc: "<<referenceMapFitFunc<<endl; errCnt++; }
      }
      if( string(openCLDeviceType) != "cpu" && string(openCLDeviceType) != "gpu"){ cerr<<"Undefined openCLDeviceType. Must be \"cpu\" or \"gpu\"."<<endl; errCnt++; }
      if( openCLMaxNumberOfCores < 0 ){           cerr<<"openCLMaxNumberOfCores"<<openCLMaxNumberOfCores<<endl; errCnt++; }
      if( maxNumberOfReco < -1 ){                 cerr<<"maxNumberOfReco: "<<maxNumberOfReco<<endl; errCnt++; }
      if( (unsigned)(constantNFilter-0) > 1){     cerr<<"constantNFilter: "<<constantNFilter<<endl; errCnt++; }
      if( surfaceCutAngle < 0.f && surfaceCutAngle > 180.f){ cerr<<"surfaceCutAngle: "<<surfaceCutAngle<<endl; errCnt++; }
      if( nchnlThreshold_anotherPol < 0){ cerr<<"nchnlThreshold_anotherPol: "<<nchnlThreshold_anotherPol<<endl; errCnt++; }

      if(errCnt > 0) return false;

   } else { cerr<<"Unable to open "<<sf<<endl; return false; }

return true;
}
