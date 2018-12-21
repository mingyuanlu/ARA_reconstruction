#ifndef ANALYSISTOOLS_H
#define ANALYSISTOOLS_H

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <algorithm>

#include "RawIcrrStationEvent.h"
#include "RawAtriStationEvent.h"
#include "UsefulAraStationEvent.h"
#include "UsefulIcrrStationEvent.h"
#include "UsefulAtriStationEvent.h"
#include "AraGeomTool.h"

#include "TChain.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TObject.h"

#include "calibrationTools.h"
#include "calibrationToolsVs3.h"
#include "evProcessTools.h"
#include "Healpix_Onion.h"
#include "recoTools.h"
#include "recoSettings.h"
#include "recoData.h"
#include "trackEngine.h"

#include "Detector.h"
#include "Trigger.h"
#include "Settings.h"
#include "Report.h"
#include "Event.h"
#include "Position.h"
#include "signal.hh"
#include "IceModel.h"

using namespace std;

int getRunType(string STATION, int runNum);

class cutParameter : public TObject
{

public:

   double val, plus, minus;

   cutParameter();
   ~cutParameter();

   ClassDef(cutParameter, 1);
};

class ARA02_cutValues : public TObject
{
private:

protected:

public:

   ARA02_cutValues();
   ~ARA02_cutValues();

   //CW impulsivity
   cutParameter cwImpCut[5];

   //Thermal box cut
   cutParameter coherenceCut_inBand[5];
   cutParameter snrCut_inBand[5];
   cutParameter coherenceCut_outOfBand[5];
   cutParameter snrCut_outOfBand[5];

   //Impulsivity cut
   cutParameter impCut;

   //Calpulser cut
   cutParameter zenMin[4];
   cutParameter zenMax[4];
   cutParameter aziMin[4];
   cutParameter aziMax[4];

   //Surface cut
   cutParameter surfaceCut_constantN;
   cutParameter surfaceCut_iterReco;

   void setValue(cutParamater& param, double _val, double _plus, double _minus);

   ClassDef(ARA02_cutValues, 1);

};

#endif
