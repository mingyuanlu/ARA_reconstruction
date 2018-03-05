#ifndef RECOSETTINGS_H
#define RECOSETTINGS_H

#include "TObject.h"
#include <string>

#define CSTRING_MAX 200

using namespace std;

class recoSettings : public TObject
{

private:

protected:

public:

   recoSettings();
   ~recoSettings();
   void initialize();
   bool readRecoSetupFile(string recoSetupFile);

   /*

   programFile:
      Default kernel_3D_analysis.c. Program file defining OpenCL kernels.
   nSideExp:
      Default 2. Number of pixels on a skympa = 12 * 2^nSideExp * 2^nSideExp
   nLayer:
      Default 1. Number of "layers" in the Healpix_Onion. I.e., number of skymaps at different radial distances.
   dataType:
      Default 0. 0: AraSim events. 1: real events.
   triggerCode:
      Default 111. This is a three digit code specifying what kind of events to reco. Only meaningful when analyzing real data.
      1st digit: RF non-calpulser trigger. 2nd digit: RF calpulser trigger. 3rd digit: software trigger. 0: do not reco. 1: do reco
   layerAllocationMode:
      Default 0. 0: nonlinear separation betweem layers. 1: linear separation between layers.
   skymapSearchMode:
      Default 0. 0: one pixelation throughout the entire search. 1: multi-step grid search with increasing pixelation - aka zoom approach.
   beamformMethod:
      Default 0. 0: cross-correlation coherence method. 1: coherently-summed waveforms (CSW) method.
   recoVertexingMode:
      Default 0. 0: reco vertex = pixel with max beam value.
   getSkymapMode:
      Default 0. 0: save skymap giving the reco vertex. 1: do not save.
   recoPolType:
      Default "vpol". Three choices are available: "vpol", "hpol", "both". This determines the polarization of antennas going into the reco. NB: only vpol or hpol is supported now.
   nSideExpStart:
      Default 2. This is only used when skymapSearchMode == 1, ie zoom search mode is used. This number gives the starting pixelization.
   nSideExpEnd:
      Default 8. This is only used when skymapSearchMode == 1, ie zoom search mode is used. This number gives the ending pixelization.
   nchnlFilter:
      Default 0. 0: no nchnl filter applied. 1: vpol nchnl filter applied. 2: hpol nchnl filter applied. 3: both pol nchnl filter applied.
   nchnlCut:
      Default 0. The least number of good channels in an event to clear the nchnl filter. Note that this may be dependent on "nchnlFilter" parameter.
   nchnlThreshold:
      Default 0. The voltage threshold in unit of noise rms used to determine if a channel is good (that is, has signal).
   nchnlThreshold_anotherPol:
      Deafult 0. The voltage threshold as above, but for another polarization if nchnlFilter == 1 or 2. In cases where nchnlFilter = 0 or 3, this parameter will not be used.
   nchnlThreshold_A1:
      Default 0. The nchnlThreshold value for A1 to reduce the noise rate to 1%, assuming nchnlFilter == 1 && nchnlCut == 3. This should be a value determined
      by examining RF data of A1.
   nchnlThreshold_A2:
      Default 0. The nchnlThreshold value for A2 to reduce the noise rate to 1%, assuming nchnlFilter == 1 && nchnlCut == 3. This should be a value determined
      by examining RF data of A2.
   nchnlThreshold_A3:
      Default 0. The nchnlThreshold value for A3 to reduce the noise rate to 1%, assuming nchnlFilter == 1 && nchnlCut == 3. This should be a value determined
      by examining RF data of A3.
   constantNFilter:
      Default 0. 0: no constant-N single-layer no-bound (at 5km) reconstruction applied. This is used to detect specifically plane-waves from above the ice surface. 1: the above reconstruction is applied and one can cut on its zen value as a surface cut. Cut value set in surfaceCutAngle parameter.
   surfaceCutAngle:
      Default 0. Surface angle cut below which an event will be considered surface event and filtered out.
   iceModel:
      Default 0. 0: depth-dependent IoR. Radiospline will be used to compute delays. 1: Bulk ice with a single IoR value of 1.76. Direct ray paths will be used
      to compute delays.
   topN:
      Default 50. Number of pixels with the most beam values to be selected.
   layerFirstRadius:
      Default 40. The radial distance in meters of the first layer in the Healpix_Onion.
   layerLastRadius:
      Default 5000. The radial distance in meters of the last layer in the Healpix_Onion.
   recordMapData:
      Default 0. 0: do not record map data histograms. 1: do record map data histograms and save to final output ROOT file
   computeLLHAndPValue:
      Default 0. 0: do not compute LLH & p-vale. 1: compute event LLH & p-value according to given reference map fit file and fit function.
   referenceMapFitFile:
      Default "". ROOT file containing the fits to reference map test statistics.
   referenceMapFitFunc:
      Default "". Function used to fit reference map test statistics distribution. Should be either "expo" or "gaus".
   openCLDeviceType:
      Default "cpu". Specifies which OpenCL devices to looks for in the platform. Currently only a single CPU is supported. In the future, should implement CPU, GPU, and a combination of both.
   openCLMaxNumberOfCores:
      Default 0. The upper limit of CPU cores to use in the OpenCL device. This uses the "Device fission extensions" of OpenCL and as I understand applies only to CPU devices. At default (0), no fission will be made and the program will in principle try to use up all given cores.
   remark:
      Default "". Any remark one wishes to add to the reco setup file. The remarks will then be carried along in the analysis output ROOT file. Note that number of characters should not exceed CSTRING_MAX defined in recoSettings.h
*/

   //string programFile;
   char programFile[CSTRING_MAX];
   int nSideExp;
   int nLayer;
   int dataType;
   char triggerCode[CSTRING_MAX];
   int layerAllocationMode;
   int skymapSearchMode;
   int beamformMethod;
   int recoVertexingMode;
   int getSkymapMode;
   //string recoPolType;
   char recoPolType[CSTRING_MAX];

   int nSideExpStart;
   int nSideExpEnd;

   int nchnlFilter;
   int nchnlCut;
   double nchnlThreshold;
   double nchnlThreshold_A1;
   double nchnlThreshold_A2;
   double nchnlThreshold_A3;

   int iceModel;
   int topN;

   double layerFirstRadius;
   double layerLastRadius;

   int recordMapData;
   //int computePValue;
   int computeLLHAndPValue;
   //string referenceMapFitFile;
   //string referenceMapFitFunc;
   char referenceMapFitFile[CSTRING_MAX];
   char referenceMapFitFunc[CSTRING_MAX];

   char openCLDeviceType[CSTRING_MAX];
   int openCLMaxNumberOfCores;

   //string remark;
   char remark[CSTRING_MAX];

   //ClassDef 4
   int maxNumberOfReco;

   //ClassDef 5
   int constantNFilter;
   float surfaceCutAngle;

   //ClassDef 6
   double nchnlThreshold_anotherPol;

   //ClassDef 7
   int dropARA02D4BH;
   int dropARA03D4;

   //ClassDef 8
   int use2ndRayReco;

   ClassDef(recoSettings, 8); //2: convert all string parameters to char
                              //3: add openCLDeviceType and openCLMaxNumberOfDevices parameters
};

#endif
