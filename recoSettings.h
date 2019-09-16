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
      Default 0. 0: linear separation betweem layers. 1: nonlinear separation between layers.
   skymapSearchMode:
      Default 0. 0: one pixelation throughout the entire search. 1: multi-step grid search with increasing pixelation - aka zoom approach.
   beamformMethod:
      Default 0. 0: cross-correlation coherence method. 1: coherently-summed waveforms (CSW) method (Deprecated! Use at your own risk).
   recoVertexingMode:
      Default 0. 0: reco vertex = pixel with max beam value.
   getSkymapMode:
      Default 0. 0: save skymap giving the reco vertex. 1: do not save.
   recoPolType:
      Default "vpol". Three choices are available: "vpol", "hpol", "both". This determines the polarization of antennas going into the reco. NB: only vpol or hpol is supported now.
   nSideExpStartG:
      Default 2. This is only used when skymapSearchMode == 1, ie zoom search mode is used. This number gives the starting pixelization.
   nSideExpEnd:
      Default 8. This is only used when skymapSearchMode == 1, ie zoom search mode is used. This number gives the ending pixelization.
   nchnlFilter:
      Default 0. 0: no nchnl filter applied. 1: vpol nchnl filter applied. 2: hpol nchnl filter applied. 3: both pol nchnl filter applied.
   nchnlCut:
      Default 0. The least number of good channels in an event to clear the nchnl filter. Note that this may be dependent on "nchnlFilter" parameter.
   nchnlThreshold:
      Default 0. The SNR threshold used to determine if a channel is good (that is, has signal).
   nchnlThreshold_anotherPol:
      Deafult 0. The SNR threshold as above, but for another polarization if nchnlFilter == 1 or 2. In cases where nchnlFilter = 0 or 3, this parameter will not be used.
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
   wInt_V:
      Default 0.4. The sample step size in nanoseconds for interpolating Vpol waveforms
   wInt_H:
      Default 0.625. The sample step size in nanoseconds for interpolating Hpol waveforms
   wInt_both:
      Default 0.5. The sample step size in nanoseconds for interpolating both pol waveforms, if for some reason one decides to do interferometry with both Vpol and Hpol waveforms
   maxPaddedSample:
      Default 2048. The number of samples in a padded waveform right before it enters the interferometor. Should be a combination of powers of 2, 3, 5, 7.
   recoEventIndex:
      Default -1. The event index to reconstruct. Use eventNumber to specify. Only the event with this index will be reconstructed.
   windowingType:
      Default 0. 0: uniform window (no window). 1: Hann window (1/4). 2: Bartlett window
   maskSubThresholdChannels:
      Default 0. 0: do not mask sub-threshold channels. 1: mask sub-threshold channels.
   maskSaturatedChannels:
      Default 1. 0: do not mask saturated channels. 1: mask saturated channels.
   saturationVoltage_mV:
      Default 600. The saturation voltage value, in mV. Saturation is checked in both +/- voltages.
   runIterativeReconstruction:
      Default 0. 0: do not run iterative reconstruction. 1: run iterative reconstruction.
   offsetBlock_threshold_V:
      Default -20. The negative threshold in mV for a Vpol channel's rolling mean to go below to count as a offset-block channel.
   offsetBlock_threshold_H:
      Default -12. The negative threshold in mV for an Hpol channel's rolling mean to go below to count as a offset-block channel.
   offsetBlock_timeRangeCut:
      Default 40. The max time range, in ns, for an event's offset-block channels' peak times to be in to count as an offset-block event. Time range is defined as max(peak times) - min(peak times).
   cwFilter:
      Default 0. 0: do not apply CW filter. 1: apply CW filter.
   minCWCoincidence:
      Default 3. The minimum number of coincident (in terms of frequency peak) channels for the event to count as a CW event.
   impulsivityFilter:
      Default 0. 0: no average impulsivity filer.1: with average impulsivity filter.
   impulsivityThreshold:
      Default 0. The filter threshold value for average impulsivity.
   flattenSaturatedAmplitude:
      Default 0. Whether the amplitudes in a waveforms that exceeds the +-saturationVoltage_mV will be flattened to be equal to saturationVoltage_mV or not. 0: no. 1: yes.
   chanMask:
      Default 1111111111111111. 16-digit code specifying which channels will be used in reconstruction. 1: use. 0: don't use. This is to replace the implementation in the body of analysis.cxx.
   powerEnvIntDuration:
      Default 25ns. Time duration to integrate when using evProcessTools::getSqrtVoltageSquaredSummedWaveform.
   snrMode:
      Default 0. 0: V_peak/RMS (crest factor). 1: sqrt((1/N)*(Sum V^2)) within powerEnvIntDuration of time. 2: E_s+n - E_n / E_n / (1ns/T)
   applyA2Ch6Correction:
      Default: 1. 0: do not apply. 1: apply.
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

   //ClassDef 9
   double wInt_V;
   double wInt_H;
   double wInt_both;
   int maxPaddedSample;

   //ClassDef 10
   int recoEventIndex;

   //ClassDef 11
   int windowingType;

   //ClassDef 12
   int maskSubThresholdChannels;
   int maskSaturatedChannels;
   double saturationVoltage_mV;

   //ClassDef 13
   int runIterativeReconstruction;
   double offsetBlock_threshold_V;
   double offsetBlock_threshold_H;
   float offsetBlock_timeRangeCut;
   int cwFilter;
   int minCWCoincidence;

   //ClassDef 14
   int impulsivityFilter;
   double impulsivityThreshold;

   //ClassDef 15
   int flattenSaturatedAmplitude;

   //ClassDef 16
   char chanMask[CSTRING_MAX];

   //ClassDef 17
   int powerEnvIntDuration;

   //ClassDef 18
   int snrMode;

   //ClassDef 19
   int applyA2Ch6Correction;


   ClassDef(recoSettings, 19); //2: convert all string parameters to char
                              //3: add openCLDeviceType and openCLMaxNumberOfDevices parameters
};

#endif
