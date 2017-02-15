/*
 * On different machine, should run get_info once to determine how many platforms, how many devices,
 * what types of devices are available. Then should decide on which platform and which devices to use.
 * This example code is based on the cobaltgpu@icecube.wisc.edu server, which has 2 platforms:
 * Platform 0: 1 AMD CPU device
 * Platform 1: 4 Nvidia GPU devices
 *
 * In this example, we will use platform 0
 *
 */
#include "evProcessTools.h"
#include "recoTools.h"
#include "TObject.h"

//ClassImp(Healpix_Onion);
ClassImp(recoSettings);
ClassImp(recoData);

using namespace std;
/*
#ifndef XCORRSUMGRAPH
TGraph *sillygr = new TGraph();
//TGraph *envelopeSum = new TGraph();
#define XCORRSUMGRAPH
#endif
*/

int setupCLRecoEnv(recoEnvData *clEnv, const char *programFile){

   cl_int i, err;
   cl_uint num_platforms=2;
   cl_uint num_devices=1;
   FILE *program_handle;
   size_t program_size, log_size;
   char *program_buffer;
   char *program_log;
   char name_data[1024];
   cl_uint ref_count;

/*
 * Get platform and devices
 */
   clEnv->platforms = (cl_platform_id*)malloc(sizeof(cl_platform_id)*num_platforms); //already know there are 2 platforms in total
   clGetPlatformIDs(num_platforms, clEnv->platforms, NULL);
   clGetPlatformInfo(clEnv->platforms[0], CL_PLATFORM_NAME, sizeof(name_data), &name_data, NULL);
   cout<<"Platform 0 name: "<<name_data<<endl;
   clEnv->devices = (cl_device_id*)malloc(sizeof(cl_device_id)*num_devices);
   err = clGetDeviceIDs(clEnv->platforms[0], CL_DEVICE_TYPE_CPU, num_devices, clEnv->devices, NULL);
   if(err<0){ cerr<<"No device found"<<endl; return -1; }

/*
 * Create context
 */
   clEnv->context = clCreateContext(NULL, num_devices, clEnv->devices, NULL, NULL, &err);
   clGetContextInfo(clEnv->context, CL_CONTEXT_REFERENCE_COUNT, sizeof(ref_count), &ref_count, NULL);
   cout<<"context reference count: "<<ref_count<<endl;

/*
 * Read in program source code
 */
   program_handle = fopen(programFile, "r");
   fseek(program_handle, 0, SEEK_END);
   program_size = ftell(program_handle);
   rewind(program_handle);
   program_buffer = (char*)malloc(program_size+1);
   program_buffer[program_size] = '\0';
   fread(program_buffer, sizeof(char), program_size, program_handle);
   fclose(program_handle);

/*
 * Create and build program
 * Program is built one device at a time in this example
 */
   clEnv->program = clCreateProgramWithSource(clEnv->context, 1 //assume only 1 program source code
                                    , (const char**)&program_buffer
                                    , &program_size, &err);
   for(i=0; i<num_devices; i++){

      err = clBuildProgram(clEnv->program, 1, &clEnv->devices[i], NULL, NULL, NULL);
      if(err<0){
      clGetProgramBuildInfo(clEnv->program, clEnv->devices[i], CL_PROGRAM_BUILD_LOG,
                            0, NULL, &log_size);
      program_log=(char*)malloc(log_size);
      clGetProgramBuildInfo(clEnv->program, clEnv->devices[i], CL_PROGRAM_BUILD_LOG,
                            log_size, program_log, NULL);
      printf("Device %u:\n%s\n", i, program_log);
      free(program_log);
   }
}

/*
 * Create kernels
 */
   cl_uint num_kernels;
   err = clCreateKernelsInProgram(clEnv->program, 0, NULL, &num_kernels);
   if(err<0){ cerr<<"No kernel found"<<endl; return -1; }
   cout<<"Number of kernels: "<<num_kernels<<endl;
   clEnv->kernels=(cl_kernel*)malloc(sizeof(cl_kernel)*num_kernels);
   clCreateKernelsInProgram(clEnv->program, num_kernels, clEnv->kernels, NULL);

/*
 * Create command queue on device 0
 */
   clEnv->queue = clCreateCommandQueue(clEnv->context, clEnv->devices[0], 0, &err);


/*
 * Example of getting kernel info
 */

   //cl_kernel shiftWf;
   //cl_kernel sumWf;
   //cl_kernel wfPwr;
   //cl_kernel getRecoDelays;

   for(int i=0; i<num_kernels; i++){
   clGetKernelInfo(clEnv->kernels[i], CL_KERNEL_FUNCTION_NAME,
                sizeof(name_data), name_data, NULL);
   cout<<"Kernel "<<i<<" name: "<<name_data<<endl;
   string kernel_name(name_data);

   if(kernel_name=="shiftWf")                      clEnv->shiftWf                =clEnv->kernels[i];
   else if (kernel_name=="sumWf")                  clEnv->sumWf                  =clEnv->kernels[i];
   else if (kernel_name=="wfPwr")                  clEnv->wfPwr                  =clEnv->kernels[i];
   else if (kernel_name=="xCorrWf")                clEnv->xCorrWf                =clEnv->kernels[i];
   else if (kernel_name=="computeXCorrCoef")       clEnv->computeXCorrCoef       =clEnv->kernels[i];
   else if (kernel_name=="computeCoherence")       clEnv->computeCoherence       =clEnv->kernels[i];
   else if (kernel_name=="bandPassFilter")         clEnv->bandPassFilter         =clEnv->kernels[i];
   else if (kernel_name=="getMaxPixInfoEachLayer") clEnv->getMaxPixInfoEachLayer =clEnv->kernels[i];
   //else if (kernel_name=="getRecoDelays") getRecoDelays=kernels[i];
   else if (kernel_name=="computeNormalizedCoherence") clEnv->computeNormalizedCoherence =clEnv->kernels[i];
   else if (kernel_name=="computeXCorrCoef_overlapCorrection") clEnv->computeXCorrCoef_overlapCorrection =clEnv->kernels[i];
   else { cerr<<"Invalid kernel name!\n"; return -1; }

   }

/*
 * Set up clFFT
 */

   err = clfftInitSetupData(&clEnv->fftSetup);
   err = clfftSetup(&clEnv->fftSetup);

   free(program_buffer);
   return 0;
}



//   recoData::recoData(){ recoData::initialize(); }
//   recoData::~recoData(){ /* default destructor */ }
/*
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

   }



   void recoData::setAllData(
     double w
   , float zen_true, float azi_true, float zen_reco, float azi_reco, float r_true, float r_reco
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

   weight = w;
   trueZen = zen_true;
   trueAzi = azi_true;
   recoZen = zen_reco;
   recoAzi = azi_reco;
   trueRadius = r_true;
   recoRadius = r_reco;

   for(int i=0; i<16; i++) recoChan[i] = usedChan[i];

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

   void recoData::setMaxPixInfo(int idx, float xCorrValue){

   maxPixIdx = idx;
   maxPixCoherence = xCorrValue;

   }
*/
   /*
   void recoData::setOnionInfo(Healpix_Onion *_onion){

   onion = _onion;

   }
   */
   /*
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

   void recoData::setFlag(float _flag){

   flag = _flag;

   }
*//*
   void recoData::duplicate(recoSettings *settings, recoData *old){

   int *_topMaxPixIdx         = (int*)calloc(old->topN, sizeof(int));
   float *_topMaxPixCoherence = (float*)calloc(old->topN, sizeof(float));
   for(int i=0; i<old->topN; i++){
      _topMaxPixIdx[i]       = old->topMaxPixIdx[i];
      _topMaxPixCoherence[i] = old->topMaxPixCoherence[i];
   }

   int *_maxPixIdxEachLayer       = (int*)calloc(settings->nLayer, sizeof(int));
   float *_maxPixCoherenceEachLayer = (float*)calloc(settings->nLayer, sizeof(float));
   for(int i=0; i<settings->nLayer; i++){
      _maxPixIdxEachLayer[i]       = old->maxPixIdxEachLayer[i];
      _maxPixCoherenceEachLayer[i] = old->maxPixCoherenceEachLayer[i];
   }

   setAllData( old->weight
             , old->trueZen,    old->trueAzi
             , old->recoZen,    old->recoAzi
             , old->trueRadius, old->recoRadius
             , old->recoChan
             , old->maxPixIdx,  old->maxPixCoherence
             //, old->onion
             , old->topN
             , _topMaxPixIdx, _topMaxPixCoherence
             , _maxPixIdxEachLayer, _maxPixCoherenceEachLayer
             , old->likelihood, old->pValue
             , old->inWindowSNR, old->unmodSNR
             , old->flag);
*/
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
*//*
   free(_topMaxPixIdx);
   free(_topMaxPixCoherence);
   free(_maxPixIdxEachLayer);
   free(_maxPixCoherenceEachLayer);

   }
*//*
   void recoData::clear(){

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
   }
*/

int reconstructCSW(unsigned int dataType, vector<TGraph *>& cleanEvent, recoEnvData *clEnv,
                float *recoDelays, float *recoDelays_V, float *recoDelays_H,
                int nDir, string pol, const int *chanMask, char *filename)
{

cout<<"Entered reconstructCSW method\n";
int nSamp;
int nAnt = (int)cleanEvent.size()/2; // Divide by 2 for only one polarization
//int totalNAnt = (int)cleanEvent.size();
int unmaskedNChan=0;
for(int ch=0; ch<2*nAnt; ch++) unmaskedNChan+=chanMask[ch];
cout<<"unmaskedNChan: "<<unmaskedNChan<<" nAnt: "<<nAnt<<endl;
float wInt;
if( dataType == 0 ) wInt = 0.5f; //AraSim event
else if( dataType == 1 ){ //real event
if( pol == "vpol" ) wInt = 0.4f;
else if ( pol == "hpol" ) wInt = 0.625f;
else { cerr<<"recoPolType undefined\n"; return -1; }
} else {
cerr<<"dataType undefined\n"; return -1; }
/*
 * Normalize wf
 */

TGraph *grNew[16];
for(int ch=0; ch<nAnt*2; ch++){
   //grNew[ch] = evProcessTools::getNormalizedGraph(cleanEvent[ch]); //For some unknown reason this line, partucularly the getNormalizedGraph function, is causing linking errors. I can't fix it.
   grNew[ch] = cleanEvent[ch];
}

cleanEvent.clear();
for(int ch=0; ch<nAnt*2; ch++){
   //if(chanMask[ch] == 1)
   cleanEvent.push_back(grNew[ch]);
}


/*
 * Loading voltsFlat array
 */

double t, v;
float *voltsFlat;

if( pol == "vpol" ){

   /* Using the 1st vpol wf for nSamp. Should make sure all vpol channels have the same nSamp. FIXME */
   nSamp = cleanEvent[0]->GetN();
   voltsFlat = (float*)calloc(nAnt*nSamp, sizeof(float));
   if(voltsFlat == NULL){ cerr<<"Null pointer to voltsFlat\n"; return -1; }

   for(int ch=0; ch<nAnt; ch++){

      if( chanMask[ch] ){
         for(int s=0; s<nSamp; s++){
            cleanEvent[ch]->GetPoint(s,t,v);
            /* Bartlett window applied in main analysis code */
            //voltsFlat[nSamp*ch + s] = static_cast<float>(v * FFTtools::bartlettWindow(s, nSamp));
            voltsFlat[nSamp*ch + s] = static_cast<float>(v);
         }
      } else {
         for(int s=0; s<nSamp; s++){
            voltsFlat[nSamp*ch + s] = 0.f;
         }
      }
   }

} else if ( pol == "hpol" ){

   /* Using the 1st hpol wf for nSamp. Should make sure all hpol channels have the same nSamp. FIXME */
   nSamp = cleanEvent[nAnt]->GetN();
   voltsFlat = (float*)calloc(nAnt*nSamp, sizeof(float));
   if(voltsFlat == NULL){ cerr<<"Null pointer to voltsFlat\n"; return -1; }

   for(int ch=0; ch<nAnt; ch++){

      if( chanMask[ch+nAnt] ){
         for(int s=0; s<nSamp; s++){
            cleanEvent[ch+nAnt]->GetPoint(s,t,v);
            /* Bartlett window applied in main analysis code */
            //voltsFlat[nSamp*ch + s] = static_cast<float>(v * FFTtools::bartlettWindow(s, nSamp));
            voltsFlat[nSamp*ch + s] = static_cast<float>(v);
         }
      } else {
         for(int s=0; s<nSamp; s++){
            voltsFlat[nSamp*ch + s] = 0.f;
         }
      }
   }

} else {
   cerr<<"recoPolType not defined\n"; return -1;
}
cout<<"voltsFlat loaded\n";

/*
 * Preparation for OUT_OF_PLACE transforms
 */

int interlvHermOutputSize = 2*(1 + nSamp/2); //Hermitian layout
int interlvOutputSize     = 2*nSamp;         //Not Hermitian layout
int planarHermOutputSize = (1 + nSamp/2);    //Hermitian planar layout

float *intensity_data_r, *intensity_data_c;
intensity_data_r = (float*)calloc(nAnt*planarHermOutputSize, sizeof(float));
intensity_data_c = (float*)calloc(nAnt*planarHermOutputSize, sizeof(float));

/*
 * FFT library related declarations
 */
//clfftPlanHandle planHandle;
clfftDim fftDim = CLFFT_1D;
size_t clLengths[1]   = {nSamp};
size_t clInStride[1]  = {1};
size_t clOutStride[1] = {1};
size_t inDist  = nSamp;
size_t outDist   = planarHermOutputSize;
size_t batchSize = nAnt;

/*
 * Set up clFFT
 */
//clfftSetupData fftSetup;
//err = clfftInitSetupData(&fftSetup);
//err = clfftSetup(&fftSetup);

/*
 * Prepare plan
 */
cout<<"Preparing plan...\n";
int err;
err = clfftCreateDefaultPlan(&clEnv->planHandle, clEnv->context, fftDim, clLengths);
err = clfftSetPlanPrecision(clEnv->planHandle, CLFFT_SINGLE);
err = clfftSetLayout(clEnv->planHandle, CLFFT_REAL, CLFFT_HERMITIAN_PLANAR);
err = clfftSetPlanScale(clEnv->planHandle, CLFFT_FORWARD, 1.f);
err = clfftSetPlanBatchSize(clEnv->planHandle, batchSize);
err = clfftSetPlanInStride(clEnv->planHandle, fftDim, clInStride);
err = clfftSetPlanOutStride(clEnv->planHandle, fftDim, clOutStride);
err = clfftSetPlanDistance(clEnv->planHandle, inDist, outDist);
err = clfftSetResultLocation(clEnv->planHandle, CLFFT_OUTOFPLACE);
err = clfftBakePlan(clEnv->planHandle, 1, &clEnv->queue, NULL, NULL);
cout<<"Plan prepared\n";
/* The plan is now ready to be executed */
cl_mem voltsFlatBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*nSamp*sizeof(float), NULL, &err);
err = clEnqueueWriteBuffer(clEnv->queue, voltsFlatBuffer, CL_TRUE, 0, nAnt*nSamp*sizeof(float), voltsFlat, 0, NULL, NULL);

cl_mem intensityRBuffer= clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*planarHermOutputSize*sizeof(float),
                         NULL, &err);
cl_mem intensityCBuffer= clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*planarHermOutputSize*sizeof(float),
                         NULL, &err);
cl_mem outBuffers[2] = {intensityRBuffer, intensityCBuffer};

cout<<"Enqueueing FFT\n";
err = clfftEnqueueTransform(clEnv->planHandle, CLFFT_FORWARD, 1, &clEnv->queue, 0, NULL, NULL,
                            &voltsFlatBuffer, outBuffers, NULL);
err = clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, outBuffers[0], CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_r,
                    0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, outBuffers[1], CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_c,
                    0, NULL, NULL);
cout<<"FFT done\n";
/*
 * Clean up CLFFT
 */
cout<<"Destroying plan...\n";
err = clfftDestroyPlan(&clEnv->planHandle);
//clfftTeardown();
cout<<"Plan destroyed\n";

/*
 * Shift FFTed wf according to reco directions
 */

cout<<"Shifting wfs...\n";
cl_mem recoDelaysBuffer;

if(pol == "vpol" ){
recoDelaysBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt*nDir,
                          recoDelays_V, &err);
/*
for(int i=0; i<nDir; i++){
   for(int j=0; j<nAnt; j++){
   cout<<recoDelays_V[i*nAnt+j]<<"\t";
   }
  cout<<endl;
} cout<<endl;
*/
}
else if (pol == "hpol" )
recoDelaysBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt*nDir,
                          recoDelays_H, &err);
else {
cerr<<"recoPolType not defined\n"; return -1;
}
cout<<"recoDelaysBuffer created\n";

float *shifted_data_r;
float *shifted_data_c;
shifted_data_r = (float*)calloc(nDir*nAnt*nSamp, sizeof(float));
shifted_data_c = (float*)calloc(nDir*nAnt*nSamp, sizeof(float));

cl_mem shiftedRBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nDir*nAnt*nSamp,
                          shifted_data_r, &err);
cl_mem shiftedCBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nDir*nAnt*nSamp,
                          shifted_data_c, &err);
cout<<"ShiftedR/CBuffer created\n";

clSetKernelArg(clEnv->shiftWf, 0, sizeof(cl_mem), &shiftedRBuffer  );
clSetKernelArg(clEnv->shiftWf, 1, sizeof(cl_mem), &shiftedCBuffer  );
clSetKernelArg(clEnv->shiftWf, 2, sizeof(cl_mem), &recoDelaysBuffer);
clSetKernelArg(clEnv->shiftWf, 3, sizeof(cl_mem), &intensityRBuffer );
clSetKernelArg(clEnv->shiftWf, 4, sizeof(cl_mem), &intensityCBuffer );
clSetKernelArg(clEnv->shiftWf, 5, sizeof(float),  &wInt);

unsigned int workDim=3;
size_t globalWorkSize[3]={(size_t)nDir, (size_t)nAnt, (size_t)nSamp};
clEnqueueNDRangeKernel(clEnv->queue, clEnv->shiftWf, workDim, NULL, globalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, shiftedRBuffer, CL_TRUE, 0, sizeof(float)*nDir*nAnt*nSamp, shifted_data_r, 0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, shiftedCBuffer, CL_TRUE, 0, sizeof(float)*nDir*nAnt*nSamp, shifted_data_c, 0, NULL, NULL);
cout<<"Shifted wfs\n";

/*
 * Make beam from summing shifted wfs
 */

cout<<"Making beam...\n";
float *beam_data_r, *beam_data_c;
beam_data_r = (float*)calloc(nDir*nSamp, sizeof(float));
beam_data_c = (float*)calloc(nDir*nSamp, sizeof(float));

cl_mem beamRBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nDir*nSamp,
                     beam_data_r, &err);
cl_mem beamCBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nDir*nSamp,
                     beam_data_c, &err);

clSetKernelArg(clEnv->sumWf, 0, sizeof(cl_mem), &beamRBuffer);
clSetKernelArg(clEnv->sumWf, 1, sizeof(cl_mem), &beamCBuffer);
clSetKernelArg(clEnv->sumWf, 2, sizeof(cl_mem), &shiftedRBuffer);
clSetKernelArg(clEnv->sumWf, 3, sizeof(cl_mem), &shiftedCBuffer);
clSetKernelArg(clEnv->sumWf, 4, sizeof(int),    &nAnt  );

workDim=2;
size_t sumGlobalWorkSize[2]={(size_t)nDir, (size_t)nSamp};
clEnqueueNDRangeKernel(clEnv->queue, clEnv->sumWf, workDim, NULL, sumGlobalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, beamRBuffer, CL_TRUE, 0, sizeof(float)*nDir*nSamp, beam_data_r, 0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, beamCBuffer, CL_TRUE, 0, sizeof(float)*nDir*nSamp, beam_data_c, 0, NULL, NULL);
cout<<"Beam made\n";

/*
 * Compute beam pwr in each reco direction
 */

cout<<"Computing beam pwr\n";
float *beamPwr;
beamPwr = (float*)calloc(nDir, sizeof(float));

cl_mem pwrBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nDir,
                   beamPwr, &err);
if(err<0) cerr<<"pwrBuffer creation fails"<<endl;

clSetKernelArg(clEnv->wfPwr, 0, sizeof(cl_mem), &pwrBuffer);
clSetKernelArg(clEnv->wfPwr, 1, sizeof(cl_mem), &beamRBuffer);
clSetKernelArg(clEnv->wfPwr, 2, sizeof(cl_mem), &beamCBuffer);
clSetKernelArg(clEnv->wfPwr, 3, sizeof(int),    &nSamp);

workDim=1;
size_t pwrGlobalWorkSize=nDir;
clEnqueueNDRangeKernel(clEnv->queue, clEnv->wfPwr, workDim, NULL, &pwrGlobalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, pwrBuffer, CL_TRUE, 0, sizeof(float)*nDir, beamPwr, 0, NULL, NULL);
cout<<"Beam pwr computed\n";

/*
 * Write FITS file
 */
cout<<"Creating Healpix map and writing to FITS....\n";
arr<float> beamPwrArr = arr<float>(&beamPwr[0], (size_t)nDir);
Healpix_Map<float> skyMap = Healpix_Map<float>(beamPwrArr, HEALPIX_ORDERING);

fitshandle fitsOut;
//char filename[] = "testCSWSkyMap.fits";
remove(filename);
fitsOut.create(filename);

write_Healpix_map_to_fits(fitsOut, skyMap, PLANCK_FLOAT32);
cout<<"Healpix map written\n";

/*
 * Deallocate memories
 */

cout<<"Deallocating memories...\n";
clReleaseMemObject(recoDelaysBuffer);
clReleaseMemObject(intensityRBuffer);
clReleaseMemObject(intensityCBuffer);
clReleaseMemObject(voltsFlatBuffer);
clReleaseMemObject(beamCBuffer);
clReleaseMemObject(beamRBuffer);
clReleaseMemObject(pwrBuffer);
clReleaseMemObject(shiftedRBuffer);
clReleaseMemObject(shiftedCBuffer);
free(voltsFlat);
//free(volts);
//free(recoDelays);
free(intensity_data_r);
free(intensity_data_c);
free(shifted_data_r);
free(shifted_data_c);
free(beam_data_r);
free(beam_data_c);
free(beamPwr);
cout<<"Memories deallocated\n";

return 0;
}

int reconstructCSW(recoSettings *settings, vector<TGraph *>& cleanEvent, recoEnvData *clEnv,
                float *recoDelays, float *recoDelays_V, float *recoDelays_H,
                int nDir, const int *chanMask, char *filename)
{

cout<<"Entered reconstructCSW method\n";
int nSamp;
int nAnt = (int)cleanEvent.size()/2; // Divide by 2 for only one polarization
//int totalNAnt = (int)cleanEvent.size();
int unmaskedNChan=0;
string pol = string(settings->recoPolType);
for(int ch=0; ch<2*nAnt; ch++) unmaskedNChan+=chanMask[ch];
cout<<"unmaskedNChan: "<<unmaskedNChan<<" nAnt: "<<nAnt<<endl;
float wInt;
int dataType = settings->dataType;
if( dataType == 0 ) wInt = 0.5f; //AraSim event
else if( dataType == 1 ){ //real event
if( pol == "vpol" ) wInt = 0.4f;
else if ( pol == "hpol" ) wInt = 0.625f;
else { cerr<<"recoPolType undefined\n"; return -1; }
} else {
cerr<<"dataType undefined\n"; return -1; }
/*
 * Normalize wf
 */

TGraph *grNew[16];
for(int ch=0; ch<nAnt*2; ch++){
   //grNew[ch] = evProcessTools::getNormalizedGraph(cleanEvent[ch]); //For some unknown reason this line, partucularly the getNormalizedGraph function, is causing linking errors. I can't fix it.
   grNew[ch] = cleanEvent[ch];
}

cleanEvent.clear();
for(int ch=0; ch<nAnt*2; ch++){
   //if(chanMask[ch] == 1)
   cleanEvent.push_back(grNew[ch]);
}


/*
 * Loading voltsFlat array
 */

double t, v;
float *voltsFlat;

if( pol == "vpol" ){

   /* Using the 1st vpol wf for nSamp. Should make sure all vpol channels have the same nSamp. FIXME */
   nSamp = cleanEvent[0]->GetN();
   voltsFlat = (float*)calloc(nAnt*nSamp, sizeof(float));
   if(voltsFlat == NULL){ cerr<<"Null pointer to voltsFlat\n"; return -1; }

   for(int ch=0; ch<nAnt; ch++){

      if( chanMask[ch] ){
         for(int s=0; s<nSamp; s++){
            cleanEvent[ch]->GetPoint(s,t,v);
            /* Bartlett window applied in main analysis code */
            //voltsFlat[nSamp*ch + s] = static_cast<float>(v * FFTtools::bartlettWindow(s, nSamp));
            voltsFlat[nSamp*ch + s] = static_cast<float>(v);
         }
      } else {
         for(int s=0; s<nSamp; s++){
            voltsFlat[nSamp*ch + s] = 0.f;
         }
      }
   }

} else if ( pol == "hpol" ){

   /* Using the 1st hpol wf for nSamp. Should make sure all hpol channels have the same nSamp. FIXME */
   nSamp = cleanEvent[nAnt]->GetN();
   voltsFlat = (float*)calloc(nAnt*nSamp, sizeof(float));
   if(voltsFlat == NULL){ cerr<<"Null pointer to voltsFlat\n"; return -1; }

   for(int ch=0; ch<nAnt; ch++){

      if( chanMask[ch+nAnt] ){
         for(int s=0; s<nSamp; s++){
            cleanEvent[ch+nAnt]->GetPoint(s,t,v);
            /* Bartlett window applied in main analysis code */
            //voltsFlat[nSamp*ch + s] = static_cast<float>(v * FFTtools::bartlettWindow(s, nSamp));
            voltsFlat[nSamp*ch + s] = static_cast<float>(v);
         }
      } else {
         for(int s=0; s<nSamp; s++){
            voltsFlat[nSamp*ch + s] = 0.f;
         }
      }
   }

} else {
   cerr<<"recoPolType not defined\n"; return -1;
}
cout<<"voltsFlat loaded\n";

/*
 * Preparation for OUT_OF_PLACE transforms
 */

int interlvHermOutputSize = 2*(1 + nSamp/2); //Hermitian layout
int interlvOutputSize     = 2*nSamp;         //Not Hermitian layout
int planarHermOutputSize = (1 + nSamp/2);    //Hermitian planar layout

float *intensity_data_r, *intensity_data_c;
intensity_data_r = (float*)calloc(nAnt*planarHermOutputSize, sizeof(float));
intensity_data_c = (float*)calloc(nAnt*planarHermOutputSize, sizeof(float));

/*
 * FFT library related declarations
 */
//clfftPlanHandle planHandle;
clfftDim fftDim = CLFFT_1D;
size_t clLengths[1]   = {nSamp};
size_t clInStride[1]  = {1};
size_t clOutStride[1] = {1};
size_t inDist  = nSamp;
size_t outDist   = planarHermOutputSize;
size_t batchSize = nAnt;

/*
 * Set up clFFT
 */
//clfftSetupData fftSetup;
//err = clfftInitSetupData(&fftSetup);
//err = clfftSetup(&fftSetup);

/*
 * Prepare plan
 */
cout<<"Preparing plan...\n";
int err;
err = clfftCreateDefaultPlan(&clEnv->planHandle, clEnv->context, fftDim, clLengths);
err = clfftSetPlanPrecision(clEnv->planHandle, CLFFT_SINGLE);
err = clfftSetLayout(clEnv->planHandle, CLFFT_REAL, CLFFT_HERMITIAN_PLANAR);
err = clfftSetPlanScale(clEnv->planHandle, CLFFT_FORWARD, 1.f);
err = clfftSetPlanBatchSize(clEnv->planHandle, batchSize);
err = clfftSetPlanInStride(clEnv->planHandle, fftDim, clInStride);
err = clfftSetPlanOutStride(clEnv->planHandle, fftDim, clOutStride);
err = clfftSetPlanDistance(clEnv->planHandle, inDist, outDist);
err = clfftSetResultLocation(clEnv->planHandle, CLFFT_OUTOFPLACE);
err = clfftBakePlan(clEnv->planHandle, 1, &clEnv->queue, NULL, NULL);
cout<<"Plan prepared\n";
/* The plan is now ready to be executed */
cl_mem voltsFlatBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*nSamp*sizeof(float), NULL, &err);
err = clEnqueueWriteBuffer(clEnv->queue, voltsFlatBuffer, CL_TRUE, 0, nAnt*nSamp*sizeof(float), voltsFlat, 0, NULL, NULL);

cl_mem intensityRBuffer= clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*planarHermOutputSize*sizeof(float),
                         NULL, &err);
cl_mem intensityCBuffer= clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*planarHermOutputSize*sizeof(float),
                         NULL, &err);
cl_mem outBuffers[2] = {intensityRBuffer, intensityCBuffer};

cout<<"Enqueueing FFT\n";
err = clfftEnqueueTransform(clEnv->planHandle, CLFFT_FORWARD, 1, &clEnv->queue, 0, NULL, NULL,
                            &voltsFlatBuffer, outBuffers, NULL);
err = clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, outBuffers[0], CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_r,
                    0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, outBuffers[1], CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_c,
                    0, NULL, NULL);
cout<<"FFT done\n";
/*
 * Clean up CLFFT
 */
cout<<"Destroying plan...\n";
err = clfftDestroyPlan(&clEnv->planHandle);
//clfftTeardown();
cout<<"Plan destroyed\n";

/*
 * Shift FFTed wf according to reco directions
 */

cout<<"Shifting wfs...\n";
cl_mem recoDelaysBuffer;

if(pol == "vpol" ){
recoDelaysBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt*nDir,
                          recoDelays_V, &err);
/*
for(int i=0; i<nDir; i++){
   for(int j=0; j<nAnt; j++){
   cout<<recoDelays_V[i*nAnt+j]<<"\t";
   }
  cout<<endl;
} cout<<endl;
*/
}
else if (pol == "hpol" )
recoDelaysBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt*nDir,
                          recoDelays_H, &err);
else {
cerr<<"recoPolType not defined\n"; return -1;
}
cout<<"recoDelaysBuffer created\n";

float *shifted_data_r;
float *shifted_data_c;
shifted_data_r = (float*)calloc(nDir*nAnt*nSamp, sizeof(float));
shifted_data_c = (float*)calloc(nDir*nAnt*nSamp, sizeof(float));

cl_mem shiftedRBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nDir*nAnt*nSamp,
                          shifted_data_r, &err);
cl_mem shiftedCBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nDir*nAnt*nSamp,
                          shifted_data_c, &err);
cout<<"ShiftedR/CBuffer created\n";

clSetKernelArg(clEnv->shiftWf, 0, sizeof(cl_mem), &shiftedRBuffer  );
clSetKernelArg(clEnv->shiftWf, 1, sizeof(cl_mem), &shiftedCBuffer  );
clSetKernelArg(clEnv->shiftWf, 2, sizeof(cl_mem), &recoDelaysBuffer);
clSetKernelArg(clEnv->shiftWf, 3, sizeof(cl_mem), &intensityRBuffer );
clSetKernelArg(clEnv->shiftWf, 4, sizeof(cl_mem), &intensityCBuffer );
clSetKernelArg(clEnv->shiftWf, 5, sizeof(float),  &wInt);

unsigned int workDim=3;
size_t globalWorkSize[3]={(size_t)nDir, (size_t)nAnt, (size_t)nSamp};
clEnqueueNDRangeKernel(clEnv->queue, clEnv->shiftWf, workDim, NULL, globalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, shiftedRBuffer, CL_TRUE, 0, sizeof(float)*nDir*nAnt*nSamp, shifted_data_r, 0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, shiftedCBuffer, CL_TRUE, 0, sizeof(float)*nDir*nAnt*nSamp, shifted_data_c, 0, NULL, NULL);
cout<<"Shifted wfs\n";

/*
 * Make beam from summing shifted wfs
 */

cout<<"Making beam...\n";
float *beam_data_r, *beam_data_c;
beam_data_r = (float*)calloc(nDir*nSamp, sizeof(float));
beam_data_c = (float*)calloc(nDir*nSamp, sizeof(float));

cl_mem beamRBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nDir*nSamp,
                     beam_data_r, &err);
cl_mem beamCBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nDir*nSamp,
                     beam_data_c, &err);

clSetKernelArg(clEnv->sumWf, 0, sizeof(cl_mem), &beamRBuffer);
clSetKernelArg(clEnv->sumWf, 1, sizeof(cl_mem), &beamCBuffer);
clSetKernelArg(clEnv->sumWf, 2, sizeof(cl_mem), &shiftedRBuffer);
clSetKernelArg(clEnv->sumWf, 3, sizeof(cl_mem), &shiftedCBuffer);
clSetKernelArg(clEnv->sumWf, 4, sizeof(int),    &nAnt  );

workDim=2;
size_t sumGlobalWorkSize[2]={(size_t)nDir, (size_t)nSamp};
clEnqueueNDRangeKernel(clEnv->queue, clEnv->sumWf, workDim, NULL, sumGlobalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, beamRBuffer, CL_TRUE, 0, sizeof(float)*nDir*nSamp, beam_data_r, 0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, beamCBuffer, CL_TRUE, 0, sizeof(float)*nDir*nSamp, beam_data_c, 0, NULL, NULL);
cout<<"Beam made\n";

/*
 * Compute beam pwr in each reco direction
 */

cout<<"Computing beam pwr\n";
float *beamPwr;
beamPwr = (float*)calloc(nDir, sizeof(float));

cl_mem pwrBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nDir,
                   beamPwr, &err);
if(err<0) cerr<<"pwrBuffer creation fails"<<endl;

clSetKernelArg(clEnv->wfPwr, 0, sizeof(cl_mem), &pwrBuffer);
clSetKernelArg(clEnv->wfPwr, 1, sizeof(cl_mem), &beamRBuffer);
clSetKernelArg(clEnv->wfPwr, 2, sizeof(cl_mem), &beamCBuffer);
clSetKernelArg(clEnv->wfPwr, 3, sizeof(int),    &nSamp);

workDim=1;
size_t pwrGlobalWorkSize=nDir;
clEnqueueNDRangeKernel(clEnv->queue, clEnv->wfPwr, workDim, NULL, &pwrGlobalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, pwrBuffer, CL_TRUE, 0, sizeof(float)*nDir, beamPwr, 0, NULL, NULL);
cout<<"Beam pwr computed\n";

/*
 * Write FITS file
 */
cout<<"Creating Healpix map and writing to FITS....\n";
arr<float> beamPwrArr = arr<float>(&beamPwr[0], (size_t)nDir);
Healpix_Map<float> skyMap = Healpix_Map<float>(beamPwrArr, RING);

fitshandle fitsOut;
//char filename[] = "testCSWSkyMap.fits";
remove(filename);
fitsOut.create(filename);

write_Healpix_map_to_fits(fitsOut, skyMap, PLANCK_FLOAT32);
cout<<"Healpix map written\n";

/*
 * Deallocate memories
 */

cout<<"Deallocating memories...\n";
clReleaseMemObject(recoDelaysBuffer);
clReleaseMemObject(intensityRBuffer);
clReleaseMemObject(intensityCBuffer);
clReleaseMemObject(voltsFlatBuffer);
clReleaseMemObject(beamCBuffer);
clReleaseMemObject(beamRBuffer);
clReleaseMemObject(pwrBuffer);
clReleaseMemObject(shiftedRBuffer);
clReleaseMemObject(shiftedCBuffer);
free(voltsFlat);
//free(volts);
//free(recoDelays);
free(intensity_data_r);
free(intensity_data_c);
free(shifted_data_r);
free(shifted_data_c);
free(beam_data_r);
free(beam_data_c);
free(beamPwr);
cout<<"Memories deallocated\n";

return 0;
}

int reconstructCSW_Serial(unsigned int dataType, vector<TGraph *>& cleanEvent,
                float *recoDelays, float *recoDelays_V, float *recoDelays_H,
                int nDir, string pol, const int *chanMask, char *filename)
                {

cout<<"Entered reconstructCSW_Serial method\n";
int nSamp;
int nAnt = (int)cleanEvent.size()/2; // Divide by 2 for only one polarization
//int totalNAnt = (int)cleanEvent.size();
int unmaskedNChan=0;
for(int ch=0; ch<2*nAnt; ch++) unmaskedNChan+=chanMask[ch];
cout<<"unmaskedNChan: "<<unmaskedNChan<<" nAnt: "<<nAnt<<endl;
float wInt;
if( dataType == 0 ) wInt = 0.5f; //AraSim event
else if( dataType == 1 ){ //real event
if( pol == "vpol" ) wInt = 0.4f;
else if ( pol == "hpol" ) wInt = 0.625f;
else { cerr<<"recoPolType undefined\n"; return -1; }
} else {
cerr<<"dataType undefined\n"; return -1; }

/*
 * Normalize wf
 */

TGraph *grNew[16];
for(int ch=0; ch<nAnt*2; ch++){
   //grNew[ch] = evProcessTools::getNormalizedGraph(cleanEvent[ch]);
   grNew[ch] = cleanEvent[ch];
}

cleanEvent.clear();
for(int ch=0; ch<nAnt*2; ch++){
   //if(chanMask[ch] == 1)
   cleanEvent.push_back(grNew[ch]);
}


/*
 * Loading voltsFlat array
 */

double t, v;
float *voltsFlat;

int aug = 2; //augmentation factor for shifted waveforms

if( pol == "vpol" ){

   /* Using the 1st vpol wf for nSamp. Should make sure all vpol channels have the same nSamp. FIXME */
   nSamp = cleanEvent[0]->GetN();
   voltsFlat = (float*)calloc(nAnt*nSamp*aug, sizeof(float));
   if(voltsFlat == NULL){ cerr<<"Null pointer to voltsFlat\n"; return -1; }

   for(int ch=0; ch<nAnt; ch++){

      if( chanMask[ch] ){
         for(int s=0; s<aug*nSamp; s++){

            if(s>nSamp*(aug-1)/2 && s<nSamp*(aug+1)/2){
            cleanEvent[ch]->GetPoint(s-nSamp*(aug-1)/2,t,v);
            /* Bartlett window applied in main analysis code */
            //voltsFlat[nSamp*ch + s] = static_cast<float>(v * FFTtools::bartlettWindow(s, nSamp));
            //cout<<"v: "<<v<<endl;
            voltsFlat[aug*nSamp*ch + s] = static_cast<float>(v);
            }
            else{
            voltsFlat[aug*nSamp*ch + s] = 0.f;
            }
         }
      } else {
         for(int s=0; s<aug*nSamp; s++){
            voltsFlat[aug*nSamp*ch + s] = 0.f;
         }
      }
   }

} else if ( pol == "hpol" ){

   /* Using the 1st hpol wf for nSamp. Should make sure all hpol channels have the same nSamp. FIXME */
   nSamp = cleanEvent[nAnt]->GetN();
   voltsFlat = (float*)calloc(nAnt*nSamp*aug, sizeof(float));
   if(voltsFlat == NULL){ cerr<<"Null pointer to voltsFlat\n"; return -1; }

   for(int ch=0; ch<nAnt; ch++){

      if( chanMask[ch+nAnt] ){
         for(int s=0; s<aug*nSamp; s++){

            if(s>nSamp*(aug-1)/2 && s<nSamp*(aug+1)/2){
            cleanEvent[ch+nAnt]->GetPoint(s-nSamp*(aug-1)/2,t,v);
            /* Bartlett window applied in main analysis code */
            //voltsFlat[nSamp*ch + s] = static_cast<float>(v * FFTtools::bartlettWindow(s, nSamp));
            voltsFlat[aug*nSamp*ch + s] = static_cast<float>(v);
            }
            else{
            voltsFlat[aug*nSamp*ch + s] = 0.f;
            }
         }
      } else {
         for(int s=0; s<aug*nSamp; s++){
            voltsFlat[aug*nSamp*ch + s] = 0.f;
         }
      }
   }

} else {
   cerr<<"recoPolType not defined\n"; return -1;
}
cout<<"voltsFlat loaded\n";

/*
 * Preparation for OUT_OF_PLACE transforms
 */
/*
int interlvHermOutputSize = 2*(1 + nSamp/2); //Hermitian layout
int interlvOutputSize     = 2*nSamp;         //Not Hermitian layout
int planarHermOutputSize = (1 + nSamp/2);    //Hermitian planar layout

float *intensity_data_r, *intensity_data_c;
intensity_data_r = (float*)calloc(nAnt*planarHermOutputSize, sizeof(float));
intensity_data_c = (float*)calloc(nAnt*planarHermOutputSize, sizeof(float));
*/
/*
 * FFT library related declarations
 */
//clfftPlanHandle planHandle;
/*
clfftDim fftDim = CLFFT_1D;
size_t clLengths[1]   = {nSamp};
size_t clInStride[1]  = {1};
size_t clOutStride[1] = {1};
size_t inDist  = nSamp;
size_t outDist   = planarHermOutputSize;
size_t batchSize = nAnt;
*/
/*
 * Set up clFFT
 */
//clfftSetupData fftSetup;
//err = clfftInitSetupData(&fftSetup);
//err = clfftSetup(&fftSetup);

/*
 * Prepare plan
 */
/*
cout<<"Preparing plan...\n";
int err;
err = clfftCreateDefaultPlan(&clEnv->planHandle, clEnv->context, fftDim, clLengths);
err = clfftSetPlanPrecision(clEnv->planHandle, CLFFT_SINGLE);
err = clfftSetLayout(clEnv->planHandle, CLFFT_REAL, CLFFT_HERMITIAN_PLANAR);
err = clfftSetPlanScale(clEnv->planHandle, CLFFT_FORWARD, 1.f);
err = clfftSetPlanBatchSize(clEnv->planHandle, batchSize);
err = clfftSetPlanInStride(clEnv->planHandle, fftDim, clInStride);
err = clfftSetPlanOutStride(clEnv->planHandle, fftDim, clOutStride);
err = clfftSetPlanDistance(clEnv->planHandle, inDist, outDist);
err = clfftSetResultLocation(clEnv->planHandle, CLFFT_OUTOFPLACE);
err = clfftBakePlan(clEnv->planHandle, 1, &clEnv->queue, NULL, NULL);
cout<<"Plan prepared\n";
*/
/* The plan is now ready to be executed */
/*
cl_mem voltsFlatBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*nSamp*sizeof(float), NULL, &err);
err = clEnqueueWriteBuffer(clEnv->queue, voltsFlatBuffer, CL_TRUE, 0, nAnt*nSamp*sizeof(float), voltsFlat, 0, NULL, NULL);

cl_mem intensityRBuffer= clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*planarHermOutputSize*sizeof(float),
                         NULL, &err);
cl_mem intensityCBuffer= clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*planarHermOutputSize*sizeof(float),
                         NULL, &err);
cl_mem outBuffers[2] = {intensityRBuffer, intensityCBuffer};

cout<<"Enqueueing FFT\n";
err = clfftEnqueueTransform(clEnv->planHandle, CLFFT_FORWARD, 1, &clEnv->queue, 0, NULL, NULL,
                            &voltsFlatBuffer, outBuffers, NULL);
err = clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, outBuffers[0], CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_r,
                    0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, outBuffers[1], CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_c,
                    0, NULL, NULL);
cout<<"FFT done\n";
*/
/*
 * Clean up CLFFT
 */
/*
cout<<"Destroying plan...\n";
err = clfftDestroyPlan(&clEnv->planHandle);
//clfftTeardown();
cout<<"Plan destroyed\n";
*/
/*
 * Shift FFTed wf according to reco directions
 */

cout<<"Shifting wfs...\n";
/*
cl_mem recoDelaysBuffer;

if(pol == "vpol" ){
recoDelaysBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt*nDir,
                          recoDelays_V, &err);
}
else if (pol == "hpol" )
recoDelaysBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt*nDir,
                          recoDelays_H, &err);
else {
cerr<<"recoPolType not defined\n"; return -1;
}
cout<<"recoDelaysBuffer created\n";

float *shifted_data_r;
float *shifted_data_c;
shifted_data_r = (float*)calloc(nDir*nAnt*nSamp, sizeof(float));
shifted_data_c = (float*)calloc(nDir*nAnt*nSamp, sizeof(float));

cl_mem shiftedRBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nDir*nAnt*nSamp,
                          shifted_data_r, &err);
cl_mem shiftedCBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nDir*nAnt*nSamp,
                          shifted_data_c, &err);
cout<<"ShiftedR/CBuffer created\n";

clSetKernelArg(clEnv->shiftWf, 0, sizeof(cl_mem), &shiftedRBuffer  );
clSetKernelArg(clEnv->shiftWf, 1, sizeof(cl_mem), &shiftedCBuffer  );
clSetKernelArg(clEnv->shiftWf, 2, sizeof(cl_mem), &recoDelaysBuffer);
clSetKernelArg(clEnv->shiftWf, 3, sizeof(cl_mem), &intensityRBuffer );
clSetKernelArg(clEnv->shiftWf, 4, sizeof(cl_mem), &intensityCBuffer );
clSetKernelArg(clEnv->shiftWf, 5, sizeof(float),  &wInt);

unsigned int workDim=3;
size_t globalWorkSize[3]={(size_t)nDir, (size_t)nAnt, (size_t)nSamp};
clEnqueueNDRangeKernel(clEnv->queue, clEnv->shiftWf, workDim, NULL, globalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, shiftedRBuffer, CL_TRUE, 0, sizeof(float)*nDir*nAnt*nSamp, shifted_data_r, 0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, shiftedCBuffer, CL_TRUE, 0, sizeof(float)*nDir*nAnt*nSamp, shifted_data_c, 0, NULL, NULL);
*/
cout<<"About to allocate shifted_data\n";
int extendedNSamp = aug*nSamp;
float *shifted_data = (float*)calloc(nDir*nAnt*extendedNSamp, sizeof(float));
cout<<"shifted_data allocated\n";
int delaysBin;
for(int i=0; i<nDir; i++){
   for(int j=0; j<nAnt; j++){

   //N.B.: here the existence of raytrace solution is _NOT_ checked!! For calpulser reconstruction this should not be a problem though
   if(pol == "vpol" )      delaysBin = (int)(recoDelays_V[i*nAnt+j]/wInt);
   else if(pol == "hpol" ) delaysBin = (int)(recoDelays_H[i*nAnt+j]/wInt);
   else { cerr<<"recoPolType not defined\n"; return -1; }

   for(int k=0; k<extendedNSamp; k++){

   if( (k+delaysBin >= 0 ) &&
       (k+delaysBin <= (extendedNSamp-1))
     )
   {
        //cout<<"voltsFlat: "<<voltsFlat[j*extendedNSamp + k +delaysBin]<<endl;
        shifted_data[i*nAnt*extendedNSamp+j*extendedNSamp+k] = voltsFlat[j*extendedNSamp + k + delaysBin];
   }
   else {
        //cout<<delaysBin<<"\t"<<k+delaysBin<<"\t"<<extendedNSamp-1<<endl;
        shifted_data[i*nAnt*extendedNSamp+j*extendedNSamp+k] = 0.f;
        //cout<<"Shifting out of range\n";
   }
   }


   }
}
cout<<"Shifted wfs\n";

/*
 * Make beam from summing shifted wfs
 */

cout<<"Making beam...\n";
/*
float *beam_data_r, *beam_data_c;
beam_data_r = (float*)calloc(nDir*nSamp, sizeof(float));
beam_data_c = (float*)calloc(nDir*nSamp, sizeof(float));

cl_mem beamRBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nDir*nSamp,
                     beam_data_r, &err);
cl_mem beamCBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nDir*nSamp,
                     beam_data_c, &err);

clSetKernelArg(clEnv->sumWf, 0, sizeof(cl_mem), &beamRBuffer);
clSetKernelArg(clEnv->sumWf, 1, sizeof(cl_mem), &beamCBuffer);
clSetKernelArg(clEnv->sumWf, 2, sizeof(cl_mem), &shiftedRBuffer);
clSetKernelArg(clEnv->sumWf, 3, sizeof(cl_mem), &shiftedCBuffer);
clSetKernelArg(clEnv->sumWf, 4, sizeof(int),    &nAnt  );

workDim=2;
size_t sumGlobalWorkSize[2]={(size_t)nDir, (size_t)nSamp};
clEnqueueNDRangeKernel(clEnv->queue, clEnv->sumWf, workDim, NULL, sumGlobalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, beamRBuffer, CL_TRUE, 0, sizeof(float)*nDir*nSamp, beam_data_r, 0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, beamCBuffer, CL_TRUE, 0, sizeof(float)*nDir*nSamp, beam_data_c, 0, NULL, NULL);
*/

float *beam_data = (float*)calloc(nDir*extendedNSamp,sizeof(float));

for(int i=0; i<nDir; i++){
   for(int j=0; j<extendedNSamp; j++){
      for(int k=0; k<nAnt; k++){

      //cout<<"shifted_data: "<<shifted_data[i*nAnt*extendedNSamp + k*extendedNSamp + j]<<endl;
      beam_data[i*extendedNSamp+j] += shifted_data[i*nAnt*extendedNSamp + k*extendedNSamp + j];

      }
      beam_data[i*extendedNSamp+j] /= (float)nAnt;
   }
}


cout<<"Beam made\n";

/*
 * Compute beam pwr in each reco direction
 */

cout<<"Computing beam pwr\n";
float *beamPwr;
beamPwr = (float*)calloc(nDir, sizeof(float));
/*
cl_mem pwrBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nDir,
                   beamPwr, &err);
if(err<0) cerr<<"pwrBuffer creation fails"<<endl;

clSetKernelArg(clEnv->wfPwr, 0, sizeof(cl_mem), &pwrBuffer);
clSetKernelArg(clEnv->wfPwr, 1, sizeof(cl_mem), &beamRBuffer);
clSetKernelArg(clEnv->wfPwr, 2, sizeof(cl_mem), &beamCBuffer);
clSetKernelArg(clEnv->wfPwr, 3, sizeof(int),    &nSamp);

workDim=1;
size_t pwrGlobalWorkSize=nDir;
clEnqueueNDRangeKernel(clEnv->queue, clEnv->wfPwr, workDim, NULL, &pwrGlobalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, pwrBuffer, CL_TRUE, 0, sizeof(float)*nDir, beamPwr, 0, NULL, NULL);
*/
for(int i=0; i<nDir; i++){

   for(int j=0; j<extendedNSamp; j++){
   //cout<<beam_data[i*extendedNSamp + j]<<endl;
   beamPwr[i] += (beam_data[i*extendedNSamp + j] * beam_data[i*extendedNSamp + j]);
   }

}

cout<<"Beam pwr computed\n";

/*
 * Write FITS file
 */
cout<<"Creating Healpix map and writing to FITS....\n";
arr<float> beamPwrArr = arr<float>(&beamPwr[0], (size_t)nDir);
Healpix_Map<float> skyMap = Healpix_Map<float>(beamPwrArr, HEALPIX_ORDERING);

fitshandle fitsOut;
//char filename[] = "testCSWSkyMap.fits";
remove(filename);
fitsOut.create(filename);

write_Healpix_map_to_fits(fitsOut, skyMap, PLANCK_FLOAT32);
cout<<"Healpix map written\n";

/*
 * Deallocate memories
 */

cout<<"Deallocating memories...\n";
//clReleaseMemObject(recoDelaysBuffer);
//clReleaseMemObject(intensityRBuffer);
//clReleaseMemObject(intensityCBuffer);
//clReleaseMemObject(voltsFlatBuffer);
//clReleaseMemObject(beamCBuffer);
//clReleaseMemObject(beamRBuffer);
//clReleaseMemObject(pwrBuffer);
//clReleaseMemObject(shiftedRBuffer);
//clReleaseMemObject(shiftedCBuffer);
free(voltsFlat);
//free(volts);
//free(recoDelays);
//free(intensity_data_r);
//free(intensity_data_c);
//free(shifted_data_r);
//free(shifted_data_c);
//free(beam_data_r);
//free(beam_data_c);
free(shifted_data);
free(beam_data);
free(beamPwr);
cout<<"Memories deallocated\n";

return 0;
}

int reconstructCSWGetMaxPix(unsigned int dataType, vector<TGraph *>& cleanEvent, recoEnvData *clEnv,
                float *recoDelays, float *recoDelays_V, float *recoDelays_H,
                int nDir, string pol, const int *chanMask, recoData *summary)
{

cout<<"Entered reconstructCSW routine\n";
int nSamp;
int nAnt = (int)cleanEvent.size()/2; // Divide by 2 for only one polarization
int unmaskedNChan=0;
for(int ch=0; ch<2*nAnt; ch++) unmaskedNChan+=chanMask[ch];
cout<<"unmaskedNChan: "<<unmaskedNChan<<" nAnt: "<<nAnt<<endl;
float wInt;
if( dataType == 0 ) wInt = 0.5f; //AraSim event
else if( dataType == 1 ){ //real event
if( pol == "vpol" ) wInt = 0.4f;
else if ( pol == "hpol" ) wInt = 0.625f;
else { cerr<<"recoPolType undefined\n"; return -1; }
} else {
cerr<<"dataType undefined\n"; return -1; }

/*
 * Loading voltsFlat array
 */

double t, v;
float *voltsFlat;

if( pol == "vpol" ){

   /* Using the 1st vpol wf for nSamp. Should make sure all vpol channels have the same nSamp. FIXME */
   nSamp = cleanEvent[0]->GetN();
   voltsFlat = (float*)calloc(nAnt*nSamp, sizeof(float));
   if(voltsFlat == NULL){ cerr<<"Null pointer to voltsFlat\n"; return -1; }

   for(int ch=0; ch<nAnt; ch++){

      if( chanMask[ch] ){
         for(int s=0; s<nSamp; s++){
            cleanEvent[ch]->GetPoint(s,t,v);
            /* Bartlett window applied in main analysis code */
            //voltsFlat[nSamp*ch + s] = static_cast<float>(v * FFTtools::bartlettWindow(s, nSamp));
            voltsFlat[nSamp*ch + s] = static_cast<float>(v);
         }
      } else {
         for(int s=0; s<nSamp; s++){
            voltsFlat[nSamp*ch + s] = 0.f;
         }
      }
   }

} else if ( pol == "hpol" ){

   /* Using the 1st hpol wf for nSamp. Should make sure all hpol channels have the same nSamp. FIXME */
   nSamp = cleanEvent[nAnt]->GetN();
   voltsFlat = (float*)calloc(nAnt*nSamp, sizeof(float));
   if(voltsFlat == NULL){ cerr<<"Null pointer to voltsFlat\n"; return -1; }

   for(int ch=0; ch<nAnt; ch++){
       if( chanMask[ch+nAnt] ){
         for(int s=0; s<nSamp; s++){
            cleanEvent[ch+nAnt]->GetPoint(s,t,v);
            /* Bartlett window applied in main analysis code */
            //voltsFlat[nSamp*ch + s] = static_cast<float>(v * FFTtools::bartlettWindow(s, nSamp));
            voltsFlat[nSamp*ch + s] = static_cast<float>(v);
         }
      } else {
         for(int s=0; s<nSamp; s++){
            voltsFlat[nSamp*ch + s] = 0.f;
         }
      }
   }

} else {
   cerr<<"recoPolType not defined\n"; return -1;
}
cout<<"voltsFlat loaded\n";

/*
 * Preparation for OUT_OF_PLACE transforms
 */

int interlvHermOutputSize = 2*(1 + nSamp/2); //Hermitian layout
int interlvOutputSize     = 2*nSamp;         //Not Hermitian layout
int planarHermOutputSize = (1 + nSamp/2);    //Hermitian planar layout

float *intensity_data_r, *intensity_data_c;
intensity_data_r = (float*)calloc(nAnt*planarHermOutputSize, sizeof(float));
intensity_data_c = (float*)calloc(nAnt*planarHermOutputSize, sizeof(float));

/*
 * FFT library related declarations
 */
//clfftPlanHandle planHandle;
clfftDim fftDim = CLFFT_1D;
size_t clLengths[1]   = {nSamp};
size_t clInStride[1]  = {1};
size_t clOutStride[1] = {1};
size_t inDist  = nSamp;
size_t outDist   = planarHermOutputSize;
size_t batchSize = nAnt;

/*
 * Set up clFFT
 */
//clfftSetupData fftSetup;
//err = clfftInitSetupData(&fftSetup);
//err = clfftSetup(&fftSetup);

/*
 * Prepare plan
 */
cout<<"Preparing plan...\n";
int err;
err = clfftCreateDefaultPlan(&clEnv->planHandle, clEnv->context, fftDim, clLengths);
err = clfftSetPlanPrecision(clEnv->planHandle, CLFFT_SINGLE);
err = clfftSetLayout(clEnv->planHandle, CLFFT_REAL, CLFFT_HERMITIAN_PLANAR);
err = clfftSetPlanScale(clEnv->planHandle, CLFFT_FORWARD, 1.f);
err = clfftSetPlanBatchSize(clEnv->planHandle, batchSize);
err = clfftSetPlanInStride(clEnv->planHandle, fftDim, clInStride);
err = clfftSetPlanOutStride(clEnv->planHandle, fftDim, clOutStride);
err = clfftSetPlanDistance(clEnv->planHandle, inDist, outDist);
err = clfftSetResultLocation(clEnv->planHandle, CLFFT_OUTOFPLACE);
err = clfftBakePlan(clEnv->planHandle, 1, &clEnv->queue, NULL, NULL);
cout<<"Plan prepared\n";
/* The plan is now ready to be executed */
cl_mem voltsFlatBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*nSamp*sizeof(float), NULL, &err);
err = clEnqueueWriteBuffer(clEnv->queue, voltsFlatBuffer, CL_TRUE, 0, nAnt*nSamp*sizeof(float), voltsFlat, 0, NULL, NULL);

cl_mem intensityRBuffer= clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*planarHermOutputSize*sizeof(float),
                         NULL, &err);
cl_mem intensityCBuffer= clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*planarHermOutputSize*sizeof(float),
                         NULL, &err);
cl_mem outBuffers[2] = {intensityRBuffer, intensityCBuffer};

cout<<"Enqueueing FFT\n";
err = clfftEnqueueTransform(clEnv->planHandle, CLFFT_FORWARD, 1, &clEnv->queue, 0, NULL, NULL,
                            &voltsFlatBuffer, outBuffers, NULL);
err = clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, outBuffers[0], CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_r,
                    0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, outBuffers[1], CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_c,
                    0, NULL, NULL);
cout<<"FFT done\n";
/*
 * Clean up CLFFT
 */
cout<<"Destroying plan...\n";
err = clfftDestroyPlan(&clEnv->planHandle);
//clfftTeardown();
cout<<"Plan destroyed\n";

/*
 * Shift FFTed wf according to reco directions
 */

cout<<"Shifting wfs...\n";
cl_mem recoDelaysBuffer;

if(pol == "vpol" ){
recoDelaysBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt*nDir,
                          recoDelays_V, &err);
/*
for(int i=0; i<nDir; i++){
   for(int j=0; j<nAnt; j++){
   cout<<recoDelays_V[i*nAnt+j]<<"\t";
   }
  cout<<endl;
} cout<<endl;
*/
}
else if (pol == "hpol" )
recoDelaysBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt*nDir,
                          recoDelays_H, &err);
else {
cerr<<"recoPolType not defined\n"; return -1;
}
cout<<"recoDelaysBuffer created\n";

float *shifted_data_r;
float *shifted_data_c;
shifted_data_r = (float*)calloc(nDir*nAnt*nSamp, sizeof(float));
shifted_data_c = (float*)calloc(nDir*nAnt*nSamp, sizeof(float));

cl_mem shiftedRBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nDir*nAnt*nSamp,
                          shifted_data_r, &err);
cl_mem shiftedCBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nDir*nAnt*nSamp,
                          shifted_data_c, &err);
cout<<"ShiftedR/CBuffer created\n";

clSetKernelArg(clEnv->shiftWf, 0, sizeof(cl_mem), &shiftedRBuffer  );
clSetKernelArg(clEnv->shiftWf, 1, sizeof(cl_mem), &shiftedCBuffer  );
clSetKernelArg(clEnv->shiftWf, 2, sizeof(cl_mem), &recoDelaysBuffer);
clSetKernelArg(clEnv->shiftWf, 3, sizeof(cl_mem), &intensityRBuffer );
clSetKernelArg(clEnv->shiftWf, 4, sizeof(cl_mem), &intensityCBuffer );
clSetKernelArg(clEnv->shiftWf, 5, sizeof(float),  &wInt);

unsigned int workDim=3;
size_t globalWorkSize[3]={(size_t)nDir, (size_t)nAnt, (size_t)nSamp};
clEnqueueNDRangeKernel(clEnv->queue, clEnv->shiftWf, workDim, NULL, globalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, shiftedRBuffer, CL_TRUE, 0, sizeof(float)*nDir*nAnt*nSamp, shifted_data_r, 0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, shiftedCBuffer, CL_TRUE, 0, sizeof(float)*nDir*nAnt*nSamp, shifted_data_c, 0, NULL, NULL);
cout<<"Shifted wfs\n";

/*
 * Make beam from summing shifted wfs
 */

cout<<"Making beam...\n";
float *beam_data_r, *beam_data_c;
beam_data_r = (float*)calloc(nDir*nSamp, sizeof(float));
beam_data_c = (float*)calloc(nDir*nSamp, sizeof(float));

cl_mem beamRBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nDir*nSamp,
                     beam_data_r, &err);
cl_mem beamCBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nDir*nSamp,
                     beam_data_c, &err);

clSetKernelArg(clEnv->sumWf, 0, sizeof(cl_mem), &beamRBuffer);
clSetKernelArg(clEnv->sumWf, 1, sizeof(cl_mem), &beamCBuffer);
clSetKernelArg(clEnv->sumWf, 2, sizeof(cl_mem), &shiftedRBuffer);
clSetKernelArg(clEnv->sumWf, 3, sizeof(cl_mem), &shiftedCBuffer);
clSetKernelArg(clEnv->sumWf, 4, sizeof(int),    &nAnt  );

workDim=2;
size_t sumGlobalWorkSize[2]={(size_t)nDir, (size_t)nSamp};
clEnqueueNDRangeKernel(clEnv->queue, clEnv->sumWf, workDim, NULL, sumGlobalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, beamRBuffer, CL_TRUE, 0, sizeof(float)*nDir*nSamp, beam_data_r, 0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, beamCBuffer, CL_TRUE, 0, sizeof(float)*nDir*nSamp, beam_data_c, 0, NULL, NULL);
cout<<"Beam made\n";

/*
 * Compute beam pwr in each reco direction
 */

cout<<"Computing beam pwr\n";
float *beamPwr;
beamPwr = (float*)calloc(nDir, sizeof(float));

cl_mem pwrBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nDir,
                   beamPwr, &err);
if(err<0) cerr<<"pwrBuffer creation fails"<<endl;

clSetKernelArg(clEnv->wfPwr, 0, sizeof(cl_mem), &pwrBuffer);
clSetKernelArg(clEnv->wfPwr, 1, sizeof(cl_mem), &beamRBuffer);
clSetKernelArg(clEnv->wfPwr, 2, sizeof(cl_mem), &beamCBuffer);
clSetKernelArg(clEnv->wfPwr, 3, sizeof(int),    &nSamp);

workDim=1;
size_t pwrGlobalWorkSize=nDir;
clEnqueueNDRangeKernel(clEnv->queue, clEnv->wfPwr, workDim, NULL, &pwrGlobalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, pwrBuffer, CL_TRUE, 0, sizeof(float)*nDir, beamPwr, 0, NULL, NULL);
cout<<"Beam pwr computed\n";

/*
 * Loop over beamPwr to find max power and its pix index
 */

float max=0.f;
int maxPixIdx;

for(int idx=0; idx<nDir; idx++){
   if(beamPwr[idx] > max){
      max = beamPwr[idx];
      maxPixIdx = idx;
   }
}

summary->setMaxPixInfo(maxPixIdx, max);

/*
 * Write FITS file
 */
/*
cout<<"Creating Halpix map and writing to FITS....\n";
arr<float> beamPwrArr = arr<float>(&beamPwr[0], (size_t)nDir);
Healpix_Map<float> skyMap = Healpix_Map<float>(beamPwrArr, RING);

fitshandle fitsOut;
char filename[] = "testSkyMap.fits";
remove(filename);
fitsOut.create(filename);

write_Healpix_map_to_fits(fitsOut, skyMap, PLANCK_FLOAT32);
cout<<"Healpix map written\n";
*/
/*
 * Deallocate memories
 */

cout<<"Deallocating memories...\n";
clReleaseMemObject(recoDelaysBuffer);
clReleaseMemObject(intensityRBuffer);
clReleaseMemObject(intensityCBuffer);
clReleaseMemObject(voltsFlatBuffer);
clReleaseMemObject(beamCBuffer);
clReleaseMemObject(beamRBuffer);
clReleaseMemObject(pwrBuffer);
clReleaseMemObject(shiftedRBuffer);
clReleaseMemObject(shiftedCBuffer);
free(voltsFlat);
//free(volts);
//free(recoDelays);
free(intensity_data_r);
free(intensity_data_c);
free(shifted_data_r);
free(shifted_data_c);
free(beam_data_r);
free(beam_data_c);
free(beamPwr);
cout<<"Memories deallocated\n";

return maxPixIdx;
}

int reconstructCSWGetMaxPix(recoSettings *settings, vector<TGraph *>& cleanEvent, recoEnvData *clEnv,
                float *recoDelays, float *recoDelays_V, float *recoDelays_H,
                int nDir, const int *chanMask, recoData *summary)
{

cout<<"Entered reconstructCSW routine\n";
int nSamp;
int nAnt = (int)cleanEvent.size()/2; // Divide by 2 for only one polarization
int unmaskedNChan=0;
string pol = string(settings->recoPolType);
for(int ch=0; ch<2*nAnt; ch++) unmaskedNChan+=chanMask[ch];
cout<<"unmaskedNChan: "<<unmaskedNChan<<" nAnt: "<<nAnt<<endl;
float wInt;
int dataType = settings->dataType;
if( dataType == 0 ) wInt = 0.5f; //AraSim event
else if( dataType == 1 ){ //real event
if( pol == "vpol" ) wInt = 0.4f;
else if ( pol == "hpol" ) wInt = 0.625f;
else { cerr<<"recoPolType undefined\n"; return -1; }
} else {
cerr<<"dataType undefined\n"; return -1; }

/*
 * Loading voltsFlat array
 */

double t, v;
float *voltsFlat;

if( pol == "vpol" ){

   /* Using the 1st vpol wf for nSamp. Should make sure all vpol channels have the same nSamp. FIXME */
   nSamp = cleanEvent[0]->GetN();
   voltsFlat = (float*)calloc(nAnt*nSamp, sizeof(float));
   if(voltsFlat == NULL){ cerr<<"Null pointer to voltsFlat\n"; return -1; }

   for(int ch=0; ch<nAnt; ch++){

      if( chanMask[ch] ){
         for(int s=0; s<nSamp; s++){
            cleanEvent[ch]->GetPoint(s,t,v);
            /* Bartlett window applied in main analysis code */
            //voltsFlat[nSamp*ch + s] = static_cast<float>(v * FFTtools::bartlettWindow(s, nSamp));
            voltsFlat[nSamp*ch + s] = static_cast<float>(v);
         }
      } else {
         for(int s=0; s<nSamp; s++){
            voltsFlat[nSamp*ch + s] = 0.f;
         }
      }
   }

} else if ( pol == "hpol" ){

   /* Using the 1st hpol wf for nSamp. Should make sure all hpol channels have the same nSamp. FIXME */
   nSamp = cleanEvent[nAnt]->GetN();
   voltsFlat = (float*)calloc(nAnt*nSamp, sizeof(float));
   if(voltsFlat == NULL){ cerr<<"Null pointer to voltsFlat\n"; return -1; }

   for(int ch=0; ch<nAnt; ch++){
       if( chanMask[ch+nAnt] ){
         for(int s=0; s<nSamp; s++){
            cleanEvent[ch+nAnt]->GetPoint(s,t,v);
            /* Bartlett window applied in main analysis code */
            //voltsFlat[nSamp*ch + s] = static_cast<float>(v * FFTtools::bartlettWindow(s, nSamp));
            voltsFlat[nSamp*ch + s] = static_cast<float>(v);
         }
      } else {
         for(int s=0; s<nSamp; s++){
            voltsFlat[nSamp*ch + s] = 0.f;
         }
      }
   }

} else {
   cerr<<"recoPolType not defined\n"; return -1;
}
cout<<"voltsFlat loaded\n";

/*
 * Preparation for OUT_OF_PLACE transforms
 */

int interlvHermOutputSize = 2*(1 + nSamp/2); //Hermitian layout
int interlvOutputSize     = 2*nSamp;         //Not Hermitian layout
int planarHermOutputSize = (1 + nSamp/2);    //Hermitian planar layout

float *intensity_data_r, *intensity_data_c;
intensity_data_r = (float*)calloc(nAnt*planarHermOutputSize, sizeof(float));
intensity_data_c = (float*)calloc(nAnt*planarHermOutputSize, sizeof(float));

/*
 * FFT library related declarations
 */
//clfftPlanHandle planHandle;
clfftDim fftDim = CLFFT_1D;
size_t clLengths[1]   = {nSamp};
size_t clInStride[1]  = {1};
size_t clOutStride[1] = {1};
size_t inDist  = nSamp;
size_t outDist   = planarHermOutputSize;
size_t batchSize = nAnt;

/*
 * Set up clFFT
 */
//clfftSetupData fftSetup;
//err = clfftInitSetupData(&fftSetup);
//err = clfftSetup(&fftSetup);

/*
 * Prepare plan
 */
cout<<"Preparing plan...\n";
int err;
err = clfftCreateDefaultPlan(&clEnv->planHandle, clEnv->context, fftDim, clLengths);
err = clfftSetPlanPrecision(clEnv->planHandle, CLFFT_SINGLE);
err = clfftSetLayout(clEnv->planHandle, CLFFT_REAL, CLFFT_HERMITIAN_PLANAR);
err = clfftSetPlanScale(clEnv->planHandle, CLFFT_FORWARD, 1.f);
err = clfftSetPlanBatchSize(clEnv->planHandle, batchSize);
err = clfftSetPlanInStride(clEnv->planHandle, fftDim, clInStride);
err = clfftSetPlanOutStride(clEnv->planHandle, fftDim, clOutStride);
err = clfftSetPlanDistance(clEnv->planHandle, inDist, outDist);
err = clfftSetResultLocation(clEnv->planHandle, CLFFT_OUTOFPLACE);
err = clfftBakePlan(clEnv->planHandle, 1, &clEnv->queue, NULL, NULL);
cout<<"Plan prepared\n";
/* The plan is now ready to be executed */
cl_mem voltsFlatBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*nSamp*sizeof(float), NULL, &err);
err = clEnqueueWriteBuffer(clEnv->queue, voltsFlatBuffer, CL_TRUE, 0, nAnt*nSamp*sizeof(float), voltsFlat, 0, NULL, NULL);

cl_mem intensityRBuffer= clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*planarHermOutputSize*sizeof(float),
                         NULL, &err);
cl_mem intensityCBuffer= clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*planarHermOutputSize*sizeof(float),
                         NULL, &err);
cl_mem outBuffers[2] = {intensityRBuffer, intensityCBuffer};

cout<<"Enqueueing FFT\n";
err = clfftEnqueueTransform(clEnv->planHandle, CLFFT_FORWARD, 1, &clEnv->queue, 0, NULL, NULL,
                            &voltsFlatBuffer, outBuffers, NULL);
err = clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, outBuffers[0], CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_r,
                    0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, outBuffers[1], CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_c,
                    0, NULL, NULL);
cout<<"FFT done\n";
/*
 * Clean up CLFFT
 */
cout<<"Destroying plan...\n";
err = clfftDestroyPlan(&clEnv->planHandle);
//clfftTeardown();
cout<<"Plan destroyed\n";

/*
 * Shift FFTed wf according to reco directions
 */

cout<<"Shifting wfs...\n";
cl_mem recoDelaysBuffer;

if(pol == "vpol" ){
recoDelaysBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt*nDir,
                          recoDelays_V, &err);
/*
for(int i=0; i<nDir; i++){
   for(int j=0; j<nAnt; j++){
   cout<<recoDelays_V[i*nAnt+j]<<"\t";
   }
  cout<<endl;
} cout<<endl;
*/
}
else if (pol == "hpol" )
recoDelaysBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt*nDir,
                          recoDelays_H, &err);
else {
cerr<<"recoPolType not defined\n"; return -1;
}
cout<<"recoDelaysBuffer created\n";

float *shifted_data_r;
float *shifted_data_c;
shifted_data_r = (float*)calloc(nDir*nAnt*nSamp, sizeof(float));
shifted_data_c = (float*)calloc(nDir*nAnt*nSamp, sizeof(float));

cl_mem shiftedRBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nDir*nAnt*nSamp,
                          shifted_data_r, &err);
cl_mem shiftedCBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nDir*nAnt*nSamp,
                          shifted_data_c, &err);
cout<<"ShiftedR/CBuffer created\n";

clSetKernelArg(clEnv->shiftWf, 0, sizeof(cl_mem), &shiftedRBuffer  );
clSetKernelArg(clEnv->shiftWf, 1, sizeof(cl_mem), &shiftedCBuffer  );
clSetKernelArg(clEnv->shiftWf, 2, sizeof(cl_mem), &recoDelaysBuffer);
clSetKernelArg(clEnv->shiftWf, 3, sizeof(cl_mem), &intensityRBuffer );
clSetKernelArg(clEnv->shiftWf, 4, sizeof(cl_mem), &intensityCBuffer );
clSetKernelArg(clEnv->shiftWf, 5, sizeof(float),  &wInt);

unsigned int workDim=3;
size_t globalWorkSize[3]={(size_t)nDir, (size_t)nAnt, (size_t)nSamp};
clEnqueueNDRangeKernel(clEnv->queue, clEnv->shiftWf, workDim, NULL, globalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, shiftedRBuffer, CL_TRUE, 0, sizeof(float)*nDir*nAnt*nSamp, shifted_data_r, 0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, shiftedCBuffer, CL_TRUE, 0, sizeof(float)*nDir*nAnt*nSamp, shifted_data_c, 0, NULL, NULL);
cout<<"Shifted wfs\n";

/*
 * Make beam from summing shifted wfs
 */

cout<<"Making beam...\n";
float *beam_data_r, *beam_data_c;
beam_data_r = (float*)calloc(nDir*nSamp, sizeof(float));
beam_data_c = (float*)calloc(nDir*nSamp, sizeof(float));

cl_mem beamRBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nDir*nSamp,
                     beam_data_r, &err);
cl_mem beamCBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nDir*nSamp,
                     beam_data_c, &err);

clSetKernelArg(clEnv->sumWf, 0, sizeof(cl_mem), &beamRBuffer);
clSetKernelArg(clEnv->sumWf, 1, sizeof(cl_mem), &beamCBuffer);
clSetKernelArg(clEnv->sumWf, 2, sizeof(cl_mem), &shiftedRBuffer);
clSetKernelArg(clEnv->sumWf, 3, sizeof(cl_mem), &shiftedCBuffer);
clSetKernelArg(clEnv->sumWf, 4, sizeof(int),    &nAnt  );

workDim=2;
size_t sumGlobalWorkSize[2]={(size_t)nDir, (size_t)nSamp};
clEnqueueNDRangeKernel(clEnv->queue, clEnv->sumWf, workDim, NULL, sumGlobalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, beamRBuffer, CL_TRUE, 0, sizeof(float)*nDir*nSamp, beam_data_r, 0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, beamCBuffer, CL_TRUE, 0, sizeof(float)*nDir*nSamp, beam_data_c, 0, NULL, NULL);
cout<<"Beam made\n";

/*
 * Compute beam pwr in each reco direction
 */

cout<<"Computing beam pwr\n";
float *beamPwr;
beamPwr = (float*)calloc(nDir, sizeof(float));

cl_mem pwrBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nDir,
                   beamPwr, &err);
if(err<0) cerr<<"pwrBuffer creation fails"<<endl;

clSetKernelArg(clEnv->wfPwr, 0, sizeof(cl_mem), &pwrBuffer);
clSetKernelArg(clEnv->wfPwr, 1, sizeof(cl_mem), &beamRBuffer);
clSetKernelArg(clEnv->wfPwr, 2, sizeof(cl_mem), &beamCBuffer);
clSetKernelArg(clEnv->wfPwr, 3, sizeof(int),    &nSamp);

workDim=1;
size_t pwrGlobalWorkSize=nDir;
clEnqueueNDRangeKernel(clEnv->queue, clEnv->wfPwr, workDim, NULL, &pwrGlobalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, pwrBuffer, CL_TRUE, 0, sizeof(float)*nDir, beamPwr, 0, NULL, NULL);
cout<<"Beam pwr computed\n";

/*
 * Loop over beamPwr to find max power and its pix index
 */

float max=0.f;
int maxPixIdx;

for(int idx=0; idx<nDir; idx++){
   if(beamPwr[idx] > max){
      max = beamPwr[idx];
      maxPixIdx = idx;
   }
}

summary->setMaxPixInfo(maxPixIdx, max);

/*
 * Write FITS file
 */
/*
cout<<"Creating Halpix map and writing to FITS....\n";
arr<float> beamPwrArr = arr<float>(&beamPwr[0], (size_t)nDir);
Healpix_Map<float> skyMap = Healpix_Map<float>(beamPwrArr, RING);

fitshandle fitsOut;
char filename[] = "testSkyMap.fits";
remove(filename);
fitsOut.create(filename);

write_Healpix_map_to_fits(fitsOut, skyMap, PLANCK_FLOAT32);
cout<<"Healpix map written\n";
*/
/*
 * Deallocate memories
 */

cout<<"Deallocating memories...\n";
clReleaseMemObject(recoDelaysBuffer);
clReleaseMemObject(intensityRBuffer);
clReleaseMemObject(intensityCBuffer);
clReleaseMemObject(voltsFlatBuffer);
clReleaseMemObject(beamCBuffer);
clReleaseMemObject(beamRBuffer);
clReleaseMemObject(pwrBuffer);
clReleaseMemObject(shiftedRBuffer);
clReleaseMemObject(shiftedCBuffer);
free(voltsFlat);
//free(volts);
//free(recoDelays);
free(intensity_data_r);
free(intensity_data_c);
free(shifted_data_r);
free(shifted_data_c);
free(beam_data_r);
free(beam_data_c);
free(beamPwr);
cout<<"Memories deallocated\n";

return maxPixIdx;
}

int reconstructXCorr(unsigned int dataType, vector<TGraph *>& cleanEvent, recoEnvData *clEnv,
                float *recoDelays, float *recoDelays_V, float *recoDelays_H,
                int nDir, string pol, const int *chanMask, char *filename)
{

cout<<"Entered reconstructXCorr method\n";
int nSamp;
int nAnt = (int)cleanEvent.size()/2; // Divide by 2 for only one polarization
int unmaskedNChan=0;
for(int ch=0; ch<2*nAnt; ch++) unmaskedNChan+=chanMask[ch];
cout<<"unmaskedNChan: "<<unmaskedNChan<<" nAnt: "<<nAnt<<endl;
float wInt;
if( dataType == 0 ) wInt = 0.5f; //AraSim event
else if( dataType == 1 ){ //real event
if( pol == "vpol" ) wInt = 0.4f;
else if ( pol == "hpol" ) wInt = 0.625f;
else { cerr<<"recoPolType undefined\n"; return -1; }
} else {
cerr<<"dataType undefined\n"; return -1; }

/*
 * Loading voltsFlat array
 */

double t, v;
float *voltsFlat;

if( pol == "vpol" ){

   /* Using the 1st vpol wf for nSamp. Should make sure all vpol channels have the same nSamp. FIXME */
   nSamp = cleanEvent[0]->GetN();
   voltsFlat = (float*)calloc(nAnt*nSamp, sizeof(float));
   if(voltsFlat == NULL){ cerr<<"Null pointer to voltsFlat\n"; return -1; }

   for(int ch=0; ch<nAnt; ch++){
      if( chanMask[ch] ){
         for(int s=0; s<nSamp; s++){
            cleanEvent[ch]->GetPoint(s,t,v);
	    /* Bartlett window applied in main analysis code */
	    //voltsFlat[nSamp*ch + s] = static_cast<float>(v * FFTtools::bartlettWindow(s, nSamp));
	    voltsFlat[nSamp*ch + s] = static_cast<float>(v);
         }
      } else {
         for(int s=0; s<nSamp; s++){
	    voltsFlat[nSamp*ch + s] = 0.f;
         }
      }
   }

} else if ( pol == "hpol" ){

   /* Using the 1st hpol wf for nSamp. Should make sure all hpol channels have the same nSamp. FIXME */
   nSamp = cleanEvent[nAnt]->GetN();
   voltsFlat = (float*)calloc(nAnt*nSamp, sizeof(float));
   if(voltsFlat == NULL){ cerr<<"Null pointer to voltsFlat\n"; return -1; }

   for(int ch=0; ch<nAnt; ch++){
       if( chanMask[ch+nAnt] ){
         for(int s=0; s<nSamp; s++){
            cleanEvent[ch+nAnt]->GetPoint(s,t,v);
            /* Bartlett window applied in main analysis code */
            //voltsFlat[nSamp*ch + s] = static_cast<float>(v * FFTtools::bartlettWindow(s, nSamp));
            voltsFlat[nSamp*ch + s] = static_cast<float>(v);
         }
      } else {
         for(int s=0; s<nSamp; s++){
            voltsFlat[nSamp*ch + s] = 0.f;
         }
      }
   }

} else {
   cerr<<"recoPolType not defined\n"; return -1;
}
cout<<"voltsFlat loaded\n";

/*
 * Preparation for OUT_OF_PLACE transforms
 */

int interlvHermOutputSize = 2*(1 + nSamp/2); //Hermitian layout
int interlvOutputSize     = 2*nSamp;         //Not Hermitian layout
int planarHermOutputSize = (1 + nSamp/2);    //Hermitian planar layout

float *intensity_data_r, *intensity_data_c;
intensity_data_r = (float*)calloc(nAnt*planarHermOutputSize, sizeof(float));
intensity_data_c = (float*)calloc(nAnt*planarHermOutputSize, sizeof(float));

/*
 * FFT library related declarations
 */
//clfftPlanHandle planHandle;
clfftDim fftDim = CLFFT_1D;
size_t clLengths[1]   = {nSamp};
size_t clInStride[1]  = {1};
size_t clOutStride[1] = {1};
size_t inDist  = nSamp;
size_t outDist   = planarHermOutputSize;
size_t batchSize = nAnt;

/*
 * Set up clFFT
 */
//clfftSetupData fftSetup;
//err = clfftInitSetupData(&fftSetup);
//err = clfftSetup(&fftSetup);

/*
 * Prepare plan
 */
cout<<"Preparing plan...\n";
int err;
err = clfftCreateDefaultPlan(&clEnv->planHandle, clEnv->context, fftDim, clLengths);
err = clfftSetPlanPrecision(clEnv->planHandle, CLFFT_SINGLE);
err = clfftSetLayout(clEnv->planHandle, CLFFT_REAL, CLFFT_HERMITIAN_PLANAR);
err = clfftSetPlanScale(clEnv->planHandle, CLFFT_FORWARD, 1.f);
err = clfftSetPlanBatchSize(clEnv->planHandle, batchSize);
err = clfftSetPlanInStride(clEnv->planHandle, fftDim, clInStride);
err = clfftSetPlanOutStride(clEnv->planHandle, fftDim, clOutStride);
err = clfftSetPlanDistance(clEnv->planHandle, inDist, outDist);
err = clfftSetResultLocation(clEnv->planHandle, CLFFT_OUTOFPLACE);
err = clfftBakePlan(clEnv->planHandle, 1, &clEnv->queue, NULL, NULL);
cout<<"Plan prepared\n";
/* The plan is now ready to be executed */
cl_mem voltsFlatBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*nSamp*sizeof(float), NULL, &err);
err = clEnqueueWriteBuffer(clEnv->queue, voltsFlatBuffer, CL_TRUE, 0, nAnt*nSamp*sizeof(float), voltsFlat, 0, NULL, NULL);

cl_mem intensityRBuffer= clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*planarHermOutputSize*sizeof(float),
                         NULL, &err);
cl_mem intensityCBuffer= clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*planarHermOutputSize*sizeof(float),
                         NULL, &err);
cl_mem outBuffers[2] = {intensityRBuffer, intensityCBuffer};

cout<<"Enqueueing FFT\n";
err = clfftEnqueueTransform(clEnv->planHandle, CLFFT_FORWARD, 1, &clEnv->queue, 0, NULL, NULL,
                            &voltsFlatBuffer, outBuffers, NULL);
err = clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, outBuffers[0], CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_r,
                    0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, outBuffers[1], CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_c,
                    0, NULL, NULL);
cout<<"FFT done\n";
/*
 * Clean up CLFFT
 */
cout<<"Destroying plan...\n";
err = clfftDestroyPlan(&clEnv->planHandle);
//clfftTeardown();
cout<<"Plan destroyed\n";

#ifdef bandpass
/*
 * Bandpass signals
 */

cout<<"Preparing bandpass filter..."<<endl;
float freqBin = 1e3 / (wInt * (float)nSamp); // wInt in ns. 1e3 for MHz
float lowFreq = 200.f;
float highFreq= 450.f;

clSetKernelArg(clEnv->bandPassFilter, 0, sizeof(cl_mem), &intensityRBuffer);
clSetKernelArg(clEnv->bandPassFilter, 1, sizeof(cl_mem), &intensityCBuffer);
clSetKernelArg(clEnv->bandPassFilter, 2, sizeof(float),  &freqBin);
clSetKernelArg(clEnv->bandPassFilter, 3, sizeof(float),  &lowFreq);
clSetKernelArg(clEnv->bandPassFilter, 4, sizeof(float),  &highFreq);

unsigned int dim = 2;
size_t bandPassWorkSize[2] = {(size_t)nAnt, (size_t)planarHermOutputSize};
clEnqueueNDRangeKernel(clEnv->queue, clEnv->bandPassFilter, dim, NULL, bandPassWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, intensityRBuffer, CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_r,
                    0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, intensityCBuffer, CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_c,
                    0, NULL, NULL);

cout<<"Bandpass filter done"<<endl;
#endif

/*
 * Cross-correlate wfs
 */

int nBaseline = nAnt*nAnt;
cout<<"nBaseline: "<<nBaseline<<endl;

float *xCorr_data_r = (float*)calloc(nBaseline*planarHermOutputSize, sizeof(float));
float *xCorr_data_c = (float*)calloc(nBaseline*planarHermOutputSize, sizeof(float));

cl_mem xCorrRBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                    sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_r, &err);

cl_mem xCorrCBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                    sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_c, &err);

clSetKernelArg(clEnv->xCorrWf, 0, sizeof(cl_mem), &xCorrRBuffer);
clSetKernelArg(clEnv->xCorrWf, 1, sizeof(cl_mem), &xCorrCBuffer);
clSetKernelArg(clEnv->xCorrWf, 2, sizeof(cl_mem), &intensityRBuffer);
clSetKernelArg(clEnv->xCorrWf, 3, sizeof(cl_mem), &intensityCBuffer);
clSetKernelArg(clEnv->xCorrWf, 4, sizeof(int),    &nAnt);

unsigned int workDim = 2;
size_t globalWorkSize[2] = {(size_t)nBaseline, (size_t)planarHermOutputSize};
clEnqueueNDRangeKernel(clEnv->queue, clEnv->xCorrWf, workDim, NULL, globalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, xCorrRBuffer, CL_TRUE, 0, sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_r, 0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, xCorrCBuffer, CL_TRUE, 0, sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_c, 0, NULL, NULL);
cout<<"Done xCorrWf\n";

/*
 * Inverse FFT cross-correlation
 */

inDist  = planarHermOutputSize;
outDist = nSamp;
batchSize = nBaseline;

/*
 * Prepare plan
 */

err = clfftCreateDefaultPlan(&clEnv->planHandle, clEnv->context, fftDim, clLengths);
err = clfftSetPlanPrecision(clEnv->planHandle, CLFFT_SINGLE);
err = clfftSetLayout(clEnv->planHandle, CLFFT_HERMITIAN_PLANAR, CLFFT_REAL);
err = clfftSetPlanScale(clEnv->planHandle, CLFFT_BACKWARD, 1.f/(float)nSamp);
err = clfftSetPlanBatchSize(clEnv->planHandle, batchSize);
err = clfftSetPlanInStride(clEnv->planHandle, fftDim, clInStride);
err = clfftSetPlanOutStride(clEnv->planHandle, fftDim, clOutStride);
err = clfftSetPlanDistance(clEnv->planHandle, inDist, outDist);
err = clfftSetResultLocation(clEnv->planHandle, CLFFT_OUTOFPLACE);
err = clfftBakePlan(clEnv->planHandle, 1, &clEnv->queue, NULL, NULL);

/* The plan is now ready to be executed */

cl_mem inBuffers[2] = {xCorrRBuffer, xCorrCBuffer};

float *xCorrTime = (float*)calloc(nBaseline*nSamp, sizeof(float));
cl_mem xCorrTimeBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nBaseline*nSamp*sizeof(float),
                           NULL, &err);
err = clfftEnqueueTransform(clEnv->planHandle, CLFFT_BACKWARD, 1, &clEnv->queue, 0, NULL, NULL, inBuffers, &xCorrTimeBuffer, NULL);
err = clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, xCorrTimeBuffer, CL_TRUE, 0, sizeof(float)*nBaseline*nSamp, xCorrTime, 0, NULL, NULL);

/*
 * Clean up FFT
 */

err = clfftDestroyPlan(&clEnv->planHandle);
cout<<"Done inverse FFT\n";

/*
 * Compute cross-correlation coefficients Cij in each direction
 */

/* Compute the square root of each channel's total wf power */

float *sqrtWfPwr = (float*)calloc(nAnt, sizeof(float));
float pwr=0.f;
for(int ant=0; ant<nAnt; ant++){
   pwr = 0.;
   for(int s=0; s<nSamp; s++){
      pwr += (voltsFlat[ant*nSamp + s] * voltsFlat[ant*nSamp + s]);
   }
   sqrtWfPwr[ant] = sqrt(pwr);
   //sqrtWfPwr[ant] = 1.f;
}
cout<<"Done sqrtWfPwr\n";

cl_mem sqrtWfPwrBuffer  = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt,
                          sqrtWfPwr, &err);
cl_mem recoDelaysBuffer;
if(pol == "vpol" ){
recoDelaysBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt*nDir,
                                  recoDelays_V, &err);
} else if (pol == "hpol" ) {
recoDelaysBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt*nDir,
                                  recoDelays_H, &err);
} else {
cerr<<"recoPolType not defined\n"; return -1;
}
cout<<"recoDelaysBuffer created\n";

float *Cij = (float*)calloc(nDir*nBaseline, sizeof(float));

cl_mem CijBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nDir*nBaseline,
                   Cij, &err);

clSetKernelArg(clEnv->computeXCorrCoef, 0, sizeof(cl_mem), &CijBuffer);
clSetKernelArg(clEnv->computeXCorrCoef, 1, sizeof(cl_mem), &xCorrTimeBuffer);
clSetKernelArg(clEnv->computeXCorrCoef, 2, sizeof(cl_mem), &recoDelaysBuffer);
clSetKernelArg(clEnv->computeXCorrCoef, 3, sizeof(cl_mem), &sqrtWfPwrBuffer);
clSetKernelArg(clEnv->computeXCorrCoef, 4, sizeof(float),  &wInt);
clSetKernelArg(clEnv->computeXCorrCoef, 5, sizeof(int),    &nAnt);
clSetKernelArg(clEnv->computeXCorrCoef, 6, sizeof(int),    &nSamp);

workDim = 2;
size_t CijGlobalWorkSize[2] = {(size_t)nDir, (size_t)nBaseline};

err = clEnqueueNDRangeKernel(clEnv->queue, clEnv->computeXCorrCoef, workDim, NULL, CijGlobalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, CijBuffer, CL_TRUE, 0, sizeof(float)*nDir*nBaseline, Cij, 0, NULL, NULL);

cout<<"Done computeXCorrCoef\n";

/*
 * Sum Cij's of all baselines in each reco direction to obtain coherence M(r-hat)
 */

float *M = (float*)calloc(nDir, sizeof(float));
cl_mem MBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nDir,
                  M, &err);

clSetKernelArg(clEnv->computeCoherence, 0, sizeof(cl_mem), &MBuffer);
clSetKernelArg(clEnv->computeCoherence, 1, sizeof(cl_mem), &CijBuffer);
clSetKernelArg(clEnv->computeCoherence, 2, sizeof(int),    &nBaseline);

workDim = 1;
size_t MGlobalWorkSize = nDir;

clEnqueueNDRangeKernel(clEnv->queue, clEnv->computeCoherence, workDim, NULL, &MGlobalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, MBuffer, CL_TRUE, 0, sizeof(float)*nDir, M, 0, NULL, NULL);

cout<<"Done computeCoherence\n";

/*
 * Write FITS file
 */
cout<<"Creating Healpix map and writing to FITS....\n";
arr<float> MArr = arr<float>(&M[0], (size_t)nDir);
Healpix_Map<float> skyMap = Healpix_Map<float>(MArr, HEALPIX_ORDERING);

fitshandle fitsOut;
//char filename[] = "testXCorrSkyMap.fits";
remove(filename);
fitsOut.create(filename);

write_Healpix_map_to_fits(fitsOut, skyMap, PLANCK_FLOAT32);
cout<<"Healpix map written\n";

/*
 * Deallocate memories
 */

cout<<"Deallocating memories...\n";
clReleaseMemObject(recoDelaysBuffer);
clReleaseMemObject(intensityRBuffer);
clReleaseMemObject(intensityCBuffer);
clReleaseMemObject(voltsFlatBuffer);
clReleaseMemObject(xCorrCBuffer);
clReleaseMemObject(xCorrRBuffer);
clReleaseMemObject(xCorrTimeBuffer);
clReleaseMemObject(sqrtWfPwrBuffer);
clReleaseMemObject(CijBuffer);
clReleaseMemObject(MBuffer);
free(voltsFlat);
//free(volts);
//free(recoDelays);
free(intensity_data_r);
free(intensity_data_c);
free(xCorr_data_r);
free(xCorr_data_c);
free(xCorrTime);
free(sqrtWfPwr);
free(Cij);
free(M);
cout<<"Memories deallocated\n";

return 0;
}


int reconstructXCorrEnvelope(unsigned int dataType, vector<TGraph *>& cleanEvent, recoEnvData *clEnv,
                float *recoDelays, float *recoDelays_V, float *recoDelays_H,
                int nDir, string pol, const int *chanMask, char *filename)
{

cout<<"Entered reconstructXCorrEnvelope method\n";
int nSamp;
int nAnt = (int)cleanEvent.size()/2; // Divide by 2 for only one polarization
int unmaskedNChan=0;
for(int ch=0; ch<2*nAnt; ch++) unmaskedNChan+=chanMask[ch];
cout<<"unmaskedNChan: "<<unmaskedNChan<<" nAnt: "<<nAnt<<endl;
float wInt;
if( dataType == 0 ) wInt = 0.5f; //AraSim event
else if( dataType == 1 ){ //real event
if( pol == "vpol" ) wInt = 0.4f;
else if ( pol == "hpol" ) wInt = 0.625f;
else { cerr<<"recoPolType undefined\n"; return -1; }
} else {
cerr<<"dataType undefined\n"; return -1; }

/*
 * Loading voltsFlat array
 */

double t, v;
float *voltsFlat;

if( pol == "vpol" ){

   /* Using the 1st vpol wf for nSamp. Should make sure all vpol channels have the same nSamp. FIXME */
   nSamp = cleanEvent[0]->GetN();
   voltsFlat = (float*)calloc(nAnt*nSamp, sizeof(float));
   if(voltsFlat == NULL){ cerr<<"Null pointer to voltsFlat\n"; return -1; }

   for(int ch=0; ch<nAnt; ch++){
      if( chanMask[ch] ){
         for(int s=0; s<nSamp; s++){
            cleanEvent[ch]->GetPoint(s,t,v);
	    /* Bartlett window applied in main analysis code */
	    //voltsFlat[nSamp*ch + s] = static_cast<float>(v * FFTtools::bartlettWindow(s, nSamp));
	    voltsFlat[nSamp*ch + s] = static_cast<float>(v);
         }
      } else {
         for(int s=0; s<nSamp; s++){
	    voltsFlat[nSamp*ch + s] = 0.f;
         }
      }
   }

} else if ( pol == "hpol" ){

   /* Using the 1st hpol wf for nSamp. Should make sure all hpol channels have the same nSamp. FIXME */
   nSamp = cleanEvent[nAnt]->GetN();
   voltsFlat = (float*)calloc(nAnt*nSamp, sizeof(float));
   if(voltsFlat == NULL){ cerr<<"Null pointer to voltsFlat\n"; return -1; }

   for(int ch=0; ch<nAnt; ch++){
       if( chanMask[ch+nAnt] ){
         for(int s=0; s<nSamp; s++){
            cleanEvent[ch+nAnt]->GetPoint(s,t,v);
            /* Bartlett window applied in main analysis code */
            //voltsFlat[nSamp*ch + s] = static_cast<float>(v * FFTtools::bartlettWindow(s, nSamp));
            voltsFlat[nSamp*ch + s] = static_cast<float>(v);
         }
      } else {
         for(int s=0; s<nSamp; s++){
            voltsFlat[nSamp*ch + s] = 0.f;
         }
      }
   }

} else {
   cerr<<"recoPolType not defined\n"; return -1;
}
cout<<"voltsFlat loaded\n";

/*
 * Preparation for OUT_OF_PLACE transforms
 */

int interlvHermOutputSize = 2*(1 + nSamp/2); //Hermitian layout
int interlvOutputSize     = 2*nSamp;         //Not Hermitian layout
int planarHermOutputSize = (1 + nSamp/2);    //Hermitian planar layout

float *intensity_data_r, *intensity_data_c;
intensity_data_r = (float*)calloc(nAnt*planarHermOutputSize, sizeof(float));
intensity_data_c = (float*)calloc(nAnt*planarHermOutputSize, sizeof(float));

/*
 * FFT library related declarations
 */
//clfftPlanHandle planHandle;
clfftDim fftDim = CLFFT_1D;
size_t clLengths[1]   = {nSamp};
size_t clInStride[1]  = {1};
size_t clOutStride[1] = {1};
size_t inDist  = nSamp;
size_t outDist   = planarHermOutputSize;
size_t batchSize = nAnt;

/*
 * Set up clFFT
 */
//clfftSetupData fftSetup;
//err = clfftInitSetupData(&fftSetup);
//err = clfftSetup(&fftSetup);

/*
 * Prepare plan
 */
cout<<"Preparing plan...\n";
int err;
err = clfftCreateDefaultPlan(&clEnv->planHandle, clEnv->context, fftDim, clLengths);
err = clfftSetPlanPrecision(clEnv->planHandle, CLFFT_SINGLE);
err = clfftSetLayout(clEnv->planHandle, CLFFT_REAL, CLFFT_HERMITIAN_PLANAR);
err = clfftSetPlanScale(clEnv->planHandle, CLFFT_FORWARD, 1.f);
err = clfftSetPlanBatchSize(clEnv->planHandle, batchSize);
err = clfftSetPlanInStride(clEnv->planHandle, fftDim, clInStride);
err = clfftSetPlanOutStride(clEnv->planHandle, fftDim, clOutStride);
err = clfftSetPlanDistance(clEnv->planHandle, inDist, outDist);
err = clfftSetResultLocation(clEnv->planHandle, CLFFT_OUTOFPLACE);
err = clfftBakePlan(clEnv->planHandle, 1, &clEnv->queue, NULL, NULL);
cout<<"Plan prepared\n";
/* The plan is now ready to be executed */
cl_mem voltsFlatBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*nSamp*sizeof(float), NULL, &err);
err = clEnqueueWriteBuffer(clEnv->queue, voltsFlatBuffer, CL_TRUE, 0, nAnt*nSamp*sizeof(float), voltsFlat, 0, NULL, NULL);

cl_mem intensityRBuffer= clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*planarHermOutputSize*sizeof(float),
                         NULL, &err);
cl_mem intensityCBuffer= clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*planarHermOutputSize*sizeof(float),
                         NULL, &err);
cl_mem outBuffers[2] = {intensityRBuffer, intensityCBuffer};

cout<<"Enqueueing FFT\n";
err = clfftEnqueueTransform(clEnv->planHandle, CLFFT_FORWARD, 1, &clEnv->queue, 0, NULL, NULL,
                            &voltsFlatBuffer, outBuffers, NULL);
err = clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, outBuffers[0], CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_r,
                    0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, outBuffers[1], CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_c,
                    0, NULL, NULL);
cout<<"FFT done\n";
/*
 * Clean up CLFFT
 */
cout<<"Destroying plan...\n";
err = clfftDestroyPlan(&clEnv->planHandle);
//clfftTeardown();
cout<<"Plan destroyed\n";

#ifdef bandpass
/*
 * Bandpass signals
 */

cout<<"Preparing bandpass filter..."<<endl;
float freqBin = 1e3 / (wInt * (float)nSamp); // wInt in ns. 1e3 for MHz
float lowFreq = 200.f;
float highFreq= 450.f;

clSetKernelArg(clEnv->bandPassFilter, 0, sizeof(cl_mem), &intensityRBuffer);
clSetKernelArg(clEnv->bandPassFilter, 1, sizeof(cl_mem), &intensityCBuffer);
clSetKernelArg(clEnv->bandPassFilter, 2, sizeof(float),  &freqBin);
clSetKernelArg(clEnv->bandPassFilter, 3, sizeof(float),  &lowFreq);
clSetKernelArg(clEnv->bandPassFilter, 4, sizeof(float),  &highFreq);

unsigned int dim = 2;
size_t bandPassWorkSize[2] = {(size_t)nAnt, (size_t)planarHermOutputSize};
clEnqueueNDRangeKernel(clEnv->queue, clEnv->bandPassFilter, dim, NULL, bandPassWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, intensityRBuffer, CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_r,
                    0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, intensityCBuffer, CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_c,
                    0, NULL, NULL);

cout<<"Bandpass filter done"<<endl;
#endif

/*
 * Cross-correlate wfs
 */

int nBaseline = nAnt*nAnt;
cout<<"nBaseline: "<<nBaseline<<endl;

float *xCorr_data_r = (float*)calloc(nBaseline*planarHermOutputSize, sizeof(float));
float *xCorr_data_c = (float*)calloc(nBaseline*planarHermOutputSize, sizeof(float));

cl_mem xCorrRBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                    sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_r, &err);

cl_mem xCorrCBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                    sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_c, &err);

clSetKernelArg(clEnv->xCorrWf, 0, sizeof(cl_mem), &xCorrRBuffer);
clSetKernelArg(clEnv->xCorrWf, 1, sizeof(cl_mem), &xCorrCBuffer);
clSetKernelArg(clEnv->xCorrWf, 2, sizeof(cl_mem), &intensityRBuffer);
clSetKernelArg(clEnv->xCorrWf, 3, sizeof(cl_mem), &intensityCBuffer);
clSetKernelArg(clEnv->xCorrWf, 4, sizeof(int),    &nAnt);

unsigned int workDim = 2;
size_t globalWorkSize[2] = {(size_t)nBaseline, (size_t)planarHermOutputSize};
clEnqueueNDRangeKernel(clEnv->queue, clEnv->xCorrWf, workDim, NULL, globalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, xCorrRBuffer, CL_TRUE, 0, sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_r, 0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, xCorrCBuffer, CL_TRUE, 0, sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_c, 0, NULL, NULL);
cout<<"Done xCorrWf\n";

/*
 * Inverse FFT cross-correlation
 */

inDist  = planarHermOutputSize;
outDist = nSamp;
batchSize = nBaseline;

/*
 * Prepare plan
 */

err = clfftCreateDefaultPlan(&clEnv->planHandle, clEnv->context, fftDim, clLengths);
err = clfftSetPlanPrecision(clEnv->planHandle, CLFFT_SINGLE);
err = clfftSetLayout(clEnv->planHandle, CLFFT_HERMITIAN_PLANAR, CLFFT_REAL);
err = clfftSetPlanScale(clEnv->planHandle, CLFFT_BACKWARD, 1.f/(float)nSamp);
err = clfftSetPlanBatchSize(clEnv->planHandle, batchSize);
err = clfftSetPlanInStride(clEnv->planHandle, fftDim, clInStride);
err = clfftSetPlanOutStride(clEnv->planHandle, fftDim, clOutStride);
err = clfftSetPlanDistance(clEnv->planHandle, inDist, outDist);
err = clfftSetResultLocation(clEnv->planHandle, CLFFT_OUTOFPLACE);
err = clfftBakePlan(clEnv->planHandle, 1, &clEnv->queue, NULL, NULL);

/* The plan is now ready to be executed */

cl_mem inBuffers[2] = {xCorrRBuffer, xCorrCBuffer};

float *xCorrTime = (float*)calloc(nBaseline*nSamp, sizeof(float));
cl_mem xCorrTimeBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nBaseline*nSamp*sizeof(float),
                           NULL, &err);
err = clfftEnqueueTransform(clEnv->planHandle, CLFFT_BACKWARD, 1, &clEnv->queue, 0, NULL, NULL, inBuffers, &xCorrTimeBuffer, NULL);
err = clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, xCorrTimeBuffer, CL_TRUE, 0, sizeof(float)*nBaseline*nSamp, xCorrTime, 0, NULL, NULL);

/*
 * Clean up FFT
 */

err = clfftDestroyPlan(&clEnv->planHandle);
clfftTeardown();
cout<<"Done inverse FFT\n";

/*
 * Get Hilbert transform of XCorr function
 */

//TCanvas cvs("cvs","cva",800,600);
//cvs.Divide(1,2);


float dt[nSamp];
float xCorrValue[nSamp];
double t_temp, v_temp;

for(int baseline=0; baseline<nBaseline; baseline++){

   for(int s=0; s<nSamp; s++){

   dt[s] = wInt*s;
   xCorrValue[s] = xCorrTime[nSamp*baseline + s];

   }

   TGraph *xCorrGraph = new TGraph(nSamp, dt, xCorrValue);

   //cvs.cd(1);
   //xCorrGraph->Draw("AL");

   TGraph* envelope = FFTtools::getHilbertEnvelope( xCorrGraph );

   //cvs.cd(2);
   //envelope->Draw("AL");

   //cvs.SaveAs("xCorrEnvelope.C");

   for(int s=0; s<nSamp; s++){

   envelope->GetPoint(s,t_temp,v_temp);
   xCorrTime[nSamp*baseline + s] = static_cast<float>(v_temp);

   }

  delete xCorrGraph;
  delete envelope;

}

cl_mem xCorrEnvBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                       sizeof(float)*nBaseline*nSamp,
                                       xCorrTime, &err);


/*
 * Compute cross-correlation coefficients Cij in each direction
 */

/* Compute the square root of each channel's total wf power */

float *sqrtWfPwr = (float*)calloc(nAnt, sizeof(float));
float pwr=0.f;
for(int ant=0; ant<nAnt; ant++){
   pwr = 0.;
   for(int s=0; s<nSamp; s++){
      pwr += (voltsFlat[ant*nSamp + s] * voltsFlat[ant*nSamp + s]);
   }
   sqrtWfPwr[ant] = sqrt(pwr);
   //sqrtWfPwr[ant] = 1.f;
}
cout<<"Done sqrtWfPwr\n";

cl_mem sqrtWfPwrBuffer  = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt,
                          sqrtWfPwr, &err);
cl_mem recoDelaysBuffer;
if(pol == "vpol" ){
recoDelaysBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt*nDir,
                                  recoDelays_V, &err);
} else if (pol == "hpol" ) {
recoDelaysBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt*nDir,
                                  recoDelays_H, &err);
} else {
cerr<<"recoPolType not defined\n"; return -1;
}
cout<<"recoDelaysBuffer created\n";

float *Cij = (float*)calloc(nDir*nBaseline, sizeof(float));

cl_mem CijBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nDir*nBaseline,
                   Cij, &err);

clSetKernelArg(clEnv->computeXCorrCoef, 0, sizeof(cl_mem), &CijBuffer);
clSetKernelArg(clEnv->computeXCorrCoef, 1, sizeof(cl_mem), &xCorrEnvBuffer);
clSetKernelArg(clEnv->computeXCorrCoef, 2, sizeof(cl_mem), &recoDelaysBuffer);
clSetKernelArg(clEnv->computeXCorrCoef, 3, sizeof(cl_mem), &sqrtWfPwrBuffer);
clSetKernelArg(clEnv->computeXCorrCoef, 4, sizeof(float),  &wInt);
clSetKernelArg(clEnv->computeXCorrCoef, 5, sizeof(int),    &nAnt);
clSetKernelArg(clEnv->computeXCorrCoef, 6, sizeof(int),    &nSamp);

workDim = 2;
size_t CijGlobalWorkSize[2] = {(size_t)nDir, (size_t)nBaseline};

err = clEnqueueNDRangeKernel(clEnv->queue, clEnv->computeXCorrCoef, workDim, NULL, CijGlobalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, CijBuffer, CL_TRUE, 0, sizeof(float)*nDir*nBaseline, Cij, 0, NULL, NULL);

cout<<"Done computeXCorrCoef\n";

/*
 * Sum Cij's of all baselines in each reco direction to obtain coherence M(r-hat)
 */

float *M = (float*)calloc(nDir, sizeof(float));
cl_mem MBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nDir,
                  M, &err);

clSetKernelArg(clEnv->computeCoherence, 0, sizeof(cl_mem), &MBuffer);
clSetKernelArg(clEnv->computeCoherence, 1, sizeof(cl_mem), &CijBuffer);
clSetKernelArg(clEnv->computeCoherence, 2, sizeof(int),    &nBaseline);

workDim = 1;
size_t MGlobalWorkSize = nDir;

clEnqueueNDRangeKernel(clEnv->queue, clEnv->computeCoherence, workDim, NULL, &MGlobalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, MBuffer, CL_TRUE, 0, sizeof(float)*nDir, M, 0, NULL, NULL);

cout<<"Done computeCoherence\n";

/*
 * Write FITS file
 */
cout<<"Creating Healpix map and writing to FITS....\n";
arr<float> MArr = arr<float>(&M[0], (size_t)nDir);
Healpix_Map<float> skyMap = Healpix_Map<float>(MArr, HEALPIX_ORDERING);

fitshandle fitsOut;
//char filename[] = "testXCorrSkyMap.fits";
remove(filename);
fitsOut.create(filename);

write_Healpix_map_to_fits(fitsOut, skyMap, PLANCK_FLOAT32);
cout<<"Healpix map written\n";

/*
 * Deallocate memories
 */

cout<<"Deallocating memories...\n";
clReleaseMemObject(recoDelaysBuffer);
clReleaseMemObject(intensityRBuffer);
clReleaseMemObject(intensityCBuffer);
clReleaseMemObject(voltsFlatBuffer);
clReleaseMemObject(xCorrCBuffer);
clReleaseMemObject(xCorrRBuffer);
clReleaseMemObject(xCorrTimeBuffer);
clReleaseMemObject(xCorrEnvBuffer);
clReleaseMemObject(sqrtWfPwrBuffer);
clReleaseMemObject(CijBuffer);
clReleaseMemObject(MBuffer);
free(voltsFlat);
//free(volts);
//free(recoDelays);
free(intensity_data_r);
free(intensity_data_c);
free(xCorr_data_r);
free(xCorr_data_c);
free(xCorrTime);
free(sqrtWfPwr);
free(Cij);
free(M);
cout<<"Memories deallocated\n";

return 0;
}

int reconstructXCorrGetMaxPix(unsigned int dataType, vector<TGraph *>& cleanEvent, recoEnvData *clEnv,
                float *recoDelays, float *recoDelays_V, float *recoDelays_H,
                int nDir, string pol, const int *chanMask, recoData *summary)
{

cout<<"Entered reconstructXCorrGetMaxPix method\n";
int nSamp;
int nAnt = (int)cleanEvent.size()/2; // Divide by 2 for only one polarization
int unmaskedNChan=0;
for(int ch=0; ch<2*nAnt; ch++) unmaskedNChan+=chanMask[ch];
cout<<"unmaskedNChan: "<<unmaskedNChan<<" nAnt: "<<nAnt<<endl;
float wInt;
if( dataType == 0 ) wInt = 0.5f; //AraSim event
else if( dataType == 1 ){ //real event
if( pol == "vpol" ) wInt = 0.4f;
else if ( pol == "hpol" ) wInt = 0.625f;
else { cerr<<"recoPolType undefined\n"; return -1; }
} else {
cerr<<"dataType undefined\n"; return -1; }

/*
 * Loading voltsFlat array
 */

double t, v;
float *voltsFlat;

if( pol == "vpol" ){

   /* Using the 1st vpol wf for nSamp. Should make sure all vpol channels have the same nSamp. FIXME */
   nSamp = cleanEvent[0]->GetN();
   voltsFlat = (float*)calloc(nAnt*nSamp, sizeof(float));
   if(voltsFlat == NULL){ cerr<<"Null pointer to voltsFlat\n"; return -1; }

   for(int ch=0; ch<nAnt; ch++){
       if( chanMask[ch] ){
         for(int s=0; s<nSamp; s++){
            cleanEvent[ch]->GetPoint(s,t,v);
            /* Bartlett window applied in main analysis code */
            //voltsFlat[nSamp*ch + s] = static_cast<float>(v * FFTtools::bartlettWindow(s, nSamp));
            voltsFlat[nSamp*ch + s] = static_cast<float>(v);
         }
      } else {
         for(int s=0; s<nSamp; s++){
            voltsFlat[nSamp*ch + s] = 0.f;
         }
      }
   }

} else if ( pol == "hpol" ){

   /* Using the 1st hpol wf for nSamp. Should make sure all hpol channels have the same nSamp. FIXME */
   nSamp = cleanEvent[nAnt]->GetN();
   voltsFlat = (float*)calloc(nAnt*nSamp, sizeof(float));
   if(voltsFlat == NULL){ cerr<<"Null pointer to voltsFlat\n"; return -1; }

   for(int ch=0; ch<nAnt; ch++){
       if( chanMask[ch+nAnt] ){
         for(int s=0; s<nSamp; s++){
            cleanEvent[ch+nAnt]->GetPoint(s,t,v);
            /* Bartlett window applied in main analysis code */
            //voltsFlat[nSamp*ch + s] = static_cast<float>(v * FFTtools::bartlettWindow(s, nSamp));
            voltsFlat[nSamp*ch + s] = static_cast<float>(v);
         }
      } else {
         for(int s=0; s<nSamp; s++){
            voltsFlat[nSamp*ch + s] = 0.f;
         }
      }
   }

} else {
   cerr<<"recoPolType not defined\n"; return -1;
}
cout<<"voltsFlat loaded\n";

/*
 * Preparation for OUT_OF_PLACE transforms
 */

int interlvHermOutputSize = 2*(1 + nSamp/2); //Hermitian layout
int interlvOutputSize     = 2*nSamp;         //Not Hermitian layout
int planarHermOutputSize = (1 + nSamp/2);    //Hermitian planar layout

float *intensity_data_r, *intensity_data_c;
intensity_data_r = (float*)calloc(nAnt*planarHermOutputSize, sizeof(float));
intensity_data_c = (float*)calloc(nAnt*planarHermOutputSize, sizeof(float));

/*
 * FFT library related declarations
 */
//clfftPlanHandle planHandle;
clfftDim fftDim = CLFFT_1D;
size_t clLengths[1]   = {nSamp};
size_t clInStride[1]  = {1};
size_t clOutStride[1] = {1};
size_t inDist  = nSamp;
size_t outDist   = planarHermOutputSize;
size_t batchSize = nAnt;

/*
 * Set up clFFT
 */
//clfftSetupData fftSetup;
//err = clfftInitSetupData(&fftSetup);
//err = clfftSetup(&fftSetup);

/*
 * Prepare plan
 */
cout<<"Preparing plan...\n";
int err;
err = clfftCreateDefaultPlan(&clEnv->planHandle, clEnv->context, fftDim, clLengths);
err = clfftSetPlanPrecision(clEnv->planHandle, CLFFT_SINGLE);
err = clfftSetLayout(clEnv->planHandle, CLFFT_REAL, CLFFT_HERMITIAN_PLANAR);
err = clfftSetPlanScale(clEnv->planHandle, CLFFT_FORWARD, 1.f);
err = clfftSetPlanBatchSize(clEnv->planHandle, batchSize);
err = clfftSetPlanInStride(clEnv->planHandle, fftDim, clInStride);
err = clfftSetPlanOutStride(clEnv->planHandle, fftDim, clOutStride);
err = clfftSetPlanDistance(clEnv->planHandle, inDist, outDist);
err = clfftSetResultLocation(clEnv->planHandle, CLFFT_OUTOFPLACE);
err = clfftBakePlan(clEnv->planHandle, 1, &clEnv->queue, NULL, NULL);
cout<<"Plan prepared\n";
/* The plan is now ready to be executed */
cl_mem voltsFlatBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*nSamp*sizeof(float), NULL, &err);
err = clEnqueueWriteBuffer(clEnv->queue, voltsFlatBuffer, CL_TRUE, 0, nAnt*nSamp*sizeof(float), voltsFlat, 0, NULL, NULL);

cl_mem intensityRBuffer= clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*planarHermOutputSize*sizeof(float),
                         NULL, &err);
cl_mem intensityCBuffer= clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*planarHermOutputSize*sizeof(float),
                         NULL, &err);
cl_mem outBuffers[2] = {intensityRBuffer, intensityCBuffer};

cout<<"Enqueueing FFT\n";
err = clfftEnqueueTransform(clEnv->planHandle, CLFFT_FORWARD, 1, &clEnv->queue, 0, NULL, NULL,
                            &voltsFlatBuffer, outBuffers, NULL);
err = clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, outBuffers[0], CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_r,
                    0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, outBuffers[1], CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_c,
                    0, NULL, NULL);
cout<<"FFT done\n";
/*
 * Clean up CLFFT
 */
cout<<"Destroying plan...\n";
err = clfftDestroyPlan(&clEnv->planHandle);
//clfftTeardown();
cout<<"Plan destroyed\n";


#ifdef bandpass
/*
 * Bandpass signals
 */
cout<<"Preparing bandpass filter..."<<endl;
float freqBin = 1e3 / (wInt * (float)nSamp); // wInt in ns. 1e3 for MHz
float lowFreq = 200.f;
float highFreq= 450.f;

clSetKernelArg(clEnv->bandPassFilter, 0, sizeof(cl_mem), &intensityRBuffer);
clSetKernelArg(clEnv->bandPassFilter, 1, sizeof(cl_mem), &intensityCBuffer);
clSetKernelArg(clEnv->bandPassFilter, 2, sizeof(float),  &freqBin);
clSetKernelArg(clEnv->bandPassFilter, 3, sizeof(float),  &lowFreq);
clSetKernelArg(clEnv->bandPassFilter, 4, sizeof(float),  &highFreq);

unsigned int dim = 2;
size_t bandPassWorkSize[2] = {(size_t)nAnt, (size_t)planarHermOutputSize};
clEnqueueNDRangeKernel(clEnv->queue, clEnv->bandPassFilter, dim, NULL, bandPassWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, intensityRBuffer, CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_r,
                    0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, intensityCBuffer, CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_c,
                    0, NULL, NULL);

cout<<"Bandpass filter done"<<endl;
#endif

/*
 * Cross-correlate wfs
 */

int nBaseline = nAnt*nAnt;
cout<<"nBaseline: "<<nBaseline<<endl;

float *xCorr_data_r = (float*)calloc(nBaseline*planarHermOutputSize, sizeof(float));
float *xCorr_data_c = (float*)calloc(nBaseline*planarHermOutputSize, sizeof(float));

cl_mem xCorrRBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                    sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_r, &err);

cl_mem xCorrCBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                    sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_c, &err);

clSetKernelArg(clEnv->xCorrWf, 0, sizeof(cl_mem), &xCorrRBuffer);
clSetKernelArg(clEnv->xCorrWf, 1, sizeof(cl_mem), &xCorrCBuffer);
clSetKernelArg(clEnv->xCorrWf, 2, sizeof(cl_mem), &intensityRBuffer);
clSetKernelArg(clEnv->xCorrWf, 3, sizeof(cl_mem), &intensityCBuffer);
clSetKernelArg(clEnv->xCorrWf, 4, sizeof(int),    &nAnt);

unsigned int workDim = 2;
size_t globalWorkSize[2] = {(size_t)nBaseline, (size_t)planarHermOutputSize};
clEnqueueNDRangeKernel(clEnv->queue, clEnv->xCorrWf, workDim, NULL, globalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, xCorrRBuffer, CL_TRUE, 0, sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_r, 0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, xCorrCBuffer, CL_TRUE, 0, sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_c, 0, NULL, NULL);
cout<<"Done xCorrWf\n";

/*
 * Inverse FFT cross-correlation
 */

inDist  = planarHermOutputSize;
outDist = nSamp;
batchSize = nBaseline;

/*
 * Prepare plan
 */

err = clfftCreateDefaultPlan(&clEnv->planHandle, clEnv->context, fftDim, clLengths);
err = clfftSetPlanPrecision(clEnv->planHandle, CLFFT_SINGLE);
err = clfftSetLayout(clEnv->planHandle, CLFFT_HERMITIAN_PLANAR, CLFFT_REAL);
err = clfftSetPlanScale(clEnv->planHandle, CLFFT_BACKWARD, 1.f/(float)nSamp);
err = clfftSetPlanBatchSize(clEnv->planHandle, batchSize);
err = clfftSetPlanInStride(clEnv->planHandle, fftDim, clInStride);
err = clfftSetPlanOutStride(clEnv->planHandle, fftDim, clOutStride);
err = clfftSetPlanDistance(clEnv->planHandle, inDist, outDist);
err = clfftSetResultLocation(clEnv->planHandle, CLFFT_OUTOFPLACE);
err = clfftBakePlan(clEnv->planHandle, 1, &clEnv->queue, NULL, NULL);

/* The plan is now ready to be executed */

cl_mem inBuffers[2] = {xCorrRBuffer, xCorrCBuffer};

float *xCorrTime = (float*)calloc(nBaseline*nSamp, sizeof(float));
cl_mem xCorrTimeBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nBaseline*nSamp*sizeof(float),
                           NULL, &err);
err = clfftEnqueueTransform(clEnv->planHandle, CLFFT_BACKWARD, 1, &clEnv->queue, 0, NULL, NULL, inBuffers, &xCorrTimeBuffer, NULL);
err = clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, xCorrTimeBuffer, CL_TRUE, 0, sizeof(float)*nBaseline*nSamp, xCorrTime, 0, NULL, NULL);

/*
 * Clean up FFT
 */

err = clfftDestroyPlan(&clEnv->planHandle);
cout<<"Done inverse FFT\n";

/*
 * Compute cross-correlation coefficients Cij in each direction
 */

/* Compute the square root of each channel's total wf power */

float *sqrtWfPwr = (float*)calloc(nAnt, sizeof(float));
float pwr=0.f;
for(int ant=0; ant<nAnt; ant++){
   pwr = 0.;
   for(int s=0; s<nSamp; s++){
      pwr += (voltsFlat[ant*nSamp + s] * voltsFlat[ant*nSamp + s]);
   }
   sqrtWfPwr[ant] = sqrt(pwr);
   //sqrtWfPwr[ant] = 1.f;
}
cout<<"Done sqrtWfPwr\n";

cl_mem sqrtWfPwrBuffer  = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt,
                          sqrtWfPwr, &err);
cl_mem recoDelaysBuffer;
if(pol == "vpol" ){
recoDelaysBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt*nDir,
                                  recoDelays_V, &err);
} else if (pol == "hpol" ) {
recoDelaysBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt*nDir,
                                  recoDelays_H, &err);
} else {
cerr<<"recoPolType not defined\n"; return -1;
}
cout<<"recoDelaysBuffer created\n";

float *Cij = (float*)calloc(nDir*nBaseline, sizeof(float));

cl_mem CijBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nDir*nBaseline,
                   Cij, &err);

clSetKernelArg(clEnv->computeXCorrCoef, 0, sizeof(cl_mem), &CijBuffer);
clSetKernelArg(clEnv->computeXCorrCoef, 1, sizeof(cl_mem), &xCorrTimeBuffer);
clSetKernelArg(clEnv->computeXCorrCoef, 2, sizeof(cl_mem), &recoDelaysBuffer);
clSetKernelArg(clEnv->computeXCorrCoef, 3, sizeof(cl_mem), &sqrtWfPwrBuffer);
clSetKernelArg(clEnv->computeXCorrCoef, 4, sizeof(float),  &wInt);
clSetKernelArg(clEnv->computeXCorrCoef, 5, sizeof(int),    &nAnt);
clSetKernelArg(clEnv->computeXCorrCoef, 6, sizeof(int),    &nSamp);

workDim = 2;
size_t CijGlobalWorkSize[2] = {(size_t)nDir, (size_t)nBaseline};

err = clEnqueueNDRangeKernel(clEnv->queue, clEnv->computeXCorrCoef, workDim, NULL, CijGlobalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, CijBuffer, CL_TRUE, 0, sizeof(float)*nDir*nBaseline, Cij, 0, NULL, NULL);

cout<<"Done computeXCorrCoef\n";

/*
 * Sum Cij's of all baselines in each reco direction to obtain coherence M(r-hat)
 */

float *M = (float*)calloc(nDir, sizeof(float));
cl_mem MBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nDir,
                  M, &err);

clSetKernelArg(clEnv->computeCoherence, 0, sizeof(cl_mem), &MBuffer);
clSetKernelArg(clEnv->computeCoherence, 1, sizeof(cl_mem), &CijBuffer);
clSetKernelArg(clEnv->computeCoherence, 2, sizeof(int),    &nBaseline);

workDim = 1;
size_t MGlobalWorkSize = nDir;

clEnqueueNDRangeKernel(clEnv->queue, clEnv->computeCoherence, workDim, NULL, &MGlobalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, MBuffer, CL_TRUE, 0, sizeof(float)*nDir, M, 0, NULL, NULL);

cout<<"Done computeCoherence\n";

/*
 * Loop over M to find the max coherence and its pix index
 */

float max=0.f;
int maxPixIdx;

for(int idx=0; idx<nDir; idx++){
   if(M[idx] > max){
      max = M[idx];
      maxPixIdx = idx;
   }
}

summary->setMaxPixInfo(maxPixIdx, max);

/*
 * Write FITS file
 */
/*
cout<<"Creating Healpix map and writing to FITS....\n";
arr<float> MArr = arr<float>(&M[0], (size_t)nDir);
Healpix_Map<float> skyMap = Healpix_Map<float>(MArr, RING);

fitshandle fitsOut;
#ifdef CSW
char filename[] = "testCSWSkyMap.fits";
#else
char filename[] = "testXCorrSkyMap.fits";
#endif
remove(filename);
fitsOut.create(filename);

write_Healpix_map_to_fits(fitsOut, skyMap, PLANCK_FLOAT32);
cout<<"Healpix map written\n";
*/
/*
 * Deallocate memories
 */

cout<<"Deallocating memories...\n";
clReleaseMemObject(recoDelaysBuffer);
clReleaseMemObject(intensityRBuffer);
clReleaseMemObject(intensityCBuffer);
clReleaseMemObject(voltsFlatBuffer);
clReleaseMemObject(xCorrCBuffer);
clReleaseMemObject(xCorrRBuffer);
clReleaseMemObject(xCorrTimeBuffer);
clReleaseMemObject(sqrtWfPwrBuffer);
clReleaseMemObject(CijBuffer);
clReleaseMemObject(MBuffer);
free(voltsFlat);
//free(volts);
//free(recoDelays);
free(intensity_data_r);
free(intensity_data_c);
free(xCorr_data_r);
free(xCorr_data_c);
free(xCorrTime);
free(sqrtWfPwr);
free(Cij);
free(M);
cout<<"Memories deallocated\n";

return maxPixIdx;
}


int reconstructXCorrEnvelopeGetMaxPix(unsigned int dataType, vector<TGraph *>& cleanEvent, recoEnvData *clEnv,
                float *recoDelays, float *recoDelays_V, float *recoDelays_H,
                int nDir, string pol, const int *chanMask, recoData *summary)
{

cout<<"Entered reconstructXCorrEnvelopeGetMaxPix method\n";
int nSamp;
int nAnt = (int)cleanEvent.size()/2; // Divide by 2 for only one polarization
int unmaskedNChan=0;
for(int ch=0; ch<2*nAnt; ch++) unmaskedNChan+=chanMask[ch];
cout<<"unmaskedNChan: "<<unmaskedNChan<<" nAnt: "<<nAnt<<endl;
float wInt;
if( dataType == 0 ) wInt = 0.5f; //AraSim event
else if( dataType == 1 ){ //real event
if( pol == "vpol" ) wInt = 0.4f;
else if ( pol == "hpol" ) wInt = 0.625f;
else { cerr<<"recoPolType undefined\n"; return -1; }
} else {
cerr<<"dataType undefined\n"; return -1; }

/*
 * Loading voltsFlat array
 */

double t, v;
float *voltsFlat;

if( pol == "vpol" ){

   /* Using the 1st vpol wf for nSamp. Should make sure all vpol channels have the same nSamp. FIXME */
   nSamp = cleanEvent[0]->GetN();
   voltsFlat = (float*)calloc(nAnt*nSamp, sizeof(float));
   if(voltsFlat == NULL){ cerr<<"Null pointer to voltsFlat\n"; return -1; }

   for(int ch=0; ch<nAnt; ch++){
       if( chanMask[ch] ){
         for(int s=0; s<nSamp; s++){
            cleanEvent[ch]->GetPoint(s,t,v);
            /* Bartlett window applied in main analysis code */
            //voltsFlat[nSamp*ch + s] = static_cast<float>(v * FFTtools::bartlettWindow(s, nSamp));
            voltsFlat[nSamp*ch + s] = static_cast<float>(v);
         }
      } else {
         for(int s=0; s<nSamp; s++){
            voltsFlat[nSamp*ch + s] = 0.f;
         }
      }
   }

} else if ( pol == "hpol" ){

   /* Using the 1st hpol wf for nSamp. Should make sure all hpol channels have the same nSamp. FIXME */
   nSamp = cleanEvent[nAnt]->GetN();
   voltsFlat = (float*)calloc(nAnt*nSamp, sizeof(float));
   if(voltsFlat == NULL){ cerr<<"Null pointer to voltsFlat\n"; return -1; }

   for(int ch=0; ch<nAnt; ch++){
       if( chanMask[ch+nAnt] ){
         for(int s=0; s<nSamp; s++){
            cleanEvent[ch+nAnt]->GetPoint(s,t,v);
            /* Bartlett window applied in main analysis code */
            //voltsFlat[nSamp*ch + s] = static_cast<float>(v * FFTtools::bartlettWindow(s, nSamp));
            voltsFlat[nSamp*ch + s] = static_cast<float>(v);
         }
      } else {
         for(int s=0; s<nSamp; s++){
            voltsFlat[nSamp*ch + s] = 0.f;
         }
      }
   }

} else {
   cerr<<"recoPolType not defined\n"; return -1;
}
cout<<"voltsFlat loaded\n";

/*
 * Preparation for OUT_OF_PLACE transforms
 */

int interlvHermOutputSize = 2*(1 + nSamp/2); //Hermitian layout
int interlvOutputSize     = 2*nSamp;         //Not Hermitian layout
int planarHermOutputSize = (1 + nSamp/2);    //Hermitian planar layout

float *intensity_data_r, *intensity_data_c;
intensity_data_r = (float*)calloc(nAnt*planarHermOutputSize, sizeof(float));
intensity_data_c = (float*)calloc(nAnt*planarHermOutputSize, sizeof(float));

/*
 * FFT library related declarations
 */
//clfftPlanHandle planHandle;
clfftDim fftDim = CLFFT_1D;
size_t clLengths[1]   = {nSamp};
size_t clInStride[1]  = {1};
size_t clOutStride[1] = {1};
size_t inDist  = nSamp;
size_t outDist   = planarHermOutputSize;
size_t batchSize = nAnt;

/*
 * Set up clFFT
 */
//clfftSetupData fftSetup;
//err = clfftInitSetupData(&fftSetup);
//err = clfftSetup(&fftSetup);

/*
 * Prepare plan
 */
cout<<"Preparing plan...\n";
int err;
err = clfftCreateDefaultPlan(&clEnv->planHandle, clEnv->context, fftDim, clLengths);
err = clfftSetPlanPrecision(clEnv->planHandle, CLFFT_SINGLE);
err = clfftSetLayout(clEnv->planHandle, CLFFT_REAL, CLFFT_HERMITIAN_PLANAR);
err = clfftSetPlanScale(clEnv->planHandle, CLFFT_FORWARD, 1.f);
err = clfftSetPlanBatchSize(clEnv->planHandle, batchSize);
err = clfftSetPlanInStride(clEnv->planHandle, fftDim, clInStride);
err = clfftSetPlanOutStride(clEnv->planHandle, fftDim, clOutStride);
err = clfftSetPlanDistance(clEnv->planHandle, inDist, outDist);
err = clfftSetResultLocation(clEnv->planHandle, CLFFT_OUTOFPLACE);
err = clfftBakePlan(clEnv->planHandle, 1, &clEnv->queue, NULL, NULL);
cout<<"Plan prepared\n";
/* The plan is now ready to be executed */
cl_mem voltsFlatBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*nSamp*sizeof(float), NULL, &err);
err = clEnqueueWriteBuffer(clEnv->queue, voltsFlatBuffer, CL_TRUE, 0, nAnt*nSamp*sizeof(float), voltsFlat, 0, NULL, NULL);

cl_mem intensityRBuffer= clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*planarHermOutputSize*sizeof(float),
                         NULL, &err);
cl_mem intensityCBuffer= clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*planarHermOutputSize*sizeof(float),
                         NULL, &err);
cl_mem outBuffers[2] = {intensityRBuffer, intensityCBuffer};

cout<<"Enqueueing FFT\n";
err = clfftEnqueueTransform(clEnv->planHandle, CLFFT_FORWARD, 1, &clEnv->queue, 0, NULL, NULL,
                            &voltsFlatBuffer, outBuffers, NULL);
err = clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, outBuffers[0], CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_r,
                    0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, outBuffers[1], CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_c,
                    0, NULL, NULL);
cout<<"FFT done\n";
/*
 * Clean up CLFFT
 */
cout<<"Destroying plan...\n";
err = clfftDestroyPlan(&clEnv->planHandle);
//clfftTeardown();
cout<<"Plan destroyed\n";


#ifdef bandpass
/*
 * Bandpass signals
 */
cout<<"Preparing bandpass filter..."<<endl;
float freqBin = 1e3 / (wInt * (float)nSamp); // wInt in ns. 1e3 for MHz
float lowFreq = 200.f;
float highFreq= 450.f;

clSetKernelArg(clEnv->bandPassFilter, 0, sizeof(cl_mem), &intensityRBuffer);
clSetKernelArg(clEnv->bandPassFilter, 1, sizeof(cl_mem), &intensityCBuffer);
clSetKernelArg(clEnv->bandPassFilter, 2, sizeof(float),  &freqBin);
clSetKernelArg(clEnv->bandPassFilter, 3, sizeof(float),  &lowFreq);
clSetKernelArg(clEnv->bandPassFilter, 4, sizeof(float),  &highFreq);

unsigned int dim = 2;
size_t bandPassWorkSize[2] = {(size_t)nAnt, (size_t)planarHermOutputSize};
clEnqueueNDRangeKernel(clEnv->queue, clEnv->bandPassFilter, dim, NULL, bandPassWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, intensityRBuffer, CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_r,
                    0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, intensityCBuffer, CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_c,
                    0, NULL, NULL);

cout<<"Bandpass filter done"<<endl;
#endif

/*
 * Cross-correlate wfs
 */

int nBaseline = nAnt*nAnt;
cout<<"nBaseline: "<<nBaseline<<endl;

float *xCorr_data_r = (float*)calloc(nBaseline*planarHermOutputSize, sizeof(float));
float *xCorr_data_c = (float*)calloc(nBaseline*planarHermOutputSize, sizeof(float));

cl_mem xCorrRBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                    sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_r, &err);

cl_mem xCorrCBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                    sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_c, &err);

clSetKernelArg(clEnv->xCorrWf, 0, sizeof(cl_mem), &xCorrRBuffer);
clSetKernelArg(clEnv->xCorrWf, 1, sizeof(cl_mem), &xCorrCBuffer);
clSetKernelArg(clEnv->xCorrWf, 2, sizeof(cl_mem), &intensityRBuffer);
clSetKernelArg(clEnv->xCorrWf, 3, sizeof(cl_mem), &intensityCBuffer);
clSetKernelArg(clEnv->xCorrWf, 4, sizeof(int),    &nAnt);

unsigned int workDim = 2;
size_t globalWorkSize[2] = {(size_t)nBaseline, (size_t)planarHermOutputSize};
clEnqueueNDRangeKernel(clEnv->queue, clEnv->xCorrWf, workDim, NULL, globalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, xCorrRBuffer, CL_TRUE, 0, sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_r, 0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, xCorrCBuffer, CL_TRUE, 0, sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_c, 0, NULL, NULL);
cout<<"Done xCorrWf\n";

/*
 * Inverse FFT cross-correlation
 */

inDist  = planarHermOutputSize;
outDist = nSamp;
batchSize = nBaseline;

/*
 * Prepare plan
 */

err = clfftCreateDefaultPlan(&clEnv->planHandle, clEnv->context, fftDim, clLengths);
err = clfftSetPlanPrecision(clEnv->planHandle, CLFFT_SINGLE);
err = clfftSetLayout(clEnv->planHandle, CLFFT_HERMITIAN_PLANAR, CLFFT_REAL);
err = clfftSetPlanScale(clEnv->planHandle, CLFFT_BACKWARD, 1.f/(float)nSamp);
err = clfftSetPlanBatchSize(clEnv->planHandle, batchSize);
err = clfftSetPlanInStride(clEnv->planHandle, fftDim, clInStride);
err = clfftSetPlanOutStride(clEnv->planHandle, fftDim, clOutStride);
err = clfftSetPlanDistance(clEnv->planHandle, inDist, outDist);
err = clfftSetResultLocation(clEnv->planHandle, CLFFT_OUTOFPLACE);
err = clfftBakePlan(clEnv->planHandle, 1, &clEnv->queue, NULL, NULL);

/* The plan is now ready to be executed */

cl_mem inBuffers[2] = {xCorrRBuffer, xCorrCBuffer};

float *xCorrTime = (float*)calloc(nBaseline*nSamp, sizeof(float));
cl_mem xCorrTimeBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nBaseline*nSamp*sizeof(float),
                           NULL, &err);
err = clfftEnqueueTransform(clEnv->planHandle, CLFFT_BACKWARD, 1, &clEnv->queue, 0, NULL, NULL, inBuffers, &xCorrTimeBuffer, NULL);
err = clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, xCorrTimeBuffer, CL_TRUE, 0, sizeof(float)*nBaseline*nSamp, xCorrTime, 0, NULL, NULL);

/*
 * Clean up FFT
 */

err = clfftDestroyPlan(&clEnv->planHandle);
clfftTeardown();
cout<<"Done inverse FFT\n";

/*
 * Get Hilbert transform of XCorr function
 */

//TCanvas cvs("cvs","cva",800,600);
//cvs.Divide(1,2);


float dt[nSamp];
float xCorrValue[nSamp];
double t_temp, v_temp;

for(int baseline=0; baseline<nBaseline; baseline++){

   for(int s=0; s<nSamp; s++){

   dt[s] = wInt*s;
   xCorrValue[s] = xCorrTime[nSamp*baseline + s];

   }

   TGraph *xCorrGraph = new TGraph(nSamp, dt, xCorrValue);

   //cvs.cd(1);
   //xCorrGraph->Draw("AL");

   TGraph* envelope = FFTtools::getHilbertEnvelope( xCorrGraph );

   //cvs.cd(2);
   //envelope->Draw("AL");

   //cvs.SaveAs("xCorrEnvelope.C");

   for(int s=0; s<nSamp; s++){

   envelope->GetPoint(s,t_temp,v_temp);
   xCorrTime[nSamp*baseline + s] = static_cast<float>(v_temp);

   }

  delete xCorrGraph;
  delete envelope;

}

cl_mem xCorrEnvBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                       sizeof(float)*nBaseline*nSamp,
                                       xCorrTime, &err);


/*
 * Compute cross-correlation coefficients Cij in each direction
 */

/* Compute the square root of each channel's total wf power */

float *sqrtWfPwr = (float*)calloc(nAnt, sizeof(float));
float pwr=0.f;
for(int ant=0; ant<nAnt; ant++){
   pwr = 0.;
   for(int s=0; s<nSamp; s++){
      pwr += (voltsFlat[ant*nSamp + s] * voltsFlat[ant*nSamp + s]);
   }
   sqrtWfPwr[ant] = sqrt(pwr);
   //sqrtWfPwr[ant] = 1.f;
}
cout<<"Done sqrtWfPwr\n";

cl_mem sqrtWfPwrBuffer  = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt,
                          sqrtWfPwr, &err);
cl_mem recoDelaysBuffer;
if(pol == "vpol" ){
recoDelaysBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt*nDir,
                                  recoDelays_V, &err);
} else if (pol == "hpol" ) {
recoDelaysBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt*nDir,
                                  recoDelays_H, &err);
} else {
cerr<<"recoPolType not defined\n"; return -1;
}
cout<<"recoDelaysBuffer created\n";

float *Cij = (float*)calloc(nDir*nBaseline, sizeof(float));

cl_mem CijBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nDir*nBaseline,
                   Cij, &err);

clSetKernelArg(clEnv->computeXCorrCoef, 0, sizeof(cl_mem), &CijBuffer);
clSetKernelArg(clEnv->computeXCorrCoef, 1, sizeof(cl_mem), &xCorrEnvBuffer);
clSetKernelArg(clEnv->computeXCorrCoef, 2, sizeof(cl_mem), &recoDelaysBuffer);
clSetKernelArg(clEnv->computeXCorrCoef, 3, sizeof(cl_mem), &sqrtWfPwrBuffer);
clSetKernelArg(clEnv->computeXCorrCoef, 4, sizeof(float),  &wInt);
clSetKernelArg(clEnv->computeXCorrCoef, 5, sizeof(int),    &nAnt);
clSetKernelArg(clEnv->computeXCorrCoef, 6, sizeof(int),    &nSamp);

workDim = 2;
size_t CijGlobalWorkSize[2] = {(size_t)nDir, (size_t)nBaseline};

err = clEnqueueNDRangeKernel(clEnv->queue, clEnv->computeXCorrCoef, workDim, NULL, CijGlobalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, CijBuffer, CL_TRUE, 0, sizeof(float)*nDir*nBaseline, Cij, 0, NULL, NULL);

cout<<"Done computeXCorrCoef\n";

/*
 * Sum Cij's of all baselines in each reco direction to obtain coherence M(r-hat)
 */

float *M = (float*)calloc(nDir, sizeof(float));
cl_mem MBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nDir,
                  M, &err);

clSetKernelArg(clEnv->computeCoherence, 0, sizeof(cl_mem), &MBuffer);
clSetKernelArg(clEnv->computeCoherence, 1, sizeof(cl_mem), &CijBuffer);
clSetKernelArg(clEnv->computeCoherence, 2, sizeof(int),    &nBaseline);

workDim = 1;
size_t MGlobalWorkSize = nDir;

clEnqueueNDRangeKernel(clEnv->queue, clEnv->computeCoherence, workDim, NULL, &MGlobalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, MBuffer, CL_TRUE, 0, sizeof(float)*nDir, M, 0, NULL, NULL);

cout<<"Done computeCoherence\n";

/*
 * Loop over M to find the max coherence and its pix index
 */

float max=0.f;
int maxPixIdx;

for(int idx=0; idx<nDir; idx++){
   if(M[idx] > max){
      max = M[idx];
      maxPixIdx = idx;
   }
}

summary->setMaxPixInfo(maxPixIdx, max);

/*
 * Write FITS file
 */
/*
cout<<"Creating Healpix map and writing to FITS....\n";
arr<float> MArr = arr<float>(&M[0], (size_t)nDir);
Healpix_Map<float> skyMap = Healpix_Map<float>(MArr, RING);

fitshandle fitsOut;
#ifdef CSW
char filename[] = "testCSWSkyMap.fits";
#else
char filename[] = "testXCorrSkyMap.fits";
#endif
remove(filename);
fitsOut.create(filename);

write_Healpix_map_to_fits(fitsOut, skyMap, PLANCK_FLOAT32);
cout<<"Healpix map written\n";
*/
/*
 * Deallocate memories
 */

cout<<"Deallocating memories...\n";
clReleaseMemObject(recoDelaysBuffer);
clReleaseMemObject(intensityRBuffer);
clReleaseMemObject(intensityCBuffer);
clReleaseMemObject(voltsFlatBuffer);
clReleaseMemObject(xCorrCBuffer);
clReleaseMemObject(xCorrRBuffer);
clReleaseMemObject(xCorrTimeBuffer);
clReleaseMemObject(xCorrEnvBuffer);
clReleaseMemObject(sqrtWfPwrBuffer);
clReleaseMemObject(CijBuffer);
clReleaseMemObject(MBuffer);
free(voltsFlat);
//free(volts);
//free(recoDelays);
free(intensity_data_r);
free(intensity_data_c);
free(xCorr_data_r);
free(xCorr_data_c);
free(xCorrTime);
free(sqrtWfPwr);
free(Cij);
free(M);
cout<<"Memories deallocated\n";

return maxPixIdx;
}

int reconstructXCorrEnvelopeGetMaxPix(recoSettings *settings, vector<TGraph *>& cleanEvent, recoEnvData *clEnv,
                float *recoDelays, float *recoDelays_V, float *recoDelays_H,
                int nDir, const int *chanMask, recoData *summary)
{

cout<<"Entered reconstructXCorrEnvelopeGetMaxPix method\n";
int nSamp;
int nAnt = (int)cleanEvent.size()/2; // Divide by 2 for only one polarization
int unmaskedNChan=0;
string pol = string(settings->recoPolType);
for(int ch=0; ch<2*nAnt; ch++) unmaskedNChan+=chanMask[ch];
cout<<"unmaskedNChan: "<<unmaskedNChan<<" nAnt: "<<nAnt<<endl;
float wInt;
int dataType = settings->dataType;
if( dataType == 0 ) wInt = 0.5f; //AraSim event
else if( dataType == 1 ){ //real event
if( pol == "vpol" ) wInt = 0.4f;
else if ( pol == "hpol" ) wInt = 0.625f;
else { cerr<<"recoPolType undefined\n"; return -1; }
} else {
cerr<<"dataType undefined\n"; return -1; }

/*
 * Loading voltsFlat array
 */

double t, v;
float *voltsFlat;

if( pol == "vpol" ){

   /* Using the 1st vpol wf for nSamp. Should make sure all vpol channels have the same nSamp. FIXME */
   nSamp = cleanEvent[0]->GetN();
   voltsFlat = (float*)calloc(nAnt*nSamp, sizeof(float));
   if(voltsFlat == NULL){ cerr<<"Null pointer to voltsFlat\n"; return -1; }

   for(int ch=0; ch<nAnt; ch++){
       if( chanMask[ch] ){
         for(int s=0; s<nSamp; s++){
            cleanEvent[ch]->GetPoint(s,t,v);
            /* Bartlett window applied in main analysis code */
            //voltsFlat[nSamp*ch + s] = static_cast<float>(v * FFTtools::bartlettWindow(s, nSamp));
            voltsFlat[nSamp*ch + s] = static_cast<float>(v);
         }
      } else {
         for(int s=0; s<nSamp; s++){
            voltsFlat[nSamp*ch + s] = 0.f;
         }
      }
   }

} else if ( pol == "hpol" ){

   /* Using the 1st hpol wf for nSamp. Should make sure all hpol channels have the same nSamp. FIXME */
   nSamp = cleanEvent[nAnt]->GetN();
   voltsFlat = (float*)calloc(nAnt*nSamp, sizeof(float));
   if(voltsFlat == NULL){ cerr<<"Null pointer to voltsFlat\n"; return -1; }

   for(int ch=0; ch<nAnt; ch++){
       if( chanMask[ch+nAnt] ){
         for(int s=0; s<nSamp; s++){
            cleanEvent[ch+nAnt]->GetPoint(s,t,v);
            /* Bartlett window applied in main analysis code */
            //voltsFlat[nSamp*ch + s] = static_cast<float>(v * FFTtools::bartlettWindow(s, nSamp));
            voltsFlat[nSamp*ch + s] = static_cast<float>(v);
         }
      } else {
         for(int s=0; s<nSamp; s++){
            voltsFlat[nSamp*ch + s] = 0.f;
         }
      }
   }

} else {
   cerr<<"recoPolType not defined\n"; return -1;
}
cout<<"voltsFlat loaded\n";

/*
 * Preparation for OUT_OF_PLACE transforms
 */

int interlvHermOutputSize = 2*(1 + nSamp/2); //Hermitian layout
int interlvOutputSize     = 2*nSamp;         //Not Hermitian layout
int planarHermOutputSize = (1 + nSamp/2);    //Hermitian planar layout

float *intensity_data_r, *intensity_data_c;
intensity_data_r = (float*)calloc(nAnt*planarHermOutputSize, sizeof(float));
intensity_data_c = (float*)calloc(nAnt*planarHermOutputSize, sizeof(float));

/*
 * FFT library related declarations
 */
//clfftPlanHandle planHandle;
clfftDim fftDim = CLFFT_1D;
size_t clLengths[1]   = {nSamp};
size_t clInStride[1]  = {1};
size_t clOutStride[1] = {1};
size_t inDist  = nSamp;
size_t outDist   = planarHermOutputSize;
size_t batchSize = nAnt;

/*
 * Set up clFFT
 */
//clfftSetupData fftSetup;
//err = clfftInitSetupData(&fftSetup);
//err = clfftSetup(&fftSetup);

/*
 * Prepare plan
 */
cout<<"Preparing plan...\n";
int err;
err = clfftCreateDefaultPlan(&clEnv->planHandle, clEnv->context, fftDim, clLengths);
err = clfftSetPlanPrecision(clEnv->planHandle, CLFFT_SINGLE);
err = clfftSetLayout(clEnv->planHandle, CLFFT_REAL, CLFFT_HERMITIAN_PLANAR);
err = clfftSetPlanScale(clEnv->planHandle, CLFFT_FORWARD, 1.f);
err = clfftSetPlanBatchSize(clEnv->planHandle, batchSize);
err = clfftSetPlanInStride(clEnv->planHandle, fftDim, clInStride);
err = clfftSetPlanOutStride(clEnv->planHandle, fftDim, clOutStride);
err = clfftSetPlanDistance(clEnv->planHandle, inDist, outDist);
err = clfftSetResultLocation(clEnv->planHandle, CLFFT_OUTOFPLACE);
err = clfftBakePlan(clEnv->planHandle, 1, &clEnv->queue, NULL, NULL);
cout<<"Plan prepared\n";
/* The plan is now ready to be executed */
cl_mem voltsFlatBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*nSamp*sizeof(float), NULL, &err);
err = clEnqueueWriteBuffer(clEnv->queue, voltsFlatBuffer, CL_TRUE, 0, nAnt*nSamp*sizeof(float), voltsFlat, 0, NULL, NULL);

cl_mem intensityRBuffer= clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*planarHermOutputSize*sizeof(float),
                         NULL, &err);
cl_mem intensityCBuffer= clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*planarHermOutputSize*sizeof(float),
                         NULL, &err);
cl_mem outBuffers[2] = {intensityRBuffer, intensityCBuffer};

cout<<"Enqueueing FFT\n";
err = clfftEnqueueTransform(clEnv->planHandle, CLFFT_FORWARD, 1, &clEnv->queue, 0, NULL, NULL,
                            &voltsFlatBuffer, outBuffers, NULL);
err = clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, outBuffers[0], CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_r,
                    0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, outBuffers[1], CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_c,
                    0, NULL, NULL);
cout<<"FFT done\n";
/*
 * Clean up CLFFT
 */
cout<<"Destroying plan...\n";
err = clfftDestroyPlan(&clEnv->planHandle);
//clfftTeardown();
cout<<"Plan destroyed\n";


#ifdef bandpass
/*
 * Bandpass signals
 */
cout<<"Preparing bandpass filter..."<<endl;
float freqBin = 1e3 / (wInt * (float)nSamp); // wInt in ns. 1e3 for MHz
float lowFreq = 200.f;
float highFreq= 450.f;

clSetKernelArg(clEnv->bandPassFilter, 0, sizeof(cl_mem), &intensityRBuffer);
clSetKernelArg(clEnv->bandPassFilter, 1, sizeof(cl_mem), &intensityCBuffer);
clSetKernelArg(clEnv->bandPassFilter, 2, sizeof(float),  &freqBin);
clSetKernelArg(clEnv->bandPassFilter, 3, sizeof(float),  &lowFreq);
clSetKernelArg(clEnv->bandPassFilter, 4, sizeof(float),  &highFreq);

unsigned int dim = 2;
size_t bandPassWorkSize[2] = {(size_t)nAnt, (size_t)planarHermOutputSize};
clEnqueueNDRangeKernel(clEnv->queue, clEnv->bandPassFilter, dim, NULL, bandPassWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, intensityRBuffer, CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_r,
                    0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, intensityCBuffer, CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_c,
                    0, NULL, NULL);

cout<<"Bandpass filter done"<<endl;
#endif

/*
 * Cross-correlate wfs
 */

int nBaseline = nAnt*nAnt;
cout<<"nBaseline: "<<nBaseline<<endl;

float *xCorr_data_r = (float*)calloc(nBaseline*planarHermOutputSize, sizeof(float));
float *xCorr_data_c = (float*)calloc(nBaseline*planarHermOutputSize, sizeof(float));

cl_mem xCorrRBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                    sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_r, &err);

cl_mem xCorrCBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                    sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_c, &err);

clSetKernelArg(clEnv->xCorrWf, 0, sizeof(cl_mem), &xCorrRBuffer);
clSetKernelArg(clEnv->xCorrWf, 1, sizeof(cl_mem), &xCorrCBuffer);
clSetKernelArg(clEnv->xCorrWf, 2, sizeof(cl_mem), &intensityRBuffer);
clSetKernelArg(clEnv->xCorrWf, 3, sizeof(cl_mem), &intensityCBuffer);
clSetKernelArg(clEnv->xCorrWf, 4, sizeof(int),    &nAnt);

unsigned int workDim = 2;
size_t globalWorkSize[2] = {(size_t)nBaseline, (size_t)planarHermOutputSize};
clEnqueueNDRangeKernel(clEnv->queue, clEnv->xCorrWf, workDim, NULL, globalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, xCorrRBuffer, CL_TRUE, 0, sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_r, 0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, xCorrCBuffer, CL_TRUE, 0, sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_c, 0, NULL, NULL);
cout<<"Done xCorrWf\n";

/*
 * Inverse FFT cross-correlation
 */

inDist  = planarHermOutputSize;
outDist = nSamp;
batchSize = nBaseline;

/*
 * Prepare plan
 */

err = clfftCreateDefaultPlan(&clEnv->planHandle, clEnv->context, fftDim, clLengths);
err = clfftSetPlanPrecision(clEnv->planHandle, CLFFT_SINGLE);
err = clfftSetLayout(clEnv->planHandle, CLFFT_HERMITIAN_PLANAR, CLFFT_REAL);
err = clfftSetPlanScale(clEnv->planHandle, CLFFT_BACKWARD, 1.f/(float)nSamp);
err = clfftSetPlanBatchSize(clEnv->planHandle, batchSize);
err = clfftSetPlanInStride(clEnv->planHandle, fftDim, clInStride);
err = clfftSetPlanOutStride(clEnv->planHandle, fftDim, clOutStride);
err = clfftSetPlanDistance(clEnv->planHandle, inDist, outDist);
err = clfftSetResultLocation(clEnv->planHandle, CLFFT_OUTOFPLACE);
err = clfftBakePlan(clEnv->planHandle, 1, &clEnv->queue, NULL, NULL);

/* The plan is now ready to be executed */

cl_mem inBuffers[2] = {xCorrRBuffer, xCorrCBuffer};

float *xCorrTime = (float*)calloc(nBaseline*nSamp, sizeof(float));
cl_mem xCorrTimeBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nBaseline*nSamp*sizeof(float),
                           NULL, &err);
err = clfftEnqueueTransform(clEnv->planHandle, CLFFT_BACKWARD, 1, &clEnv->queue, 0, NULL, NULL, inBuffers, &xCorrTimeBuffer, NULL);
err = clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, xCorrTimeBuffer, CL_TRUE, 0, sizeof(float)*nBaseline*nSamp, xCorrTime, 0, NULL, NULL);

/*
 * Clean up FFT
 */

err = clfftDestroyPlan(&clEnv->planHandle);
clfftTeardown();
cout<<"Done inverse FFT\n";

/*
 * Get Hilbert transform of XCorr function
 */

//TCanvas cvs("cvs","cva",800,600);
//cvs.Divide(1,2);


float dt[nSamp];
float xCorrValue[nSamp];
double t_temp, v_temp;

for(int baseline=0; baseline<nBaseline; baseline++){

   for(int s=0; s<nSamp; s++){

   dt[s] = wInt*s;
   xCorrValue[s] = xCorrTime[nSamp*baseline + s];

   }

   TGraph *xCorrGraph = new TGraph(nSamp, dt, xCorrValue);

   //cvs.cd(1);
   //xCorrGraph->Draw("AL");

   TGraph* envelope = FFTtools::getHilbertEnvelope( xCorrGraph );

   //cvs.cd(2);
   //envelope->Draw("AL");

   //cvs.SaveAs("xCorrEnvelope.C");

   for(int s=0; s<nSamp; s++){

   envelope->GetPoint(s,t_temp,v_temp);
   xCorrTime[nSamp*baseline + s] = static_cast<float>(v_temp);

   }

  delete xCorrGraph;
  delete envelope;

}

cl_mem xCorrEnvBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                       sizeof(float)*nBaseline*nSamp,
                                       xCorrTime, &err);


/*
 * Compute cross-correlation coefficients Cij in each direction
 */

/* Compute the square root of each channel's total wf power */

float *sqrtWfPwr = (float*)calloc(nAnt, sizeof(float));
float pwr=0.f;
for(int ant=0; ant<nAnt; ant++){
   pwr = 0.;
   for(int s=0; s<nSamp; s++){
      pwr += (voltsFlat[ant*nSamp + s] * voltsFlat[ant*nSamp + s]);
   }
   sqrtWfPwr[ant] = sqrt(pwr);
   //sqrtWfPwr[ant] = 1.f;
}
cout<<"Done sqrtWfPwr\n";

cl_mem sqrtWfPwrBuffer  = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt,
                          sqrtWfPwr, &err);
cl_mem recoDelaysBuffer;
if(pol == "vpol" ){
recoDelaysBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt*nDir,
                                  recoDelays_V, &err);
} else if (pol == "hpol" ) {
recoDelaysBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt*nDir,
                                  recoDelays_H, &err);
} else {
cerr<<"recoPolType not defined\n"; return -1;
}
cout<<"recoDelaysBuffer created\n";

float *Cij = (float*)calloc(nDir*nBaseline, sizeof(float));

cl_mem CijBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nDir*nBaseline,
                   Cij, &err);

clSetKernelArg(clEnv->computeXCorrCoef, 0, sizeof(cl_mem), &CijBuffer);
clSetKernelArg(clEnv->computeXCorrCoef, 1, sizeof(cl_mem), &xCorrEnvBuffer);
clSetKernelArg(clEnv->computeXCorrCoef, 2, sizeof(cl_mem), &recoDelaysBuffer);
clSetKernelArg(clEnv->computeXCorrCoef, 3, sizeof(cl_mem), &sqrtWfPwrBuffer);
clSetKernelArg(clEnv->computeXCorrCoef, 4, sizeof(float),  &wInt);
clSetKernelArg(clEnv->computeXCorrCoef, 5, sizeof(int),    &nAnt);
clSetKernelArg(clEnv->computeXCorrCoef, 6, sizeof(int),    &nSamp);

workDim = 2;
size_t CijGlobalWorkSize[2] = {(size_t)nDir, (size_t)nBaseline};

err = clEnqueueNDRangeKernel(clEnv->queue, clEnv->computeXCorrCoef, workDim, NULL, CijGlobalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, CijBuffer, CL_TRUE, 0, sizeof(float)*nDir*nBaseline, Cij, 0, NULL, NULL);

cout<<"Done computeXCorrCoef\n";

/*
 * Sum Cij's of all baselines in each reco direction to obtain coherence M(r-hat)
 */

float *M = (float*)calloc(nDir, sizeof(float));
cl_mem MBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nDir,
                  M, &err);

clSetKernelArg(clEnv->computeCoherence, 0, sizeof(cl_mem), &MBuffer);
clSetKernelArg(clEnv->computeCoherence, 1, sizeof(cl_mem), &CijBuffer);
clSetKernelArg(clEnv->computeCoherence, 2, sizeof(int),    &nBaseline);

workDim = 1;
size_t MGlobalWorkSize = nDir;

clEnqueueNDRangeKernel(clEnv->queue, clEnv->computeCoherence, workDim, NULL, &MGlobalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, MBuffer, CL_TRUE, 0, sizeof(float)*nDir, M, 0, NULL, NULL);

cout<<"Done computeCoherence\n";

/*
 * Loop over M to find the max coherence and its pix index
 */

float max=0.f;
int maxPixIdx;

for(int idx=0; idx<nDir; idx++){
   if(M[idx] > max){
      max = M[idx];
      maxPixIdx = idx;
   }
}

summary->setMaxPixInfo(maxPixIdx, max);

/*
 * Write FITS file
 */
/*
cout<<"Creating Healpix map and writing to FITS....\n";
arr<float> MArr = arr<float>(&M[0], (size_t)nDir);
Healpix_Map<float> skyMap = Healpix_Map<float>(MArr, RING);

fitshandle fitsOut;
#ifdef CSW
char filename[] = "testCSWSkyMap.fits";
#else
char filename[] = "testXCorrSkyMap.fits";
#endif
remove(filename);
fitsOut.create(filename);

write_Healpix_map_to_fits(fitsOut, skyMap, PLANCK_FLOAT32);
cout<<"Healpix map written\n";
*/
/*
 * Deallocate memories
 */

cout<<"Deallocating memories...\n";
clReleaseMemObject(recoDelaysBuffer);
clReleaseMemObject(intensityRBuffer);
clReleaseMemObject(intensityCBuffer);
clReleaseMemObject(voltsFlatBuffer);
clReleaseMemObject(xCorrCBuffer);
clReleaseMemObject(xCorrRBuffer);
clReleaseMemObject(xCorrTimeBuffer);
clReleaseMemObject(xCorrEnvBuffer);
clReleaseMemObject(sqrtWfPwrBuffer);
clReleaseMemObject(CijBuffer);
clReleaseMemObject(MBuffer);
free(voltsFlat);
//free(volts);
//free(recoDelays);
free(intensity_data_r);
free(intensity_data_c);
free(xCorr_data_r);
free(xCorr_data_c);
free(xCorrTime);
free(sqrtWfPwr);
free(Cij);
free(M);
cout<<"Memories deallocated\n";

return maxPixIdx;
}

int reconstructXCorrGetMaxPixAndMap(unsigned int dataType, vector<TGraph *>& cleanEvent, recoEnvData *clEnv,
                float *recoDelays, float *recoDelays_V, float *recoDelays_H,
                int nDir, string pol, const int *chanMask, recoData *summary, char *filename)
{

cout<<"Entered reconstructXCorrGetMaxPixAndMap method\n";
int nSamp;
int nAnt = (int)cleanEvent.size()/2; // Divide by 2 for only one polarization
int unmaskedNChan=0;
for(int ch=0; ch<2*nAnt; ch++) unmaskedNChan+=chanMask[ch];
cout<<"unmaskedNChan: "<<unmaskedNChan<<" nAnt: "<<nAnt<<endl;
float wInt;
if( dataType == 0 ) wInt = 0.5f; //AraSim event
else if( dataType == 1 ){ //real event
if( pol == "vpol" ) wInt = 0.4f;
else if ( pol == "hpol" ) wInt = 0.625f;
else { cerr<<"recoPolType undefined\n"; return -1; }
} else {
cerr<<"dataType undefined\n"; return -1; }

/*
 * Loading voltsFlat array
 */

double t, v;
float *voltsFlat;

if( pol == "vpol" ){

   /* Using the 1st vpol wf for nSamp. Should make sure all vpol channels have the same nSamp. FIXME */
   nSamp = cleanEvent[0]->GetN();
   voltsFlat = (float*)calloc(nAnt*nSamp, sizeof(float));
   if(voltsFlat == NULL){ cerr<<"Null pointer to voltsFlat\n"; return -1; }

   for(int ch=0; ch<nAnt; ch++){
       if( chanMask[ch] ){
         for(int s=0; s<nSamp; s++){
            cleanEvent[ch]->GetPoint(s,t,v);
            /* Bartlett window applied in main analysis code */
            //voltsFlat[nSamp*ch + s] = static_cast<float>(v * FFTtools::bartlettWindow(s, nSamp));
            voltsFlat[nSamp*ch + s] = static_cast<float>(v);
         }
      } else {
         for(int s=0; s<nSamp; s++){
            voltsFlat[nSamp*ch + s] = 0.f;
         }
      }
   }

} else if ( pol == "hpol" ){

   /* Using the 1st hpol wf for nSamp. Should make sure all hpol channels have the same nSamp. FIXME */
   nSamp = cleanEvent[nAnt]->GetN();
   voltsFlat = (float*)calloc(nAnt*nSamp, sizeof(float));
   if(voltsFlat == NULL){ cerr<<"Null pointer to voltsFlat\n"; return -1; }

   for(int ch=0; ch<nAnt; ch++){
       if( chanMask[ch+nAnt] ){
         for(int s=0; s<nSamp; s++){
            cleanEvent[ch+nAnt]->GetPoint(s,t,v);
            /* Bartlett window applied in main analysis code */
            //voltsFlat[nSamp*ch + s] = static_cast<float>(v * FFTtools::bartlettWindow(s, nSamp));
            voltsFlat[nSamp*ch + s] = static_cast<float>(v);
         }
      } else {
         for(int s=0; s<nSamp; s++){
            voltsFlat[nSamp*ch + s] = 0.f;
         }
      }
   }

} else {
   cerr<<"recoPolType not defined\n"; return -1;
}
cout<<"voltsFlat loaded\n";

/*
 * Preparation for OUT_OF_PLACE transforms
 */

int interlvHermOutputSize = 2*(1 + nSamp/2); //Hermitian layout
int interlvOutputSize     = 2*nSamp;         //Not Hermitian layout
int planarHermOutputSize = (1 + nSamp/2);    //Hermitian planar layout

float *intensity_data_r, *intensity_data_c;
intensity_data_r = (float*)calloc(nAnt*planarHermOutputSize, sizeof(float));
intensity_data_c = (float*)calloc(nAnt*planarHermOutputSize, sizeof(float));

/*
 * FFT library related declarations
 */
//clfftPlanHandle planHandle;
clfftDim fftDim = CLFFT_1D;
size_t clLengths[1]   = {nSamp};
size_t clInStride[1]  = {1};
size_t clOutStride[1] = {1};
size_t inDist  = nSamp;
size_t outDist   = planarHermOutputSize;
size_t batchSize = nAnt;

/*
 * Set up clFFT
 */
//clfftSetupData fftSetup;
//err = clfftInitSetupData(&fftSetup);
//err = clfftSetup(&fftSetup);

/*
 * Prepare plan
 */
cout<<"Preparing plan...\n";
int err;
err = clfftCreateDefaultPlan(&clEnv->planHandle, clEnv->context, fftDim, clLengths);
err = clfftSetPlanPrecision(clEnv->planHandle, CLFFT_SINGLE);
err = clfftSetLayout(clEnv->planHandle, CLFFT_REAL, CLFFT_HERMITIAN_PLANAR);
err = clfftSetPlanScale(clEnv->planHandle, CLFFT_FORWARD, 1.f);
err = clfftSetPlanBatchSize(clEnv->planHandle, batchSize);
err = clfftSetPlanInStride(clEnv->planHandle, fftDim, clInStride);
err = clfftSetPlanOutStride(clEnv->planHandle, fftDim, clOutStride);
err = clfftSetPlanDistance(clEnv->planHandle, inDist, outDist);
err = clfftSetResultLocation(clEnv->planHandle, CLFFT_OUTOFPLACE);
err = clfftBakePlan(clEnv->planHandle, 1, &clEnv->queue, NULL, NULL);
cout<<"Plan prepared\n";
/* The plan is now ready to be executed */
cl_mem voltsFlatBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*nSamp*sizeof(float), NULL, &err);
err = clEnqueueWriteBuffer(clEnv->queue, voltsFlatBuffer, CL_TRUE, 0, nAnt*nSamp*sizeof(float), voltsFlat, 0, NULL, NULL);

cl_mem intensityRBuffer= clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*planarHermOutputSize*sizeof(float),
                         NULL, &err);
cl_mem intensityCBuffer= clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*planarHermOutputSize*sizeof(float),
                         NULL, &err);
cl_mem outBuffers[2] = {intensityRBuffer, intensityCBuffer};

cout<<"Enqueueing FFT\n";
err = clfftEnqueueTransform(clEnv->planHandle, CLFFT_FORWARD, 1, &clEnv->queue, 0, NULL, NULL,
                            &voltsFlatBuffer, outBuffers, NULL);
err = clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, outBuffers[0], CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_r,
                    0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, outBuffers[1], CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_c,
                    0, NULL, NULL);
cout<<"FFT done\n";
/*
 * Clean up CLFFT
 */
cout<<"Destroying plan...\n";
err = clfftDestroyPlan(&clEnv->planHandle);
//clfftTeardown();
cout<<"Plan destroyed\n";


#ifdef bandpass
/*
 * Bandpass signals
 */
cout<<"Preparing bandpass filter..."<<endl;
float freqBin = 1e3 / (wInt * (float)nSamp); // wInt in ns. 1e3 for MHz
float lowFreq = 200.f;
float highFreq= 450.f;

clSetKernelArg(clEnv->bandPassFilter, 0, sizeof(cl_mem), &intensityRBuffer);
clSetKernelArg(clEnv->bandPassFilter, 1, sizeof(cl_mem), &intensityCBuffer);
clSetKernelArg(clEnv->bandPassFilter, 2, sizeof(float),  &freqBin);
clSetKernelArg(clEnv->bandPassFilter, 3, sizeof(float),  &lowFreq);
clSetKernelArg(clEnv->bandPassFilter, 4, sizeof(float),  &highFreq);

unsigned int dim = 2;
size_t bandPassWorkSize[2] = {(size_t)nAnt, (size_t)planarHermOutputSize};
clEnqueueNDRangeKernel(clEnv->queue, clEnv->bandPassFilter, dim, NULL, bandPassWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, intensityRBuffer, CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_r,
                    0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, intensityCBuffer, CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_c,
                    0, NULL, NULL);

cout<<"Bandpass filter done"<<endl;
#endif

/*
 * Cross-correlate wfs
 */

int nBaseline = nAnt*nAnt;
cout<<"nBaseline: "<<nBaseline<<endl;

float *xCorr_data_r = (float*)calloc(nBaseline*planarHermOutputSize, sizeof(float));
float *xCorr_data_c = (float*)calloc(nBaseline*planarHermOutputSize, sizeof(float));

cl_mem xCorrRBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                    sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_r, &err);

cl_mem xCorrCBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                    sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_c, &err);

clSetKernelArg(clEnv->xCorrWf, 0, sizeof(cl_mem), &xCorrRBuffer);
clSetKernelArg(clEnv->xCorrWf, 1, sizeof(cl_mem), &xCorrCBuffer);
clSetKernelArg(clEnv->xCorrWf, 2, sizeof(cl_mem), &intensityRBuffer);
clSetKernelArg(clEnv->xCorrWf, 3, sizeof(cl_mem), &intensityCBuffer);
clSetKernelArg(clEnv->xCorrWf, 4, sizeof(int),    &nAnt);

unsigned int workDim = 2;
size_t globalWorkSize[2] = {(size_t)nBaseline, (size_t)planarHermOutputSize};
clEnqueueNDRangeKernel(clEnv->queue, clEnv->xCorrWf, workDim, NULL, globalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, xCorrRBuffer, CL_TRUE, 0, sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_r, 0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, xCorrCBuffer, CL_TRUE, 0, sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_c, 0, NULL, NULL);
cout<<"Done xCorrWf\n";

/*
 * Inverse FFT cross-correlation
 */

inDist  = planarHermOutputSize;
outDist = nSamp;
batchSize = nBaseline;

/*
 * Prepare plan
 */

err = clfftCreateDefaultPlan(&clEnv->planHandle, clEnv->context, fftDim, clLengths);
err = clfftSetPlanPrecision(clEnv->planHandle, CLFFT_SINGLE);
err = clfftSetLayout(clEnv->planHandle, CLFFT_HERMITIAN_PLANAR, CLFFT_REAL);
err = clfftSetPlanScale(clEnv->planHandle, CLFFT_BACKWARD, 1.f/(float)nSamp);
err = clfftSetPlanBatchSize(clEnv->planHandle, batchSize);
err = clfftSetPlanInStride(clEnv->planHandle, fftDim, clInStride);
err = clfftSetPlanOutStride(clEnv->planHandle, fftDim, clOutStride);
err = clfftSetPlanDistance(clEnv->planHandle, inDist, outDist);
err = clfftSetResultLocation(clEnv->planHandle, CLFFT_OUTOFPLACE);
err = clfftBakePlan(clEnv->planHandle, 1, &clEnv->queue, NULL, NULL);

/* The plan is now ready to be executed */

cl_mem inBuffers[2] = {xCorrRBuffer, xCorrCBuffer};

float *xCorrTime = (float*)calloc(nBaseline*nSamp, sizeof(float));
cl_mem xCorrTimeBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nBaseline*nSamp*sizeof(float),
                           NULL, &err);
err = clfftEnqueueTransform(clEnv->planHandle, CLFFT_BACKWARD, 1, &clEnv->queue, 0, NULL, NULL, inBuffers, &xCorrTimeBuffer, NULL);
err = clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, xCorrTimeBuffer, CL_TRUE, 0, sizeof(float)*nBaseline*nSamp, xCorrTime, 0, NULL, NULL);

/*
 * Clean up FFT
 */

err = clfftDestroyPlan(&clEnv->planHandle);
cout<<"Done inverse FFT\n";

/*
 * Compute cross-correlation coefficients Cij in each direction
 */

/* Compute the square root of each channel's total wf power */

float *sqrtWfPwr = (float*)calloc(nAnt, sizeof(float));
float pwr=0.f;
for(int ant=0; ant<nAnt; ant++){
   pwr = 0.;
   for(int s=0; s<nSamp; s++){
      pwr += (voltsFlat[ant*nSamp + s] * voltsFlat[ant*nSamp + s]);
   }
   sqrtWfPwr[ant] = sqrt(pwr);
   //sqrtWfPwr[ant] = 1.f;
}
cout<<"Done sqrtWfPwr\n";

cl_mem sqrtWfPwrBuffer  = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt,
                          sqrtWfPwr, &err);
cl_mem recoDelaysBuffer;
if(pol == "vpol" ){
recoDelaysBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt*nDir,
                                  recoDelays_V, &err);
} else if (pol == "hpol" ) {
recoDelaysBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt*nDir,
                                  recoDelays_H, &err);
} else {
cerr<<"recoPolType not defined\n"; return -1;
}
cout<<"recoDelaysBuffer created\n";

float *Cij = (float*)calloc(nDir*nBaseline, sizeof(float));

cl_mem CijBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nDir*nBaseline,
                   Cij, &err);

clSetKernelArg(clEnv->computeXCorrCoef, 0, sizeof(cl_mem), &CijBuffer);
clSetKernelArg(clEnv->computeXCorrCoef, 1, sizeof(cl_mem), &xCorrTimeBuffer);
clSetKernelArg(clEnv->computeXCorrCoef, 2, sizeof(cl_mem), &recoDelaysBuffer);
clSetKernelArg(clEnv->computeXCorrCoef, 3, sizeof(cl_mem), &sqrtWfPwrBuffer);
clSetKernelArg(clEnv->computeXCorrCoef, 4, sizeof(float),  &wInt);
clSetKernelArg(clEnv->computeXCorrCoef, 5, sizeof(int),    &nAnt);
clSetKernelArg(clEnv->computeXCorrCoef, 6, sizeof(int),    &nSamp);

workDim = 2;
size_t CijGlobalWorkSize[2] = {(size_t)nDir, (size_t)nBaseline};

err = clEnqueueNDRangeKernel(clEnv->queue, clEnv->computeXCorrCoef, workDim, NULL, CijGlobalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, CijBuffer, CL_TRUE, 0, sizeof(float)*nDir*nBaseline, Cij, 0, NULL, NULL);

cout<<"Done computeXCorrCoef\n";

/*
 * Sum Cij's of all baselines in each reco direction to obtain coherence M(r-hat)
 */

float *M = (float*)calloc(nDir, sizeof(float));
cl_mem MBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nDir,
                  M, &err);

clSetKernelArg(clEnv->computeCoherence, 0, sizeof(cl_mem), &MBuffer);
clSetKernelArg(clEnv->computeCoherence, 1, sizeof(cl_mem), &CijBuffer);
clSetKernelArg(clEnv->computeCoherence, 2, sizeof(int),    &nBaseline);

workDim = 1;
size_t MGlobalWorkSize = nDir;

clEnqueueNDRangeKernel(clEnv->queue, clEnv->computeCoherence, workDim, NULL, &MGlobalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, MBuffer, CL_TRUE, 0, sizeof(float)*nDir, M, 0, NULL, NULL);

cout<<"Done computeCoherence\n";

/*
 * Loop over M to find the max coherence and its pix index
 */

float max=0.f;
int maxPixIdx;

for(int idx=0; idx<nDir; idx++){
   if(M[idx] > max){
      max = M[idx];
      maxPixIdx = idx;
   }
}

summary->setMaxPixInfo(maxPixIdx, max);

/*
 * Write FITS file
 */

cout<<"Creating Healpix map and writing to FITS....\n";
arr<float> MArr = arr<float>(&M[0], (size_t)nDir);
Healpix_Map<float> skyMap = Healpix_Map<float>(MArr, HEALPIX_ORDERING);

fitshandle fitsOut;
//#ifdef CSW
//char filename[] = "testCSWSkyMap.fits";
//#else
//char filename[] = "testXCorrSkyMap.fits";
//#endif
remove(filename);
fitsOut.create(filename);

write_Healpix_map_to_fits(fitsOut, skyMap, PLANCK_FLOAT32);
cout<<"Healpix map written\n";

/*
 * Deallocate memories
 */

cout<<"Deallocating memories...\n";
clReleaseMemObject(recoDelaysBuffer);
clReleaseMemObject(intensityRBuffer);
clReleaseMemObject(intensityCBuffer);
clReleaseMemObject(voltsFlatBuffer);
clReleaseMemObject(xCorrCBuffer);
clReleaseMemObject(xCorrRBuffer);
clReleaseMemObject(xCorrTimeBuffer);
clReleaseMemObject(sqrtWfPwrBuffer);
clReleaseMemObject(CijBuffer);
clReleaseMemObject(MBuffer);
free(voltsFlat);
//free(volts);
//free(recoDelays);
free(intensity_data_r);
free(intensity_data_c);
free(xCorr_data_r);
free(xCorr_data_c);
free(xCorrTime);
free(sqrtWfPwr);
free(Cij);
free(M);
cout<<"Memories deallocated\n";

return maxPixIdx;
}


int reconstructXCorrEnvelopeGetMaxPixAndMap(unsigned int dataType, vector<TGraph *>& cleanEvent, recoEnvData *clEnv,
                float *recoDelays, float *recoDelays_V, float *recoDelays_H,
                int nDir, string pol, const int *chanMask, recoData *summary, char *filename)
{

cout<<"Entered reconstructXCorrEnvelopeGetMaxPixAndMap method\n";
int nSamp;
int nAnt = (int)cleanEvent.size()/2; // Divide by 2 for only one polarization
int unmaskedNChan=0;
for(int ch=0; ch<2*nAnt; ch++) unmaskedNChan+=chanMask[ch];
cout<<"unmaskedNChan: "<<unmaskedNChan<<" nAnt: "<<nAnt<<endl;
float wInt;
if( dataType == 0 ) wInt = 0.5f; //AraSim event
else if( dataType == 1 ){ //real event
if( pol == "vpol" ) wInt = 0.4f;
else if ( pol == "hpol" ) wInt = 0.625f;
else { cerr<<"recoPolType undefined\n"; return -1; }
} else {
cerr<<"dataType undefined\n"; return -1; }

/*
 * Loading voltsFlat array
 */

double t, v;
float *voltsFlat;

if( pol == "vpol" ){

   /* Using the 1st vpol wf for nSamp. Should make sure all vpol channels have the same nSamp. FIXME */
   nSamp = cleanEvent[0]->GetN();
   voltsFlat = (float*)calloc(nAnt*nSamp, sizeof(float));
   if(voltsFlat == NULL){ cerr<<"Null pointer to voltsFlat\n"; return -1; }

   for(int ch=0; ch<nAnt; ch++){
       if( chanMask[ch] ){
         for(int s=0; s<nSamp; s++){
            cleanEvent[ch]->GetPoint(s,t,v);
            /* Bartlett window applied in main analysis code */
            //voltsFlat[nSamp*ch + s] = static_cast<float>(v * FFTtools::bartlettWindow(s, nSamp));
            voltsFlat[nSamp*ch + s] = static_cast<float>(v);
         }
      } else {
         for(int s=0; s<nSamp; s++){
            voltsFlat[nSamp*ch + s] = 0.f;
         }
      }
   }

} else if ( pol == "hpol" ){

   /* Using the 1st hpol wf for nSamp. Should make sure all hpol channels have the same nSamp. FIXME */
   nSamp = cleanEvent[nAnt]->GetN();
   voltsFlat = (float*)calloc(nAnt*nSamp, sizeof(float));
   if(voltsFlat == NULL){ cerr<<"Null pointer to voltsFlat\n"; return -1; }

   for(int ch=0; ch<nAnt; ch++){
       if( chanMask[ch+nAnt] ){
         for(int s=0; s<nSamp; s++){
            cleanEvent[ch+nAnt]->GetPoint(s,t,v);
            /* Bartlett window applied in main analysis code */
            //voltsFlat[nSamp*ch + s] = static_cast<float>(v * FFTtools::bartlettWindow(s, nSamp));
            voltsFlat[nSamp*ch + s] = static_cast<float>(v);
         }
      } else {
         for(int s=0; s<nSamp; s++){
            voltsFlat[nSamp*ch + s] = 0.f;
         }
      }
   }

} else {
   cerr<<"recoPolType not defined\n"; return -1;
}
cout<<"voltsFlat loaded\n";

/*
 * Preparation for OUT_OF_PLACE transforms
 */

int interlvHermOutputSize = 2*(1 + nSamp/2); //Hermitian layout
int interlvOutputSize     = 2*nSamp;         //Not Hermitian layout
int planarHermOutputSize = (1 + nSamp/2);    //Hermitian planar layout

float *intensity_data_r, *intensity_data_c;
intensity_data_r = (float*)calloc(nAnt*planarHermOutputSize, sizeof(float));
intensity_data_c = (float*)calloc(nAnt*planarHermOutputSize, sizeof(float));

/*
 * FFT library related declarations
 */
//clfftPlanHandle planHandle;
clfftDim fftDim = CLFFT_1D;
size_t clLengths[1]   = {nSamp};
size_t clInStride[1]  = {1};
size_t clOutStride[1] = {1};
size_t inDist  = nSamp;
size_t outDist   = planarHermOutputSize;
size_t batchSize = nAnt;

/*
 * Set up clFFT
 */
//clfftSetupData fftSetup;
//err = clfftInitSetupData(&fftSetup);
//err = clfftSetup(&fftSetup);

/*
 * Prepare plan
 */
cout<<"Preparing plan...\n";
int err;
err = clfftCreateDefaultPlan(&clEnv->planHandle, clEnv->context, fftDim, clLengths);
err = clfftSetPlanPrecision(clEnv->planHandle, CLFFT_SINGLE);
err = clfftSetLayout(clEnv->planHandle, CLFFT_REAL, CLFFT_HERMITIAN_PLANAR);
err = clfftSetPlanScale(clEnv->planHandle, CLFFT_FORWARD, 1.f);
err = clfftSetPlanBatchSize(clEnv->planHandle, batchSize);
err = clfftSetPlanInStride(clEnv->planHandle, fftDim, clInStride);
err = clfftSetPlanOutStride(clEnv->planHandle, fftDim, clOutStride);
err = clfftSetPlanDistance(clEnv->planHandle, inDist, outDist);
err = clfftSetResultLocation(clEnv->planHandle, CLFFT_OUTOFPLACE);
err = clfftBakePlan(clEnv->planHandle, 1, &clEnv->queue, NULL, NULL);
cout<<"Plan prepared\n";
/* The plan is now ready to be executed */
cl_mem voltsFlatBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*nSamp*sizeof(float), NULL, &err);
err = clEnqueueWriteBuffer(clEnv->queue, voltsFlatBuffer, CL_TRUE, 0, nAnt*nSamp*sizeof(float), voltsFlat, 0, NULL, NULL);

cl_mem intensityRBuffer= clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*planarHermOutputSize*sizeof(float),
                         NULL, &err);
cl_mem intensityCBuffer= clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*planarHermOutputSize*sizeof(float),
                         NULL, &err);
cl_mem outBuffers[2] = {intensityRBuffer, intensityCBuffer};

cout<<"Enqueueing FFT\n";
err = clfftEnqueueTransform(clEnv->planHandle, CLFFT_FORWARD, 1, &clEnv->queue, 0, NULL, NULL,
                            &voltsFlatBuffer, outBuffers, NULL);
err = clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, outBuffers[0], CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_r,
                    0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, outBuffers[1], CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_c,
                    0, NULL, NULL);
cout<<"FFT done\n";
/*
 * Clean up CLFFT
 */
cout<<"Destroying plan...\n";
err = clfftDestroyPlan(&clEnv->planHandle);
//clfftTeardown();
cout<<"Plan destroyed\n";


#ifdef bandpass
/*
 * Bandpass signals
 */
cout<<"Preparing bandpass filter..."<<endl;
float freqBin = 1e3 / (wInt * (float)nSamp); // wInt in ns. 1e3 for MHz
float lowFreq = 200.f;
float highFreq= 450.f;

clSetKernelArg(clEnv->bandPassFilter, 0, sizeof(cl_mem), &intensityRBuffer);
clSetKernelArg(clEnv->bandPassFilter, 1, sizeof(cl_mem), &intensityCBuffer);
clSetKernelArg(clEnv->bandPassFilter, 2, sizeof(float),  &freqBin);
clSetKernelArg(clEnv->bandPassFilter, 3, sizeof(float),  &lowFreq);
clSetKernelArg(clEnv->bandPassFilter, 4, sizeof(float),  &highFreq);

unsigned int dim = 2;
size_t bandPassWorkSize[2] = {(size_t)nAnt, (size_t)planarHermOutputSize};
clEnqueueNDRangeKernel(clEnv->queue, clEnv->bandPassFilter, dim, NULL, bandPassWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, intensityRBuffer, CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_r,
                    0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, intensityCBuffer, CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_c,
                    0, NULL, NULL);

cout<<"Bandpass filter done"<<endl;
#endif

/*
 * Cross-correlate wfs
 */

int nBaseline = nAnt*nAnt;
cout<<"nBaseline: "<<nBaseline<<endl;

float *xCorr_data_r = (float*)calloc(nBaseline*planarHermOutputSize, sizeof(float));
float *xCorr_data_c = (float*)calloc(nBaseline*planarHermOutputSize, sizeof(float));

cl_mem xCorrRBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                    sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_r, &err);

cl_mem xCorrCBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                    sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_c, &err);

clSetKernelArg(clEnv->xCorrWf, 0, sizeof(cl_mem), &xCorrRBuffer);
clSetKernelArg(clEnv->xCorrWf, 1, sizeof(cl_mem), &xCorrCBuffer);
clSetKernelArg(clEnv->xCorrWf, 2, sizeof(cl_mem), &intensityRBuffer);
clSetKernelArg(clEnv->xCorrWf, 3, sizeof(cl_mem), &intensityCBuffer);
clSetKernelArg(clEnv->xCorrWf, 4, sizeof(int),    &nAnt);

unsigned int workDim = 2;
size_t globalWorkSize[2] = {(size_t)nBaseline, (size_t)planarHermOutputSize};
clEnqueueNDRangeKernel(clEnv->queue, clEnv->xCorrWf, workDim, NULL, globalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, xCorrRBuffer, CL_TRUE, 0, sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_r, 0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, xCorrCBuffer, CL_TRUE, 0, sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_c, 0, NULL, NULL);
cout<<"Done xCorrWf\n";

/*
 * Inverse FFT cross-correlation
 */

inDist  = planarHermOutputSize;
outDist = nSamp;
batchSize = nBaseline;

/*
 * Prepare plan
 */

err = clfftCreateDefaultPlan(&clEnv->planHandle, clEnv->context, fftDim, clLengths);
err = clfftSetPlanPrecision(clEnv->planHandle, CLFFT_SINGLE);
err = clfftSetLayout(clEnv->planHandle, CLFFT_HERMITIAN_PLANAR, CLFFT_REAL);
err = clfftSetPlanScale(clEnv->planHandle, CLFFT_BACKWARD, 1.f/(float)nSamp);
err = clfftSetPlanBatchSize(clEnv->planHandle, batchSize);
err = clfftSetPlanInStride(clEnv->planHandle, fftDim, clInStride);
err = clfftSetPlanOutStride(clEnv->planHandle, fftDim, clOutStride);
err = clfftSetPlanDistance(clEnv->planHandle, inDist, outDist);
err = clfftSetResultLocation(clEnv->planHandle, CLFFT_OUTOFPLACE);
err = clfftBakePlan(clEnv->planHandle, 1, &clEnv->queue, NULL, NULL);

/* The plan is now ready to be executed */

cl_mem inBuffers[2] = {xCorrRBuffer, xCorrCBuffer};

float *xCorrTime = (float*)calloc(nBaseline*nSamp, sizeof(float));
cl_mem xCorrTimeBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nBaseline*nSamp*sizeof(float),
                           NULL, &err);
err = clfftEnqueueTransform(clEnv->planHandle, CLFFT_BACKWARD, 1, &clEnv->queue, 0, NULL, NULL, inBuffers, &xCorrTimeBuffer, NULL);
err = clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, xCorrTimeBuffer, CL_TRUE, 0, sizeof(float)*nBaseline*nSamp, xCorrTime, 0, NULL, NULL);

/*
 * Clean up FFT
 */

err = clfftDestroyPlan(&clEnv->planHandle);
clfftTeardown();
cout<<"Done inverse FFT\n";

/*
 * Get Hilbert transform of XCorr function
 */

TCanvas cvs("cvs","cva",800,600);
cvs.Divide(1,2);


float dt[nSamp];
float xCorrValue[nSamp];
double t_temp, v_temp;
int ant1, ant2;

for(int baseline=0; baseline<nBaseline; baseline++){


   ant1 = baseline / nAnt;
   ant2 = baseline % nAnt;

   for(int s=0; s<nSamp; s++){

   dt[s] = wInt*s;
   xCorrValue[s] = xCorrTime[nSamp*baseline + s];

   }

   TGraph *xCorrGraph = new TGraph(nSamp, dt, xCorrValue);

   //cvs.cd(1);
   //xCorrGraph->Draw("AL");

   TGraph* envelope = FFTtools::getHilbertEnvelope( xCorrGraph );
/*
   if( ant1 == 0 && ant2 == 3){
   cvs.cd(1);
   envelope->Draw("AL");
   //xCorrGraph->Draw("AL");
   cvs.SaveAs("xCorrEnvelope_chan0_3.C");
   }
*/ /*
   if( ant1 == 3 && ant2 == 0){
   cvs.cd(2);
   envelope->Draw("AL");
   }
   */
   //cvs.SaveAs("xCorrEnvelope_chan0_3.C");

   for(int s=0; s<nSamp; s++){

   envelope->GetPoint(s,t_temp,v_temp);
   xCorrTime[nSamp*baseline + s] = static_cast<float>(v_temp);

   }

  delete xCorrGraph;
  delete envelope;

}

cl_mem xCorrEnvBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                       sizeof(float)*nBaseline*nSamp,
                                       xCorrTime, &err);


/*
 * Compute cross-correlation coefficients Cij in each direction
 */

/* Compute the square root of each channel's total wf power */

float *sqrtWfPwr = (float*)calloc(nAnt, sizeof(float));
float pwr=0.f;
for(int ant=0; ant<nAnt; ant++){
   pwr = 0.;
   for(int s=0; s<nSamp; s++){
      pwr += (voltsFlat[ant*nSamp + s] * voltsFlat[ant*nSamp + s]);
   }
   sqrtWfPwr[ant] = sqrt(pwr);
   //sqrtWfPwr[ant] = 1.f;
}
cout<<"Done sqrtWfPwr\n";

cl_mem sqrtWfPwrBuffer  = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt,
                          sqrtWfPwr, &err);
cl_mem recoDelaysBuffer;
if(pol == "vpol" ){
recoDelaysBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt*nDir,
                                  recoDelays_V, &err);
} else if (pol == "hpol" ) {
recoDelaysBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt*nDir,
                                  recoDelays_H, &err);
} else {
cerr<<"recoPolType not defined\n"; return -1;
}
cout<<"recoDelaysBuffer created\n";

float *Cij = (float*)calloc(nDir*nBaseline, sizeof(float));

cl_mem CijBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nDir*nBaseline,
                   Cij, &err);

clSetKernelArg(clEnv->computeXCorrCoef, 0, sizeof(cl_mem), &CijBuffer);
clSetKernelArg(clEnv->computeXCorrCoef, 1, sizeof(cl_mem), &xCorrEnvBuffer);
clSetKernelArg(clEnv->computeXCorrCoef, 2, sizeof(cl_mem), &recoDelaysBuffer);
clSetKernelArg(clEnv->computeXCorrCoef, 3, sizeof(cl_mem), &sqrtWfPwrBuffer);
clSetKernelArg(clEnv->computeXCorrCoef, 4, sizeof(float),  &wInt);
clSetKernelArg(clEnv->computeXCorrCoef, 5, sizeof(int),    &nAnt);
clSetKernelArg(clEnv->computeXCorrCoef, 6, sizeof(int),    &nSamp);

workDim = 2;
size_t CijGlobalWorkSize[2] = {(size_t)nDir, (size_t)nBaseline};

err = clEnqueueNDRangeKernel(clEnv->queue, clEnv->computeXCorrCoef, workDim, NULL, CijGlobalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, CijBuffer, CL_TRUE, 0, sizeof(float)*nDir*nBaseline, Cij, 0, NULL, NULL);

cout<<"Done computeXCorrCoef\n";

/*
 * Sum Cij's of all baselines in each reco direction to obtain coherence M(r-hat)
 */

float *M = (float*)calloc(nDir, sizeof(float));
cl_mem MBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nDir,
                  M, &err);

clSetKernelArg(clEnv->computeCoherence, 0, sizeof(cl_mem), &MBuffer);
clSetKernelArg(clEnv->computeCoherence, 1, sizeof(cl_mem), &CijBuffer);
clSetKernelArg(clEnv->computeCoherence, 2, sizeof(int),    &nBaseline);

workDim = 1;
size_t MGlobalWorkSize = nDir;

clEnqueueNDRangeKernel(clEnv->queue, clEnv->computeCoherence, workDim, NULL, &MGlobalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, MBuffer, CL_TRUE, 0, sizeof(float)*nDir, M, 0, NULL, NULL);

cout<<"Done computeCoherence\n";

/*
 * Loop over M to find the max coherence and its pix index
 */

float max=0.f;
int maxPixIdx;

for(int idx=0; idx<nDir; idx++){
   if(M[idx] > max){
      max = M[idx];
      maxPixIdx = idx;
   }
}

cout<<"max: "<<max<<endl;
summary->setMaxPixInfo(maxPixIdx, max);

/*
 * Write FITS file
 */

cout<<"Creating Healpix map and writing to FITS....\n";
arr<float> MArr = arr<float>(&M[0], (size_t)nDir);
Healpix_Map<float> skyMap = Healpix_Map<float>(MArr, HEALPIX_ORDERING);

fitshandle fitsOut;
//#ifdef CSW
//char filename[] = "testCSWSkyMap.fits";
//#else
//char filename[] = "testXCorrSkyMap.fits";
//#endif
remove(filename);
fitsOut.create(filename);

write_Healpix_map_to_fits(fitsOut, skyMap, PLANCK_FLOAT32);
cout<<"Healpix map written\n";

/*
 * Deallocate memories
 */

cout<<"Deallocating memories...\n";
clReleaseMemObject(recoDelaysBuffer);
clReleaseMemObject(intensityRBuffer);
clReleaseMemObject(intensityCBuffer);
clReleaseMemObject(voltsFlatBuffer);
clReleaseMemObject(xCorrCBuffer);
clReleaseMemObject(xCorrRBuffer);
clReleaseMemObject(xCorrTimeBuffer);
clReleaseMemObject(xCorrEnvBuffer);
clReleaseMemObject(sqrtWfPwrBuffer);
clReleaseMemObject(CijBuffer);
clReleaseMemObject(MBuffer);
free(voltsFlat);
//free(volts);
//free(recoDelays);
free(intensity_data_r);
free(intensity_data_c);
free(xCorr_data_r);
free(xCorr_data_c);
free(xCorrTime);
free(sqrtWfPwr);
free(Cij);
free(M);
cout<<"Memories deallocated\n";

return maxPixIdx;
}
/*
int reconstruct3DXCorrEnvelopeGetMaxPixAndMap(unsigned int dataType, vector<TGraph *>& cleanEvent, recoEnvData *clEnv,
                float *recoDelays, float *recoDelays_V, float *recoDelays_H,
                string pol, const int *chanMask, recoData *summary, char *filename,
                TH1F *xCorrAroundPeakHist[] )
{

cout<<"Entered reconstruct3DXCorrEnvelopeGetMaxPixAndMap method\n";
int nSamp;
int nAnt = (int)cleanEvent.size()/2; // Divide by 2 for only one polarization
int unmaskedNChan=0;
for(int ch=0; ch<2*nAnt; ch++) unmaskedNChan+=chanMask[ch];
cout<<"unmaskedNChan: "<<unmaskedNChan<<" nAnt: "<<nAnt<<endl;
float wInt;
if( dataType == 0 ) wInt = 0.5f; //AraSim event
else if( dataType == 1 ){ //real event
if( pol == "vpol" ) wInt = 0.4f;
else if ( pol == "hpol" ) wInt = 0.625f;
else { cerr<<"recoPolType undefined\n"; return -1; }
} else {
cerr<<"dataType undefined\n"; return -1; }

int nDir   = summary->onion->nDir;
int nLayer = summary->onion->nLayer;
printf("nDir: %d nLayer: %d\n",nDir,nLayer);
*//*
 * Loading voltsFlat array
 */
/*
double t, v;
float *voltsFlat;

if( pol == "vpol" ){
*/
   /* Using the 1st vpol wf for nSamp. Should make sure all vpol channels have the same nSamp. FIXME */
/*   nSamp = cleanEvent[0]->GetN();
   voltsFlat = (float*)calloc(nAnt*nSamp, sizeof(float));
   if(voltsFlat == NULL){ cerr<<"Null pointer to voltsFlat\n"; return -1; }

   for(int ch=0; ch<nAnt; ch++){
       if( chanMask[ch] ){
         for(int s=0; s<nSamp; s++){
            cleanEvent[ch]->GetPoint(s,t,v);
  */          /* Bartlett window applied in main analysis code */
            //voltsFlat[nSamp*ch + s] = static_cast<float>(v * FFTtools::bartlettWindow(s, nSamp));
 /*           voltsFlat[nSamp*ch + s] = static_cast<float>(v);
         }
      } else {
         for(int s=0; s<nSamp; s++){
            voltsFlat[nSamp*ch + s] = 0.f;
         }
      }
   }

} else if ( pol == "hpol" ){

*/   /* Using the 1st hpol wf for nSamp. Should make sure all hpol channels have the same nSamp. FIXME */
 /*  nSamp = cleanEvent[nAnt]->GetN();
   voltsFlat = (float*)calloc(nAnt*nSamp, sizeof(float));
   if(voltsFlat == NULL){ cerr<<"Null pointer to voltsFlat\n"; return -1; }

   for(int ch=0; ch<nAnt; ch++){
       if( chanMask[ch+nAnt] ){
         for(int s=0; s<nSamp; s++){
            cleanEvent[ch+nAnt]->GetPoint(s,t,v);
*/           /* Bartlett window applied in main analysis code */
            //voltsFlat[nSamp*ch + s] = static_cast<float>(v * FFTtools::bartlettWindow(s, nSamp));
/*            voltsFlat[nSamp*ch + s] = static_cast<float>(v);
         }
      } else {
         for(int s=0; s<nSamp; s++){
            voltsFlat[nSamp*ch + s] = 0.f;
         }
      }
   }

} else {
   cerr<<"recoPolType not defined\n"; return -1;
}
cout<<"voltsFlat loaded\n";
*/
/*
 * Preparation for OUT_OF_PLACE transforms
 */
/*
int interlvHermOutputSize = 2*(1 + nSamp/2); //Hermitian layout
int interlvOutputSize     = 2*nSamp;         //Not Hermitian layout
int planarHermOutputSize = (1 + nSamp/2);    //Hermitian planar layout

float *intensity_data_r, *intensity_data_c;
intensity_data_r = (float*)calloc(nAnt*planarHermOutputSize, sizeof(float));
intensity_data_c = (float*)calloc(nAnt*planarHermOutputSize, sizeof(float));
*/
/*
 * FFT library related declarations
 *//*
//clfftPlanHandle planHandle;
clfftDim fftDim = CLFFT_1D;
size_t clLengths[1]   = {nSamp};
size_t clInStride[1]  = {1};
size_t clOutStride[1] = {1};
size_t inDist  = nSamp;
size_t outDist   = planarHermOutputSize;
size_t batchSize = nAnt;
*/
/*
 * Set up clFFT
 */
//clfftSetupData fftSetup;
//err = clfftInitSetupData(&fftSetup);
//err = clfftSetup(&fftSetup);

/*
 * Prepare plan
 *//*
cout<<"Preparing plan...\n";
int err;
err = clfftCreateDefaultPlan(&clEnv->planHandle, clEnv->context, fftDim, clLengths);
err = clfftSetPlanPrecision(clEnv->planHandle, CLFFT_SINGLE);
err = clfftSetLayout(clEnv->planHandle, CLFFT_REAL, CLFFT_HERMITIAN_PLANAR);
err = clfftSetPlanScale(clEnv->planHandle, CLFFT_FORWARD, 1.f);
err = clfftSetPlanBatchSize(clEnv->planHandle, batchSize);
err = clfftSetPlanInStride(clEnv->planHandle, fftDim, clInStride);
err = clfftSetPlanOutStride(clEnv->planHandle, fftDim, clOutStride);
err = clfftSetPlanDistance(clEnv->planHandle, inDist, outDist);
err = clfftSetResultLocation(clEnv->planHandle, CLFFT_OUTOFPLACE);
err = clfftBakePlan(clEnv->planHandle, 1, &clEnv->queue, NULL, NULL);
cout<<"Plan prepared\n";
*//* The plan is now ready to be executed */
/*cl_mem voltsFlatBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*nSamp*sizeof(float), NULL, &err);
err = clEnqueueWriteBuffer(clEnv->queue, voltsFlatBuffer, CL_TRUE, 0, nAnt*nSamp*sizeof(float), voltsFlat, 0, NULL, NULL);

cl_mem intensityRBuffer= clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*planarHermOutputSize*sizeof(float),
                         NULL, &err);
cl_mem intensityCBuffer= clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*planarHermOutputSize*sizeof(float),
                         NULL, &err);
cl_mem outBuffers[2] = {intensityRBuffer, intensityCBuffer};

cout<<"Enqueueing FFT\n";
err = clfftEnqueueTransform(clEnv->planHandle, CLFFT_FORWARD, 1, &clEnv->queue, 0, NULL, NULL,
                            &voltsFlatBuffer, outBuffers, NULL);
err = clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, outBuffers[0], CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_r,
                    0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, outBuffers[1], CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_c,
                    0, NULL, NULL);
cout<<"FFT done\n";
*//*
 * Clean up CLFFT
 *//*
cout<<"Destroying plan...\n";
err = clfftDestroyPlan(&clEnv->planHandle);
//clfftTeardown();
cout<<"Plan destroyed\n";


#ifdef bandpass
*//*
 * Bandpass signals
 *//*
cout<<"Preparing bandpass filter..."<<endl;
float freqBin = 1e3 / (wInt * (float)nSamp); // wInt in ns. 1e3 for MHz
float lowFreq = 200.f;
float highFreq= 450.f;

clSetKernelArg(clEnv->bandPassFilter, 0, sizeof(cl_mem), &intensityRBuffer);
clSetKernelArg(clEnv->bandPassFilter, 1, sizeof(cl_mem), &intensityCBuffer);
clSetKernelArg(clEnv->bandPassFilter, 2, sizeof(float),  &freqBin);
clSetKernelArg(clEnv->bandPassFilter, 3, sizeof(float),  &lowFreq);
clSetKernelArg(clEnv->bandPassFilter, 4, sizeof(float),  &highFreq);

unsigned int dim = 2;
size_t bandPassWorkSize[2] = {(size_t)nAnt, (size_t)planarHermOutputSize};
clEnqueueNDRangeKernel(clEnv->queue, clEnv->bandPassFilter, dim, NULL, bandPassWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, intensityRBuffer, CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_r,
                    0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, intensityCBuffer, CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_c,
                    0, NULL, NULL);

cout<<"Bandpass filter done"<<endl;
#endif
*/
/*
 * Cross-correlate wfs
 */
/*
int nBaseline = nAnt*nAnt;
cout<<"nBaseline: "<<nBaseline<<endl;

float *xCorr_data_r = (float*)calloc(nBaseline*planarHermOutputSize, sizeof(float));
float *xCorr_data_c = (float*)calloc(nBaseline*planarHermOutputSize, sizeof(float));

cl_mem xCorrRBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                    sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_r, &err);

cl_mem xCorrCBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                    sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_c, &err);

clSetKernelArg(clEnv->xCorrWf, 0, sizeof(cl_mem), &xCorrRBuffer);
clSetKernelArg(clEnv->xCorrWf, 1, sizeof(cl_mem), &xCorrCBuffer);
clSetKernelArg(clEnv->xCorrWf, 2, sizeof(cl_mem), &intensityRBuffer);
clSetKernelArg(clEnv->xCorrWf, 3, sizeof(cl_mem), &intensityCBuffer);
clSetKernelArg(clEnv->xCorrWf, 4, sizeof(int),    &nAnt);

unsigned int workDim = 2;
size_t globalWorkSize[2] = {(size_t)nBaseline, (size_t)planarHermOutputSize};
clEnqueueNDRangeKernel(clEnv->queue, clEnv->xCorrWf, workDim, NULL, globalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, xCorrRBuffer, CL_TRUE, 0, sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_r, 0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, xCorrCBuffer, CL_TRUE, 0, sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_c, 0, NULL, NULL);
cout<<"Done xCorrWf\n";
*/
/*
 * Inverse FFT cross-correlation
 */
/*
inDist  = planarHermOutputSize;
outDist = nSamp;
batchSize = nBaseline;
*/
/*
 * Prepare plan
 */
/*
err = clfftCreateDefaultPlan(&clEnv->planHandle, clEnv->context, fftDim, clLengths);
err = clfftSetPlanPrecision(clEnv->planHandle, CLFFT_SINGLE);
err = clfftSetLayout(clEnv->planHandle, CLFFT_HERMITIAN_PLANAR, CLFFT_REAL);
err = clfftSetPlanScale(clEnv->planHandle, CLFFT_BACKWARD, 1.f/(float)nSamp);
err = clfftSetPlanBatchSize(clEnv->planHandle, batchSize);
err = clfftSetPlanInStride(clEnv->planHandle, fftDim, clInStride);
err = clfftSetPlanOutStride(clEnv->planHandle, fftDim, clOutStride);
err = clfftSetPlanDistance(clEnv->planHandle, inDist, outDist);
err = clfftSetResultLocation(clEnv->planHandle, CLFFT_OUTOFPLACE);
err = clfftBakePlan(clEnv->planHandle, 1, &clEnv->queue, NULL, NULL);
*/
/* The plan is now ready to be executed */
/*
cl_mem inBuffers[2] = {xCorrRBuffer, xCorrCBuffer};

float *xCorrTime = (float*)calloc(nBaseline*nSamp, sizeof(float));
cl_mem xCorrTimeBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nBaseline*nSamp*sizeof(float),
                           NULL, &err);
err = clfftEnqueueTransform(clEnv->planHandle, CLFFT_BACKWARD, 1, &clEnv->queue, 0, NULL, NULL, inBuffers, &xCorrTimeBuffer, NULL);
err = clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, xCorrTimeBuffer, CL_TRUE, 0, sizeof(float)*nBaseline*nSamp, xCorrTime, 0, NULL, NULL);
*/
/*
 * Clean up FFT
 */
/*
err = clfftDestroyPlan(&clEnv->planHandle);
clfftTeardown();
cout<<"Done inverse FFT\n";
*/
/*
 * Get Hilbert transform of XCorr function
 */
/*
TCanvas cvs("cvs","cva",800,600);
cvs.Divide(1,2);


float dt[nSamp];
float xCorrValue[nSamp];
double t_temp, v_temp;
double t_temp2, v_temp2;
int ant1, ant2;
float plusMinusTime = 25.f; //record +-25 ns around peak in XCorr function
//int plusMinusRange = (int)(plusMinusTime / wInt);
//#ifndef XCORRSUMGRAPH
static TGraph *sillygr = new TGraph();
//TGraph *envelopeSum = new TGraph();
//#define XCORRSUMGRAPH
//#endif

for(int baseline=0; baseline<nBaseline; baseline++){


   ant1 = baseline / nAnt;
   ant2 = baseline % nAnt;

   for(int s=0; s<nSamp; s++){

   dt[s] = wInt*s;
   xCorrValue[s] = xCorrTime[nSamp*baseline + s];

   }

   TGraph *xCorrGraph = new TGraph(nSamp, dt, xCorrValue);

   //cvs.cd(1);
   //xCorrGraph->Draw("AL");

   TGraph* envelope = FFTtools::getHilbertEnvelope( xCorrGraph );
*//*
   if( ant1 == 0 && ant2 == 3){
   cvs.cd(1);
   envelope->Draw("AL");
   //xCorrGraph->Draw("AL");
   cvs.SaveAs("xCorrEnvelope_chan0_3.C");
   }
*/
/*
   if( ant1 == 3 && ant2 == 0){
   cvs.cd(2);
   envelope->Draw("AL");
   }

   cvs.SaveAs("xCorrEnvelope_chan0_3.C");
*/
/*
   if( ant1 == 0 && chanMask[ant1] ){
      cout<<"recording XCorr values around peak of each baseline\n";

      if( ant2 == 1 && chanMask[ant2] ){
      stackXCorrAroundPeak(envelope, xCorrAroundPeakHist[0], plusMinusTime);
      } else if
      ( ant2 == 3 && chanMask[ant2] ){
      stackXCorrAroundPeak(envelope, xCorrAroundPeakHist[1], plusMinusTime);
      } else if
      ( ant2 == 7 && chanMask[ant2] ){
      stackXCorrAroundPeak(envelope, xCorrAroundPeakHist[2], plusMinusTime);
      }
   }

   for(int s=0; s<nSamp; s++){

   if(ant1 == 0  && ant2 == 4 ){
   xCorrGraph->GetPoint(s,t_temp,v_temp);
   if(sillygr->GetN() == 0){
      //if(envelopeSum->GetN() != 0 ) cerr<<"Both TGraphs should be uninitialized!"<<endl;

      sillygr->SetPoint(s, t_temp, v_temp);

   } else {

      sillygr->GetPoint(s, t_temp2, v_temp2);
      sillygr->SetPoint(s, t_temp+t_temp2, v_temp+v_temp2);

   }
   }
   envelope->GetPoint(s,t_temp,v_temp);
   xCorrTime[nSamp*baseline + s] = static_cast<float>(v_temp);

   }

  cvs.cd();
  sillygr->Draw("AL");
  cvs.SaveAs("xCorrSumGraph_chan0_4.C");

  delete xCorrGraph;
  delete envelope;

}

cl_mem xCorrEnvBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                       sizeof(float)*nBaseline*nSamp,
                                       xCorrTime, &err);

*/
/*
 * Compute cross-correlation coefficients Cij in each direction
 */

/* Compute the square root of each channel's total wf power */
/*
float *sqrtWfPwr = (float*)calloc(nAnt, sizeof(float));
float pwr=0.f;
for(int ant=0; ant<nAnt; ant++){
   pwr = 0.;
   for(int s=0; s<nSamp; s++){
      pwr += (voltsFlat[ant*nSamp + s] * voltsFlat[ant*nSamp + s]);
   }
   sqrtWfPwr[ant] = sqrt(pwr);
   //sqrtWfPwr[ant] = 1.f;
}
cout<<"Done sqrtWfPwr\n";

cl_mem sqrtWfPwrBuffer  = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt,
                          sqrtWfPwr, &err);
cl_mem recoDelaysBuffer;
if(pol == "vpol" ){
recoDelaysBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt*nDir*nLayer,
                                  recoDelays_V, &err);
} else if (pol == "hpol" ) {
recoDelaysBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt*nDir*nLayer,
                                  recoDelays_H, &err);
} else {
cerr<<"recoPolType not defined\n"; return -1;
}
cout<<"recoDelaysBuffer created\n";

float *Cij = (float*)calloc(nLayer*nDir*nBaseline, sizeof(float));

cl_mem CijBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nLayer*nDir*nBaseline,
                   Cij, &err);

clSetKernelArg(clEnv->computeXCorrCoef, 0, sizeof(cl_mem), &CijBuffer);
clSetKernelArg(clEnv->computeXCorrCoef, 1, sizeof(cl_mem), &xCorrEnvBuffer);
clSetKernelArg(clEnv->computeXCorrCoef, 2, sizeof(cl_mem), &recoDelaysBuffer);
clSetKernelArg(clEnv->computeXCorrCoef, 3, sizeof(cl_mem), &sqrtWfPwrBuffer);
clSetKernelArg(clEnv->computeXCorrCoef, 4, sizeof(float),  &wInt);
clSetKernelArg(clEnv->computeXCorrCoef, 5, sizeof(int),    &nAnt);
clSetKernelArg(clEnv->computeXCorrCoef, 6, sizeof(int),    &nSamp);

workDim = 3;
size_t CijGlobalWorkSize[3] = {(size_t)nLayer, (size_t)nDir, (size_t)nBaseline};

err = clEnqueueNDRangeKernel(clEnv->queue, clEnv->computeXCorrCoef, workDim, NULL, CijGlobalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, CijBuffer, CL_TRUE, 0, sizeof(float)*nLayer*nDir*nBaseline, Cij, 0, NULL, NULL);

cout<<"Done computeXCorrCoef\n";
*/
/*
 * Sum Cij's of all baselines in each reco direction to obtain coherence M(r-hat)
 */
/*
float *M = (float*)calloc(nLayer*nDir, sizeof(float));
cl_mem MBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nLayer*nDir,
                  M, &err);

clSetKernelArg(clEnv->computeCoherence, 0, sizeof(cl_mem), &MBuffer);
clSetKernelArg(clEnv->computeCoherence, 1, sizeof(cl_mem), &CijBuffer);
clSetKernelArg(clEnv->computeCoherence, 2, sizeof(int),    &nBaseline);

workDim = 2;
size_t MGlobalWorkSize[2] = {(size_t)nLayer, (size_t)nDir};

clEnqueueNDRangeKernel(clEnv->queue, clEnv->computeCoherence, workDim, NULL, MGlobalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, MBuffer, CL_TRUE, 0, sizeof(float)*nLayer*nDir, M, 0, NULL, NULL);

cout<<"Done computeCoherence\n";
*/
/*
 * Loop over M to find the max coherence and its pix index
 */
/*
float max=0.f;
int maxPixIdx;
int *rank = (int*)calloc(nLayer*nDir, sizeof(int));

TMath::Sort(nLayer*nDir, M, rank);

int *topMaxPixIdx = (int*)calloc(summary->topN, sizeof(int));
float *topMaxPixCoherence = (float*)calloc(summary->topN, sizeof(float));

for(int i=0; i<summary->topN; i++){
   topMaxPixIdx[i] = rank[i];
   topMaxPixCoherence[i] = M[rank[i]];
}

summary->setTopMaxPixInfo(topMaxPixIdx, topMaxPixCoherence);

//int *rankEachLayer = (int*)calloc(nDir, sizeof(int));
int *maxPixIdxEachLayer = (int*)calloc(nLayer, sizeof(int));
float *maxPixCoherenceEachLayer = (float*)calloc(nLayer, sizeof(float));

cl_mem maxPixIdxEachLayerBuffer       = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(int)*nLayer,
                                        maxPixIdxEachLayer, &err);
cl_mem maxPixCoherenceEachLayerBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nLayer,
                                        maxPixCoherenceEachLayer, &err);

clSetKernelArg(clEnv->getMaxPixInfoEachLayer, 0, sizeof(cl_mem), &maxPixIdxEachLayerBuffer);
clSetKernelArg(clEnv->getMaxPixInfoEachLayer, 1, sizeof(cl_mem), &maxPixCoherenceEachLayerBuffer);
clSetKernelArg(clEnv->getMaxPixInfoEachLayer, 2, sizeof(cl_mem), &MBuffer);
clSetKernelArg(clEnv->getMaxPixInfoEachLayer, 3, sizeof(int),    &nDir);

workDim = 1;
size_t maxPixInfoGlobalWorkSize = nLayer;

clEnqueueNDRangeKernel(clEnv->queue, clEnv->getMaxPixInfoEachLayer, workDim, NULL, &maxPixInfoGlobalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, maxPixIdxEachLayerBuffer, CL_TRUE, 0, sizeof(int)*nLayer, maxPixIdxEachLayer, 0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, maxPixCoherenceEachLayerBuffer, CL_TRUE, 0, sizeof(float)*nLayer, maxPixCoherenceEachLayer, 0, NULL, NULL);

summary->setMaxPixInfoEachLayer(maxPixIdxEachLayer, maxPixCoherenceEachLayer);
*/
/*
for(int idx=0; idx<nDir; idx++){
   if(M[idx] > max){
      max = M[idx];
      maxPixIdx = idx;
   }
}
*/
/*
maxPixIdx = rank[0];
cout<<"max: "<<M[rank[0]]<<endl;
summary->setMaxPixInfo(rank[0], M[rank[0]]);
*/
/*
 * Compute likelihood and p value of skymap compared to reference map fits
 */
/*
cout<<"Computing map likelihood and p-value w.r.t. the referenc map..."<<endl;

double pValue, likelihood;
err = computeMapLikelihoodAndPValue(summary->onion->nDir, summary->onion->nLayer, REFERENCE_MAP_FIT_FILE, M, likelihood, pValue);

TH1F *likelihoodDist =  new TH1F("likelihoodDist","likelihoodDist",1000, -1000,1000);
TH1F *pValueDist     =  new TH1F("pValueDist","pValueDist",1000, -1000,1000);

likelihoodDist->Fill(likelihood);
pValueDist->Fill(pValue);

likelihoodDist->Draw();
cvs.SaveAs("testLikelihoodDist.C");
pValueDist->Draw();
cvs.SaveAs("testPValueDist.C");
*/
/*
 * Write FITS file
 */
/*
cout<<"Creating Healpix map and writing to FITS....\n";
//arr<float> MArr = arr<float>(&M[0], (size_t)nDir);
arr<float> MArr = arr<float>(&M[rank[0] / nDir], (size_t)nDir);
Healpix_Map<float> skyMap = Healpix_Map<float>(MArr, RING);

fitshandle fitsOut;
//#ifdef CSW
//char filename[] = "testCSWSkyMap.fits";
//#else
//char filename[] = "testXCorrSkyMap.fits";
//#endif
remove(filename);
fitsOut.create(filename);

write_Healpix_map_to_fits(fitsOut, skyMap, PLANCK_FLOAT32);
cout<<"Healpix map written\n";
*/
/*
arr<float> MEachLayerArr;

for(int i=0; i<nLayer; i++){

   MEachLayerArr = arr<float>(&M[i*nDir], (size_t)nDir);
   skyMap = Healpix_Map<float>(MEachLayerArr, RING);

   char layerMap[200];
   sprintf(layerMap,"layer_%d_skymap.fits",i);
   remove(layerMap);
   fitsOut.create(layerMap);
   write_Healpix_map_to_fits(fitsOut, skyMap, PLANCK_FLOAT32);

}
*/

/*
 * Deallocate memories
 */
/*
cout<<"Deallocating memories...\n";
clReleaseMemObject(recoDelaysBuffer);
clReleaseMemObject(intensityRBuffer);
clReleaseMemObject(intensityCBuffer);
clReleaseMemObject(voltsFlatBuffer);
clReleaseMemObject(xCorrCBuffer);
clReleaseMemObject(xCorrRBuffer);
clReleaseMemObject(xCorrTimeBuffer);
clReleaseMemObject(xCorrEnvBuffer);
clReleaseMemObject(sqrtWfPwrBuffer);
clReleaseMemObject(CijBuffer);
clReleaseMemObject(MBuffer);
clReleaseMemObject(maxPixIdxEachLayerBuffer);
clReleaseMemObject(maxPixCoherenceEachLayerBuffer);
free(voltsFlat);
//free(volts);
//free(recoDelays);
free(intensity_data_r);
free(intensity_data_c);
free(xCorr_data_r);
free(xCorr_data_c);
free(xCorrTime);
free(sqrtWfPwr);
free(Cij);
free(M);
free(rank);
free(topMaxPixIdx);
free(topMaxPixCoherence);
free(maxPixIdxEachLayer);
free(maxPixCoherenceEachLayer);
cout<<"Memories deallocated\n";

return maxPixIdx;
}
*/
/*
int reconstruct3DXCorrEnvelopeGetMaxPixAndMapData(unsigned int dataType, vector<TGraph *>& cleanEvent, recoEnvData *clEnv,
                float *recoDelays, float *recoDelays_V, float *recoDelays_H,
                string pol, const int *chanMask, recoData *summary, char *filename, float *mapData,
                TH1F *xCorrAroundPeakHist[], TGraph *sillygr[] )
{

cout<<"Entered reconstruct3DXCorrEnvelopeGetMaxPixAndMapData method\n";
int nSamp;
//int nAnt = (int)cleanEvent.size()/2; // Divide by 2 for only one polarization
int nAnt;// = (int)cleanEvent.size();
if( pol == "vpol" || pol == "hpol" ) nAnt = (int)cleanEvent.size()/2; else  nAnt = (int)cleanEvent.size();
int unmaskedNChan=0;
for(int ch=0; ch<(int)cleanEvent.size(); ch++) unmaskedNChan+=chanMask[ch];
cout<<"unmaskedNChan: "<<unmaskedNChan<<" nAnt: "<<nAnt<<endl;
float wInt;
if( dataType == 0 ) wInt = 0.5f; //AraSim event
else if( dataType == 1 ){ //real event
if( pol == "vpol" ) wInt = 0.4f;
else if ( pol == "hpol" ) wInt = 0.625f;
else if ( pol == "both" ) wInt = 0.5f;
else { cerr<<"recoPolType undefined\n"; return -1; }
} else {
cerr<<"dataType undefined\n"; return -1; }

int nDir   = summary->onion->nDir;
int nLayer = summary->onion->nLayer;
printf("nDir: %d nLayer: %d\n",nDir,nLayer);
*//*
 * Loading voltsFlat array
 */
/*
double t, v;
float *voltsFlat;

if( pol == "vpol" ){
*/
   /* Using the 1st vpol wf for nSamp. Should make sure all vpol channels have the same nSamp. FIXME */
/*   nSamp = cleanEvent[0]->GetN();
   voltsFlat = (float*)calloc(nAnt*nSamp, sizeof(float));
   if(voltsFlat == NULL){ cerr<<"Null pointer to voltsFlat\n"; return -1; }

   for(int ch=0; ch<nAnt; ch++){
       if( chanMask[ch] ){
         for(int s=0; s<nSamp; s++){
            cleanEvent[ch]->GetPoint(s,t,v);
*/            /* Bartlett window applied in main analysis code */
            //voltsFlat[nSamp*ch + s] = static_cast<float>(v * FFTtools::bartlettWindow(s, nSamp));
/*            voltsFlat[nSamp*ch + s] = static_cast<float>(v);
         }
      } else {
         for(int s=0; s<nSamp; s++){
            voltsFlat[nSamp*ch + s] = 0.f;
         }
      }
   }

} else if ( pol == "hpol" ){

  */ /* Using the 1st hpol wf for nSamp. Should make sure all hpol channels have the same nSamp. FIXME */
/*   nSamp = cleanEvent[nAnt]->GetN();
   voltsFlat = (float*)calloc(nAnt*nSamp, sizeof(float));
   if(voltsFlat == NULL){ cerr<<"Null pointer to voltsFlat\n"; return -1; }

   for(int ch=0; ch<nAnt; ch++){
       if( chanMask[ch+nAnt] ){
         for(int s=0; s<nSamp; s++){
            cleanEvent[ch+nAnt]->GetPoint(s,t,v);
*/          /* Bartlett window applied in main analysis code */
            //voltsFlat[nSamp*ch + s] = static_cast<float>(v * FFTtools::bartlettWindow(s, nSamp));
/*            voltsFlat[nSamp*ch + s] = static_cast<float>(v);
         }
      } else {
         for(int s=0; s<nSamp; s++){
            voltsFlat[nSamp*ch + s] = 0.f;
         }
      }
   }

} else if ( pol == "both" ){

  */ /* Using the 1st vpol wf for nSamp. Should make sure all hpol channels have the same nSamp. FIXME */
/*   nSamp = cleanEvent[0]->GetN();
   voltsFlat = (float*)calloc(nAnt*nSamp, sizeof(float));
   if(voltsFlat == NULL){ cerr<<"Null pointer to voltsFlat\n"; return -1; }

   for(int ch=0; ch<nAnt; ch++){
       if( chanMask[ch] ){
         for(int s=0; s<nSamp; s++){
            cleanEvent[ch]->GetPoint(s,t,v);
*/            /* Bartlett window applied in main analysis code */
            //voltsFlat[nSamp*ch + s] = static_cast<float>(v * FFTtools::bartlettWindow(s, nSamp));
/*            voltsFlat[nSamp*ch + s] = static_cast<float>(v);
         }
      } else {
         for(int s=0; s<nSamp; s++){
            voltsFlat[nSamp*ch + s] = 0.f;
         }
      }
   }

} else {
   cerr<<"recoPolType not defined\n"; return -1;
}
cout<<"voltsFlat loaded\n";
*/
/*
 * Preparation for OUT_OF_PLACE transforms
 */
/*
int interlvHermOutputSize = 2*(1 + nSamp/2); //Hermitian layout
int interlvOutputSize     = 2*nSamp;         //Not Hermitian layout
int planarHermOutputSize = (1 + nSamp/2);    //Hermitian planar layout

float *intensity_data_r, *intensity_data_c;
intensity_data_r = (float*)calloc(nAnt*planarHermOutputSize, sizeof(float));
intensity_data_c = (float*)calloc(nAnt*planarHermOutputSize, sizeof(float));
*/
/*
 * FFT library related declarations
 *//*
//clfftPlanHandle planHandle;
clfftDim fftDim = CLFFT_1D;
size_t clLengths[1]   = {nSamp};
size_t clInStride[1]  = {1};
size_t clOutStride[1] = {1};
size_t inDist  = nSamp;
size_t outDist   = planarHermOutputSize;
size_t batchSize = nAnt;
*/
/*
 * Set up clFFT
 */
//clfftSetupData fftSetup;
//err = clfftInitSetupData(&fftSetup);
//err = clfftSetup(&fftSetup);

/*
 * Prepare plan
 *//*
cout<<"Preparing plan...\n";
int err;
err = clfftCreateDefaultPlan(&clEnv->planHandle, clEnv->context, fftDim, clLengths);
err = clfftSetPlanPrecision(clEnv->planHandle, CLFFT_SINGLE);
err = clfftSetLayout(clEnv->planHandle, CLFFT_REAL, CLFFT_HERMITIAN_PLANAR);
err = clfftSetPlanScale(clEnv->planHandle, CLFFT_FORWARD, 1.f);
err = clfftSetPlanBatchSize(clEnv->planHandle, batchSize);
err = clfftSetPlanInStride(clEnv->planHandle, fftDim, clInStride);
err = clfftSetPlanOutStride(clEnv->planHandle, fftDim, clOutStride);
err = clfftSetPlanDistance(clEnv->planHandle, inDist, outDist);
err = clfftSetResultLocation(clEnv->planHandle, CLFFT_OUTOFPLACE);
err = clfftBakePlan(clEnv->planHandle, 1, &clEnv->queue, NULL, NULL);
cout<<"Plan prepared\n";
*//* The plan is now ready to be executed */
/*
cl_mem voltsFlatBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*nSamp*sizeof(float), NULL, &err);
err = clEnqueueWriteBuffer(clEnv->queue, voltsFlatBuffer, CL_TRUE, 0, nAnt*nSamp*sizeof(float), voltsFlat, 0, NULL, NULL);

cl_mem intensityRBuffer= clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*planarHermOutputSize*sizeof(float),
                         NULL, &err);
cl_mem intensityCBuffer= clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*planarHermOutputSize*sizeof(float),
                         NULL, &err);
cl_mem outBuffers[2] = {intensityRBuffer, intensityCBuffer};

cout<<"Enqueueing FFT\n";
err = clfftEnqueueTransform(clEnv->planHandle, CLFFT_FORWARD, 1, &clEnv->queue, 0, NULL, NULL,
                            &voltsFlatBuffer, outBuffers, NULL);
err = clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, outBuffers[0], CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_r,
                    0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, outBuffers[1], CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_c,
                    0, NULL, NULL);
cout<<"FFT done\n";
*/
/*
 * Clean up CLFFT
 *//*
cout<<"Destroying plan...\n";
err = clfftDestroyPlan(&clEnv->planHandle);
//clfftTeardown();
cout<<"Plan destroyed\n";


#ifdef bandpass
*//*
 * Bandpass signals
 *//*
cout<<"Preparing bandpass filter..."<<endl;
float freqBin = 1e3 / (wInt * (float)nSamp); // wInt in ns. 1e3 for MHz
float lowFreq = 200.f;
float highFreq= 450.f;

clSetKernelArg(clEnv->bandPassFilter, 0, sizeof(cl_mem), &intensityRBuffer);
clSetKernelArg(clEnv->bandPassFilter, 1, sizeof(cl_mem), &intensityCBuffer);
clSetKernelArg(clEnv->bandPassFilter, 2, sizeof(float),  &freqBin);
clSetKernelArg(clEnv->bandPassFilter, 3, sizeof(float),  &lowFreq);
clSetKernelArg(clEnv->bandPassFilter, 4, sizeof(float),  &highFreq);

unsigned int dim = 2;
size_t bandPassWorkSize[2] = {(size_t)nAnt, (size_t)planarHermOutputSize};
clEnqueueNDRangeKernel(clEnv->queue, clEnv->bandPassFilter, dim, NULL, bandPassWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, intensityRBuffer, CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_r,
                    0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, intensityCBuffer, CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_c,
                    0, NULL, NULL);

cout<<"Bandpass filter done"<<endl;
#endif
*/
/*
 * Cross-correlate wfs
 */
/*
const int nBaseline = nAnt*nAnt;
cout<<"nBaseline: "<<nBaseline<<endl;

float *xCorr_data_r = (float*)calloc(nBaseline*planarHermOutputSize, sizeof(float));
float *xCorr_data_c = (float*)calloc(nBaseline*planarHermOutputSize, sizeof(float));

cl_mem xCorrRBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                    sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_r, &err);

cl_mem xCorrCBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                    sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_c, &err);

clSetKernelArg(clEnv->xCorrWf, 0, sizeof(cl_mem), &xCorrRBuffer);
clSetKernelArg(clEnv->xCorrWf, 1, sizeof(cl_mem), &xCorrCBuffer);
clSetKernelArg(clEnv->xCorrWf, 2, sizeof(cl_mem), &intensityRBuffer);
clSetKernelArg(clEnv->xCorrWf, 3, sizeof(cl_mem), &intensityCBuffer);
clSetKernelArg(clEnv->xCorrWf, 4, sizeof(int),    &nAnt);

unsigned int workDim = 2;
size_t globalWorkSize[2] = {(size_t)nBaseline, (size_t)planarHermOutputSize};
clEnqueueNDRangeKernel(clEnv->queue, clEnv->xCorrWf, workDim, NULL, globalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, xCorrRBuffer, CL_TRUE, 0, sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_r, 0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, xCorrCBuffer, CL_TRUE, 0, sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_c, 0, NULL, NULL);
cout<<"Done xCorrWf\n";
*/
/*
 * Inverse FFT cross-correlation
 */
/*
inDist  = planarHermOutputSize;
outDist = nSamp;
batchSize = nBaseline;
*/
/*
 * Prepare plan
 */
/*
err = clfftCreateDefaultPlan(&clEnv->planHandle, clEnv->context, fftDim, clLengths);
err = clfftSetPlanPrecision(clEnv->planHandle, CLFFT_SINGLE);
err = clfftSetLayout(clEnv->planHandle, CLFFT_HERMITIAN_PLANAR, CLFFT_REAL);
err = clfftSetPlanScale(clEnv->planHandle, CLFFT_BACKWARD, 1.f/(float)nSamp);
err = clfftSetPlanBatchSize(clEnv->planHandle, batchSize);
err = clfftSetPlanInStride(clEnv->planHandle, fftDim, clInStride);
err = clfftSetPlanOutStride(clEnv->planHandle, fftDim, clOutStride);
err = clfftSetPlanDistance(clEnv->planHandle, inDist, outDist);
err = clfftSetResultLocation(clEnv->planHandle, CLFFT_OUTOFPLACE);
err = clfftBakePlan(clEnv->planHandle, 1, &clEnv->queue, NULL, NULL);
*/
/* The plan is now ready to be executed */
/*
cl_mem inBuffers[2] = {xCorrRBuffer, xCorrCBuffer};

float *xCorrTime = (float*)calloc(nBaseline*nSamp, sizeof(float));
cl_mem xCorrTimeBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nBaseline*nSamp*sizeof(float),
                           NULL, &err);
err = clfftEnqueueTransform(clEnv->planHandle, CLFFT_BACKWARD, 1, &clEnv->queue, 0, NULL, NULL, inBuffers, &xCorrTimeBuffer, NULL);
err = clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, xCorrTimeBuffer, CL_TRUE, 0, sizeof(float)*nBaseline*nSamp, xCorrTime, 0, NULL, NULL);
*/
/*
 * Clean up FFT
 */
/*
err = clfftDestroyPlan(&clEnv->planHandle);
clfftTeardown();
cout<<"Done inverse FFT\n";
*/
/*
 * Get Hilbert transform of XCorr function
 */
/*
TCanvas cvs("cvs","cvs",800,600);
//cvs.Divide(1,2);
char envelopename[200];

float dt[nSamp];
float xCorrValue[nSamp];
double t_temp, v_temp;
double t_temp2, v_temp2;
int ant1, ant2;
float plusMinusTime = 25.f; //record +-25 ns around peak in XCorr function
//int plusMinusRange = (int)(plusMinusTime / wInt);
//static TGraph *sillygr[64];// = new TGraph();
//static TGraph *sillygr = (TGraph*)malloc(nBaseline*sizeof(TGraph));
*//*
for(int i=0; i<nBaseline; i++){

   sillygr[i]=new TGraph();

}
*/
//static int lock;
//char sillygrname[200];
/*
for(int baseline=0; baseline<nBaseline; baseline++){


   ant1 = baseline / nAnt;
   ant2 = baseline % nAnt;

   for(int s=0; s<nSamp; s++){

   dt[s] = wInt*s;
   xCorrValue[s] = xCorrTime[nSamp*baseline + s];

   }

   TGraph *xCorrGraph = new TGraph(nSamp, dt, xCorrValue);

   //cvs.cd(1);
   //xCorrGraph->Draw("AL");

   TGraph* envelope = FFTtools::getHilbertEnvelope( xCorrGraph );
*//*
   sprintf(envelopename,"xCorrEnvelope_2014_A3_burn_RF_chan%d_%d.C", ant1, ant2);
   cvs.cd();
   //envelope->Draw("AL");
   envelope->SetLineColor(kRed);
   xCorrGraph->Draw("AL");
   envelope->Draw("Lsame");
   cvs.SaveAs(envelopename);
*/
/*
   if( ant1 == 0 && ant2 == 3){
   cvs.cd(1);
   envelope->Draw("AL");
   //xCorrGraph->Draw("AL");
   cvs.SaveAs("xCorrEnvelope_chan0_3.C");
   }
*/
/*
   if( ant1 == 3 && ant2 == 0){
   cvs.cd(2);
   envelope->Draw("AL");
   }

   cvs.SaveAs("xCorrEnvelope_chan0_3.C");
*/
/*
   if( ant1 == 0 && chanMask[ant1] ){
      cout<<"recording XCorr values around peak of each baseline\n";

      if( ant2 == 1 && chanMask[ant2] ){
      stackXCorrAroundPeak(envelope, xCorrAroundPeakHist[0], plusMinusTime);
      } else if
      ( ant2 == 3 && chanMask[ant2] ){
      stackXCorrAroundPeak(envelope, xCorrAroundPeakHist[1], plusMinusTime);
      } else if
      ( ant2 == 7 && chanMask[ant2] ){
      stackXCorrAroundPeak(envelope, xCorrAroundPeakHist[2], plusMinusTime);
      }
   }

   //if(ant1==1 && ant2==5){
   //if(ant1 < ant2 ){
   //for(int s=0; s<nSamp; s++){

   //if(ant1 == 0  && ant2 == 4 ){
   //xCorrGraph->GetPoint(s,t_temp,v_temp);
   //if(lock == 0){
   //if(envelopeSum->GetN() != 0 ) cerr<<"Both TGraphs should be uninitialized!"<<endl;
*//*
   if(sillygr[baseline]->GetN() == 0){

      for(int s=0; s<nSamp; s++){
      xCorrGraph->GetPoint(s,t_temp,v_temp);
      sillygr[baseline]->SetPoint(s, t_temp, v_temp);
      }

   } else {

      for(int s=0; s<nSamp; s++){
      xCorrGraph->GetPoint(s,t_temp,v_temp);
      sillygr[baseline]->GetPoint(s, t_temp2, v_temp2);
      if(t_temp != t_temp2) cerr<<"t_temp != t_temp2 !!\n";
      sillygr[baseline]->SetPoint(s, t_temp2, v_temp+v_temp2);
      }
   }
*/
   //}
/*
   for(int s=0; s<nSamp; s++){

   envelope->GetPoint(s,t_temp,v_temp);
   xCorrTime[nSamp*baseline + s] = static_cast<float>(v_temp);

   }
   //if(lock == 0) lock=1;
*//*
  cvs.cd();
  sprintf(sillygrname,"xCorrSumGraph_chan%d_%d.C",ant1,ant2);
  sillygr[baseline]->Draw("AL");
  cvs.SaveAs(sillygrname);
*/
/*
  delete xCorrGraph;
  delete envelope;

}

cl_mem xCorrEnvBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                       sizeof(float)*nBaseline*nSamp,
                                       xCorrTime, &err);

*/
/*
 * Compute cross-correlation coefficients Cij in each direction
 */

/* Compute the square root of each channel's total wf power */
/*
float *sqrtWfPwr = (float*)calloc(nAnt, sizeof(float));
float pwr=0.f;
for(int ant=0; ant<nAnt; ant++){
   pwr = 0.;
   for(int s=0; s<nSamp; s++){
      pwr += (voltsFlat[ant*nSamp + s] * voltsFlat[ant*nSamp + s]);
   }
   sqrtWfPwr[ant] = sqrt(pwr);
   //sqrtWfPwr[ant] = 1.f;
}
cout<<"Done sqrtWfPwr\n";

cl_mem sqrtWfPwrBuffer  = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt,
                          sqrtWfPwr, &err);
cl_mem recoDelaysBuffer;
if(pol == "vpol" ){
recoDelaysBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt*nDir*nLayer,
                                  recoDelays_V, &err);
} else if (pol == "hpol" ) {
recoDelaysBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt*nDir*nLayer,
                                  recoDelays_H, &err);
} else if (pol == "both" ) {
recoDelaysBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt*nDir*nLayer,
                                  recoDelays, &err);
} else {
cerr<<"recoPolType not defined\n"; return -1;
}
cout<<"recoDelaysBuffer created\n";

float *Cij = (float*)calloc(nLayer*nDir*nBaseline, sizeof(float));

cl_mem CijBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nLayer*nDir*nBaseline,
                   Cij, &err);

clSetKernelArg(clEnv->computeXCorrCoef, 0, sizeof(cl_mem), &CijBuffer);
clSetKernelArg(clEnv->computeXCorrCoef, 1, sizeof(cl_mem), &xCorrEnvBuffer);
clSetKernelArg(clEnv->computeXCorrCoef, 2, sizeof(cl_mem), &recoDelaysBuffer);
clSetKernelArg(clEnv->computeXCorrCoef, 3, sizeof(cl_mem), &sqrtWfPwrBuffer);
clSetKernelArg(clEnv->computeXCorrCoef, 4, sizeof(float),  &wInt);
clSetKernelArg(clEnv->computeXCorrCoef, 5, sizeof(int),    &nAnt);
clSetKernelArg(clEnv->computeXCorrCoef, 6, sizeof(int),    &nSamp);

workDim = 3;
size_t CijGlobalWorkSize[3] = {(size_t)nLayer, (size_t)nDir, (size_t)nBaseline};

err = clEnqueueNDRangeKernel(clEnv->queue, clEnv->computeXCorrCoef, workDim, NULL, CijGlobalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, CijBuffer, CL_TRUE, 0, sizeof(float)*nLayer*nDir*nBaseline, Cij, 0, NULL, NULL);

cout<<"Done computeXCorrCoef\n";
*/
/*
 * Sum Cij's of all baselines in each reco direction to obtain coherence M(r-hat)
 */
/*
float *M = (float*)calloc(nLayer*nDir, sizeof(float));
cl_mem MBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nLayer*nDir,
                  M, &err);
*//*
clSetKernelArg(clEnv->computeCoherence, 0, sizeof(cl_mem), &MBuffer);
clSetKernelArg(clEnv->computeCoherence, 1, sizeof(cl_mem), &CijBuffer);
clSetKernelArg(clEnv->computeCoherence, 2, sizeof(int),    &nBaseline);
*//*
clSetKernelArg(clEnv->computeNormalizedCoherence, 0, sizeof(cl_mem), &MBuffer);
clSetKernelArg(clEnv->computeNormalizedCoherence, 1, sizeof(cl_mem), &CijBuffer);
clSetKernelArg(clEnv->computeNormalizedCoherence, 2, sizeof(int),    &nBaseline);

workDim = 2;
size_t MGlobalWorkSize[2] = {(size_t)nLayer, (size_t)nDir};

//clEnqueueNDRangeKernel(clEnv->queue, clEnv->computeCoherence, workDim, NULL, MGlobalWorkSize, NULL, 0, NULL, NULL);
clEnqueueNDRangeKernel(clEnv->queue, clEnv->computeNormalizedCoherence, workDim, NULL, MGlobalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, MBuffer, CL_TRUE, 0, sizeof(float)*nLayer*nDir, M, 0, NULL, NULL);

cout<<"Done computeNormalizedCoherence\n";
*/
/*
 * Loop over M to find the max coherence and its pix index
 */
/*
float max=0.f;
int maxPixIdx;
int *rank = (int*)calloc(nLayer*nDir, sizeof(int));

TMath::Sort(nLayer*nDir, M, rank);

int *topMaxPixIdx = (int*)calloc(summary->topN, sizeof(int));
float *topMaxPixCoherence = (float*)calloc(summary->topN, sizeof(float));

for(int i=0; i<summary->topN; i++){
   topMaxPixIdx[i] = rank[i];
   topMaxPixCoherence[i] = M[rank[i]];
}

summary->setTopMaxPixInfo(topMaxPixIdx, topMaxPixCoherence);

//int *rankEachLayer = (int*)calloc(nDir, sizeof(int));
int *maxPixIdxEachLayer = (int*)calloc(nLayer, sizeof(int));
float *maxPixCoherenceEachLayer = (float*)calloc(nLayer, sizeof(float));

cl_mem maxPixIdxEachLayerBuffer       = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(int)*nLayer,
                                        maxPixIdxEachLayer, &err);
cl_mem maxPixCoherenceEachLayerBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nLayer,
                                        maxPixCoherenceEachLayer, &err);

clSetKernelArg(clEnv->getMaxPixInfoEachLayer, 0, sizeof(cl_mem), &maxPixIdxEachLayerBuffer);
clSetKernelArg(clEnv->getMaxPixInfoEachLayer, 1, sizeof(cl_mem), &maxPixCoherenceEachLayerBuffer);
clSetKernelArg(clEnv->getMaxPixInfoEachLayer, 2, sizeof(cl_mem), &MBuffer);
clSetKernelArg(clEnv->getMaxPixInfoEachLayer, 3, sizeof(int),    &nDir);

workDim = 1;
size_t maxPixInfoGlobalWorkSize = nLayer;

clEnqueueNDRangeKernel(clEnv->queue, clEnv->getMaxPixInfoEachLayer, workDim, NULL, &maxPixInfoGlobalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, maxPixIdxEachLayerBuffer, CL_TRUE, 0, sizeof(int)*nLayer, maxPixIdxEachLayer, 0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, maxPixCoherenceEachLayerBuffer, CL_TRUE, 0, sizeof(float)*nLayer, maxPixCoherenceEachLayer, 0, NULL, NULL);

summary->setMaxPixInfoEachLayer(maxPixIdxEachLayer, maxPixCoherenceEachLayer);
*/
/*
for(int idx=0; idx<nDir; idx++){
   if(M[idx] > max){
      max = M[idx];
      maxPixIdx = idx;
   }
}
*//*
maxPixIdx = rank[0];
cout<<"max: "<<M[rank[0]]<<endl;
summary->setMaxPixInfo(rank[0], M[rank[0]]);
*/
/*
 * Transferring map data
 */
 /*
cout<<"Transferring maxPixIdx map data...\n";
//mapData[0] = &M[rank[0] / nDir];
for(int dir=0; dir<nDir; dir++) mapData[dir] = M[rank[0] / nDir + dir];
*/
/*
 * Compute likelihood and p value of skymap compared to reference map fits
 */
/*
cout<<"Computing map likelihood and p-value w.r.t. the referenc map..."<<endl;

double pValue, likelihood;
err = computeMapLikelihoodAndPValue(summary->onion->nDir, summary->onion->nLayer, REFERENCE_MAP_FIT_FILE, M, likelihood, pValue);
summary->setLikelihoodAndPValue(likelihood, pValue);
cout<<"LLH: "<<likelihood<<" P Value: "<<pValue<<endl;
*/
/*
TH1F *likelihoodDist =  new TH1F("likelihoodDist","likelihoodDist",1000, -5000,50000);
TH1F *pValueDist     =  new TH1F("pValueDist","pValueDist",1000, -1000,1000);

likelihoodDist->Fill(likelihood);
pValueDist->Fill(pValue);

likelihoodDist->Draw();
cvs.SaveAs("testLikelihoodDist.C");
pValueDist->Draw();
cvs.SaveAs("testPValueDist.C");
*/
/*
 * Write FITS file
 */
/*
cout<<"Creating Healpix map and writing to FITS....\n";
//arr<float> MArr = arr<float>(&M[0], (size_t)nDir);
arr<float> MArr = arr<float>(&M[rank[0] / nDir], (size_t)nDir);
Healpix_Map<float> skyMap = Healpix_Map<float>(MArr, RING);

fitshandle fitsOut;
//#ifdef CSW
//char filename[] = "testCSWSkyMap.fits";
//#else
//char filename[] = "testXCorrSkyMap.fits";
//#endif
remove(filename);
fitsOut.create(filename);

write_Healpix_map_to_fits(fitsOut, skyMap, PLANCK_FLOAT32);
cout<<"Healpix map written\n";
*//*
arr<float> MEachLayerArr;

for(int i=0; i<nLayer; i++){

   MEachLayerArr = arr<float>(&M[i*nDir], (size_t)nDir);
   skyMap = Healpix_Map<float>(MEachLayerArr, RING);

   char layerMap[200];
   sprintf(layerMap,"layer_%d_skymap.fits",i);
   remove(layerMap);
   fitsOut.create(layerMap);
   write_Healpix_map_to_fits(fitsOut, skyMap, PLANCK_FLOAT32);

}
*/

/*
 * Deallocate memories
 */
/*
cout<<"Deallocating memories...\n";
clReleaseMemObject(recoDelaysBuffer);
clReleaseMemObject(intensityRBuffer);
clReleaseMemObject(intensityCBuffer);
clReleaseMemObject(voltsFlatBuffer);
clReleaseMemObject(xCorrCBuffer);
clReleaseMemObject(xCorrRBuffer);
clReleaseMemObject(xCorrTimeBuffer);
clReleaseMemObject(xCorrEnvBuffer);
clReleaseMemObject(sqrtWfPwrBuffer);
clReleaseMemObject(CijBuffer);
clReleaseMemObject(MBuffer);
clReleaseMemObject(maxPixIdxEachLayerBuffer);
clReleaseMemObject(maxPixCoherenceEachLayerBuffer);
free(voltsFlat);
//free(volts);
//free(recoDelays);
free(intensity_data_r);
free(intensity_data_c);
free(xCorr_data_r);
free(xCorr_data_c);
free(xCorrTime);
free(sqrtWfPwr);
free(Cij);
free(M);
free(rank);
free(topMaxPixIdx);
free(topMaxPixCoherence);
free(maxPixIdxEachLayer);
free(maxPixCoherenceEachLayer);
cout<<"Memories deallocated\n";

return maxPixIdx;
}
*/
int reconstruct3DXCorrEnvelopeGetMaxPixAndMapData(recoSettings *settings, vector<TGraph *>& cleanEvent, recoEnvData *clEnv,
                float *recoDelays, float *recoDelays_V, float *recoDelays_H,
                const int *chanMask, recoData *summary, char *filename, float *mapData)
{

cout<<"Entered reconstruct3DXCorrEnvelopeGetMaxPixAndMapData method\n";
int nSamp;
//int nAnt = (int)cleanEvent.size()/2; // Divide by 2 for only one polarization
int nAnt;// = (int)cleanEvent.size();
string pol = string(settings->recoPolType);
if( pol == "vpol" || pol == "hpol" ) nAnt = (int)cleanEvent.size()/2; else  nAnt = (int)cleanEvent.size();
int unmaskedNChan=0;
for(int ch=0; ch<(int)cleanEvent.size(); ch++) unmaskedNChan+=chanMask[ch];
cout<<"unmaskedNChan: "<<unmaskedNChan<<" nAnt: "<<nAnt<<endl;
float wInt;
int dataType = settings->dataType;
if( dataType == 0 ) wInt = 0.5f; //AraSim event
else if( dataType == 1 ){ //real event
if( pol == "vpol" ) wInt = 0.4f;
else if ( pol == "hpol" ) wInt = 0.625f;
else if ( pol == "both" ) wInt = 0.5f;
else { cerr<<"recoPolType undefined\n"; return -1; }
} else {
cerr<<"dataType undefined\n"; return -1; }

int nSideExp = settings->nSideExp;
int nDir   = /*summary->onion*/12 * pow(2,nSideExp) * pow(2,nSideExp);;
int nLayer = /*summary->onion*/settings->nLayer;
printf("nDir: %d nLayer: %d\n",nDir,nLayer);
/*
 * Loading voltsFlat array
 */

double t, v;
float *voltsFlat;

if( pol == "vpol" ){

   /* Using the 1st vpol wf for nSamp. Should make sure all vpol channels have the same nSamp. FIXME */
   nSamp = cleanEvent[0]->GetN();
   voltsFlat = (float*)calloc(nAnt*nSamp, sizeof(float));
   if(voltsFlat == NULL){ cerr<<"Null pointer to voltsFlat\n"; return -1; }

   for(int ch=0; ch<nAnt; ch++){
       if( chanMask[ch] ){
         for(int s=0; s<nSamp; s++){
            cleanEvent[ch]->GetPoint(s,t,v);
            /* Bartlett window applied in main analysis code */
            //voltsFlat[nSamp*ch + s] = static_cast<float>(v * FFTtools::bartlettWindow(s, nSamp));
            voltsFlat[nSamp*ch + s] = static_cast<float>(v);
         }
      } else {
         for(int s=0; s<nSamp; s++){
            voltsFlat[nSamp*ch + s] = 0.f;
         }
      }
   }

} else if ( pol == "hpol" ){

   /* Using the 1st hpol wf for nSamp. Should make sure all hpol channels have the same nSamp. FIXME */
   nSamp = cleanEvent[nAnt]->GetN();
   voltsFlat = (float*)calloc(nAnt*nSamp, sizeof(float));
   if(voltsFlat == NULL){ cerr<<"Null pointer to voltsFlat\n"; return -1; }

   for(int ch=0; ch<nAnt; ch++){
       if( chanMask[ch+nAnt] ){
         for(int s=0; s<nSamp; s++){
            cleanEvent[ch+nAnt]->GetPoint(s,t,v);
            /* Bartlett window applied in main analysis code */
            //voltsFlat[nSamp*ch + s] = static_cast<float>(v * FFTtools::bartlettWindow(s, nSamp));
            voltsFlat[nSamp*ch + s] = static_cast<float>(v);
         }
      } else {
         for(int s=0; s<nSamp; s++){
            voltsFlat[nSamp*ch + s] = 0.f;
         }
      }
   }

} else if ( pol == "both" ){

   /* Using the 1st vpol wf for nSamp. Should make sure all hpol channels have the same nSamp. FIXME */
   nSamp = cleanEvent[0]->GetN();
   voltsFlat = (float*)calloc(nAnt*nSamp, sizeof(float));
   if(voltsFlat == NULL){ cerr<<"Null pointer to voltsFlat\n"; return -1; }

   for(int ch=0; ch<nAnt; ch++){
       if( chanMask[ch] ){
         for(int s=0; s<nSamp; s++){
            cleanEvent[ch]->GetPoint(s,t,v);
            /* Bartlett window applied in main analysis code */
            //voltsFlat[nSamp*ch + s] = static_cast<float>(v * FFTtools::bartlettWindow(s, nSamp));
            voltsFlat[nSamp*ch + s] = static_cast<float>(v);
         }
      } else {
         for(int s=0; s<nSamp; s++){
            voltsFlat[nSamp*ch + s] = 0.f;
         }
      }
   }

} else {
   cerr<<"recoPolType not defined\n"; return -1;
}
cout<<"voltsFlat loaded\n";

/*
 * Preparation for OUT_OF_PLACE transforms
 */

int interlvHermOutputSize = 2*(1 + nSamp/2); //Hermitian layout
int interlvOutputSize     = 2*nSamp;         //Not Hermitian layout
int planarHermOutputSize = (1 + nSamp/2);    //Hermitian planar layout

float *intensity_data_r, *intensity_data_c;
intensity_data_r = (float*)calloc(nAnt*planarHermOutputSize, sizeof(float));
intensity_data_c = (float*)calloc(nAnt*planarHermOutputSize, sizeof(float));

/*
 * FFT library related declarations
 */
//clfftPlanHandle planHandle;
clfftDim fftDim = CLFFT_1D;
size_t clLengths[1]   = {nSamp};
size_t clInStride[1]  = {1};
size_t clOutStride[1] = {1};
size_t inDist  = nSamp;
size_t outDist   = planarHermOutputSize;
size_t batchSize = nAnt;

/*
 * Set up clFFT
 */
//clfftSetupData fftSetup;
//err = clfftInitSetupData(&fftSetup);
//err = clfftSetup(&fftSetup);

/*
 * Prepare plan
 */
cout<<"Preparing plan...\n";
int err;
err = clfftCreateDefaultPlan(&clEnv->planHandle, clEnv->context, fftDim, clLengths);
err = clfftSetPlanPrecision(clEnv->planHandle, CLFFT_SINGLE);
err = clfftSetLayout(clEnv->planHandle, CLFFT_REAL, CLFFT_HERMITIAN_PLANAR);
err = clfftSetPlanScale(clEnv->planHandle, CLFFT_FORWARD, 1.f);
err = clfftSetPlanBatchSize(clEnv->planHandle, batchSize);
err = clfftSetPlanInStride(clEnv->planHandle, fftDim, clInStride);
err = clfftSetPlanOutStride(clEnv->planHandle, fftDim, clOutStride);
err = clfftSetPlanDistance(clEnv->planHandle, inDist, outDist);
err = clfftSetResultLocation(clEnv->planHandle, CLFFT_OUTOFPLACE);
err = clfftBakePlan(clEnv->planHandle, 1, &clEnv->queue, NULL, NULL);
cout<<"Plan prepared\n";
/* The plan is now ready to be executed */
cl_mem voltsFlatBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*nSamp*sizeof(float), NULL, &err);
err = clEnqueueWriteBuffer(clEnv->queue, voltsFlatBuffer, CL_TRUE, 0, nAnt*nSamp*sizeof(float), voltsFlat, 0, NULL, NULL);

cl_mem intensityRBuffer= clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*planarHermOutputSize*sizeof(float),
                         NULL, &err);
cl_mem intensityCBuffer= clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*planarHermOutputSize*sizeof(float),
                         NULL, &err);
cl_mem outBuffers[2] = {intensityRBuffer, intensityCBuffer};

cout<<"Enqueueing FFT\n";
err = clfftEnqueueTransform(clEnv->planHandle, CLFFT_FORWARD, 1, &clEnv->queue, 0, NULL, NULL,
                            &voltsFlatBuffer, outBuffers, NULL);
err = clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, outBuffers[0], CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_r,
                    0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, outBuffers[1], CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_c,
                    0, NULL, NULL);
cout<<"FFT done\n";
/*
 * Clean up CLFFT
 */
cout<<"Destroying plan...\n";
err = clfftDestroyPlan(&clEnv->planHandle);
//clfftTeardown();
cout<<"Plan destroyed\n";


#ifdef bandpass
/*
 * Bandpass signals
 */
cout<<"Preparing bandpass filter..."<<endl;
float freqBin = 1e3 / (wInt * (float)nSamp); // wInt in ns. 1e3 for MHz
float lowFreq = 200.f;
float highFreq= 450.f;

clSetKernelArg(clEnv->bandPassFilter, 0, sizeof(cl_mem), &intensityRBuffer);
clSetKernelArg(clEnv->bandPassFilter, 1, sizeof(cl_mem), &intensityCBuffer);
clSetKernelArg(clEnv->bandPassFilter, 2, sizeof(float),  &freqBin);
clSetKernelArg(clEnv->bandPassFilter, 3, sizeof(float),  &lowFreq);
clSetKernelArg(clEnv->bandPassFilter, 4, sizeof(float),  &highFreq);

unsigned int dim = 2;
size_t bandPassWorkSize[2] = {(size_t)nAnt, (size_t)planarHermOutputSize};
clEnqueueNDRangeKernel(clEnv->queue, clEnv->bandPassFilter, dim, NULL, bandPassWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, intensityRBuffer, CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_r,
                    0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, intensityCBuffer, CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_c,
                    0, NULL, NULL);

cout<<"Bandpass filter done"<<endl;
#endif

/*
 * Cross-correlate wfs
 */

const int nBaseline = nAnt*nAnt;
cout<<"nBaseline: "<<nBaseline<<endl;

float *xCorr_data_r = (float*)calloc(nBaseline*planarHermOutputSize, sizeof(float));
float *xCorr_data_c = (float*)calloc(nBaseline*planarHermOutputSize, sizeof(float));

cl_mem xCorrRBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                    sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_r, &err);

cl_mem xCorrCBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                    sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_c, &err);

clSetKernelArg(clEnv->xCorrWf, 0, sizeof(cl_mem), &xCorrRBuffer);
clSetKernelArg(clEnv->xCorrWf, 1, sizeof(cl_mem), &xCorrCBuffer);
clSetKernelArg(clEnv->xCorrWf, 2, sizeof(cl_mem), &intensityRBuffer);
clSetKernelArg(clEnv->xCorrWf, 3, sizeof(cl_mem), &intensityCBuffer);
clSetKernelArg(clEnv->xCorrWf, 4, sizeof(int),    &nAnt);

unsigned int workDim = 2;
size_t globalWorkSize[2] = {(size_t)nBaseline, (size_t)planarHermOutputSize};
clEnqueueNDRangeKernel(clEnv->queue, clEnv->xCorrWf, workDim, NULL, globalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, xCorrRBuffer, CL_TRUE, 0, sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_r, 0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, xCorrCBuffer, CL_TRUE, 0, sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_c, 0, NULL, NULL);
cout<<"Done xCorrWf\n";

/*
 * Inverse FFT cross-correlation
 */

inDist  = planarHermOutputSize;
outDist = nSamp;
batchSize = nBaseline;

/*
 * Prepare plan
 */

err = clfftCreateDefaultPlan(&clEnv->planHandle, clEnv->context, fftDim, clLengths);
err = clfftSetPlanPrecision(clEnv->planHandle, CLFFT_SINGLE);
err = clfftSetLayout(clEnv->planHandle, CLFFT_HERMITIAN_PLANAR, CLFFT_REAL);
err = clfftSetPlanScale(clEnv->planHandle, CLFFT_BACKWARD, 1.f/(float)nSamp);
err = clfftSetPlanBatchSize(clEnv->planHandle, batchSize);
err = clfftSetPlanInStride(clEnv->planHandle, fftDim, clInStride);
err = clfftSetPlanOutStride(clEnv->planHandle, fftDim, clOutStride);
err = clfftSetPlanDistance(clEnv->planHandle, inDist, outDist);
err = clfftSetResultLocation(clEnv->planHandle, CLFFT_OUTOFPLACE);
err = clfftBakePlan(clEnv->planHandle, 1, &clEnv->queue, NULL, NULL);

/* The plan is now ready to be executed */

cl_mem inBuffers[2] = {xCorrRBuffer, xCorrCBuffer};

float *xCorrTime = (float*)calloc(nBaseline*nSamp, sizeof(float));
cl_mem xCorrTimeBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nBaseline*nSamp*sizeof(float),
                           NULL, &err);
err = clfftEnqueueTransform(clEnv->planHandle, CLFFT_BACKWARD, 1, &clEnv->queue, 0, NULL, NULL, inBuffers, &xCorrTimeBuffer, NULL);
err = clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, xCorrTimeBuffer, CL_TRUE, 0, sizeof(float)*nBaseline*nSamp, xCorrTime, 0, NULL, NULL);

/*
 * Clean up FFT
 */

err = clfftDestroyPlan(&clEnv->planHandle);
clfftTeardown();
cout<<"Done inverse FFT\n";

/*
 * Get Hilbert transform of XCorr function
 */

TCanvas cvs("cvs","cvs",800,600);
//cvs.Divide(1,2);
char envelopename[200];

float dt[nSamp];
float xCorrValue[nSamp];
double t_temp, v_temp;
double t_temp2, v_temp2;
int ant1, ant2;
float plusMinusTime = 25.f; //record +-25 ns around peak in XCorr function
//int plusMinusRange = (int)(plusMinusTime / wInt);
//static TGraph *sillygr[64];// = new TGraph();
//static TGraph *sillygr = (TGraph*)malloc(nBaseline*sizeof(TGraph));
/*
for(int i=0; i<nBaseline; i++){

   sillygr[i]=new TGraph();

}
*/
//static int lock;
//char sillygrname[200];

for(int baseline=0; baseline<nBaseline; baseline++){


   ant1 = baseline / nAnt;
   ant2 = baseline % nAnt;

   for(int s=0; s<nSamp; s++){

   dt[s] = wInt*s;
   xCorrValue[s] = xCorrTime[nSamp*baseline + s];

   }

   TGraph *xCorrGraph = new TGraph(nSamp, dt, xCorrValue);

   //cvs.cd(1);
   //xCorrGraph->Draw("AL");

   TGraph* envelope = FFTtools::getHilbertEnvelope( xCorrGraph );
/*
   sprintf(envelopename,"xCorrEnvelope_2014_A3_burn_RF_chan%d_%d.C", ant1, ant2);
   cvs.cd();
   //envelope->Draw("AL");
   envelope->SetLineColor(kRed);
   xCorrGraph->Draw("AL");
   envelope->Draw("Lsame");
   cvs.SaveAs(envelopename);
*/
/*
   if( ant1 == 0 && ant2 == 3){
   cvs.cd(1);
   envelope->Draw("AL");
   //xCorrGraph->Draw("AL");
   cvs.SaveAs("xCorrEnvelope_chan0_3.C");
   }
*/
/*
   if( ant1 == 3 && ant2 == 0){
   cvs.cd(2);
   envelope->Draw("AL");
   }

   cvs.SaveAs("xCorrEnvelope_chan0_3.C");
*/
/*
   if( ant1 == 0 && chanMask[ant1] ){
      cout<<"recording XCorr values around peak of each baseline\n";

      if( ant2 == 1 && chanMask[ant2] ){
      stackXCorrAroundPeak(envelope, xCorrAroundPeakHist[0], plusMinusTime);
      } else if
      ( ant2 == 3 && chanMask[ant2] ){
      stackXCorrAroundPeak(envelope, xCorrAroundPeakHist[1], plusMinusTime);
      } else if
      ( ant2 == 7 && chanMask[ant2] ){
      stackXCorrAroundPeak(envelope, xCorrAroundPeakHist[2], plusMinusTime);
      }
   }
*/
   //if(ant1==1 && ant2==5){
   //if(ant1 < ant2 ){
   //for(int s=0; s<nSamp; s++){

   //if(ant1 == 0  && ant2 == 4 ){
   //xCorrGraph->GetPoint(s,t_temp,v_temp);
   //if(lock == 0){
   //if(envelopeSum->GetN() != 0 ) cerr<<"Both TGraphs should be uninitialized!"<<endl;
/*
   if(sillygr[baseline]->GetN() == 0){

      for(int s=0; s<nSamp; s++){
      xCorrGraph->GetPoint(s,t_temp,v_temp);
      sillygr[baseline]->SetPoint(s, t_temp, v_temp);
      }

   } else {

      for(int s=0; s<nSamp; s++){
      xCorrGraph->GetPoint(s,t_temp,v_temp);
      sillygr[baseline]->GetPoint(s, t_temp2, v_temp2);
      if(t_temp != t_temp2) cerr<<"t_temp != t_temp2 !!\n";
      sillygr[baseline]->SetPoint(s, t_temp2, v_temp+v_temp2);
      }
   }
*/
   //}

   for(int s=0; s<nSamp; s++){

   envelope->GetPoint(s,t_temp,v_temp);
   xCorrTime[nSamp*baseline + s] = static_cast<float>(v_temp);

   }
   //if(lock == 0) lock=1;
/*
  cvs.cd();
  sprintf(sillygrname,"xCorrSumGraph_chan%d_%d.C",ant1,ant2);
  sillygr[baseline]->Draw("AL");
  cvs.SaveAs(sillygrname);
*/

  delete xCorrGraph;
  delete envelope;

}

cl_mem xCorrEnvBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                       sizeof(float)*nBaseline*nSamp,
                                       xCorrTime, &err);


/*
 * Compute cross-correlation coefficients Cij in each direction
 */

/* Compute the square root of each channel's total wf power */

float *sqrtWfPwr = (float*)calloc(nAnt, sizeof(float));
float pwr=0.f;
for(int ant=0; ant<nAnt; ant++){
   pwr = 0.;
   for(int s=0; s<nSamp; s++){
      pwr += (voltsFlat[ant*nSamp + s] * voltsFlat[ant*nSamp + s]);
   }
   sqrtWfPwr[ant] = sqrt(pwr);
   //sqrtWfPwr[ant] = 1.f;
}
cout<<"Done sqrtWfPwr\n";

cl_mem sqrtWfPwrBuffer  = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt,
                          sqrtWfPwr, &err);
cl_mem recoDelaysBuffer;
if(pol == "vpol" ){
recoDelaysBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt*nDir*nLayer,
                                  recoDelays_V, &err);
} else if (pol == "hpol" ) {
recoDelaysBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt*nDir*nLayer,
                                  recoDelays_H, &err);
} else if (pol == "both" ) {
recoDelaysBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt*nDir*nLayer,
                                  recoDelays, &err);
} else {
cerr<<"recoPolType not defined\n"; return -1;
}
cout<<"recoDelaysBuffer created\n";

float *Cij = (float*)calloc(nLayer*nDir*nBaseline, sizeof(float));

cl_mem CijBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nLayer*nDir*nBaseline,
                   Cij, &err);

clSetKernelArg(clEnv->computeXCorrCoef, 0, sizeof(cl_mem), &CijBuffer);
clSetKernelArg(clEnv->computeXCorrCoef, 1, sizeof(cl_mem), &xCorrEnvBuffer);
clSetKernelArg(clEnv->computeXCorrCoef, 2, sizeof(cl_mem), &recoDelaysBuffer);
clSetKernelArg(clEnv->computeXCorrCoef, 3, sizeof(cl_mem), &sqrtWfPwrBuffer);
clSetKernelArg(clEnv->computeXCorrCoef, 4, sizeof(float),  &wInt);
clSetKernelArg(clEnv->computeXCorrCoef, 5, sizeof(int),    &nAnt);
clSetKernelArg(clEnv->computeXCorrCoef, 6, sizeof(int),    &nSamp);

workDim = 3;
size_t CijGlobalWorkSize[3] = {(size_t)nLayer, (size_t)nDir, (size_t)nBaseline};

err = clEnqueueNDRangeKernel(clEnv->queue, clEnv->computeXCorrCoef, workDim, NULL, CijGlobalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, CijBuffer, CL_TRUE, 0, sizeof(float)*nLayer*nDir*nBaseline, Cij, 0, NULL, NULL);

cout<<"Done computeXCorrCoef\n";

/*
 * Sum Cij's of all baselines in each reco direction to obtain coherence M(r-hat)
 */

float *M = (float*)calloc(nLayer*nDir, sizeof(float));
cl_mem MBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nLayer*nDir,
                  M, &err);
/*
clSetKernelArg(clEnv->computeCoherence, 0, sizeof(cl_mem), &MBuffer);
clSetKernelArg(clEnv->computeCoherence, 1, sizeof(cl_mem), &CijBuffer);
clSetKernelArg(clEnv->computeCoherence, 2, sizeof(int),    &nBaseline);
*/
clSetKernelArg(clEnv->computeNormalizedCoherence, 0, sizeof(cl_mem), &MBuffer);
clSetKernelArg(clEnv->computeNormalizedCoherence, 1, sizeof(cl_mem), &CijBuffer);
clSetKernelArg(clEnv->computeNormalizedCoherence, 2, sizeof(int),    &nBaseline);

workDim = 2;
size_t MGlobalWorkSize[2] = {(size_t)nLayer, (size_t)nDir};

//clEnqueueNDRangeKernel(clEnv->queue, clEnv->computeCoherence, workDim, NULL, MGlobalWorkSize, NULL, 0, NULL, NULL);
clEnqueueNDRangeKernel(clEnv->queue, clEnv->computeNormalizedCoherence, workDim, NULL, MGlobalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, MBuffer, CL_TRUE, 0, sizeof(float)*nLayer*nDir, M, 0, NULL, NULL);

cout<<"Done computeNormalizedCoherence\n";

/*
 * Loop over M to find the max coherence and its pix index
 */

float max=0.f;
int maxPixIdx;
int *rank = (int*)calloc(nLayer*nDir, sizeof(int));

TMath::Sort(nLayer*nDir, M, rank);

int *topMaxPixIdx = (int*)calloc(summary->topN, sizeof(int));
float *topMaxPixCoherence = (float*)calloc(summary->topN, sizeof(float));

for(int i=0; i<summary->topN; i++){
   topMaxPixIdx[i] = rank[i];
   topMaxPixCoherence[i] = M[rank[i]];
}

summary->setTopMaxPixInfo(topMaxPixIdx, topMaxPixCoherence);

//int *rankEachLayer = (int*)calloc(nDir, sizeof(int));
int *maxPixIdxEachLayer = (int*)calloc(nLayer, sizeof(int));
float *maxPixCoherenceEachLayer = (float*)calloc(nLayer, sizeof(float));

cl_mem maxPixIdxEachLayerBuffer       = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(int)*nLayer,
                                        maxPixIdxEachLayer, &err);
cl_mem maxPixCoherenceEachLayerBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nLayer,
                                        maxPixCoherenceEachLayer, &err);

clSetKernelArg(clEnv->getMaxPixInfoEachLayer, 0, sizeof(cl_mem), &maxPixIdxEachLayerBuffer);
clSetKernelArg(clEnv->getMaxPixInfoEachLayer, 1, sizeof(cl_mem), &maxPixCoherenceEachLayerBuffer);
clSetKernelArg(clEnv->getMaxPixInfoEachLayer, 2, sizeof(cl_mem), &MBuffer);
clSetKernelArg(clEnv->getMaxPixInfoEachLayer, 3, sizeof(int),    &nDir);

workDim = 1;
size_t maxPixInfoGlobalWorkSize = nLayer;

clEnqueueNDRangeKernel(clEnv->queue, clEnv->getMaxPixInfoEachLayer, workDim, NULL, &maxPixInfoGlobalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, maxPixIdxEachLayerBuffer, CL_TRUE, 0, sizeof(int)*nLayer, maxPixIdxEachLayer, 0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, maxPixCoherenceEachLayerBuffer, CL_TRUE, 0, sizeof(float)*nLayer, maxPixCoherenceEachLayer, 0, NULL, NULL);

summary->setMaxPixInfoEachLayer(settings, maxPixIdxEachLayer, maxPixCoherenceEachLayer);

/*
for(int idx=0; idx<nDir; idx++){
   if(M[idx] > max){
      max = M[idx];
      maxPixIdx = idx;
   }
}
*/
maxPixIdx = rank[0];
cout<<"max: "<<M[rank[0]]<<endl;
summary->setMaxPixInfo(rank[0], M[rank[0]]);

/*
 * Transferring map data
 */

cout<<"Transferring maxPixIdx map data...\n";
//mapData[0] = &M[rank[0] / nDir];
for(int dir=0; dir<nDir; dir++) mapData[dir] = M[rank[0] / nDir + dir];

/*
 * Compute likelihood and p value of skymap compared to reference map fits
 */

if(settings->computeLLHAndPValue == 1){
cout<<"Computing map likelihood and p-value w.r.t. the referenc map..."<<endl;

double pValue, likelihood;
err = computeMapLikelihoodAndPValue(/*summary->onion*/nDir, /*summary->onion*/nLayer, settings->referenceMapFitFunc/*.c_str()*/, settings->referenceMapFitFile/*.c_str()*/, M, likelihood, pValue);
summary->setLikelihoodAndPValue(likelihood, pValue);
cout<<"LLH: "<<likelihood<<" P Value: "<<pValue<<endl;
}
/*
TH1F *likelihoodDist =  new TH1F("likelihoodDist","likelihoodDist",1000, -5000,50000);
TH1F *pValueDist     =  new TH1F("pValueDist","pValueDist",1000, -1000,1000);

likelihoodDist->Fill(likelihood);
pValueDist->Fill(pValue);

likelihoodDist->Draw();
cvs.SaveAs("testLikelihoodDist.C");
pValueDist->Draw();
cvs.SaveAs("testPValueDist.C");
*/
/*
 * Write FITS file
 */

cout<<"Creating Healpix map and writing to FITS....\n";
//arr<float> MArr = arr<float>(&M[0], (size_t)nDir);
arr<float> MArr = arr<float>(&M[rank[0] / nDir], (size_t)nDir);
Healpix_Map<float> skyMap = Healpix_Map<float>(MArr, HEALPIX_ORDERING);

fitshandle fitsOut;
//#ifdef CSW
//char filename[] = "testCSWSkyMap.fits";
//#else
//char filename[] = "testXCorrSkyMap.fits";
//#endif
remove(filename);
fitsOut.create(filename);

write_Healpix_map_to_fits(fitsOut, skyMap, PLANCK_FLOAT32);
cout<<"Healpix map written\n";
/*
arr<float> MEachLayerArr;

for(int i=0; i<nLayer; i++){

   MEachLayerArr = arr<float>(&M[i*nDir], (size_t)nDir);
   skyMap = Healpix_Map<float>(MEachLayerArr, RING);

   char layerMap[200];
   sprintf(layerMap,"layer_%d_skymap.fits",i);
   remove(layerMap);
   fitsOut.create(layerMap);
   write_Healpix_map_to_fits(fitsOut, skyMap, PLANCK_FLOAT32);

}
*/

/*
 * Deallocate memories
 */

cout<<"Deallocating memories...\n";
clReleaseMemObject(recoDelaysBuffer);
clReleaseMemObject(intensityRBuffer);
clReleaseMemObject(intensityCBuffer);
clReleaseMemObject(voltsFlatBuffer);
clReleaseMemObject(xCorrCBuffer);
clReleaseMemObject(xCorrRBuffer);
clReleaseMemObject(xCorrTimeBuffer);
clReleaseMemObject(xCorrEnvBuffer);
clReleaseMemObject(sqrtWfPwrBuffer);
clReleaseMemObject(CijBuffer);
clReleaseMemObject(MBuffer);
clReleaseMemObject(maxPixIdxEachLayerBuffer);
clReleaseMemObject(maxPixCoherenceEachLayerBuffer);
free(voltsFlat);
//free(volts);
//free(recoDelays);
free(intensity_data_r);
free(intensity_data_c);
free(xCorr_data_r);
free(xCorr_data_c);
free(xCorrTime);
free(sqrtWfPwr);
free(Cij);
free(M);
free(rank);
free(topMaxPixIdx);
free(topMaxPixCoherence);
free(maxPixIdxEachLayer);
free(maxPixCoherenceEachLayer);
cout<<"Memories deallocated\n";

return maxPixIdx;
}

int reconstruct3DXCorrEnvelopeGetMaxPixAndMapData_overlapCorrection(recoSettings *settings, vector<TGraph *>& cleanEvent, recoEnvData *clEnv,
                float *recoDelays, float *recoDelays_V, float *recoDelays_H, const double *beginTimeByChannel, const int *wfNBins,
                const int *chanMask, recoData *summary, char *filename, float *mapData)
{

cout<<"Entered reconstruct3DXCorrEnvelopeGetMaxPixAndMapData_overlapCorrection method\n";
int nSamp;
//int nAnt = (int)cleanEvent.size()/2; // Divide by 2 for only one polarization
int nAnt;// = (int)cleanEvent.size();
string pol = string(settings->recoPolType);
if( pol == "vpol" || pol == "hpol" ) nAnt = (int)cleanEvent.size()/2; else  nAnt = (int)cleanEvent.size();
int unmaskedNChan=0;
for(int ch=0; ch<(int)cleanEvent.size(); ch++) unmaskedNChan+=chanMask[ch];
cout<<"unmaskedNChan: "<<unmaskedNChan<<" nAnt: "<<nAnt<<endl;
float wInt;
int dataType = settings->dataType;
if( dataType == 0 ) wInt = 0.5f; //AraSim event
else if( dataType == 1 ){ //real event
if( pol == "vpol" ) wInt = 0.4f;
else if ( pol == "hpol" ) wInt = 0.625f;
else if ( pol == "both" ) wInt = 0.5f;
else { cerr<<"recoPolType undefined\n"; return -1; }
} else {
cerr<<"dataType undefined\n"; return -1; }

int nSideExp = settings->nSideExp;
int nDir   = /*summary->onion*/12 * pow(2,nSideExp) * pow(2,nSideExp);;
int nLayer = /*summary->onion*/settings->nLayer;
printf("nDir: %d nLayer: %d\n",nDir,nLayer);
/*
 * Loading voltsFlat array
 */

double t, v;
float *voltsFlat;

if( pol == "vpol" ){

   /* Using the 1st vpol wf for nSamp. Should make sure all vpol channels have the same nSamp. FIXME */
   nSamp = cleanEvent[0]->GetN();
   voltsFlat = (float*)calloc(nAnt*nSamp, sizeof(float));
   if(voltsFlat == NULL){ cerr<<"Null pointer to voltsFlat\n"; return -1; }

   for(int ch=0; ch<nAnt; ch++){
       if( chanMask[ch] ){
         for(int s=0; s<nSamp; s++){
            cleanEvent[ch]->GetPoint(s,t,v);
            /* Bartlett window applied in main analysis code */
            //voltsFlat[nSamp*ch + s] = static_cast<float>(v * FFTtools::bartlettWindow(s, nSamp));
            voltsFlat[nSamp*ch + s] = static_cast<float>(v);
         }
      } else {
         for(int s=0; s<nSamp; s++){
            voltsFlat[nSamp*ch + s] = 0.f;
         }
      }
   }

} else if ( pol == "hpol" ){

   /* Using the 1st hpol wf for nSamp. Should make sure all hpol channels have the same nSamp. FIXME */
   nSamp = cleanEvent[nAnt]->GetN();
   voltsFlat = (float*)calloc(nAnt*nSamp, sizeof(float));
   if(voltsFlat == NULL){ cerr<<"Null pointer to voltsFlat\n"; return -1; }

   for(int ch=0; ch<nAnt; ch++){
       if( chanMask[ch+nAnt] ){
         for(int s=0; s<nSamp; s++){
            cleanEvent[ch+nAnt]->GetPoint(s,t,v);
            /* Bartlett window applied in main analysis code */
            //voltsFlat[nSamp*ch + s] = static_cast<float>(v * FFTtools::bartlettWindow(s, nSamp));
            voltsFlat[nSamp*ch + s] = static_cast<float>(v);
         }
      } else {
         for(int s=0; s<nSamp; s++){
            voltsFlat[nSamp*ch + s] = 0.f;
         }
      }
   }

} else if ( pol == "both" ){

   /* Using the 1st vpol wf for nSamp. Should make sure all hpol channels have the same nSamp. FIXME */
   nSamp = cleanEvent[0]->GetN();
   voltsFlat = (float*)calloc(nAnt*nSamp, sizeof(float));
   if(voltsFlat == NULL){ cerr<<"Null pointer to voltsFlat\n"; return -1; }

   for(int ch=0; ch<nAnt; ch++){
       if( chanMask[ch] ){
         for(int s=0; s<nSamp; s++){
            cleanEvent[ch]->GetPoint(s,t,v);
            /* Bartlett window applied in main analysis code */
            //voltsFlat[nSamp*ch + s] = static_cast<float>(v * FFTtools::bartlettWindow(s, nSamp));
            voltsFlat[nSamp*ch + s] = static_cast<float>(v);
         }
      } else {
         for(int s=0; s<nSamp; s++){
            voltsFlat[nSamp*ch + s] = 0.f;
         }
      }
   }

} else {
   cerr<<"recoPolType not defined\n"; return -1;
}
cout<<"voltsFlat loaded\n";

/*
 * Preparation for OUT_OF_PLACE transforms
 */

int interlvHermOutputSize = 2*(1 + nSamp/2); //Hermitian layout
int interlvOutputSize     = 2*nSamp;         //Not Hermitian layout
int planarHermOutputSize = (1 + nSamp/2);    //Hermitian planar layout

float *intensity_data_r, *intensity_data_c;
intensity_data_r = (float*)calloc(nAnt*planarHermOutputSize, sizeof(float));
intensity_data_c = (float*)calloc(nAnt*planarHermOutputSize, sizeof(float));

/*
 * FFT library related declarations
 */
//clfftPlanHandle planHandle;
clfftDim fftDim = CLFFT_1D;
size_t clLengths[1]   = {nSamp};
size_t clInStride[1]  = {1};
size_t clOutStride[1] = {1};
size_t inDist  = nSamp;
size_t outDist   = planarHermOutputSize;
size_t batchSize = nAnt;

/*
 * Set up clFFT
 */
//clfftSetupData fftSetup;
//err = clfftInitSetupData(&fftSetup);
//err = clfftSetup(&fftSetup);

/*
 * Prepare plan
 */
cout<<"Preparing plan...\n";
int err;
err = clfftCreateDefaultPlan(&clEnv->planHandle, clEnv->context, fftDim, clLengths);
err = clfftSetPlanPrecision(clEnv->planHandle, CLFFT_SINGLE);
err = clfftSetLayout(clEnv->planHandle, CLFFT_REAL, CLFFT_HERMITIAN_PLANAR);
err = clfftSetPlanScale(clEnv->planHandle, CLFFT_FORWARD, 1.f);
err = clfftSetPlanBatchSize(clEnv->planHandle, batchSize);
err = clfftSetPlanInStride(clEnv->planHandle, fftDim, clInStride);
err = clfftSetPlanOutStride(clEnv->planHandle, fftDim, clOutStride);
err = clfftSetPlanDistance(clEnv->planHandle, inDist, outDist);
err = clfftSetResultLocation(clEnv->planHandle, CLFFT_OUTOFPLACE);
err = clfftBakePlan(clEnv->planHandle, 1, &clEnv->queue, NULL, NULL);
cout<<"Plan prepared\n";
/* The plan is now ready to be executed */
cl_mem voltsFlatBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*nSamp*sizeof(float), NULL, &err);
err = clEnqueueWriteBuffer(clEnv->queue, voltsFlatBuffer, CL_TRUE, 0, nAnt*nSamp*sizeof(float), voltsFlat, 0, NULL, NULL);

cl_mem intensityRBuffer= clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*planarHermOutputSize*sizeof(float),
                         NULL, &err);
cl_mem intensityCBuffer= clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*planarHermOutputSize*sizeof(float),
                         NULL, &err);
cl_mem outBuffers[2] = {intensityRBuffer, intensityCBuffer};

cout<<"Enqueueing FFT\n";
err = clfftEnqueueTransform(clEnv->planHandle, CLFFT_FORWARD, 1, &clEnv->queue, 0, NULL, NULL,
                            &voltsFlatBuffer, outBuffers, NULL);
err = clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, outBuffers[0], CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_r,
                    0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, outBuffers[1], CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_c,
                    0, NULL, NULL);
cout<<"FFT done\n";
/*
 * Clean up CLFFT
 */
cout<<"Destroying plan...\n";
err = clfftDestroyPlan(&clEnv->planHandle);
//clfftTeardown();
cout<<"Plan destroyed\n";


#ifdef bandpass
/*
 * Bandpass signals
 */
cout<<"Preparing bandpass filter..."<<endl;
float freqBin = 1e3 / (wInt * (float)nSamp); // wInt in ns. 1e3 for MHz
float lowFreq = 200.f;
float highFreq= 450.f;

clSetKernelArg(clEnv->bandPassFilter, 0, sizeof(cl_mem), &intensityRBuffer);
clSetKernelArg(clEnv->bandPassFilter, 1, sizeof(cl_mem), &intensityCBuffer);
clSetKernelArg(clEnv->bandPassFilter, 2, sizeof(float),  &freqBin);
clSetKernelArg(clEnv->bandPassFilter, 3, sizeof(float),  &lowFreq);
clSetKernelArg(clEnv->bandPassFilter, 4, sizeof(float),  &highFreq);

unsigned int dim = 2;
size_t bandPassWorkSize[2] = {(size_t)nAnt, (size_t)planarHermOutputSize};
clEnqueueNDRangeKernel(clEnv->queue, clEnv->bandPassFilter, dim, NULL, bandPassWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, intensityRBuffer, CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_r,
                    0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, intensityCBuffer, CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_c,
                    0, NULL, NULL);

cout<<"Bandpass filter done"<<endl;
#endif

/*
 * Cross-correlate wfs
 */

const int nBaseline = nAnt*nAnt;
cout<<"nBaseline: "<<nBaseline<<endl;

float *xCorr_data_r = (float*)calloc(nBaseline*planarHermOutputSize, sizeof(float));
float *xCorr_data_c = (float*)calloc(nBaseline*planarHermOutputSize, sizeof(float));

cl_mem xCorrRBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                    sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_r, &err);

cl_mem xCorrCBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                    sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_c, &err);

clSetKernelArg(clEnv->xCorrWf, 0, sizeof(cl_mem), &xCorrRBuffer);
clSetKernelArg(clEnv->xCorrWf, 1, sizeof(cl_mem), &xCorrCBuffer);
clSetKernelArg(clEnv->xCorrWf, 2, sizeof(cl_mem), &intensityRBuffer);
clSetKernelArg(clEnv->xCorrWf, 3, sizeof(cl_mem), &intensityCBuffer);
clSetKernelArg(clEnv->xCorrWf, 4, sizeof(int),    &nAnt);

unsigned int workDim = 2;
size_t globalWorkSize[2] = {(size_t)nBaseline, (size_t)planarHermOutputSize};
clEnqueueNDRangeKernel(clEnv->queue, clEnv->xCorrWf, workDim, NULL, globalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, xCorrRBuffer, CL_TRUE, 0, sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_r, 0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, xCorrCBuffer, CL_TRUE, 0, sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_c, 0, NULL, NULL);
cout<<"Done xCorrWf\n";

/*
 * Inverse FFT cross-correlation
 */

inDist  = planarHermOutputSize;
outDist = nSamp;
batchSize = nBaseline;

/*
 * Prepare plan
 */

err = clfftCreateDefaultPlan(&clEnv->planHandle, clEnv->context, fftDim, clLengths);
err = clfftSetPlanPrecision(clEnv->planHandle, CLFFT_SINGLE);
err = clfftSetLayout(clEnv->planHandle, CLFFT_HERMITIAN_PLANAR, CLFFT_REAL);
err = clfftSetPlanScale(clEnv->planHandle, CLFFT_BACKWARD, 1.f/(float)nSamp);
err = clfftSetPlanBatchSize(clEnv->planHandle, batchSize);
err = clfftSetPlanInStride(clEnv->planHandle, fftDim, clInStride);
err = clfftSetPlanOutStride(clEnv->planHandle, fftDim, clOutStride);
err = clfftSetPlanDistance(clEnv->planHandle, inDist, outDist);
err = clfftSetResultLocation(clEnv->planHandle, CLFFT_OUTOFPLACE);
err = clfftBakePlan(clEnv->planHandle, 1, &clEnv->queue, NULL, NULL);

/* The plan is now ready to be executed */

cl_mem inBuffers[2] = {xCorrRBuffer, xCorrCBuffer};

float *xCorrTime = (float*)calloc(nBaseline*nSamp, sizeof(float));
cl_mem xCorrTimeBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nBaseline*nSamp*sizeof(float),
                           NULL, &err);
err = clfftEnqueueTransform(clEnv->planHandle, CLFFT_BACKWARD, 1, &clEnv->queue, 0, NULL, NULL, inBuffers, &xCorrTimeBuffer, NULL);
err = clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, xCorrTimeBuffer, CL_TRUE, 0, sizeof(float)*nBaseline*nSamp, xCorrTime, 0, NULL, NULL);

/*
 * Clean up FFT
 */

err = clfftDestroyPlan(&clEnv->planHandle);
clfftTeardown();
cout<<"Done inverse FFT\n";

/*
 * Get Hilbert transform of XCorr function
 */

TCanvas cvs("cvs","cvs",800,600);
//cvs.Divide(1,2);
char envelopename[200];

float dt[nSamp];
float xCorrValue[nSamp];
double t_temp, v_temp;
double t_temp2, v_temp2;
int ant1, ant2;
float plusMinusTime = 25.f; //record +-25 ns around peak in XCorr function
//int plusMinusRange = (int)(plusMinusTime / wInt);
//static TGraph *sillygr[64];// = new TGraph();
//static TGraph *sillygr = (TGraph*)malloc(nBaseline*sizeof(TGraph));
/*
for(int i=0; i<nBaseline; i++){

   sillygr[i]=new TGraph();

}
*/
//static int lock;
//char sillygrname[200];

for(int baseline=0; baseline<nBaseline; baseline++){


   ant1 = baseline / nAnt;
   ant2 = baseline % nAnt;

   for(int s=0; s<nSamp; s++){

   dt[s] = wInt*s;
   xCorrValue[s] = xCorrTime[nSamp*baseline + s];

   }

   TGraph *xCorrGraph = new TGraph(nSamp, dt, xCorrValue);

   //cvs.cd(1);
   //xCorrGraph->Draw("AL");

   TGraph* envelope = FFTtools::getHilbertEnvelope( xCorrGraph );
/*
   sprintf(envelopename,"xCorrEnvelope_2014_A3_burn_RF_chan%d_%d.C", ant1, ant2);
   cvs.cd();
   //envelope->Draw("AL");
   envelope->SetLineColor(kRed);
   xCorrGraph->Draw("AL");
   envelope->Draw("Lsame");
   cvs.SaveAs(envelopename);
*/
/*
   if( ant1 == 0 && ant2 == 3){
   cvs.cd(1);
   envelope->Draw("AL");
   //xCorrGraph->Draw("AL");
   cvs.SaveAs("xCorrEnvelope_chan0_3.C");
   }
*/
/*
   if( ant1 == 3 && ant2 == 0){
   cvs.cd(2);
   envelope->Draw("AL");
   }

   cvs.SaveAs("xCorrEnvelope_chan0_3.C");
*/
/*
   if( ant1 == 0 && chanMask[ant1] ){
      cout<<"recording XCorr values around peak of each baseline\n";

      if( ant2 == 1 && chanMask[ant2] ){
      stackXCorrAroundPeak(envelope, xCorrAroundPeakHist[0], plusMinusTime);
      } else if
      ( ant2 == 3 && chanMask[ant2] ){
      stackXCorrAroundPeak(envelope, xCorrAroundPeakHist[1], plusMinusTime);
      } else if
      ( ant2 == 7 && chanMask[ant2] ){
      stackXCorrAroundPeak(envelope, xCorrAroundPeakHist[2], plusMinusTime);
      }
   }
*/
   //if(ant1==1 && ant2==5){
   //if(ant1 < ant2 ){
   //for(int s=0; s<nSamp; s++){

   //if(ant1 == 0  && ant2 == 4 ){
   //xCorrGraph->GetPoint(s,t_temp,v_temp);
   //if(lock == 0){
   //if(envelopeSum->GetN() != 0 ) cerr<<"Both TGraphs should be uninitialized!"<<endl;
/*
   if(sillygr[baseline]->GetN() == 0){

      for(int s=0; s<nSamp; s++){
      xCorrGraph->GetPoint(s,t_temp,v_temp);
      sillygr[baseline]->SetPoint(s, t_temp, v_temp);
      }

   } else {

      for(int s=0; s<nSamp; s++){
      xCorrGraph->GetPoint(s,t_temp,v_temp);
      sillygr[baseline]->GetPoint(s, t_temp2, v_temp2);
      if(t_temp != t_temp2) cerr<<"t_temp != t_temp2 !!\n";
      sillygr[baseline]->SetPoint(s, t_temp2, v_temp+v_temp2);
      }
   }
*/
   //}

   for(int s=0; s<nSamp; s++){

   envelope->GetPoint(s,t_temp,v_temp);
   xCorrTime[nSamp*baseline + s] = static_cast<float>(v_temp);

   }
   //if(lock == 0) lock=1;
/*
  cvs.cd();
  sprintf(sillygrname,"xCorrSumGraph_chan%d_%d.C",ant1,ant2);
  sillygr[baseline]->Draw("AL");
  cvs.SaveAs(sillygrname);
*/

  delete xCorrGraph;
  delete envelope;

}

cl_mem xCorrEnvBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                       sizeof(float)*nBaseline*nSamp,
                                       xCorrTime, &err);


/*
 * Compute cross-correlation coefficients Cij in each direction
 */

/* Compute the square root of each channel's total wf power */

float *sqrtWfPwr = (float*)calloc(nAnt, sizeof(float));
float pwr=0.f;
for(int ant=0; ant<nAnt; ant++){
   pwr = 0.;
   for(int s=0; s<nSamp; s++){
      pwr += (voltsFlat[ant*nSamp + s] * voltsFlat[ant*nSamp + s]);
   }
   sqrtWfPwr[ant] = sqrt(pwr);
   //sqrtWfPwr[ant] = 1.f;
}
cout<<"Done sqrtWfPwr\n";

cl_mem sqrtWfPwrBuffer  = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt,
                          sqrtWfPwr, &err);
cl_mem recoDelaysBuffer;
if(pol == "vpol" ){
recoDelaysBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt*nDir*nLayer,
                                  recoDelays_V, &err);
} else if (pol == "hpol" ) {
recoDelaysBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt*nDir*nLayer,
                                  recoDelays_H, &err);
} else if (pol == "both" ) {
recoDelaysBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt*nDir*nLayer,
                                  recoDelays, &err);
} else {
cerr<<"recoPolType not defined\n"; return -1;
}
cout<<"recoDelaysBuffer created\n";

float *Cij = (float*)calloc(nLayer*nDir*nBaseline, sizeof(float));

cl_mem CijBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nLayer*nDir*nBaseline,
                   Cij, &err);

double *beginTimeByChannel_temp = (double*)calloc(nAnt, sizeof(double));
int *wfNBins_temp = (int*)calloc(nAnt, sizeof(int));

for(int ch=0; ch<nAnt; ch++){
  if      (pol=="vpol"){ beginTimeByChannel_temp[ch] = beginTimeByChannel[ch]; wfNBins_temp[ch] = wfNBins[ch]; }
  else if (pol=="hpol"){ beginTimeByChannel_temp[ch] = beginTimeByChannel[nAnt+ch]; wfNBins_temp[ch] = wfNBins[nAnt+ch]; }
  else if (pol=="both"){ beginTimeByChannel_temp[ch] = beginTimeByChannel[ch]; wfNBins_temp[ch] = wfNBins[ch]; }
  else    { cerr<<"recoPolType not defined\n"; return -1; }
}

cl_mem beginTimeByChannelBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(double)*nAnt,
                                  beginTimeByChannel_temp, &err);

cl_mem wfNBinsBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(int)*nAnt,
                       wfNBins_temp, &err);

clSetKernelArg(clEnv->computeXCorrCoef_overlapCorrection, 0, sizeof(cl_mem), &CijBuffer);
clSetKernelArg(clEnv->computeXCorrCoef_overlapCorrection, 1, sizeof(cl_mem), &xCorrEnvBuffer);
clSetKernelArg(clEnv->computeXCorrCoef_overlapCorrection, 2, sizeof(cl_mem), &recoDelaysBuffer);
clSetKernelArg(clEnv->computeXCorrCoef_overlapCorrection, 3, sizeof(cl_mem), &sqrtWfPwrBuffer);
clSetKernelArg(clEnv->computeXCorrCoef_overlapCorrection, 4, sizeof(cl_mem), &beginTimeByChannelBuffer);
clSetKernelArg(clEnv->computeXCorrCoef_overlapCorrection, 5, sizeof(cl_mem), &wfNBinsBuffer);
clSetKernelArg(clEnv->computeXCorrCoef_overlapCorrection, 6, sizeof(float),  &wInt);
clSetKernelArg(clEnv->computeXCorrCoef_overlapCorrection, 7, sizeof(int),    &nAnt);
clSetKernelArg(clEnv->computeXCorrCoef_overlapCorrection, 8, sizeof(int),    &nSamp);

workDim = 3;
size_t CijGlobalWorkSize[3] = {(size_t)nLayer, (size_t)nDir, (size_t)nBaseline};

err = clEnqueueNDRangeKernel(clEnv->queue, clEnv->computeXCorrCoef_overlapCorrection, workDim, NULL, CijGlobalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, CijBuffer, CL_TRUE, 0, sizeof(float)*nLayer*nDir*nBaseline, Cij, 0, NULL, NULL);

cout<<"Done computeXCorrCoef_overlapCorrection\n";

/*
 * Sum Cij's of all baselines in each reco direction to obtain coherence M(r-hat)
 */

float *M = (float*)calloc(nLayer*nDir, sizeof(float));
cl_mem MBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nLayer*nDir,
                  M, &err);
/*
clSetKernelArg(clEnv->computeCoherence, 0, sizeof(cl_mem), &MBuffer);
clSetKernelArg(clEnv->computeCoherence, 1, sizeof(cl_mem), &CijBuffer);
clSetKernelArg(clEnv->computeCoherence, 2, sizeof(int),    &nBaseline);
*/
clSetKernelArg(clEnv->computeNormalizedCoherence, 0, sizeof(cl_mem), &MBuffer);
clSetKernelArg(clEnv->computeNormalizedCoherence, 1, sizeof(cl_mem), &CijBuffer);
clSetKernelArg(clEnv->computeNormalizedCoherence, 2, sizeof(int),    &nBaseline);

workDim = 2;
size_t MGlobalWorkSize[2] = {(size_t)nLayer, (size_t)nDir};

//clEnqueueNDRangeKernel(clEnv->queue, clEnv->computeCoherence, workDim, NULL, MGlobalWorkSize, NULL, 0, NULL, NULL);
clEnqueueNDRangeKernel(clEnv->queue, clEnv->computeNormalizedCoherence, workDim, NULL, MGlobalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, MBuffer, CL_TRUE, 0, sizeof(float)*nLayer*nDir, M, 0, NULL, NULL);

cout<<"Done computeNormalizedCoherence\n";

/*
 * Loop over M to find the max coherence and its pix index
 */

float max=0.f;
int maxPixIdx;
int *rank = (int*)calloc(nLayer*nDir, sizeof(int));

TMath::Sort(nLayer*nDir, M, rank);

int *topMaxPixIdx = (int*)calloc(summary->topN, sizeof(int));
float *topMaxPixCoherence = (float*)calloc(summary->topN, sizeof(float));

for(int i=0; i<summary->topN; i++){
   topMaxPixIdx[i] = rank[i];
   topMaxPixCoherence[i] = M[rank[i]];
}

summary->setTopMaxPixInfo(topMaxPixIdx, topMaxPixCoherence);

//int *rankEachLayer = (int*)calloc(nDir, sizeof(int));
int *maxPixIdxEachLayer = (int*)calloc(nLayer, sizeof(int));
float *maxPixCoherenceEachLayer = (float*)calloc(nLayer, sizeof(float));

cl_mem maxPixIdxEachLayerBuffer       = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(int)*nLayer,
                                        maxPixIdxEachLayer, &err);
cl_mem maxPixCoherenceEachLayerBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nLayer,
                                        maxPixCoherenceEachLayer, &err);

clSetKernelArg(clEnv->getMaxPixInfoEachLayer, 0, sizeof(cl_mem), &maxPixIdxEachLayerBuffer);
clSetKernelArg(clEnv->getMaxPixInfoEachLayer, 1, sizeof(cl_mem), &maxPixCoherenceEachLayerBuffer);
clSetKernelArg(clEnv->getMaxPixInfoEachLayer, 2, sizeof(cl_mem), &MBuffer);
clSetKernelArg(clEnv->getMaxPixInfoEachLayer, 3, sizeof(int),    &nDir);

workDim = 1;
size_t maxPixInfoGlobalWorkSize = nLayer;

clEnqueueNDRangeKernel(clEnv->queue, clEnv->getMaxPixInfoEachLayer, workDim, NULL, &maxPixInfoGlobalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, maxPixIdxEachLayerBuffer, CL_TRUE, 0, sizeof(int)*nLayer, maxPixIdxEachLayer, 0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, maxPixCoherenceEachLayerBuffer, CL_TRUE, 0, sizeof(float)*nLayer, maxPixCoherenceEachLayer, 0, NULL, NULL);

summary->setMaxPixInfoEachLayer(settings, maxPixIdxEachLayer, maxPixCoherenceEachLayer);

/*
for(int idx=0; idx<nDir; idx++){
   if(M[idx] > max){
      max = M[idx];
      maxPixIdx = idx;
   }
}
*/
maxPixIdx = rank[0];
cout<<"max: "<<M[rank[0]]<<endl;
summary->setMaxPixInfo(rank[0], M[rank[0]]);

/*
 * Transferring map data
 */

cout<<"Transferring maxPixIdx map data...\n";
//mapData[0] = &M[rank[0] / nDir];
for(int dir=0; dir<nDir; dir++) mapData[dir] = M[rank[0] / nDir + dir];

/*
 * Compute likelihood and p value of skymap compared to reference map fits
 */

if(settings->computeLLHAndPValue == 1){
cout<<"Computing map likelihood and p-value w.r.t. the referenc map..."<<endl;

double pValue, likelihood;
err = computeMapLikelihoodAndPValue(/*summary->onion*/nDir, /*summary->onion*/nLayer, settings->referenceMapFitFunc/*.c_str()*/, settings->referenceMapFitFile/*.c_str()*/, M, likelihood, pValue);
summary->setLikelihoodAndPValue(likelihood, pValue);
cout<<"LLH: "<<likelihood<<" P Value: "<<pValue<<endl;
}
/*
TH1F *likelihoodDist =  new TH1F("likelihoodDist","likelihoodDist",1000, -5000,50000);
TH1F *pValueDist     =  new TH1F("pValueDist","pValueDist",1000, -1000,1000);

likelihoodDist->Fill(likelihood);
pValueDist->Fill(pValue);

likelihoodDist->Draw();
cvs.SaveAs("testLikelihoodDist.C");
pValueDist->Draw();
cvs.SaveAs("testPValueDist.C");
*/
/*
 * Write FITS file
 */

cout<<"Creating Healpix map and writing to FITS....\n";
//arr<float> MArr = arr<float>(&M[0], (size_t)nDir);
arr<float> MArr = arr<float>(&M[rank[0] / nDir], (size_t)nDir);
Healpix_Map<float> skyMap = Healpix_Map<float>(MArr, HEALPIX_ORDERING);

fitshandle fitsOut;
//#ifdef CSW
//char filename[] = "testCSWSkyMap.fits";
//#else
//char filename[] = "testXCorrSkyMap.fits";
//#endif
remove(filename);
fitsOut.create(filename);

write_Healpix_map_to_fits(fitsOut, skyMap, PLANCK_FLOAT32);
cout<<"Healpix map written\n";
/*
arr<float> MEachLayerArr;

for(int i=0; i<nLayer; i++){

   MEachLayerArr = arr<float>(&M[i*nDir], (size_t)nDir);
   skyMap = Healpix_Map<float>(MEachLayerArr, RING);

   char layerMap[200];
   sprintf(layerMap,"layer_%d_skymap.fits",i);
   remove(layerMap);
   fitsOut.create(layerMap);
   write_Healpix_map_to_fits(fitsOut, skyMap, PLANCK_FLOAT32);

}
*/

/*
 * Deallocate memories
 */

cout<<"Deallocating memories...\n";
clReleaseMemObject(recoDelaysBuffer);
clReleaseMemObject(intensityRBuffer);
clReleaseMemObject(intensityCBuffer);
clReleaseMemObject(voltsFlatBuffer);
clReleaseMemObject(xCorrCBuffer);
clReleaseMemObject(xCorrRBuffer);
clReleaseMemObject(xCorrTimeBuffer);
clReleaseMemObject(xCorrEnvBuffer);
clReleaseMemObject(sqrtWfPwrBuffer);
clReleaseMemObject(CijBuffer);
clReleaseMemObject(MBuffer);
clReleaseMemObject(maxPixIdxEachLayerBuffer);
clReleaseMemObject(maxPixCoherenceEachLayerBuffer);
clReleaseMemObject(beginTimeByChannelBuffer);
clReleaseMemObject(wfNBinsBuffer);
free(voltsFlat);
//free(volts);
//free(recoDelays);
free(intensity_data_r);
free(intensity_data_c);
free(xCorr_data_r);
free(xCorr_data_c);
free(xCorrTime);
free(sqrtWfPwr);
free(Cij);
free(M);
free(rank);
free(topMaxPixIdx);
free(topMaxPixCoherence);
free(maxPixIdxEachLayer);
free(maxPixCoherenceEachLayer);
free(beginTimeByChannel_temp);
free(wfNBins_temp);
cout<<"Memories deallocated\n";

return maxPixIdx;
}

int reconstruct3DXCorrEnvelopeGetMaxPix_ZoomMode(recoSettings *settings, vector<TGraph *>& cleanEvent, recoEnvData *clEnv,
                const float stationCenterDepth, const vector<vector<double> >& antLocation,
                float *recoDelays, float *recoDelays_V, float *recoDelays_H,
                const int *chanMask, recoData *summary, char *filename)
{

cout<<"Entered reconstruct3DXCorrEnvelopeGetMaxPix_ZoomMode method\n";
int nSamp;
//int nAnt = (int)cleanEvent.size()/2; // Divide by 2 for only one polarization
int nAnt;// = (int)cleanEvent.size();
string pol = string(settings->recoPolType);
if( pol == "vpol" || pol == "hpol" ) nAnt = (int)cleanEvent.size()/2; else  nAnt = (int)cleanEvent.size();
int unmaskedNChan=0;
for(int ch=0; ch<(int)cleanEvent.size(); ch++) unmaskedNChan+=chanMask[ch];
cout<<"unmaskedNChan: "<<unmaskedNChan<<" nAnt: "<<nAnt<<endl;
float wInt;
int dataType = settings->dataType;
if( dataType == 0 ) wInt = 0.5f; //AraSim event
else if( dataType == 1 ){ //real event
if( pol == "vpol" ) wInt = 0.4f;
else if ( pol == "hpol" ) wInt = 0.625f;
else if ( pol == "both" ) wInt = 0.5f;
else { cerr<<"recoPolType undefined\n"; return -1; }
} else {
cerr<<"dataType undefined\n"; return -1; }

//int nSideExp = settings->nSideExp;
//int nDir   = /*summary->onion*/12 * pow(2,nSideExp) * pow(2,nSideExp);;
//int nLayer = /*summary->onion*/settings->nLayer;
//printf("nDir: %d nLayer: %d\n",nDir,nLayer);
int nSideExp, nDir, nLayer;

/*
 * Loading voltsFlat array
 */

double t, v;
float *voltsFlat;

if( pol == "vpol" ){

   /* Using the 1st vpol wf for nSamp. Should make sure all vpol channels have the same nSamp. FIXME */
   nSamp = cleanEvent[0]->GetN();
   voltsFlat = (float*)calloc(nAnt*nSamp, sizeof(float));
   if(voltsFlat == NULL){ cerr<<"Null pointer to voltsFlat\n"; return -1; }

   for(int ch=0; ch<nAnt; ch++){
       if( chanMask[ch] ){
         for(int s=0; s<nSamp; s++){
            cleanEvent[ch]->GetPoint(s,t,v);
            /* Bartlett window applied in main analysis code */
            //voltsFlat[nSamp*ch + s] = static_cast<float>(v * FFTtools::bartlettWindow(s, nSamp));
            voltsFlat[nSamp*ch + s] = static_cast<float>(v);
         }
      } else {
         for(int s=0; s<nSamp; s++){
            voltsFlat[nSamp*ch + s] = 0.f;
         }
      }
   }

} else if ( pol == "hpol" ){

   /* Using the 1st hpol wf for nSamp. Should make sure all hpol channels have the same nSamp. FIXME */
   nSamp = cleanEvent[nAnt]->GetN();
   voltsFlat = (float*)calloc(nAnt*nSamp, sizeof(float));
   if(voltsFlat == NULL){ cerr<<"Null pointer to voltsFlat\n"; return -1; }

   for(int ch=0; ch<nAnt; ch++){
       if( chanMask[ch+nAnt] ){
         for(int s=0; s<nSamp; s++){
            cleanEvent[ch+nAnt]->GetPoint(s,t,v);
            /* Bartlett window applied in main analysis code */
            //voltsFlat[nSamp*ch + s] = static_cast<float>(v * FFTtools::bartlettWindow(s, nSamp));
            voltsFlat[nSamp*ch + s] = static_cast<float>(v);
         }
      } else {
         for(int s=0; s<nSamp; s++){
            voltsFlat[nSamp*ch + s] = 0.f;
         }
      }
   }

} else if ( pol == "both" ){

   /* Using the 1st vpol wf for nSamp. Should make sure all hpol channels have the same nSamp. FIXME */
   nSamp = cleanEvent[0]->GetN();
   voltsFlat = (float*)calloc(nAnt*nSamp, sizeof(float));
   if(voltsFlat == NULL){ cerr<<"Null pointer to voltsFlat\n"; return -1; }

   for(int ch=0; ch<nAnt; ch++){
       if( chanMask[ch] ){
         for(int s=0; s<nSamp; s++){
            cleanEvent[ch]->GetPoint(s,t,v);
            /* Bartlett window applied in main analysis code */
            //voltsFlat[nSamp*ch + s] = static_cast<float>(v * FFTtools::bartlettWindow(s, nSamp));
            voltsFlat[nSamp*ch + s] = static_cast<float>(v);
         }
      } else {
         for(int s=0; s<nSamp; s++){
            voltsFlat[nSamp*ch + s] = 0.f;
         }
      }
   }

} else {
   cerr<<"recoPolType not defined\n"; return -1;
}
cout<<"voltsFlat loaded\n";

/*
 * Preparation for OUT_OF_PLACE transforms
 */

int interlvHermOutputSize = 2*(1 + nSamp/2); //Hermitian layout
int interlvOutputSize     = 2*nSamp;         //Not Hermitian layout
int planarHermOutputSize = (1 + nSamp/2);    //Hermitian planar layout

float *intensity_data_r, *intensity_data_c;
intensity_data_r = (float*)calloc(nAnt*planarHermOutputSize, sizeof(float));
intensity_data_c = (float*)calloc(nAnt*planarHermOutputSize, sizeof(float));

/*
 * FFT library related declarations
 */
//clfftPlanHandle planHandle;
clfftDim fftDim = CLFFT_1D;
size_t clLengths[1]   = {nSamp};
size_t clInStride[1]  = {1};
size_t clOutStride[1] = {1};
size_t inDist  = nSamp;
size_t outDist   = planarHermOutputSize;
size_t batchSize = nAnt;

/*
 * Set up clFFT
 */
//clfftSetupData fftSetup;
//err = clfftInitSetupData(&fftSetup);
//err = clfftSetup(&fftSetup);

/*
 * Prepare plan
 */
cout<<"Preparing plan...\n";
int err;
err = clfftCreateDefaultPlan(&clEnv->planHandle, clEnv->context, fftDim, clLengths);
err = clfftSetPlanPrecision(clEnv->planHandle, CLFFT_SINGLE);
err = clfftSetLayout(clEnv->planHandle, CLFFT_REAL, CLFFT_HERMITIAN_PLANAR);
err = clfftSetPlanScale(clEnv->planHandle, CLFFT_FORWARD, 1.f);
err = clfftSetPlanBatchSize(clEnv->planHandle, batchSize);
err = clfftSetPlanInStride(clEnv->planHandle, fftDim, clInStride);
err = clfftSetPlanOutStride(clEnv->planHandle, fftDim, clOutStride);
err = clfftSetPlanDistance(clEnv->planHandle, inDist, outDist);
err = clfftSetResultLocation(clEnv->planHandle, CLFFT_OUTOFPLACE);
err = clfftBakePlan(clEnv->planHandle, 1, &clEnv->queue, NULL, NULL);
cout<<"Plan prepared\n";
/* The plan is now ready to be executed */
cl_mem voltsFlatBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*nSamp*sizeof(float), NULL, &err);
err = clEnqueueWriteBuffer(clEnv->queue, voltsFlatBuffer, CL_TRUE, 0, nAnt*nSamp*sizeof(float), voltsFlat, 0, NULL, NULL);

cl_mem intensityRBuffer= clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*planarHermOutputSize*sizeof(float),
                         NULL, &err);
cl_mem intensityCBuffer= clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nAnt*planarHermOutputSize*sizeof(float),
                         NULL, &err);
cl_mem outBuffers[2] = {intensityRBuffer, intensityCBuffer};

cout<<"Enqueueing FFT\n";
err = clfftEnqueueTransform(clEnv->planHandle, CLFFT_FORWARD, 1, &clEnv->queue, 0, NULL, NULL,
                            &voltsFlatBuffer, outBuffers, NULL);
err = clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, outBuffers[0], CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_r,
                    0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, outBuffers[1], CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_c,
                    0, NULL, NULL);
cout<<"FFT done\n";
/*
 * Clean up CLFFT
 */
cout<<"Destroying plan...\n";
err = clfftDestroyPlan(&clEnv->planHandle);
//clfftTeardown();
cout<<"Plan destroyed\n";


#ifdef bandpass
/*
 * Bandpass signals
 */
cout<<"Preparing bandpass filter..."<<endl;
float freqBin = 1e3 / (wInt * (float)nSamp); // wInt in ns. 1e3 for MHz
float lowFreq = 200.f;
float highFreq= 450.f;

clSetKernelArg(clEnv->bandPassFilter, 0, sizeof(cl_mem), &intensityRBuffer);
clSetKernelArg(clEnv->bandPassFilter, 1, sizeof(cl_mem), &intensityCBuffer);
clSetKernelArg(clEnv->bandPassFilter, 2, sizeof(float),  &freqBin);
clSetKernelArg(clEnv->bandPassFilter, 3, sizeof(float),  &lowFreq);
clSetKernelArg(clEnv->bandPassFilter, 4, sizeof(float),  &highFreq);

unsigned int dim = 2;
size_t bandPassWorkSize[2] = {(size_t)nAnt, (size_t)planarHermOutputSize};
clEnqueueNDRangeKernel(clEnv->queue, clEnv->bandPassFilter, dim, NULL, bandPassWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, intensityRBuffer, CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_r,
                    0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, intensityCBuffer, CL_TRUE, 0, sizeof(float)*nAnt*planarHermOutputSize, intensity_data_c,
                    0, NULL, NULL);

cout<<"Bandpass filter done"<<endl;
#endif

/*
 * Cross-correlate wfs
 */

const int nBaseline = nAnt*nAnt;
cout<<"nBaseline: "<<nBaseline<<endl;

float *xCorr_data_r = (float*)calloc(nBaseline*planarHermOutputSize, sizeof(float));
float *xCorr_data_c = (float*)calloc(nBaseline*planarHermOutputSize, sizeof(float));

cl_mem xCorrRBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                    sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_r, &err);

cl_mem xCorrCBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                    sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_c, &err);

clSetKernelArg(clEnv->xCorrWf, 0, sizeof(cl_mem), &xCorrRBuffer);
clSetKernelArg(clEnv->xCorrWf, 1, sizeof(cl_mem), &xCorrCBuffer);
clSetKernelArg(clEnv->xCorrWf, 2, sizeof(cl_mem), &intensityRBuffer);
clSetKernelArg(clEnv->xCorrWf, 3, sizeof(cl_mem), &intensityCBuffer);
clSetKernelArg(clEnv->xCorrWf, 4, sizeof(int),    &nAnt);

unsigned int workDim = 2;
size_t globalWorkSize[2] = {(size_t)nBaseline, (size_t)planarHermOutputSize};
clEnqueueNDRangeKernel(clEnv->queue, clEnv->xCorrWf, workDim, NULL, globalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, xCorrRBuffer, CL_TRUE, 0, sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_r, 0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, xCorrCBuffer, CL_TRUE, 0, sizeof(float)*nBaseline*planarHermOutputSize, xCorr_data_c, 0, NULL, NULL);
cout<<"Done xCorrWf\n";

/*
 * Inverse FFT cross-correlation
 */

inDist  = planarHermOutputSize;
outDist = nSamp;
batchSize = nBaseline;

/*
 * Prepare plan
 */

err = clfftCreateDefaultPlan(&clEnv->planHandle, clEnv->context, fftDim, clLengths);
err = clfftSetPlanPrecision(clEnv->planHandle, CLFFT_SINGLE);
err = clfftSetLayout(clEnv->planHandle, CLFFT_HERMITIAN_PLANAR, CLFFT_REAL);
err = clfftSetPlanScale(clEnv->planHandle, CLFFT_BACKWARD, 1.f/(float)nSamp);
err = clfftSetPlanBatchSize(clEnv->planHandle, batchSize);
err = clfftSetPlanInStride(clEnv->planHandle, fftDim, clInStride);
err = clfftSetPlanOutStride(clEnv->planHandle, fftDim, clOutStride);
err = clfftSetPlanDistance(clEnv->planHandle, inDist, outDist);
err = clfftSetResultLocation(clEnv->planHandle, CLFFT_OUTOFPLACE);
err = clfftBakePlan(clEnv->planHandle, 1, &clEnv->queue, NULL, NULL);

/* The plan is now ready to be executed */

cl_mem inBuffers[2] = {xCorrRBuffer, xCorrCBuffer};

float *xCorrTime = (float*)calloc(nBaseline*nSamp, sizeof(float));
cl_mem xCorrTimeBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE, nBaseline*nSamp*sizeof(float),
                           NULL, &err);
err = clfftEnqueueTransform(clEnv->planHandle, CLFFT_BACKWARD, 1, &clEnv->queue, 0, NULL, NULL, inBuffers, &xCorrTimeBuffer, NULL);
err = clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, xCorrTimeBuffer, CL_TRUE, 0, sizeof(float)*nBaseline*nSamp, xCorrTime, 0, NULL, NULL);

/*
 * Clean up FFT
 */

err = clfftDestroyPlan(&clEnv->planHandle);
//clfftTeardown();
cout<<"Done inverse FFT\n";

/*
 * Get Hilbert transform of XCorr function
 */

TCanvas cvs("cvs","cvs",800,600);
//cvs.Divide(1,2);
char envelopename[200];

float dt[nSamp];
float xCorrValue[nSamp];
double t_temp, v_temp;
double t_temp2, v_temp2;
int ant1, ant2;
float plusMinusTime = 25.f; //record +-25 ns around peak in XCorr function
//int plusMinusRange = (int)(plusMinusTime / wInt);
//static TGraph *sillygr[64];// = new TGraph();
//static TGraph *sillygr = (TGraph*)malloc(nBaseline*sizeof(TGraph));
/*
for(int i=0; i<nBaseline; i++){

   sillygr[i]=new TGraph();

}
*/
//static int lock;
//char sillygrname[200];

for(int baseline=0; baseline<nBaseline; baseline++){


   ant1 = baseline / nAnt;
   ant2 = baseline % nAnt;

   for(int s=0; s<nSamp; s++){

   dt[s] = wInt*s;
   xCorrValue[s] = xCorrTime[nSamp*baseline + s];

   }

   TGraph *xCorrGraph = new TGraph(nSamp, dt, xCorrValue);

   //cvs.cd(1);
   //xCorrGraph->Draw("AL");

   TGraph* envelope = FFTtools::getHilbertEnvelope( xCorrGraph );
/*
   sprintf(envelopename,"xCorrEnvelope_2014_A3_burn_RF_chan%d_%d.C", ant1, ant2);
   cvs.cd();
   //envelope->Draw("AL");
   envelope->SetLineColor(kRed);
   xCorrGraph->Draw("AL");
   envelope->Draw("Lsame");
   cvs.SaveAs(envelopename);
*/
/*
   if( ant1 == 0 && ant2 == 3){
   cvs.cd(1);
   envelope->Draw("AL");
   //xCorrGraph->Draw("AL");
   cvs.SaveAs("xCorrEnvelope_chan0_3.C");
   }
*/
/*
   if( ant1 == 3 && ant2 == 0){
   cvs.cd(2);
   envelope->Draw("AL");
   }

   cvs.SaveAs("xCorrEnvelope_chan0_3.C");
*/
/*
   if( ant1 == 0 && chanMask[ant1] ){
      cout<<"recording XCorr values around peak of each baseline\n";

      if( ant2 == 1 && chanMask[ant2] ){
      stackXCorrAroundPeak(envelope, xCorrAroundPeakHist[0], plusMinusTime);
      } else if
      ( ant2 == 3 && chanMask[ant2] ){
      stackXCorrAroundPeak(envelope, xCorrAroundPeakHist[1], plusMinusTime);
      } else if
      ( ant2 == 7 && chanMask[ant2] ){
      stackXCorrAroundPeak(envelope, xCorrAroundPeakHist[2], plusMinusTime);
      }
   }
*/
   //if(ant1==1 && ant2==5){
   //if(ant1 < ant2 ){
   //for(int s=0; s<nSamp; s++){

   //if(ant1 == 0  && ant2 == 4 ){
   //xCorrGraph->GetPoint(s,t_temp,v_temp);
   //if(lock == 0){
   //if(envelopeSum->GetN() != 0 ) cerr<<"Both TGraphs should be uninitialized!"<<endl;
/*
   if(sillygr[baseline]->GetN() == 0){

      for(int s=0; s<nSamp; s++){
      xCorrGraph->GetPoint(s,t_temp,v_temp);
      sillygr[baseline]->SetPoint(s, t_temp, v_temp);
      }

   } else {

      for(int s=0; s<nSamp; s++){
      xCorrGraph->GetPoint(s,t_temp,v_temp);
      sillygr[baseline]->GetPoint(s, t_temp2, v_temp2);
      if(t_temp != t_temp2) cerr<<"t_temp != t_temp2 !!\n";
      sillygr[baseline]->SetPoint(s, t_temp2, v_temp+v_temp2);
      }
   }
*/
   //}

   for(int s=0; s<nSamp; s++){

   envelope->GetPoint(s,t_temp,v_temp);
   xCorrTime[nSamp*baseline + s] = static_cast<float>(v_temp);

   }
   //if(lock == 0) lock=1;
/*
  cvs.cd();
  sprintf(sillygrname,"xCorrSumGraph_chan%d_%d.C",ant1,ant2);
  sillygr[baseline]->Draw("AL");
  cvs.SaveAs(sillygrname);
*/

  delete xCorrGraph;
  delete envelope;

}

cl_mem xCorrEnvBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                       sizeof(float)*nBaseline*nSamp,
                                       xCorrTime, &err);

/* Compute the square root of each channel's total wf power */

float *sqrtWfPwr = (float*)calloc(nAnt, sizeof(float));
float pwr=0.f;
for(int ant=0; ant<nAnt; ant++){
   pwr = 0.;
   for(int s=0; s<nSamp; s++){
      pwr += (voltsFlat[ant*nSamp + s] * voltsFlat[ant*nSamp + s]);
   }
   sqrtWfPwr[ant] = sqrt(pwr);
   //sqrtWfPwr[ant] = 1.f;
}
cout<<"Done sqrtWfPwr\n";

cl_mem sqrtWfPwrBuffer  = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt,
                          sqrtWfPwr, &err);

/* At this point all correlation coefficients are ready to be read out according to the delays. For zoom search mode, we will loop over
 * skymaps (whole or patches) from low to high resolutions, and read out the relevant coefficients.
 */

cl_mem recoDelaysBuffer;
float max=0.f;
int maxPixIdx;
int *rank;
nSideExp = settings->nSideExpStart;
nDir = 12 * pow(2, nSideExp) * pow(2, nSideExp);
nLayer = settings->nLayer;
Healpix_Onion *onion = new Healpix_Onion(nSideExp, nLayer);
//Healpix_Onion onion(nSideExp, nLayer);

if(pol == "both" ){
recoDelays  = (float*)malloc(nLayer*nDir*nAnt*sizeof(float));
recoDelays_V= (float*)malloc(nLayer*nDir*(nAnt/2)*sizeof(float));
recoDelays_H= (float*)malloc(nLayer*nDir*(nAnt/2)*sizeof(float));
} else {
recoDelays  = (float*)malloc(nLayer*nDir*nAnt*2*sizeof(float));
recoDelays_V= (float*)malloc(nLayer*nDir*nAnt*sizeof(float));
recoDelays_H= (float*)malloc(nLayer*nDir*nAnt*sizeof(float));
}

/* compute reco delays for skymap of nSideExpStart */
// N.B. nAnt could be 8 or 16 depending on recoPolType

if(settings->iceModel == 1){
   if( pol == "both" ){
   err = computeRecoDelaysWithConstantN(nAnt, -1.f*stationCenterDepth, antLocation,
                                        //radius, nSideExp,
                                        onion, recoDelays, recoDelays_V, recoDelays_H);
   } else {
   err = computeRecoDelaysWithConstantN(nAnt*2, -1.f*stationCenterDepth, antLocation,
                                        //radius, nSideExp,
                                        onion, recoDelays, recoDelays_V, recoDelays_H);
   }
} else if(settings->iceModel == 0){
   if( pol == "both" ){
   err = compute3DRecoDelaysWithRadioSpline(nAnt, -1.f*stationCenterDepth, antLocation,
                                            onion, recoDelays, recoDelays_V, recoDelays_H);
   } else {
   err = compute3DRecoDelaysWithRadioSpline(nAnt*2, -1.f*stationCenterDepth, antLocation,
                                            onion, recoDelays, recoDelays_V, recoDelays_H);
   }
} else { cerr<<"Undefined iceModel parameter\n"; return -1; }
if( err<0 ){ cerr<<"Error computing reco delays\n"; return -1; }


if(pol == "vpol" ){
recoDelaysBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt*nDir*nLayer,
                                  recoDelays_V, &err);
} else if (pol == "hpol" ) {
recoDelaysBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt*nDir*nLayer,
                                  recoDelays_H, &err);
} else if (pol == "both" ) {
recoDelaysBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt*nDir*nLayer,
                                  recoDelays, &err);
} else {
cerr<<"recoPolType not defined\n"; return -1;
}
cout<<"recoDelaysBuffer created\n";

/*
 * Compute cross-correlation coefficients Cij in each direction
 */

float *Cij = (float*)calloc(nLayer*nDir*nBaseline, sizeof(float));

cl_mem CijBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nLayer*nDir*nBaseline,
                   Cij, &err);

clSetKernelArg(clEnv->computeXCorrCoef, 0, sizeof(cl_mem), &CijBuffer);
clSetKernelArg(clEnv->computeXCorrCoef, 1, sizeof(cl_mem), &xCorrEnvBuffer);
clSetKernelArg(clEnv->computeXCorrCoef, 2, sizeof(cl_mem), &recoDelaysBuffer);
clSetKernelArg(clEnv->computeXCorrCoef, 3, sizeof(cl_mem), &sqrtWfPwrBuffer);
clSetKernelArg(clEnv->computeXCorrCoef, 4, sizeof(float),  &wInt);
clSetKernelArg(clEnv->computeXCorrCoef, 5, sizeof(int),    &nAnt);
clSetKernelArg(clEnv->computeXCorrCoef, 6, sizeof(int),    &nSamp);
workDim = 3;
size_t CijGlobalWorkSize[3] = {(size_t)nLayer, (size_t)nDir, (size_t)nBaseline};

err = clEnqueueNDRangeKernel(clEnv->queue, clEnv->computeXCorrCoef, workDim, NULL, CijGlobalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, CijBuffer, CL_TRUE, 0, sizeof(float)*nLayer*nDir*nBaseline, Cij, 0, NULL, NULL);

cout<<"Done computeXCorrCoef\n";

/*
 * Sum Cij's of all baselines in each reco direction to obtain coherence M(r-hat)
 */

float *M = (float*)calloc(nLayer*nDir, sizeof(float));
cl_mem MBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nLayer*nDir,
                  M, &err);
/*
clSetKernelArg(clEnv->computeCoherence, 0, sizeof(cl_mem), &MBuffer);
clSetKernelArg(clEnv->computeCoherence, 1, sizeof(cl_mem), &CijBuffer);
clSetKernelArg(clEnv->computeCoherence, 2, sizeof(int),    &nBaseline);
*/
clSetKernelArg(clEnv->computeNormalizedCoherence, 0, sizeof(cl_mem), &MBuffer);
clSetKernelArg(clEnv->computeNormalizedCoherence, 1, sizeof(cl_mem), &CijBuffer);
clSetKernelArg(clEnv->computeNormalizedCoherence, 2, sizeof(int),    &nBaseline);

workDim = 2;
size_t MGlobalWorkSize[2] = {(size_t)nLayer, (size_t)nDir};

//clEnqueueNDRangeKernel(clEnv->queue, clEnv->computeCoherence, workDim, NULL, MGlobalWorkSize, NULL, 0, NULL, NULL);
clEnqueueNDRangeKernel(clEnv->queue, clEnv->computeNormalizedCoherence, workDim, NULL, MGlobalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, MBuffer, CL_TRUE, 0, sizeof(float)*nLayer*nDir, M, 0, NULL, NULL);

cout<<"Done computeNormalizedCoherence\n";

/*
 * Loop over M to find the max coherence and its pix index
 */

//float max=0.f;
//int maxPixIdx;
rank = (int*)calloc(nLayer*nDir, sizeof(int));

TMath::Sort(nLayer*nDir, M, rank);

maxPixIdx = rank[0];
int last2DMaxPixIdx = maxPixIdx % nDir;
//int finalMaxPixIdx = maxPixIdx;
/*
   delete onion;
   clReleaseMemObject(recoDelaysBuffer);
   clReleaseMemObject(CijBuffer);
   clReleaseMemObject(MBuffer);
   free(recoDelays);
   free(recoDelays_V);
   free(recoDelays_H);
   free(Cij);
   free(M);
   free(rank);
*/

float *delays;


for(nSideExp = settings->nSideExpStart+1; nSideExp <= settings->nSideExpEnd; nSideExp++){

   cout<<"nSideExp: "<<nSideExp<<endl;

   delete onion;
   clReleaseMemObject(recoDelaysBuffer);
   clReleaseMemObject(CijBuffer);
   clReleaseMemObject(MBuffer);
   //clFinish(clEnv->queue);
   free(recoDelays);
   free(recoDelays_V);
   free(recoDelays_H);
   free(Cij);
   free(M);
   free(rank);
   recoDelays = NULL;
   recoDelays_V = NULL;
   recoDelays_H = NULL;
   Cij = NULL;
   M = NULL;
   rank = NULL;

   onion = new Healpix_Onion(nSideExp, nLayer);
   nDir = NPIX_NESTED; //defined in recoTools.h


   //finalMaxPixIdx *= NPIX_NESTED;

   if(pol == "both" ){
   recoDelays  = (float*)malloc(nLayer*nDir*nAnt*sizeof(float));
   recoDelays_V= (float*)malloc(nLayer*nDir*(nAnt/2)*sizeof(float));
   recoDelays_H= (float*)malloc(nLayer*nDir*(nAnt/2)*sizeof(float));
   } else {
   recoDelays  = (float*)malloc(nLayer*nDir*nAnt*2*sizeof(float));
   recoDelays_V= (float*)malloc(nLayer*nDir*nAnt*sizeof(float));
   recoDelays_H= (float*)malloc(nLayer*nDir*nAnt*sizeof(float));
   }



   /* compute reco delays for zoomed patch of skymap */



   if(settings->iceModel == 1){
      if(pol == "both"){
      err = computeZoomedRecoDelaysWithConstantN(nAnt, -1.f*stationCenterDepth, antLocation,
                                           //radius, nSideExp,
                                           onion, recoDelays, recoDelays_V, recoDelays_H,
                                           last2DMaxPixIdx, nDir);
      } else {
      err = computeZoomedRecoDelaysWithConstantN(nAnt*2, -1.f*stationCenterDepth, antLocation,
                                           //radius, nSideExp,
                                           onion, recoDelays, recoDelays_V, recoDelays_H,
                                           last2DMaxPixIdx, nDir);
      }
   } else if(settings->iceModel == 0){
      if(pol == "both"){
      err = compute3DZoomedRecoDelaysWithRadioSpline(nAnt, -1.f*stationCenterDepth, antLocation,
                                               onion, recoDelays, recoDelays_V, recoDelays_H,
                                               last2DMaxPixIdx, nDir);
      } else {
      err = compute3DZoomedRecoDelaysWithRadioSpline(nAnt*2, -1.f*stationCenterDepth, antLocation,
                                               onion, recoDelays, recoDelays_V, recoDelays_H,
                                               last2DMaxPixIdx, nDir);
      }
   } else { cerr<<"Undefined iceModel parameter\n"; return -1; }
   if( err<0 ){ cerr<<"Error computing reco delays\n"; return -1; }

   if(pol == "vpol" ){
   //delays = recoDelays_V;
   recoDelaysBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt*nDir*nLayer,
                                     recoDelays_V, &err);
   } else if (pol == "hpol" ) {
   //delays = recoDelays_H;
   recoDelaysBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt*nDir*nLayer,
                                     recoDelays_H, &err);
   } else if (pol == "both" ) {
   //delays = recoDelays;
   recoDelaysBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nAnt*nDir*nLayer,
                                     recoDelays, &err);
   } else {
   cerr<<"recoPolType not defined\n"; return -1;
   }
   //clFinish(clEnv->queue);
   //cout<<"recoDelaysBuffer created\n";



   /*
    * Compute cross-correlation coefficients Cij in each direction
    */



   Cij = (float*)calloc(nLayer*nDir*nBaseline, sizeof(float));
/*
   int ant1, ant2;
   for(int layer=0; layer<nLayer; layer++){
      for(int dir=0; dir<nDir; dir++){
         for(int baseline=0; baseline<nBaseline; baseline++){

         ant1 = baseline / nAnt;
         ant2 = baseline % nAnt;
         if( (layer*nDir*nAnt + dir*nAnt + ant1) >= nLayer*nDir*nBaseline || (layer*nDir*nAnt + dir*nAnt + ant2) >= nLayer*nDir*nBaseline ) printf("****** Warning! out of range!!! ****\n");

         if( sqrtWfPwr[ant1] != 0.f && sqrtWfPwr[ant2] != 0.f ){
         if(delays[layer*nDir*nAnt + dir*nAnt + ant1] > -1e9 && delays[layer*nDir*nAnt + dir*nAnt + ant2] > -1e9 ){
         int shiftBin = (delays[layer*nDir*nAnt + dir*nAnt + ant1] - delays[layer*nDir*nAnt + dir*nAnt + ant2]) / wInt;
  */       /* The default circular  FFT of cross-correlation outputs in wrap-around order */
         /* The correlation at -i is in r_N-i */
/*         if( shiftBin < 0 ) shiftBin += nSamp;
         Cij[layer*nDir*nBaseline + dir*nBaseline + baseline] = xCorrTime[baseline*nSamp + shiftBin] / (sqrtWfPwr[ant1] * sqrtWfPwr[ant2]);
         } else {
         //printf("Not both delays exists!\n");
         //if(delays[gid0*nAnt + ant1] != delays[gid0*nAnt + ant2] )
         //   if(delays[gid0*nAnt + ant1] > 0 || delays[gid0*nAnt + ant2] >0) printf("ant1 delay: %f ant2 delay: %f\n",delays[gid0*nAnt + ant1],delays[gid0*nAnt + ant2]    );
         Cij[layer*nDir*nBaseline + dir*nBaseline + baseline] = 0.f;
         }
         } else {
         Cij[layer*nDir*nBaseline + dir*nBaseline + baseline] = 0.f;
         }
         }
      }
   }
*/
   CijBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nLayer*nDir*nBaseline,
                      Cij, &err);
   clSetKernelArg(clEnv->computeXCorrCoef, 0, sizeof(cl_mem), &CijBuffer);
   clSetKernelArg(clEnv->computeXCorrCoef, 1, sizeof(cl_mem), &xCorrEnvBuffer);
   clSetKernelArg(clEnv->computeXCorrCoef, 2, sizeof(cl_mem), &recoDelaysBuffer);
   clSetKernelArg(clEnv->computeXCorrCoef, 3, sizeof(cl_mem), &sqrtWfPwrBuffer);
   clSetKernelArg(clEnv->computeXCorrCoef, 4, sizeof(float),  &wInt);
   clSetKernelArg(clEnv->computeXCorrCoef, 5, sizeof(int),    &nAnt);
   clSetKernelArg(clEnv->computeXCorrCoef, 6, sizeof(int),    &nSamp);
   workDim = 3;
   //size_t CijGlobalWorkSize[3] = {(size_t)nLayer, (size_t)nDir, (size_t)nBaseline};
   CijGlobalWorkSize[1] = (size_t)nDir; //reset to new nDir value
   cout<<CijGlobalWorkSize[0]<<" "<<CijGlobalWorkSize[1]<<" "<<CijGlobalWorkSize[2]<<endl;
   err = clEnqueueNDRangeKernel(clEnv->queue, clEnv->computeXCorrCoef, workDim, NULL, CijGlobalWorkSize, NULL, 0, NULL, NULL);
   clFinish(clEnv->queue);
   clEnqueueReadBuffer(clEnv->queue, CijBuffer, CL_TRUE, 0, sizeof(float)*nLayer*nDir*nBaseline, Cij, 0, NULL, NULL);

   cout<<"Done computeXCorrCoef\n";



   /*
    * Sum Cij's of all baselines in each reco direction to obtain coherence M(r-hat)
    */

   M = (float*)calloc(nLayer*nDir, sizeof(float));
/*
    M[gid0*nDir + gid1] = 0.f;
   for(int i=0; i<nAnt; i++){
      for(int j=i+1; j<nAnt; j++){
       M[gid0*nDir + gid1] += Cij[gid0*nDir*nBaseline + gid1*nBaseline + i*nAnt + j];
       }
    }
    //M[gid0] = sqrt(M[gid0]*M[gid0]);
    M[gid0*nDir + gid1] /= (float)((nAnt*(nAnt-1))/2); //Normalize M by the number of baselines summed.
                                                       //This makes the comparison among events with different numbers of reco channel
*/                                                       //more natural. 16.15.16


   MBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nLayer*nDir,
                  M, &err);


   /*
   clSetKernelArg(clEnv->computeCoherence, 0, sizeof(cl_mem), &MBuffer);
   clSetKernelArg(clEnv->computeCoherence, 1, sizeof(cl_mem), &CijBuffer);
   clSetKernelArg(clEnv->computeCoherence, 2, sizeof(int),    &nBaseline);
   */



   clSetKernelArg(clEnv->computeNormalizedCoherence, 0, sizeof(cl_mem), &MBuffer);
   clSetKernelArg(clEnv->computeNormalizedCoherence, 1, sizeof(cl_mem), &CijBuffer);
   clSetKernelArg(clEnv->computeNormalizedCoherence, 2, sizeof(int),    &nBaseline);

   workDim = 2;
   //size_t MGlobalWorkSize[2] = {(size_t)nLayer, (size_t)nDir};
   MGlobalWorkSize[1] = (size_t)nDir; //reset to new nDir value

   //clEnqueueNDRangeKernel(clEnv->queue, clEnv->computeCoherence, workDim, NULL, MGlobalWorkSize, NULL, 0, NULL, NULL);
   clEnqueueNDRangeKernel(clEnv->queue, clEnv->computeNormalizedCoherence, workDim, NULL, MGlobalWorkSize, NULL, 0, NULL, NULL);
   clFinish(clEnv->queue);
   clEnqueueReadBuffer(clEnv->queue, MBuffer, CL_TRUE, 0, sizeof(float)*nLayer*nDir, M, 0, NULL, NULL);

   cout<<"Done computeNormalizedCoherence\n";



   /*
    * Loop over M to find the max coherence and its pix index
    */



   rank = (int*)calloc(nLayer*nDir, sizeof(int));

   TMath::Sort(nLayer*nDir, M, rank);

   //maxPixIdx = rank[0]%nDir + last2DMaxPidxIdx*NPIX_NESTED + onion->nDir*(rank[0]/nDir);
   last2DMaxPixIdx =  rank[0]%nDir + last2DMaxPixIdx*NPIX_NESTED;
   maxPixIdx = last2DMaxPixIdx +  onion->nDir*(rank[0]/nDir);



/*
   delete onion;
   clReleaseMemObject(recoDelaysBuffer);
   clReleaseMemObject(CijBuffer);
   clReleaseMemObject(MBuffer);
   free(recoDelays);
   free(recoDelays_V);
   free(recoDelays_H);
   free(Cij);
   free(M);
   free(rank);
 */



  }



/*
int *topMaxPixIdx = (int*)calloc(summary->topN, sizeof(int));
float *topMaxPixCoherence = (float*)calloc(summary->topN, sizeof(float));

for(int i=0; i<summary->topN; i++){
   topMaxPixIdx[i] = rank[i];
   topMaxPixCoherence[i] = M[rank[i]];
}

summary->setTopMaxPixInfo(topMaxPixIdx, topMaxPixCoherence);

//int *rankEachLayer = (int*)calloc(nDir, sizeof(int));
int *maxPixIdxEachLayer = (int*)calloc(nLayer, sizeof(int));
float *maxPixCoherenceEachLayer = (float*)calloc(nLayer, sizeof(float));

cl_mem maxPixIdxEachLayerBuffer       = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(int)*nLayer,
                                        maxPixIdxEachLayer, &err);
cl_mem maxPixCoherenceEachLayerBuffer = clCreateBuffer(clEnv->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float)*nLayer,
                                        maxPixCoherenceEachLayer, &err);

clSetKernelArg(clEnv->getMaxPixInfoEachLayer, 0, sizeof(cl_mem), &maxPixIdxEachLayerBuffer);
clSetKernelArg(clEnv->getMaxPixInfoEachLayer, 1, sizeof(cl_mem), &maxPixCoherenceEachLayerBuffer);
clSetKernelArg(clEnv->getMaxPixInfoEachLayer, 2, sizeof(cl_mem), &MBuffer);
clSetKernelArg(clEnv->getMaxPixInfoEachLayer, 3, sizeof(int),    &nDir);

workDim = 1;
size_t maxPixInfoGlobalWorkSize = nLayer;

clEnqueueNDRangeKernel(clEnv->queue, clEnv->getMaxPixInfoEachLayer, workDim, NULL, &maxPixInfoGlobalWorkSize, NULL, 0, NULL, NULL);
clFinish(clEnv->queue);
clEnqueueReadBuffer(clEnv->queue, maxPixIdxEachLayerBuffer, CL_TRUE, 0, sizeof(int)*nLayer, maxPixIdxEachLayer, 0, NULL, NULL);
clEnqueueReadBuffer(clEnv->queue, maxPixCoherenceEachLayerBuffer, CL_TRUE, 0, sizeof(float)*nLayer, maxPixCoherenceEachLayer, 0, NULL, NULL);

summary->setMaxPixInfoEachLayer(settings, maxPixIdxEachLayer, maxPixCoherenceEachLayer);
*/
/*
for(int idx=0; idx<nDir; idx++){
   if(M[idx] > max){
      max = M[idx];
      maxPixIdx = idx;
   }
}
*/
//maxPixIdx = rank[0];
cout<<"max: "<<M[rank[0]]<<endl;
//summary->setMaxPixInfo(rank[0], M[rank[0]]);
summary->setMaxPixInfo(maxPixIdx, M[rank[0]]);

/*
 * Transferring map data
 */

//cout<<"Transferring maxPixIdx map data...\n";
//mapData[0] = &M[rank[0] / nDir];
//for(int dir=0; dir<nDir; dir++) mapData[dir] = M[rank[0] / nDir + dir];

/*
 * Compute likelihood and p value of skymap compared to reference map fits
 */

//if(settings->computeLLHAndPValue == 1){
//cout<<"Computing map likelihood and p-value w.r.t. the referenc map..."<<endl;

//double pValue, likelihood;
//err = computeMapLikelihoodAndPValue(/*summary->onion*/nDir, /*summary->onion*/nLayer, settings->referenceMapFitFunc.c_str(), settings->referenceMapFitFile.c_str(), M, likelihood, pValue);
//summary->setLikelihoodAndPValue(likelihood, pValue);
//cout<<"LLH: "<<likelihood<<" P Value: "<<pValue<<endl;
//}

/*
TH1F *likelihoodDist =  new TH1F("likelihoodDist","likelihoodDist",1000, -5000,50000);
TH1F *pValueDist     =  new TH1F("pValueDist","pValueDist",1000, -1000,1000);

likelihoodDist->Fill(likelihood);
pValueDist->Fill(pValue);

likelihoodDist->Draw();
cvs.SaveAs("testLikelihoodDist.C");
pValueDist->Draw();
cvs.SaveAs("testPValueDist.C");
*/
/*
 * Write FITS file
 */
/*
cout<<"Creating Healpix map and writing to FITS....\n";
//arr<float> MArr = arr<float>(&M[0], (size_t)nDir);
arr<float> MArr = arr<float>(&M[rank[0] / nDir], (size_t)nDir);
Healpix_Map<float> skyMap = Healpix_Map<float>(MArr, RING);

fitshandle fitsOut;
//#ifdef CSW
//char filename[] = "testCSWSkyMap.fits";
//#else
//char filename[] = "testXCorrSkyMap.fits";
//#endif
remove(filename);
fitsOut.create(filename);

write_Healpix_map_to_fits(fitsOut, skyMap, PLANCK_FLOAT32);
cout<<"Healpix map written\n";
*//*
arr<float> MEachLayerArr;

for(int i=0; i<nLayer; i++){

   MEachLayerArr = arr<float>(&M[i*nDir], (size_t)nDir);
   skyMap = Healpix_Map<float>(MEachLayerArr, RING);

   char layerMap[200];
   sprintf(layerMap,"layer_%d_skymap.fits",i);
   remove(layerMap);
   fitsOut.create(layerMap);
   write_Healpix_map_to_fits(fitsOut, skyMap, PLANCK_FLOAT32);

}
*/

/*
 * Deallocate memories
 */

cout<<"Deallocating memories...\n";
clReleaseMemObject(recoDelaysBuffer);
clReleaseMemObject(intensityRBuffer);
clReleaseMemObject(intensityCBuffer);
clReleaseMemObject(voltsFlatBuffer);
clReleaseMemObject(xCorrCBuffer);
clReleaseMemObject(xCorrRBuffer);
clReleaseMemObject(xCorrTimeBuffer);
clReleaseMemObject(xCorrEnvBuffer);
clReleaseMemObject(sqrtWfPwrBuffer);
clReleaseMemObject(CijBuffer);
clReleaseMemObject(MBuffer);
//clReleaseMemObject(maxPixIdxEachLayerBuffer);
//clReleaseMemObject(maxPixCoherenceEachLayerBuffer);
free(voltsFlat);
//free(volts);
free(recoDelays);
free(recoDelays_V);
free(recoDelays_H);
free(intensity_data_r);
free(intensity_data_c);
free(xCorr_data_r);
free(xCorr_data_c);
free(xCorrTime);
free(sqrtWfPwr);
free(Cij);
free(M);
free(rank);
//free(topMaxPixIdx);
//free(topMaxPixCoherence);
//free(maxPixIdxEachLayer);
//free(maxPixCoherenceEachLayer);
cout<<"Memories deallocated\n";

return maxPixIdx;
}

int plotMaxPix(const int nDir, int *maxPix, char *filename){

   arr<int> maxPixArr = arr<int>(&maxPix[0], (size_t)nDir);
   Healpix_Map<int> maxPixMap = Healpix_Map<int>(maxPixArr, HEALPIX_ORDERING);

   fitshandle fitsOut;
   //char filename[] = "testMaxPixMap.fits";
   remove(filename);
   fitsOut.create(filename);

   write_Healpix_map_to_fits(fitsOut, maxPixMap, PLANCK_INT8);
   cout<<"Healpix maxPixMap written\n";

   return 0;
}

int plotMaxPixZenAzi(const int nSideExp, int *maxPix, char *rootFilename){

   ifstream rootfile(rootFilename);
   TH1F *zenDist;
   TH1F *aziDist;

   if( !rootfile ){
   TFile fp(rootFilename, "NEW");
   zenDist = new TH1F("recoZenDist", "recoZenDist", 180, 0, 180);
   aziDist = new TH1F("recoAziDist", "recoAziDist", 360, 0, 360);
   zenDist->Write();
   aziDist->Write();
   fp.Close();
   }
   //else{
   TFile fp_2(rootFilename, "READ");
   zenDist = (TH1F*)fp_2.Get("recoZenDist");
   aziDist = (TH1F*)fp_2.Get("recoAziDist");
   //}

   int nSide = pow(2, nSideExp);
   Healpix_Base hpBase = Healpix_Base(nSide, HEALPIX_ORDERING, SET_NSIDE);
   int nDir  = hpBase.Npix();
   pointing pt;

   float zen, azi;

   for(int pix=0; pix<nDir; pix++){

   pt = hpBase.pix2ang( pix );
   zen = pt.theta*180.f / M_PI;;
   azi = pt.phi * 180.f / M_PI;

   if( maxPix[pix] != 0 ){
   zenDist->Fill(zen, (float)maxPix[pix]);
   aziDist->Fill(azi, (float)maxPix[pix]);
   }
   }

   TFile fp_3(rootFilename, "RECREATE");
   zenDist->Write();
   aziDist->Write();
   fp_2.Close();
   fp_3.Close();

   return 0;
}


int tearDown(recoEnvData *clEnv){
/*
   clReleaseKernel(clEnv->shiftWf);
   clReleaseKernel(clEnv->sumWf);
   clReleaseKernel(clEnv->wfPwr);
   clReleaseKernel(clEnv->xCorrWf);
   clReleaseKernel(clEnv->computeXCorrCoef);
   clReleaseKernel(clEnv->computeCoherence);
*/
   free(clEnv->kernels);
   clReleaseCommandQueue(clEnv->queue);
   clReleaseProgram(clEnv->program);
   clReleaseContext(clEnv->context);
   free(clEnv->devices);
   free(clEnv->platforms);

   return 0;
}

float getMeanDelay( vector<float>& solvedDelay ){

   //cout<<"solvedDelay size: "<<solvedDelay.size()<<endl;
   float mean = 0.f;
   for(int i=0; i<(int)solvedDelay.size(); i++){
   //cout<<"solvedDelay: "<<solvedDelay[i]<<endl;
   mean+=solvedDelay[i];
   }
   mean /= (float)solvedDelay.size();

   solvedDelay.clear(); //clear the vector once it's used

   return mean;
}

float getMeanDelay_passByValue( vector<float> solvedDelay ){

   //cout<<"solvedDelay size: "<<solvedDelay.size()<<endl;
   float mean = 0.f;
   for(int i=0; i<(int)solvedDelay.size(); i++){
   //cout<<"solvedDelay: "<<solvedDelay[i]<<endl;
   mean+=solvedDelay[i];
   }
   mean /= (float)solvedDelay.size();

   //solvedDelay.clear(); //clear the vector once it's used

   return mean;
}

int computeRecoDelaysWithRadioSpline(const int nAnt, const float zCenter, const vector<vector<double> >& antLoc,
                                     const float radius, const int nSideExp,
                                     float *recoDelays, float *recoDelays_V, float *recoDelays_H)
{

if(zCenter > 0 ) cerr<<"zCenter should be negative under the ice surfac\n";
/*
 * Initializing Healpix base
 */
if(nSideExp < 0 || nSideExp > 8){ cerr<<"Invalid nSideExp\n"; return -1; }
int nSide = pow(2, nSideExp);
Healpix_Base hpBase = Healpix_Base(nSide, HEALPIX_ORDERING, SET_NSIDE);
int nDir = hpBase.Npix();
pointing pt;
cout<<"Healpix base initialization success. nDir: "<<nDir<<endl;

/*
 * Get delay from each pixel direction
 */

//float *recoDelays;
//recoDelays  = (float*)malloc(nDir*nAnt*sizeof(float));
//recoDelays_V= (float*)malloc(nDir*(nAnt/2)*sizeof(float));
//recoDelays_H= (float*)malloc(nDir*(nAnt/2)*sizeof(float));
float test_zenith, test_azimuth;

/* Using radiospline to delay times */
char * tablePath = getenv("RADIOSPLINE_TABLE_DIR");
if (tablePath == NULL) {
    std::cout << "ERROR: please point the RADIOSPLINE_TABLE_DIR environment variable to" << std::endl;
    std::cout << " the spline .fits table directory." << std::endl;
    return -1;
}
std::string tablePathStr(tablePath);
RayDelay ray(tablePathStr+"/"+ICE_FILE,
             tablePathStr+"/"+AIR_FILE,
             tablePathStr+"/"+SHADOW_FILE);
cout<<"RayDelay object created\n";

double r, zRec, zSrc;
float tempDelay, meanDelay=0.f;
vector<float> solvedDelay;
double coordSrc[3], coordTrg[3];
//cout<<"recoDelays:\n";

   for(int pix=0; pix<nDir; pix++){

      pt = hpBase.pix2ang( pix );
      test_zenith  = pt.theta; //  in radians
      test_azimuth = pt.phi  ;

      coordSrc[0] = radius*sin(test_zenith)*cos(test_azimuth);
      coordSrc[1] = radius*sin(test_zenith)*sin(test_azimuth);
      coordSrc[2] = radius*cos(test_zenith);

      //cout<<"coordSrc: "<<coordSrc[0]<<"\t"<<coordSrc[1]<<"\t"<<coordSrc[2]<<endl;
      //cout<<"nAnt: "<<nAnt<<endl;
      //cout<<"tempDelay:\n";
      for(int k=0; k<nAnt; k++){
      //cout<<"k: "<<k<<endl;
      coordTrg[0] = (antLoc[k][0]);
      coordTrg[1] = (antLoc[k][1]);
      coordTrg[2] = (antLoc[k][2]);
      if (Detector2Cylinder(coordSrc, coordTrg, zCenter, &r, &zRec, &zSrc) != 0)
      std::cout << "ERROR: couldn't convert to cylindrical coordinates." << std::endl;

      tempDelay = static_cast<float>(ray.GetPropagationTime(r, zRec, zSrc));
      //cout<<tempDelay<<" ";
      if( tempDelay > 1.f )
         //if( k<8 || k>11 )
            solvedDelay.push_back(tempDelay);

      recoDelays[pix*nAnt + k] = tempDelay;

      }
      //cout<<endl;
      meanDelay = getMeanDelay( solvedDelay );
      //cout<<"meanDelay = "<<meanDelay<<endl;

      for(int k=0; k<nAnt; k++){

      if(recoDelays[pix*nAnt + k] > 1.f ){
         recoDelays[pix*nAnt + k] -= meanDelay;
         //cout<<recoDelays[pix*nAnt + k]<<" ";
         if(k<8) recoDelays_V[pix*nAnt/2 + k]   = recoDelays[pix*nAnt + k];
         else    recoDelays_H[pix*nAnt/2 + k-8] = recoDelays[pix*nAnt + k];
      } else {
          recoDelays[pix*nAnt + k] = -1e10; // no solution!
          if(k<8) recoDelays_V[pix*nAnt/2 + k]   = recoDelays[pix*nAnt + k];
          else    recoDelays_H[pix*nAnt/2 + k-8] = recoDelays[pix*nAnt + k];
          //cout<<"recoDelays: "<<recoDelays[pix*nAnt + k]<<"\t";
      }
      //cout<<"End of assigning delays\n";
      }//end of nAnt
      //cout<<endl;
   }//end of pix

/*
 * Write FITS file
 */
/*
cout<<"Creating Healpix map of channel 0 & 3 delays, and writing to FITS....\n";

float *delays = (float*)calloc(nDir, sizeof(float));

for(int pix=0; pix<nDir; pix++){

   delays[pix] = recoDelays[pix*nAnt + 3] - recoDelays[pix*nAnt + 0];

}
arr<float> delaysArr = arr<float>(&delays[0], (size_t)nDir);
Healpix_Map<float> delaysSkyMap = Healpix_Map<float>(delaysArr, RING);

fitshandle fitsOut;
remove("delaysSkyMap_chan0_3.fits");
fitsOut.create("delaysSkyMap_chan0_3.fits");

write_Healpix_map_to_fits(fitsOut, delaysSkyMap, PLANCK_FLOAT32);
cout<<"Healpix map written\n";

free(delays);
*/
   return 0;
}

int compute3DRecoDelaysWithRadioSpline(const int nAnt, const float zCenter, const vector<vector<double> >& antLoc,
                                     //const float radius, const int nSideExp,
                                     Healpix_Onion *onion,
                                     float *recoDelays, float *recoDelays_V, float *recoDelays_H)
{

if(zCenter > 0 ) cerr<<"zCenter should be negative under the ice surfac\n";
/*
 * Initializing Healpix base
 */
/*
if(nSideExp < 0 || nSideExp > 8){ cerr<<"Invalid nSideExp\n"; return -1; }
int nSide = pow(2, nSideExp);
Healpix_Base hpBase = Healpix_Base(nSide, RING, SET_NSIDE);
int nDir = hpBase.Npix();
pointing pt;
cout<<"Healpix base initialization success. nDir: "<<nDir<<endl;
*/

/*
 * Initialize variables related to Healpix_Onion
 */

int nDir   = onion->nDir;
int nLayer = onion->nLayer;
cout<<"Healpix_Onion info obtained. nDir: "<<nDir<<" nLayer: "<<nLayer<<endl;

/*
 * Get delay from each layer and pixel direction
 */

//float *recoDelays;
//recoDelays  = (float*)malloc(nDir*nAnt*sizeof(float));
//recoDelays_V= (float*)malloc(nDir*(nAnt/2)*sizeof(float));
//recoDelays_H= (float*)malloc(nDir*(nAnt/2)*sizeof(float));
float test_r, test_zenith, test_azimuth;

/* Using radiospline to delay times */
char * tablePath = getenv("RADIOSPLINE_TABLE_DIR");
if (tablePath == NULL) {
    std::cout << "ERROR: please point the RADIOSPLINE_TABLE_DIR environment variable to" << std::endl;
    std::cout << " the spline .fits table directory." << std::endl;
    return -1;
}
std::string tablePathStr(tablePath);
RayDelay ray(tablePathStr+"/"+ICE_FILE,
             tablePathStr+"/"+AIR_FILE,
             tablePathStr+"/"+SHADOW_FILE);
cout<<"RayDelay object created\n";

double r, zRec, zSrc;
float tempDelay, meanDelay=0.f;
vector<float> solvedDelay;
double coordSrc[3], coordTrg[3];
//cout<<"recoDelays:\n";
for(int layer=0; layer<nLayer; layer++){

   test_r = onion->layerRadii[layer]; //in meters
   //cout<<"test_r at layer "<<layer<<" is: "<<test_r<<endl;
   for(int pix=0; pix<nDir; pix++){

      //pt = hpBase.pix2ang( pix );
      //test_zenith  = pt.theta; //  in radians
      //test_azimuth = pt.phi  ;
      test_zenith  = onion->getPointing( pix ).theta;  //in radians
      test_azimuth = onion->getPointing( pix ).phi  ;
      /*
      coordSrc[0] = radius*sin(test_zenith)*cos(test_azimuth);
      coordSrc[1] = radius*sin(test_zenith)*sin(test_azimuth);
      coordSrc[2] = radius*cos(test_zenith);
      */
      coordSrc[0] = test_r*sin(test_zenith)*cos(test_azimuth);
      coordSrc[1] = test_r*sin(test_zenith)*sin(test_azimuth);
      coordSrc[2] = test_r*cos(test_zenith);


      //cout<<"coordSrc: "<<coordSrc[0]<<"\t"<<coordSrc[1]<<"\t"<<coordSrc[2]<<endl;
      //cout<<"nAnt: "<<nAnt<<endl;
      //cout<<"tempDelay:\n";
      for(int k=0; k<nAnt; k++){
      //cout<<"k: "<<k<<endl;
      coordTrg[0] = (antLoc[k][0]);
      coordTrg[1] = (antLoc[k][1]);
      coordTrg[2] = (antLoc[k][2]);
      if (Detector2Cylinder(coordSrc, coordTrg, zCenter, &r, &zRec, &zSrc) != 0)
      std::cout << "ERROR: couldn't convert to cylindrical coordinates." << std::endl;

      tempDelay = static_cast<float>(ray.GetPropagationTime(r, zRec, zSrc));
      //cout<<tempDelay<<" ";
      if( tempDelay > 1.f )
         //if( k<8 || k>11 )
            solvedDelay.push_back(tempDelay);
      if( tempDelay > 1e9 ) cout<<"Unreasonbaly large delay\n";

      recoDelays[layer*nDir*nAnt + pix*nAnt + k] = tempDelay;

      }//end of k
      //cout<<endl;
      meanDelay = getMeanDelay( solvedDelay );
      //cout<<"meanDelay = "<<meanDelay<<endl;

      for(int k=0; k<nAnt; k++){

      if(recoDelays[layer*nDir*nAnt + pix*nAnt + k] > 1.f ){
         recoDelays[layer*nDir*nAnt + pix*nAnt + k] -= meanDelay;
         //cout<<recoDelays[layer*nDir*nAnt + pix*nAnt + k]<<" ";
         if(k<8) recoDelays_V[layer*nDir*nAnt/2 + pix*nAnt/2 + k]   = recoDelays[layer*nDir*nAnt + pix*nAnt + k];
         else    recoDelays_H[layer*nDir*nAnt/2 + pix*nAnt/2 + k-8] = recoDelays[layer*nDir*nAnt + pix*nAnt + k];
      } else {
          recoDelays[layer*nDir*nAnt + pix*nAnt + k] = -1e10; // no solution!
          if(k<8) recoDelays_V[layer*nDir*nAnt/2 + pix*nAnt/2 + k]   = recoDelays[layer*nDir*nAnt + pix*nAnt + k];
          else    recoDelays_H[layer*nDir*nAnt/2 + pix*nAnt/2 + k-8] = recoDelays[layer*nDir*nAnt + pix*nAnt + k];
          //cout<<"recoDelays: "<<recoDelays[layer*nDir*nAnt + pix*nAnt + k]<<"\t";
      }
      //cout<<"End of assigning delays\n";
      }//end of nAnt
      //cout<<endl;
   }//end of pix
}//end of layer
/*
 * Write FITS file
 */
/*
cout<<"Creating Healpix map of channel 0 & 4 delays, and writing to FITS....\n";

float *delays = (float*)calloc(nDir, sizeof(float));
arr<float> delaysArr;
Healpix_Map<float> delaysSkyMap;
fitshandle fitsOut;
char filename[200];

for(int layer=0; layer<nLayer; layer++){
   for(int pix=0; pix<nDir; pix++){

   delays[pix] = recoDelays[layer*nDir*nAnt + pix*nAnt + 0] - recoDelays[layer*nDir*nAnt + pix*nAnt + 4];
   if(delays[pix] > 1e9 || delays[pix] < -1e9 ) delays[pix] = 0.f;

   }
   cout<<endl;
delaysArr = arr<float>(&delays[0], (size_t)nDir);
delaysSkyMap = Healpix_Map<float>(delaysArr, RING);

//fitshandle fitsOut;
sprintf(filename, "delaysSkyMap_A3_layer%d_chan0_4.fits", layer);
remove(filename);
fitsOut.create(filename);

write_Healpix_map_to_fits(fitsOut, delaysSkyMap, PLANCK_FLOAT32);
cout<<"Healpix map written\n";

}

free(delays);
*/
   return 0;
}

int compute3DZoomedRecoDelaysWithRadioSpline(const int nAnt, const float zCenter, const vector<vector<double> >& antLoc,
                                     //const float radius, const int nSideExp,
                                     Healpix_Onion *onion,
                                     float *recoDelays, float *recoDelays_V, float *recoDelays_H,
                                     const int last2DMaxPixIdx, const int nPix_nested)
{

if(zCenter > 0 ) cerr<<"zCenter should be negative under the ice surfac\n";
/*
 * Initializing Healpix base
 */
/*
if(nSideExp < 0 || nSideExp > 8){ cerr<<"Invalid nSideExp\n"; return -1; }
int nSide = pow(2, nSideExp);
Healpix_Base hpBase = Healpix_Base(nSide, RING, SET_NSIDE);
int nDir = hpBase.Npix();
pointing pt;
cout<<"Healpix base initialization success. nDir: "<<nDir<<endl;
*/

/*
 * Initialize variables related to Healpix_Onion
 */

int nDir   = onion->nDir;
int nLayer = onion->nLayer;
cout<<"Healpix_Onion info obtained. nDir: "<<nDir<<" nLayer: "<<nLayer<<endl;

/*
 * Get delay from each layer and pixel direction
 */

//float *recoDelays;
//recoDelays  = (float*)malloc(nDir*nAnt*sizeof(float));
//recoDelays_V= (float*)malloc(nDir*(nAnt/2)*sizeof(float));
//recoDelays_H= (float*)malloc(nDir*(nAnt/2)*sizeof(float));
float test_r, test_zenith, test_azimuth;

/* Using radiospline to delay times */
char * tablePath = getenv("RADIOSPLINE_TABLE_DIR");
if (tablePath == NULL) {
    std::cout << "ERROR: please point the RADIOSPLINE_TABLE_DIR environment variable to" << std::endl;
    std::cout << " the spline .fits table directory." << std::endl;
    return -1;
}
std::string tablePathStr(tablePath);
RayDelay ray(tablePathStr+"/"+ICE_FILE,
             tablePathStr+"/"+AIR_FILE,
             tablePathStr+"/"+SHADOW_FILE);
cout<<"RayDelay object created\n";

double r, zRec, zSrc;
float tempDelay, meanDelay=0.f;
vector<float> solvedDelay;
double coordSrc[3], coordTrg[3];
//cout<<"recoDelays:\n";
int zoomedPix;

for(int layer=0; layer<nLayer; layer++){

   test_r = onion->layerRadii[layer]; //in meters
   //cout<<"test_r at layer "<<layer<<" is: "<<test_r<<endl;
   //for(int pix=0; pix<nDir; pix++){
   for(int pix=0; pix<nPix_nested; pix++){

      zoomedPix = last2DMaxPixIdx * nPix_nested + pix;
      //pt = hpBase.pix2ang( pix );
      //test_zenith  = pt.theta; //  in radians
      //test_azimuth = pt.phi  ;
      test_zenith  = onion->getPointing( zoomedPix + nDir*layer ).theta;  //in radians
      test_azimuth = onion->getPointing( zoomedPix + nDir*layer ).phi  ;
      /*
      coordSrc[0] = radius*sin(test_zenith)*cos(test_azimuth);
      coordSrc[1] = radius*sin(test_zenith)*sin(test_azimuth);
      coordSrc[2] = radius*cos(test_zenith);
      */
      coordSrc[0] = test_r*sin(test_zenith)*cos(test_azimuth);
      coordSrc[1] = test_r*sin(test_zenith)*sin(test_azimuth);
      coordSrc[2] = test_r*cos(test_zenith);


      //cout<<"coordSrc: "<<coordSrc[0]<<"\t"<<coordSrc[1]<<"\t"<<coordSrc[2]<<endl;
      //cout<<"nAnt: "<<nAnt<<endl;
      //cout<<"tempDelay:\n";
      for(int k=0; k<nAnt; k++){
      //cout<<"k: "<<k<<endl;
      coordTrg[0] = (antLoc[k][0]);
      coordTrg[1] = (antLoc[k][1]);
      coordTrg[2] = (antLoc[k][2]);
      if (Detector2Cylinder(coordSrc, coordTrg, zCenter, &r, &zRec, &zSrc) != 0)
      std::cout << "ERROR: couldn't convert to cylindrical coordinates." << std::endl;

      tempDelay = static_cast<float>(ray.GetPropagationTime(r, zRec, zSrc));
      //cout<<tempDelay<<" ";
      if( tempDelay > 1.f )
         //if( k<8 || k>11 )
            solvedDelay.push_back(tempDelay);
      if( tempDelay > 1e9 ) cout<<"Unreasonbaly large delay\n";

      //recoDelays[layer*nDir*nAnt + pix*nAnt + k] = tempDelay;
      recoDelays[layer*nPix_nested*nAnt + pix*nAnt + k] = tempDelay;

      }//end of k
      //cout<<endl;
      meanDelay = getMeanDelay( solvedDelay );
      //cout<<"meanDelay = "<<meanDelay<<endl;

      for(int k=0; k<nAnt; k++){

      if(recoDelays[layer*nPix_nested*nAnt + pix*nAnt + k] > 1.f ){
         recoDelays[layer*nPix_nested*nAnt + pix*nAnt + k] -= meanDelay;
         //cout<<recoDelays[layer*nDir*nAnt + pix*nAnt + k]<<" ";
         if(k<8) recoDelays_V[layer*nPix_nested*nAnt/2 + pix*nAnt/2 + k]   = recoDelays[layer*nPix_nested*nAnt + pix*nAnt + k];
         else    recoDelays_H[layer*nPix_nested*nAnt/2 + pix*nAnt/2 + k-8] = recoDelays[layer*nPix_nested*nAnt + pix*nAnt + k];
      } else {
          recoDelays[layer*nPix_nested*nAnt + pix*nAnt + k] = -1e10; // no solution!
          if(k<8) recoDelays_V[layer*nPix_nested*nAnt/2 + pix*nAnt/2 + k]   = recoDelays[layer*nPix_nested*nAnt + pix*nAnt + k];
          else    recoDelays_H[layer*nPix_nested*nAnt/2 + pix*nAnt/2 + k-8] = recoDelays[layer*nPix_nested*nAnt + pix*nAnt + k];
          //cout<<"recoDelays: "<<recoDelays[layer*nDir*nAnt + pix*nAnt + k]<<"\t";
      }
      //cout<<"End of assigning delays\n";
      }//end of nAnt
      //cout<<endl;
   }//end of pix
}//end of layer
/*
 * Write FITS file
 */
/*
cout<<"Creating Healpix map of channel 0 & 4 delays, and writing to FITS....\n";

float *delays = (float*)calloc(nDir, sizeof(float));
arr<float> delaysArr;
Healpix_Map<float> delaysSkyMap;
fitshandle fitsOut;
char filename[200];

for(int layer=0; layer<nLayer; layer++){
   for(int pix=0; pix<nDir; pix++){

   delays[pix] = recoDelays[layer*nDir*nAnt + pix*nAnt + 0] - recoDelays[layer*nDir*nAnt + pix*nAnt + 4];
   if(delays[pix] > 1e9 || delays[pix] < -1e9 ) delays[pix] = 0.f;

   }
   cout<<endl;
delaysArr = arr<float>(&delays[0], (size_t)nDir);
delaysSkyMap = Healpix_Map<float>(delaysArr, RING);

//fitshandle fitsOut;
sprintf(filename, "delaysSkyMap_A3_layer%d_chan0_4.fits", layer);
remove(filename);
fitsOut.create(filename);

write_Healpix_map_to_fits(fitsOut, delaysSkyMap, PLANCK_FLOAT32);
cout<<"Healpix map written\n";

}

free(delays);
*/
   return 0;
}

int computeRecoDelaysWithConstantN(const int nAnt, const float zCenter, const vector<vector<double> >& antLoc,
                                     //const float radius, const int nSideExp,
                                     Healpix_Onion *onion,
                                     float *recoDelays, float *recoDelays_V, float *recoDelays_H)
{

if(zCenter > 0 ) cerr<<"zCenter should be negative under the ice surfac\n";
/*
 * Initializing Healpix base
 */
//if(nSideExp < 0 || nSideExp > 7){ cerr<<"Invalid nSideExp\n"; return -1; }
//int nSide = pow(2, nSideExp);
//Healpix_Base hpBase = Healpix_Base(nSide, RING, SET_NSIDE);
//int nDir = hpBase.Npix();
//pointing pt;
//cout<<"Healpix base initialization success. nDir: "<<nDir<<endl;

/*
 * Initialize variables related to Healpix_Onion
 */

int nDir   = onion->nDir;
int nLayer = onion->nLayer;
cout<<"Healpix_Onion info obtained. nDir: "<<nDir<<" nLayer: "<<nLayer<<endl;


/*
 * Get delay from each pixel direction
 */

//float *recoDelays;
//recoDelays  = (float*)malloc(nDir*nAnt*sizeof(float));
//recoDelays_V= (float*)malloc(nDir*(nAnt/2)*sizeof(float));
//recoDelays_H= (float*)malloc(nDir*(nAnt/2)*sizeof(float));
float test_r, test_zenith, test_azimuth;

float tempDelay, meanDelay=0.f;
vector<float> solvedDelay;
double coordSrc[3], coordTrg[3];
//cout<<"recoDelays:\n";
for(int layer=0; layer<nLayer; layer++){

   test_r = onion->layerRadii[layer];

   for(int pix=0; pix<nDir; pix++){

      //pt = hpBase.pix2ang( pix );
      //test_zenith  = pt.theta; //  in radians
      //test_azimuth = pt.phi  ;

      test_zenith  = onion->getPointing( pix ).theta; // in radians
      test_azimuth = onion->getPointing( pix ).phi  ;

      coordSrc[0] = test_r*sin(test_zenith)*cos(test_azimuth);
      coordSrc[1] = test_r*sin(test_zenith)*sin(test_azimuth);
      coordSrc[2] = test_r*cos(test_zenith);

      if(coordSrc[2] > fabs(zCenter)){

      cerr<<"Warning!! Source above ice surface. The different index of refraction in air is _NOT_ implemented!"<<endl;
      for(int k=0; k<nAnt; k++){
        recoDelays[layer*nDir*nAnt + pix*nAnt + k] = -1e10;;
        if(k<8) recoDelays_V[layer*nDir*nAnt/2 + pix*nAnt/2 + k]   = recoDelays[layer*nDir*nAnt + pix*nAnt + k];
        else    recoDelays_H[layer*nDir*nAnt/2 + pix*nAnt/2 + k-8] = recoDelays[layer*nDir*nAnt + pix*nAnt + k];
      }
      } else if ( (zCenter + coordSrc[2]) <= -3000){

      cerr<<"Warning!! Source is in bedrock"<<endl;
      for(int k=0; k<nAnt; k++){
        recoDelays[layer*nDir*nAnt + pix*nAnt + k] = -1e10;;
        if(k<8) recoDelays_V[layer*nDir*nAnt/2 + pix*nAnt/2 + k]   = recoDelays[layer*nDir*nAnt + pix*nAnt + k];
        else    recoDelays_H[layer*nDir*nAnt/2 + pix*nAnt/2 + k-8] = recoDelays[layer*nDir*nAnt + pix*nAnt + k];
      }
      } else {
      //cout<<"coordSrc: "<<coordSrc[0]<<"\t"<<coordSrc[1]<<"\t"<<coordSrc[2]<<endl;

      for(int k=0; k<nAnt; k++){

      coordTrg[0] = (antLoc[k][0]);
      coordTrg[1] = (antLoc[k][1]);
      coordTrg[2] = (antLoc[k][2]);


      /* Assume simple spherical wave propagation with homogeneous isotropic ice */
      tempDelay = nIce * sqrt( (coordTrg[0] - coordSrc[0])*(coordTrg[0] - coordSrc[0])
                             + (coordTrg[1] - coordSrc[1])*(coordTrg[1] - coordSrc[1])
                             + (coordTrg[2] - coordSrc[2])*(coordTrg[2] - coordSrc[2])
                             ) / speedOfLight;
      solvedDelay.push_back(tempDelay);
      //cout<<tempDelay<<" ";
      recoDelays[layer*nDir*nAnt + pix*nAnt + k] = tempDelay;

      }
      //cout<<endl;
      meanDelay = getMeanDelay( solvedDelay );
      //cout<<"meanDelay = "<<meanDelay<<endl;

      for(int k=0; k<nAnt; k++){
         recoDelays[layer*nDir*nAnt + pix*nAnt + k] -= meanDelay;
         //cout<<recoDelays[pix*nAnt + k]<<" ";
         if(k<8) recoDelays_V[layer*nDir*nAnt/2 + pix*nAnt/2 + k]   = recoDelays[layer*nDir*nAnt + pix*nAnt + k];
         else    recoDelays_H[layer*nDir*nAnt/2 + pix*nAnt/2 + k-8] = recoDelays[layer*nDir*nAnt + pix*nAnt + k];
      //cout<<"End of assigning delays\n";
      }//end of nAnt
      //cout<<endl;
      }//end of else
   }//end of pix
}//end of layer

   return 0;
}

int computeZoomedRecoDelaysWithConstantN(const int nAnt, const float zCenter, const vector<vector<double> >& antLoc,
                                     //const float radius, const int nSideExp,
                                     Healpix_Onion *onion,
                                     float *recoDelays, float *recoDelays_V, float *recoDelays_H,
                                     const int last2DMaxPixIdx, const int nPix_nested)
{

if(zCenter > 0 ) cerr<<"zCenter should be negative under the ice surfac\n";
/*
 * Initializing Healpix base
 */
//if(nSideExp < 0 || nSideExp > 7){ cerr<<"Invalid nSideExp\n"; return -1; }
//int nSide = pow(2, nSideExp);
//Healpix_Base hpBase = Healpix_Base(nSide, RING, SET_NSIDE);
//int nDir = hpBase.Npix();
//pointing pt;
//cout<<"Healpix base initialization success. nDir: "<<nDir<<endl;

/*
 * Initialize variables related to Healpix_Onion
 */

int nDir   = onion->nDir;
int nLayer = onion->nLayer;
cout<<"Healpix_Onion info obtained. nDir: "<<nDir<<" nLayer: "<<nLayer<<endl;


/*
 * Get delay from each pixel direction
 */

//float *recoDelays;
//recoDelays  = (float*)malloc(nDir*nAnt*sizeof(float));
//recoDelays_V= (float*)malloc(nDir*(nAnt/2)*sizeof(float));
//recoDelays_H= (float*)malloc(nDir*(nAnt/2)*sizeof(float));
float test_r, test_zenith, test_azimuth;

float tempDelay, meanDelay=0.f;
vector<float> solvedDelay;
double coordSrc[3], coordTrg[3];
//cout<<"recoDelays:\n";
int zoomedPix;

for(int layer=0; layer<nLayer; layer++){

   test_r = onion->layerRadii[layer];

   //for(int pix=0; pix<nDir; pix++){
   for(int pix=0; pix<nPix_nested; pix++){

      zoomedPix = last2DMaxPixIdx * nPix_nested + pix;
      //pt = hpBase.pix2ang( pix );
      //test_zenith  = pt.theta; //  in radians
      //test_azimuth = pt.phi  ;

      test_zenith  = onion->getPointing( zoomedPix + nDir*layer ).theta; // in radians
      test_azimuth = onion->getPointing( zoomedPix + nDir*layer ).phi  ;

      coordSrc[0] = test_r*sin(test_zenith)*cos(test_azimuth);
      coordSrc[1] = test_r*sin(test_zenith)*sin(test_azimuth);
      coordSrc[2] = test_r*cos(test_zenith);

      if(coordSrc[2] > fabs(zCenter)){

      cerr<<"Warning!! Source above ice surface. The different index of refraction in air is _NOT_ implemented!"<<endl;
      for(int k=0; k<nAnt; k++){
        recoDelays[layer*nPix_nested*nAnt + pix*nAnt + k] = -1e10;;
        if(k<8) recoDelays_V[layer*nPix_nested*nAnt/2 + pix*nAnt/2 + k]   = recoDelays[layer*nPix_nested*nAnt + pix*nAnt + k];
        else    recoDelays_H[layer*nPix_nested*nAnt/2 + pix*nAnt/2 + k-8] = recoDelays[layer*nPix_nested*nAnt + pix*nAnt + k];
      }
      } else if ( (zCenter + coordSrc[2]) <= -3000){

      cerr<<"Warning!! Source is in bedrock"<<endl;
      for(int k=0; k<nAnt; k++){
        recoDelays[layer*nPix_nested*nAnt + pix*nAnt + k] = -1e10;;
        if(k<8) recoDelays_V[layer*nPix_nested*nAnt/2 + pix*nAnt/2 + k]   = recoDelays[layer*nPix_nested*nAnt + pix*nAnt + k];
        else    recoDelays_H[layer*nPix_nested*nAnt/2 + pix*nAnt/2 + k-8] = recoDelays[layer*nPix_nested*nAnt + pix*nAnt + k];
      }
      } else {
      //cout<<"coordSrc: "<<coordSrc[0]<<"\t"<<coordSrc[1]<<"\t"<<coordSrc[2]<<endl;

      for(int k=0; k<nAnt; k++){

      coordTrg[0] = (antLoc[k][0]);
      coordTrg[1] = (antLoc[k][1]);
      coordTrg[2] = (antLoc[k][2]);


      /* Assume simple spherical wave propagation with homogeneous isotropic ice */
      tempDelay = nIce * sqrt( (coordTrg[0] - coordSrc[0])*(coordTrg[0] - coordSrc[0])
                             + (coordTrg[1] - coordSrc[1])*(coordTrg[1] - coordSrc[1])
                             + (coordTrg[2] - coordSrc[2])*(coordTrg[2] - coordSrc[2])
                             ) / speedOfLight;
      solvedDelay.push_back(tempDelay);
      //cout<<tempDelay<<" ";
      recoDelays[layer*nPix_nested*nAnt + pix*nAnt + k] = tempDelay;

      }
      //cout<<endl;
      meanDelay = getMeanDelay( solvedDelay );
      //cout<<"meanDelay = "<<meanDelay<<endl;

      for(int k=0; k<nAnt; k++){
         recoDelays[layer*nPix_nested*nAnt + pix*nAnt + k] -= meanDelay;
         //cout<<recoDelays[pix*nAnt + k]<<" ";
         if(k<8) recoDelays_V[layer*nPix_nested*nAnt/2 + pix*nAnt/2 + k]   = recoDelays[layer*nPix_nested*nAnt + pix*nAnt + k];
         else    recoDelays_H[layer*nPix_nested*nAnt/2 + pix*nAnt/2 + k-8] = recoDelays[layer*nPix_nested*nAnt + pix*nAnt + k];
      //cout<<"End of assigning delays\n";
      }//end of nAnt
      //cout<<endl;
      }//end of else
   }//end of pix
}//end of layer

   return 0;
}

int getMaxBin(TGraph *gr){

   double t, v, max;
   max = 0.;
   int maxBin = 0;
   //cout<<"gr->GetN(): "<<gr->GetN()<<endl;

   for(int s=0; s<gr->GetN(); s++){

   gr->GetPoint(s, t, v);
   if( fabs(v) > max ){
      max = fabs(v);
      maxBin = s;
   }
   }
   return maxBin;
}

void setMeanAndSigmaInNoMax(TGraph *gr, double *stats){

   int bin = gr->GetN();
   int MaxBin = getMaxBin( gr );
   //cout<<"MaxBin: "<<MaxBin<<endl;
   int binCounter=0;

   double mean =0;
   double sigma=0;
   double t, v;

   if( MaxBin <= bin/4 ){

      for (int i=MaxBin+bin/4; i<bin; i++){
      gr->GetPoint(i, t, v);
      mean  += v;
      sigma += v * v;
      binCounter++;
      }
   }

   else if( MaxBin >= 3*bin/4 ){

      for (int i=0; i<MaxBin-bin/4; i++){
      gr->GetPoint(i, t, v);
      mean  += v;
      sigma += v * v;
      binCounter++;
      }
   }

   else{

      for (int i=0; i<MaxBin-bin/4; i++){
      gr->GetPoint(i, t, v);
      mean  += v;
      sigma += v * v;
      binCounter++;
      }

      for (int i=MaxBin+bin/4; i<bin; i++){
      gr->GetPoint(i, t, v);
      mean  += v;
      sigma += v * v;
      binCounter++;
      }
   }

   mean  = mean / (double)binCounter;
   sigma = TMath::Sqrt( ( sigma - ((double)binCounter * mean * mean )) / (double)(binCounter - 1) );
   //cout<<"mean="<<mean<<"\tsigma="<<sigma<<endl;
   //delete gr;

   stats[0] = mean;
   stats[1] = sigma;

}
void getNchnl(const vector<TGraph *>& cleanEvent, double threshold, int *nchnlArray){

   double sigma;
   double mean;
   double statsArray[2]={0};
   //TGraph *gr;
   //Double_t *volts;
   double t, v;
   int bin;
   int passThreshold[16] = {0}; //IF THE CHANNEL PASSED THRESHOLD*SIGMA, THIS NUMBER IS SET TO ONE
   int totalPassedChnl, totalPassedVpol, totalPassedHpol;
   totalPassedChnl = totalPassedVpol = totalPassedHpol = 0;
   //AraGeomTool *tempGeom = AraGeomTool::Instance();
   //Int_t stationId = theEvent->stationId;

   for ( int ch=0; ch<(int)cleanEvent.size(); ch++){

         //AraAntPol::AraAntPol_t polType = tempGeom->getStationInfo(stationId)->getAntennaInfo(ch)->polType;

         //if ( polType == AraAntPol::kVertical ){
         //          //cout<<"mean="<<mean<<"\tsigma="<<sigma<<endl;
         //gr = theEvent->getGraphFromRFChan(ch);
         //gr = cleanEvent[ch];
         //volts = gr->GetY();
         bin   = cleanEvent[ch]->GetN();
         //setMeanAndSigmaInNoMax(gr,statsArray);
         setMeanAndSigmaInNoMax(cleanEvent[ch], statsArray);

         mean  = statsArray[0];
         sigma = statsArray[1];

         for (int binCounter=0; binCounter<bin; binCounter++){

            cleanEvent[ch]->GetPoint(binCounter, t, v);
            if ( fabs(v - mean) > threshold * sigma ){

               totalPassedChnl += 1;

               if( ch<8 ) totalPassedVpol += 1;
               else       totalPassedHpol += 1;
               //if ( polType == AraAntPol::kVertical){ totalPassedVpol += 1;
               //} else if ( polType == AraAntPol::kHorizontal){ totalPassedHpol += 1;
               //} else { cerr<<"********************* polType not vpol or hpol !!! ***********************"<<endl;
               //}
               break;
            }
         }

      //delete gr;


   }//end of ch

nchnlArray[0] = totalPassedChnl;
nchnlArray[1] = totalPassedVpol;
nchnlArray[2] = totalPassedHpol;

}

void getNchnlMask(const vector<TGraph *>& cleanEvent, double threshold, int *nchnlArray, const int *chanMask, int *goodChan){

   double sigma;
   double mean;
   double statsArray[2]={0};
   //TGraph *gr;
   //Double_t *volts;
   double t, v;
   int bin;
   int passThreshold[cleanEvent.size()];  //IF THE CHANNEL PASSED THRESHOLD*SIGMA, THIS NUMBER IS SET TO ONE
   int totalPassedChnl, totalPassedVpol, totalPassedHpol;
   totalPassedChnl = totalPassedVpol = totalPassedHpol = 0;
   //AraGeomTool *tempGeom = AraGeomTool::Instance();
   //Int_t stationId = theEvent->stationId;

   for ( int ch=0; ch<(int)cleanEvent.size(); ch++){

         passThreshold[ch] = 0;
         //AraAntPol::AraAntPol_t polType = tempGeom->getStationInfo(stationId)->getAntennaInfo(ch)->polType;

         //if ( polType == AraAntPol::kVertical ){
         //          //cout<<"mean="<<mean<<"\tsigma="<<sigma<<endl;
         //gr = theEvent->getGraphFromRFChan(ch);
         //gr = cleanEvent[ch];
         //volts = gr->GetY();
         bin   = cleanEvent[ch]->GetN();
         //setMeanAndSigmaInNoMax(gr,statsArray);
         setMeanAndSigmaInNoMax(cleanEvent[ch], statsArray);

         mean  = statsArray[0];
         sigma = statsArray[1];

         for (int binCounter=0; binCounter<bin; binCounter++){

            cleanEvent[ch]->GetPoint(binCounter, t, v);
            if ( fabs(v - mean) > threshold * sigma ){

               //totalPassedChnl += 1;
               passThreshold[ch] = 1;

               //if( ch<8 ) totalPassedVpol += 1;
               //else       totalPassedHpol += 1;
               //if ( polType == AraAntPol::kVertical){ totalPassedVpol += 1;
               //} else if ( polType == AraAntPol::kHorizontal){ totalPassedHpol += 1;
               //} else { cerr<<"********************* polType not vpol or hpol !!! ***********************"<<endl;
               //}
               break;
            }
         }

      //delete gr;


   }//end of ch

   for(int ch=0; ch<16; ch++){

   goodChan[ch] = (chanMask[ch] && passThreshold[ch] );

   totalPassedChnl += goodChan[ch];
   if(ch<8) totalPassedVpol += goodChan[ch];
   else     totalPassedHpol += goodChan[ch];

   }

nchnlArray[0] = totalPassedChnl;
nchnlArray[1] = totalPassedVpol;
nchnlArray[2] = totalPassedHpol;

}

/*
 * Version of getNchnlMask taking into account saturation. Saturated channels are not used
 */
void getNchnlMaskSat(const vector<TGraph *>& cleanEvent, double threshold, int *nchnlArray, const int *chanMask, int *goodChan, int& numSatChan){

   double sigma;
   double mean;
   double statsArray[2]={0};
   //TGraph *gr;
   //Double_t *volts;
   double t, v;
   int bin;
   int passThreshold[cleanEvent.size()];  //IF THE CHANNEL PASSED THRESHOLD*SIGMA, THIS NUMBER IS SET TO ONE
   int saturated[cleanEvent.size()];
   int totalPassedChnl, totalPassedVpol, totalPassedHpol;
   int totalSatChnl, totalSatVpol, totalSatHpol;
   totalPassedChnl = totalPassedVpol = totalPassedHpol = 0;
   totalSatChnl = totalSatVpol = totalSatHpol = 0;
   //AraGeomTool *tempGeom = AraGeomTool::Instance();
   //Int_t stationId = theEvent->stationId;

   for ( int ch=0; ch<(int)cleanEvent.size(); ch++){

         passThreshold[ch] = 0;
         saturated[ch]= 0;
         //AraAntPol::AraAntPol_t polType = tempGeom->getStationInfo(stationId)->getAntennaInfo(ch)->polType;

         //if ( polType == AraAntPol::kVertical ){
         //          //cout<<"mean="<<mean<<"\tsigma="<<sigma<<endl;
         //gr = theEvent->getGraphFromRFChan(ch);
         //gr = cleanEvent[ch];
         //volts = gr->GetY();
         bin   = cleanEvent[ch]->GetN();
         //setMeanAndSigmaInNoMax(gr,statsArray);
         setMeanAndSigmaInNoMax(cleanEvent[ch], statsArray);

         mean  = statsArray[0];
         sigma = statsArray[1];

         for (int binCounter=0; binCounter<bin; binCounter++){

            cleanEvent[ch]->GetPoint(binCounter, t, v);

            /* check if SNR > threshold*sigma */
            if ( fabs(v - mean) > threshold * sigma ){

               //totalPassedChnl += 1;
               passThreshold[ch] = 1;

               //if( ch<8 ) totalPassedVpol += 1;
               //else       totalPassedHpol += 1;
               //if ( polType == AraAntPol::kVertical){ totalPassedVpol += 1;
               //} else if ( polType == AraAntPol::kHorizontal){ totalPassedHpol += 1;
               //} else { cerr<<"********************* polType not vpol or hpol !!! ***********************"<<endl;
               //}
               //break;
            }

            /* check if saturated at +/- 1000mV */
            if( fabs( fabs(v) - 1000. ) < 0.5 ){

               saturated[ch] = 1;
            }


         }

      //delete gr;


   }//end of ch

   for(int ch=0; ch<16; ch++){

   goodChan[ch] = (chanMask[ch] && passThreshold[ch] );
   if( saturated[ch] ) goodChan[ch] = 0;

   totalPassedChnl += goodChan[ch];
   totalSatChnl    += saturated[ch];
   if(ch<8){ totalPassedVpol += goodChan[ch]; totalSatVpol += saturated[ch]; }
   else    { totalPassedHpol += goodChan[ch]; totalSatHpol += saturated[ch]; }

   }

/* only look at Vpols now */
numSatChan = totalSatVpol;

nchnlArray[0] = totalPassedChnl;
nchnlArray[1] = totalPassedVpol;
nchnlArray[2] = totalPassedHpol;

}

/*
 * Returns the SNR as the number of sigmas between the waveform peak and mean for all channels. SNR is defined to be always positive (absolute value)
 */
void getChannelSNR(const vector<TGraph *>& cleanEvent, float *snrArray){

   double sigma;
   double mean;
   double absPeak;
   double statsArray[2]={0};
   //TGraph *gr;
   //Double_t *volts;
   double t, v;
   int bin;
   //int passThreshold[cleanEvent.size()];  //IF THE CHANNEL PASSED THRESHOLD*SIGMA, THIS NUMBER IS SET TO ONE
   //int saturated[cleanEvent.size()];
   //int totalPassedChnl, totalPassedVpol, totalPassedHpol;
   //int totalSatChnl, totalSatVpol, totalSatHpol;
   //totalPassedChnl = totalPassedVpol = totalPassedHpol = 0;
   //totalSatChnl = totalSatVpol = totalSatHpol = 0;
   //AraGeomTool *tempGeom = AraGeomTool::Instance();
   //Int_t stationId = theEvent->stationId;

   for ( int ch=0; ch<(int)cleanEvent.size(); ch++){

         absPeak = 0.;
         //passThreshold[ch] = 0;
         //saturated[ch]= 0;
         //AraAntPol::AraAntPol_t polType = tempGeom->getStationInfo(stationId)->getAntennaInfo(ch)->polType;

         //if ( polType == AraAntPol::kVertical ){
         //          //cout<<"mean="<<mean<<"\tsigma="<<sigma<<endl;
         //gr = theEvent->getGraphFromRFChan(ch);
         //gr = cleanEvent[ch];
         //volts = gr->GetY();
         bin   = cleanEvent[ch]->GetN();
         //setMeanAndSigmaInNoMax(gr,statsArray);
         setMeanAndSigmaInNoMax(cleanEvent[ch], statsArray);

         mean  = statsArray[0];
         sigma = statsArray[1];

         for (int binCounter=0; binCounter<bin; binCounter++){

            cleanEvent[ch]->GetPoint(binCounter, t, v);

            /* check if SNR > threshold*sigma */
            //if ( fabs(v - mean) > threshold * sigma ){
            if( fabs(v-mean) > absPeak ){
               //totalPassedChnl += 1;
               //passThreshold[ch] = 1;
               absPeak = fabs(v-mean);
               //if( ch<8 ) totalPassedVpol += 1;
               //else       totalPassedHpol += 1;
               //if ( polType == AraAntPol::kVertical){ totalPassedVpol += 1;
               //} else if ( polType == AraAntPol::kHorizontal){ totalPassedHpol += 1;
               //} else { cerr<<"********************* polType not vpol or hpol !!! ***********************"<<endl;
               //}
               //break;
            }

            /* check if saturated at +/- 1000mV */
            //if( fabs( fabs(v) - 1000. ) < 0.5 ){

            //   saturated[ch] = 1;
            //}


         }//end of binCounter

      //delete gr;
      snrArray[ch] = static_cast<float>(absPeak / sigma);

   }//end of ch
   /*
   for(int ch=0; ch<16; ch++){

   goodChan[ch] = (chanMask[ch] && passThreshold[ch] );
   if( saturated[ch] ) goodChan[ch] = 0;

   totalPassedChnl += goodChan[ch];
   totalSatChnl    += saturated[ch];
   if(ch<8){ totalPassedVpol += goodChan[ch]; totalSatVpol += saturated[ch]; }
   else    { totalPassedHpol += goodChan[ch]; totalSatHpol += saturated[ch]; }

   }
   */
/* only look at Vpols now */
/*
numSatChan = totalSatVpol;

nchnlArray[0] = totalPassedChnl;
nchnlArray[1] = totalPassedVpol;
nchnlArray[2] = totalPassedHpol;
*/
}

/*
 * Returns the unmodified SNR as the number of sigmas between the waveform peak and mean for all channels. SNR is defined to be always positive (absolute value)
 * Unmodified: in the case of the Hann windowing function being used, this means in the center 1/2 of the waveform where the amplitude is not modified.
 */
void getChannelUnmodifiedSNR(const vector<TGraph *>& cleanEvent, float *snrArray){

   double sigma;
   double mean;
   double absPeak;
   double statsArray[2]={0};
   //TGraph *gr;
   //Double_t *volts;
   double t, v;
   int bin;
   //int passThreshold[cleanEvent.size()];  //IF THE CHANNEL PASSED THRESHOLD*SIGMA, THIS NUMBER IS SET TO ONE
   //int saturated[cleanEvent.size()];
   //int totalPassedChnl, totalPassedVpol, totalPassedHpol;
   //int totalSatChnl, totalSatVpol, totalSatHpol;
   //totalPassedChnl = totalPassedVpol = totalPassedHpol = 0;
   //totalSatChnl = totalSatVpol = totalSatHpol = 0;
   //AraGeomTool *tempGeom = AraGeomTool::Instance();
   //Int_t stationId = theEvent->stationId;
   int mod = 4;//for modified Hann window

   for ( int ch=0; ch<(int)cleanEvent.size(); ch++){

         absPeak = 0.;
         //passThreshold[ch] = 0;
         //saturated[ch]= 0;
         //AraAntPol::AraAntPol_t polType = tempGeom->getStationInfo(stationId)->getAntennaInfo(ch)->polType;

         //if ( polType == AraAntPol::kVertical ){
         //          //cout<<"mean="<<mean<<"\tsigma="<<sigma<<endl;
         //gr = theEvent->getGraphFromRFChan(ch);
         //gr = cleanEvent[ch];
         //volts = gr->GetY();
         bin   = cleanEvent[ch]->GetN();
         //setMeanAndSigmaInNoMax(gr,statsArray);
         setMeanAndSigmaInNoMax(cleanEvent[ch], statsArray);

         mean  = statsArray[0];
         sigma = statsArray[1];

         for (int binCounter=(bin/mod); binCounter<((bin*(mod-1))/mod); binCounter++){

            cleanEvent[ch]->GetPoint(binCounter, t, v);

            /* check if SNR > threshold*sigma */
            //if ( fabs(v - mean) > threshold * sigma ){
            if( fabs(v-mean) > absPeak ){
               //totalPassedChnl += 1;
               //passThreshold[ch] = 1;
               absPeak = fabs(v-mean);
               //if( ch<8 ) totalPassedVpol += 1;
               //else       totalPassedHpol += 1;
               //if ( polType == AraAntPol::kVertical){ totalPassedVpol += 1;
               //} else if ( polType == AraAntPol::kHorizontal){ totalPassedHpol += 1;
               //} else { cerr<<"********************* polType not vpol or hpol !!! ***********************"<<endl;
               //}
               //break;
            }

            /* check if saturated at +/- 1000mV */
            //if( fabs( fabs(v) - 1000. ) < 0.5 ){

            //   saturated[ch] = 1;
            //}


         }//end of binCounter

      //delete gr;
      snrArray[ch] = static_cast<float>(absPeak / sigma);

   }//end of ch
   /*
   for(int ch=0; ch<16; ch++){

   goodChan[ch] = (chanMask[ch] && passThreshold[ch] );
   if( saturated[ch] ) goodChan[ch] = 0;

   totalPassedChnl += goodChan[ch];
   totalSatChnl    += saturated[ch];
   if(ch<8){ totalPassedVpol += goodChan[ch]; totalSatVpol += saturated[ch]; }
   else    { totalPassedHpol += goodChan[ch]; totalSatHpol += saturated[ch]; }

   }
   */
/* only look at Vpols now */
/*
numSatChan = totalSatVpol;

nchnlArray[0] = totalPassedChnl;
nchnlArray[1] = totalPassedVpol;
nchnlArray[2] = totalPassedHpol;
*/
}


//int recordDiff(int nSideExp, int maxPixIdx, float maxPixValue, double weight, float zen_true, float azi_true, float r_true, int *usedChan, char *rootFilename){
/*
int recordDiff(int nSideExp, recoData *summary, char *rootFilename){


   if(nSideExp < 0 || nSideExp > 7){ cerr<<"Invalid nSideExp\n"; return -1; }
   int nSide = pow(2, nSideExp);
   Healpix_Base hpBase = Healpix_Base(nSide, RING, SET_NSIDE);
   int nDir = hpBase.Npix();
   pointing pt;

   pt = hpBase.pix2ang( summary->maxPixIdx );
   //float dZen = (pt.theta - zen_true) * 180.f / M_PI;
   //float dAzi = (pt.phi   - azi_true) * 180.f / M_PI;
   float dZen = pt.theta * 180.f / M_PI - summary->trueZen;
   float dAzi = pt.phi   * 180.f / M_PI - summary->trueAzi;


   ifstream rootfile(rootFilename);
   TH1F *dZenDist;
   TH1F *dAziDist;
   TH2F *recoTrueZenDist;
   TH2F *recoTrueAziDist;
   TTree *dataTree;

   //double w = weight;
   //double zen = zen_true * 180.f / M_PI;
   //double azi = azi_true * 180.f / M_PI;
   double recZen = pt.theta * 180.f / M_PI;
   double recAzi = pt.phi   * 180.f / M_PI;
   summary->setRecoDir(recZen, recAzi);


   recoData dummyData;
   dummyData.duplicate(summary);
   //dummyData.getData( weight
   //                 , zen_true * 180.f / M_PI, azi_true * 180.f / M_PI
   //                 , pt.theta * 180.f / M_PI, pt.phi   * 180.f / M_PI
   //                 , r_true, 0.f
   //                 , usedChan
   //                 , maxPixValue
   //                 );
*//*
   printf("old weight: %f new weight: %f\nold trueRadius: %f new trueRadius: %f\nold recoRadius: %f new recoRadius: %f\nold trueZen: %f new trueZen: %f\nold trueAzi: %f new trueAzi: %f\nold maxPixIdx: %d new maxPixIdx: %d\nold maxPixCoherence: %f new maxPixCoherence: %f\n"
           , summary->weight, dummyData.weight
           , summary->trueRadius, dummyData.trueRadius
           , summary->recoRadius, dummyData.recoRadius
           , summary->trueZen, dummyData.trueZen
           , summary->trueAzi, dummyData.trueAzi
           , summary->maxPixIdx, dummyData.maxPixIdx
           , summary->maxPixCoherence, dummyData.maxPixCoherence
          );
   for(int ch=0; ch<16; ch++) printf("recoChan %d: old %d new %d", ch, summary->recoChan[ch], dummyData.recoChan[ch]);
*/
 /*
   if( !rootfile ){
   TFile fp(rootFilename, "NEW");
   dZenDist = new TH1F("recoZenDiff", "recoZenDiff", 360, -180, 180);
   dAziDist = new TH1F("recoAziDiff", "recoAziDiff", 720, -360, 360);
   recoTrueZenDist = new TH2F("recoTrueZenDist", "recoTrueZenDist", 180, 0, 180, 180, 0, 180);
   recoTrueAziDist = new TH2F("recoTrueAziDist", "recoTrueAziDist", 360, 0, 360, 360, 0, 360);

   dataTree = new TTree("dataTree", "dataTree");
   dataTree->Branch("recoData", &dummyData, "weight/D:trueZen/F:trueAzi/F:recoZen/F:recoAzi/F:trueRadius/F:recoRadius/F:recoChan[16]/I:maxPixIdx/I:coherence/F");

   dZenDist->Write();
   dAziDist->Write();
   recoTrueZenDist->Write();
   recoTrueAziDist->Write();
   dataTree->Write();
   fp.Close();
   }
   //else{
   TFile fp_2(rootFilename, "UPDATE");
   dZenDist = (TH1F*)fp_2.Get("recoZenDiff");
   dAziDist = (TH1F*)fp_2.Get("recoAziDiff");
   recoTrueZenDist = (TH2F*)fp_2.Get("recoTrueZenDist");
   recoTrueAziDist = (TH2F*)fp_2.Get("recoTrueAziDist");
   dataTree = (TTree*)fp_2.Get("dataTree");
   //dataTree->ResetBranchAddresses();
   dataTree->SetBranchAddress("recoData",&dummyData);

   cout<<"trueRadius: "<<dummyData.trueRadius<<endl;
   dZenDist->Fill(dZen);
   dAziDist->Fill(dAzi);
   //recoTrueZenDist->Fill(zen, recZen);
   //recoTrueAziDist->Fill(azi, recAzi);
   recoTrueZenDist->Fill(dummyData.trueZen, recZen);
   recoTrueAziDist->Fill(dummyData.trueAzi, recAzi);
   dataTree->Fill();

   //TFile fp_3(rootFilename, "RECREATE");
   dZenDist->Write();
   dAziDist->Write();
   recoTrueZenDist->Write();
   recoTrueAziDist->Write();
   dataTree->Write();
   fp_2.Close();
   //fp_3.Close();

   return 0;
}
*/
//int recordDiffGetFlag(int nSideExp, int maxPixIdx, float maxPixValue, double weight, float zen_true, float azi_true, float r_true, int *usedChan, char *rootFilename){
int recordDiffGetFlag(int nSideExp, recoData *summary, char *rootFilename){

/*
 * flag signals whether the event loop should break and save the current skymap, depending on whether the flag condition is met
 * flag 0: do not break
 * flag 1: break and save skymap
 */

   int flag = 0;

   if(nSideExp < 0 || nSideExp > 7){ cerr<<"Invalid nSideExp\n"; return -1; }
   int nSide = pow(2, nSideExp);
   Healpix_Base hpBase = Healpix_Base(nSide, HEALPIX_ORDERING, SET_NSIDE);
   int nDir = hpBase.Npix();
   pointing pt;

   pt = hpBase.pix2ang( summary->maxPixIdx );
   //float dZen = (pt.theta - zen_true) * 180.f / M_PI;
   //float dAzi = (pt.phi   - azi_true) * 180.f / M_PI;
   float dZen = pt.theta * 180.f / M_PI - summary->trueZen;
   float dAzi = pt.phi   * 180.f / M_PI - summary->trueAzi;

   ifstream rootfile(rootFilename);
   TH1F *dZenDist;
   TH1F *dAziDist;
   TH2F *recoTrueZenDist;
   TH2F *recoTrueAziDist;
   TTree *dataTree;

   //double w = weight;
   //double zen = zen_true * 180.f / M_PI;
   //double azi = azi_true * 180.f / M_PI;
   double recZen = pt.theta * 180.f / M_PI;
   double recAzi = pt.phi   * 180.f / M_PI;
   summary->setRecoDir(recZen, recAzi);

   recoData dummyData;
   //dummyData.duplicate(summary);
   dummyData = *summary;
   //dummyData.getData( weight
   //                , zen_true * 180.f / M_PI, azi_true * 180.f / M_PI
   //                 , pt.theta * 180.f / M_PI, pt.phi   * 180.f / M_PI
   //                 , r_true, 0.f
   //                 , usedChan
   //                 , maxPixValue
   //                 );
/*
   printf("old weight: %f new weight %f\nold trueRadius: %f new trueRadius: %f\nold recoRadius: %f new recoRadius: %f\n           old trueZen: %f new trueZen: %f\nold trueAzi: %f new trueAzi: %f\nold maxPixIdx: %d new maxPixIdx: %d\nold maxPixCoherence: %f new maxPixCoherence: %f\n"
           , summary->weight, dummyData.weight
           , summary->trueRadius, dummyData.trueRadius
           , summary->recoRadius, dummyData.recoRadius
           , summary->trueZen, dummyData.trueZen
           , summary->trueAzi, dummyData.trueAzi
           , summary->maxPixIdx, dummyData.maxPixIdx
           , summary->maxPixCoherence, dummyData.maxPixCoherence
          );
   for(int ch=0; ch<16; ch++) printf("recoChan %d: old %d new %d", ch, summary->recoChan[ch], dummyData.recoChan[ch]);
*/


   if( !rootfile ){
   TFile fp(rootFilename, "NEW");
   dZenDist = new TH1F("recoZenDiff", "recoZenDiff", 360, -180, 180);
   dAziDist = new TH1F("recoAziDiff", "recoAziDiff", 720, -360, 360);
   recoTrueZenDist = new TH2F("recoTrueZenDist", "recoTrueZenDist", 180, 0, 180, 180, 0, 180);
   recoTrueAziDist = new TH2F("recoTrueAziDist", "recoTrueAziDist", 360, 0, 360, 360, 0, 360);

   dataTree = new TTree("dataTree", "dataTree");
   //dataTree->Branch("recoData", &dummyData, "weight/D:trueZen/F:trueAzi/F:recoZen/F:recoAzi/F:trueRadius/F:recoRadius/F:recoChan[16]/I:maxPixIdx/I:coherence/F");

   dZenDist->Write();
   dAziDist->Write();
   recoTrueZenDist->Write();
   recoTrueAziDist->Write();
   dataTree->Write();
   fp.Close();
   }
   //else{
   TFile fp_2(rootFilename, "UPDATE");
   dZenDist = (TH1F*)fp_2.Get("recoZenDiff");
   dAziDist = (TH1F*)fp_2.Get("recoAziDiff");
   recoTrueZenDist = (TH2F*)fp_2.Get("recoTrueZenDist");
   recoTrueAziDist = (TH2F*)fp_2.Get("recoTrueAziDist");
   dataTree = (TTree*)fp_2.Get("dataTree");
   //dataTree->ResetBranchAddresses();
   //dataTree->SetBranchAddress("recoData",&dummyData);

   cout<<"trueRadius: "<<dummyData.trueRadius<<endl;
   cout<<"trueZen: "<<dummyData.trueZen<<" recoZen: "<<dummyData.recoZen<<endl;
   cout<<"trueAzi: "<<dummyData.trueAzi<<" recoAzi: "<<dummyData.recoAzi<<endl;
   dZenDist->Fill(dZen);
   dAziDist->Fill(dAzi);
   //recoTrueZenDist->Fill(zen, recZen);
   //recoTrueAziDist->Fill(azi, recAzi);
   recoTrueZenDist->Fill(dummyData.trueZen, recZen);
   recoTrueAziDist->Fill(dummyData.trueAzi, recAzi);
   dataTree->Fill();

   //TFile fp_3(rootFilename, "RECREATE");
   dZenDist->Write();
   dAziDist->Write();
   recoTrueZenDist->Write();
   recoTrueAziDist->Write();
   //dataTree->Write();
   fp_2.Close();
   //fp_3.Close();

/*
 * Set flag condition HERE !!!!
 */

   //if( fabs(dAzi) > 10.f ) flag = 1;
   if( dZen > 15.f ) flag = 1;

   return flag;
}

//int record3DDiffGetFlag(recoData *summary, char *rootFilename){

/*
 * flag signals whether the event loop should break and save the current skymap, depending on whether the flag condition is met
 * flag 0: do not break
 * flag 1: break and save skymap
 */
 /*
   int flag = 0;

   int nSideExp = summary->onion->nSideExp;
   int nLayer   = summary->onion->nLayer;

   if(nSideExp < 0 || nSideExp > 7){ cerr<<"Invalid nSideExp\n"; return -1; }
   int nSide = pow(2, nSideExp);
   Healpix_Base hpBase = Healpix_Base(nSide, RING, SET_NSIDE);
   int nDir = hpBase.Npix();
   pointing pt;

   cout<<"maxPixIdx: "<<summary->maxPixIdx<<" maxPixIdx on the skymap: "<<summary->maxPixIdx%nDir<<endl;
   pt = hpBase.pix2ang( summary->maxPixIdx % nDir );
   //float dZen = (pt.theta - zen_true) * 180.f / M_PI;
   //float dAzi = (pt.phi   - azi_true) * 180.f / M_PI;
   float dZen = pt.theta * 180.f / M_PI - summary->trueZen;
   float dAzi = pt.phi   * 180.f / M_PI - summary->trueAzi;

   ifstream rootfile(rootFilename);
   TH1F *dZenDist;
   TH1F *dAziDist;
   TH2F *recoTrueZenDist;
   TH2F *recoTrueAziDist;
   TTree *onionTree;
   TTree *dataTree;

   //double w = weight;
   //double zen = zen_true * 180.f / M_PI;
   //double azi = azi_true * 180.f / M_PI;
   double recZen = pt.theta * 180.f / M_PI;
   double recAzi = pt.phi   * 180.f / M_PI;
   summary->setRecoDir(recZen, recAzi);
   summary->setRecoRadius( summary->onion->getLayerRadius(summary->maxPixIdx) );

   recoData *dummyData = new recoData;
   dummyData->duplicate(summary);
   vector<int>   * topMaxPixIdx             = &dummyData->topMaxPixIdx;
   vector<float> * topMaxPixCoherence       = &dummyData->topMaxPixCoherence;
   vector<int>   * maxPixIdxEachLayer       = &dummyData->maxPixIdxEachLayer;
   vector<float> * maxPixCoherenceEachLayer = &dummyData->maxPixCoherenceEachLayer;
   //dummyData.getData( weight
   //                , zen_true * 180.f / M_PI, azi_true * 180.f / M_PI
   //                 , pt.theta * 180.f / M_PI, pt.phi   * 180.f / M_PI
   //                 , r_true, 0.f
   //                 , usedChan
   //                 , maxPixValue
   //                 );
*//*
   printf("old weight: %f new weight %f\nold trueRadius: %f new trueRadius: %f\nold recoRadius: %f new recoRadius: %f\n           old trueZen: %f new trueZen: %f\nold trueAzi: %f new trueAzi: %f\nold maxPixIdx: %d new maxPixIdx: %d\nold maxPixCoherence: %f new maxPixCoherence: %f\n"
           , summary->weight, dummyData->weight
           , summary->trueRadius, dummyData->trueRadius
           , summary->recoRadius, dummyData->recoRadius
           , summary->trueZen, dummyData->trueZen
           , summary->trueAzi, dummyData->trueAzi
           , summary->maxPixIdx, dummyData->maxPixIdx
           , summary->maxPixCoherence, dummyData->maxPixCoherence
          );
   for(int ch=0; ch<16; ch++) printf("recoChan %d: old %d new %d", ch, summary->recoChan[ch], dummyData->recoChan[ch]);
   cout<<endl;
*/
/*
   if( !rootfile ){
   TFile fp(rootFilename, "NEW");
   dZenDist = new TH1F("recoZenDiff", "recoZenDiff", 360, -180, 180);
   dAziDist = new TH1F("recoAziDiff", "recoAziDiff", 720, -360, 360);
   recoTrueZenDist = new TH2F("recoTrueZenDist", "recoTrueZenDist", 180, 0, 180, 180, 0, 180);
   recoTrueAziDist = new TH2F("recoTrueAziDist", "recoTrueAziDist", 360, 0, 360, 360, 0, 360);

   onionTree = new TTree("onionTree", "onionTree");
   onionTree->Branch("nSideExp", &nSideExp, "nSideExp/I");
   onionTree->Branch("nLayer",   &nLayer, "nLayer/I");
   dataTree  = new TTree("dataTree", "dataTree");
   //dataTree->Branch("recoData", &dummyData, "weight/D:trueZen/F:trueAzi/F:recoZen/F:recoAzi/F:trueRadius/F:recoRadius/F:recoChan[16]/I:maxPixIdx/I:coherence/F:topN/I");
   dataTree->Branch("weight",     &dummyData->weight,          "weight/D");
   dataTree->Branch("trueZen",    &dummyData->trueZen,         "trueZen/F");
   dataTree->Branch("trueAzi",    &dummyData->trueAzi,         "trueAzi/F");
   dataTree->Branch("recoZen",    &dummyData->recoZen,         "recoZen/F");
   dataTree->Branch("recoAzi",    &dummyData->recoAzi,         "recoAzi/F");
   dataTree->Branch("trueRadius", &dummyData->trueRadius,      "trueRadius/F");
   dataTree->Branch("recoRadius", &dummyData->recoRadius,      "recoRadius/F");
   dataTree->Branch("recoChan",   &dummyData->recoChan,        "recoChan[16]/I");
   dataTree->Branch("maxPixIdx",  &dummyData->maxPixIdx,       "maxPixIdx/I");
   dataTree->Branch("coherence",  &dummyData->maxPixCoherence, "coherence/F");
   dataTree->Branch("topN",       &dummyData->topN,            "topN/I");
   dataTree->Branch("topMaxPixIdx",             &dummyData->topMaxPixIdx);
   dataTree->Branch("topMaxPixCoherence",       &dummyData->topMaxPixCoherence);
   dataTree->Branch("maxPixIdxEachLayer",       &dummyData->maxPixIdxEachLayer);
   dataTree->Branch("maxPixCoherenceEachLayer", &dummyData->maxPixCoherenceEachLayer);
   dataTree->Branch("likelihood", &dummyData->likelihood,      "likelihood/D");
   dataTree->Branch("pValue",     &dummyData->pValue,          "pValue/D");
   dataTree->Branch("inWindowSNR", &dummyData->inWindowSNR,    "inWindowSNR/F");
   dataTree->Branch("unmodSNR",    &dummyData->unmodSNR,       "unmodSNR/F");

   //cout<<"4765\n";
   //dataTree->Branch("recoData","recoData",&dummyData,128000,0);
*/
   /* onionTree only needs to be filled once, at creation */
/*   onionTree->Fill();

   dZenDist->Write();
   dAziDist->Write();
   recoTrueZenDist->Write();
   recoTrueAziDist->Write();
   onionTree->Write();
   dataTree->Write();
   fp.Close();
   }
   //else{
   TFile fp_2(rootFilename, "UPDATE");
   dZenDist = (TH1F*)fp_2.Get("recoZenDiff");
   dAziDist = (TH1F*)fp_2.Get("recoAziDiff");
   recoTrueZenDist = (TH2F*)fp_2.Get("recoTrueZenDist");
   recoTrueAziDist = (TH2F*)fp_2.Get("recoTrueAziDist");
   dataTree = (TTree*)fp_2.Get("dataTree");
   //dataTree->ResetBranchAddresses();
   //cout<<"4784\n";
   //dataTree->SetBranchAddress("recoData",                 &dummyData);
   dataTree->SetBranchAddress("weight",     &dummyData->weight);
   dataTree->SetBranchAddress("trueZen",    &dummyData->trueZen);
   dataTree->SetBranchAddress("trueAzi",    &dummyData->trueAzi);
   dataTree->SetBranchAddress("recoZen",    &dummyData->recoZen);
   dataTree->SetBranchAddress("recoAzi",    &dummyData->recoAzi);
   dataTree->SetBranchAddress("trueRadius", &dummyData->trueRadius);
   dataTree->SetBranchAddress("recoRadius", &dummyData->recoRadius);
   dataTree->SetBranchAddress("recoChan",   &dummyData->recoChan);
   dataTree->SetBranchAddress("maxPixIdx",  &dummyData->maxPixIdx);
   dataTree->SetBranchAddress("coherence",  &dummyData->maxPixCoherence);
   dataTree->SetBranchAddress("topN",       &dummyData->topN);
   dataTree->SetBranchAddress("topMaxPixIdx",             &topMaxPixIdx);
   dataTree->SetBranchAddress("topMaxPixCoherence",       &topMaxPixCoherence);
   dataTree->SetBranchAddress("maxPixIdxEachLayer",       &maxPixIdxEachLayer);
   dataTree->SetBranchAddress("maxPixCoherenceEachLayer", &maxPixCoherenceEachLayer);
   dataTree->SetBranchAddress("likelihood", &dummyData->likelihood);
   dataTree->SetBranchAddress("pValue",     &dummyData->pValue);
   dataTree->SetBranchAddress("inWindowSNR",&dummyData->inWindowSNR);
   dataTree->SetBranchAddress("unmodSNR",   &dummyData->unmodSNR);

   cout<<"trueRadius: "<<dummyData->trueRadius<<" recoRadius: "<<dummyData->recoRadius<<endl;
   cout<<"trueZen: "<<dummyData->trueZen<<" recoZen: "<<dummyData->recoZen<<endl;
   cout<<"trueAzi: "<<dummyData->trueAzi<<" recoAzi: "<<dummyData->recoAzi<<endl;
   dZenDist->Fill(dZen);
   dAziDist->Fill(dAzi);
   //recoTrueZenDist->Fill(zen, recZen);
   //recioTrueAziDist->Fill(azi, recAzi);
   recoTrueZenDist->Fill(dummyData->trueZen, recZen);
   recoTrueAziDist->Fill(dummyData->trueAzi, recAzi);
   dataTree->Fill();

   //TFile fp_3(rootFilename, "RECREATE");
   dZenDist->Write();
   dAziDist->Write();
   recoTrueZenDist->Write();
   recoTrueAziDist->Write();
   dataTree->Write();
   fp_2.Close();
   //}//end of else
   //fp_3.Close();
   delete dummyData;
*//*
 * Set flag condition HERE !!!!
 */

   //if( fabs(dAzi) > 10.f ) flag = 1;
/*   if( dZen > 15.f ) flag = 1;

   return flag;
}
*/
int record3DDiffGetFlag(recoSettings *settings, recoData *summary, TH1F *dZenDist, TH1F *dAziDist, TH2F *recoTrueZenDist, TH2F *recoTrueAziDist){

   cout<<"*********** record3DDiffGetFlag ***********\n";

   int flag = 0;

   int nSideExp = /*summary->onion*/settings->nSideExp;
   int nLayer   = /*summary->onion*/settings->nLayer;

   //if(nSideExp < 0 || nSideExp > 7){ cerr<<"Invalid nSideExp\n"; return -1; }
   int nSide = pow(2, nSideExp);
   Healpix_Base hpBase = Healpix_Base(nSide, HEALPIX_ORDERING, SET_NSIDE);
   int nDir = hpBase.Npix();
   pointing pt;

   Healpix_Onion onion(nSideExp, nLayer);

   cout<<"maxPixIdx: "<<summary->maxPixIdx<<" maxPixIdx on the skymap: "<<summary->maxPixIdx%nDir<<endl;
   pt = hpBase.pix2ang( summary->maxPixIdx % nDir );
   //float dZen = (pt.theta - zen_true) * 180.f / M_PI;
   //float dAzi = (pt.phi   - azi_true) * 180.f / M_PI;

   double recZen = pt.theta * 180.f / M_PI;
   double recAzi = pt.phi   * 180.f / M_PI;
   float dZen = recZen - summary->trueZen;
   float dAzi = recAzi - summary->trueAzi;

   //double w = weight;
   //double zen = zen_true * 180.f / M_PI;
   //double azi = azi_true * 180.f / M_PI;
   //double recZen = pt.theta * 180.f / M_PI;
   //double recAzi = pt.phi   * 180.f / M_PI;
   summary->setRecoDir(recZen, recAzi);
   summary->setRecoRadius( /*summary->onion->*/onion.getLayerRadius(summary->maxPixIdx) );

   dZenDist->Fill(dZen);
   dAziDist->Fill(dAzi);
   //recoTrueZenDist->Fill(zen, recZen);
   //recioTrueAziDist->Fill(azi, recAzi);
   recoTrueZenDist->Fill(summary->trueZen, recZen);
   recoTrueAziDist->Fill(summary->trueAzi, recAzi);

/*
 * Set flag condition here!!!
 */

   if( fabs(dZen) > 15.f ) flag = 1;

return flag;
}

int record3DZoomedDiffGetFlag(recoSettings *settings, recoData *summary, TH1F *dZenDist, TH1F *dAziDist, TH2F *recoTrueZenDist, TH2F *recoTrueAziDist){

   cout<<"******** record3DZoomedDiffGetFlag **********\n";

   int flag = 0;

   int nSideExp = /*summary->onion*/settings->/*nSideExpStart*/nSideExpEnd;
   int nLayer   = /*summary->onion*/settings->nLayer;

   //if(nSideExp < 0 || nSideExp > 7){ cerr<<"Invalid nSideExp\n"; return -1; }
   int nSide = pow(2, nSideExp);
   Healpix_Base hpBase = Healpix_Base(nSide, HEALPIX_ORDERING, SET_NSIDE);
   int nDir = hpBase.Npix();
   pointing pt;

   Healpix_Onion onion(nSideExp, nLayer);

   cout<<"maxPixIdx: "<<summary->maxPixIdx<<" maxPixIdx on the skymap: "<<summary->maxPixIdx%nDir<<endl;
   pt = hpBase.pix2ang( summary->maxPixIdx % nDir );
   //float dZen = (pt.theta - zen_true) * 180.f / M_PI;
   //float dAzi = (pt.phi   - azi_true) * 180.f / M_PI;

   double recZen = pt.theta * 180.f / M_PI;
   double recAzi = pt.phi   * 180.f / M_PI;
   float dZen = recZen - summary->trueZen;
   float dAzi = recAzi - summary->trueAzi;

   //double w = weight;
   //double zen = zen_true * 180.f / M_PI;
   //double azi = azi_true * 180.f / M_PI;
   //double recZen = pt.theta * 180.f / M_PI;
   //double recAzi = pt.phi   * 180.f / M_PI;
   summary->setRecoDir(recZen, recAzi);
   summary->setRecoRadius( /*summary->onion->*/onion.getLayerRadius(summary->maxPixIdx) );

   dZenDist->Fill(dZen);
   dAziDist->Fill(dAzi);
   //recoTrueZenDist->Fill(zen, recZen);
   //recioTrueAziDist->Fill(azi, recAzi);
   recoTrueZenDist->Fill(summary->trueZen, recZen);
   recoTrueAziDist->Fill(summary->trueAzi, recAzi);

/*
 * Set flag condition here!!!
 */

   if( fabs(dZen) > 15.f ) flag = 1;

return flag;
}

float getSpaceAngle(float theta1, float phi1, float theta2, float phi2){

return acos( sin(theta1)*cos(phi1)*sin(theta2)*cos(phi2) +
             sin(theta1)*sin(phi1)*sin(theta2)*sin(phi2) +
             cos(theta1)*cos(theta2)
           );
}

void stackXCorrAroundPeak(const TGraph *gr, TH1F *hist, float plusMinusTime){

   double t, v;
   double max = -1e10;
   double maxTime;
   for(int s=0; s<gr->GetN(); s++){

      gr->GetPoint(s,t,v);
      if( v > max ) { max = v; maxTime = t; }
   }

   for(int s=0; s<gr->GetN(); s++){
      gr->GetPoint(s,t,v);
      if( fabs(t - maxTime) < plusMinusTime ) hist->Fill(static_cast<float>(t-maxTime), static_cast<float>(v / max));
   }

}

int doNchnlScan(const double eventWeight, const vector<TGraph *>& cleanEvent, TH2F *mnMap, int *nchnlArray, const int *chanMask, int *goodChan, const int nThresStep, const double minThres, const double maxThres){

double nchnl_threshold;
//TH2F *mnMap = new TH2F("mnMap","mnMap",17,-0.5,16.5,nThresStep,minThres-((maxThres-minThres)/(double)nThresStep)*0.5,maxThres-((maxThres-minThres)/(double)nThresStep)*0.5);
if(mnMap == NULL){ cerr<<"Need to pass *TH2F to doNchnlScan\n"; return -1; }

for(int step=0; step<nThresStep; step++){

   nchnl_threshold = minThres + (double)step * (maxThres - minThres) / (double)nThresStep ;

   getNchnlMask(cleanEvent, nchnl_threshold, nchnlArray, chanMask, goodChan);

   for(int n=0; n<17; n++){
   if( nchnlArray[0] >= n ) mnMap->Fill(n, nchnl_threshold, eventWeight);
   }

}

return 0;
}

int doNchnlScan(const double eventWeight, const vector<TGraph *>& cleanEvent, TH2F *mnMap, TH2F *mnMap_V, TH2F *mnMap_H, int *nchnlArray, const int *chanMask, int *goodChan, const int nThresStep, const double minThres, const double maxThres){

double nchnl_threshold;
//TH2F *mnMap = new TH2F("mnMap","mnMap",17,-0.5,16.5,nThresStep,minThres-((maxThres-minThres)/(double)nThresStep)*0.5,maxThres-((maxThres-minThres)/(double)nThresStep)*0.5);
if(mnMap == NULL || mnMap_V == NULL || mnMap_H == NULL){ cerr<<"Need to pass 3 *TH2F to doNchnlScan\n"; return -1; }

for(int step=0; step<nThresStep; step++){

   nchnl_threshold = minThres + (double)step * (maxThres - minThres) / (double)nThresStep ;

   getNchnlMask(cleanEvent, nchnl_threshold, nchnlArray, chanMask, goodChan);

   for(int n=0; n<17; n++){
   if( nchnlArray[0] >= n ) mnMap->Fill(n, nchnl_threshold, eventWeight);
   }
   for(int n=0; n<9; n++){
   if( nchnlArray[1] >= n ) mnMap_V->Fill(n, nchnl_threshold, eventWeight);
   if( nchnlArray[2] >= n ) mnMap_H->Fill(n, nchnl_threshold, eventWeight);
   }

}

return 0;
}


int computeMapLikelihoodAndPValue(const int nDir, const int nLayer, const char *fitFunc, const char *fitFuncFile, float *mapData, double& likelihood, double& pValue){

   //int nDir = 12 * pow(2, nSideExp) * pow(2, nSideExp);
   char pixname[200];

   //ifstream ifs(fitFuncFile);
   //if(!ifs) { cerr<<"No fitFuncFile!!\n"; return -1; }
   TFile ref(fitFuncFile);
   TF1 /**expo,*/ *normexpo; //expo: fitted exponential fx. normexpo: normalized expo. A probability distribution
   TF1 *normgaus;

   //double likelihood, pValue,
   likelihood = pValue = 0.;
   int pixCnt = 0;
   double llh;
   double a; //e^(b+ax), ie. parameter no. 2
   double mean, sigma;

   //for(int pix=0; pix<nDir*nLayer; pix++){

      //sprintf(pixname, "expo_pix_%d", pix);
      //if(ref.GetListOfKeys()->Contains(pixname)){

      //pixCnt++;

      //expo = (TF1*)ref.Get(pixname);
   if( fitFunc == "expo" ){

      for(int pix=0; pix<nDir*nLayer; pix++){
      sprintf(pixname, "normexpo_pix_%d", pix);
      if(ref.GetListOfKeys()->Contains(pixname)){

      pixCnt++;
      normexpo = (TF1*)ref.Get(pixname);
      a = normexpo->GetParameter(2);
      //likelihood = (likelihood*(float)(pixCnt-1) + log(normexpo->Eval(mapData[pix])))/(float)pixCnt;
      //llh = log(normexpo->Eval(mapData[pix]));
      //likelihood += log(normexpo->Eval(mapData[pix]));
      llh = log(-1. * a) + a * mapData[pix];
      likelihood += llh;
      if( !isfinite(llh) ) cout<<"pix: "<<pix<<" llh: "<<llh<<" normexpo eval: "<<normexpo->Eval(mapData[pix])<<" mapData: "<<mapData[pix]<<endl;
      //cout<<"likelihood: "<<likelihood<<endl;
      //pValue = (pValue*(float)(pixCnt-1) + log(exp(normexpo->GetParameter(2)*mapData[pix])))/(float)pixCnt;
      pValue += a * mapData[pix];
      //cout<<"pValue: "<<pValue<<endl;
      //delete expo;
      delete normexpo;
      }
      }//end of pix
   }

   else if( fitFunc == "gaus"){

      for(int pix=0; pix<nDir*nLayer; pix++){
      sprintf(pixname, "normgaus_pix_%d", pix);
      if(ref.GetListOfKeys()->Contains(pixname)){

      pixCnt++;
      normgaus = (TF1*)ref.Get(pixname);
      mean  = normgaus->GetParameter(1);
      sigma = normgaus->GetParameter(2);

      llh = -log(sqrt(2*TMath::Pi())*sigma) - 0.5 * pow((mapData[pix] - mean)/sigma, 2.);
      //llh = log(normgaus->Eval(mapData[pix]));
      likelihood += llh;
      if(!isfinite(llh) ) cout<<"pix: "<<pix<<" llh: "<<llh<<" normgaus eval: "<<normgaus->Eval(mapData[pix])<<" mapData: "<<mapData[pix]<<endl;

      //pValue += log((1. - TMath::Erf((mapData[pix]-mean)/(sigma*sqrt(2.))))/2.);

      delete normgaus;
      }
      }//end of pix
   }
   //}//end of pix

   likelihood /= (double)pixCnt;
   pValue     /= (double)pixCnt;
   //ifs.close();
   ref.Close();
   //delete expo;
   //delete normexpo;
//
return 0;
}
