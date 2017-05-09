#ifndef HEALPIX_ONION_H
#define HEALPIX_ONION_H
/*
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <cmath>
*/
/* OpenCL/CLFFT includes */
//#include <CL/opencl.h>
//#include <clFFT.h>

/* radiospline includes */
//#include "radiospline/IceGeometry.h"
//#include "radiospline/RayDelay.h"
//#include "radiospline/FirnShadow.h"

//#define AIR_FILE "InAir1.fits"
//#define ICE_FILE "InIce1.fits"
//#define AIR_FILE "InAirFinal_mod.fits"
//#define ICE_FILE "InIceFinal_mod.fits"
//#define SHADOW_FILE "shadow_5km.fits"
//#define AIR_FILE "delay_inair.fits"
//#define ICE_FILE "delay_inice.fits"
//#define SHADOW_FILE "firn_shadow.fits"

/* CFITSIO include */
//#include "fitsio.h"

/* Healpix includes */
#include "healpix_base.h"
#include "healpix_data_io.h"
#include "healpix_map_fitsio.h"
#include "healpix_map.h"
#include "fitshandle.h"
#include "pointing.h"

#include "vec3.h"
#include "arr.h"

/* ROOT includes */
//#include "TGraph.h"
//#include "TCanvas.h"
//#include "TGaxis.h"
#include "TObject.h"

//#include "FFTtools.h"
//#include "evProcessTools.h"
//#include "recoData.h"
//#include "recoSettings.h"
/*
#define C_INV 3.34
#define speedOfLight 0.3 // m/ns
#define nIce 1.76 // ~-180m
*/
//#define REFERENCE_MAP_FIT_FILE "testFitFuntFile.root"
//#define REFERENCE_MAP_FIT_FILE "testFitFuncFile_2013_A3_pulserSweep_run410-430.root"
#define HEALPIX_ORDERING NEST //could be RING or NEST

//#define REFERENCE_MAP_FIT_FILE "testFitFuncFile_2014_ARA03_run3623.root"
//#define FIT_FUNC "gaus"
/* If bandpass defined, waveforms will be bandpassed before being cross-correlated in the freq domain */
//#define bandpass

using namespace std;
/*
#ifndef XCORRSUMGRAPH
TGraph *sillygr = new TGraph();
//TGraph *envelopeSum = new TGraph();
#define XCORRSUMGRAPH
#endif
*/

class Healpix_Onion : public TObject
{
private:

protected:

public:



   int nSideExp;
   int nLayer;
   int nDir;
   Healpix_Base hpBase;
   vector<float> layerRadii;

   Healpix_Onion();//{ initialize(); }
   ~Healpix_Onion();//{ }

   void initialize();/*{
   nSideExp = 2;
   nLayer   = 1;
   hpBase = Healpix_Base(pow(2,nSideExp), RING, SET_NSIDE);
   nDir = hpBase.Npix();
   layerRadii.clear();
   layerRadii.push_back(3000);

   }
  */
   //~Healpix_Onion();

   Healpix_Onion(int _nSideExp, int _nLayer, float _layerFirstRadius, float _layerLastRadius);/*{

   nSideExp = _nSideExp;
   nLayer   = _nLayer;
   hpBase = Healpix_Base(pow(2,nSideExp), RING, SET_NSIDE);
   nDir = hpBase.Npix();
   layerRadii.clear();
   float r=0.f;

   //while( r<= 5000.f ){

   //   r += 5000.f / (float)nLayer;
   //   layerRadii.push_back( r );

   //}


   for(int i=1; i<=nLayer; i++){

   //layerRadii.push_back( (float)i * 5000.f / (float)nLayer);
   layerRadii.push_back( (float)i * 3000.f / (float)nLayer); //for testing 3D calpulser reco

   }

   if( (int)layerRadii.size() != nLayer ) cerr<<"Warning!! layerRadii.size(): "<<layerRadii.size()<<" nLayer: "<<nLayer<<endl;

   }
*/
   pointing getPointing(int pixNum);/*{

   return hpBase.pix2ang( pixNum % nDir );

   }
*/
   int getLayerNumber(int pixNum);/*{

   return pixNum / nDir ;

   }
*/
   float getLayerRadius(int pixNum);/*{

   return layerRadii[ getLayerNumber(pixNum) ];

   }
*/
   //ClassDef(Healpix_Onion, 1);
};

#endif
