#include "Healpix_Onion.h"

//ClassImp(Healpix_Onion);

using namespace std;

   Healpix_Onion::Healpix_Onion(){ Healpix_Onion::initialize(); }
   Healpix_Onion::~Healpix_Onion(){ }

   void Healpix_Onion::initialize(){
   nSideExp = 2;
   nLayer   = 1;
   hpBase = Healpix_Base(pow(2,nSideExp), HEALPIX_ORDERING/*RING*/, SET_NSIDE);
   nDir = hpBase.Npix();
   layerRadii.clear();
   layerRadii.push_back(3000);

   }

   Healpix_Onion::Healpix_Onion(int _nSideExp, int _nLayer){

   nSideExp = _nSideExp;
   nLayer   = _nLayer;
   hpBase = Healpix_Base(pow(2,nSideExp), HEALPIX_ORDERING/*RING*/, SET_NSIDE);
   nDir = hpBase.Npix();
   layerRadii.clear();
   float r=0.f;

  for(int i=1; i<=nLayer; i++){

   layerRadii.push_back( (float)i * /*5000.f*/42.f / (float)nLayer); //for testing 3D calpulser reco

   }

   if( (int)layerRadii.size() != nLayer ) cerr<<"Warning!! layerRadii.size(): "<<layerRadii.size()<<" nLayer: "<<nLayer<<endl;

   }

   pointing Healpix_Onion::getPointing(int pixNum){

   return hpBase.pix2ang( pixNum % nDir );

   }

   int Healpix_Onion::getLayerNumber(int pixNum){

   return pixNum / nDir ;

   }

  float Healpix_Onion::getLayerRadius(int pixNum){

   return layerRadii[ getLayerNumber(pixNum) ];

   }
