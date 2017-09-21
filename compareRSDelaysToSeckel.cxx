#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TObject.h"
#include "TPolyMarker3D.h"
#include "TH3D.h"

#include "calibrationTools.h"
#include "calibrationToolsVs3.h"
#include "evProcessTools.h"
#include "Healpix_Onion.h"
#include "recoTools.h"
#include "recoSettings.h"
#include "recoData.h"

int main(int argc, char** argv){

/* Read in Seckel geometry (ant & src locations) and delays (1st and 2nd rays) */

vector<double> srcPosVec;
vector<vector<float> > recoDelaysVec, recoDelaysVec_V, recoDelaysVec_H, recoDelaysVec_ctr;

int err = getRecoDelaysFromSeckel("raygrid.txt", srcPosVec, recoDelaysVec_ctr, recoDelaysVec, recoDelaysVec_V, recoDelaysVec_H);
if( err<0 ){ cerr<<"Error computing reco delays\n"; return -1; }

vector<vector<double> > antLocation;
err = getSeckelStationGeometry(antLocation);
if( err<0 ){ cerr<<"Error calibrating geometry and delays\n"; return -1; }

/* Call RayTrace object, feed Seckel ant & src locations to it, and obtain delays */

int nAnt = 16;
float *recoDelays = (float*)malloc(11*11*11*nAnt*sizeof(float));
float *recoRefracDelays = (float*)malloc(11*11*11*nAnt*sizeof(float));
err = compute3DRecoBothDelaysWithRadioSplineWithSeckelGeom(srcPosVec, antLocation,
                                          recoDelays, recoRefracDelays);

/* Compare delays */

TCanvas *cvs;
TH2F *d1stDelay[11], *d2ndDelay[11];
char histname[200];
char cvsname[200];

int ch = 0;

for(int layer=0; layer<11; layer++){
  //for(int ch=0; ch<nAnt; ch++){

    snprintf(histname,sizeof(histname),"direct_layer%d_ch%d",layer, ch);
    d1stDelay[layer] = new TH2F(histname, histname,11,-0.5,10.5,11,-0.5,10.5);
    snprintf(histname,sizeof(histname),"refracted_layer%d_ch%d",layer, ch);
    d2ndDelay[layer] = new TH2F(histname, histname,11,-0.5,10.5,11,-0.5,10.5);

    for(int theta=0; theta<11; theta++){
      for(int phi=0; phi<11; phi++){

        d1stDelay[layer]->Fill(phi,theta,recoDelaysVec[0][(layer*11*11+theta*11+phi)*nAnt+ch] - recoDelays[(layer*11*11+theta*11+phi)*nAnt+ch]);
        d2ndDelay[layer]->Fill(phi,theta,recoDelaysVec[1][(layer*11*11+theta*11+phi)*nAnt+ch] - recoRefracDelays[(layer*11*11+theta*11+phi)*nAnt+ch]);

      }
    }

    cvs = new TCanvas("cvs","cvs",800,600);
    snprintf(cvsname,sizeof(cvsname),"layer_%d.C",layer);
    cvs->Divide(2,1);
    cvs->cd(1);
    d1stDelay[layer]->Draw("colz");
    cvs->cd(2);
    d2ndDelay[layer]->Draw("colz");
    cvs->SaveAs(cvsname);
    delete cvs;
  //}
}

return 0;
}
