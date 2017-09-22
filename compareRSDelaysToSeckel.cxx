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
#include "TRandom3.h"

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

TH1F *hist=new TH1F("hist","hist",500,-20,20);
TH1F *directHist[16*16], *refractHist[16*16];
//directHist = new TH1F("directHist","directHist",100,-5,5);
//refractHist = new TH1F("refractHist","refractHist",100,-5,5);
int ch1 /*= 0*/;
int ch2 /*= 4*/;

cvs = new TCanvas("cvs","cvs",800,600);
cvs->Divide(2,1);
TRandom3 *rnd = new TRandom3();
cvs->cd(1); hist->Draw(); hist->SetStats(0); hist->GetYaxis()->SetRangeUser(0,1400);
cvs->cd(2); hist->Draw(); hist->SetStats(0); hist->GetYaxis()->SetRangeUser(0,1400);
int ci; //color index

for(ch1=0; ch1<16;  ch1++){
for(ch2=0; ch2<16; ch2++){

snprintf(histname,sizeof(histname),"direct_%d_%d",ch1,ch2);
directHist[ch1*16+ch2]=new TH1F(histname, histname, 500,-20,20);
snprintf(histname,sizeof(histname),"refract_%d_%d",ch1,ch2);
refractHist[ch1*16+ch2]=new TH1F(histname, histname, 500,-20,20);

for(int layer=0; layer<11; layer++){
  //for(int ch=0; ch<nAnt; ch++){
    /*
    snprintf(histname,sizeof(histname),"direct_layer%d_ch%d_ch%d",layer, ch1, ch2);
    d1stDelay[layer] = new TH2F(histname, histname,11,-0.5,10.5,11,-0.5,10.5);
    snprintf(histname,sizeof(histname),"refracted_layer%d_ch%d_ch%d",layer, ch1, ch2);
    d2ndDelay[layer] = new TH2F(histname, histname,11,-0.5,10.5,11,-0.5,10.5);
    */
    for(int theta=0; theta<11; theta++){
      for(int phi=0; phi<11; phi++){

        //d1stDelay[layer]->Fill(phi,theta,(recoDelaysVec[0][(layer*11*11+theta*11+phi)*nAnt+ch1] - recoDelaysVec[0][(layer*11*11+theta*11+phi)*nAnt+ch2]) -
        //                                 (recoDelays[(layer*11*11+theta*11+phi)*nAnt+ch1] - recoDelays[(layer*11*11+theta*11+phi)*nAnt+ch2])
        //                      );
        directHist[ch1*16+ch2]->Fill((recoDelaysVec[0][(layer*11*11+theta*11+phi)*nAnt+ch1] - recoDelaysVec[0][(layer*11*11+theta*11+phi)*nAnt+ch2]) -
                                     (recoDelays[(layer*11*11+theta*11+phi)*nAnt+ch1] - recoDelays[(layer*11*11+theta*11+phi)*nAnt+ch2])
                                    );
        //cout<<recoDelays[(layer*11*11+theta*11+phi)*nAnt+ch]<<endl;
        //d2ndDelay[layer]->Fill(phi,theta,(recoDelaysVec[1][(layer*11*11+theta*11+phi)*nAnt+ch1] - recoDelaysVec[1][(layer*11*11+theta*11+phi)*nAnt+ch2]) -
        //                                 (recoRefracDelays[(layer*11*11+theta*11+phi)*nAnt+ch1] - recoRefracDelays[(layer*11*11+theta*11+phi)*nAnt+ch2])
        //                      );
        refractHist[ch1*16+ch2]->Fill((recoDelaysVec[1][(layer*11*11+theta*11+phi)*nAnt+ch1] - recoDelaysVec[1][(layer*11*11+theta*11+phi)*nAnt+ch2]) -
                                      (recoRefracDelays[(layer*11*11+theta*11+phi)*nAnt+ch1] - recoRefracDelays[(layer*11*11+theta*11+phi)*nAnt+ch2])
                                     );
        //cout<<recoRefracDelays[(layer*11*11+theta*11+phi)*nAnt+ch]<<endl;

      }
    }
    /*
    cvs = new TCanvas("cvs","cvs",800,600);
    snprintf(cvsname,sizeof(cvsname),"layer_%d.C",layer);
    cvs->Divide(2,1);
    cvs->cd(1);
    d1stDelay[layer]->Draw("colz");
    cvs->cd(2);
    d2ndDelay[layer]->Draw("colz");
    cvs->SaveAs(cvsname);
    delete cvs;
    */
  //}
}

if(ch1!=ch2){

ci = rnd->Rndm()*29+20;

cvs->cd(1);
directHist[ch1*16+ch2]->SetLineColor(ci);
directHist[ch1*16+ch2]->Draw("same");

cvs->cd(2);
refractHist[ch1*16+ch2]->SetLineColor(ci);
refractHist[ch1*16+ch2]->Draw("same");

}

}
}

cvs->SaveAs("RS_Seckel_comparison.C");

return 0;
}
