#include <iostream>
#include "evProcessTools.h"
#include "TGraph.h"
#include "TFile.h"
#include "TRandom3.h"

int main(int argc, char** argv){

  if(argc<5){ cerr<<"Usage: 1.wfLen 2.wInt 3.range 4.ROOT filename"<<endl; return -1; }

  TFile fp(argv[4],"RECREATE");
  const int nchan = 16;
  int wfLen = atoi(argv[1]); //number of samples
  double wInt = atof(argv[2]); // sample bin size
  double range = atof(argv[3]); // +-range
  TGraph *gr[nchan];

  TRandom3 *rnd = new TRandom3(100);
  char grName[200];

  for(int ch=0; ch<nchan; ch++){

    gr[ch] = evProcessTools::getRandomVoltageGraph(wfLen, wInt, range, rnd);
    snprintf(grName,sizeof(char)*200,"gr_%d",ch);
    gr[ch]->SetName(grName);
    gr[ch]->Write();
  }

  fp.Close();

}
