////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////  exampleLoop.cxx
////      Just a very simple example that loops over RawAraEvent objects
////      calibrating them to make a UsefulAraEvent
////
////    Feb 2011,  rjn@hep.ucl.ac.uk
////////////////////////////////////////////////////////////////////////////////

//Includes
#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <dirent.h>


//AraRoot Includes
#include "RawIcrrStationEvent.h"

#include "RawAtriStationEvent.h"
#include "UsefulAraStationEvent.h"
#include "UsefulIcrrStationEvent.h"
#include "UsefulAtriStationEvent.h"

//Include FFTtools.h if you want to ask the correlation, etc. tools
#include "FFTtools.h"
#include "AraGeomTool.h"

//ROOT Includes
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TSystem.h"
#include "TMath.h"
#include "TF1.h"
#include "TChain.h"
#include "TH2D.h"
#include "TLine.h"
#include "TFitResult.h"
#include "TVirtualFitter.h"
#include "TArc.h"

//selfmade includes
///#include "filter_alg.h"
/////#include "WaveProcessing.h"
#include "calibrationToolsVs3.h"
//#include "interFerometryTools.h"
//#include "timeSequenceFilter.h"
//#include "filterTools.h"
//#include "analysisTools.h"
//For pdestal file search
#include <sys/types.h>
#include <dirent.h>
#include <vector>

#include "evProcessTools.h"

using namespace std;


Double_t evProcessTools::integratePower(TH1D* histo)
{

  Double_t integral=0;
  Double_t power;
  Int_t firstBin=1;
  Int_t lastBin=histo->GetEntries();
  Double_t df = histo->GetBinWidth(firstBin);
  for(int i=firstBin;i<=lastBin;i++) {
    power =  histo->GetBinContent(i);
    integral+=power*df;
  }
  return integral;
}

TGraph *evProcessTools::getfftfromgraph(TGraph *gr, Double_t intSample, const int maxSample)
{
   //   static AraGeomTool *fGeomTool = AraGeomTool::Instance();
   if(!gr) return NULL;
   Double_t newX[maxSample],newY[maxSample];
   Int_t maxSamps=maxSample;

   //Added by Lu. This step is in getFFTForRFChan routine
   //TGraph *grInt = FFTtools::getInterpolatedGraph(gr, intSample);
   //TGraph *grInt = (TGraph*)gr->Clone();

   //Note that if the graph is not interpolated in this routine, it should be interpolated _before_ being passed to this function

   Int_t numSamps  = gr->GetN();
   Double_t xVals;
   Double_t yVals;
   for(int i=0;i<maxSamps;i++) {
      if(i<numSamps) {
	 gr->GetPoint(i,xVals,yVals);
	 newX[i]=xVals;
	 newY[i]=yVals*FFTtools::bartlettWindow(i,numSamps);
      }
      else {
	 newX[i]=newX[i-1]+intSample;
	 newY[i]=0;
      }
  }
   TGraph *grNew = new TGraph(maxSamps,newX,newY);
//   TGraph *grFFT = FFTtools::makePowerSpectrumVoltsSecondsBartlett(grNew);
//   makePowerSpectrumMilliVoltsNanoSeconds
   TGraph *grFFT = FFTtools::makePowerSpectrumMilliVoltsNanoSeconds(grNew); //dB is the one used in getFFTForRFChan
                                                                              //non-dB is used in Thomas' original code
//   TGraph *grFFT = FFTtools::makePowerSpectrumVoltsSeconds(grNew);
   delete gr;
   delete grNew;
  // delete xVals;
//   delete yVals;
   return grFFT;
}

TGraph *evProcessTools::getfftdbfromgraph(TGraph *gr, Double_t intSample, const int maxSample)
{
   //   static AraGeomTool *fGeomTool = AraGeomTool::Instance();
   if(!gr) return NULL;
   Double_t newX[maxSample],newY[maxSample];
   Int_t maxSamps=maxSample;

   //Added by Lu. This step is in getFFTForRFChan routine
   //TGraph *grInt = FFTtools::getInterpolatedGraph(gr, intSample);
   //TGraph *grInt = (TGraph*)gr->Clone();

   //Note that if the graph is not interpolated in this routine, it should be interpolated _before_ being passed to this function

   Int_t numSamps  = gr->GetN();
   Double_t xVals;
   Double_t yVals;
   for(int i=0;i<maxSamps;i++) {
      if(i<numSamps) {
	 gr->GetPoint(i,xVals,yVals);
	 newX[i]=xVals;
	 newY[i]=yVals*FFTtools::bartlettWindow(i,numSamps);
      }
      else {
	 newX[i]=newX[i-1]+intSample;
	 newY[i]=0;
      }
  }
   TGraph *grNew = new TGraph(maxSamps,newX,newY);
//   TGraph *grFFT = FFTtools::makePowerSpectrumVoltsSecondsBartlett(grNew);
//   makePowerSpectrumMilliVoltsNanoSeconds
   TGraph *grFFT = FFTtools::makePowerSpectrumMilliVoltsNanoSecondsdB(grNew); //dB is the one used in getFFTForRFChan
                                                                              //non-dB is used in Thomas' original code
//   TGraph *grFFT = FFTtools::makePowerSpectrumVoltsSeconds(grNew);
   //delete gr;
   delete grNew;
  // delete xVals;
//   delete yVals;
   return grFFT;
}

TGraph *evProcessTools::getBartlettAndPaddedGraph(TGraph *gr, const int maxSample)
{

   //   static AraGeomTool *fGeomTool = AraGeomTool::Instance();
   if(!gr) return NULL;
   Double_t newX[maxSample],newY[maxSample];
   Int_t maxSamps=maxSample;

   Double_t time1=0;
   Double_t time2=0;
   Double_t xVals;
   Double_t yVals;

   gr->GetPoint(0,time1,yVals);
   gr->GetPoint(1,time2,yVals);

  Double_t intSample = time2 - time1;
   cerr << "The interpolation: " << intSample << endl;
   Int_t numSamps  = gr->GetN();
   cerr << "The number of points: " << numSamps << endl;
   for(int i=0;i<maxSamps;i++) {
      if(i<numSamps) {
	 gr->GetPoint(i,xVals,yVals);
	 newX[i]=xVals;
	 newY[i]=yVals*FFTtools::bartlettWindow(i,numSamps);
      }
      else {
	 newX[i]=newX[i-1]+intSample;
	 newY[i]=0;
      }
  }
   TGraph *grNew = new TGraph(maxSamps,newX,newY);
   return grNew;
}


/* Make sure the wfs all start with t=beginTime. Zero-pad wf from t=beginTime to the 1st bin of the original wf */
TGraph *evProcessTools::getWindowedAndPaddedEqualBeginGraph(TGraph *gr, const int maxSample, const double beginTime)
{

   //   static AraGeomTool *fGeomTool = AraGeomTool::Instance();
   if(!gr) return NULL;
   Double_t newX[maxSample],newY[maxSample];
   Int_t maxSamps=maxSample;

   Double_t time1=0;
   Double_t time2=0;
   Double_t xVals;
   Double_t yVals;

   gr->GetPoint(0,time1,yVals);
   gr->GetPoint(1,time2,yVals);

   Double_t intSample = time2 - time1;
   //cerr << "The interpolation: " << intSample << endl;
   Int_t numSamps  = gr->GetN();
   //cerr << "The number of points: " << numSamps << endl;

   Int_t time1_bin = (time1-beginTime)/intSample; //Number of bins from t=beginTime to t=time1

   for(int i=0;i<maxSamps;i++) {

      if(i < time1_bin ){
         newX[i] = time1 - (time1_bin-i)*intSample;
         newY[i] = 0;
      }

      else if( i<(time1_bin + numSamps)) {
	 gr->GetPoint(i-time1_bin,xVals,yVals);
	 newX[i]=xVals;
         /* if want to use Bartlett window */
	 //newY[i]=yVals*FFTtools::bartlettWindow(i-time1_bin,numSamps);
         /* if want to use modified Hann window */
         int modFrac = 4;
         newY[i]=yVals*evProcessTools::modifiedHannWindow(i-time1_bin, numSamps, modFrac);
      }
      else {
	 newX[i]=newX[i-1]+intSample;
	 newY[i]=0;
      }
  }
   TGraph *grNew = new TGraph(maxSamps,newX,newY);
   return grNew;
}

TGraph *evProcessTools::getScaledGraph(TGraph *gr)
{
  double statsArray[2] = {0.};
  double sigma, mean;

  setMeanAndSigmaInNoMax(gr, statsArray);
  mean  = statsArray[0];
  sigma = statsArray[1];

  //if( sigma == 0 ) return gr; //No waveform

  int nSamp = gr->GetN();
  double times, volts;
  double newX[nSamp], newY[nSamp];

  for(int i=0; i<nSamp; i++){
    gr->GetPoint(i, times, volts);
    newX[i] = times;
    if( sigma==0 ) newY[i] = 0.;
    else newY[i] = (volts - mean) / sigma;
  }

  TGraph *grNew = new TGraph(nSamp, newX, newY);
  return grNew;

}

TGraph *evProcessTools::getNormalizedGraph(TGraph *gr)
{

   int nSamp = gr->GetN();
   double times, volts;
   double pwr = 0.;
   double newX[nSamp], newY[nSamp];

   for(int i=0; i<gr->GetN(); i++){

      gr->GetPoint(i, times, volts);
      pwr += (volts*volts);
      newX[i] = times;
      newY[i] = volts;

   }

   pwr = sqrt(pwr);

   for(int i=0; i<nSamp; i++) newY[i] /= pwr;

   TGraph *grNew = new TGraph(nSamp, newX, newY);
   return grNew;

}

TGraph *evProcessTools::getRandomVoltageGraph(const int nSamp, const double wInt, const double range, TRandom3 *rnd){

  //TRandom3 *rnd = new TRandom3();
  //rnd->SetSeed(seed);
  double x[nSamp], y[nSamp];

  for(int i=0; i<nSamp; i++){
    y[i] = (2.*rnd->Rndm() - 1) * range;
    if(nSamp==0)cout<<y[i]<<endl;
    x[i] = wInt*i;
  }

  TGraph *grRnd = new TGraph(nSamp, x, y);

  return grRnd;
}

TH1D *evProcessTools::makeffthisto(TGraph *gr, char* name, Double_t intSample)
{

	TGraph *gr_new = getfftfromgraph(gr, intSample, 2000);

	int bin_number = gr_new->GetN();
	double bin_min = 0;
	double bin_max = 0;
	double dummy_y = 0;
	double temp_x;
	double temp_y;

	gr_new->GetPoint(0, bin_min, dummy_y);
	gr_new->GetPoint(bin_number-1, bin_max, dummy_y);


	TH1D *h1 = new TH1D(name, name, bin_number, bin_min, bin_max);

	for(int i=0;i<bin_number;i++)
	{
		gr_new->GetPoint(i, temp_x, temp_y);
		h1->Fill(temp_x,temp_y);
	}
	delete gr_new;
	return h1;
}


double evProcessTools::getTotalPower(std::vector<TGraph*> gr, std::vector<std::vector<double> > ant_loc)
{

	TH1D *histo[16] = {0};
	char histName[20];
//	std::vector<TGraph *> gr;
	double TCNewEl = 0;
//   	for(int a=0;a<16;a++)
//   	{
//		if(a<8){
//   		   gr.push_back( FFTtools::getInterpolatedGraph(gr_in[a],0.4));
//		}
//		else{
//   		   gr.push_back( FFTtools::getInterpolatedGraph(gr_in[a],0.625));
//		}
//   	}
	double totalPower =1;
	for(int a=0;a<16;a++){
	sprintf(histName, "fftHisto%d", a);
	if(a<8)	histo[a] = makeffthisto(gr[a], histName, 0.4);
	else histo[a] = makeffthisto(gr[a], histName, 0.625);
	totalPower = totalPower*integratePower(histo[a]);
	}


	for(int i=0;i<16;i++){
		delete histo[i];
	}
	return totalPower;
}


std::vector < double > evProcessTools::vecSum(std::vector< double > v1, std::vector< double > v2, int factor)
{
  std::vector < double > sumV;
  if(v1.size()!=v2.size())return sumV;
  for(int i = 0; i<v1.size();i++)
  {
    sumV.push_back(v1[i] + factor*v2[i]);
  }
  return sumV;
}



std::vector < double > evProcessTools::vecCrossProduct3D(std::vector< double > v1, std::vector< double > v2)
{
  std::vector < double > productV;
  if(v1.size()!=v2.size())return productV;
  double tempP = 0;
  int epsilon = 0;
  for(int i = 0; i<3;i++)
  {
    for(int j=0;j<3;j++){
      for(int k=0;k<3;k++){
	 if((i==j)||(i==k)||(j==k))epsilon =0;
    	 else epsilon = ( (i-j)*(i-k)*(j-k)/2 );
	 tempP+=epsilon*v1[j]*v2[k];
      }
    }
    productV.push_back(tempP);
  }
  return productV;
}



double evProcessTools::vecScalarProduct(std::vector< double > v1, std::vector< double > v2)
{
  if(v1.size()!=v2.size())return NULL;
  double product = 0;;
  for(int i = 0; i<v1.size();i++)
  {
    product+=v1[i]*v2[i];
  }
  return product;
}





double evProcessTools::vecLength(std::vector< double > v1)
{
  double squareSum = 0;
  for(int i=0;i<v1.size();i++)
  {
    squareSum+=v1[i]*v1[i];
  }
  return TMath::Sqrt(squareSum);
}





std::vector< double > evProcessTools::getXYProjection(std::vector< double > v)
{
  std::vector< double > projV;
  return projV;
}






std::vector<double> evProcessTools::getCrudeReconstruction(vector<vector< double > > delayMap, std::vector<std::vector<double> > ant_loc)
{
  std::vector<double> recoVec;
  double theta = 0;
  double weightSum = 0;
  double distance = 0;
  double clustering[360] = {0};
  for(int i=0;i<28;i++){
    if(TMath::Abs(delayMap[i][0] - delayMap[i][1])==4.0)
    {
	distance = TMath::Sqrt( pow(ant_loc[(int)delayMap[i][0]][0] - ant_loc[(int)delayMap[i][1]][0], 2)
				+ pow(ant_loc[(int)delayMap[i][0]][1] - ant_loc[(int)delayMap[i][1]][1], 2)
				+ pow(ant_loc[(int)delayMap[i][0]][2] - ant_loc[(int)delayMap[i][1]][2], 2)   );
	if(TMath::Abs(delayMap[i][2])<distance/0.171){
	  theta+=delayMap[i][3]*delayMap[i][3]*TMath::ACos( 0.171/distance*(-delayMap[i][2]) );
	  weightSum+=delayMap[i][3]*delayMap[i][3];
//	  cerr << "Theta is: " << theta << " for channels: " << (int)delayMap[i][0] << "and channel: " <<(int)delayMap[i][1] << " with delay " << delayMap[i][2] << " and distance: " << distance << endl;
	}
    }
  }
  theta = theta/weightSum;
//  cerr << "Final Theta: " << theta << endl;
 double alpha = 0;
 std::vector < double> xVec;
 xVec.push_back(1.0);
 xVec.push_back(0.0);
 xVec.push_back(0.0);
 weightSum = 0;
 for(int j=0;j<4;j++)
 {
 int level = j;
  for(int i=0;i<28;i++){
    if(delayMap[i][0]<4*(level+1) && delayMap[i][1]<4*(level+1) && delayMap[i][0]>=4*(level) && delayMap[i][1]>=4*(level)  )
    {
	std::vector<double> pair = evProcessTools::vecSum(ant_loc[(int)delayMap[i][0]], ant_loc[(int)delayMap[i][1]], -1.0);
	distance = evProcessTools::vecLength( pair );
	if(TMath::Abs(delayMap[i][2])<distance/0.171){
	  if(pair[1]>0) alpha = TMath::ACos( evProcessTools::vecScalarProduct(pair,xVec)/distance );
	  else alpha = 2.0*TMath::Pi() - TMath::ACos( evProcessTools::vecScalarProduct(pair,xVec)/distance );
//	  cerr << "Alpha is: " << alpha << " for channels: " << (int)delayMap[i][0] << "and channel: " <<(int)delayMap[i][1] << endl;
//	  cerr << "Writing to: " <<  0.171/distance*(-delayMap[i][2])*TMath::Sin(theta) << " With theta: " << theta << endl;
//	  cerr << "Writing to: " << (TMath::Nint(180/TMath::Pi()*( alpha +  TMath::ACos( 0.171/distance*(-delayMap[i][2])*TMath::Sin(theta) ) )) + 360)%360 << " With weight: " << delayMap[i][3]*delayMap[i][3] << endl;
	  clustering[(TMath::Nint(180/TMath::Pi()*( alpha -  TMath::ACos( 0.171/distance*(-delayMap[i][2])*TMath::Sin(theta) ) ))+360)%360]+=delayMap[i][3]*delayMap[i][3];
	  clustering[(TMath::Nint(180/TMath::Pi()*( alpha +  TMath::ACos( 0.171/distance*(-delayMap[i][2])*TMath::Sin(theta) ) ))+360)%360]+=delayMap[i][3]*delayMap[i][3];
//	cerr << "Final angles are: " << (TMath::Nint(180/TMath::Pi()*( alpha -  TMath::ACos( 0.171/distance*(-delayMap[i][2])*TMath::Sin(theta) ) ))+360)%360 << " and " << (TMath::Nint(180/TMath::Pi()*( alpha + TMath::ACos( 0.171/distance*(-delayMap[i][2])*TMath::Sin(theta) ) ))+360)%360;
//	cerr << "Without INT angles are: " << (180/TMath::Pi()*( alpha -  TMath::ACos( 0.171/distance*(-delayMap[i][2])*TMath::Sin(theta) ) ))+360 << " and " << (180/TMath::Pi()*( alpha + TMath::ACos( 0.171/distance*(-delayMap[i][2])*TMath::Sin(theta) ) ))+360 << endl;
	}
    }
  }
 }
  double tempWeight = 0;
  double maxAngle = 0;
  double maxCluster = 0;
  double tempAngle = 0;
//  cerr << "This was ok!" << endl;
  for(int a=0;a<360;a++)
  {
     tempWeight=0;
     tempAngle = 0;
     for(int dd=0;dd<30;dd++){
	tempWeight+=clustering[(a+dd)%360];
	tempAngle+=(a+dd)%360*clustering[(a+dd)%360];
     }
//     cerr << "Here we get: " << tempWeight << " for an angle of: " << a << endl;
     if(tempWeight>maxCluster){ maxCluster=tempWeight; maxAngle=tempAngle/tempWeight;}
  }

  recoVec.push_back(theta);
  recoVec.push_back(maxAngle/180.0*TMath::Pi());
  return recoVec;
}





double evProcessTools::getCorehentlySummedWaveAmp(std::vector<TGraph *> gr, vector<vector< double > >delayMap, double intSample, std::vector<std::vector< double > > *cohSumWave )
{
  int iterator[28];
  int countSmaller = 0;
  TGraph *gr1 = 0;
  TGraph *grSum = 0;
  for(int i=0;i<28;i++){
    countSmaller=0;
    for(int j=0;j<28;j++){
	if(i!=j && delayMap[i][3]<delayMap[j][3]) countSmaller++;
    }
    iterator[countSmaller]=i;
  }

  std::vector<double> vecVolts;
  std::vector<double> vecTimes;


  double times2, volts2;
  double times, volts;
  gr[0]->GetPoint(0, times, volts);
  gr[0]->GetPoint(1, times2, volts2);
  double currentInt = times2-times;
  int refChan = (int)delayMap[iterator[0]][0];
  int binOffset = 0;
//  gr[refChan]->SetPoint(0,0.0,0.0);
  grSum = FFTtools::getInterpolatedGraph(gr[refChan],intSample);
   Double_t *xSum = grSum->GetX();
   Double_t *ySum = grSum->GetY();

  double waveMax = 0;
//  cout << "The ref channel is: " << refChan << "  With: " << grSum->GetN() << " points." << endl;
  double graphOffset = 0;
  for(int i=0;i<28;i++){
    if((int)delayMap[i][0]==refChan){
	for(int p=0;p<gr[(int)delayMap[i][1]]->GetN();p++){
	  gr[(int)delayMap[i][1]]->GetPoint(p,times,volts);
	  if(p==0)graphOffset = xSum[0] - times;
	  gr[(int)delayMap[i][1]]->SetPoint(p,times+delayMap[i][2],volts);
//	  if(p==0)gr[(int)delayMap[i][1]]->SetPoint(p,0.0 + delayMap[i][2],0.0);
	}
	gr1 = FFTtools::getInterpolatedGraph(gr[(int)delayMap[i][1]],intSample);
        binOffset = -TMath::Nint((delayMap[i][2] - graphOffset)/intSample);
//	cout << "Delays: " << delayMap[i][2] << " binoffset: " << binOffset << " With partner: Ch" <<(int)delayMap[i][1] << endl;
    }
    else if((int)delayMap[i][1]==refChan){
	for(int p=0;p<gr[(int)delayMap[i][0]]->GetN();p++){
	  gr[(int)delayMap[i][0]]->GetPoint(p,times,volts);
	  if(p==0)graphOffset = xSum[0] - times;
	  gr[(int)delayMap[i][0]]->SetPoint(p,times - delayMap[i][2],volts);
//	  if(p==0)gr[int(delayMap[i][0])]->SetPoint(p,0.0 - delayMap[i][2],0.0);
	}
	gr1 = FFTtools::getInterpolatedGraph(gr[(int)delayMap[i][0]],intSample);
        binOffset = TMath::Nint((delayMap[i][2] + graphOffset)/intSample);
//	cout << "Delays: " << delayMap[i][2] << " binoffset: " << binOffset << " With partner: Ch" <<(int)delayMap[i][0]  << endl;
    }
    else{ continue;}
    waveMax = 0;
    for(int g = 0;g<grSum->GetN();g++){
	  if(g<gr1->GetN()-binOffset && g+binOffset>=0){
	    gr1->GetPoint(g+binOffset, times, volts);
	    ySum[g]+=volts;
//	    if(g%200==0)cout << "Compare: " << xSum[g] << " and " << times << endl;
	}
	if(ySum[g]>waveMax)waveMax=ySum[g];
	if(i==28){vecVolts.push_back(ySum[g]);vecTimes.push_back(xSum[g]);}
    }
    delete gr1;
  }

  cohSumWave->push_back(vecTimes);
  cohSumWave->push_back(vecVolts);



//  cout << "The max waveAmp is: " << waveMax << endl;
  delete grSum;
  return waveMax;
}







//*******************************************************************************************************************//
//*** This module calculates the correlation between signal waveforms and writes vectors with************************//
//***  delay times and correlation amplitudes************************************************************************//
//*******************************************************************************************************************//
std::vector<double> evProcessTools::getMaximumCorrelationSum(
		std::vector<TGraph*> gr, 				//** The vector of input waveforms to be correlated
		std::vector<std::vector<double> > ant_loc, 		//** FIXME: This is not used anymore, can be deleted
		Double_t *corrSum, 					//** Array to store correlation sum, if needed
		std::vector<std::vector<double> > *delayMaps, 		//** Vector to store delay times and correlation amps.
		std::vector<std::vector< double > > *cohSumWave,	//** The coherently summed wave, just for testing
		int stationId, 						//** We need the station ID because of possibly different treatment of stations
		int plot						//** In case we want to plot correlations, we set this "1"
		)
{

	double maxNorm[2] ={0.0};
	double t1,v1;
	double t2,v2;
	double tg1,vg1;
	double tg2,vg2;
	double tg1a,vg1a;
	double tg2a,vg2a;
	double integral = 0;
	double graphIntegral1 = 0;
	double graphIntegral2 = 0;
	double finalDelay=0;
//	std::vector<TGraph *> gr;
	double tempA;
	TGraph *grCorr = 0;
	TGraph *grEnv = 0;
	TGraph *grEnv2 = 0;
	TGraph *grCorrEnv = 0;

	corrSum[0]=0;
	corrSum[1]=0;
	corrSum[2]=0;
	corrSum[3]=0;
	corrSum[4]=0;
	corrSum[5]=0;
	vector< double > delays;
	vector<vector< double > > delayMapV;
	vector<vector< double > > delayMapH;
	double meanGr = 0;
	double meanCorr = 0;
	std::vector< double > rangeGr;
	std::vector< double > rangeCorr;

	double peaksCorr[2] = {0};
	double peaksGr[2] = {0};
	double start1 = 0;
	double start2 = 0;
	double length1 = 0;
	double length2 = 0;
	double end1 = 0;
	double end2 = 0;
	double center = 0;
	double holdEnv1 =0;
	double holdEnv2 =0;
	int ovCount1 = 0;
	int ovCount2 = 0;
	double intSamp1 = 0;
	double intSamp2 = 0;
	for(int ch1=0;ch1<15;ch1++){
	  for(int ch2=ch1+1;ch2<16;ch2++){
		if(ch1<8)intSamp1=0.4;
		else intSamp1=0.625;
		if(ch2<8)intSamp2=0.4;
		else intSamp2=0.625;
		ovCount1 = 0;
		ovCount2 = 0;
		integral = 0;
		graphIntegral1 = 0;
		graphIntegral2 = 0;
		tempA=0;
		meanGr = 0;
		meanCorr = 0;
		double overlap = 0;
		int delayBin = 0;
		finalDelay = 500.0;
	      if(!(stationId==2 && ch2==15) && !(ch2==ch1+8) ){
		grCorrEnv = new TGraph();//FFTtools::getHilbertTransform(grCorr);
		grEnv = FFTtools::getHilbertEnvelope(gr[ch1]);
		grEnv2 = FFTtools::getHilbertEnvelope(gr[ch2]);
		grCorr = FFTtools::getInterpolatedCorrelationGraph(grEnv, grEnv2, 0.1);
		gr[ch1]->GetPoint(gr[ch1]->GetN()-1, tg1, vg1);
		gr[ch2]->GetPoint(gr[ch2]->GetN()-1, tg2, vg2);
		grCorr->GetPoint(grCorr->GetN()/2,t1,v1);
		center = t1;
		length1 = tg1;
		length2 = tg2;
		end1 = tg1;
		end2 = tg2;
		gr[ch1]->GetPoint(0,t1,v1);
		gr[ch2]->GetPoint(0,t2,v2);
//		cerr << "Offset, gr1: " << t1 << " gr2: " << t2 << " center: " << center << endl;
		double graphOffset = t1 - t2;
		holdEnv1 =0;
		holdEnv2 =0;
		int envLength = 20;
		int grLength = 0;
		double ovH1 = 0;
		double ovH2 = 0;
		if(gr[ch1]->GetN()>gr[ch2]->GetN() ) grLength = gr[ch1]->GetN();
		else grLength = gr[ch2]->GetN();
		for(int s=0;s<envLength;s++){
			if( s<gr[ch1]->GetN() ){gr[ch1]->GetPoint(s,t1,v1);}
			if( s<gr[ch2]->GetN() ){gr[ch2]->GetPoint(s,t2,v2);}
			if(s<envLength){
				holdEnv1+=v1*v1;
				holdEnv2+=v2*v2;
			}
		}
//		cout << "The offset is: " << graphOffset << endl;
		rangeGr.clear();
		rangeCorr.clear();
		for(int g=0;g<grCorr->GetN();g++)
		{
			if(g<gr[ch1]->GetN()){
			   gr[ch1]->GetPoint(g, tg1, vg1);
			   graphIntegral1+=vg1*vg1;
			   if(g>envLength-1){ gr[ch1]->GetPoint(g-envLength,tg1a,vg1a);  holdEnv1 = holdEnv1 - vg1a*vg1a + vg1*vg1;}
		//	   grEnv->SetPoint(g,tg1,TMath::Sqrt(holdEnv1));  //This replaces a hilbert envelope, which takes a lot of time to calculate!
			}
			if(g<gr[ch2]->GetN()){
			   gr[ch2]->GetPoint(g, tg2, vg2);
			   graphIntegral2+=vg2*vg2;
			   if(g>envLength-1){ gr[ch2]->GetPoint(g-envLength,tg2a,vg2a);  holdEnv2 = holdEnv2 - vg2a*vg2a + vg2*vg2;}
		//	   grEnv2->SetPoint(g,tg2,TMath::Sqrt(holdEnv2));  //This replaces a hilbert envelope, which takes a lot of time to calculate!
			}
			if(g==0){
				start1 = tg1;
				start2 = tg2;
				length1 = length1 - start1;
				length2 = length2 - start2;
			}
			grCorr->GetPoint(g,t1,v1);
			if(ovCount1<2*gr[ch1]->GetN() && ovCount2<2*gr[ch2]->GetN() ){
				if( (t1 - center)<=0.0 && (length1 + (t1 - center) )>intSamp1*ovCount1){gr[ch1]->GetPoint(ovCount1,tg1,vg1);ovH1+=vg1*vg1;ovCount1++;}
				if( (t1 - center)<=0.0 && (length2 + (t1 - center) )>intSamp2*ovCount2){gr[ch2]->GetPoint(gr[ch2]->GetN()-1-ovCount2,tg2,vg2);ovH2+=vg2*vg2;ovCount2++;}
				if( (t1 - center)>0.0 && (length1 + (t1 - center) )>intSamp1*ovCount1){gr[ch1]->GetPoint(ovCount1 - gr[ch1]->GetN(),tg1,vg1);ovH1-=vg1*vg1;ovCount1++;}
				if( (t1 - center)>0.0 && (length2 + (t1 - center) )>intSamp2*ovCount2){gr[ch2]->GetPoint(2*gr[ch2]->GetN()-1-ovCount2,tg2,vg2);ovH2-=vg2*vg2;ovCount2++;}
			}
			if(ovH1<0.0) ovH1=1.0;
			if(ovH2<0.0) ovH2=1.0;
			if(ovH1*ovH2<0.0) overlap=1.0;
			overlap = 1.0;// TMath::Sqrt(ovH1*ovH2);

			if(overlap>0.1 && /*TMath::Abs*/(v1)/TMath::Sqrt(overlap)>tempA ){finalDelay = t1;tempA=(v1)/TMath::Sqrt(overlap);delayBin=g;}
			integral+=TMath::Sqrt(v1*v1)/grCorr->GetN();
			if(overlap>0.1){grCorrEnv->SetPoint(g, t1, v1/TMath::Sqrt(overlap) );}
			else grCorrEnv->SetPoint(g, t1, 0.0 );

		}
//		cerr << "The first final delay: " << finalDelay << "   For ch: " << ch1 << "  " << ch2 << endl;
/*
		if(1){
//		cerr << "The single delay is: " << finalDelay << endl;
//		grCorr->GetPoint(grCorr->GetN()/2 -1000,t1,v1);
//		cerr << "Corr: " << t1 << " " << length1 << "  " << length2 << " start1: " << start1 << " start2: " << start2 << " overlap: " << abs(t1 - center) << endl;
//		cerr <<abs( t1 - (start1 - start2)) << endl;
		double waveMean = getMeanBaseLine(grEnv);
		double waveMean2 = getMeanBaseLine(grEnv2);
//		double corrMean = getMeanBaseLine(grCorrEnv);
		double wavePeakPos[2] = {0};
		double wavePeakAmp[2] = {0};
		double wavePeakPos2[2] = {0};
		double wavePeakAmp2[2] = {0};
		double corrPeakPos[3] = {0};
		double corrPeakAmp[3] ={0};
		int wPeakCount = 1;
		int wPeakCount2 = 1;
		int cPeakCount = 1;
		wPeakCount = peakFinder(grEnv, 2, 200, &wavePeakPos[0], &wavePeakAmp[0], 2.5*waveMean, 2);
		wPeakCount2 = peakFinder(grEnv2, 2, 200, &wavePeakPos2[0], &wavePeakAmp2[0], 2.5*waveMean2, 2);
		if((end1-start1)>300.0 && (end2-start2)>300.0 && ( (wavePeakPos[0]-start1)<15.0 || (wavePeakPos2[0]-start2)<15.0 || (end1 - wavePeakPos[0])<15.0 || (end2 - wavePeakPos2[0])<15.0 ) ){
			tempA = 0.0;
		}
		else
		{
		if(wPeakCount>1){
			if(wavePeakAmp[1]>wavePeakAmp[0]*2.0/4.0) wPeakCount = 2;
			else wPeakCount = 1;
		}
		if(wPeakCount2>1){
			if(wavePeakAmp2[1]>wavePeakAmp2[0]*2.0/4.0) wPeakCount2 = 2;
			else wPeakCount2 = 1;
		}
		if(cPeakCount>1){
			int holdPeak = 1;
			if(corrPeakAmp[1]>corrPeakAmp[0]*2.0/4.0) holdPeak+=1;
			if(corrPeakAmp[2]>corrPeakAmp[0]*2.0/4.0) holdPeak+=1;
			cPeakCount=holdPeak;
		}
//		cerr << "We found: " << wPeakCount << "  " << wPeakCount2 << " peaks in the wave and "	<< cPeakCount << " peaks in the correlation-graph." << endl;
		//Here we go through different situations:
		//Ex.: If the reference ch1 has one signal-peak: earlier hits on ch2 will give positive correlation delays.
		//We want the first hit, because the second one is the reflection.
		double tempWPos1, tempWPos2;
		if(wPeakCount==2 || wPeakCount2==2){
		  if(wPeakCount==1 && wPeakCount2==2){
//			cerr << "Double peak2: " << wavePeakPos2[0] << "  " << wavePeakPos2[1] << endl;
			tempWPos1 = wavePeakPos[0];
			if(wavePeakPos2[1]>wavePeakPos2[0]) tempWPos2 = wavePeakPos2[0];
			else tempWPos2 = wavePeakPos2[1];
			finalDelay = tempWPos1 - tempWPos2;
		  }
		  if(wPeakCount==2 && wPeakCount2==1){
//			cerr << "Double peak1: " << wavePeakPos[0] << "  " << wavePeakPos[1] << endl;
			if(wavePeakPos[1]>wavePeakPos[0]) tempWPos1 = wavePeakPos[0];
			else tempWPos1 = wavePeakPos[1];
			tempWPos2 = wavePeakPos2[0];
			finalDelay = tempWPos1 - tempWPos2;
		  }
		  if(wPeakCount==2 && wPeakCount2==2){
//			cerr << "Double peak1: " << wavePeakPos[0] << "  " << wavePeakPos[1] << endl;
//			cerr << "Double peak2: " << wavePeakPos2[0] << "  " << wavePeakPos2[1] << endl;
			if(wavePeakPos[1]>wavePeakPos[0]) tempWPos1 = wavePeakPos[0];
			else tempWPos1 = wavePeakPos[1];
			if(wavePeakPos2[1]>wavePeakPos2[0]) tempWPos2 = wavePeakPos2[0];
			else tempWPos2 = wavePeakPos2[1];
				finalDelay = tempWPos1 - tempWPos2;
//			cerr << tempWPos1 << "  " << tempWPos2 << "   " << finalDelay << endl;
		  }
		  int pp = int( (finalDelay - center)*10.0) + grCorr->GetN()/2;
		  tempA=0;
		  for(int g=pp-200;g<pp+200;g++)
		  {
			grCorrEnv->GetPoint(g,t1,v1);
			if((v1)>tempA){finalDelay = t1;tempA=TMath::Abs(v1);}
		  }
		}

		}
	}//added for exclusion! FIXME
*/
		if(plot==1){
		  TGraph *grShift = new TGraph();
		  for(int ff=0;ff<gr[ch1]->GetN();ff++){
			grEnv->GetPoint(ff,t1,v1);
			grShift->SetPoint(ff, t1 - finalDelay, v1);
		  }
		  char title[40];
		  sprintf(title, "Correlation: Ch%d, Ch%d", ch1, ch2);
		  char cname[40];
		  sprintf(cname, "CorrelationCh%dCh%d", ch1, ch2);
		  char g1title[40];
		  sprintf(g1title, "Wave Ch%d", ch1);
		  char g2title[40];
		  sprintf(g2title, "Wave Ch%d", ch2);
		  TCanvas *c2 = new TCanvas(cname);
		  c2->Divide(1,4);
		  c2->cd(1);
			grEnv->SetTitle(g1title);
			grEnv->Draw("ALP");
		  c2->cd(2);
			grEnv2->SetTitle(g2title);
			grEnv2->Draw("ALP");
		  c2->cd(3);
//		  	grCorrEnv->SetTitle(title);
//		  	grCorr->Draw("ALP");
		  	grCorr->SetLineColor(2);
		  	grCorr->Draw("ALP");
		  c2->cd(4);
			grEnv2->Draw("ALP");
			grShift->SetLineColor(2);
			grShift->Draw("same");
		  TFile *ff = new TFile("corrFile.root", "UPDATE");
		  c2->Write();
		  ff->Close();
		  delete grShift;
		  delete ff;
		  delete c2;
		}
//		cerr << "Final delay after new voodoo algorithm: "<< wPeakCount << " and " << wPeakCount2 << " delay: " << finalDelay << endl;
		delete grCorr;
		delete grEnv;
		delete grCorrEnv;
		delete grEnv2;
	      }
		double sampFactor = 0;
		if(ch1<8 && ch2<8) sampFactor = 0.4/0.625;
		else if(ch1<8 || ch2<8) sampFactor = TMath::Sqrt(0.4/0.625);
		else sampFactor = 1.0;
		if( ((stationId==2) && (ch1==15 || ch2==15 ))||(ch2==ch1+8) ){//Excluding bad channels from ARA02 and Hpol-Vpol pairs right next to each other
		delays.push_back(ch1);
		delays.push_back(ch2);
		delays.push_back(0.0);
		delays.push_back(0.0);
		delays.push_back(0.0);
		}
		else{
		delays.push_back(ch1);
		delays.push_back(ch2);
		delays.push_back(finalDelay);
		delays.push_back(tempA/integral);
		delays.push_back(tempA*tempA/TMath::Sqrt((graphIntegral1)*(graphIntegral2)));
		}


		if(ch2<8){
		  delayMapV.push_back(delays);
		  corrSum[0]+=tempA/integral;
		  corrSum[1]+=tempA*tempA/TMath::Sqrt(graphIntegral1*graphIntegral2);
		  corrSum[4]+=tempA;
		}
		if(ch1>7){
		  delayMapH.push_back(delays);
		  corrSum[2]+=tempA/integral;
		  corrSum[3]+=tempA*tempA/TMath::Sqrt(graphIntegral1*graphIntegral2);
		  corrSum[5]+=tempA;
		}
		delayMaps->push_back(delays);
		delays.clear();
		if(ch1==0 && ch2==1)corrSum[6] = grCorr->GetN();
		if(ch1==8 && ch2==9)corrSum[7] = grCorr->GetN();
	  }
	}

	vector <double> VpolReco;// = getCrudeReconstruction(delayMapV, ant_loc);
//	double waveAmpV = 0;
	double waveAmpV = 0;//getCorehentlySummedWaveAmp(gr, delayMapV, 0.1, cohSumWave);
	VpolReco.push_back(waveAmpV);
	delayMapV.clear();

	vector <double> HpolReco;// = getCrudeReconstruction(delayMapH, ant_loc);
	double waveAmpH = 0;// getCorehentlySummedWaveAmp(gr, delayMapH, 0.1, cohSumWave);
	HpolReco.push_back(waveAmpH);
	delayMapH.clear();

	if(waveAmpH>waveAmpV)return HpolReco;
	else return VpolReco;
}










std::vector< vector< double> > evProcessTools::getFFT(std::vector<TGraph*> gr, std::vector<std::vector<double> > ant_loc, vector<vector< double > > * FFTpower, vector<vector<double> > * FFTfreq, vector<vector< double > > * normFFT,  vector<vector< double > > refArray, double *cwValueV, double *cwValueH, vector<vector< double > > * rayleighV, vector<vector<double> > * rayleighF)
{

	TH1D *histo[16] = {0};
	char histName[20];
//	std::vector<TGraph *> gr;
//   	for(int a=0;a<16;a++)
//   	{
//		if(a<8){
//   		   gr.push_back( FFTtools::getInterpolatedGraph(gr_in[a],0.4));// actually 0.4
//		}
//		else{
//   		   gr.push_back( FFTtools::getInterpolatedGraph(gr_in[a],0.625));// actually 0.625
//		}
//   	}
	int bin_number = 0;
	double temp_x;
	double temp_y;
	std::vector<double> freqs;
	std::vector<double> powers;
	std::vector<double> normPowers;
	std::vector<double> rayleighFreq;
	std::vector<double> rayleighVolts;
	std::vector<double> cwSingleVecIn;
	std::vector<vector<double> > cwSingleCh;

	TGraph * gr_new = 0;
	int rebinCount =0;
	double fftP = 0;
	double fftF = 0;
	int subBins = 0;
	double cwSingle[200] = {0};
	double cwMultV[200] = {1.0};
	double cwMultH[200] = {1.0};
	double peakBase = 0;
	double wInt = 0.0;
	for(int a=0;a<16;a++){
	  rebinCount=0;
	  fftP=0;
	  subBins=0;
//	if(gr[a]->GetN()<550){
//	  cerr << "Run FFT on ch" << a << endl;
	  if(a<8)wInt=0.4;
	  else wInt=0.625;
	  gr_new = getfftfromgraph(gr[a], wInt, 2000);
//	}
	  bin_number = gr_new->GetN();
	  for(int i=0;i<bin_number;i++)
	  {
		gr_new->GetPoint(i, temp_x, temp_y);
	    	if(temp_x<=1.0E3){
		  if(!(temp_x<(rebinCount+1)*5.0)){
		    fftF=rebinCount*5.0;
		    normPowers.push_back(fftP/subBins/refArray[a][rebinCount]);
		    freqs.push_back(fftF);
		    powers.push_back(fftP/subBins);
//		  cerr << "The number of subbins is: " << subBins << endl;
		    fftP=0;
		    subBins=0;
		    rebinCount++;
		    if(rebinCount==200) break;
		  }
	       	  if(temp_x>=rebinCount*5.0 && temp_x<(rebinCount+1)*5.0){fftP+=/*TMath::Sqrt*/(temp_y)/gr[a]->GetN()/*/wInt*/;subBins++;}
//		 if(temp_x>249 && temp_x<251) cerr << temp_x << "   " << temp_y << endl;
		}
		rayleighFreq.push_back(temp_x);
		rayleighVolts.push_back(TMath::Sqrt(temp_y)/wInt);
	  }//edn loop bin_number
//	  cerr << "Rebinning done!" << endl;
	  cwSingleVecIn.push_back(0.0);
	  cwSingleVecIn.push_back(0.0);
	  for(int i=2;i<198;i++){
	    peakBase = (normPowers[i-2] + normPowers[i-1] + normPowers[i+1] + normPowers[i+2])/4.0;
	    cwSingle[i] = normPowers[i];
	    cwSingleVecIn.push_back(cwSingle[i]);
	    if(a<8){
		if(a==2){
//		    if(cwSingle[i]>1.0)
			cwMultV[i]=cwSingle[i];
//		    else cwMultV[i]=1.0;
		}
	//	else /*if(cwSingle[i]>1.0)*/ cwMultV[i] = cwMultV[i]*cwSingle[i];
	    }
	    else{
		if(a==8){
//		    if(cwSingle[i]>1.0)
			cwMultH[i]=cwSingle[i];
//		    else cwMultH[i]=1.0;
		}
	//	else /*if(cwSingle[i]>1.0)*/ cwMultH[i] = cwMultH[i]*cwSingle[i];
	    }
	  }
	  cwSingleVecIn.push_back(0.0);
	  cwSingleVecIn.push_back(0.0);
	  delete gr_new;
//	  cerr << "Is this the problem, with size: " << cwSingleVecIn.size() << endl;
	  cwSingleCh.push_back(cwSingleVecIn);
//	  cerr << "Not anymore!" << endl;
	  FFTpower->push_back(powers);
	  FFTfreq->push_back(freqs);
	  normFFT->push_back(normPowers);
	  rayleighF->push_back(rayleighFreq);
	  rayleighV->push_back(rayleighVolts);
	  rayleighVolts.clear();
	  rayleighFreq.clear();
	  normPowers.clear();
	  freqs.clear();
	  powers.clear();
	  cwSingleVecIn.clear();
	}
//	cerr << "CW-Max found!" << endl;
	double cwMax = 0;
	for(int i=0;i<200;i++){
	  cwValueV[i]=cwMultV[i];
	  cwValueH[i]=cwMultH[i];
	  if(i*5.0>100.0 && i*5.0<850.0)
	  {
		if(cwValueV[i]>cwMax)cwMax = cwValueV[i];
		if(cwValueH[i]>cwMax)cwMax = cwValueH[i];
	  }
	}

//	for(int i=0;i<16;i++){
//		delete gr[i];
//	}
	return cwSingleCh;
}






TGraph *evProcessTools::cwCleaner(TGraph *grWaveIn, Double_t minFreq, Double_t maxFreq, double reduction)
{

//cerr << "Entered the cleaning module, with maxFreq: " << maxFreq << " and minFreq: " << minFreq << endl;
    TGraph *grWave = FFTtools::padWaveToLength(grWaveIn, 2000);
    int numPoints = grWaveIn->GetN();
//    TGraph *grWave = evProcessTools::getBartlettAndPaddedGraph(grWaveIn, 2000);
//    cerr << "Number of dataPoints: " << grWave->GetN() << endl;
//    TGraph *grWave =
    double *oldY = grWave->GetY();
    double *oldX = grWave->GetX();
    double deltaT=oldX[1]-oldX[0];
    int length=grWave->GetN();
    FFTWComplex *theFFT=FFTtools::doFFT(length,oldY);

    int newLength=(length/2)+1;

    //    double fMax = 1/(2*deltaT);  // In Hz
    double deltaF=1/(deltaT*length); //Hz
    deltaF*=1e3; //MHz
    //    std::cout << fMax << "\t" << deltaF << "\t" << deltaT << "\t" << length << std::endl;

    double tempF=0;
    double ratio[10];
    double sum = 0;
    int count=0;
//cerr << "Produced the FFT!" << endl;
    for(int i=0;i<newLength;i++) {
      //      std::cout << tempF << "\t" << theFFT[i].re << "\t" << theFFT[i].im << "\n";

      if(tempF>minFreq && tempF<maxFreq) {
	sum+= theFFT[i].re*theFFT[i].re + theFFT[i].im*theFFT[i].im;
	count++;
//            std::cout << tempF << "\t" << theFFT[i].re << "\t" << theFFT[i].im << "\n";
      }
      tempF+=deltaF;
    }
    sum = sum/count;
  //  cerr << "Calculated the normSum: " << sum << endl;
    tempF=0;
    count=0;
    for(int i=0;i<newLength;i++) {
      //      std::cout << tempF << "\t" << theFFT[i].re << "\t" << theFFT[i].im << "\n";

      if(tempF>minFreq && tempF<maxFreq) {
	ratio[count] = (theFFT[i].re*theFFT[i].re + theFFT[i].im*theFFT[i].im)/sum;
	theFFT[i].re=1.0/TMath::Sqrt(ratio[count]*reduction)*theFFT[i].re;
	theFFT[i].im=1.0/TMath::Sqrt(ratio[count]*reduction)*theFFT[i].im;
//	std::cout << "The ratio: " << ratio[count] << endl;
	count++;
      }
      //      std::cout << tempF << "\t" << theFFT[i].re << "\t" << theFFT[i].im << "\n";
      tempF+=deltaF;
    }
//cerr << "Did the cleaning!" << endl;
    double *filteredVals = FFTtools::doInvFFT(length,theFFT);
//cerr << "performed the back transform!" << endl;
   int presamps = length - numPoints;
   TGraph *grFiltered = new TGraph();
   for(int i=0;i<numPoints;i++){
	grFiltered->SetPoint(i, oldX[i+presamps/2], filteredVals[i+presamps/2]);
//    filteredVals[i-presamps/2] = filteredVals[i];///FFTtools::bartlettWindow(i,numPoints);
//    oldX[i-presamps/2] = oldX[i];
  }


//    TGraph *grFiltered = new TGraph(numPoints,oldX,filteredVals);
//    TGraph *grFiltered = new TGraph(length,oldX,filteredVals);
//    cerr << "Number of dataPoints in filtered graph: " << grFiltered->GetN() << endl;
//cerr << "The new graph has: " << grFiltered->GetN() << " datapoints. " << endl;
    delete [] theFFT;
    delete [] filteredVals;
    return grFiltered;
}

















TGraph *evProcessTools::universalCleaner(TGraph *grWaveIn, std::vector< double> cwSingleVec )
{

//cerr << "Entered the cleaning module, with maxFreq: " << maxFreq << " and minFreq: " << minFreq << endl;

    int cwFreq=0;
    TGraph *grWave = FFTtools::padWaveToLength(grWaveIn, 2000);
    int numPoints = grWaveIn->GetN();


//    TGraph *grWave = evProcessTools::getBartlettAndPaddedGraph(grWaveIn, 2000);
//    cerr << "Number of dataPoints: " << grWave->GetN() << endl;
//    TGraph *grWave =
    double *oldY = grWave->GetY();
    double *oldX = grWave->GetX();
    double deltaT=oldX[1]-oldX[0];
    int length=grWave->GetN();
    FFTWComplex *theFFT=FFTtools::doFFT(length,oldY);

    int newLength=(length/2)+1;

    //    double fMax = 1/(2*deltaT);  // In Hz
    double deltaF=1/(deltaT*length); //Hz
    deltaF*=1e3; //MHz
    //    std::cout << fMax << "\t" << deltaF << "\t" << deltaT << "\t" << length << std::endl;
    double tempF=0;
for(int dd=0;dd<200;dd++){
    cwFreq=dd;
    Double_t minFreq = cwFreq*5.0;
    Double_t maxFreq =(cwFreq + 1.0)*5.0;
    double reduction = cwSingleVec[cwFreq];
    if(reduction>2.0){
    tempF=0;
    double ratio[10];
    double sum = 0;
    int count=0;
//cerr << "Produced the FFT!" << endl;
    for(int i=0;i<newLength;i++) {
      //      std::cout << tempF << "\t" << theFFT[i].re << "\t" << theFFT[i].im << "\n";

      if(tempF>minFreq && tempF<maxFreq) {
	sum+= theFFT[i].re*theFFT[i].re + theFFT[i].im*theFFT[i].im;
	count++;
//            std::cout << tempF << "\t" << theFFT[i].re << "\t" << theFFT[i].im << "\n";
      }
      tempF+=deltaF;
    }
    sum = sum/count;
  //  cerr << "Calculated the normSum: " << sum << endl;
    tempF=0;
    count=0;
    for(int i=0;i<newLength;i++) {
      //      std::cout << tempF << "\t" << theFFT[i].re << "\t" << theFFT[i].im << "\n";

      if(tempF>minFreq && tempF<maxFreq) {
	ratio[count] = (theFFT[i].re*theFFT[i].re + theFFT[i].im*theFFT[i].im)/sum;
	theFFT[i].re=1.0/TMath::Sqrt(ratio[count]*reduction)*theFFT[i].re;
	theFFT[i].im=1.0/TMath::Sqrt(ratio[count]*reduction)*theFFT[i].im;
//	std::cout << "The ratio: " << ratio[count] << endl;
	count++;
      }
      //      std::cout << tempF << "\t" << theFFT[i].re << "\t" << theFFT[i].im << "\n";
      tempF+=deltaF;
    }
   }//if reduction
}//loop dd
//cerr << "Did the cleaning!" << endl;
    for(int i=0;i<newLength;i++) {
      if(tempF<150.0 && tempF>850.0) {
	theFFT[i].re=0.0;
	theFFT[i].im=0.0;
      }
      tempF+=deltaF;
    }
    double *filteredVals = FFTtools::doInvFFT(length,theFFT);
//cerr << "performed the back transform!" << endl;
   int presamps = length - numPoints;
   TGraph *grFiltered = new TGraph();
   for(int i=0;i<numPoints;i++){
	grFiltered->SetPoint(i, oldX[i+presamps/2], filteredVals[i+presamps/2]);
//    filteredVals[i-presamps/2] = filteredVals[i];///FFTtools::bartlettWindow(i,numPoints);
//    oldX[i-presamps/2] = oldX[i];
  }


//    TGraph *grFiltered = new TGraph(numPoints,oldX,filteredVals);
//    TGraph *grFiltered = new TGraph(length,oldX,filteredVals);
//    cerr << "Number of dataPoints in filtered graph: " << grFiltered->GetN() << endl;
//cerr << "The new graph has: " << grFiltered->GetN() << " datapoints. " << endl;
    delete [] theFFT;
    delete [] filteredVals;
    return grFiltered;
}


















double evProcessTools::getSimpleFFT(std::vector<TGraph*> gr, vector<vector< double > > * FFTpower, vector<vector<double> > * FFTfreq)
{

	TH1D *histo[16] = {0};
	char histName[20];
//	std::vector<TGraph *> gr;
//   	for(int a=0;a<16;a++)
//  	{
//		if(a<8){
//   		   gr.push_back( FFTtools::getInterpolatedGraph(gr_in[a],0.4));// actually 0.4
//		}
//		else{
//   		   gr.push_back( FFTtools::getInterpolatedGraph(gr_in[a],0.625));// actually 0.625
//		}
//   	}
	int bin_number = 0;
	double temp_x;
	double temp_y;
	std::vector<double> freqs;
	std::vector<double> powers;

	TGraph * gr_new = 0;
	int rebinCount =0;
	double fftP = 0;
	double fftF = 0;
	int subBins = 0;
	double peakBase = 0;
	double wInt = 0.0;
	for(int a=0;a<16;a++){
	  rebinCount=0;
	  fftP=0;
	  subBins=0;
//	if(gr[a]->GetN()<550){
	  if(a<8)wInt=0.4;
	  else wInt=0.625;
	  gr_new = getfftfromgraph(gr[a], wInt, 2000);
//	}
	  bin_number = gr_new->GetN();
	  for(int i=0;i<bin_number;i++)
	  {
		gr_new->GetPoint(i, temp_x, temp_y);
	    	if(temp_x<=1.0E3){
		  if(!(temp_x<(rebinCount+1)*5.0)){
		    fftF=rebinCount*5.0;
		    freqs.push_back(fftF);
		    powers.push_back(fftP/subBins);
//		  cerr << "The number of subbins is: " << subBins << endl;
		    fftP=0;
		    subBins=0;
		    rebinCount++;
		    if(rebinCount==200) break;
		  }
	       	  if(temp_x>=rebinCount*5.0 && temp_x<(rebinCount+1)*5.0){fftP+=/*TMath::Sqrt*/(temp_y)/gr[a]->GetN()/*/wInt*/;subBins++;}
//		 if(temp_x>249 && temp_x<251) cerr << temp_x << "   " << temp_y << endl;
		}
	  }//edn loop bin_number
//	  cerr << "Rebinning done!" << endl;
	  delete gr_new;
	  FFTpower->push_back(powers);
	  FFTfreq->push_back(freqs);
	  freqs.clear();
	  powers.clear();
	}
//	cerr << "CW-Max found!" << endl;
	for(int i=0;i<16;i++){
		delete gr[i];
	}
	return 0;
}


//This calculates the mean value of an envelope function. To account for peaks (to exclude them!), the mean waveform is divided into four parts
//and the two minimum values are taken to calculate the mean.
double evProcessTools::getMeanBaseLine(TGraph *gr)
{
	double mean[4] = {0};
	double t1, v1;
	for(int i=0;i<gr->GetN()/5;i++){
		gr->GetPoint(i*5,t1,v1);
		if(i*5<1.0/4.0*gr->GetN()) mean[0]+=v1/( 1.0/4.0*gr->GetN() );
		if(i*5>=1.0/4.0*gr->GetN() && i*5<2.0/4.0*gr->GetN()) mean[1]+=v1/( 1.0/4.0*gr->GetN() );
		if(i*5>=2.0/4.0*gr->GetN() && i*5<3.0/4.0*gr->GetN()) mean[2]+=v1/( 1.0/4.0*gr->GetN() );
		if(i*5>=3.0/4.0*gr->GetN()) mean[3]+=v1/( 1.0/4.0*gr->GetN() );
	}
//This would help to find the two lowest mean values and to exclude signal. Unfortunately it is problematic for the correlation graph.
	double temp1 = mean[0];
	double temp2 = 300.0;
	for(int i=0;i<4;i++){
//		cerr << "Temp-mean: " << mean[i] << endl;
		if(mean[i]<temp2){
			if(mean[i]<temp1){
				temp2 = temp1;
				temp1 = mean[i];
			}
			else{
				temp2 = mean[i];
			}
		}
	}
//	cerr << temp1 << "  " << temp2 << endl;
	return (mean[0]+mean[1]+mean[2]+mean[3])/4.0*5.0;
}

//This peakfinder can be used for envelope functions, which are constructed with the Hilbert transform.
//gr = Envelope graph. Generated through hilbert transform.
//numberOfpeaks = number of peaks we want starting from the strongest.
//range = specified range in which we search for one maximum. The full waveform is scanned in steps of half the range.
//peakPos = an array-address to store the peak positions
//peakHight = an array-address to store the peak hight. Not sure what I can use this for.
//threshold = a minimum threshold for the peaks hight to cross, to be considered as a peak.
int evProcessTools::peakFinder(TGraph *gr, int numberOfPeaks, int range, double *peakPos, double *peakHight, double threshold, int precision)
{
	double absMax = 0;
	double absMaxPos = 0;
	double tempMax1 = 0;
	double tempMaxPos1 = 0;
	double tempMax2 = 0;
	double tempMaxPos2 = 0;
	double t1,v1;
	int peakIt1=0;
	int peakIt2=0;
	int it1 = 0;
	int it2 = 0;
	std::vector< double > allPeaks;
	std::vector< double > allPos;
	int hit1 = 0;
	int hit2 = 0;
	for(int i=0;i<gr->GetN(); i+=precision){
		gr->GetPoint(i,t1,v1);
		if(v1>absMax){
			absMax = v1;
			absMaxPos = t1;
		}
		if(v1>threshold){
			if(i>=range*it1 && i<range*(it1+1)){
				if(v1>tempMax1){
					tempMax1 = v1;
					tempMaxPos1 = t1;//t1;
					peakIt1 = i;
					hit1 = 1;
				}
			}
			if(i>=range*it2+range/2 && i<range*(it2+1)+range/2){
				if(v1>tempMax2){
					tempMax2 = v1;
					tempMaxPos2 = t1;//t1;
					peakIt2 = i;
					hit2 =1;
				}
			}

		}
		if(range*(it1+1)-1-i<=precision){
//			cerr << "Range1: " << i << "    " << range*(it1+1)-1 << endl;
			if(hit1==1 ){
				if(peakIt1>range*it1+precision && peakIt1<range*(it1+1)-1-precision){
					if(allPeaks.size()>0){
					if(allPos[allPeaks.size()-1]!=tempMaxPos1){
					allPeaks.push_back(tempMax1);
					allPos.push_back(tempMaxPos1);
					}
					}
					else{
					allPeaks.push_back(tempMax1);
					allPos.push_back(tempMaxPos1);
					}
//					cerr << "Is this happening?" << endl;
				}
			}
			it1++;
			tempMax1 = 0;
			hit1=0;
		}
		if(range*(it2+1)+range/2-1-i<=precision){
//			cerr << "Range2: " << i << "    " << range*(it2+1)-1 << endl;
			if(hit2==1 ){
				if(peakIt2>range*it2+range/2+precision && peakIt2<range*(it2+1)+range/2-1-precision){
					if(allPeaks.size()>0){
					if(allPos[allPeaks.size()-1]!=tempMaxPos2){
					allPeaks.push_back(tempMax2);
					allPos.push_back(tempMaxPos2);
					}
					}
					else{
					allPeaks.push_back(tempMax2);
					allPos.push_back(tempMaxPos2);
					}
//					cerr << "Is this happening?" << endl;
				}
			}
			it2++;
			tempMax2 = 0;
			hit2=0;
		}

	}
	int finalCount = 0;
	if(allPeaks.size()==0){
		peakPos[0] = absMaxPos;
		peakHight[0] = absMax;
		return 1;
	}
	else{
		  int iterator[100] = {0};
		  int countSmaller = 0;
		  for(int i=0;i<allPeaks.size();i++){
//		    cerr << "Show: " << allPeaks[i] << "  " << allPos[i] << endl;
		    countSmaller=0;
		    for(int j=0;j<allPeaks.size();j++){
			if(i!=j && allPeaks[i]<allPeaks[j]) countSmaller++;
		    }
		    iterator[countSmaller]=i;
		  }
		  for(int p=0;p<allPeaks.size();p++){
			  if(p<numberOfPeaks){
				  peakPos[p] = allPos[iterator[p]];
				  peakHight[p] = allPeaks[iterator[p]];
				  finalCount++;
			  }
		  }
		  return finalCount;
	}
}





int evProcessTools::spikingString(std::vector<TGraph*> gr, int stationId)
{
	double graphMax[16] = {0};
	double stringAverage[4];
	double stringGradient[4];
	double totalAverage=0;
	double position=0;
	double tv, vv;
//	cerr << "the spikes: ";
	for(int i=0;i<16;i++)
	{
		graphMax[i]=0.0;
		for(int g=0;g<gr[i]->GetN();g++){
			gr[i]->GetPoint(g,tv,vv);
			if(TMath::Abs(vv)>graphMax[i]){graphMax[i] = TMath::Abs(vv);position=tv;}
		}
		if(stationId==2 && i%4==3)stringAverage[i%4]+=graphMax[i]/3.0;
		else stringAverage[i%4]+=graphMax[i]/4.0;
		totalAverage+=graphMax[i]/16.0;
//		cerr << graphMax[i] << " at " << position << "     ";
	}
//	cerr << endl;
	int stringSpike = 0;
	double tempMax = 0.0;
	for(int j=0;j<4;j++){
//		for(int k=0;k<4;k++){
//			if(stringAverage[j]>4.0/3.0*(stringAverage[(j+1)%4] +stringAverage[(j+2)%4] + stringAverage[(j+3)%4])){cerr << "Discrepancy: " << stringAverage[j] << " " << 1.0/3.0*(stringAverage[(j+1)%4] +stringAverage[(j+2)%4] + stringAverage[(j+3)%4]) << endl; stringSpike++;}
//		}
			if(100.0*stringAverage[j]*3.0/(stringAverage[(j+1)%4] +stringAverage[(j+2)%4] + stringAverage[(j+3)%4])>tempMax) tempMax = 100.0*stringAverage[j]*3.0/(stringAverage[(j+1)%4] +stringAverage[(j+2)%4] + stringAverage[(j+3)%4]);
	}
	stringSpike = tempMax;
	return stringSpike;
}

/*
 * Windowing function. Increases to 1.0 after numSample/modFrac
 * Standard Hann: modFrac = 2; LOPES: modFrac = 4
 */
double evProcessTools::modifiedHannWindow(int idx, int numSample, int modFrac){

   int winLen = 2*(numSample/modFrac);
   if( idx < numSample/modFrac )                   return 0.5*(1-cos(2*M_PI*idx/(double)(winLen-1)));
   else if ( idx > (modFrac-1)*numSample/modFrac)  return 0.5*(1-cos(2*M_PI*(winLen+idx-numSample)/(double)(winLen-1)));
   else                                            return 1.;

}
