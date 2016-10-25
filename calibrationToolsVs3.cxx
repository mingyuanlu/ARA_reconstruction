
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////  calibrateTime.cxx 
////      CUse previously determind calibration parameters, to calibrate the sample timing of a given channel/event.
////
////    Sep. 2013, tmeures@ulb.ac.be
////////////////////////////////////////////////////////////////////////////////

//Includes
#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <fstream>

using namespace std;

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
#include "TSystem.h"
#include "TMath.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"


double calConst[4][6][512][64][9] = {{{{{0}}}}};
double calConstErr[4][6][512][64][9] = {{{{{0}}}}};
double calConstHigh[4][6][512][64][5] = {{{{{0}}}}};
double baseCalConst[128][13] = {{0}};
double cableDelays[4][6] = {{0}};


//string calFilesDir = "/home/meures/AraCalibration/testFolder/calFiles";
string calFilesDir = "/home/meures/analysis/analysisCode/calModules/calFiles";


















TGraph *calibrateTimeS(int elchan, UsefulAtriStationEvent * owncalEvPtr, double *tempCorr, int stationId)
{
	double samplecorr[128] = {0};
	double sampleRms[128] = {0};
	double meanV1[40];
	double freqBevSev, freqBevSodd, freqBoddSev, freqBoddSodd, shift_odd, shift_ev, shift_evodd;
	TGraph *grOwnCal = owncalEvPtr->getGraphFromElecChan(elchan);
	int Nsamples = grOwnCal->GetN();

	ifstream calRead;
	char calReadchar[60];
	sprintf(calReadchar,(calFilesDir + "/calFileVs2ARA%2.2dCh%dSamples0.txt").c_str(), stationId, elchan);
	calRead.open(calReadchar);
	calRead >> shift_evodd >> shift_ev >> freqBevSev >> freqBoddSev;
	for(int d=0;d<64;d++)
	{
		calRead >> samplecorr[d*2] >> sampleRms[d*2];
	}
	calRead.close();
	sprintf(calReadchar,(calFilesDir+ "/calFileVs2ARA%2.2dCh%dSamples1.txt").c_str(), stationId, elchan);
	calRead.open(calReadchar);
	calRead >> shift_evodd >> shift_odd >> freqBevSodd >> freqBoddSodd;
	for(int d=0;d<64;d++)
	{
		calRead >> samplecorr[d*2+1] >> sampleRms[d*2+1];
	}
	calRead.close();

	double grtime, grvolt;
	double corrTime = 0;

	double voltCalValues[128] ={0};
	ifstream in;
	char baseVCF[200];
	sprintf(baseVCF, (calFilesDir + "/baseVoltCorrARA%2.2dCH%d.txt").c_str(),stationId, elchan);
	in.open(baseVCF);
	if(in.good())
	{
		for(int i=0;i<128;i++)
		{
//			in >> voltCalValues[i];
			voltCalValues[i] = 0;
		}
		in.close();
	}
	else{
		in.close();
		for(int i=0;i<128;i++)
		{
			voltCalValues[i] = 0;
		}
	}
			
//	for(int s=0;s<Nsamples;s++)
//	{
//		if(s%64==0)meanV1[s/64]=0;
//		grOwnCal->GetPoint(s,grtime,grvolt);
//		meanV1[s/128]+=grvolt/128.0;
//	}

	double tempt, tempv;
	double sort_times[2048];
	double sort_volts[2048];
	double sort_voltCalValues[2048];
	double sortRms[2048];
	double pointRms = 0;
	double block_time = 0.0;
	double blockTarray[40][64] = {{0}};
	double blockVarray[40][64] = {{0}};
	int block_count = 0;
	int grBlockN  = 0;
	for(int j= 0;j<Nsamples;j++)
	{
		grOwnCal->GetPoint(j,grtime,grvolt);
		grBlockN = owncalEvPtr->blockVec[j/16].getBlock();
		blockTarray[block_count][(j)%64]=grtime/3.2;
		blockVarray[block_count][(j)%64]=grvolt;
		if((j)%64==63)block_count++;
	}

	double temp_max =  1.0E22;
	double temp_values[3];
	double max_values[3];

	int grBlock0 = 0;
	double sampleTime[128] = {0};
	for(int s=0;s<128;s++)
	{
		if(s<2){
			sampleTime[0] = 0.0;
			sampleTime[1] = 1.0/3.2 + shift_evodd;
		}
		else{
			sampleTime[s]= sampleTime[s-2] + (1.0/1.6*samplecorr[s])*tempCorr[s%2];
		}
	}

	for(int j= 0;j<Nsamples;j++)
	{
		grBlockN = owncalEvPtr->blockVec[j/16].getBlock();
		if(j==0)grBlock0 = owncalEvPtr->blockVec[0].getBlock();

		if(grBlock0%2==0)
		{
			if(j%2==0){
				 corrTime = (j/128)*40.0 + sampleTime[j%128];
				 pointRms = sampleRms[j%128];
				 sort_voltCalValues[j] = voltCalValues[j%128];
			}
			else{
				 corrTime = (j/128)*40.0 + sampleTime[j%128]; 
				 pointRms = sampleRms[j%128];
				 sort_voltCalValues[j] = voltCalValues[j%128];
			}
		}
		else
		{
			if(j%2==0){
				 corrTime =  ((j+64)/128)*40.0 - 20.0 + sampleTime[(j+64)%128];
//				 corrTime =  ((j+64)/128)*40.0 + sampleTime[(j+64)%128] - sampleTime[64];
				 pointRms = sampleRms[(j+64)%128];
				 sort_voltCalValues[j] = voltCalValues[(j+64)%128];
			}
			else{
				 corrTime = ((j+64)/128)*40.0 - 20.0 + sampleTime[(j+64)%128]; 
//				 corrTime = ((j+64)/128)*40.0 + sampleTime[(j+64)%128] - sampleTime[64]; 
				 pointRms = sampleRms[(j+64)%128];
				 sort_voltCalValues[j] = voltCalValues[(j+64)%128];
			}
		}



		sort_times[j] = corrTime - cableDelays[elchan/8][elchan%8];
		sort_volts[j] = blockVarray[j/64][j%64];
		sortRms[j] =pointRms;
	}
	
	TGraph *NgrOwnCal = new TGraph(Nsamples);
	for(int go=0;go<Nsamples;go++)
	{
		NgrOwnCal->SetPoint(go,sort_times[go],sort_volts[go] + sort_voltCalValues[go]);
//		NgrOwnCal->SetPointError(go, sortRms[go], 0.0);
	}
	return NgrOwnCal;
}








TGraphErrors *calibrateTime(int elchan, UsefulAtriStationEvent * owncalEvPtr, double *tempCorr, int stationId)
{
	double samplecorr[128] = {0};
	double sampleRms[128] = {0};
	double meanV1[40];
	double freqBevSev, freqBevSodd, freqBoddSev, freqBoddSodd, shift_odd, shift_ev, shift_evodd;
	TGraph *grOwnCal = owncalEvPtr->getGraphFromElecChan(elchan);
	int Nsamples = grOwnCal->GetN();

	ifstream calRead;
	char calReadchar[60];
	sprintf(calReadchar,(calFilesDir + "/calFileVs2ARA%2.2dCh%dSamples0.txt").c_str(), stationId, elchan);
	calRead.open(calReadchar);
	calRead >> shift_evodd >> shift_ev >> freqBevSev >> freqBoddSev;
	for(int d=0;d<64;d++)
	{
		calRead >> samplecorr[d*2] >> sampleRms[d*2];
	}
	calRead.close();
	sprintf(calReadchar,(calFilesDir + "/calFileVs2ARA%2.2dCh%dSamples1.txt").c_str(), stationId, elchan);
	calRead.open(calReadchar);
	calRead >> shift_evodd >> shift_odd >> freqBevSodd >> freqBoddSodd;
	for(int d=0;d<64;d++)
	{
		calRead >> samplecorr[d*2+1] >> sampleRms[d*2+1];
	}
	calRead.close();

	double grtime, grvolt;
	double corrTime = 0;

	double voltCalValues[128] ={0};
	ifstream in;
	char baseVCF[200];
	sprintf(baseVCF, (calFilesDir + "/baseVoltCorrARA%2.2dCH%d.txt").c_str(), stationId, elchan);
	in.open(baseVCF);
	if(in.good())
	{
		for(int i=0;i<128;i++)
		{
//			in >> voltCalValues[i];
			voltCalValues[i] = 0;
		}
		in.close();
	}
	else{
		in.close();
		for(int i=0;i<128;i++)
		{
			voltCalValues[i] = 0;
		}
	}
			
//	for(int s=0;s<Nsamples;s++)
//	{
//		if(s%64==0)meanV1[s/64]=0;
//		grOwnCal->GetPoint(s,grtime,grvolt);
//		meanV1[s/128]+=grvolt/128.0;
//	}

	double tempt, tempv;
	double sort_times[2048];
	double sort_volts[2048];
	double sort_voltCalValues[2048];
	double sortRms[2048];
	double pointRms = 0;
	double block_time = 0.0;
	double blockTarray[40][64] = {{0}};
	double blockVarray[40][64] = {{0}};
	int block_count = 0;
	int grBlockN  = 0;
	for(int j= 0;j<Nsamples;j++)
	{
		grOwnCal->GetPoint(j,grtime,grvolt);
		grBlockN = owncalEvPtr->blockVec[j/16].getBlock();
		blockTarray[block_count][(j)%64]=grtime/3.2;
		blockVarray[block_count][(j)%64]=grvolt;
		if((j)%64==63)block_count++;
	}

	double temp_max =  1.0E22;
	double temp_values[3];
	double max_values[3];

	int grBlock0 = 0;
	double sampleTime[128] = {0};
	for(int s=0;s<128;s++)
	{
		if(s<2){
			sampleTime[0] = 0.0;
			sampleTime[1] = 1.0/3.2 + shift_evodd;
		}
		else{
			sampleTime[s]= sampleTime[s-2] + (1.0/1.6*samplecorr[s])*tempCorr[s%2];
		}
	}

	for(int j= 0;j<Nsamples;j++)
	{
		grBlockN = owncalEvPtr->blockVec[j/16].getBlock();
		if(j==0)grBlock0 = owncalEvPtr->blockVec[0].getBlock();

		if(grBlock0%2==0)
		{
			if(j%2==0){
				 corrTime = (j/128)*40.0 + sampleTime[j%128];
				 pointRms = sampleRms[j%128];
				 sort_voltCalValues[j] = voltCalValues[j%128];
			}
			else{
				 corrTime = (j/128)*40.0 + sampleTime[j%128]; 
				 pointRms = sampleRms[j%128];
				 sort_voltCalValues[j] = voltCalValues[j%128];
			}
		}
		else
		{
			if(j%2==0){
				 corrTime =  ((j+64)/128)*40.0 - 20.0 + sampleTime[(j+64)%128];
//				 corrTime =  ((j+64)/128)*40.0 + sampleTime[(j+64)%128] - sampleTime[64];
				 pointRms = sampleRms[(j+64)%128];
				 sort_voltCalValues[j] = voltCalValues[(j+64)%128];
			}
			else{
				 corrTime =  ((j+64)/128)*40.0 - 20.0 + sampleTime[(j+64)%128];
//				 corrTime =  ((j+64)/128)*40.0 + sampleTime[(j+64)%128] - sampleTime[64];
				 pointRms = sampleRms[(j+64)%128];
				 sort_voltCalValues[j] = voltCalValues[(j+64)%128];
			}
		}



		sort_times[j] = corrTime - cableDelays[elchan/8][elchan%8];
		sort_volts[j] = blockVarray[j/64][j%64];
		sortRms[j] =pointRms;
	}
	
	TGraphErrors *NgrOwnCal = new TGraphErrors(Nsamples);
	for(int go=0;go<Nsamples;go++)
	{
		NgrOwnCal->SetPoint(go,sort_times[go],sort_volts[go] + sort_voltCalValues[go]);
		NgrOwnCal->SetPointError(go, sortRms[go], 0.0);
	}
	return NgrOwnCal;
}






TGraph * sortGraphS(TGraph * gr)
{
	const int Nsamples = gr->GetN();
	double sort_times[Nsamples], sort_volts[Nsamples];
	double sortXErr[Nsamples];
	double sortYErr[Nsamples];
	double times, volts;
	double tErr, vErr;
	for(int i=0;i<gr->GetN();i++)
	{
		gr->GetPoint(i, times, volts);
		tErr = gr->GetErrorX(i);
		vErr = gr->GetErrorY(i);
		sort_times[i] = times;
		sort_volts[i] = volts;
		sortXErr[i] = tErr;
		sortYErr[i] = vErr;
	}

	double sorted_times[2048];
	double sorted_volts[2048];
	double sortedXErr[2048];
	double sortedYErr[2048];
	double tempT, tempV, tempEX, tempEY;
	int it = 0;
	int err_it = 0;
	for(int i=0;i<Nsamples;i++)
	{
		sorted_times[i] = -1;
		sorted_volts[i] = 0;
		sortedXErr[i] = 0;
		sortedYErr[i] = 0;
	}
	for(int so=0;so<Nsamples;so++)
	{
		it = 0;
		for(int to=0;to<Nsamples;to++)
		{
			if(sort_times[so]>sort_times[to]) it++;
		}
		if(it<Nsamples){
			sorted_times[it] = sort_times[so];
			sorted_volts[it] = sort_volts[so];
			sortedXErr[it] = sortXErr[so];
			sortedYErr[it] = sortYErr[so];
		}
		
	}
	int mover = 0;
	
	TGraph *NgrOwnCal = new TGraph(Nsamples);
	for(int go=0;go<Nsamples;go++)
	{
		if(go+mover<Nsamples){	
			if(sorted_times[go+mover]==-1) mover++;
			NgrOwnCal->SetPoint(go,sorted_times[go+mover],sorted_volts[go+mover]);
//			NgrOwnCal->SetPointError(go,sortedXErr[go+mover],sortedYErr[go+mover]);
		}
		else{				
			NgrOwnCal->SetPoint(go,sorted_times[Nsamples]+1.0/1000.0*go,0);
//			NgrOwnCal->SetPointError(go,0,0);
		}	
	}
	return NgrOwnCal;
}













TGraphErrors * sortGraph(TGraphErrors * gr)
{
	const int Nsamples = gr->GetN();
	double sort_times[Nsamples], sort_volts[Nsamples];
	double sortXErr[Nsamples];
	double sortYErr[Nsamples];
	double times, volts;
	double tErr, vErr;
	for(int i=0;i<gr->GetN();i++)
	{
		gr->GetPoint(i, times, volts);
		tErr = gr->GetErrorX(i);
		vErr = gr->GetErrorY(i);
		sort_times[i] = times;
		sort_volts[i] = volts;
		sortXErr[i] = tErr;
		sortYErr[i] = vErr;
	}

	double sorted_times[2048];
	double sorted_volts[2048];
	double sortedXErr[2048];
	double sortedYErr[2048];
	double tempT, tempV, tempEX, tempEY;
	int it = 0;
	int err_it = 0;
	for(int i=0;i<Nsamples;i++)
	{
		sorted_times[i] = -1;
		sorted_volts[i] = 0;
		sortedXErr[i] = 0;
		sortedYErr[i] = 0;
	}
	for(int so=0;so<Nsamples;so++)
	{
		it = 0;
		for(int to=0;to<Nsamples;to++)
		{
			if(sort_times[so]>sort_times[to]) it++;
		}
		if(it<Nsamples){
			sorted_times[it] = sort_times[so];
			sorted_volts[it] = sort_volts[so];
			sortedXErr[it] = sortXErr[so];
			sortedYErr[it] = sortYErr[so];
		}
		
	}
	int mover = 0;
	
	TGraphErrors *NgrOwnCal = new TGraphErrors(Nsamples);
	for(int go=0;go<Nsamples;go++)
	{
		if(go+mover<Nsamples){	
			if(sorted_times[go+mover]==-1) mover++;
			NgrOwnCal->SetPoint(go,sorted_times[go+mover],sorted_volts[go+mover]);
			NgrOwnCal->SetPointError(go,sortedXErr[go+mover],sortedYErr[go+mover]);
		}
		else{				
			NgrOwnCal->SetPoint(go,sorted_times[Nsamples]+1.0/1000.0*go,0);
			NgrOwnCal->SetPointError(go,0,0);
		}	
	}
	return NgrOwnCal;
}


TGraphErrors * calibrateVolts(UsefulAtriStationEvent * owncalEvPtr, int chip, int elchan, TGraphErrors *calGrin, int stationId)
{
	TGraphErrors * calGr = new TGraphErrors(calGrin->GetN());
	double newVolts = 0;
	double voltError = 0;
	double overEstimate = 0;
	double times, volts, XErr, YErr;
	for(int s=0;s<calGrin->GetN();s++)
	{
		int grBlockN = owncalEvPtr->blockVec[s/16].getBlock();
		calGrin->GetPoint(s,times,volts);
		XErr = calGrin->GetErrorX(s);
//		YErr = calGrin->GetErrorY(s);
		while(calConst[chip][elchan][grBlockN][s%64][8]>1.0) grBlockN = (grBlockN - 2 + 512)%512;		
		if(TMath::Abs(volts)<400)
		{
			if(volts>0)
			{
				 newVolts = calConst[chip][elchan][grBlockN][s%64][6] 
						+ (volts - calConst[chip][elchan][grBlockN][s%64][7])*calConst[chip][elchan][grBlockN][s%64][0] 
						+ pow((volts - calConst[chip][elchan][grBlockN][s%64][7]), 2)*calConst[chip][elchan][grBlockN][s%64][1] 
						+ pow((volts - calConst[chip][elchan][grBlockN][s%64][7]), 3)*calConst[chip][elchan][grBlockN][s%64][2]; 
			}
			else
			{
				 newVolts = calConst[chip][elchan][grBlockN][s%64][6] 
						+ (volts - calConst[chip][elchan][grBlockN][s%64][7])*calConst[chip][elchan][grBlockN][s%64][3] 
						+ pow((volts - calConst[chip][elchan][grBlockN][s%64][7]), 2)*calConst[chip][elchan][grBlockN][s%64][4]
						+  pow((volts - calConst[chip][elchan][grBlockN][s%64][7]), 3)*calConst[chip][elchan][grBlockN][s%64][5];
			}
		}
		else
		{
			if(volts>0)
			{
				 newVolts = calConstHigh[chip][elchan][grBlockN][s%64][0] + volts*calConstHigh[chip][elchan][grBlockN][s%64][1];
			}
			else
			{
				 newVolts = calConstHigh[chip][elchan][grBlockN][s%64][2] + volts*calConstHigh[chip][elchan][grBlockN][s%64][3];
			}
		}

//		overEstimate = TMath::Sqrt(calConstErr[chip][elchan][grBlockN][s%64][8] );
			
/*		if(volts>0)
		{
			 voltError = TMath::Sqrt(  pow(calConstErr[chip][elchan][grBlockN][s%64][6]*overEstimate,2)
					+ pow((calConstErr[chip][elchan][grBlockN][s%64][7])*overEstimate*(
						-calConst[chip][elchan][grBlockN][s%64][0] 
					- 2.0*pow((volts - calConst[chip][elchan][grBlockN][s%64][7]), 1)*calConst[chip][elchan][grBlockN][s%64][1] 
					- 3.0*pow((volts - calConst[chip][elchan][grBlockN][s%64][7]), 2)*calConst[chip][elchan][grBlockN][s%64][2]),2)
					+ pow( (volts - calConst[chip][elchan][grBlockN][s%64][7])*calConstErr[chip][elchan][grBlockN][s%64][0]*overEstimate,2 )
					+ pow( pow((volts - calConst[chip][elchan][grBlockN][s%64][7]), 2)*calConstErr[chip][elchan][grBlockN][s%64][1]*overEstimate,2  )
					+ pow( pow((volts - calConst[chip][elchan][grBlockN][s%64][7]), 3)*calConstErr[chip][elchan][grBlockN][s%64][2]*overEstimate  ,2) 
					   ); 
		}
		else
		{
			 voltError = TMath::Sqrt(  pow(calConstErr[chip][elchan][grBlockN][s%64][6]*overEstimate,2)
					+ pow((calConstErr[chip][elchan][grBlockN][s%64][7])*overEstimate*(
						-calConst[chip][elchan][grBlockN][s%64][3] 
					- 2.0*pow((volts - calConst[chip][elchan][grBlockN][s%64][7]), 1)*calConst[chip][elchan][grBlockN][s%64][4] 
					- 3.0*pow((volts - calConst[chip][elchan][grBlockN][s%64][7]), 2)*calConst[chip][elchan][grBlockN][s%64][5]),2)
					+ pow( (volts - calConst[chip][elchan][grBlockN][s%64][7])*calConstErr[chip][elchan][grBlockN][s%64][3]*overEstimate,2 )
					+ pow( pow((volts - calConst[chip][elchan][grBlockN][s%64][7]), 2)*calConstErr[chip][elchan][grBlockN][s%64][4]*overEstimate,2  )
					+ pow( pow((volts - calConst[chip][elchan][grBlockN][s%64][7]), 3)*calConstErr[chip][elchan][grBlockN][s%64][5]*overEstimate  ,2) 
					   ); 
		}
*/
		calGr->SetPoint(s, times, newVolts);
//		calGr->SetPointError(s, XErr, calConst[chip][elchan][grBlockN][s%64][8]);
		calGr->SetPointError(s, XErr, 0/*voltError*/);
	}
	return calGr;
}










TGraphErrors * getVoltageError(UsefulAtriStationEvent * owncalEvPtr, int chip, int elchan, TGraphErrors *calGrin, int stationId)
{
	TGraphErrors * calGr = new TGraphErrors(calGrin->GetN());
	double voltError = 0;
	double times, volts, XErr, YErr;
	int firstBlock = 0;
	firstBlock = owncalEvPtr->blockVec[0].getBlock();
//	int division = calGrin->GetN()/owncalEvPtr->blockVec.size();
	for(int s=0;s<calGrin->GetN();s++)
	{
		int grBlockN = owncalEvPtr->blockVec[s/16].getBlock();
		calGrin->GetPoint(s,times,volts);
		XErr = calGrin->GetErrorX(s);
		YErr = calGrin->GetErrorY(s);
		while(calConst[chip][elchan][grBlockN][s%64][8]>1.0) grBlockN = (grBlockN - 2 + 512)%512;		
//		if(s%64==1 || s%64==3) grBlockN = 0;		
//		if(TMath::Abs(volts)<600)
		{
			if(volts>0)
			{
				 voltError = TMath::Sqrt(  pow(calConstErr[chip][elchan][grBlockN][s%64][6],2)
						+ pow((calConstErr[chip][elchan][grBlockN][s%64][7])*(
							-calConst[chip][elchan][grBlockN][s%64][0] 
						- 2.0*pow((volts - calConst[chip][elchan][grBlockN][s%64][7]), 1)*calConst[chip][elchan][grBlockN][s%64][1] 
						- 3.0*pow((volts - calConst[chip][elchan][grBlockN][s%64][7]), 2)*calConst[chip][elchan][grBlockN][s%64][2]),2)
						+ pow( (volts - calConst[chip][elchan][grBlockN][s%64][7])*calConstErr[chip][elchan][grBlockN][s%64][0],2 )
						+ pow( pow((volts - calConst[chip][elchan][grBlockN][s%64][7]), 2)*calConstErr[chip][elchan][grBlockN][s%64][1],2  )
						+ pow( pow((volts - calConst[chip][elchan][grBlockN][s%64][7]), 3)*calConstErr[chip][elchan][grBlockN][s%64][2]  ,2)    ); 
			}
			else
			{
				 voltError = TMath::Sqrt(  pow(calConstErr[chip][elchan][grBlockN][s%64][6],2)
						+ pow((calConstErr[chip][elchan][grBlockN][s%64][7])*(
							-calConst[chip][elchan][grBlockN][s%64][3] 
						- 2.0*pow((volts - calConst[chip][elchan][grBlockN][s%64][7]), 1)*calConst[chip][elchan][grBlockN][s%64][4] 
						- 3.0*pow((volts - calConst[chip][elchan][grBlockN][s%64][7]), 2)*calConst[chip][elchan][grBlockN][s%64][5]),2)
						+ pow( (volts - calConst[chip][elchan][grBlockN][s%64][7])*calConstErr[chip][elchan][grBlockN][s%64][3],2 )
						+ pow( pow((volts - calConst[chip][elchan][grBlockN][s%64][7]), 2)*calConstErr[chip][elchan][grBlockN][s%64][4],2  )
						+ pow( pow((volts - calConst[chip][elchan][grBlockN][s%64][7]), 3)*calConstErr[chip][elchan][grBlockN][s%64][5]  ,2)    ); 
			}
		}
//		else
//		{
//			if(volts>0)
//			{
//				 newVolts = calConstHigh[chip][elchan][grBlockN][s%64][0] + volts*calConstHigh[chip][elchan][grBlockN][s%64][1];
//			}
//			else
//			{
//				 newVolts = calConstHigh[chip][elchan][grBlockN][s%64][2] + volts*calConstHigh[chip][elchan][grBlockN][s%64][3];
//			}
//		}
//		calGr->SetPoint(s, times, newVolts);
		calGr->SetPointError(s, XErr, voltError);
	}
	return calGr;
}






int cutOverTimes(TGraphErrors * inGr, TGraphErrors * outGr, int halfFul=0)
{
	if(halfFul==0){
	int SdtCount = 0;
	double XErr, YErr, dt, times1, volts1;
	double Ltimes = -5.0;
	for(int i=0;i<inGr->GetN();i++)
	{
		inGr->GetPoint(i,times1,volts1);
		XErr = inGr->GetErrorX(i);
		YErr = inGr->GetErrorY(i);
		if(i%64>50 && times1>( (i/64)+1)*20.0)
		{
			SdtCount++;
		}
		else
		{
			outGr->SetPoint(i - SdtCount,times1,volts1);
			outGr->SetPointError(i - SdtCount, XErr, YErr);
			
		}
	}	
	return SdtCount;
	}
	else{
	int SdtCount = 0;
	double XErr, YErr, dt, times1, volts1;
	double Ltimes = -5.0;
	for(int i=0;i<inGr->GetN();i++)
	{
		inGr->GetPoint(i,times1,volts1);
		XErr = inGr->GetErrorX(i);
		YErr = inGr->GetErrorY(i);
		if(i%32>25 && times1>( (i/32)+1)*20.0)
		{
			SdtCount++;
		}
		else
		{
			outGr->SetPoint(i - SdtCount,times1,volts1);
			outGr->SetPointError(i - SdtCount, XErr, YErr);
			
		}
	}	
	return SdtCount;
	}
}





int selectSamples(TGraphErrors * SgrOwnCal, TGraphErrors * SgrOwnCal2)
{
	int SdtCount = 0;
	double XErr, XErr2, dt, times1, volts1, YErr;
	double Ltimes = -5.0;
	for(int i=0;i<SgrOwnCal->GetN();i++)
	{
		SgrOwnCal->GetPoint(i,times1,volts1);
		XErr = SgrOwnCal->GetErrorX(i);
		YErr = SgrOwnCal->GetErrorY(i);
		double XErr2 = SgrOwnCal->GetErrorX(i-1);
		dt = times1 - Ltimes;
		if(XErr>0.2)
		{
			SdtCount++;
		}
		else
		{
			if(dt<0.12)
			{
				if(XErr2<=0.2){
					if(XErr>XErr2 )	SdtCount++;
					else {SdtCount++; SgrOwnCal2->SetPoint(i - SdtCount,times1,volts1);SgrOwnCal2->SetPointError(i - SdtCount, XErr, YErr);}
				}
				else {SgrOwnCal2->SetPoint(i - SdtCount,times1,volts1);SgrOwnCal2->SetPointError(i - SdtCount, XErr, YErr);}
			}
			else
			{
				SgrOwnCal2->SetPoint(i - SdtCount,times1,volts1);
				SgrOwnCal2->SetPointError(i - SdtCount, XErr, YErr);
			}
		}		
		Ltimes = times1;	
	}
	return SdtCount;
}











TGraph * calibrateVoltsS(UsefulAtriStationEvent * owncalEvPtr, int chip, int elchan, TGraph *calGrin, int stationId)
{
	TGraph * calGr = new TGraph(calGrin->GetN());
	double newVolts = 0;
	double times, volts, XErr, YErr;
	int firstBlock = 0;
//	int division = calGrin->GetN()/owncalEvPtr->blockVec.size();
	for(int s=0;s<calGrin->GetN();s++)
	{
		int grBlockN = owncalEvPtr->blockVec[s/16].getBlock();
		calGrin->GetPoint(s,times,volts);
//		XErr = calGrin->GetErrorX(s);
//		YErr = calGrin->GetErrorY(s);
		while(calConst[chip][elchan][grBlockN][s%64][8]>1.0) grBlockN = (grBlockN - 2 + 512)%512;		
		if(TMath::Abs(volts)<400)
		{
			if(volts>0)
			{
				 newVolts = calConst[chip][elchan][grBlockN][s%64][6] 
						+ (volts - calConst[chip][elchan][grBlockN][s%64][7])*calConst[chip][elchan][grBlockN][s%64][0] 
						+ pow((volts - calConst[chip][elchan][grBlockN][s%64][7]), 2)*calConst[chip][elchan][grBlockN][s%64][1] 
						+ pow((volts - calConst[chip][elchan][grBlockN][s%64][7]), 3)*calConst[chip][elchan][grBlockN][s%64][2]; 
			}
			else
			{
				 newVolts = calConst[chip][elchan][grBlockN][s%64][6] 
						+ (volts - calConst[chip][elchan][grBlockN][s%64][7])*calConst[chip][elchan][grBlockN][s%64][3] 
						+ pow((volts - calConst[chip][elchan][grBlockN][s%64][7]), 2)*calConst[chip][elchan][grBlockN][s%64][4]
						+  pow((volts - calConst[chip][elchan][grBlockN][s%64][7]), 3)*calConst[chip][elchan][grBlockN][s%64][5];
			}
		}
		else
		{
			if(volts>0)
			{
				 newVolts = calConstHigh[chip][elchan][grBlockN][s%64][0] + volts*calConstHigh[chip][elchan][grBlockN][s%64][1];
				if(calConstHigh[chip][elchan][grBlockN][s%64][1]<0.0)cerr <<"CALIBRATION ALERT!!" << endl;
			}
			else
			{
				 newVolts = calConstHigh[chip][elchan][grBlockN][s%64][2] + volts*calConstHigh[chip][elchan][grBlockN][s%64][3];
				if(calConstHigh[chip][elchan][grBlockN][s%64][3]<0.0)cerr <<"CALIBRATION ALERT!!" << endl;
			}
		}
//		else
//		{
//			 if(firstBlock%2==0)newVolts = baseCalConst[s%128][1] + volts*baseCalConst[s%128][2]; 
//			 else newVolts = baseCalConst[(s+64)%128][1] + volts*baseCalConst[(s+64)%128][2]; 
//		}
//		if(TMath::Abs(newVolts) >550 && TMath::Abs(volts) - TMath::Abs(newVolts)<80 )
//		{
//			 if(firstBlock%2==0)newVolts = baseCalConst[s%128][1] + volts*baseCalConst[s%128][2]; 
//			 else newVolts = baseCalConst[(s+64)%128][1] + volts*baseCalConst[(s+64)%128][2]; 
//		}
	
		calGr->SetPoint(s, times, newVolts);
//		calGr->SetPointError(s, XErr, calConst[chip][elchan][grBlockN][s%64][8]);
	}
	return calGr;
}



















TGraph * calibrateBasicVoltsS(UsefulAtriStationEvent * owncalEvPtr, int chip, int elchan, TGraph *calGrin, int stationId)
{
	TGraph * calGr = new TGraph(calGrin->GetN());
	double newVolts = 0;
	double times, volts, XErr, YErr;
	int firstBlock = 0;
	firstBlock = owncalEvPtr->blockVec[0].getBlock();
//	int division = calGrin->GetN()/owncalEvPtr->blockVec.size();
	for(int s=0;s<calGrin->GetN();s++)
	{
		int grBlockN = owncalEvPtr->blockVec[s/16].getBlock();
		calGrin->GetPoint(s,times,volts);
//		XErr = calGrin->GetErrorX(s);
//		YErr = calGrin->GetErrorY(s);
		newVolts = 0;
		if(firstBlock%2==0){
			if(TMath::Abs(volts)<510)
			{
				 newVolts = baseCalConst[s%128][3] + volts*baseCalConst[s%128][4] 
								+ volts*volts*baseCalConst[s%128][5]
								+ volts*volts*volts*baseCalConst[s%128][6]
								+ volts*volts*volts*volts*baseCalConst[s%128][7]
								+ volts*volts*volts*volts*volts*baseCalConst[s%128][8]
								+ volts*volts*volts*volts*volts*volts*baseCalConst[s%128][9]
								+ volts*volts*volts*volts*volts*volts*volts*baseCalConst[s%128][10]
								+ volts*volts*volts*volts*volts*volts*volts*volts*baseCalConst[s%128][11]
								+ volts*volts*volts*volts*volts*volts*volts*volts*volts*baseCalConst[s%128][12];
			}
			else
			{
				 newVolts = baseCalConst[s%128][1] + volts*baseCalConst[s%128][2]; 
			}
		}
		else{
			if(TMath::Abs(volts)<510)
			{
				 newVolts = baseCalConst[(s+64)%128][3] + volts*baseCalConst[(s+64)%128][4] 
								+ volts*volts*baseCalConst[(s+64)%128][5]
								+ volts*volts*volts*baseCalConst[(s+64)%128][6]
								+ volts*volts*volts*volts*baseCalConst[(s+64)%128][7]
								+ volts*volts*volts*volts*volts*baseCalConst[(s+64)%128][8]
								+ volts*volts*volts*volts*volts*volts*baseCalConst[(s+64)%128][9]
								+ volts*volts*volts*volts*volts*volts*volts*baseCalConst[(s+64)%128][10]
								+ volts*volts*volts*volts*volts*volts*volts*volts*baseCalConst[(s+64)%128][11]
								+ volts*volts*volts*volts*volts*volts*volts*volts*volts*baseCalConst[(s+64)%128][12];
			}
			else
			{
				 newVolts = baseCalConst[(s+64)%128][1] + volts*baseCalConst[(s+64)%128][2]; 
			}

		}
		calGr->SetPoint(s, times, newVolts);
//		calGr->SetPointError(s, XErr, calConst[chip][elchan][grBlockN][s%64][8]);
	}
	return calGr;
}

void loadBasicVoltCalValues( int stationID, int elchan )
{
	int maxSample = 128;
	ifstream in;
	double temp;
	char infile[200];
	sprintf(infile,"/home/meures/AraCalibration/ARA%2.2d/finalVCalOutFileCH%d.txt",stationID, elchan);
	in.open(infile);
	

		    for(int l=0;l<maxSample;l++)
		    {
			for(int m=0;m<13;m++)
			{
			     if(in.good()){
				in >> temp;
				baseCalConst[l][m] = temp;
//				if(m==12) cerr << baseCalConst[l][m] << endl;
			     }
			}
		    }
	in.close();
}










void loadHighVoltCalValues( int stationID )
{
	int maxChip = 4;  //FIXME: So far only values for one chip
	int maxChan = 6;  //FIXME: Just one calibrated channel so far
	int maxBlock = 512;
	int maxSample = 64;
	ifstream in;
	int tChip, tChannel, tBlock;
	double temp;
	char infile[200];
	sprintf(infile, (calFilesDir +  "/voltCalibHighARA%2.2d.txt").c_str(), stationID);
	in.open(infile);
	if(in.good()){	

	for(int i = 0;i<maxChip;i++)
	{
	   if(i==0) maxChan = 6;
	   else if(i==3) maxChan = 6;
	   else maxChan=4;
  	    for(int j = 0;j<maxChan;j++)
	    {
		for(int k=0;k<maxBlock;k++)
		{
		    in >> tChip >> tChannel >> tBlock;
		    for(int l=0;l<maxSample;l++)
		    {
			for(int m=0;m<5;m++)
			{
			     if(in.good()){
				in >> temp;
				calConstHigh[tChip][tChannel][k][l][m] = temp;
			     }
			}
		    }
		}
	    }
	}
	}
	in.close();
}













void loadVoltCalValues( int stationID )
{
	int maxChip = 4;  //FIXME: So far only values for one chip
	int maxChan = 6;  //FIXME: Just one calibrated channel so far
	int maxBlock = 512;
	int maxSample = 64;
	ifstream in;
	int tChip, tChannel, tBlock;
	double temp;
	char infile[200];
	sprintf(infile, (calFilesDir + "/voltCalibVppCorrectDerivativeVs3ARA%2.2d.txt").c_str(), stationID);
	in.open(infile);
	
	if(in.good()){
	for(int i = 0;i<maxChip;i++)
	{
	   if(i==0) maxChan = 6;
	   else if(i==3) maxChan = 6;
	   else maxChan=4;
  	    for(int j = 0;j<maxChan;j++)
	    {
		for(int k=0;k<maxBlock;k++)
		{
		    in >> tChip >> tChannel >> tBlock;
		    for(int l=0;l<maxSample;l++)
		    {
			for(int m=0;m<9;m++)
			{
			     if(in.good()){
				in >> temp;
				calConst[tChip][tChannel][k][l][m] = temp;
			     }
			}
		    }
		}
	    }
	}
	}
	in.close();
}



void loadCableDelays(int stationId)
{
	int maxChip = 4;  //FIXME: So far only values for one chip
	int maxChan = 6;  //FIXME: Just one calibrated channel so far
	ifstream in;
	int tChip, tChannel;
	double temp;
	char infile[200];
	sprintf(infile, (calFilesDir + "/cableDelaysARA%2.2d.txt").c_str(), stationId);
	in.open(infile);
	

	for(int i = 0;i<maxChip;i++)
	{
	   if(i==0) maxChan = 6;
	   else if(i==3) maxChan = 6;
	   else maxChan=4;
  	    for(int j = 0;j<maxChan;j++)
	    {
		in >> tChip >> tChannel;
	     	if(in.good())
		{
			in >> temp;
			cableDelays[tChip][tChannel] = temp;
		}
	    }
	}
	in.close();
}








void loadVoltCalErrors( int stationID )
{
	int maxChip = 4;  //FIXME: So far only values for one chip
	int maxChan = 6;  //FIXME: Just one calibrated channel so far
	int maxBlock = 512;
	int maxSample = 64;
	ifstream in;
	int tChip, tChannel, tBlock;
	double temp;
	char infile[200];
	sprintf(infile, (calFilesDir + "/voltCalibErrorVs3ARA%2.2d.txt").c_str(), stationID);
	in.open(infile);
	

	for(int i = 0;i<maxChip;i++)
	{
	   if(i==0) maxChan = 6;
	   else if(i==3) maxChan = 6;
	   else maxChan=4;
  	    for(int j = 0;j<maxChan;j++)
	    {
		for(int k=0;k<maxBlock;k++)
		{
		    in >> tChip >> tChannel >> tBlock;
		    for(int l=0;l<maxSample;l++)
		    {
			for(int m=0;m<9;m++)
			{
			     if(in.good()){
				in >> temp;
				calConstErr[tChip][tChannel][k][l][m] = temp;
			     }
			}
		    }
		}
	    }
	}
	in.close();
}



















TGraph * halfSamplesS(TGraph *ingr, int evodd)
{
//evodd can be 0 = even, 1=odd.
	TGraph * grNew = new TGraph(ingr->GetN()/2);
	double times, volts;
	for(int i = 0;i<ingr->GetN();i++)
	{
		ingr->GetPoint(i, times, volts);
		if(i%2==evodd){grNew->SetPoint(i/2,times,volts);}
	}
	return grNew;
}





TGraphErrors * halfSamples(TGraphErrors *ingr, int evodd)
{
//evodd can be 0 = even, 1=odd.
	TGraphErrors * grNew = new TGraphErrors(ingr->GetN()/2);
	double times, volts, XErr, YErr;
	for(int i = 0;i<ingr->GetN();i++)
	{
		ingr->GetPoint(i, times, volts);
		XErr = ingr->GetErrorX(i);
		YErr = ingr->GetErrorY(i);
		if(i%2==evodd){grNew->SetPoint(i/2,times,volts);grNew->SetPointError(i/2,XErr,YErr);}
	}
	return grNew;
}



void writeGraphToFile(TGraph * ingr, string fileName)
{
	double temp1, temp2;
	double times;
	double volts;
	fstream out;
	out.open(fileName.c_str(), ios::out);
	for(int i=0;i<ingr->GetN();i++)
	{
		ingr->GetPoint(i,temp1, temp2);
		times = temp1;
		volts = temp2;
//		cout << times << "\t" << volts << endl;

		out << times << "\t" << volts << endl;
	}
	out.close();
}



void writeHistoToFile(TH1D * inh, string fileName)
{
	double times, volts;
	fstream out;
	out.open(fileName.c_str(), ios::out);
	for(int i=0;i<inh->GetNbinsX();i++)
	{
		volts = inh->GetBinContent(i);
		times = inh->GetBinCenter(i);
		out << times << "\t" << volts << endl;
	}
	out.close();
}



int getTempCorr(TTree *eventTree, int runNumber, int stationID, int chipNumber, double *meanSqCfreq, AraEventCalibrator *calib, string rootBaseDir )
{
	RawAtriStationEvent *rawAtriEvPtr=0;
	RawAraStationEvent *rawEvPtr=0;
	UsefulAtriStationEvent *realAtriEvPtr=0;
	eventTree->ResetBranchAddresses();
	eventTree->SetBranchAddress("event", &rawAtriEvPtr);

   Long64_t numEntries=eventTree->GetEntries();
   cout << "Total number of events in the read files: " << numEntries << endl;
   
   int elchan = chipNumber*8 + 7;
   
	char run_char[30];
	sprintf(run_char, "run%d", runNumber);
	string run_string = string(run_char);

	char station_char[30];
	sprintf(station_char, "ARA%02d", stationID);
	string station_string = string(station_char);

//	char dir_char[200];
//	int ped_run_no = 367;
//	sprintf(dir_char,"/home/meures/AraCalibration/uhen/ara/data/calibration/ARA03/raw_data/run_%06d/pedestalValues.run%06d.dat", ped_run_no, ped_run_no );
//	AraEventCalibrator *calib = AraEventCalibrator::Instance();
//	calib->setAtriPedFile(dir_char, stationID);


  	int pass_evts = 0;
	char elchan_char[30];
	sprintf(elchan_char, "elChan%d", elchan);
	string elchan_string = string(elchan_char);

	Long64_t eventNumber = 0;
	double first_value = 0.0;	

	cerr << "Start loop!" << endl;
	
	int count_blocks = 0;
	double times = 0;
	double volts = 0;
	double block_time = 0;

	TGraph *grSQNoCalc = 0; 
	TCanvas *csq = new TCanvas("csq");
	TCanvas *csqO = new TCanvas("csqO");

	int block_number; 
	double meanSqc = 0.0;
	TH1D *hzc = new TH1D("hzc", "mean crossings", 2000, 0, 3);
	TH1D *hffc = new TH1D("hffc", "first crossings", 10000, 0, 200);
	TH1D *hffcO = new TH1D("hffcO", "first crossings", 10000, 0, 200);
	TH1D *hmc = new TH1D("hmc", "Time vs tempsqc", 10000, -100, 100);
	TH1D *hmcO = new TH1D("hmcO", "Time vs tempsqc", 10000, -100, 100);
	TH1D *hzcO = new TH1D("hzcO", "mean crossings odd samples", 2000, 0, 3);
	double tempSqVc=0;
	double tempSqTc=0;
	double tempSqVc2=0;
	double tempSqVc3=0;
	double dummyT = 0;
	double firstCrossingc = 0.0;
	int firstc = 0;
	int firstfirstc= 0;
	int countCrossingc = 0;
	int countCrossingcO = 0;
	double crossingPosc = 0.0;
	double meanCrossingsBevSevc[2] = {0};
	double meanCrossingsBoddSevc[2] = {0};
	double meanCrossingsBevSoddc[2] = {0};
	double meanCrossingsBoddSoddc[2] = {0};
	double sqtimesc, sqvoltsc;

  	TFile *fout = new TFile((rootBaseDir + "temp_cal_data_"+elchan_string + "_" +run_string+".root").c_str(),"RECREATE");
	double lastCrossingc =0.0;
	int secondCross = 0;	
	int thirdCross = 0;	
	int Cross4 = 0;	
	int Cross5 = 0;	
	int Cross6 = 0;	
	int Cross7 = 0;	
	int Cross8 = 0;
	cerr << "We are running on channel: " << elchan << endl;
	
	for(Long64_t event=0;event<numEntries;event++) 
	{//loop all events
     
     	//This line gets the RawIcrr or RawAtri Event
     	eventTree->GetEntry(event);
     	eventNumber=event;
     
     	realAtriEvPtr = new UsefulAtriStationEvent(rawAtriEvPtr, AraCalType::kJustPed);

		
		TGraph *gr1 = realAtriEvPtr->getGraphFromElecChan(elchan);
		meanSqc = 0.0;
		grSQNoCalc = realAtriEvPtr->getGraphFromElecChan(elchan);
		int blockN = realAtriEvPtr->blockVec[0].getBlock();
		TGraph *grSqc = new TGraph(grSQNoCalc->GetN()/2);
		TGraph *grSqcO = new TGraph(grSQNoCalc->GetN()/2);
		if(blockN%2==0){
		for(int s=0;s<grSQNoCalc->GetN();s++)
		{
			grSQNoCalc->GetPoint(s,sqtimesc,sqvoltsc);
			if(s%2==0)grSqc->SetPoint(s/2,s/3.2,sqvoltsc);
			if(s%2==1)grSqcO->SetPoint(s/2,s/3.2,sqvoltsc);
			meanSqc +=sqvoltsc/grSQNoCalc->GetN();				
		}
		csq->cd();
		grSqc->Draw("ALP");
		csqO->cd();
		grSqcO->Draw("ALP");
		if(event<5)csq->Write();
		if(event<5)csqO->Write();
//		cerr << "The mean is: " << meanSqc << "   ";
		
		firstc = 0;
		firstfirstc=0;
		lastCrossingc =0.0;
		secondCross=0;
		thirdCross = 0;	
		Cross4 = 0;	
		Cross5 = 0;	
		Cross6 = 0;	
		Cross7 = 0;	
		Cross8 = 0;	

		for(int z=64;z<grSqc->GetN();z++)
		{
			grSqc->GetPoint(z,sqtimesc,sqvoltsc);
			if(z<grSqc->GetN()-1)grSqc->GetPoint(z+1,dummyT,tempSqVc3);
			if(z==64){tempSqVc2=sqvoltsc; tempSqVc=sqvoltsc;tempSqTc=sqtimesc;}	
			else
			{
				if(/* sqtimesc>(0.1 + 40.0*(z/64)) &&*/z>65 && ( (meanSqc<tempSqVc3 && meanSqc>tempSqVc2 && meanSqc>tempSqVc && meanSqc<sqvoltsc) || (meanSqc>tempSqVc3 && meanSqc<tempSqVc2 && meanSqc<tempSqVc && meanSqc>sqvoltsc) ) )
				{
					if(firstfirstc==0)
					{
//						if(meanSqc>tempSqVc && meanSqc<sqvoltsc){
						hffc->Fill( tempSqTc + 1.0/(-tempSqVc + sqvoltsc)*(-tempSqTc + sqtimesc)*(meanSqc - tempSqVc));
						firstfirstc=1;
//						}
					}
	
					if(firstc==0)
					{
						if(stationID==2 || (meanSqc>tempSqVc && meanSqc<sqvoltsc)){
						firstCrossingc = tempSqTc + 1.0/(-tempSqVc + sqvoltsc)*(-tempSqTc + sqtimesc)*(meanSqc - tempSqVc);
						hffc->Fill( tempSqTc + 1.0/(-tempSqVc + sqvoltsc)*(-tempSqTc + sqtimesc)*(meanSqc - tempSqVc));
//						hzc->Fill(0.0);
						meanCrossingsBevSevc[0] = firstCrossingc;
						lastCrossingc = 0.0;
						firstc=1;
						}
//						cerr << " first: " << firstCrossingc;
					}		
					else if(secondCross==0)
					{
						crossingPosc = tempSqTc + 1.0/(-tempSqVc + sqvoltsc)*(-tempSqTc + sqtimesc)*(meanSqc - tempSqVc) - firstCrossingc; 
						if(crossingPosc>2.0 && crossingPosc<8.0){ 
			//				meanCrossingsBevSevc[1] += 5.0/(crossingPosc - lastCrossingc);
			//				countCrossingc++;
			//				hzc->Fill(5.0/(crossingPosc - lastCrossingc));
						}
						lastCrossingc = crossingPosc;
						secondCross=1;
					}
					else if(thirdCross==0)
					{
						crossingPosc = tempSqTc + 1.0/(-tempSqVc + sqvoltsc)*(-tempSqTc + sqtimesc)*(meanSqc - tempSqVc) - firstCrossingc; 
						if(crossingPosc>7.0 && crossingPosc<13.0){ 
			//				meanCrossingsBevSevc[1] += 10.0/(crossingPosc);
			//				countCrossingc++;
			//				hzc->Fill(10.0/(crossingPosc));
						}
						lastCrossingc = crossingPosc;
						thirdCross=1;
					}
					else if(Cross4==0)
					{
						crossingPosc = tempSqTc + 1.0/(-tempSqVc + sqvoltsc)*(-tempSqTc + sqtimesc)*(meanSqc - tempSqVc) - firstCrossingc; 
						if(crossingPosc>12.0 && crossingPosc<18.0){ 
			//				meanCrossingsBevSevc[1] += 15.0/(crossingPosc);
			//				countCrossingc++;
			//				hzc->Fill(15.0/(crossingPosc));
						}
						lastCrossingc = crossingPosc;
						Cross4=1;
					}
					else if(Cross5==0)
					{
						crossingPosc = tempSqTc + 1.0/(-tempSqVc + sqvoltsc)*(-tempSqTc + sqtimesc)*(meanSqc - tempSqVc) - firstCrossingc; 
						if(crossingPosc>17.0 && crossingPosc<23.0){ 
			//				meanCrossingsBevSevc[1] += 20.0/(crossingPosc);
			//				countCrossingc++;
			//				hzc->Fill(20.0/(crossingPosc));
						}
						lastCrossingc = crossingPosc;
						Cross5=1;
					}
					else if(Cross6==0)
					{
						crossingPosc = tempSqTc + 1.0/(-tempSqVc + sqvoltsc)*(-tempSqTc + sqtimesc)*(meanSqc - tempSqVc) - firstCrossingc; 
						if(crossingPosc>22.0 && crossingPosc<28.0){ 
			//				meanCrossingsBevSevc[1] += 25.0/(crossingPosc);
			//				countCrossingc++;
			//				hzc->Fill(25.0/(crossingPosc));
						}
						lastCrossingc = crossingPosc;
						Cross6=1;
					}
					else if(Cross7==0)
					{
						crossingPosc = tempSqTc  + 1.0/(-tempSqVc + sqvoltsc)*(-tempSqTc + sqtimesc)*(meanSqc - tempSqVc) - firstCrossingc;
						if(crossingPosc>27.0 && crossingPosc<33.0){ 
						hmc->Fill( 1.0/(-tempSqVc + sqvoltsc)*(-tempSqTc + sqtimesc)*(meanSqc - tempSqVc));
							meanCrossingsBevSevc[1] += 30.0/(crossingPosc);
							countCrossingc++;
							hzc->Fill(30.0/crossingPosc);
						}
						lastCrossingc = crossingPosc;
						Cross7=1;
//						cerr << " Last: " << crossingPosc;
					}
					else if(Cross8==0)
					{
						crossingPosc = tempSqTc + 1.0/(-tempSqVc + sqvoltsc)*(-tempSqTc + sqtimesc)*(meanSqc - tempSqVc) - firstCrossingc; 
//						meanCrossingsBevSevc[1] += 35.0/(crossingPosc);
//						countCrossingc++;
//						hzc->Fill(35.0/(crossingPosc));
						lastCrossingc = crossingPosc;
						Cross8=1;
					}
				}
				tempSqVc2 = tempSqVc;
				tempSqVc=sqvoltsc;
				tempSqTc=sqtimesc;
				if(z%64==63){firstc=0;secondCross=0;thirdCross=0;Cross4=0;Cross5=0;Cross6=0;Cross7=0;Cross8=0;}
			}
		}
		firstc = 0;
		firstfirstc=0;
		lastCrossingc =0.0;
		secondCross=0;
		thirdCross = 0;	
		Cross4 = 0;	
		Cross5 = 0;	
		Cross6 = 0;	
		Cross7 = 0;	
		Cross8 = 0;	

//		cerr << endl;

		for(int z=64;z<grSqcO->GetN();z++)
		{
			grSqcO->GetPoint(z,sqtimesc,sqvoltsc);
			if(z<grSqc->GetN()-1)grSqc->GetPoint(z+1,dummyT,tempSqVc3);
			if(z==64){ tempSqVc2=sqvoltsc;tempSqVc=sqvoltsc;tempSqTc=sqtimesc;}	
			else
			{
				if( /*sqtimesc>(0.1 + 40.0*(z/64)) &&*/ z>65 && ( (meanSqc<tempSqVc3 && meanSqc>tempSqVc2 && meanSqc>tempSqVc && meanSqc<sqvoltsc) || (meanSqc>tempSqVc3 && meanSqc<tempSqVc2 && meanSqc<tempSqVc && meanSqc>sqvoltsc) ) )
				{
					if(firstfirstc==0)
					{
//						if(meanSqc>tempSqVc && meanSqc<sqvoltsc){
						hffcO->Fill( tempSqTc + 1.0/(-tempSqVc + sqvoltsc)*(-tempSqTc + sqtimesc)*(meanSqc - tempSqVc));
						firstfirstc=1;
//						}
					}
	
					if(firstc==0)
					{
						if(stationID==2||(meanSqc>tempSqVc && meanSqc<sqvoltsc)){
						firstCrossingc = tempSqTc + 1.0/(-tempSqVc + sqvoltsc)*(-tempSqTc + sqtimesc)*(meanSqc - tempSqVc);
						hffcO->Fill( tempSqTc + 1.0/(-tempSqVc + sqvoltsc)*(-tempSqTc + sqtimesc)*(meanSqc - tempSqVc));
//						hzcO->Fill(0.0);
						meanCrossingsBevSoddc[0] = firstCrossingc;
						lastCrossingc = 0.0;
						firstc=1;
						}
					}		
					else if(secondCross==0)
					{
						crossingPosc = tempSqTc + 1.0/(-tempSqVc + sqvoltsc)*(-tempSqTc + sqtimesc)*(meanSqc - tempSqVc) - firstCrossingc; 
						if(crossingPosc>2.0 && crossingPosc<8.0){ 
			//				meanCrossingsBevSoddc[1] += 5.0/(crossingPosc - lastCrossingc);
			//				countCrossingcO++;
			//				hzcO->Fill(5.0/(crossingPosc - lastCrossingc));
						}
						lastCrossingc = crossingPosc;
						secondCross = 1;
					}
					else if(thirdCross==0)
					{
						crossingPosc = tempSqTc + 1.0/(-tempSqVc + sqvoltsc)*(-tempSqTc + sqtimesc)*(meanSqc - tempSqVc) - firstCrossingc; 
						if(crossingPosc>7.0 && crossingPosc<13.0){ 
			//				meanCrossingsBevSoddc[1] += 10.0/(crossingPosc);
			//				countCrossingcO++;
			//				hzcO->Fill(10.0/(crossingPosc));
						}
						lastCrossingc = crossingPosc;
						thirdCross=1;
					}
					else if(Cross4==0)
					{
						crossingPosc = tempSqTc + 1.0/(-tempSqVc + sqvoltsc)*(-tempSqTc + sqtimesc)*(meanSqc - tempSqVc) - firstCrossingc; 
						if(crossingPosc>12.0 && crossingPosc<18.0){ 
			//				meanCrossingsBevSoddc[1] += 15.0/(crossingPosc);
			//				countCrossingcO++;
			//				hzcO->Fill(15.0/(crossingPosc));
						}
						lastCrossingc = crossingPosc;
						Cross4=1;
					}
					else if(Cross5==0)
					{
						crossingPosc = tempSqTc + 1.0/(-tempSqVc + sqvoltsc)*(-tempSqTc + sqtimesc)*(meanSqc - tempSqVc) - firstCrossingc; 
						if(crossingPosc>17.0 && crossingPosc<23.0){ 
			//				meanCrossingsBevSoddc[1] += 20.0/(crossingPosc);
			//				countCrossingcO++;
			//				hzcO->Fill(20.0/(crossingPosc));
						}
						lastCrossingc = crossingPosc;
						Cross5=1;
					}
					else if(Cross6==0)
					{
						crossingPosc = tempSqTc + 1.0/(-tempSqVc + sqvoltsc)*(-tempSqTc + sqtimesc)*(meanSqc - tempSqVc) - firstCrossingc; 
						if(crossingPosc>22.0 && crossingPosc<28.0){ 
			//				meanCrossingsBevSoddc[1] += 25.0/(crossingPosc);
			//				countCrossingcO++;
			//				hzcO->Fill(25.0/(crossingPosc));
						}
						lastCrossingc = crossingPosc;
						Cross6=1;
					}
					else if(Cross7==0)
					{
						crossingPosc = tempSqTc + 1.0/(-tempSqVc + sqvoltsc)*(-tempSqTc + sqtimesc)*(meanSqc - tempSqVc) - firstCrossingc; 
						if(crossingPosc>27.0 && crossingPosc<33.0){ 
						hmcO->Fill( 1.0/(-tempSqVc + sqvoltsc)*(-tempSqTc + sqtimesc)*(meanSqc - tempSqVc));
							meanCrossingsBevSoddc[1] += 30.0/(crossingPosc);
							countCrossingcO++;
							hzcO->Fill(30.0/crossingPosc);
						}
						lastCrossingc = crossingPosc;
						Cross7=1;
					}
					else if(Cross8==0)
					{
						crossingPosc = tempSqTc + 1.0/(-tempSqVc + sqvoltsc)*(-tempSqTc + sqtimesc)*(meanSqc - tempSqVc) - firstCrossingc; 
//						meanCrossingsBevSoddc[1] += 35.0/(crossingPosc);
//						countCrossingcO++;
//						hzcO->Fill(35.0/(crossingPosc));
						lastCrossingc = crossingPosc;
						Cross8=1;
					}
				}
				tempSqVc2 = tempSqVc;
				tempSqVc=sqvoltsc;
				tempSqTc=sqtimesc;
				if(z%64==63){firstc=0;secondCross=0;thirdCross=0;Cross4=0;Cross5=0;Cross6=0;Cross7=0;Cross8=0;}
			}
		}
		}//ifblockN%2==0
		delete grSqc;
		delete grSqcO;
		delete realAtriEvPtr;

//		loadBar(event, numEntries, 100, 50);
  	    pass_evts++;

	}//end loop all events

	meanSqCfreq[0] = (meanCrossingsBevSevc[1]/countCrossingc) ;
	meanSqCfreq[1] = (meanCrossingsBevSoddc[1]/countCrossingcO) ;

	cerr << "Lets see: " << countCrossingc << "   " << countCrossingcO << endl;
	cerr << "Lets see: " << meanCrossingsBevSevc[1] << "   " << meanCrossingsBevSoddc[1] << endl;
	cerr << "Results: " << meanSqCfreq[0] << "    " << meanSqCfreq[1] << endl;

	TCanvas *cffc = new TCanvas("cffc");
	hffc->Draw("");
	hffcO->SetLineColor(2);
	hffcO->Draw("same");
	cffc->Write();
	TCanvas *cmc = new TCanvas("cmc");
	hmc->Draw("");
	hmcO->SetLineColor(2);
	hmcO->Draw("same");
	cmc->Write();
	TCanvas *czc = new TCanvas("czc");
	hzc->Draw("");
	hzcO->SetLineColor(2);
	hzcO->Draw("same");
	czc->Write();
	fout->Close();

	cout << "All crossing stuff C: " << meanSqCfreq[0] << "    "  << meanSqCfreq[1] << endl;
	return 0;
}


//Due to some features of the waveforms, this module is specifically designed for waveforms of uncalibrated times, voltages and all 64 samples per block.
//It subtracts the average as well as some fixed offset, which appears due to the non-linearity of the waveforms.
double clearOffset(TGraph *gr, int stationId)
{
	double times, volts;
	double meanV = 0;
	double meanV2 = 0;
	double meanVfirst = 0;
	int meanCountFirst = 0;
	int meanCount1 = 0;
	int meanCount2 = 0;
	double rms = 0;
	double rms2 = 0;
	double fixedOffset = 0;
	double blockOffset[80] = {0};
	double blockRMS[80] = {0};
	int blockCount[80] = {0};

	if(stationId==2)fixedOffset=15;
	else fixedOffset=10;
	if(gr->GetN()>600){
	  for(int i=0;i<gr->GetN();i++)
	  {
	     gr->GetPoint(i,times,volts);
//	     if(i==0) cerr << "The first time is: " << times << endl;
	     if(times>64.0){
	        blockOffset[(int)(times-1)/64]+=volts;
	        blockRMS[(int)(times-1)/64]+=volts*volts;
	        blockCount[(int)(times-1)/64]++;
		if(i>=gr->GetN()/3*2){meanV+=volts;meanCount1++;rms+=volts*volts;}
	  	if(i<gr->GetN()/3){meanV2+=volts;meanCount2++;rms2+=volts*volts;}
	     }
	     else {meanVfirst+=volts;meanCountFirst++;}
	  }
	}
	else{
	  for(int i=0;i<gr->GetN();i++)
	  {
	     gr->GetPoint(i,times,volts);
//	     if(i==0) cerr << "The first time 2 is: " << times << endl;
	     if(times>64){
	        blockOffset[(int)(times-1)/64]+=volts;
	        blockRMS[(int)(times-1)/64]+=(volts*volts);
	        blockCount[(int)(times-1)/64]++;
		meanV+=volts;meanCount1++;rms+=volts*volts;
	  	meanV2+=volts;meanCount2++;rms2+=volts*volts;
	     }
	     else {meanVfirst+=volts;meanCountFirst++;}
	  }

	}

	rms = rms/meanCount1;
	rms2 = rms2/meanCount2;
	meanV = meanV/meanCount1;
	meanV2 = meanV2/meanCount2;
	meanVfirst = meanVfirst/meanCountFirst;	
	if(rms>rms2) meanV = meanV2;
	else meanV = meanV;

	for(int i=0;i<80;i++){
//		cerr << "The offset: " << blockOffset[i] << "  " <<blockCount[i] << " the average: " << meanV << endl;
		if(blockCount[i]>5){
			if(TMath::Abs(blockOffset[i]/blockCount[i] - meanV)>25.0 &&TMath::Sqrt(blockRMS[i]/blockCount[i] - blockOffset[i]*blockOffset[i]/(blockCount[i]*blockCount[i]))<90.0){
				//cerr << "The offset: " << blockOffset[i]/blockCount[i] << "The RMS: " << TMath::Sqrt(blockRMS[i]/blockCount[i] -  blockOffset[i]*blockOffset[i]/(blockCount[i]*blockCount[i])) << " the average: " << meanV << endl;
				 blockOffset[i] = blockOffset[i]/blockCount[i] - meanV; 
				//cerr << "This happens!" << endl;
				}
			else blockOffset[i] = 0.0;
		}
		else blockOffset[i] = 0.0;
	}



	for(int i=0;i<gr->GetN();i++)
	{
		gr->GetPoint(i,times,volts);
		if(times<=64)gr->SetPoint(i,times,volts - meanVfirst - fixedOffset);
		else gr->SetPoint(i,times,volts - meanV - fixedOffset - blockOffset[int(times)/64]);
	}
	return meanV;
}



/*
double clearOffsetErr(TGraphErrors *gr)
{
	double times, volts;
	double meanV = 0;
	double meanV2 = 0;
	double rms = 0;
	double rms2 = 0;
	for(int i=0;i<gr->GetN();i++)
	{
		gr->GetPoint(i,times,volts);
		if(i>=gr->GetN()/4*3){meanV+=volts/(gr->GetN()/4);rms+=volts*volts;}
		if(i<gr->GetN()/4){meanV2+=volts/(gr->GetN()/4);rms2+=volts*volts;}
	}

	if(rms>rms2) meanV = meanV2;
	else meanV = meanV;

	for(int i=0;i<gr->GetN();i++)
	{
		gr->GetPoint(i,times,volts);
		gr->SetPoint(i,times,volts - meanV- 11);
	}
	return meanV;
}
*/



void invertGraph(TGraph *gr)
{
	double t1, v1;
	for(int i=0;i<gr->GetN();i++)
	{
		gr->GetPoint(i,t1,v1);
		gr->SetPoint(i,t1,-v1);
	}
}




/*
TGraphErrors * calibrateElecChan(int elchan, UsefulAtriStationEvent * owncalEvPtr, double *TCorr, int stationId)
{
	TGraphErrors * NgrOwnCal = 0;
	TGraphErrors * VgrOwnCal2 = 0;
	TGraphErrors * VgrOwnCal = new TGraphErrors(100);
	TGraphErrors * SgrOwnCal = 0;
	TGraphErrors * SgrOwnCal3 = new TGraphErrors(100);
	TGraphErrors * halfS = 0;
	TGraphErrors * halfScut = new TGraphErrors(100);
	TGraphErrors * halfSSorted = 0;

	int SdtCount=0;
	NgrOwnCal = calibrateTime(elchan, owncalEvPtr, &TCorr[0], stationId);
	double meanV2 = clearOffset(NgrOwnCal, stationId);  
	VgrOwnCal2 = calibrateVolts(owncalEvPtr, elchan/8, elchan%8, NgrOwnCal, stationId);
	if(elchan%8<2){
		int nonumber = cutOverTimes(VgrOwnCal2, VgrOwnCal);
		SgrOwnCal = sortGraph(VgrOwnCal);
		SdtCount = selectSamples(SgrOwnCal, SgrOwnCal3);
		if(stationId==3 && elchan<2){invertGraph(SgrOwnCal3);}
		return SgrOwnCal3;
	}
	else{
		halfS = halfSamples(VgrOwnCal2, 1);
		cutOverTimes(halfS,halfScut,1);
		halfSSorted = sortGraph(halfScut);
		if(stationId==3 && elchan==2){invertGraph(halfSSorted);}
		return halfSSorted;
	}
		

}
*/


int timeSampleNumbers[4][8][2] = {{{0}}};
int timeIndex[4][8][2][64] = {{{{0}}}};
double timeConst[4][8][2][64] = {{{{0}}}};


void loadTimeConstants(int stationId)
{

	double tTime;
	int tIndex;
	int tNumber;
	int tChip, tChan, tCap;

	ifstream calRead;
	char calReadchar[60];
	sprintf(calReadchar,(calFilesDir + "/AraStation%dSampleTimingNew.txt").c_str(), stationId);
	cout << "Opening: " << calReadchar << endl;
	calRead.open(calReadchar);
	if(calRead.good()){
	for(int chi=0;chi<4;chi++)
	{
		for(int cha=0;cha<8;cha++)
		{
			calRead >> tChip >> tChan >> tCap >> tNumber;
			timeSampleNumbers[tChip][tChan][tCap] = tNumber;
			for(int ind=0;ind<tNumber;ind++)
			{
				calRead >> tIndex;
				timeIndex[tChip][tChan][tCap][ind] = tIndex;
			}
			calRead >> tChip >> tChan >> tCap >> tNumber;
			timeSampleNumbers[tChip][tChan][tCap] = tNumber;
			for(int ind=0;ind<tNumber;ind++)
			{
				calRead >> tTime;
				timeConst[tChip][tChan][tCap][ind] = tTime;
			}
			calRead >> tChip >> tChan >> tCap >> tNumber;
			timeSampleNumbers[tChip][tChan][tCap] = tNumber;
			for(int ind=0;ind<tNumber;ind++)
			{
				calRead >> tIndex;
				timeIndex[tChip][tChan][tCap][ind] = tIndex;
			}
			calRead >> tChip >> tChan >> tCap >> tNumber;
			timeSampleNumbers[tChip][tChan][tCap] = tNumber;
			for(int ind=0;ind<tNumber;ind++)
			{
				calRead >> tTime;
				timeConst[tChip][tChan][tCap][ind] = tTime;
			}
		}
	}
	}else{
		cout << "Couldn't open timing calibration file!" << endl;
	}
	calRead.close();
}











TGraph *calibrateTimeAltern(int elchan, UsefulAtriStationEvent * owncalEvPtr, double tempCorr, int stationId)
{

	int fullSamples = 0;
	double sort_times[4096];
	double sort_volts[4096];
	int chip = elchan/8;
	int grBlockN = 0;
	int block = 0;
	double times, volts;
	double newVolts;
	double corrTime = 0;
//	cerr << "Start calibration!" << endl;
	TGraph *grOwnCal = owncalEvPtr->getGraphFromElecChan(elchan);
//	cerr << "Read graph with: " << grOwnCal->GetN() << " points." << endl;
	
	clearOffset(grOwnCal, stationId);
//	cerr << "Offset cleared!" << endl;
	int Nsamples = grOwnCal->GetN();
	double sampleTime[128] = {0};
	int Nsamples2 = 0;
	int Nsamples1 = 0;

//	cerr << "We are calibrating channel " << elchan << " sitting on chip: " << elchan/8 << " with the channel number: " << elchan%8 << endl;

	for(int bl=0;bl<grOwnCal->GetN()/64;bl++)
	{
		block = owncalEvPtr->blockVec[bl*4].getBlock();
//		cerr << "the blocks in new cal: " << block;
		if(block%2==0)
		{
			Nsamples2 = timeSampleNumbers[elchan/8][elchan%8][0];	
//			cerr << "  With samples: " << Nsamples2;
			for(int j= 0;j<Nsamples2;j++)
			{
				corrTime = bl*20.0 + timeConst[elchan/8][elchan%8][0][j]*tempCorr;
				sort_times[fullSamples] = corrTime;
				grOwnCal->GetPoint(timeIndex[elchan/8][elchan%8][0][j]+bl*64,times,volts);
//				cerr << " " <<  timeIndex[elchan/8][elchan%8][0][j];
				grBlockN = block;
				while(calConst[chip][elchan%8][grBlockN][timeIndex[elchan/8][elchan%8][0][j]][8]>1.0) grBlockN = (grBlockN - 2 + 512)%512;		
				if(TMath::Abs(volts)<400)
				{
					if(volts>0)
					{
						 newVolts = calConst[chip][elchan%8][grBlockN][timeIndex[elchan/8][elchan%8][0][j]][6] 
						+ (volts - calConst[chip][elchan%8][grBlockN][timeIndex[elchan/8][elchan%8][0][j]][7])
							*calConst[chip][elchan%8][grBlockN][timeIndex[elchan/8][elchan%8][0][j]][0] 
						+ pow((volts - calConst[chip][elchan%8][grBlockN][timeIndex[elchan/8][elchan%8][0][j]][7]), 2)
							*calConst[chip][elchan%8][grBlockN][timeIndex[elchan/8][elchan%8][0][j]][1] 
						+ pow((volts - calConst[chip][elchan%8][grBlockN][timeIndex[elchan/8][elchan%8][0][j]][7]), 3)
							*calConst[chip][elchan%8][grBlockN][timeIndex[elchan/8][elchan%8][0][j]][2]; 
					}
					else
					{
						 newVolts = calConst[chip][elchan%8][grBlockN][timeIndex[elchan/8][elchan%8][0][j]][6] 
						+ (volts - calConst[chip][elchan%8][grBlockN][timeIndex[elchan/8][elchan%8][0][j]][7])
							*calConst[chip][elchan%8][grBlockN][timeIndex[elchan/8][elchan%8][0][j]][3] 
						+ pow((volts - calConst[chip][elchan%8][grBlockN][timeIndex[elchan/8][elchan%8][0][j]][7]), 2)
							*calConst[chip][elchan%8][grBlockN][timeIndex[elchan/8][elchan%8][0][j]][4] 
						+ pow((volts - calConst[chip][elchan%8][grBlockN][timeIndex[elchan/8][elchan%8][0][j]][7]), 3)
							*calConst[chip][elchan%8][grBlockN][timeIndex[elchan/8][elchan%8][0][j]][5]; 
					}
				}
				else
				{
					if(volts>0)
					{
						 newVolts = calConstHigh[chip][elchan%8][grBlockN][timeIndex[elchan/8][elchan%8][0][j]][0] 
							+ volts*calConstHigh[chip][elchan%8][grBlockN][timeIndex[elchan/8][elchan%8][0][j]][1];
					}
					else
					{
					 	newVolts = calConstHigh[chip][elchan%8][grBlockN][timeIndex[elchan/8][elchan%8][0][j]][2] 
							+ volts*calConstHigh[chip][elchan%8][grBlockN][timeIndex[elchan/8][elchan%8][0][j]][3];
					}
				}
				sort_volts[fullSamples] = newVolts;
				fullSamples++;
			}
		}
		else
		{
			Nsamples2 = timeSampleNumbers[elchan/8][elchan%8][1];	
//			cerr << "  With samples: " << Nsamples2;
			for(int j=0;j<Nsamples2;j++)
			{
				corrTime = bl*20.0 + timeConst[elchan/8][elchan%8][1][j]*tempCorr - 20.0;
				sort_times[fullSamples] = corrTime;
				grOwnCal->GetPoint(timeIndex[elchan/8][elchan%8][1][j]+bl*64,times,volts);
//				cerr << " " <<  timeIndex[elchan/8][elchan%8][1][j];
				grBlockN = block;
				while(calConst[chip][elchan%8][grBlockN][timeIndex[elchan/8][elchan%8][1][j]][8]>1.0) grBlockN = (grBlockN - 2 + 512)%512;		
				if(TMath::Abs(volts)<400)
				{
					if(volts>0)
					{
						 newVolts = calConst[chip][elchan%8][grBlockN][timeIndex[elchan/8][elchan%8][1][j]][6] 
						+ (volts - calConst[chip][elchan%8][grBlockN][timeIndex[elchan/8][elchan%8][1][j]][7])
							*calConst[chip][elchan%8][grBlockN][timeIndex[elchan/8][elchan%8][1][j]][0] 
						+ pow((volts - calConst[chip][elchan%8][grBlockN][timeIndex[elchan/8][elchan%8][1][j]][7]), 2)
							*calConst[chip][elchan%8][grBlockN][timeIndex[elchan/8][elchan%8][1][j]][1] 
						+ pow((volts - calConst[chip][elchan%8][grBlockN][timeIndex[elchan/8][elchan%8][1][j]][7]), 3)
							*calConst[chip][elchan%8][grBlockN][timeIndex[elchan/8][elchan%8][1][j]][2]; 
					}
					else
					{
						 newVolts = calConst[chip][elchan%8][grBlockN][timeIndex[elchan/8][elchan%8][1][j]][6] 
						+ (volts - calConst[chip][elchan%8][grBlockN][timeIndex[elchan/8][elchan%8][1][j]][7])
							*calConst[chip][elchan%8][grBlockN][timeIndex[elchan/8][elchan%8][1][j]][3] 
						+ pow((volts - calConst[chip][elchan%8][grBlockN][timeIndex[elchan/8][elchan%8][1][j]][7]), 2)
							*calConst[chip][elchan%8][grBlockN][timeIndex[elchan/8][elchan%8][1][j]][4] 
						+ pow((volts - calConst[chip][elchan%8][grBlockN][timeIndex[elchan/8][elchan%8][1][j]][7]), 3)
							*calConst[chip][elchan%8][grBlockN][timeIndex[elchan/8][elchan%8][1][j]][5]; 
					}
				}
				else
				{
					if(volts>0)
					{
						 newVolts = calConstHigh[chip][elchan%8][grBlockN][timeIndex[elchan/8][elchan%8][1][j]][0] 
							+ volts*calConstHigh[chip][elchan%8][grBlockN][timeIndex[elchan/8][elchan%8][1][j]][1];
					}
					else
					{
					 	newVolts = calConstHigh[chip][elchan%8][grBlockN][timeIndex[elchan/8][elchan%8][1][j]][2] 
							+ volts*calConstHigh[chip][elchan%8][grBlockN][timeIndex[elchan/8][elchan%8][1][j]][3];
					}
				}
				sort_volts[fullSamples] = newVolts;
				fullSamples++;
			}
		}
//		cerr << endl;
	}

//	cerr << "The number of samples is: " << fullSamples << endl;

	delete grOwnCal;	
	TGraph * NgrOwnCal = new TGraph();
	int countGr=0;
	for(int k=0;k<fullSamples;k++)
	{
		if(sort_times[k]>20.0){
		NgrOwnCal->SetPoint(countGr,sort_times[k],sort_volts[k]);
		countGr++;
		}
		
	}
	NgrOwnCal->Sort();
	return NgrOwnCal;
}
