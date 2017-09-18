/*
__kernel void arithmetic( __global float *input ){

   float4 result;
   float4 x = vload4(0, input);
   float4 y = (float4)(3.f, 3.f, 3.f,3.f);
   //result = fmod(x, y); //x%y for float
   result = remainder(x, y); // x - n * y, where n is integer closest to x/y
   vstore4( result, 0, input);

}
*/
#define C_INV 3.34
/*
__kernel void getRecoDelays( __global float *delays, __global float *antLoc, __global float *zen, __global float *azi){

   int gid0=get_global_id(0);
   int gid1=get_global_id(1);
   int gid2=get_global_id(2);

   int nAzi=get_global_size(1);
   int nAnt=get_global_size(2);

   float zenRad=zen[gid0]*M_PI/180.f;
   float aziRad=azi[gid1]*M_PI/180.f;
   float dir_x = sin(zenRad)*cos(aziRad);
   float dir_y = sin(zenRad)*sin(aziRad);
   float dir_z = cos(zenRad);
   delays[gid0*nAzi*nAnt + gid1*nAnt + gid2] = -C_INV*( dir_x * antLoc[3*gid2]
                                                      + dir_y * antLoc[3*gid2+1]
                                                      + dir_z * antLoc[3*gid2+2]
                                                      ) ;
}
*/
__kernel void wfPwr( __global float *pwr, __global float *beam_r, __global float *beam_c, int nSamp){

   int gid0=get_global_id(0);

   pwr[gid0]=0.f;
   for(int i=0; i<nSamp; i++) pwr[gid0] += (beam_r[gid0*nSamp + i] * beam_r[gid0*nSamp + i] +
                                            beam_c[gid0*nSamp + i] * beam_c[gid0*nSamp + i] );
   if(pwr[gid0]<0.f) printf("pwr: %f\t", pwr[gid0]);
}

__kernel void sumWf( __global float *beam_r, __global float *beam_c,
                     __global float *shiftedRWf, __global float *shiftedCWf,
                     int nAnt){

   int gid0 = get_global_id(0);
   int gid1 = get_global_id(1);
   //size_t dim_0 = get_global_size(0); //number of waveforms collapsed together
   //int nSamp=768;
   int nDir  = get_global_size(0);
   int nSamp = get_global_size(1);

   beam_r[gid0*nSamp+gid1]=0.f;
   beam_c[gid0*nSamp+gid1]=0.f;
   for(int i=0; i<nAnt; i++){
   //beam[gid_1] = (intensity[gid_0*(nSamp) + gid_1] + intensity[(!gid_0)*(nSamp) + gid_1])/2;
   //beam_r[gid0*nSamp+gid1] += shiftedWf[gid0*nAnt*nSamp*2+i*nSamp*2+2*gid1];
   //beam_c[gid0*nSamp+gid1] += shiftedWf[gid0*nAnt*nSamp*2+i*nSamp*2+2*gid1+1];
   beam_r[gid0*nSamp+gid1] += shiftedRWf[gid0*nAnt*nSamp + i*nSamp + gid1];
   beam_c[gid0*nSamp+gid1] += shiftedCWf[gid0*nAnt*nSamp + i*nSamp + gid1];
  //atomic_add(beam[gid_1], intensity[gid_0*nSamp+gid_1]);
   }
   beam_r[gid0*nSamp+gid1] /= (float)nAnt;
   beam_c[gid0*nSamp+gid1] /= (float)nAnt;
   //beam[gid0*nSamp+gid1] *= beam[gid0*exNSamp+gid1];
}

__kernel void shiftWf( __global float *shiftedRWf, __global float *shiftedCWf,
                       __global float *delays,
                       __global float *intensity_r, __global float *intensity_c,
                       float wInt){

   int gid0 = get_global_id(0);
   int gid1 = get_global_id(1);
   int gid2 = get_global_id(2);
   int nDir = get_global_size(0);
   int nAnt = get_global_size(1);
   int nSamp = get_global_size(2);

   if( delays[gid0*nAnt+gid1] > -1e9 ) //have raytrace solution
   {
   float delaysBin = delays[gid0*nAnt+gid1]/wInt;
   //for(int i=0; i<nDir*nAnt; i++) delaysBin[i]=(int)(delays[gid0*nAnt+i]/wInt);
   //delaysBin[0] = delays[0]/2;
   //delaysBin[1] = delays[1]/2;

   /* F{ g(t+a) } = exp(i*2pi*f*a) * G(f) */
   float phi = 2.f*M_PI*(float)gid2*delaysBin/(float)(nSamp);
   float rot[2] = {cos(phi), sin(phi)};
   float intensity2[2];

   if( gid2 < (nSamp/2 + 1) ){

      intensity2[0] = intensity_r[gid1*(nSamp/2+1)+gid2];
      intensity2[1] = intensity_c[gid1*(nSamp/2+1)+gid2];

   } else {

      intensity2[0] =  intensity_r[gid1*(nSamp/2+1)+nSamp-gid2];
      intensity2[1] =  intensity_c[gid1*(nSamp/2+1)+nSamp-gid2] * -1.f;

   }

   shiftedRWf[gid0*nAnt*nSamp + gid1*nSamp + gid2] = intensity2[0]*rot[0] - intensity2[1]*rot[1];
   shiftedCWf[gid0*nAnt*nSamp + gid1*nSamp + gid2] = intensity2[1]*rot[0] + intensity2[0]*rot[1];

   } else //No raytrace solution
   {

   shiftedRWf[gid0*nAnt*nSamp + gid1*nSamp + gid2] = 0.f;
   shiftedCWf[gid0*nAnt*nSamp + gid1*nSamp + gid2] = 0.f;

   }
   //shiftedRWf[gid0*nAnt*nSamp + gid1*nSamp + gid2] = 0.1;
   //shiftedCWf[gid0*nAnt*nSamp + gid1*nSamp + gid2] = 0.1;

   /*
   if( gid2%2 == 0 ) //stores real part
   shiftedWf[gid0*nAnt*interlvOutputSize+gid1*interlvOutputSize+gid2] = intensity2[0]*rot[0] - intensity2[1]*rot[1];
   else              //stores imaginary part
   shiftedWf[gid0*nAnt*interlvOutputSize+gid1*interlvOutputSize+gid2] = intensity2[1]*rot[0] + intensity2[0]*rot[1];
   */
   /*
   if( (gid2 + delaysBin >= 0)  &&
       (gid2 + delaysBin <=(exNSamp-1))
   )
        shiftedWf[gid0*nAnt*exNSamp+gid1*exNSamp+gid2] = intensity[gid1*exNSamp + gid2 + delaysBin];
   else shiftedWf[gid0*nAnt*exNSamp+gid1*exNSamp+gid2] = 0.f;
   */
   //shift[gid_0*(nSamp)+gid_1] = intensity[gid_0*(nSamp) + gid_1 ];

}


__kernel void xCorrWf( __global float *xCorr_r, __global float *xCorr_c,
                       __global float *intensity_r, __global float *intensity_c,
                       int nAnt){

   int gid0 = get_global_id(0); //which baseline
   int gid1 = get_global_id(1); //which freq bin
   int planarHermOutputSize = get_global_size(1);

   int ant1 = gid0 / nAnt;
   int ant2 = gid0 % nAnt;

   xCorr_r[gid0*planarHermOutputSize + gid1] = intensity_r[ant1*planarHermOutputSize + gid1] * intensity_r[ant2*planarHermOutputSize + gid1]
                                             + intensity_c[ant1*planarHermOutputSize + gid1] * intensity_c[ant2*planarHermOutputSize + gid1] ;
   xCorr_c[gid0*planarHermOutputSize + gid1] = intensity_r[ant1*planarHermOutputSize + gid1] * intensity_c[ant2*planarHermOutputSize + gid1] * -1.f
                                             + intensity_c[ant1*planarHermOutputSize + gid1] * intensity_r[ant2*planarHermOutputSize + gid1] ;

}

/*
 * This kernel is used to x-corr FFT wavefoms with no Hermitian symmetry, which may be due to operations in the frequency domain. Eg, bandpass filtering
 */
/*
__kernel void xCorrAsymWf( __global float *xCorr_r, __global float *xCorr_c,
                                  __global float *full_intensity_r, __global float *full_intensity_c,
                                  int nAnt){

   int gid0 = get_global_id(0); //which baseline
   int gid1 = get_global_id(1); //which freq bin
   int planarHermOutputSize = get_global_size(1);

   int ant1 = gid0 / nAnt;
   int ant2 = gid0 % nAnt;

   xCorr_r[gid0*planarHermOutputSize + gid1] = intensity_r[ant1*planarHermOutputSize + gid1] * intensity_r[ant2*planarHermOutputSize + gid1]
                                             + intensity_c[ant1*planarHermOutputSize + gid1] * intensity_c[ant2*planarHermOutputSize + gid1] ;
   xCorr_c[gid0*planarHermOutputSize + gid1] = intensity_r[ant1*planarHermOutputSize + gid1] * intensity_c[ant2*planarHermOutputSize + gid1] * -1.f
                                             + intensity_c[ant1*planarHermOutputSize + gid1] * intensity_r[ant2*planarHermOutputSize + gid1] ;

}
*/

__kernel void computeXCorrCoef( __global float *Cij, __global float *xCorrTime,
                                __global float *delays, __global float *sqrtWfPwr,
                                float wInt, int nAnt, int nSamp){

   int gid0 = get_global_id(0); //which layer
   int gid1 = get_global_id(1); //which direction
   int gid2 = get_global_id(2); //which baseline
   int nLayer    = get_global_size(0);
   int nDir      = get_global_size(1);
   int nBaseline = get_global_size(2);
   int ant1 = gid2 / nAnt;
   int ant2 = gid2 % nAnt;

   if( (gid0*nDir*nAnt + gid1*nAnt + ant1) >= nLayer*nDir*nBaseline || (gid0*nDir*nAnt + gid1*nAnt + ant2) >= nLayer*nDir*nBaseline ) printf("****** Warning! out of range!!! **** nLayer: %d nDir: %d nBaseline: %d ant1: %d ant2: %d\n",nLayer, nDir, nBaseline, gid0*nDir*nAnt + gid1*nAnt + ant1 -  nLayer*nDir*nBaseline, gid0*nDir*nAnt + gid1*nAnt + ant2 -  nLayer*nDir*nBaseline);
   if( sqrtWfPwr[ant1] != 0.f && sqrtWfPwr[ant2] != 0.f ){
   if(delays[gid0*nDir*nAnt + gid1*nAnt + ant1] > -1e9 && delays[gid0*nDir*nAnt + gid1*nAnt + ant2] > -1e9 ){
   int shiftBin = (delays[gid0*nDir*nAnt + gid1*nAnt + ant1] - delays[gid0*nDir*nAnt + gid1*nAnt + ant2]) / wInt;
   /* The default circular  FFT of cross-correlation outputs in wrap-around order */
   /* The correlation at -i is in r_N-i */
   if( shiftBin < 0 ) shiftBin += nSamp;
   Cij[gid0*nDir*nBaseline + gid1*nBaseline + gid2] = xCorrTime[gid2*nSamp + shiftBin] / (sqrtWfPwr[ant1] * sqrtWfPwr[ant2]);
   } else {
   //printf("Not both delays exists!\n");
   //if(delays[gid0*nAnt + ant1] != delays[gid0*nAnt + ant2] )
   //   if(delays[gid0*nAnt + ant1] > 0 || delays[gid0*nAnt + ant2] >0) printf("ant1 delay: %f ant2 delay: %f\n",delays[gid0*nAnt + ant1],delays[gid0*nAnt + ant2]);
   Cij[gid0*nDir*nBaseline + gid1*nBaseline + gid2] = 0.f;
   }
   } else {
   Cij[gid0*nDir*nBaseline + gid1*nBaseline + gid2] = 0.f;
   }

}

__kernel void computeCoherence( __global float *M, __global float *Cij,
                                int nBaseline){

   int gid0 = get_global_id(0); //which layer
   int gid1 = get_global_id(1); //which direction
   int nDir = get_global_size(1);
   int nAnt = (int)sqrt((float)nBaseline);
   //printf("nAnt in computeCoherence: %d\n", nAnt);

   M[gid0*nDir + gid1] = 0.f;
  for(int i=0; i<nAnt; i++){
     for(int j=i+1; j<nAnt; j++){
      M[gid0*nDir + gid1] += Cij[gid0*nDir*nBaseline + gid1*nBaseline + i*nAnt + j];
      }
   }
   //M[gid0] = sqrt(M[gid0]*M[gid0]);

}

__kernel void computeNormalizedCoherence( __global float *M, __global float *Cij,
                                          int nBaseline){

   int gid0 = get_global_id(0); //which layer
   int gid1 = get_global_id(1); //which direction
   int nDir = get_global_size(1);
   int nAnt = (int)sqrt((float)nBaseline);
   //printf("nAnt in computeCoherence: %d\n", nAnt);

   M[gid0*nDir + gid1] = 0.f;
  for(int i=0; i<nAnt; i++){
     for(int j=i+1; j<nAnt; j++){
      M[gid0*nDir + gid1] += Cij[gid0*nDir*nBaseline + gid1*nBaseline + i*nAnt + j];
      }
   }
   //M[gid0] = sqrt(M[gid0]*M[gid0]);
   M[gid0*nDir + gid1] /= (float)((nAnt*(nAnt-1))/2); //Normalize M by the number of baselines summed.
                                                      //This makes the comparison among events with different numbers of reco channel
                                                      //more natural. 16.15.16
}
/*
__kernel void bandStopFilter( __global float *full_intensity_r, __global float *full_intensity_c,
                              __global float *intensity_r, __global float *intensity_c,
                              float freqBin, float lowFreq, float highFreq){

   int gid0 = get_global_id(0);
   int gid1 = get_global_id(1);
   int nAnt = get_global_size(0);
   int nSamp = get_global_size(1);
   //int planarHermOutputSize = get_global_size(1);

   int planarHermOutputSize = 1 +  nSamp/2;
   //int nSamp = 2 * (planarHermOutputSize - 1);
   float freq = freqBin * gid1;

   if( gid1 < planarHermOutputSize ){

      if( freq > lowFreq && freq < highFreq )
      full_intensity_r[gid0*nSamp + gid1] = full_intensity_c[gid0*nSamp + gid1] = 0.f;
      //intensity_r[gid0*planarHermOutputSize + gid1] = intensity_c[gid0*planarHermOutputSize + gid1]= 0.f;
      else{
      full_intensity_r[gid0*nSamp + gid1] = intensity_r[gid0*planarHermOutputSize + gid1];
      full_intensity_c[gid0*nSamp + gid1] = intensity_c[gid0*planarHermOutputSize + gid1];
      }

   } else {

      if( freq > lowFreq && freq < highFreq )
      full_intensity_r[gid0*nSamp + gid1] = full_intensity_c[gid0*nSamp + gid1] = 0.f;
      else{
      full_intensity_r[gid0*nSamp + gid1] = intensity_r[gid0*planarHermOutputSize + nSamp - gid1];
      full_intensity_c[gid0*nSamp + gid1] = intensity_c[gid0*planarHermOutputSize + nSamp - gid1]*-1.f;
      }

   }//end of else

}
*/
__kernel void bandPassFilter( //__global float *full_intensity_r, __global float *full_intensity_c,
                              __global float *intensity_r, __global float *intensity_c,
                              float freqBin, float lowFreq, float highFreq){

   int gid0 = get_global_id(0);
   int gid1 = get_global_id(1);
   int nAnt = get_global_size(0);
   //int nSamp = get_global_size(1);
   int planarHermOutputSize = get_global_size(1);

   //int planarHermOutputSize = 1 +  nSamp/2;
   //int nSamp = 2 * (planarHermOutputSize - 1);
   float freq = freqBin * gid1;

//   if( gid1 < planarHermOutputSize ){

      if( freq < lowFreq || freq > highFreq )
      intensity_r[gid0*planarHermOutputSize + gid1] = intensity_c[gid0*planarHermOutputSize + gid1] = 0.f;
      //intensity_r[gid0*planarHermOutputSize + gid1] = intensity_c[gid0*planarHermOutputSize + gid1]= 0.f;
/*
      else{
      full_intensity_r[gid0*nSamp + gid1] = intensity_r[gid0*planarHermOutputSize + gid1];
      full_intensity_c[gid0*nSamp + gid1] = intensity_c[gid0*planarHermOutputSize + gid1];
      }
*/
/*
   } else {

      if( freq < lowFreq || freq > highFreq )
      full_intensity_r[gid0*nSamp + gid1] = full_intensity_c[gid0*nSamp + gid1] = 0.f;
      else{
      full_intensity_r[gid0*nSamp + gid1] = intensity_r[gid0*planarHermOutputSize + nSamp - gid1];
      full_intensity_c[gid0*nSamp + gid1] = intensity_c[gid0*planarHermOutputSize + nSamp - gid1]*-1.f;
      }

   }//end of else
*/
}

__kernel void getMaxPixInfoEachLayer( __global int *maxPixIdx, __global float *maxPixCoherence,
                                      __global float *M,
                                      int nDir){

   int gid0 = get_global_id(0); //which layer
   //int nLayer = get_global_size(0);

   int maxIdx;
   float max = 0.f;

   for(int i=0; i<nDir; i++){

   if(M[gid0*nDir + i] > max){
      max = M[gid0*nDir + i];
      maxIdx = gid0*nDir + i;
    }
    }

   maxPixIdx[gid0] = maxIdx;
   maxPixCoherence[gid0] = max;

}
