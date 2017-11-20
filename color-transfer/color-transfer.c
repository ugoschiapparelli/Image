/**
 * @file color-transfert
 * @brief transfert color from source image to target image.
 *        Method from Reinhard et al. :
 *        Erik Reinhard, Michael Ashikhmin, Bruce Gooch and Peter Shirley,
 *        'Color Transfer between Images', IEEE CGA special issue on
 *        Applied Perception, Vol 21, No 5, pp 34-41, September - October 2001
 */

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <bcl.h>

#define D 3

static float RGB2LMS[D][D] = {
  {0.3811, 0.5783, 0.0402},
  {0.1967, 0.7244, 0.0782},
  {0.0241, 0.1288, 0.8444}
};

static float LMS2RGB[D][D] = {
  {4.4679, -3.5873, 0.1193},
  {-1.2186, 2.3809, -0.1624},
  {0.0497, -0.2439, 1.2045}
};

static float LMS2LAB[D][D] = {
  {0.5774, 0.5774, 0.5774},
  {0.4082, 0.4082, -0.8165},
  {0.7071, -0.7071, 0.f}
};

static float LAB2LMS[D][D] = {
  {0.5774, 0.4082, 0.7071},
  {0.5774, 0.4082, -0.7071},
  {0.5774, -0.8165, 0.f}
};

struct Stats {
  float meanLambda;
  float meanAlpha;
  float meanBeta;

  float sigmaLambda;
  float sigmaAlpha;
  float sigmaBeta;
};


void multD(float matrix[D][D], float *color) {
  float result[D];
  for (size_t i = 0; i < D; i++) {
    float r = 0;
    for (size_t k = 0; k < D; k++) {
      r += color[k]*matrix[i][k];
    }
    result[i] = r;
  }
  for (size_t i = 0; i < D; i++) {
    color[i] = result[i];
  }
}

void rgb2lab(float *color) {
  for (size_t i = 0; i < D; i++)
    color[i] = 1.f*color[i] / 255;
  multD(RGB2LMS,color);
  //Log10
  for (size_t i = 0; i < D; i++) {
    if (color[i]!=0) {
      color[i] = log10(color[i]);
      }
    else {
      color[i] = log10(0.0001);
    }
  }
  multD(LMS2LAB,color);
}

void lab2rgb(float *color) {
  multD(LAB2LMS,color);
  //Pow10
  for (size_t i = 0; i < D; i++) {
    color[i] = pow(10,color[i]);
  }
  multD(LMS2RGB,color);
  for (size_t i = 0; i < D; i++) {
    if (color[i]>1) {
      color[i] = 255;
    } else if (color[i]<0) {
      color[i] = 0;
    } else {
      color[i] = color[i]*255;
    }
  }
}

void loadLAB(float * imgLab, unsigned short *img, int size) {
  for(int i = 0; i < size; i++){
    float color[D];
    for (int c = 0; c < D; c++)
      color[c] = *img++;
    rgb2lab(color);
    for (int c = 0; c < D; c++)
      *imgLab++ = color[c];
  }
}

void saveLAB(float * imgLab, unsigned short *img, int size) {
  for(int i = 0; i < size; i++){
    float color[D];
    for (int c = 0; c < D; c++)
      color[c] = *imgLab++;
    lab2rgb(color);
    for (int c = 0; c < D; c++)
      *img++ = color[c];
  }
}

struct Stats getStats(float *img, int size) {
  struct Stats stats;

  float sumLambda = 0.f;
  float sumAlpha = 0.f;
  float sumBeta = 0.f;
  int cpt = 0;
  //On fait la somme des l, a, b
  for(int i=0; i<size; i++){
    float color[D];
    for (int c = 0; c < D; c++)
      color[c] = *img++;

    sumLambda += color[0];
    sumAlpha += color[1];
    sumBeta += color[2];
    cpt += 1;
  }
  img -= size*3; //Going back to start
  stats.meanLambda = sumLambda/cpt;
  stats.meanAlpha = sumAlpha/cpt;
  stats.meanBeta = sumBeta/cpt;

  float tmpSigmaLambda = 0;
  float tmpSigmaAlpha = 0;
  float tmpSigmaBeta = 0;

  //On calcule l'ecart-type
  for(int i=0; i<size; i++){
    float color[D];
    for (int c = 0; c < D; c++)
      color[c] = *img++;

    float tmpLambda = color[0]-stats.meanLambda;
    float tmpAlpha = color[1]-stats.meanAlpha;
    float tmpBeta = color[2]-stats.meanBeta;

    tmpSigmaLambda += tmpLambda * tmpLambda;
    tmpSigmaAlpha += tmpAlpha * tmpAlpha;
    tmpSigmaBeta += tmpBeta * tmpBeta;
  }

  stats.sigmaLambda = sqrt(tmpSigmaLambda/cpt);
  stats.sigmaAlpha = sqrt(tmpSigmaAlpha/cpt);
  stats.sigmaBeta = sqrt(tmpSigmaBeta/cpt);

  return stats;
}

void applyChange(float * img, int size, struct Stats statsImt, struct Stats statsIms) {
  for(int i=0; i<size; i++){
    float color[D];
    for (int c = 0; c < D; c++)
      color[c] = *img++;
    img-=D; //Going back in array

    //l* = l - <l>
    color[0] = color[0] - statsImt.meanLambda;
    color[1] = color[1] - statsImt.meanAlpha;
    color[2] = color[2] - statsImt.meanBeta;

    //l'=slt/sls * l*
    color[0] = (statsIms.sigmaLambda/statsImt.sigmaLambda) * color[0];
    color[1] = (statsIms.sigmaAlpha/statsImt.sigmaAlpha) * color[1];
    color[2] = (statsIms.sigmaBeta/statsImt.sigmaBeta) * color[2];

    //l = l' + <l>
    color[0] = color[0] + statsIms.meanLambda;
    color[1] = color[1] + statsIms.meanAlpha;
    color[2] = color[2] + statsIms.meanBeta;

    for (int c = 0; c < D; c++) {
      *img++ = color[c];
    }
  }
}


static void process(char *ims, char *imt, char* imd){
  //Get stats for source images
  struct Stats statsIms;
  pnm p_ims = pnm_load(ims);
  int s_ims = pnm_get_height(p_ims)*pnm_get_width(p_ims);
  unsigned short *u_ims = pnm_get_image(p_ims);
  float *f_ims = (float *) malloc(3*s_ims*sizeof(float));
  loadLAB(f_ims, u_ims, s_ims);
  statsIms = getStats(f_ims, s_ims);

  //Get stats for target images
  struct Stats statsImt;
  pnm p_imt = pnm_load(imt);
  int s_imt = pnm_get_height(p_imt)*pnm_get_width(p_imt);
  unsigned short *u_imt = pnm_get_image(p_imt);
  float *f_imt = (float *) malloc(3*s_imt*sizeof(float));
  loadLAB(f_imt, u_imt, s_imt);
  statsImt = getStats(f_imt, s_imt);

  //Apply stats to target unsigned short array
  applyChange(f_imt, s_imt, statsImt, statsIms);

  //Create a new picture of target size
  pnm p_imd = pnm_new(pnm_get_width(p_imt), pnm_get_height(p_imt), PnmRawPpm);
  unsigned short *u_imd = pnm_get_image(p_imd);
  //Convert LAB array to RGB array
  saveLAB(f_imt, u_imd, s_imt);
  //Save picture
  pnm_save(p_imd, PnmRawPpm, imd);

  pnm_free(p_imd);
  pnm_free(p_imt);
  pnm_free(p_ims);
  free(f_imt);
  free(f_ims);
}

void usage (char *s){
  fprintf(stderr, "Usage: %s <ims> <imt> <imd> \n", s);
  exit(EXIT_FAILURE);
}

#define param 3
int main(int argc, char *argv[]){
  if (argc != param+1)
    usage(argv[0]);
  process(argv[1], argv[2], argv[3]);
  return EXIT_SUCCESS;
}
