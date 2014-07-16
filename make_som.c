#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#define MAXLINE 4096

#define ELEM_SWAP(a,b) { register double t=(a);(a)=(b);(b)=t; }

double quick_median(double *arr, int n_elements) 
{
    int low, high ;
    int median;
    int middle, ll, hh;

    low = 0 ; high = n_elements-1 ; median = (low + high) / 2;
    for (;;) {
        if (high <= low) /* One element only */
            return arr[median] ;

        if (high == low + 1) {  /* Two elements only */
            if (arr[low] > arr[high])
                ELEM_SWAP(arr[low], arr[high]) ;
            return arr[median] ;
        }

    /* Find median of low, middle and high items; swap into position low */
    middle = (low + high) / 2;
    if (arr[middle] > arr[high])    ELEM_SWAP(arr[middle], arr[high]) ;
    if (arr[low] > arr[high])       ELEM_SWAP(arr[low], arr[high]) ;
    if (arr[middle] > arr[low])     ELEM_SWAP(arr[middle], arr[low]) ;

    /* Swap low item (now in position middle) into position (low+1) */
    ELEM_SWAP(arr[middle], arr[low+1]) ;

    /* Nibble from each end towards middle, swapping items when stuck */
    ll = low + 1;
    hh = high;
    for (;;) {
        do ll++; while (arr[low] > arr[ll]) ;
        do hh--; while (arr[hh]  > arr[low]) ;

        if (hh < ll)
        break;

        ELEM_SWAP(arr[ll], arr[hh]) ;
    }

    /* Swap middle item (in position low) back into correct position */
    ELEM_SWAP(arr[low], arr[hh]) ;

    /* Re-set active partition */
    if (hh <= median)
        low = ll;
        if (hh >= median)
        high = hh - 1;
    }

    return arr[median];
}

#undef ELEM_SWAP

double randn (double mu, double sigma){
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;
  if (call == 1){
    call = !call;
    return (mu + sigma*(double) X2);
  }
  do {
    U1 = -1 + ((double) rand() / RAND_MAX) * 2;
    U2 = -1 + ((double) rand() / RAND_MAX) * 2;
    W = pow(U1, 2) + pow(U2, 2);
  }
  while (W >= 1 || W == 0);
  mult = sqrt((-2*log(W))/W);
  X1 = U1*mult;
  X2 = U2*mult;

  call = !call;
  return (mu + sigma*(double) X1);
}

double learningFunction (double iter, double Niterations){
  //Return the current value of learning rate function a(t)
  return pow(0.01, (iter/Niterations));
}

double sigmaFunction (double iter, double Niterations, int gridside){
  //Return the value of sigma for the neighborhood function at this time
  double result;
  result = (double) gridside * pow((1./(double)gridside),(iter/Niterations));
  return result;
}

double neighborhoodFunction (double x1, double y1, double x2, double y2, double iter, double Niterations, double sigma){
  //Return neighborhood learning value at this iteration
  double sqdist;
  sqdist = pow((x1 - x2),2.) + pow((y1 - y2),2.);
  return exp(-sqdist/(sigma*sigma));
}

double calcMean (double *array, int n_elements){
  double total=0;
  int i=0;
  for(i=0;i<n_elements;i++){
    total += array[i];
  }
  return total/(double)n_elements;
}

double calcStdDev (double *array, int n_elements){
  double mean, total=0;
  int i = 0;
  mean = calcMean(array, n_elements);
  for(i=0;i<n_elements;i++){
    total += (array[i]-mean)*(array[i]-mean);
} 
return sqrt(total/(double)n_elements);
}

double calcSeparation(double *vec, double *cell, int Ndimensions){
  int i=0;
  double dist = 0.;
  for(i=0;i<Ndimensions;i++){
    dist += (vec[i] - cell[i])*(vec[i] - cell[i]);
  }
  return dist;  
}

//---------------------------------------------------------------------------------
int main(int argc, char *argv[]){
  // SOM implementation, on a rectangular grid.

  FILE *f_data, *f_test;
  double *fwhm, *ra, *dec, *i_auto, *redshift, *ebv, *Nf, *u, *g, *r, *iband, *z, *yg, *Y, *J, *H;
  int Ninput, *sed, *rlaw;
  long *ID;
  int gridside = 100; //size of rectangular grid on a side
  int Ncells = gridside*gridside;
  int Ndimensions = 8;
  int i = 0, j = 0;
  double iter = 0, Niterations = 100; 

  char line[MAXLINE], cmd[MAXLINE];

  if(argc != 4) {
    fprintf(stderr,"Use: %s nThreads input_data output \n", argv[0]);
    return -1;
  }

  f_data = fopen(argv[2],"r");
  if(f_data == NULL){
    fprintf(stderr,"Cannot open input file %s !\n",argv[2]);
    exit(1);
  }else{

    fprintf(stderr,"Reading data catalog %s.\n",argv[2]);
    fclose(f_data);
    sprintf(cmd,"wc -l < %s\n",argv[2]);
    f_data = popen(cmd,"r");
    fgets(line, MAXLINE, f_data); 
    sscanf(line, "%d", &Ninput);
        
    //Read in data
    f_data = fopen(argv[2],"r");

    ID = malloc(sizeof(long)*(Ninput));
    sed = malloc(sizeof(int)*(Ninput));
    rlaw = malloc(sizeof(int)*(Ninput));
    fwhm = malloc(sizeof(double)*Ninput);
    ra = malloc(sizeof(double)*(Ninput));
    dec = malloc(sizeof(double)*(Ninput));
    i_auto = malloc(sizeof(double)*(Ninput));
    redshift = malloc(sizeof(double)*(Ninput));
    ebv = malloc(sizeof(double)*(Ninput));
    Nf = malloc(sizeof(double)*(Ninput));
    u = malloc(sizeof(double)*(Ninput));
    g = malloc(sizeof(double)*(Ninput));
    r = malloc(sizeof(double)*(Ninput));
    iband = malloc(sizeof(double)*(Ninput));
    z = malloc(sizeof(double)*(Ninput));
    yg = malloc(sizeof(double)*(Ninput));
    Y = malloc(sizeof(double)*(Ninput));
    J = malloc(sizeof(double)*(Ninput));
    H = malloc(sizeof(double)*(Ninput));
   
    for(i=0;i<Ninput;i++) {
      fscanf(f_data,"%ld ",&ID[i]);
      fscanf(f_data,"%lf ",&fwhm[i]);
      fscanf(f_data,"%lf ",&ra[i]);
      fscanf(f_data,"%lf ",&dec[i]);
      fscanf(f_data,"%lf ",&i_auto[i]);
      fscanf(f_data,"%lf ",&redshift[i]);
      fscanf(f_data,"%d ",&sed[i]);
      fscanf(f_data,"%lf ",&ebv[i]);
      fscanf(f_data,"%d ",&rlaw[i]);
      fscanf(f_data,"%lf ",&Nf[i]);
      fscanf(f_data,"%lf ",&u[i]);
      fscanf(f_data,"%lf ",&g[i]);
      fscanf(f_data,"%lf ",&r[i]);
      fscanf(f_data,"%lf ",&iband[i]);
      fscanf(f_data,"%lf ",&z[i]);
      fscanf(f_data,"%lf ",&yg[i]);
      fscanf(f_data,"%lf ",&Y[i]);
      fscanf(f_data,"%lf ",&J[i]);
      fscanf(f_data,"%lf ",&H[i]);
    }
    
    fclose(f_data);

  } //end else

  //---------------------------------------------------------------------------------
  // Define colors, normalize them, and define the errors in the normalized values.

  int ngood=0;
  for(i=0;i<Ninput;i++){
    if(u[i] > 0.35 && u[i] < 1e4 && g[i] > 0.35 && g[i] < 1e4 && r[i] > 0.35 && r[i] < 1e4 && iband[i] > 0.35 && iband[i] < 1e4 && z[i] > 0.35 && z[i] < 1e4 && yg[i] > 0.35 && yg[i] < 1e4 && Y[i] > 0.35 && Y[i] < 1e4 && J[i] > 0.35 && J[i] < 1e4 && H[i] > 0.35 && H[i] < 1e4){
      ngood++;
    }
  }

  double *ug, *gr, *ri, *iz, *zyg, *ygY, *YJ, *JH;

  ug = malloc(sizeof(double)*ngood);
  gr = malloc(sizeof(double)*ngood);
  ri = malloc(sizeof(double)*ngood);
  iz = malloc(sizeof(double)*ngood);
  zyg = malloc(sizeof(double)*ngood);
  ygY = malloc(sizeof(double)*ngood);
  YJ = malloc(sizeof(double)*ngood);
  JH = malloc(sizeof(double)*ngood);
  
  int idx = 0;
  for(i=0;i<Ninput;i++){
    if(u[i] > 0.35 && u[i] < 1e4 && g[i] > 0.35 && g[i] < 1e4 && r[i] > 0.35 && r[i] < 1e4 && iband[i] > 0.35 && iband[i] < 1e4 && z[i] > 0.35 && z[i] < 1e4 && yg[i] > 0.35 && yg[i] < 1e4 && Y[i] > 0.35 && Y[i] < 1e4 && J[i] > 0.35 && J[i] < 1e4 && H[i] > 0.35 && H[i] < 1e4){
      ug[idx] = -2.5*(log10(u[i]) - log10(g[i]));   
      gr[idx] = -2.5*(log10(g[i]) - log10(r[i]));
      ri[idx] = -2.5*(log10(r[i]) - log10(iband[i]));
      iz[idx] = -2.5*(log10(iband[i]) - log10(z[i]));
      zyg[idx] = -2.5*(log10(z[i]) - log10(yg[i]));
      ygY[idx] = -2.5*(log10(yg[i]) - log10(Y[i]));
      YJ[idx] = -2.5*(log10(Y[i]) - log10(J[i]));
      JH[idx] = -2.5*(log10(J[i]) - log10(H[i]));
      idx++;
    }
  }    

  double med_ug, med_gr, med_ri, med_iz, med_zyg, med_ygY, med_YJ, med_JH;
  double stddev_ug, stddev_gr, stddev_ri, stddev_iz, stddev_zyg, stddev_ygY, stddev_YJ, stddev_JH;
  double ug_scaled[ngood], gr_scaled[ngood], ri_scaled[ngood], iz_scaled[ngood], zyg_scaled[ngood], ygY_scaled[ngood], YJ_scaled[ngood], JH_scaled[ngood];
  
  double temp[ngood];
  for(i=0;i<ngood;i++){
    temp[i] = ug[i];
  }
  med_ug = quick_median(temp,ngood);
  for(i=0;i<ngood;i++){
    temp[i] = gr[i];
  }
  med_gr = quick_median(temp,ngood);
  for(i=0;i<ngood;i++){
    temp[i] = ri[i];
  }
  med_ri = quick_median(temp,ngood);
  for(i=0;i<ngood;i++){
    temp[i] = iz[i];
  }
  med_iz = quick_median(temp,ngood);
  for(i=0;i<ngood;i++){
    temp[i] = zyg[i];
  }
  med_zyg = quick_median(temp,ngood);
  for(i=0;i<ngood;i++){
    temp[i] = ygY[i];
  }
  med_ygY = quick_median(temp,ngood);
  for(i=0;i<ngood;i++){
    temp[i] = YJ[i];
  }
  med_YJ = quick_median(temp,ngood);
  for(i=0;i<ngood;i++){
    temp[i] = JH[i];
  }
  med_JH = quick_median(temp,ngood);

  //Get data standard deviation values
  stddev_ug = calcStdDev(ug,ngood);
  stddev_gr = calcStdDev(gr,ngood);
  stddev_ri = calcStdDev(ri,ngood);
  stddev_iz = calcStdDev(iz,ngood);
  stddev_zyg = calcStdDev(zyg,ngood);
  stddev_ygY = calcStdDev(ygY,ngood);  
  stddev_YJ = calcStdDev(YJ,ngood);
  stddev_JH = calcStdDev(JH,ngood);
  
  // Scale the data using the median and standard deviation
  for(i=0;i<ngood;i++){
    ug_scaled[i] = (ug[i]-med_ug)/stddev_ug;
    gr_scaled[i] = (gr[i]-med_ug)/stddev_gr;
    ri_scaled[i] = (ri[i]-med_ug)/stddev_ri;
    iz_scaled[i] = (iz[i]-med_ug)/stddev_iz;
    zyg_scaled[i] = (zyg[i]-med_ug)/stddev_zyg;
    ygY_scaled[i] = (ygY[i]-med_ug)/stddev_ygY;
    YJ_scaled[i] = (YJ[i]-med_ug)/stddev_YJ;
    JH_scaled[i] = (JH[i]-med_ug)/stddev_JH;
  }

  //---------------------------------------------------------------------------------
  //Define and initialize the map
  struct SOM {
    double x, y; 
    double cellweights[Ndimensions]; 
  };

  struct SOM mySOM[Ncells]; 

  srand (time(NULL));
  
  //Initialize cell weights with normally distributed random numbers and assign cell positions
  for(i=0;i<Ncells;i++){
    mySOM[i].x = i % gridside;
    mySOM[i].y = floor((double) i / (double) gridside);
    for(j=0;j<Ndimensions;j++){			
      mySOM[i].cellweights[j] = randn(0,1);
    }
  }
 
  //---------------------------------------------------------------------------------
  // Train the map

  int randomIdx, idx_bmu;
  double sep, minsep, a_t, H_t, sigma_t, sigma_t_squared, trainingVector[Ndimensions];

  Niterations = 100000;

  for(iter=0;iter<Niterations;iter++){
    // Select at random from the training data
    // May want to change this to shuffling etc.
    randomIdx = rand() % ngood;
    trainingVector[0] = ug_scaled[randomIdx];
    trainingVector[1] = gr_scaled[randomIdx];
    trainingVector[2] = ri_scaled[randomIdx];
    trainingVector[3] = iz_scaled[randomIdx];
    trainingVector[4] = zyg_scaled[randomIdx];
    trainingVector[5] = ygY_scaled[randomIdx];
    trainingVector[6] = YJ_scaled[randomIdx];
    trainingVector[7] = JH_scaled[randomIdx];

    // Find the best-matching unit. 
    idx_bmu = 0;
    minsep = calcSeparation(trainingVector,mySOM[0].cellweights,Ndimensions);
    for(j=1;j<Ncells;j++){
      sep = calcSeparation(trainingVector,mySOM[j].cellweights,Ndimensions); //return squared Euclidean distance
      if (sep < minsep){
	minsep = sep;
	idx_bmu = j;
      }      
    }
    //fprintf(stderr,"This untrained minsep %lf\n",minsep);
    
    // We have the BMU, now update the map

    // Calculate a(t)
    a_t = learningFunction(iter,Niterations);
    sigma_t = sigmaFunction(iter,Niterations,gridside);
    sigma_t_squared = sigma_t * sigma_t;
    

    for(j=0;j<Ncells;j++){

      if( ((mySOM[j].x-mySOM[idx_bmu].x)*(mySOM[j].x-mySOM[idx_bmu].x) + (mySOM[j].y-mySOM[idx_bmu].y)*(mySOM[j].y-mySOM[idx_bmu].y)) <= 9*sigma_t_squared ){

	// Calculate the neighborhood function of this cell with the BMU
	H_t = neighborhoodFunction(mySOM[j].x, mySOM[j].y, mySOM[idx_bmu].x, mySOM[idx_bmu].y, iter, Niterations, sigma_t);
	
	// Update the cell's weights
	for(i=0;i<Ndimensions;i++){
	  mySOM[j].cellweights[i] = mySOM[j].cellweights[i] + a_t * H_t * (trainingVector[i] - mySOM[j].cellweights[i]);
	}
      } // end if cell is worth updating
    } // end loop over cells            
  } //end loop over iterations


  //---------------------------------------------------------------------------------
  // Calculate the quantization error in the map

  double mean_err, median_err, quantization_errors[ngood];

  struct SOM finalSOM[Ncells]; 
 
  /*
  for(i=0;i<Ncells;i++){
    finalSOM[i].x = i % gridside;
    finalSOM[i].y = floor((double) i / (double) gridside);		
    finalSOM[i].cellweights[0] = mySOM[i].cellweights[0]*stddev_ug + med_ug;
    finalSOM[i].cellweights[1] = mySOM[i].cellweights[1]*stddev_gr + med_gr;
    finalSOM[i].cellweights[2] = mySOM[i].cellweights[2]*stddev_ri + med_ri;
    finalSOM[i].cellweights[3] = mySOM[i].cellweights[3]*stddev_iz + med_iz;
    finalSOM[i].cellweights[4] = mySOM[i].cellweights[4]*stddev_zyg + med_zyg;
    finalSOM[i].cellweights[5] = mySOM[i].cellweights[5]*stddev_ygY + med_ygY;
    finalSOM[i].cellweights[6] = mySOM[i].cellweights[6]*stddev_YJ + med_YJ;
    finalSOM[i].cellweights[7] = mySOM[i].cellweights[7]*stddev_JH + med_JH;
  }
  
  for(i=0;i<ngood;i++){
    // Compute the quantization error of the map in the original (unscaled) data space
    trainingVector[0] = ug[i];
    trainingVector[1] = gr[i];
    trainingVector[2] = ri[i];
    trainingVector[3] = iz[i];
    trainingVector[4] = zyg[i];
    trainingVector[5] = ygY[i];
    trainingVector[6] = YJ[i];
    trainingVector[7] = JH[i];

    // Find the best-matching unit. 
    idx_bmu = 0;
    minsep = calcSeparation(trainingVector,finalSOM[0].cellweights,Ndimensions);
    for(j=1;j<Ncells;j++){
      sep = calcSeparation(trainingVector,finalSOM[j].cellweights,Ndimensions); //return squared Euclidean distance
      if (sep < minsep){
	minsep = sep;
	idx_bmu = j;
      }      
    }
    quantization_errors[i] = sqrt(minsep);
  }
  */
  
  for(i=0;i<ngood;i++){
    // Compute the quantization error of the map in the original (unscaled) data space
    trainingVector[0] = ug_scaled[i];
    trainingVector[1] = gr_scaled[i];
    trainingVector[2] = ri_scaled[i];
    trainingVector[3] = iz_scaled[i];
    trainingVector[4] = zyg_scaled[i];
    trainingVector[5] = ygY_scaled[i];
    trainingVector[6] = YJ_scaled[i];
    trainingVector[7] = JH_scaled[i];

    // Find the best-matching unit. 
    idx_bmu = 0;
    minsep = calcSeparation(trainingVector,mySOM[0].cellweights,Ndimensions);
    for(j=1;j<Ncells;j++){
      sep = calcSeparation(trainingVector,mySOM[j].cellweights,Ndimensions); //return squared Euclidean distance
      if (sep < minsep){
	minsep = sep;
	idx_bmu = j;
      }      
    }
    quantization_errors[i] = sqrt(minsep);
  }

  mean_err = calcMean(quantization_errors, ngood);
  median_err = quick_median(quantization_errors, ngood);
  fprintf(stderr,"The mean quantization error of this map is %lf\n",mean_err);
  fprintf(stderr,"The median quantization error of this map is %lf\n",median_err);
   
} // End SOM program


