#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <time.h>
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
  // Return the current value of learning rate function a(t)
  return pow(0.01, (iter/Niterations));
}

double sigmaFunction (double iter, double Niterations, double sigma_start){
  // Return the value of sigma for the neighborhood function at this time
  double result;
  result = sigma_start * pow((1./sigma_start),(iter/Niterations));
  return result;
}

double neighborhoodFunction (double x1, double y1, double x2, double y2, double iter, double Niterations, double sigma){
  // Return neighborhood learning value at this iteration
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
  // Return squared Euclidean distance between two vectors
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

  FILE *f_data, *f_test, *f_out;
  double *fwhm, *ra, *dec, *i_auto, *redshift, *ebv, *Nf, *u, *g, *r, *iband, *z, *yg, *Y, *J, *H;
  int Ninput, *sed, *rlaw;
  long *ID;
  int xgridside = 10; //size of x dimension of map
  int ygridside = 20; 
  int Ncells = xgridside*ygridside;
  double sigma_start = fminf((double)xgridside, (double)ygridside);
  fprintf(stderr,"sigma_start is %lf\n",sigma_start);
  int Ndimensions = 8;
  int i = 0, j = 0, k = 0;
  double iter = 0, Niterations = 100, elapsed; 
  clock_t start, end;

  char line[MAXLINE], cmd[MAXLINE];

  start = clock();

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

  double *ug, *gr, *ri, *iz, *zyg, *ygY, *YJ, *JH, *redshift_good;

  ug = malloc(sizeof(double)*ngood);
  gr = malloc(sizeof(double)*ngood);
  ri = malloc(sizeof(double)*ngood);
  iz = malloc(sizeof(double)*ngood);
  zyg = malloc(sizeof(double)*ngood);
  ygY = malloc(sizeof(double)*ngood);
  YJ = malloc(sizeof(double)*ngood);
  JH = malloc(sizeof(double)*ngood);
  redshift_good = malloc(sizeof(double)*ngood);
  
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
      redshift_good[idx] = redshift[i];
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
    gr_scaled[i] = (gr[i]-med_gr)/stddev_gr;
    ri_scaled[i] = (ri[i]-med_ri)/stddev_ri;
    iz_scaled[i] = (iz[i]-med_iz)/stddev_iz;
    zyg_scaled[i] = (zyg[i]-med_zyg)/stddev_zyg;
    ygY_scaled[i] = (ygY[i]-med_ygY)/stddev_ygY;
    YJ_scaled[i] = (YJ[i]-med_YJ)/stddev_YJ;
    JH_scaled[i] = (JH[i]-med_JH)/stddev_JH;
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
    mySOM[i].x = i % xgridside;
    mySOM[i].y = floor((double) i / (double) xgridside);
    for(j=0;j<Ndimensions;j++){			
      mySOM[i].cellweights[j] = randn(0,1);
    }
  }
 
  //---------------------------------------------------------------------------------
  // Train the map

  int randomIdx, idx_bmu;
  double sep, minsep, a_t, H_t, sigma_t, sigma_t_squared, trainingVector[Ndimensions];

  Niterations = 10000;
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
    
    // We found the best matching unit, now update the map

    // Calculate a(t)
    a_t = learningFunction(iter,Niterations);
    sigma_t = sigmaFunction(iter,Niterations,sigma_start);
    sigma_t_squared = sigma_t * sigma_t;
 
    // Loop through cells for update, only look at ones reasonably close to BMU

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
  // Define final map and compute its properties

  double mean_err, median_err, quantization_errors[ngood];

  struct finalSOM {
    int cell_index;
    double x, y; 
    double cellweights[Ndimensions];
    double cell_occupation; // total number of objects in training sample associating with this cell
    double U_value; // measure of separation from neighbor cells
    double cell_dispersion; // how similar are the objects that associate with this cell
    double cell_redshift; // average redshift of galaxies in cell, for use on training for now
  };
  struct finalSOM finalSOM[Ncells]; 
  
  for(i=0;i<Ncells;i++){
    finalSOM[i].cell_index = i;
    finalSOM[i].x = i % xgridside;
    finalSOM[i].y = floor((double) i / (double) xgridside);		
    finalSOM[i].cellweights[0] = mySOM[i].cellweights[0]*stddev_ug + med_ug;
    finalSOM[i].cellweights[1] = mySOM[i].cellweights[1]*stddev_gr + med_gr;
    finalSOM[i].cellweights[2] = mySOM[i].cellweights[2]*stddev_ri + med_ri;
    finalSOM[i].cellweights[3] = mySOM[i].cellweights[3]*stddev_iz + med_iz;
    finalSOM[i].cellweights[4] = mySOM[i].cellweights[4]*stddev_zyg + med_zyg;
    finalSOM[i].cellweights[5] = mySOM[i].cellweights[5]*stddev_ygY + med_ygY;
    finalSOM[i].cellweights[6] = mySOM[i].cellweights[6]*stddev_YJ + med_YJ;
    finalSOM[i].cellweights[7] = mySOM[i].cellweights[7]*stddev_JH + med_JH;
  }
  
  // Calculate the occupation of each cell and average redshift of galaxies in each cell
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
    // Update the occupation of the winning cell
    finalSOM[idx_bmu].cell_occupation += 1;
    finalSOM[idx_bmu].cell_redshift += redshift_good[i];
    finalSOM[idx_bmu].cell_dispersion += minsep;
  } 

  for(i=0;i<Ncells;i++){
    if(finalSOM[i].cell_occupation != 0.){
      finalSOM[i].cell_redshift /= finalSOM[i].cell_occupation;
      finalSOM[i].cell_dispersion /= finalSOM[i].cell_occupation;
      finalSOM[i].cell_dispersion = sqrt(finalSOM[i].cell_dispersion);
    }
  }

  // Compute cell U-matrix values, using largest separation with neighboring cell (Kraaijveld)
  double minx, maxx, miny, maxy, max_sep;
  int tempidx;
  
  for(i=0;i<Ncells;i++){
    max_sep = 0;
    minx = finalSOM[i].x - 1;
    if(minx < 0){minx = 0;}
    maxx = finalSOM[i].x + 1;
    if(maxx == xgridside){maxx = xgridside-1;}
    miny = finalSOM[i].y - 1;
    if(miny < 0){miny = 0;}
    maxy = finalSOM[i].y + 1;   
    if(maxy == ygridside){maxy = ygridside-1;}

    //fprintf(stderr,"This minx is %lf\n", minx);
    //fprintf(stderr,"This maxx is %lf\n", maxx);
    //fprintf(stderr,"This miny is %lf\n", miny);
    //fprintf(stderr,"This maxy is %lf\n", maxy);

    sep = 0;
    max_sep = 0;
    for(j=minx;j<=maxx;j++){
      for(k=miny;k<=maxy;k++){
	tempidx = k*xgridside + j;
	//fprintf(stderr,"This j is %d\n", j);
	//fprintf(stderr,"This k is %d\n", k);
	//fprintf(stderr,"This tempidx is %d\n", tempidx);	
	//fprintf(stderr,"This finalSOM[tempidx].x is %lf\n", finalSOM[tempidx].x);	
	//fprintf(stderr,"This finalSOM[tempidx].y is %lf\n", finalSOM[tempidx].y);
	sep = calcSeparation(finalSOM[tempidx].cellweights, finalSOM[i].cellweights, Ndimensions);
	//fprintf(stderr, "This separation is %lf\n", sqrt(sep));
	if(sep > max_sep){max_sep = sep;}
      }
    } // end loop over neighbors
    
    finalSOM[i].U_value = sqrt(max_sep);
  } // end loop over map cells 

  // Error checking, print average U value
  sep = 0;
  for(i=0;i<Ncells;i++){
    sep += finalSOM[i].U_value;
  }
  fprintf(stderr,"Average U value is %lf\n", sep/Ncells);
  

  //---------------------------------------------------------------------------------
  // Calculate the quantization error in the map

  //double total_error=0;
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
    
    /*
    total_error=0;
    //Let's check how well things match
    fprintf(stderr,"New object\n");
    for(j=0;j<Ndimensions;j++){
      fprintf(stderr, "Training vector[j] = %lf, Cell weight is %lf\n",trainingVector[j],finalSOM[idx_bmu].cellweights[j]);
      total_error += (trainingVector[j]-finalSOM[idx_bmu].cellweights[j])*(trainingVector[j]-finalSOM[idx_bmu].cellweights[j]);
    }
    fprintf(stderr,"Total error of %lf\n", sqrt(total_error));
    */
      
  }

  mean_err = calcMean(quantization_errors, ngood);
  median_err = quick_median(quantization_errors, ngood);
  fprintf(stderr,"The mean quantization error of this map is %lf\n",mean_err);
  fprintf(stderr,"The median quantization error of this map is %lf\n",median_err);

  
  //---------------------------------------------------------------------------------
  // Save the map to file
  // Format: index x  y  weight1 weight2 ... weightN U_value cell_dispersion cell_redshift

  f_out = fopen(argv[3],"w");
  if(f_out == NULL){
    fprintf(stderr,"Cannot open output file %s !\n",argv[3]);
    exit(1);
  }

  fprintf(f_out,"# Index X Y Weight1 Weight2 Weight3 Weight4 Weight5 Weight6 Weight7 Weight8 Cell_Occupation U_value Cell_Dispersion Cell_Average Redshift\n"); 
  for(i=0;i<Ncells;i++){
    fprintf(f_out,"%d %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",finalSOM[i].cell_index,finalSOM[i].x,finalSOM[i].y,finalSOM[i],finalSOM[i].cellweights[0],finalSOM[i].cellweights[1],finalSOM[i].cellweights[2],finalSOM[i].cellweights[3],finalSOM[i].cellweights[4],finalSOM[i].cellweights[5],finalSOM[i].cellweights[6],finalSOM[i].cellweights[7], finalSOM[i].cell_occupation, finalSOM[i].U_value, finalSOM[i].cell_dispersion, finalSOM[i].cell_redshift);
  }
  fclose(f_out);

  end = clock();
  elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
  fprintf(stderr,"It took %lf seconds.\n",elapsed); 

} // End SOM program


