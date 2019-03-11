#include <stdio.h>
#include "hemi/parallel_for.h"
#include "hemi/array.h"
#include "hemi/device_api.h"
#include "arrayHandle.h"
#include "boost/math/distributions/normal.hpp"
#include "stats-master/include/stats.hpp"

using namespace hemi;

/*quantili <- vector(length=0)
    for (i in 1:(ncx-1)) {
      quantili <- c(quantili, apply(abs(x[,-c(1:i), drop=FALSE]-as.vector(x[,i])), 2, use))
    }
    
  res <- quantile(quantili, probs)  
  if (size) {
    res <- list(quantile=res, stats=quantili, call=match.call())
  }
  return(res)
    */
__device__ ddouble * my_ddcpy(ddouble *dest, const ddouble *src){
  int i = 0;
  do {
    dest[i] = src[i];
	}while (src[i++].val != 0);
  return dest;
}

__device__ ddouble * my_ddcat(ddouble *dest, const ddouble *src){
  int i = 0;
  while (dest[i].val != 0) i++;
  my_ddcpy(dest+i, src);
  return dest;
}


void fillExample(ddouble * v, int n)
{ 
	double x;
	for (auto i : hemi::grid_stride_range(0, n)) 
	{
		x = (i+1)*2;
		v[i].val = x;
	}
}

double max(ddouble * vet, int n)
{
	double max = vet[0].val;
	for (auto i : hemi::grid_stride_range(0, n)) 
	{
			if(vet[i].val > max)
			{
				max = vet[i].val;
			}
	}
	return max;
}

ddouble * asVector(double ** mat, int nR, int nC)
{
	ddouble * res;
	int c = 0;
	for (auto i : hemi::grid_stride_range(0, nR)) 
	{
		for (auto j : hemi::grid_stride_range(0, nC)) 
		{
			res[c].val = mat[i][j];
			c++;
		}
	}
	return res;
}

double quantil(double **matrix, double probs, int nR, int nC)
{
	double res;

	hemi::Array<ddouble> transx(nC, true);
	hemi::Array<ddouble> supp(nR, true);

    hemi::Array<ddouble> quantili(20, true);
	hemi::Array<ddouble> arr2(20, true);
	hemi::Array<ddouble> result(40, true);

	fillExample(transx.writeOnlyHostPtr(), 20);
	fillExample(quantili.writeOnlyHostPtr(), 20);
	fillExample(arr2.writeOnlyHostPtr(), 40);
	
    ddouble a = quantili.writeOnlyHostPtr()[1].val;
	//trans(transx.writeOnlyPtr(), supp.writeOnlyPtr(), x, nR, nC);
	
    ddouble *q1 = new ddouble[20];
    ddouble *q2 = new ddouble[20];
    ddouble *ris = new ddouble[40];
    for(int i=0; i<20; i++)
    {
    	q1[i] = quantili.writeOnlyHostPtr()[i];
    	q2[i] = arr2.writeOnlyHostPtr()[i];
    	printf("Q2: %f\n", q2[i].val);
    }



	parallel_for(0, 4, [q1, q2, matrix] HEMI_LAMBDA (int f)  {
			
			my_ddcat(q1, q2);
		});

	for(int f=0; f<40; f++)
	{
		printf("Res1: %f\n", q1[f].val);
	}
	deviceSynchronize();
	boost::math::normal dist(0.0, probs);

	res = quantile(dist, 0.95);
	return res;
}

/*
__global__ void computeQuantilesKernel(double *matIn, int nRows, int nCols, int nQuantiles, double *outsideValues, double *quantilesAve, int param2)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    __shared__ double values[nCols];//big enough for 100 columns
    int keys[100];
    int nQuant[100];//big enough for 100 quantiles (percentiles)
    float thisQuantile[100];
    int quant;

    if (idx >= nRows) return;

    //read matIn from global memory
    for (int i = 0; i < nCols; i++)
    {
        values[i] = matIn[idx * nCols + i + param2 * nCols * nRows];
        keys[i] = i;
    }

    __syncthreads();
    //bubble Sort:
    for (int i = 0; i < nCols/2; i++)
	{
    int j = threadIdx.x;
    if (j % 2 == 0 && j<nCols-1)
        if (values[j+1] < values[j])
        {
            //swap(values[j+1], values[j]);
        	double appoggio;
        	appoggio = values[j+1];
        	values[j+1] = values[j];
        	values[j] = appoggio;
        }
    __syncthreads();
    if (j % 2 == 1 && j<nCols-1)
        if (values[j+1] < values[j])
        {
            //swap(values[j+1], values[j]);
            double appoggio;
        	appoggio = values[j+1];
        	values[j+1] = values[j];
        	values[j] = appoggio;
        }
    __syncthreads();
	}
    //end of bubble sort

    //reset nQuant and thisQuantile
    for (int iQuant = 0; iQuant < nQuantiles; iQuant++)
    {
        nQuant[iQuant] = 0;
        thisQuantile[iQuant] = 0;
    }

    //Compute sum of outsideValues for each quantile
    for (int i = 0; i < nCols; i++)
    {
        quant = (int)(((double)i + 0.5) / ((double)nCols / (double)nQuantiles));//quantile like Matlab
        nQuant[quant]++;
        thisQuantile[quant] += outsideValues[idx * nCols + keys[i]];
    }

    //Divide by the size of each quantile to get averages
    for (int iQuant = 0; iQuant < nQuantiles; iQuant++)
    {
        quantilesAve[idx + nRows * iQuant + param2 * nQuantiles * nRows] = thisQuantile[iQuant] / (float)nQuant[iQuant];
    }
}*/

extern "C"
int quantMain(double **matrix, int nR, int nC)
{

	boost::math::normal dist(0.0, 1.0);

	//ddouble * dist1 = asVector(matrix,nR,nC);
	// 95% of distribution is below q:
	double q = quantile(dist, 0.95);

	//double **app = stats::qbeta(matrix, 0.95, 0.9);

	double xxx = quantil(matrix, 0.03, 20, 20);

	printf("Quantile 1: %f\n", q); 

	deviceSynchronize();

	return 0;
}

