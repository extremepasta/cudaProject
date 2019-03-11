#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <algorithm>

#include "hemi/parallel_for.h"
#include "hemi/array.h"
#include "hemi/device_api.h"
#include "timer.h"

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


__host__ __device__ void swap(char *x, char *y)
{
    char temp;
    temp = *x;
    *x = *y;
    *y = temp;
}

void combn(int rangeInf, int rangeSup, int *com, int *pos, int *rInf)
{
	int limit = (rangeSup*rangeSup-1)/2;

	int OPT_SZ = limit * 2;

	printf("Dentro\n");

	printf("Dentro1 %d\n", limit);

	cudaMemset(com, 0, OPT_SZ);

	printf("Dentro4 %d\n", limit);

	hemi::parallel_for(0, limit, [=] HEMI_LAMBDA (int i) {
		if(i==0)
		{
			pos[0] = 1;
			rInf[0] = 1;
		}
		com[i] = rInf[0];
		com[i+1] = rInf[0]+pos[0];
		i++;
		if(rInf[0]+pos[0] == rangeSup)
		{
			pos[0] = rInf[0]+2;
			rInf[0]=rInf[0]+1;
		}
		else
				pos[0]=pos[0]+1; 
	});

}

double max(double * vet, int n)
{
	double max = vet[0];
	for (auto i : hemi::grid_stride_range(1, n-1)) 
	{
			if(vet[i] > max)
			{
				max = vet[i];
			}
	}
	return max;
}

double min(double * vet, int n)
{
	double min = vet[0];
	for (auto i : hemi::grid_stride_range(0, n-1)) 
	{
			if(vet[i] < min)
			{
				min = vet[i];
			}
	}
	return min;
}

double *asVecColumn(double ** mat, int col, int nR)
{
	double *res = new double[nR];
	for (auto i : hemi::grid_stride_range(0, nR-1)) 
	{
		res[i] = mat[i][col];
	}
	return res;
}

bool *banddepthforonecurve(int *x, double **xdata, double *ydata, double tau, int nR, int nC)
{
	//double **supporto;
	printf("Quaaa\n");
	int colUno = x[1];
	int colDue = x[2];
	double *sottrazione = new double[nR];
	double *envsup = new double[nR];
	double *envinf = new double[nR];
	bool *res = new bool[2];
	bool *inenvsup = new bool[nR];
	bool *inenvinf = new bool[nR];

	//SET ENVSUP & ENVINF
	for (auto i : hemi::grid_stride_range(0, nR-1)) 
	{
		if(xdata[i][colUno] >= xdata[i][colDue])
		{
			envsup[i] = xdata[i][colUno];
			envinf[i] = xdata[i][colDue];
		}
		else if(xdata[i][colUno] <= xdata[i][colDue])
		{
			envsup[i] = xdata[i][colDue];
			envinf[i] = xdata[i][colUno];
		}
		sottrazione[i] = envsup[i]-envinf[i];
	}

	int maxEnv = max(sottrazione, nR);

	//SET inenvsup & inenvinf
	for (auto j : hemi::grid_stride_range(0, nR-1)) 
	{
		if(ydata[j] <= envsup[j])
			inenvsup[j] = true;
		else
			inenvsup[j] = false;

		if(ydata[j] >= envinf[j])
			inenvinf[j] = true;
		else
			inenvinf[j] = false;
	}

	bool depth = true;
	//Set depth
	for (auto d : hemi::grid_stride_range(0, nR-1)) 
	{
		if(!inenvsup[d])
		{
			depth = false;
			d = nR;
		}
		else if(!inenvinf[d])
		{
			depth = false;
			d = nR;
		}
	}

	bool localdepth = depth & maxEnv <= tau;
	res[0] = depth;
	res[1] = localdepth;

	return res;
}


extern "C"
int secondMain(double **matrix, int nR, int nC, double tau)
{
	deviceSynchronize();

	StartTimer();

	for (auto i : hemi::grid_stride_range(0, nR)) 
	{
		for (auto j : hemi::grid_stride_range(0, nC)) 
		{
			printf(" %f", matrix[i][j]);
			if(j == nC-1)
				printf("\n");
		}
	}

	int rangeSup = 4;

	int limit = (rangeSup*rangeSup-1)/2;
	int OPT_SZ = limit * 2;
	
	bool *res = new bool[2];
    hemi::Array<int> com(OPT_SZ, true);
    hemi::Array<int> pos(1, true);
	hemi::Array<int> rInf(1, true);

	combn(1, rangeSup, com.writeOnlyPtr(), pos.writeOnlyPtr(), rInf.writeOnlyPtr());

	double ms = GetTimer();

	printf("\tSecondmain() time    : %f msec\n", ms);
	printf("Qui\n");

	//printf(" %i -", com.readOnlyPtr(hemi::host)[0]);
	//printf(" %i ", com.readOnlyPtr(hemi::host)[1]);
	printf("\n");

	double *ydata = asVecColumn(matrix, 1, nR);

	res = banddepthforonecurve(com.writeOnlyPtr(), matrix, ydata, tau, nR, nC);

	printf("Quantile 2: %f\n", tau); 

	deviceSynchronize();

	return 0;
}

