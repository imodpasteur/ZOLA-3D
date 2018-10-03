

//#include "cublas_v2.h"
//#include "cublas.h"
 
#include <thrust/sort.h>





extern "C"
__global__ void vec_test1(int n,float *d_A,int size) {
	
    //float sum = thrust::reduce(thrust::seq, d_A, d_A + size);
	
	//thrust::sort_by_key(thrust::device,d_A, d_A + size, index);
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	if (id ==0)
    {
		thrust::sort(thrust::seq,d_A, d_A + size);
		printf("max side result = %f     %d\n", *(d_A+size-1),n);
	}

    

}




extern "C"
__global__ void vec_initIndex(int n, int *index)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        index[id] = id;
    }
}






extern "C"
__global__ void vec_computeLocalMaxima(int n, int *localMaxima,int sizeImage,float *input,int sizeFilt,int edgeSizeReject)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
		localMaxima[id]=0;
		int sizeSquare=sizeImage*sizeImage;
		int z=id/sizeSquare;
		int p=id%sizeSquare;
		int x=p/sizeImage;
		int y=p%sizeImage;
		int sizeReject=max(sizeFilt,edgeSizeReject);
		if ((x-sizeReject>=0)&&(y-sizeReject>=0)&&(x+sizeReject<sizeImage)&&(y+sizeReject<sizeImage)){//not edges
			int a,aa,neighbor;
			int ok=1;
			for (a=-sizeFilt;a<=sizeFilt;a++){
				for (a=-sizeFilt;a<=sizeFilt;a++){
					neighbor=z*sizeSquare+(x+a)*sizeImage+(y+aa);
					if (input[id]<input[neighbor]){
						ok=0;
						goto next;
					}
						
				}
			}
			next:
			if (ok==1){
				localMaxima[id]=1;
			}
		}
		
    }
}





extern "C"
__global__ void vec_eraseNonLocalMaxima(int n, float *input,int *localMaxima)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
		if (localMaxima[id]==0){
			input[id]=-1;
		}
		
    }
}





extern "C"
__global__ void vec_sortRows(int n,float *d_A,int * index,int size) {
	
    //float sum = thrust::reduce(thrust::seq, d_A, d_A + size);
	
	//thrust::sort_by_key(thrust::device,d_A, d_A + size, index);
	
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	if (id ==0)
    {
		thrust::stable_sort_by_key(thrust::seq,d_A, d_A + size, index,thrust::greater<float>());

	}
}




extern "C"
__global__ void vec_set (int n, double *result, double  value)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = value;
    }
}


//=== Vector arithmetic ======================================================

extern "C"
__global__ void vec_add (int n, double *result, double  *x, double  *y)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = x[id] + y[id];
    }
}


extern "C"
__global__ void vec_sub (int n, double *result, double  *x, double  *y)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = x[id] - y[id];
    }
}


extern "C"
__global__ void vec_mul (int n, double *result, double  *x, double  *y)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = x[id] * y[id];
    }
}



extern "C"
__global__ void vec_mul_fl (int n, float *result, float  *x, float  *y)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = x[id] * y[id];
    }
}



extern "C"
__global__ void vec_div (int n, double *result, double  *x, double  *y)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = x[id] / y[id];
    }
}

extern "C"
__global__ void vec_negate (int n, double *result, double  *x)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = -x[id];
    }
}




//=== Vector-and-scalar arithmetic ===========================================

extern "C"
__global__ void vec_addScalar (int n, double *result, double  *x, double  y)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = x[id] + y;
    }
}


extern "C"
__global__ void vec_subScalar (int n, double *result, double  *x, double  y)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = x[id] - y;
    }
}


extern "C"
__global__ void vec_mulScalar (int n, double *result, double  *x, double  y)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = x[id] * y;
    }
}


extern "C"
__global__ void vec_divScalar (int n, double *result, double  *x, double  y)
{

    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = x[id] / y;
    }
}




extern "C"
__global__ void vec_scalarAdd (int n, double *result, double  x, double  *y)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = x + y[id];
    }
}


extern "C"
__global__ void vec_scalarSub (int n, double *result, double  x, double  *y)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = x - y[id];
    }
}


extern "C"
__global__ void vec_scalarMul (int n, double *result, double  x, double  *y)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = x * y[id];
    }
}


extern "C"
__global__ void vec_scalarDiv (int n, double *result, double  x, double  *y)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = x / y[id];
    }
}









//=== Vector comparison ======================================================

extern "C"
__global__ void vec_lt (int n, double *result, double  *x, double  *y)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = (x[id] < y[id])?1.0:0.0;
    }
}


extern "C"
__global__ void vec_lte (int n, double *result, double  *x, double  *y)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = (x[id] <= y[id])?1.0:0.0;
    }
}


extern "C"
__global__ void vec_eq (int n, double *result, double  *x, double  *y)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = (x[id] == y[id])?1.0:0.0;
    }
}


extern "C"
__global__ void vec_gte (int n, double *result, double  *x, double  *y)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = (x[id] >= y[id])?1.0:0.0;
    }
}


extern "C"
__global__ void vec_gt (int n, double *result, double  *x, double  *y)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = (x[id] > y[id])?1.0:0.0;
    }
}



extern "C"
__global__ void vec_ne (int n, double *result, double  *x, double  *y)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = (x[id] != y[id])?1.0:0.0;
    }
}




//=== Vector-and-scalar comparison ===========================================

extern "C"
__global__ void vec_ltScalar (int n, double *result, double  *x, double  y)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = (x[id] < y)?1.0:0.0;
    }
}


extern "C"
__global__ void vec_lteScalar (int n, double *result, double  *x, double  y)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = (x[id] <= y)?1.0:0.0;
    }
}


extern "C"
__global__ void vec_eqScalar (int n, double *result, double  *x, double  y)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = (x[id] == y)?1.0:0.0;
    }
}


extern "C"
__global__ void vec_gteScalar (int n, double *result, double  *x, double  y)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = (x[id] >= y)?1.0:0.0;
    }
}


extern "C"
__global__ void vec_gtScalar (int n, double *result, double  *x, double  y)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = (x[id] > y)?1.0:0.0;
    }
}


extern "C"
__global__ void vec_neScalar (int n, double *result, double  *x, double  y)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = (x[id] != y)?1.0:0.0;
    }
}











//=== Vector math (one argument) =============================================


// Calculate the arc cosine of the input argument.
extern "C"
__global__ void vec_acos (int n, double *result, double  *x)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = acos(x[id]);
    }
}


// Calculate the nonnegative arc hyperbolic cosine of the input argument.
extern "C"
__global__ void vec_acosh (int n, double *result, double  *x)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = acosh(x[id]);
    }
}


// Calculate the arc sine of the input argument.
extern "C"
__global__ void vec_asin (int n, double *result, double  *x)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = asin(x[id]);
    }
}


// Calculate the arc hyperbolic sine of the input argument.
extern "C"
__global__ void vec_asinh (int n, double *result, double  *x)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = asinh(x[id]);
    }
}


// Calculate the arc tangent of the input argument.
extern "C"
__global__ void vec_atan (int n, double *result, double  *x)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = atan(x[id]);
    }
}


// Calculate the arc hyperbolic tangent of the input argument.
extern "C"
__global__ void vec_atanh (int n, double *result, double  *x)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = atanh(x[id]);
    }
}


// Calculate the cube root of the input argument.
extern "C"
__global__ void vec_cbrt (int n, double *result, double  *x)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = cbrt(x[id]);
    }
}


// Calculate ceiling of the input argument.
extern "C"
__global__ void vec_ceil (int n, double *result, double  *x)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = ceil(x[id]);
    }
}


// Calculate the cosine of the input argument.
extern "C"
__global__ void vec_cos (int n, double *result, double  *x)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = cos(x[id]);
    }
}


// Calculate the hyperbolic cosine of the input argument.
extern "C"
__global__ void vec_cosh (int n, double *result, double  *x)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = cosh(x[id]);
    }
}


// Calculate the cosine of the input argument × p .
extern "C"
__global__ void vec_cospi (int n, double *result, double  *x)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = cospi(x[id]);
    }
}


// Calculate the complementary error function of the input argument.
extern "C"
__global__ void vec_erfc (int n, double *result, double  *x)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = erfc(x[id]);
    }
}


// Calculate the inverse complementary error function of the input argument.
extern "C"
__global__ void vec_erfcinv (int n, double *result, double  *y)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = erfcinv(y[id]);
    }
}


// Calculate the scaled complementary error function of the input argument.
extern "C"
__global__ void vec_erfcx (int n, double *result, double  *x)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = erfcx(x[id]);
    }
}


// Calculate the error function of the input argument.
extern "C"
__global__ void vec_erf (int n, double *result, double  *x)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = erf(x[id]);
    }
}


// Calculate the inverse error function of the input argument.
extern "C"
__global__ void vec_erfinv (int n, double *result, double  *y)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = erfinv(y[id]);
    }
}


// Calculate the base 10 exponential of the input argument.
extern "C"
__global__ void vec_exp10 (int n, double *result, double  *x)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = exp10(x[id]);
    }
}


// Calculate the base 2 exponential of the input argument.
extern "C"
__global__ void vec_exp2 (int n, double *result, double  *x)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = exp2(x[id]);
    }
}


// Calculate the base e exponential of the input argument.
extern "C"
__global__ void vec_exp (int n, double *result, double  *x)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = exp(x[id]);
    }
}


// Calculate the base e exponential of the input argument, minus 1.
extern "C"
__global__ void vec_expm1 (int n, double *result, double  *x)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = expm1(x[id]);
    }
}


// Calculate the absolute value of its argument.
extern "C"
__global__ void vec_fabs (int n, double *result, double  *x)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = fabs(x[id]);
    }
}


// Calculate the largest integer less than or equal to x.
extern "C"
__global__ void vec_floor (int n, double *result, double  *x)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = floor(x[id]);
    }
}


// Calculate the value of the Bessel function of the first kind of order 0 for the input argument.
extern "C"
__global__ void vec_j0 (int n, double *result, double  *x)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = j0(x[id]);
    }
}


// Calculate the value of the Bessel function of the first kind of order 1 for the input argument.
extern "C"
__global__ void vec_j1 (int n, double *result, double  *x)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = j1(x[id]);
    }
}


// Calculate the natural logarithm of the absolute value of the gamma function of the input argument.
extern "C"
__global__ void vec_lgamma (int n, double *result, double  *x)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = lgamma(x[id]);
    }
}


// Calculate the base 10 logarithm of the input argument.
extern "C"
__global__ void vec_log10 (int n, double *result, double  *x)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = log10(x[id]);
    }
}


// Calculate the value of l o g e ( 1 + x ) .
extern "C"
__global__ void vec_log1p (int n, double *result, double  *x)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = log1p(x[id]);
    }
}


// Calculate the base 2 logarithm of the input argument.
extern "C"
__global__ void vec_log2 (int n, double *result, double  *x)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = log2(x[id]);
    }
}


// Calculate the doubleing point representation of the exponent of the input argument.
extern "C"
__global__ void vec_logb (int n, double *result, double  *x)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = logb(x[id]);
    }
}


// Calculate the natural logarithm of the input argument.
extern "C"
__global__ void vec_log (int n, double *result, double  *x)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = log(x[id]);
    }
}


// Calculate the standard normal cumulative distribution function.
extern "C"
__global__ void vec_normcdf (int n, double *result, double  *y)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = normcdf(y[id]);
    }
}


// Calculate the inverse of the standard normal cumulative distribution function.
extern "C"
__global__ void vec_normcdfinv (int n, double *result, double  *y)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = normcdfinv(y[id]);
    }
}


// Calculate reciprocal cube root function.
extern "C"
__global__ void vec_rcbrt (int n, double *result, double  *x)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = rcbrt(x[id]);
    }
}


// Round input to nearest integer value in doubleing-point.
extern "C"
__global__ void vec_rint (int n, double *result, double  *x)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = rint(x[id]);
    }
}


// Round to nearest integer value in doubleing-point.
extern "C"
__global__ void vec_round (int n, double *result, double  *x)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = round(x[id]);
    }
}


// Calculate the reciprocal of the square root of the input argument.
extern "C"
__global__ void vec_rsqrt (int n, double *result, double  *x)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = rsqrt(x[id]);
    }
}


// Calculate the sine of the input argument.
extern "C"
__global__ void vec_sin (int n, double *result, double  *x)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = sin(x[id]);
    }
}


// Calculate the hyperbolic sine of the input argument.
extern "C"
__global__ void vec_sinh (int n, double *result, double  *x)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = sinh(x[id]);
    }
}


// Calculate the sine of the input argument × p .
extern "C"
__global__ void vec_sinpi (int n, double *result, double  *x)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = sinpi(x[id]);
    }
}


// Calculate the square root of the input argument.
extern "C"
__global__ void vec_sqrt (int n, double *result, double  *x)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = sqrt(x[id]);
    }
}


// Calculate the tangent of the input argument.
extern "C"
__global__ void vec_tan (int n, double *result, double  *x)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = tan(x[id]);
    }
}


// Calculate the hyperbolic tangent of the input argument.
extern "C"
__global__ void vec_tanh (int n, double *result, double  *x)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = tanh(x[id]);
    }
}


// Calculate the gamma function of the input argument.
extern "C"
__global__ void vec_tgamma (int n, double *result, double  *x)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = tgamma(x[id]);
    }
}


// Truncate input argument to the integral part.
extern "C"
__global__ void vec_trunc (int n, double *result, double  *x)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = trunc(x[id]);
    }
}


// Calculate the value of the Bessel function of the second kind of order 0 for the input argument.
extern "C"
__global__ void vec_y0 (int n, double *result, double  *x)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = y0(x[id]);
    }
}


// Calculate the value of the Bessel function of the second kind of order 1 for the input argument.
extern "C"
__global__ void vec_y1 (int n, double *result, double  *x)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = y1(x[id]);
    }
}











//=== Vector math (two arguments) ============================================





// Create value with given magnitude, copying sign of second value.
extern "C"
__global__ void vec_copysign (int n, double *result, double  *x, double  *y)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = copysign(x[id], y[id]);
    }
}

// Compute the positive difference between x and y.
extern "C"
__global__ void vec_fdim (int n, double *result, double  *x, double  *y)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = fdim(x[id], y[id]);
    }
}

// Divide two doubleing point values.
extern "C"
__global__ void vec_fdivide (int n, double *result, double  *x, double  *y)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = fdivide(x[id], y[id]);
    }
}

// Determine the maximum numeric value of the arguments.
extern "C"
__global__ void vec_fmax (int n, double *result, double  *x, double  *y)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = fmax(x[id], y[id]);
    }
}

// Determine the minimum numeric value of the arguments.
extern "C"
__global__ void vec_fmin (int n, double *result, double  *x, double  *y)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = fmin(x[id], y[id]);
    }
}

// Calculate the doubleing-point remainder of x / y.
extern "C"
__global__ void vec_fmod (int n, double *result, double  *x, double  *y)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = fmod(x[id], y[id]);
    }
}

// Calculate the square root of the sum of squares of two arguments.
extern "C"
__global__ void vec_hypot (int n, double *result, double  *x, double  *y)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = hypot(x[id], y[id]);
    }
}

// Return next representable single-precision doubleing-point value afer argument.
extern "C"
__global__ void vec_nextafter (int n, double *result, double  *x, double  *y)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = nextafter(x[id], y[id]);
    }
}

// Calculate the value of first argument to the power of second argument.
extern "C"
__global__ void vec_pow (int n, double *result, double  *x, double  *y)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = pow(x[id], y[id]);
    }
}

// Compute single-precision doubleing-point remainder.
extern "C"
__global__ void vec_remainder (int n, double *result, double  *x, double  *y)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = remainder(x[id], y[id]);
    }
}















extern "C"
__global__ void vec_testkernel (int n, double *result, double  *x, double  *y)
{
    
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	for (int j=0;j<100;j++){//stupid...just to test speed
	    result[id] = x[id] * y[id]+sqrt(pow(3.14159,id));;;
	}



}









extern "C"
__global__ void vec_computePSF_phase (int n, double *realOutput, double *imagOutput, double  *kx, double  *ky, double  *kz, double  *pupil, double  *phase,double dx, double dy, double dz)
{
    double x;
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	if (id < n)
    {
        x= kx[id]*dx + ky[id]*dy + kz[id]*dz + phase[id];
		realOutput[id]=pupil[id]*cos(x);
		imagOutput[id]=pupil[id]*sin(x);
		//realOutput[id]=phase[id];
		//imagOutput[id]=sin(phase[id]);
    }



}



extern "C"
__global__ void vec_computePSF_phaseN (int n, double  *kx, double  *ky, double  *kz, double  *pupil, double  *phase,double dx, double dy, double dz, int *sparseIndexEvenDisk, int *sparseIndexOddDisk, double *fft)
{
    double x;

	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	if (id < n)
    {
        x= kx[id]*dx + ky[id]*dy + kz[id]*dz + phase[id];
		
		fft[sparseIndexEvenDisk[id]]=pupil[id]*cos(x);
		fft[sparseIndexOddDisk[id]]=pupil[id]*sin(x);


		
    }



}





extern "C"
__global__ void vec_computePSF_phaseNwithOil (int n, double  *kx, double  *ky, double  *kz,double  *kz_is_imag, double  *kz_oil,double  *kz_oil_is_imag, double  *pupil, double  *phase,double dx, double dy, double dz, double dz_oil, int *sparseIndexEvenDisk, int *sparseIndexOddDisk, double *fft)
{
    double x,y,z;

	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	if (id < n)
    {
		z= kx[id]*dx + ky[id]*dy + phase[id];
        x= z + kz[id]*dz - kz_oil[id]*dz_oil;
		y= z + kz[id]*dz*kz_is_imag[id] - kz_oil[id]*dz_oil*kz_oil_is_imag[id];
		fft[sparseIndexEvenDisk[id]]=pupil[id]*cos(x);
		fft[sparseIndexOddDisk[id]]=pupil[id]*sin(y);


		
    }



}




extern "C"         
__global__ void vec_computePSF_phaseNMany (int n,int sizePart,int sizeTot, double  *kx, double  *ky, double  *kz, double  *pupil, double  *phase,double* position, int *sparseIndexEvenDisk, int *sparseIndexOddDisk, double *fft,int many)
{
    double x;
	
	int u,p;
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	if (id < n)
    {
			u=id/sizePart;
			p=id%sizePart;
		    x= kx[p]*position[u] + ky[p]*position[u+many] + kz[p]*position[u+2*many] + phase[p];
			fft[sparseIndexEvenDisk[p]+u*sizeTot]=pupil[p]*cos(x);
			fft[sparseIndexOddDisk[p]+u*sizeTot]=pupil[p]*sin(x);
			
			

    }



}




extern "C"         
__global__ void vec_computePSF_phaseNMany_f (int n,int sizePart,int sizeTot, float  *kx, float  *ky, float  *kz, float  *pupil, float  *phase,float* position, int *sparseIndexEvenDisk, int *sparseIndexOddDisk, float *fft,int many)
{
    float x;
	
	int u,p;
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	if (id < n)
    {
			u=id/sizePart;
			p=id%sizePart;
		    x= kx[p]*position[u] + ky[p]*position[u+many] + kz[p]*position[u+2*many] + phase[p];
			fft[sparseIndexEvenDisk[p]+u*sizeTot]=pupil[p]*cos(x);
			fft[sparseIndexOddDisk[p]+u*sizeTot]=pupil[p]*sin(x);
			
			

    }



}




extern "C"         
__global__ void vec_computePSF_phaseNManywithOil_f (int n,int sizePart,int sizeTot, float  *kx, float  *ky, float  *kz,float  *kz_is_imag, float  *kz_oil,float  *kz_oil_is_imag, float  *pupil, float  *phase,float* position, int *sparseIndexEvenDisk, int *sparseIndexOddDisk, float *fft,int many)
{
    float x,y,z;
	//float x;
	
	int u,p;
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	if (id < n)
    {
			u=id/sizePart;
			p=id%sizePart;

			//x= kx[p]*position[u] + ky[p]*position[u+many] + phase[p] + kz[p]*position[u+2*many]*kz_is_imag[p] - kz_oil[p]*position[u+3*many]*kz_oil_is_imag[p];
			//fft[sparseIndexEvenDisk[p]+u*sizeTot]=pupil[p]*cos(x);
			//fft[sparseIndexOddDisk[p]+u*sizeTot]=pupil[p]*sin(x);
			
			
			z= kx[p]*position[u] + ky[p]*position[u+many] + phase[p];
		    x= z + kz[p]*position[u+2*many] - kz_oil[p]*position[u+3*many];
			y= z + kz[p]*position[u+2*many]*kz_is_imag[p] - kz_oil[p]*position[u+3*many]*kz_oil_is_imag[p];
			fft[sparseIndexEvenDisk[p]+u*sizeTot]=pupil[p]*cos(x);
			fft[sparseIndexOddDisk[p]+u*sizeTot]=pupil[p]*sin(y);
			
			

    }



}



extern "C"         
__global__ void vec_thetest(int n,float  *kz_is_imag)
{
    
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	if (id < n)
    {
			
printf("kz -> %d  %f  \n",id,kz_is_imag[id]);
			
			

    }



}



//WARNING : device_sum size should be gridDim.x
__device__ void sumTheBlocks (double *device_sum)
{

	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	if (id < 1)//lose few time... i agree
    {
		for (int i=1;i<gridDim.x;i++){
			device_sum[0]+=device_sum[i];
		}
		
	}
}


__device__ int truc;
__device__ int barrier;
__device__ int barrier2;
__device__ void initSynchroBlocks(){
	if (threadIdx.x + blockIdx.x * blockDim.x==0){
		barrier=gridDim.x;
		barrier2=gridDim.x;
	}
}

__device__ void synchroBlocks(){
	__syncthreads();
	if (threadIdx.x==0){
		atomicSub( &barrier , 1 );
	}
	//atomicAdd( &truc , 1 );
	//if ( threadIdx.x == 0 )
        //while ( atomicCAS( &barrier , 0 , 0 ) != 0 );//does not work: infinite loop...

    __syncthreads();
}
__device__ void synchroBlocks2(){
	__syncthreads();
	if (threadIdx.x==0){
		atomicSub( &barrier2 , 1 );
	}
	//atomicAdd( &truc , 1 );
	//if ( threadIdx.x == 0 )
        //while ( atomicCAS( &barrier2 , 0 , 0 ) != 0 );//does not work: infinite loop...

    __syncthreads();
}








extern "C"
__device__ void vec_dense2Sparse (int n, double *device_input, int *device_sparse, double *device_output)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	if (id < n)
    {
		device_output[id]=device_input[device_sparse[id]];
    }
}


extern "C"
__device__ void vec_sparse2Dense (int n, double *device_input, int *device_sparse, double *device_output)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	if (id < n)
    {
		device_output[device_sparse[id]]=device_input[id];
    }
}



extern "C"
__global__ void vec_computePSF_signal (int n, double *result, double *real, double  *imag,double divide)
{
    double x,y;
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	if (id < n)
    {
        x=real[id]/divide;
		y=imag[id]/divide;
		result[id]=x*x+y*y;
    }



}




extern "C"
__global__ void vec_computePSF_signalN (int n, double *result, double divide, int *sparseIndexEvenShiftOutput, int *sparseIndexOddShiftOutput, double *fft)
{
    double x,y;
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	if (id < n)
    {
		x=fft[sparseIndexEvenShiftOutput[id]]/divide;
		y=fft[sparseIndexOddShiftOutput[id]]/divide;
		result[id]=x*x+y*y;
    }



}


extern "C"
__global__ void vec_computePSF_signalN2 (int n, double divide, int *sparseIndexEvenShiftOutput, int *sparseIndexOddShiftOutput, double *fft,int *sparseIndexEvenShiftOutputNext,int *sparseIndexOddShiftOutputNext, double *psffft)
{
    double x,y;
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	if (id < n)
    {
		x=fft[sparseIndexEvenShiftOutput[id]]/divide;
		y=fft[sparseIndexOddShiftOutput[id]]/divide;
		psffft[sparseIndexEvenShiftOutputNext[id]]=x*x+y*y;
		psffft[sparseIndexOddShiftOutputNext[id]]=0;

		//psffft[id]=x*x+y*y;
    }



}





extern "C"
__global__ void vec_computePSF_signalN2Many (int n,int sizePart,int sizeTot, double divide, int *sparseIndexEvenShiftOutput, int *sparseIndexOddShiftOutput, double *fft,int *sparseIndexEvenShiftOutputNext,int *sparseIndexOddShiftOutputNext, double *psffft)
{
	int u,p;
	
	
    double x,y;
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	if (id < n)
    {
		u=id/sizePart;
		p=id%sizePart;
		x=fft[sparseIndexEvenShiftOutput[p]+sizeTot*u]/divide;
		y=fft[sparseIndexOddShiftOutput[p]+sizeTot*u]/divide;
		psffft[sparseIndexEvenShiftOutputNext[p]+2*sizePart*u]=x*x+y*y;
		psffft[sparseIndexOddShiftOutputNext[p]+2*sizePart*u]=0;

		//psffft[id]=x*x+y*y;
    }



}






extern "C"
__global__ void vec_computePSF_signalN2Many_f (int n,int sizePart,int sizeTot, float divide, int *sparseIndexEvenShiftOutput, int *sparseIndexOddShiftOutput, float *fft,int *sparseIndexEvenShiftOutputNext,int *sparseIndexOddShiftOutputNext, float *psffft)
{
	int u,p;
	
	
    float x,y;
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	if (id < n)
    {
		u=id/sizePart;
		p=id%sizePart;
		x=fft[sparseIndexEvenShiftOutput[p]+sizeTot*u]/divide;
		y=fft[sparseIndexOddShiftOutput[p]+sizeTot*u]/divide;
		psffft[sparseIndexEvenShiftOutputNext[p]+2*sizePart*u]=x*x+y*y;
		psffft[sparseIndexOddShiftOutputNext[p]+2*sizePart*u]=0;

		//psffft[id]=x*x+y*y;
    }



}




extern "C"
__global__ void vec_computePSF_signalsqrt (int n, double *result, double *real, double  *imag,double divide)
{
    double x,y;
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	if (id < n)
    {
        x=real[id]/divide;
		y=imag[id]/divide;
		result[id]=sqrt(x*x+y*y);
    }



}





extern "C"
__global__ void vec_computePSF_signalNsqrt (int n, double *result, double *fft,double divide, int *sparseIndexEvenShift2D, int *sparseIndexOddShift2D)
{
    double x,y;
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	if (id < n)
    {
        x=fft[sparseIndexEvenShift2D[id]]/divide;
		y=fft[sparseIndexOddShift2D[id]]/divide;
		result[id]=sqrt(x*x+y*y);
    }



}




extern "C"
__global__ void vec_computePSF_signalNsqrtMany (int n, int sizeSubImage,double *result, double *fft,double divide, int *sparseIndexEvenShift2D, int *sparseIndexOddShift2D)
{
    double x,y;
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
 	int id2=id%sizeSubImage;
	int id3=id/sizeSubImage;
	if (id < n)
    {
        x=fft[sparseIndexEvenShift2D[id2]+id3*sizeSubImage*2]/divide;
		y=fft[sparseIndexOddShift2D[id2]+id3*sizeSubImage*2]/divide;
		result[id]=sqrt(x*x+y*y);
    }



}




extern "C"
__global__ void vec_computePSF_signalNsqrtMany_f (int n, int sizeSubImage,float *result, float *fft,float divide, int *sparseIndexEvenShift2D, int *sparseIndexOddShift2D)
{
    float x,y;
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
 	int id2=id%sizeSubImage;
	int id3=id/sizeSubImage;
	if (id < n)
    {
        x=fft[sparseIndexEvenShift2D[id2]+id3*sizeSubImage*2]/divide;
		y=fft[sparseIndexOddShift2D[id2]+id3*sizeSubImage*2]/divide;
		result[id]=sqrt(x*x+y*y);
    }



}




extern "C"
__global__ void vec_computePSF_signalNsqrtMany_fcrop (int n, int sizeSubImage, int sizeSubImageFull,float *result, float *fft,float divide, int *sparseIndexEvenShift2D, int *sparseIndexOddShift2D)
{
    float x,y;
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
 	int id2=id%sizeSubImage;
	int id3=id/sizeSubImage;
	int id4=id3*sizeSubImageFull;
	if (id < n)
    {
        x=fft[sparseIndexEvenShift2D[id2]+id4*2]/divide;
		y=fft[sparseIndexOddShift2D[id2]+id4*2]/divide;
		result[id]=sqrt(x*x+y*y);
    }



}

__device__ void divideBySum (int n, double *result, double *tmpsum)
{
    
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	synchroBlocks2();
	if (id < n)
    {	
		result[id]/=tmpsum[0];
	}

}







extern "C"
__global__ void vec_mulMany (int n, int sizeKernel, double *result, double  *x, double  *kernel)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	int id2=id%sizeKernel;
    if (id < n)
    {
        result[id] = x[id] * kernel[id2];
    }
}



extern "C"
__global__ void vec_divScalarMany (int n,  int sizeSubImage,double *result, double  *x, double  *div)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	int id2=id/sizeSubImage;
    if (id < n)
    {
		if (div[id2]!=0)
        	result[id] = x[id] / div[id2];
    }
}




extern "C"
__global__ void vec_mulMany_f (int n, int sizeKernel, float *result, float  *x, float  *kernel)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	int id2=id%sizeKernel;
    if (id < n)
    {
        result[id] = x[id] * kernel[id2];
    }
}



extern "C"
__global__ void vec_computeModelMany1 (int n, int sizeImage, double *result, double  *x, double  *amplitude,double background)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	int id2=id/sizeImage;
    if (id < n)
    {
        result[id] = x[id] * amplitude[id2] + background;
    }
}


extern "C"
__global__ void vec_computeModelMany2 (int n, int sizeImage, double *result, double  *x, double  *amplitude,double *background)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	int id2=id/sizeImage;
    if (id < n)
    {
        result[id] = x[id] * amplitude[id2] + background[id2];
    }
}



//here, background is 2D
extern "C"
__global__ void vec_computeModelMany3 (int n, int sizeImage, double *result, double  *x, double  *amplitude,double *background)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	int id2=id/sizeImage;
	int id3=id%sizeImage;
    if (id < n)
    {

        result[id] = x[id] * amplitude[id2] + background[id3];
		
    }
}



extern "C"
__global__ void vec_computeModelMany1_scmos (int n, int sizeImage, double *result, double  *x, double  *amplitude,double background,double  *scmos)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	int id2=id/sizeImage;
	int id3=id%sizeImage;
    if (id < n)
    {
        result[id] = x[id] * amplitude[id2] + background + scmos[id3];
    }
}


extern "C"
__global__ void vec_computeModelMany2_scmos (int n, int sizeImage, double *result, double  *x, double  *amplitude,double *background,double  *scmos)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	int id2=id/sizeImage;
	int id3=id%sizeImage;
    if (id < n)
    {
        result[id] = x[id] * amplitude[id2] + background[id2] + scmos[id3];
    }
}



//here, background is 2D
extern "C"
__global__ void vec_computeModelMany3_scmos (int n, int sizeImage, double *result, double  *x, double  *amplitude,double *background,double  *scmos)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	int id2=id/sizeImage;
	int id3=id%sizeImage;
    if (id < n)
    {

        result[id] = x[id] * amplitude[id2] + background[id3] + scmos[id3];
		
    }
}



extern "C"
__global__ void vec_divScalarMany_f (int n,  int sizeSubImage,float *result,double *resultDouble, float  *x, float  *div)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	int id2=id/sizeSubImage;
    if (id < n)
    {
		if (div[id2]!=0){
        	result[id] = x[id] / div[id2];
			resultDouble[id] =0;
        	resultDouble[id] =(double)(result[id]);
		}
    }
}








extern "C"
__global__ void vec_computePoissonLikelihood (int n, double *result, double *image, double  *model)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	if (id < n)
    {
		if (model[id]>0)
			result[id]=model[id]-image[id]*log(model[id]);
		else
			result[id]=10000000;
			
    }



}






//reshuffle: 
//exemple 4 PSF to merge in 2 model
//->>> PSF=1,2,3 merged with PSF=4,5,6
//->>> PSF=7,8,9 merged with PSF=10,11,12
//we need to reshuffle here to use then Dgemv for matrix operation
//1,2,3,4,5,6,7,8,9,10,11,12 -> 1,2,3,7,8,9,4,5,6,10,11,12
extern "C"
__global__ void vec_addPhotonsAndBackgroundManyReshuffle (int n, int sizeSubImage,int numberPSFperModel,double *output, double *input, double *photonAndBackground)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;

	if (id < n)
	{
		int idPSF=id/sizeSubImage;
		int idModel=idPSF/numberPSFperModel;
		int idoffset=id%sizeSubImage;
		int idposit=idPSF%numberPSFperModel;
		int idreshuffle=idModel*sizeSubImage +idposit*sizeSubImage*(n/(sizeSubImage*numberPSFperModel))       +idoffset;
		output[idreshuffle]=input[id]*photonAndBackground[idPSF*2]+photonAndBackground[idPSF*2+1];

	}

}






extern "C"
__global__ void vec_addPhotonsAndBackgroundManyReshuffle_scmos (int n, int sizeSubImage,int numberPSFperModel,double *output, double *input, double *photonAndBackground, double * scmos)
{


//print("to do as previous function");



	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	
	if (id < n)
    {
	int idPSF=id/sizeSubImage;
	int idModel=idPSF/numberPSFperModel;
	int idoffset=id%sizeSubImage;
	int idposit=idPSF%numberPSFperModel;
	int idreshuffle=idModel*sizeSubImage +idposit*sizeSubImage*(n/(sizeSubImage*numberPSFperModel))       +idoffset;
	output[idreshuffle]=input[id]*photonAndBackground[idPSF*2]+photonAndBackground[idPSF*2+1]+scmos[id];

    }

}








extern "C"
__global__ void vec_addPhotonsAndBackgroundMany (int n, int sizeSubImage,double *output, double *input, double *photonAndBackground)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	int id2=id/sizeSubImage;
	if (id < n)
    {
		output[id]=input[id]*photonAndBackground[id2*2]+photonAndBackground[id2*2+1];

    }

}






extern "C"
__global__ void vec_addPhotonsAndBackgroundMany_scmos (int n, int sizeSubImage,double *output, double *input, double *photonAndBackground, double * scmos)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	int id2=id/sizeSubImage;
	if (id < n)
    {
		output[id]=input[id]*photonAndBackground[id2*2]+photonAndBackground[id2*2+1]+scmos[id];

    }

}







extern "C"
__global__ void vec_addPhotonsAndBackgroundMany_f (int n, int sizeSubImage,float *output, float *input, float *photonAndBackground)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	int id2=id/sizeSubImage;
	if (id < n)
    {
		output[id]=input[id]*photonAndBackground[id2*2]+photonAndBackground[id2*2+1];

    }

}






extern "C"
__global__ void vec_double2float (int n, float *output, double *input)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	if (id < n)
    {
		output[id]=(float)input[id];

    }

}


extern "C"
__global__ void vec_float2double (int n, double *output, float *input)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	if (id < n)
    {
		output[id]=(double)input[id];

    }

}




//multi kernel complexe conjugate
//*2 because real and imag parts
//n is total size (complex) of kernel divided by 2
//sizeInput is total size (complex) of image divided by 2
extern "C"
__global__ void vec_complexeConjugateKernel (int n,  int sizeInput, float *output, float *input, float *inputKernel)
{
	//n size 
	//int id = 2*(threadIdx.x + blockIdx.x * blockDim.x);
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = 2*(idy * gridDim.x * blockDim.x + idx);
	int id2=id%(sizeInput*2);
	float real;
	float imag;
	float tmp;
	if (id < n*2)
    {
		real=input[id2]/sqrt((float)sizeInput);
		imag=input[id2+1]/sqrt((float)sizeInput);
		//id : real
		//id+1 : imaginary
		tmp=imag*inputKernel[id+1]+real*inputKernel[id];
		output[id+1]=imag*inputKernel[id]-real*inputKernel[id+1];
		output[id]=tmp;
		

    }

}


//multi kernel correlation result
//n is total size (complex) of kernel divided by 2
//sizeInput is total size (complex) of image divided by 2
extern "C"
__global__ void vec_makeResultCorrelation(int n, int sizeInput, int sizeFullPadded,float *output, float *input, int *sparse)
{
	//n size 
	//int id = (threadIdx.x + blockIdx.x * blockDim.x);
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	int id2=(id)%(sizeInput);
	int id3=(id)/(sizeInput);
	if (id < n)
    {
		output[id]=input[sparse[id2]+sizeFullPadded*id3]/sqrt((float)sizeFullPadded/2.);

    }

}





extern "C"
__global__ void vec_divScalarFloat ( int n, float *result, float  *x, float  y)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int offset = idy * gridDim.x * blockDim.x + idx;

	//int offset = threadIdx.x + blockIdx.x * blockDim.x;
    if (offset < n)
    {
        result[offset] = x[offset] / y;
    }
}





//11 images as input
//25 images as output
//n=3sizesubimages
extern "C"
__global__ void vec_computeCRLB (int n,int sizeMatrix,double *output, double *input,double h)
{
	
	int sizeSubImage=n/(sizeMatrix*sizeMatrix);
	
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	if (id < n)
    {
		int p=(id/sizeSubImage);//p=0..24
		int positpix=id%sizeSubImage;//0..sizeSubImage-1
		int line=p/(sizeMatrix);//0..4
		int column=p%(sizeMatrix);//0..4
		double d1=(input[positpix+sizeSubImage*(line*2+2)]-input[positpix+sizeSubImage*(line*2+1)])/(2*h);
		double d2=(input[positpix+sizeSubImage*(column*2+2)]-input[positpix+sizeSubImage*(column*2+1)])/(2*h);



		if (input[positpix]>0){
			output[id]=(d1*d2)/(input[positpix]);
		}
		else{
			output[id]=100000000;
		}
		

    }

}







extern "C"
__global__ void vec_divCorrelation (int n, float  *x,int sizeImage, float  *varImage,float  *varPSF)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
		int lengthImage=sizeImage*sizeImage;
		int positInImage=id%lengthImage;
		int zposit=id/lengthImage;
		float div=(varImage[positInImage]*varPSF[zposit]);
		if (div>0){
			x[id] = x[id] / sqrt(div);
		}
        else{
			x[id] = -1;
		}
    }
}





/*
#include <stdio.h>

int main() {

	//vec_initFFT<<<1,1>>>(16);

  printf("toto%f\n",5.);  


}*/

































