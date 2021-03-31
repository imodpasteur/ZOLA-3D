

//#include "cublas_v2.h"
//#include "cublas.h"
 
//#include <thrust/sort.h>


//#include     "math.h"



#define ACC 40.0
#define BIGNO 1.0e10
#define BIGNI 1.0e-10




__device__ double bessi0( double x )
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function In(x) and n=0.  */
/*------------------------------------------------------------*/
{
   double ax,ans;
   double y;


   if ((ax=fabs(x)) < 3.75) {
      y=x/3.75,y=y*y;
      ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
         +y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
   } else {
      y=3.75/ax;
      ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
         +y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
         +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
         +y*0.392377e-2))))))));
   }
   return ans;
}




__device__ double bessi1( double x)
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function In(x) and n=1.  */
/*------------------------------------------------------------*/
{
   double ax,ans;
   double y;


   if ((ax=fabs(x)) < 3.75) {
      y=x/3.75,y=y*y;
      ans=ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
         +y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
   } else {
      y=3.75/ax;
      ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
         -y*0.420059e-2));
      ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
         +y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
      ans *= (exp(ax)/sqrt(ax));
   }
   return x < 0.0 ? -ans : ans;
}




__device__ double bessi( int n, double x)
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function In(x) for n >= 0*/
/*------------------------------------------------------------*/
{

   int j;
   double bi,bim,bip,tox,ans;
    
   //I added this line because when n is integer --> In(x)=I|n|(x)
   n=abs(n);
   
   if (n<0)
    return(0);
   
   if (n == 0)
      return( bessi0(x) );
   if (n == 1)
      return( bessi1(x) );


   if (x == 0.0)
      return 0.0;
   else {
      tox=2.0/fabs(x);
      bip=ans=0.0;
      bi=1.0;
      for (j=2*(n+(int) sqrt(ACC*n));j>0;j--) {
         bim=bip+j*tox*bi;
         bip=bi;
         bi=bim;
         if (fabs(bi) > BIGNO) {
            ans *= BIGNI;
            bi *= BIGNI;
            bip *= BIGNI;
         }
         if (j == n) ans=bip;
      }
      ans *= bessi0(x)/bi;
      return  x < 0.0 && n%2 == 1 ? -ans : ans;
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
__global__ void vec_subFloat (int n, float *result, float  *x, float  *y)
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
__global__ void vec_mul_fl_pow (int n, float *result, float  *x, float  *y_p,float power)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        if (y_p[id]>0){
            int count=0;
            do{
                result[id] = (float)(pow((double)y_p[id],(double)power) * (double)x[id]);
                power/=2;
                count++;
            }while((isnan(result[id]))&&(count<20));
            
            if (isnan(result[id])){
                result[id] = x[id];
                //printf("NAN value %f\n",result[id]);
            }
        }
        else{
            result[id] = x[id];
        }
        
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
__global__ void vec_addScalarFloat (int n, float *result, float  *x, float  y)
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
__global__ void vec_mulScalarFloat (int n, float *result, float  *x, float  y)
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






//put psf into image
//n=widthImage*heightImage
extern "C"         
__global__ void vec_computeModelAndLikelihood(int n,int widthImage,int heightImage,int widthPSF,int heightPSF, int numberPSF,double  *likelihood,double  *model,float  *image, float  *psf, float  *parameters,int index_x_parameter, int index_y_parameter,int index_photon_parameter,float  *background)
{

	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	
	int x=id/heightImage;
	int y=id%heightImage;
	if (id < n)
    {
        int x_start,y_start;
        int X,Y;
        model[id]=(double)background[id];//we should put bckg here
        for (int index_psf=0;index_psf<numberPSF;index_psf++){
			x_start=parameters[index_x_parameter*numberPSF+index_psf];
			X=x-x_start;
			if ((X>=0)&&(X<widthPSF)){
			    y_start=parameters[index_y_parameter*numberPSF+index_psf];
			    Y=y-y_start;
			    if ((Y>=0)&&(Y<heightPSF)){
			        model[id]+=(double)psf[index_psf*widthPSF*heightPSF+X*heightPSF+Y]*(double)parameters[index_photon_parameter*numberPSF+index_psf];
		        }
		        
		    }
		}
		
		if (model[id]>0)
			likelihood[id]=model[id]-(double)image[id]*log(model[id]);
		else
			likelihood[id]=0;

    }



}





//put psf into image
//n=widthImage*heightImage
extern "C"         
__global__ void vec_computeLikelihoodAndModelwithPhotonNumberAndBackground(int n, int numberPSF,double  *likelihood,float  *image,double  *model, float  *psf, float  *parameters,int index_photon_parameter,double *bckg)
{

	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	int index_psf=id/(n/numberPSF);
	if (id < n)
    {
        //printf("%d    %d    %d \n",index_psf,(index_photon_parameter*numberPSF+index_psf));
        model[id]=(double)psf[id]*(double)parameters[index_photon_parameter*numberPSF+index_psf]+bckg[id];
        if (model[id]>0){
		    likelihood[id]=model[id]-(double)image[id]*log(model[id]);
	    }
	    else{
	        likelihood[id]=0;
        }
    }
}



//put psf into image
//n=widthImage*heightImage
extern "C"         
__global__ void vec_computeModelwithPhotonNumber(int n, int numberPSF,double  *model, float  *psf, float  *parameters,int index_photon_parameter)
{

	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	int index_psf=id/(n/numberPSF);
	if (id < n)
    {
        //printf("%d    %d    %d \n",index_psf,(index_photon_parameter*numberPSF+index_psf));
        model[id]=(double)psf[id]*(double)parameters[index_photon_parameter*numberPSF+index_psf];
		
    }
}





//put psf into image
//n=widthImage*heightImage
extern "C"         
__global__ void vec_subtractModelwithPhotonNumber(int n, int numberPSF,double  *modelsub,double  *model, float  *psf, float  *parameters,int index_photon_parameter)
{

	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	int index_psf=id/(n/numberPSF);
	//int index_psf=id%numberPSF;
	if (id < n)
    {
        
        modelsub[id]=model[id]-((double)psf[id]*(double)parameters[index_photon_parameter*numberPSF+index_psf]);
		     
    }
}






//put psf into image
//n=widthImage*heightImage
extern "C"         
__global__ void vec_partialModel(int n,double  *result,double  *a, double  *b, float  h)
{

	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	if (id < n)
    {
        
        result[id]=(a[id]-b[id])/(double)(2*h);
		     
    }
}






//put psf into image
//n=widthImage*heightImage
extern "C"         
__global__ void vec_chi2(long n,int numPixelsPSF,double  *result,double  *model, float  *image)
{

	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	if (id < n)
    {
        if (model[id]>0){
            result[id]=(model[id]-(double)image[id])*(model[id]-(double)image[id])/(model[id]*(double)numPixelsPSF);
        }
        else{
            result[id]=1./(double)numPixelsPSF;
        }
		     
    }
}





//put psf into image
//n=widthImage*heightImage
extern "C"         
__global__ void vec_computeFisherMatrix(int n,int psfsize_square, int numberPSF,float  *fisher,double  *model, double  *modelx,double  *modely,double  *modelz,double  *modelphoton)
{

	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	
	if (id < n)
    {
        int idPSF=id/(psfsize_square*16);
	    int xyf=id%(psfsize_square*16);
	    int idFisher=xyf/psfsize_square;
	    int idPixel=xyf%psfsize_square;
	    int u=idFisher/4;
	    int v=idFisher%4;
	    int position=idPSF*psfsize_square+idPixel;
	    //printf("id  %d    idPSF %d   xyf %d     uv %d %d    idPixel %d   position %d\n",id,idPSF,xyf,u,v,idPixel,position);
	    if (model[position]>0){
	    
            fisher[id]=(float)(1/model[position]);
            if (u==0){
                fisher[id]*=(float)modelx[position];
            }
            else if (u==1){
                fisher[id]*=(float)modely[position];
            }
            else if (u==2){
                fisher[id]*=(float)modelz[position];
            }
            else if (u==3){
                fisher[id]*=(float)modelphoton[position];
            }
            
            if (v==0){
                fisher[id]*=(float)modelx[position];
            }
            else if (v==1){
                fisher[id]*=(float)modely[position];
            }
            else if (v==2){
                fisher[id]*=(float)modelz[position];
            }
            else if (v==3){
                fisher[id]*=(float)modelphoton[position];
            }
            
            
            
        }
        else{
            fisher[id]=0;
        }
		     
    }
}





//does the reverse of vec_computeModelAndLikelihood: put image in PSF stack (not useful)
extern "C"         
__global__ void vec_cropFromImage(int n,int widthPSF,int heightPSF, int numberPSF,double  *result,int widthImage,int heightImage, double  *image, float  *parameters,int index_x_parameter, int index_y_parameter)
{

	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	int index_psf=id%numberPSF;
	int xy=id/numberPSF;
	int x=xy/heightPSF;
	int y=xy%heightPSF;
	if (id < n)
    {
			int x_start=(int)parameters[index_x_parameter*numberPSF+index_psf];
			int y_start=(int)parameters[index_y_parameter*numberPSF+index_psf];
			int X=(x_start+x);
			int Y=y_start+y;
			//printf("%d  %d  / X%d   Y%d     x %d     y %d      xs%d     ys %d   xy %d     index_psf %d\n",id,(xy+index_psf*widthPSF*heightPSF),X,Y,x,y,x_start,y_start,xy,index_psf);
			if ((X>=0)&&(Y>=0)&&(X<widthImage)&&(Y<heightImage)){
			    result[xy+index_psf*widthPSF*heightPSF]=image[X*heightImage+Y];
		    }
		    else{
		        result[xy+index_psf*widthPSF*heightPSF]=0;
		    }
			

    }



}





//does the reverse of vec_computeModelAndLikelihood: put image in PSF stack (not useful)
extern "C"         
__global__ void vec_cropFromImageFloat(int n,int widthPSF,int heightPSF, int numberPSF,float  *result,int widthImage,int heightImage, float  *image, float  *parameters,int index_x_parameter, int index_y_parameter)
{

	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	int index_psf=id%numberPSF;
	int xy=id/numberPSF;
	int x=xy/heightPSF;
	int y=xy%heightPSF;
	if (id < n)
    {
			int x_start=parameters[index_x_parameter*numberPSF+index_psf];
			int y_start=parameters[index_y_parameter*numberPSF+index_psf];
			int X=(x_start+x);
			int Y=y_start+y;
			//printf("%d  %d  / X%d   Y%d     x %d     y %d      xs%d     ys %d   xy %d     index_psf %d\n",id,(xy+index_psf*widthPSF*heightPSF),X,Y,x,y,x_start,y_start,xy,index_psf);
			if ((X>=0)&&(Y>=0)&&(X<widthImage)&&(Y<heightImage)){
			    result[xy+index_psf*widthPSF*heightPSF]=image[X*heightImage+Y];
		    }
		    else{
		        result[xy+index_psf*widthPSF*heightPSF]=0;
		    }
			

    }



}





//does the reverse of vec_computeModelAndLikelihood: put image in PSF stack
extern "C"         
__global__ void vec_shiftParameter(int n,int indexParameter,float h,float *parameters)
{

	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	if (id < n)
    {
			parameters[indexParameter*n+id]+=h;
    }



}






extern "C"         
__global__ void vec_updateParameter(int n,int indexParameter,float h,double *lik1,double *lik2,double *lik3,float *parameters,float *parameterSave,float *gamma_weight,float minJump,float maxJump)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	if (id < n)
    {
        if (gamma_weight[indexParameter*n+id]>0){
			double hh=(double)h;
            double grad;
            if (abs((lik3[id]+lik1[id]-2.*lik2[id])/(hh*hh))==0){
                grad=((lik3[id]-lik1[id])*(2.*hh));
            }
            else{
                grad=((lik3[id]-lik1[id])/(2.*hh))/abs((lik3[id]+lik1[id]-2.*lik2[id])/(hh*hh));
            }
            if (grad>0){
                grad=min(abs((double)maxJump),grad);
                grad=max(abs((double)minJump),grad);
            }
            else{
                grad=max(-abs((double)maxJump),grad);
                grad=min(-abs((double)minJump),grad);
            }
            //if ((indexParameter==2)&&(id==0))
            //    printf("grad   %d     %d    %f     lik  %f    %f    %f   %f\n",indexParameter,id,grad,lik1[id],lik2[id],lik3[id],gamma_weight[indexParameter*n+id]);
            parameterSave[id]=parameters[indexParameter*n+id];
            parameters[indexParameter*n+id]-=(float)grad*gamma_weight[indexParameter*n+id];
        }
            
    }



}




extern "C"         
__global__ void vec_checkLikelihood(int n,int indexParameter,double *likold,double *liknew,float *parameters,float * parameterSave,float *gamma_weight)
{

	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	if (id < n)
    {
			if (liknew[id]>likold[id]){//if new likelihood not better: change weight and put back weight
                parameters[indexParameter*n+id]=parameterSave[id];
                gamma_weight[indexParameter*n+id]/=10.;
            }
            else{//if new likelihood better: no change anymore
                likold[id]=liknew[id];
                gamma_weight[indexParameter*n+id]=0;
                parameterSave[id]=parameters[indexParameter*n+id];
            }
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
__global__ void vec_addFloat(int n, float *result, float  *a, float *b)
{

    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = a[id] + b[id];
    }
}



extern "C"
__global__ void vec_addanddivide(int n, float *result, float  *num, float *div, float *added2div)
{//perform num/(div+added2div)

    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        div[id]+=added2div[id];
        result[id] = num[id] / (div[id]);
    }
}








extern "C"
__global__ void vec_computeLikelihoodDeconvolution(int n, float *result, float  *I, float *M)
{//perform num/(div+added2div)

    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = M[id]-I[id]*log(M[id]);
    }
}







__device__ double bessiRatio( int n, double x)
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function In-1(x)/In(x) for n >= 0*/
/*------------------------------------------------------------*/
{
    n=abs(n);
    double Y=(x/2)*(x/2);
    double B_up=1;
    double BB_up=1;
    double B_down=1;
    double BB_down=1;
    double k=0;
    double n_up=n-1;
    double n_down=n;
    double ratio=1;
    double iter=sqrt(ACC*n_down);
    int nbIter=0;
    iter=5;
    for (k=1;k<iter;k++){
        B_up*=Y/(k*(k+n_up));
        BB_up+=B_up;
        B_down*=Y/(k*(k+n_down));
        BB_down+=B_down;
        if (B_up>1e+300){
            //break;
        }
        ratio=BB_up/BB_down;
        nbIter++;
    }
    //printf("BB %f  %f   /   %f  %f  / Y %f   ratio : %f  it:%d\n",n_up,n_down,B_up,B_down,Y,ratio,nbIter);
    printf("BB %f  %f   /   %f         bu:%f bd:%f r:%f \n",n_down,x,(n_down*2/x)*ratio,BB_up,BB_down,ratio);
    return (n_down*2/x)*ratio;

}





//See: On the Computation of Modified Bessel Function Ratios, Mathematics of Computation, September 1978
__device__ double R( double * n, double * x,double k,double * K){
    
    if (k>=*K){
        return 0.;
    }
    else{
        return 1./( (2./ *x)*(*n+k) + R(n,x,k+1,K));
    }
    
}

__device__ double myBessiRatio( double n, double x)
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function In-1(x)/In(x) for n >= 0*/
/*------------------------------------------------------------*/
{
    n=abs(n);
    
    double iterNumber=20;
    
    return 1./R(&n,&x,0,&iterNumber);

}


__device__ double myBessiRatioNonRecusrsive( double n, double x)
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function In-1(x)/In(x) for n >= 0*/
/*------------------------------------------------------------*/
{
    
    
    double som=0;
    for (double k=50;k>=0;k--){
        som+=((2./x)*(n+k));
        som=1./som;
    }
    return 1./som;
}


//Even if it is float, we assume I is photon count here --> we cast it to integer
//Perform Ib-Ia skellam
// n should be full size
// minReplacement = 1 by default: When I<<M --> b1==inf  --> inf value So we replace "inf" values by "minReplacement" in result
// maxReplacement = 1 by default: When I>>M --> b1==0  --> division Nan So we replace "Nan" values by "maxReplacement" in result
extern "C"
__global__ void vec_skellam_order1(long n, long imageSize,float *result, float  *I, float *M)
{//perform num/(div+added2div)

    long idx = threadIdx.x + blockIdx.x * blockDim.x;
	long idy = threadIdx.y + blockIdx.y * blockDim.y;
	long id = idy * gridDim.x * blockDim.x + idx;
	int indexFrame=id/imageSize;
	long indexImage=id%imageSize;
	float sqrtMaMb,tmp,Itmp,tt;
	//float Mb,Ma;
	//float Ib,Ia;
	//double b0,b1;
	
    if (id < n)
    {
    
        if (id<imageSize){
            result[id+n*2]=0;
            /*Itmp=I[indexImage+imageSize]-I[indexImage];
            result[id+n*2]=Itmp/(2*M[indexImage+imageSize]);*/
        }
        if (id>=(n-imageSize)){
            result[id+n]=0;
            /*Itmp=I[indexImage+imageSize]-I[indexImage];
            result[id+n]=-Itmp/(2*M[indexImage]);*/
        }
        //Ma=M[indexImage];
        //Mb=M[indexImage+imageSize];
        //Ia=I[indexImage];
        //Ib=I[indexImage+imageSize];
        int index0=indexImage+imageSize*(indexFrame);//current frame
        int index1=indexImage+imageSize*(indexFrame+1);//next frame
        if (id<(n-imageSize)){// because we do it nbFrame-1 times
            Itmp=I[index1]-I[index0];
            sqrtMaMb=(float)sqrt(M[index0]*M[index1]);
            tt=sqrtMaMb*2;
            //printf("test %f   %f   \n",Itmp,tt);
            
            tmp= myBessiRatioNonRecusrsive((double)abs(Itmp),(double)tt)    ;
            
            /*//here, we compute skellam :
            b1=bessi( (int)Itmp, (double)tmp);
            b0=bessi( (int)(Itmp-1), (double)tmp);
            //printf("b0 %f   b1 %f   tmp %f Itmp %f     sqrtMaMb %f     %f   bool: %d\n",b0,b1,tmp,Itmp,sqrtMaMb,1e+30,(b1>1e+30));
            if (b1<1e-307){//10E-30
                //result[offset+id+n*2]=maxReplacement;
                //result[offset+id+n*3]=maxReplacement;
                result[id+n]=maxReplacement;
                result[id+n*2+imageSize]=maxReplacement;
            }
            else if (b0>1e+307){//10E+36
                //result[offset+id+n*2]=minReplacement;
                //result[offset+id+n*3]=minReplacement;
                result[id+n]=minReplacement;
                result[id+n*2+imageSize]=minReplacement;
            }
            else{
                tmp=(double)((double)1./(double)sqrtMaMb)*( b0  -  ( (double)(Itmp/tmp) * b1 ) )/b1;*/
                
                //result[offset+id+n*2]=(float)(M[indexImage+imageSize]*tmp);
                //result[offset+id+n*3]=(float)(M[indexImage]*tmp);
                result[id+n]=(float)(M[index1]*( tmp- (double)(abs(Itmp)/tt) ))/(double)sqrtMaMb;
                result[id+n*2+imageSize]=(float)(M[index0]*( tmp- (double)(abs(Itmp)/tt) )/(double)sqrtMaMb);
            //}
            
            //result[offset+id]=-Itmp/(2*M[indexImage]);//  (Ia-Ib)/2Ma
            //result[offset+id+n]=Itmp/(2*M[indexImage+imageSize]);//  (Ib-Ia)/2Mb
            result[id+n]+=-Itmp/(2*M[index0]);//  (Ia-Ib)/2Ma
            result[id+n*2+imageSize]+=Itmp/(2*M[index1]);//  (Ib-Ia)/2Mb
            
            if (index0/imageSize!=0){//first image
                result[id+n]/=2;
                
            }
            if (index1/imageSize!=((n/imageSize)-1)){
                result[id+n*2+imageSize]/=2;
            }
            
        }
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





extern "C"
__global__ void vec_computeGaussianLikelihood (int n, double *result, double *image, double  *model)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	if (id < n)
    {
		if (model[id]>0)
			result[id]=(model[id]-image[id])*(model[id]-image[id]);
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
__global__ void vec_shrink (int n,float  *output, float  *input, float  threshold)
{
    
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	if (id < n)
    {
			if (input[id]<threshold){
			    output[id]=0;
			}
			else{
			    output[id]=input[id];
			}
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



extern "C"
__global__ void vec_complexeConjugateKernelSubtract (int n,  int sizeInput, float *output, float *input, float *inputKernel)
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





//multi kernel complexe conjugate
//*2 because real and imag parts
//n is total size (complex) of kernel divided by 2
//sizeInput is total size (complex) of image divided by 2
extern "C"
__global__ void vec_complexeMulKernel (int n,  int sizeInput, float *output, float *input, float *inputKernel)
{
	//n size 
	//int id = 2*(threadIdx.x + blockIdx.x * blockDim.x);
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = 2*(idy * gridDim.x * blockDim.x + idx);
	int id2=id%(n*2);//strange...id2 might be = to id
	float real;
	float imag;
	float tmp;
	if (id < n*2)
    {
		real=input[id2]/sqrt((float)sizeInput);
		imag=input[id2+1]/sqrt((float)sizeInput);
		//id : real
		//id+1 : imaginary
		tmp=real*inputKernel[id]-imag*inputKernel[id+1];
		output[id+1]=imag*inputKernel[id]+real*inputKernel[id+1];
		output[id]=tmp;
		

    }

}





//multi kernel complexe conjugate
//*2 because real and imag parts
//n is total size (complex) of kernel divided by 2
//sizeInput is total size (complex) of image divided by 2
extern "C"
__global__ void vec_complexeMulKernelMany (long n,  long depth,  long size, float *output, float *input, float *inputKernel)
{
	//n size 
	//int id = 2*(threadIdx.x + blockIdx.x * blockDim.x);
	long idx = threadIdx.x + blockIdx.x * blockDim.x;
	long idy = threadIdx.y + blockIdx.y * blockDim.y;
	long idd = 2*(idy * gridDim.x * blockDim.x + idx);
	long id=idd%(2*size*depth);
	
	//int id2=id%(n*2);
	float real;
	float imag;
	float tmp;
	if (idd < n*2)
    {
        
		real=input[idd]/size;
		imag=input[idd+1]/size;
		//id : real
		//id+1 : imaginary
		tmp=real*inputKernel[id]-imag*inputKernel[id+1];
		output[idd+1]=imag*inputKernel[id]+real*inputKernel[id+1];
		output[idd]=tmp;
		

    }

}





extern "C"
__global__ void vec_copyMany(long n, long sizeInput,long depth,long nbFrame, float *output, float *input)
{
	//n size 
	//int id = (threadIdx.x + blockIdx.x * blockDim.x);
	long idx = threadIdx.x + blockIdx.x * blockDim.x;
	long idy = threadIdx.y + blockIdx.y * blockDim.y;
	long idd = idy * gridDim.x * blockDim.x + idx;
	long id=idd%(sizeInput*depth);
	long idFrame=idd/(sizeInput*depth);
	long idPixel=(id)%(sizeInput);
	//long idZ=(id)/(sizeInput);
	
	if (idd < n)
    {
        long posit=idPixel+idFrame*sizeInput;
		output[idd]=input[posit];
        
        
    }

}




//multi kernel correlation result
//n is total size (complex) of kernel divided by 2
//sizeInput is total size (complex) of image divided by 2
extern "C"
__global__ void vec_mycusparsemoduloSsctrMany(long n, long nbFrame,long depth,long sizeShift, long width,long height,float *output, float *input, int *sparse)
{
	//n size 
	//int id = (threadIdx.x + blockIdx.x * blockDim.x);
	long idx = threadIdx.x + blockIdx.x * blockDim.x;
	long idy = threadIdx.y + blockIdx.y * blockDim.y;
	long idd = idy * gridDim.x * blockDim.x + idx;
	long sizeSparse=width*height;
	long idFrame=idd/(depth*width*height);
	long id=idd%(depth*width*height);
	long id2=(id)%(sizeSparse);
	long id3=(id)/(sizeSparse);
	if (idd < n)
    {
        
		output[(long)sparse[id2]+sizeShift*id3+idFrame*sizeShift*depth]=input[idd];

    }

}





//multi kernel correlation result
//n is total size (complex) of kernel divided by 2
//sizeInput is total size (complex) of image divided by 2
//process many images at the same times
extern "C"
__global__ void vec_mycusparsemoduloSsctr(int n, int sizeShift, int sizeSparse,float *output, float *input, int *sparse)
{
	//n size 
	//int id = (threadIdx.x + blockIdx.x * blockDim.x);
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	int id2=(id)%(sizeSparse);
	int id3=(id)/(sizeSparse);
	if (id < n)
    {
    
		output[sparse[id2]+sizeShift*id3]=input[id];

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


//multi kernel correlation result
//n is total size (complex) of kernel divided by 2
//sizeInput is total size (complex) of image divided by 2
extern "C"
__global__ void vec_normalizeCorrelation(long n, long nbFrame, long depth, long sizeInput,float *output, float *input, float *divide)
{
	//n size 
	//int id = (threadIdx.x + blockIdx.x * blockDim.x);
	long idx = threadIdx.x + blockIdx.x * blockDim.x;
	long idy = threadIdx.y + blockIdx.y * blockDim.y;
	long idd = idy * gridDim.x * blockDim.x + idx;
	long id=idd%(sizeInput*depth);
	//long idFrame=idd/(sizeInput*depth);
	//long id2=(id)%(sizeInput);
	//long id3=(id)/(sizeInput);
	if (idd < n)
    {
        
		output[idd]=input[idd]/divide[id];

    }

}




//multi kernel correlation result
//n is total size (complex) of kernel divided by 2
//sizeInput is total size (complex) of image divided by 2
extern "C"
__global__ void vec_makeResultCorrelationMany(long n, long nbFrame, long depth, long sizeInput, long sizeFullPadded,float *output, float *input, int *sparse)
{
	//n size 
	//int id = (threadIdx.x + blockIdx.x * blockDim.x);
	long idx = threadIdx.x + blockIdx.x * blockDim.x;
	long idy = threadIdx.y + blockIdx.y * blockDim.y;
	long idd = idy * gridDim.x * blockDim.x + idx;
	long id=idd%(sizeInput*depth);
	long idFrame=idd/(sizeInput*depth);
	long id2=(id)%(sizeInput);
	long id3=(id)/(sizeInput);
	if (idd < n)
    {
        
		output[idd]=input[sparse[id2]+sizeFullPadded*id3+idFrame*sizeFullPadded*depth];

    }

}




//multi kernel correlation result
//n is total size (complex) of kernel divided by 2
//sizeInput is total size (complex) of image divided by 2
extern "C"
__global__ void vec_turnMatrixMany(long n, long nbFrame, long depth, long sizeInput,float *output, float *input)
{
	//n size 
	//int id = (threadIdx.x + blockIdx.x * blockDim.x);
	long idx = threadIdx.x + blockIdx.x * blockDim.x;
	long idy = threadIdx.y + blockIdx.y * blockDim.y;
	long idd = idy * gridDim.x * blockDim.x + idx;
	long id=idd%(sizeInput*depth);
	long idFrame=idd/(sizeInput*depth);
	long idPixel=(id)%(sizeInput);
	long idZ=(id)/(sizeInput);
	
	if (idd < n)
    {
        long posit=idPixel+idFrame*sizeInput+idZ*sizeInput*nbFrame;
		output[posit]=input[idd];
        
        
    }

}


extern "C"
__global__ void vec_makeResultCorrelationNormalized(int n, int sizeInput, int sizeFullPadded,float *output, float *input, int *sparse,float divide,float* device_divide,float minValue)
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
        float tmp=(input[sparse[id2]+sizeFullPadded*id3]/(divide*sqrt((float)sizeFullPadded/2.)))-device_divide[id2];
        if (tmp>0){
		    output[id]=minValue;//tmp;
	    }
	    else{
	        output[id]=minValue;
	    }

    }

}



extern "C"
__global__ void vec_initializeDeconvolution(int n,int nbpsf ,float *o, float *op, float *m,float value)
{
	//n size 
	//int id = (threadIdx.x + blockIdx.x * blockDim.x);
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	if (id < n)
    {
        o[id]=value;
        op[id]=value;
        if (id<n/nbpsf){
            m[id]=value*(float)nbpsf;
        }

    }

}




extern "C"
__global__ void vec_initializeVectorToValue(int n ,float *v, float value)
{
	//n size 
	//int id = (threadIdx.x + blockIdx.x * blockDim.x);
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	if (id < n)
    {
        v[id]=value;
        

    }

}

extern "C"
__global__ void vec_chiScore (int n, float *result, float  *image, float  *model)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
    if (id < n)
    {
        result[id] = (image[id] - model[id])*(image[id] - model[id])/model[id];
    }
}



extern "C"
__global__ void vec_max(int n ,float *v, float value)
{
	//n size 
	//int id = (threadIdx.x + blockIdx.x * blockDim.x);
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	if (id < n)
    {
        v[id]=max(v[id],value);
        

    }

}




extern "C"
__global__ void vec_subtractMeanWithSumAsInputWithPositiveConstraint(int n, float *output, float *input, float *sum,float minValue)
{
	//n size 
	//int id = (threadIdx.x + blockIdx.x * blockDim.x);
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	
	if (id < n)
    {
        float tmp=input[id]-(sum[0]/(float)n);
        if (tmp>0){
		    output[id]=tmp;
	    }
	    else{
	        output[id]=minValue;
	    }

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







extern "C"
__global__ void vec_updateMandOP ( int n, float *m,float *op, float  *x, float  div)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int offset = idy * gridDim.x * blockDim.x + idx;
    float tmp;
	//int offset = threadIdx.x + blockIdx.x * blockDim.x;
    if (offset < n)
    {
        tmp  = m[offset] -  op[offset];
        op[offset] = (x[offset] / div);
        m[offset] = tmp + op[offset];
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






//n should be size(input)/2
/*extern "C"
__global__ void vec_sortRows(int n,float * value,int *index) {
	
	
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = 2*(idy * gridDim.x * blockDim.x + idx);
	float valtmp;
	int indtmp;
	int i;
	if (id <n*2)
    {
        
        for (i=0;i<1+(n/2);i++){//TO DO
		    if (value[id]>value[id+1]){//if yes -> change 
		        valtmp=value[id];
		        indtmp=index[id];
		        value[id]=value[id+1];
		        index[id]=index[id+1];
		        value[id+1]=valtmp;
		        index[id+1]=indtmp;
		    }
	    }

	}
}*/




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








//perform M=\sum_z(o_z*psf_z)+bckg, where * corresponds to convolution operator
extern "C"
__global__ void vec_manualFilteringTest (long n, int sizeImageX, int sizeImageY, int sizePSF,int nbPSF,float *m, float *o,  float *psf, float  *bckg,float *tmp)
{
	long idx = threadIdx.x + blockIdx.x * blockDim.x;
	long idy = threadIdx.y + blockIdx.y * blockDim.y;
	long id = idy * (long)gridDim.x * (long)blockDim.x + idx;
	
	long p_psf=(long)(id%(long)(sizePSF*sizePSF));//reste
	//long x_psf=p_psf/sizePSF;
	//long y_psf=p_psf%sizePSF;
	
	long p_image=(long)(id/(long)(sizePSF*sizePSF));//position
	long x_image=p_image/sizeImageY;
	long y_image=p_image%sizeImageY;
	
	
	if (id < n)
    {
        //initialize m to 0
		if (p_psf==0){
		    m[p_image]=0;
		}
		
        int t=0;//loop t among nb_psf
        
		
		
		//compute convolution for one single image
		tmp[id]=o[p_image+t*sizeImageX*sizeImageY]*psf[p_psf];
		//printf("0\n");
		
		////////////////////////////////////////////IMPOSSIBLE: __syncthreads does not synchronise blocks/grids
		__syncthreads();
		
		
		//printf("1\n");
		//p_image == x_image*sizeImage+y_image
		if (p_image!=x_image*sizeImageY+y_image){
		printf("ZUT %d  %d  %d\n",p_image,x_image,y_image,sizeImageY);
		}
		if (p_psf==0){
		    for (long i=0;i<sizePSF;i++){
		        for (long ii=0;ii<sizePSF;ii++){
		            //p_image==((x_image)*sizeImageY+y_image);
		            long nextY=(1+sizePSF*sizePSF)*ii+sizePSF*i;
		            long nextX=sizeImageY*sizePSF*sizePSF*i;
		            long index=p_image*sizePSF*sizePSF + nextY + nextX;
		            if (((x_image+i)<sizeImageX)&&((y_image+ii)<sizeImageY)){
		                if (index<n){
		                    if (p_image==0){
		                        //printf("%d  %d  %d  %d      x_  %d  %d   size: %d  %d  \n",i,ii,index,p_image,x_image,y_image,sizeImageY,sizePSF);
	                        }
	                        
	                        m[p_image]+=tmp[index];
                        }
                        else{
                             //printf("ZUT %d  %d  %d %d  %d\n",p_image,x_image,y_image,sizeImageY,index,n);
                        }
                    }
                    else{m[p_image]=-.1;}
                    
                }
	        }
		}
		__syncthreads();
		//for (int i=0;i<sizePSF;i++){
		/*for (int i=0;i<1;i++){
		    int power=1;

		    for (int ii=0;ii<i;ii++){

		        power*=2;

		    }

		    if ((id)%(2*power)==0){//if middle of convol kernel is even
		        int pp=power;
		        if (power

		        int index=id+power*sizePSF*sizePSF+power;

		        if (index<n){

		            //tmp[id]+=tmp[index];
		            tmp[id]=(float)index;

	            }

	        }

		}*/
			
    }



}




//convolution: apply local min
extern "C"
__global__ void vec_localMinimum(int n, int sizeImageX, int sizeImageY, int sizePSF,float *result, float *image)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	
	
	
	int x_image=id/sizeImageY;
	int y_image=id%sizeImageY;
	
	int px;
	int py;
	int i;
	int ii;
	int dist;
	if (id < n)
    {
        //initialize m to 0
        result[id]=image[id];
	    
		
		
		//search min
	    for (i=0,px=x_image-sizePSF/2;i<sizePSF;i++,px++){
	        for (ii=0,py=y_image-sizePSF/2;ii<sizePSF;ii++,py++){
	            
	            if ((px<sizeImageX)&&(py<sizeImageY)&&(px>=0)&&(py>=0)){
	                
                    dist=sqrt((float)(((sizePSF/2)-i)*((sizePSF/2)-i)+((sizePSF/2)-ii)*((sizePSF/2)-ii)));
	                if (dist<=sizePSF/2){
                        result[id]=fmin(image[px*sizeImageY+py],result[id]);
                    }
                    
                }
            }
        }
        
        
        
    }



}






//convolution: apply local max
extern "C"
__global__ void vec_localMaximum(int n, int sizeImageX, int sizeImageY, int sizePSF,float *result, float *image)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	
	
	
	int x_image=id/sizeImageY;
	int y_image=id%sizeImageY;
	
	int px;
	int py;
	int i;
	int ii;
	int dist;
	
	if (id < n)
    {
        //initialize m to 0
        result[id]=image[id];
	    
		
		
		//search max
	    for (i=0,px=x_image-sizePSF/2;i<sizePSF;i++,px++){
	        for (ii=0,py=y_image-sizePSF/2;ii<sizePSF;ii++,py++){
	            
	            if ((px<sizeImageX)&&(py<sizeImageY)&&(px>=0)&&(py>=0)){
	                dist=sqrt((float)(((sizePSF/2)-i)*((sizePSF/2)-i)+((sizePSF/2)-ii)*((sizePSF/2)-ii)));
	                if (dist<sizePSF/2){
                        result[id]=fmax(image[px*sizeImageY+py],result[id]);
                    }
                    
                }
            }
        }
        
        
        
    }



}




//convolution: apply local mean
extern "C"
__global__ void vec_localMean(int n, int sizeImageX, int sizeImageY, int sizePSF,float *result, float *image)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	
	
	
	int x_image=id/sizeImageY;
	int y_image=id%sizeImageY;
	
	int px;
	int py;
	int i;
	int ii;
	int dist;
	float mean;
	float count;
	if (id < n)
    {
        //initialize m to 0
        mean=0;
	    count=0;
		
		
		//search max
	    for (i=0,px=x_image-sizePSF/2;i<sizePSF;i++,px++){
	        for (ii=0,py=y_image-sizePSF/2;ii<sizePSF;ii++,py++){
	            
	            if ((px<sizeImageX)&&(py<sizeImageY)&&(px>=0)&&(py>=0)){
	                dist=sqrt((float)(((sizePSF/2)-i)*((sizePSF/2)-i)+((sizePSF/2)-ii)*((sizePSF/2)-ii)));
	                if (dist<=sizePSF/2){
                        mean+=image[px*sizeImageY+py];
                        count++;
                    }
                    
                }
            }
        }
        result[id]=mean/count;
        
        
    }



}




//convolution: fast if psfsize roughly < 6x6
extern "C"
__global__ void vec_manualFiltering (int n, int sizeImageX, int sizeImageY, int sizePSF,float *m, float *o,  float *psf)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int id = idy * gridDim.x * blockDim.x + idx;
	
	
	
	int x_image=id/sizeImageY;
	int y_image=id%sizeImageY;
	
	int count;
	int px;
	int py;
	int i;
	int ii;
	if (id < n)
    {
        //initialize m to 0
	    m[id]=0;
		
        
		
		
		//compute convolution for one single image
		
		
		count=0;
	    for (i=0,px=x_image-sizePSF/2;i<sizePSF;i++,px++){
	        for (ii=0,py=y_image-sizePSF/2;ii<sizePSF;ii++,py++){
	            
	            if ((px<sizeImageX)&&(py<sizeImageY)&&(px>=0)&&(py>=0)){
                    m[id]+=o[px*sizeImageY+py]*psf[i*sizePSF+ii];
                    count++;
                }
            }
        }
        m[id]/=(float)count;
    }



}







//perform M=\sum_z(o_z*psf_z)+bckg, where * corresponds to convolution operator
extern "C"
__global__ void vec_manualFilteringStacked (int n, int sizeImageX, int sizeImageY, int sizePSF,int nbPSF,float *res, float *o,  float *psf)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int idy = threadIdx.y + blockIdx.y * blockDim.y;
	int ids = idy * gridDim.x * blockDim.x + idx;
	
	int s=ids%nbPSF;
	int id=ids/nbPSF;
	int x_image=id/sizeImageY;
	int y_image=id%sizeImageY;
	
	int count;
	int px;
	int py;
	int i;
	int ii;
	int shift=s*sizeImageY*sizeImageX;
    
	if (id < n)
    {
        
        //initialize m to 0
        
	    res[ids]=0;
		
		//compute convolution for one single image
		
		count=0;
	    for (i=0,px=x_image-sizePSF/2;i<sizePSF;i++,px++){
	        for (ii=0,py=y_image-sizePSF/2;ii<sizePSF;ii++,py++){
	            
	            if ((px<sizeImageX)&&(py<sizeImageY)&&(px>=0)&&(py>=0)){
                    res[ids]+=o[shift+px*sizeImageY+py]*psf[i*sizePSF+ii];
                    count++;
                }
            }
        }
        res[ids]/=(float)count;
    }
}










//convolution
extern "C"
__global__ void vec_manualFilteringFast (long n, long sizeImageX, long sizeImageY, long sizePSF,float *m, float *o,  float *psf)
{
	long idx = threadIdx.x + blockIdx.x * blockDim.x;
	long idy = threadIdx.y + blockIdx.y * blockDim.y;
	long id = idy * gridDim.x * blockDim.x + idx;
	
	
	
	long p_psf=(long)(id%(long)(sizePSF*sizePSF));//reste
	//long x_psf=p_psf/sizePSF;
	//long y_psf=p_psf%sizePSF;
	
	long p_image=(long)(id/(long)(sizePSF*sizePSF));//position
	long x_image=p_image/sizeImageY;
	long y_image=p_image%sizeImageY;
	
	
	//long p_image=id/(sizePSF*sizePSF);
	//long p_psf=id%(sizePSF*sizePSF);
	
	
	//int x_image=p_image/sizeImageY;
	//int y_image=p_image%sizeImageY;
	int x_psf=(int)p_psf/sizePSF;
	int y_psf=(int)p_psf%sizePSF;
	
	
	
	if (id < n)
    {   
        
        int xx_im=x_image+x_psf-sizePSF/2;
        int yy_im=y_image+y_psf-sizePSF/2;
        long o_posit=(xx_im)*sizeImageY+(yy_im);
        //long m_posit=(p_image)*sizePSF*sizePSF+p_psf   ==   id ;
        
        if ((xx_im<sizeImageX)&&(yy_im<sizeImageY)&&(xx_im>=0)&&(yy_im>=0)){
            
            printf("m_posit:%ld   p_im:%ld   p_psf:%ld  x_psf:%d    y_psf:%d        o_posit:%ld      n:%ld    v_im:%f   v_psf:%f\n",id,p_image,p_psf,x_psf,y_psf,o_posit,n,o[o_posit],psf[p_psf]);
            
            //initialize m to 0
	        m[id]=o[o_posit]*psf[sizePSF*sizePSF - (int)p_psf - 1];
	    }
		
    }



}







//convolution
extern "C"
__global__ void vec_manualFilteringStackedFast (long n, long sizeImageX, long sizeImageY, long sizePSF, long nbPSF,float *m, float *o,  float *psf)
{
	long idx = threadIdx.x + blockIdx.x * blockDim.x;
	long idy = threadIdx.y + blockIdx.y * blockDim.y;
	long ids = idy * gridDim.x * blockDim.x + idx;
	
	long id = ids/nbPSF;
	long s = ids%nbPSF;
	
	long p_psf=(long)(id%(long)(sizePSF*sizePSF));//reste
	//long x_psf=p_psf/sizePSF;
	//long y_psf=p_psf%sizePSF;
	
	long p_image=(long)(id/(long)(sizePSF*sizePSF));//position
	long x_image=p_image/sizeImageY;
	long y_image=p_image%sizeImageY;
	
	
	//long p_image=id/(sizePSF*sizePSF);
	//long p_psf=id%(sizePSF*sizePSF);
	
	
	//int x_image=p_image/sizeImageY;
	//int y_image=p_image%sizeImageY;
	int x_psf=(int)p_psf/sizePSF;
	int y_psf=(int)p_psf%sizePSF;
	
	
	
	if (id < n)
    {   
        
        int xx_im=x_image+x_psf-sizePSF/2;
        int yy_im=y_image+y_psf-sizePSF/2;
        long o_posit=s*sizeImageX*sizeImageY + (xx_im)*sizeImageY+(yy_im);
        //long m_posit=(p_image)*sizePSF*sizePSF+p_psf   ==   id ;
        
        if ((xx_im<sizeImageX)&&(yy_im<sizeImageY)&&(xx_im>=0)&&(yy_im>=0)){
            
            //printf("m_posit:%ld   p_im:%ld   p_psf:%ld  x_psf:%d    y_psf:%d        o_posit:%ld      n:%ld    v_im:%f   v_psf:%f\n",id,p_image,p_psf,x_psf,y_psf,o_posit,n,o[o_posit],psf[p_psf]);
            
            //initialize m to 0
	        m[ids]=o[o_posit]*psf[sizePSF*sizePSF - (int)p_psf - 1];
	    }
		
    }



}









































