/*
 * code with Licence MIT
 */
package org.pasteur.imagej.cuda;
import ij.IJ;
import jcuda.runtime.cudaStream_t;
import jcuda.CudaException;
import jcuda.Pointer;
import jcuda.driver.CUstream;
import jcuda.jcufft.cufftHandle;

/**
 * Vector operations for double-precision floating point vectors. <br>
 * <br>
 * Before any of the vector methods of this class can be called, it
 * has to be initialized by calling {@link #init()}. In order to 
 * release all resources allocated by this class, {@link #shutdown()}
 * has to be called.
 * <br>
 * NOTE: This class only forms a thin layer around the actual CUDA 
 * kernel calls. It can not verify the function parameters. The caller
 * is responsible for giving appropriate pointers and vector sizes.
 */
public class MyVecDouble 
{
    /**
     * The {@link VecKernels} that maintains the kernels that can be
     * called via the methods of this class.
     */
    private static MyVecKernels vecKernels = null;
    
    /**
     * Private constructor to prevent instantiation
     */
    private MyVecDouble()
    {
        // Private constructor to prevent instantiation
    }
    
    /**
     * Initialize this class. This method has to be called before any
     * of the vector operation methods of this class can be called.
     * The resources that are allocated by this call may be freed
     * by calling {@link #shutdown()}.
     * 
     * @throws CudaException If the class can not be initialized. 
     * Reasons for this may be (but are not limited to) : <br>
     * <ul>
     *   <li> 
     *     It is neither possible to attach to an existing context, 
     *     nor to create a new context.
     *   </li>
     *   <li> 
     *     The resource that contains the kernels (for example, a
     *     PTX file) can not be loaded
     *   </li>
     * </ul>
     */
    public static void init()
    {
        //if (vecKernels == null){
            //IJ.log("init");
            shutdown();
            String kernelNamePrefix = "vec_";
            String kernelNameType = "double";
            String kernelNameSuffix = "";
            vecKernels = new MyDefaultVecKernels(
                kernelNameType, kernelNamePrefix, kernelNameSuffix);
        //}
    }
    
    
    
    /**
     * Perform a shutdown and release all resources allocated by this class.
     */
    public static void shutdown()
    {
        
        if (vecKernels != null)
        {
            //IJ.log("shut");
            vecKernels.shutdown();
            vecKernels = null;
        }
    }
    
    
    
    /**
     * Passes the call to the {@link VecKernels#call(String, Object...)}
     * 
     * @param name The kernel name
     * @param workSize The work size for the kernel call
     * @param arguments The kernel arguments (including the vector size)
     */
    private static void callSingle(String name, CUstream idStream, long workSize, Object ... arguments)
    {
        if (vecKernels == null)
        {
            throw new CudaException(
                "Kernels not initialized. Call init() first.");
        }
        vecKernels.callSingle(name, workSize, idStream,arguments);
    }
    
    /**
     * Passes the call to the {@link VecKernels#call(String, Object...)}
     * 
     * @param name The kernel name
     * @param workSize The work size for the kernel call
     * @param arguments The kernel arguments (including the vector size)
     */
    private static void call(String name, CUstream idStream, long workSize, Object ... arguments)
    {
        if (vecKernels == null)
        {
            throw new CudaException(
                "Kernels not initialized. Call init() first.");
        }
        vecKernels.call(name, workSize, idStream,arguments);
    }
    
    
    
    
    
    /**
     * Passes the call to the {@link VecKernels#call(String, Object...)}
     * 
     * @param name The kernel name
     * @param workSize The work size for the kernel call
     * @param arguments The kernel arguments (including the vector size)
     */
    private static void call2DGrid(String name, CUstream idStream, long workSize, Object ... arguments)
    {
        if (vecKernels == null)
        {
            throw new CudaException(
                "Kernels not initialized. Call init() first.");
        }
        vecKernels.call2DGrid(name, workSize, idStream,arguments);
    }
    
    
    /**
     * Passes the call to the {@link VecKernels#call(String, Object...)}
     * 
     * @param name The kernel name
     * @param workSize The work size for the kernel call
     * @param arguments The kernel arguments (including the vector size)
     */
    private static void callWithSharedMemory(String name, CUstream idStream, long workSize, Object ... arguments)
    {
        if (vecKernels == null)
        {
            throw new CudaException(
                "Kernels not initialized. Call init() first.");
        }
        vecKernels.callWithSharedMemory(name, workSize, idStream,arguments);
    }
    
    
    
    /**
     * Set all elements of the given vector to the given value
     * 
     * @param n The size of the vector
     * @param result The vector that will store the result
     * @param value The value to set
     */
    public static void set(CUstream idStream,long n, Pointer result, double value)
    {
        call2DGrid("set", idStream,n, n, result, value);
    }
    
    //=== Vector arithmetic ==================================================
    
    /**
     * Add the given vectors.
     * 
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     * @param y The second vector
     */
    public static void add(CUstream idStream,long n, Pointer result, Pointer x, Pointer y)
    {
        call2DGrid("add", idStream, n, n, result, x, y);
    }

    /**
     * Subtract the given vectors.
     * 
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     * @param y The second vector
     */
    public static void sub(CUstream idStream,long n, Pointer result, Pointer x, Pointer y)
    {
        call2DGrid("sub", idStream, n, n, result, x, y);
    }

    /**
     * Multiply the given vectors.
     * 
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     * @param y The second vector
     */
    public static void mul(CUstream idStream,long n, Pointer result, Pointer x, Pointer y)
    {
        call2DGrid("mul", idStream, n, n, result, x, y);
    }
    
    
    /**
     * Multiply the given vectors.
     * 
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     * @param y The second vector
     */
    public static void mul_fl(CUstream idStream,long n, Pointer result, Pointer x, Pointer y)
    {
        call2DGrid("mul_fl", idStream, n, n, result, x, y);
    }

    /**
     * Divide the given vectors.
     * 
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     * @param y The second vector
     */
    public static void div(CUstream idStream,long n, Pointer result, Pointer x, Pointer y)
    {
        call2DGrid("div", idStream, n, n, result, x, y);
    }

    /**
     * Negate the given vector.
     * @param n The size of the vectors
     * 
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void negate(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("negate", idStream, n, n, result, x);
    }
    
    //=== Vector-and-scalar arithmetic =======================================
    
    /**
     * Add the given scalar to the given vector.
     * 
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     * @param y The scalar
     */
    public static void addScalar(CUstream idStream,long n, Pointer result, Pointer x, double y)
    {
        call2DGrid("addScalar", idStream, n, n, result, x, y);
    }

    /**
     * Subtract the given scalar from the given vector.
     * 
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     * @param y The scalar
     */
    public static void subScalar(CUstream idStream,long n, Pointer result, Pointer x, double y)
    {
        call2DGrid("subScalar", idStream, n, n, result, x, y);
    }

    /**
     * Multiply the given vector with the given scalar.
     * 
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     * @param y The scalar
     */
    public static void mulScalar(CUstream idStream,long n, Pointer result, Pointer x, double y)
    {
        call2DGrid("mulScalar", idStream, n, n, result, x, y);
    }

    /**
     * Divide the given vector by the given scalar.
     * 
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     * @param y The scalar
     */
    public static void divScalar(CUstream idStream,long n, Pointer result, Pointer x, double y)
    {
        call2DGrid("divScalar", idStream, n, n, result, x, y);
    }
    
    
    /**
     * Divide the given vector by the given scalar.
     * 
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     * @param y The scalar
     */
    public static void divScalarFloat(CUstream idStream,long n,Pointer result, Pointer x, float y)
    {
        call2DGrid("divScalarFloat", idStream, n, n,result, x, y);
    }
    
    /**
     * Add the given vector to the given scalar.
     * 
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The scalar
     * @param y The vector
     */
    public static void scalarAdd(CUstream idStream,long n, Pointer result, double x, Pointer y)
    {
        call2DGrid("scalarAdd", idStream, n, n, result, x, y);
    }

    /**
     * Subtract the given vector from the given scalar.
     * 
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The scalar
     * @param y The vector
     */
    public static void scalarSub(CUstream idStream,long n, Pointer result, double x, Pointer y)
    {
        call2DGrid("scalarSub", idStream, n, n, result, x, y);
    }

    /**
     * Multiply the given scalar with the given vector.
     * 
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The scalar
     * @param y The vector
     */
    public static void scalarMul(CUstream idStream,long n, Pointer result, double x, Pointer y)
    {
        call2DGrid("scalarMul", idStream, n, n, result, x, y);
    }

    /**
     * Divide the given scalar by the given vector.
     * 
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The scalar
     * @param y The vector
     */
    public static void scalarDiv(CUstream idStream,long n, Pointer result, double x, Pointer y)
    {
        call2DGrid("scalarDiv", idStream, n, n, result, x, y);
    }
    
    
    //=== Vector comparison ==================================================
    
    /**
     * Perform a '&lt;' comparison of the given vectors. 
     * The result will be set to <code>1.0f</code> where the comparison yields 
     * <code>true</code>, and to <code>0.0f</code> where the comparison yields 
     * <code>false</code>.
     *  
     * @param n The size of the vectors
     * @param result The vector that will store the result.
     * @param x The first vector
     * @param y The second vector
     */
    public static void lt(CUstream idStream,long n, Pointer result, Pointer x, Pointer y)
    {
        call2DGrid("lt", idStream, n, n, result, x, y);
    }

    /**
     * Perform a '&lt;=' comparison of the given vectors. 
     * The result will be set to <code>1.0f</code> where the comparison yields 
     * <code>true</code>, and to <code>0.0f</code> where the comparison yields 
     * <code>false</code>.
     *  
     * @param n The size of the vectors
     * @param result The vector that will store the result.
     * @param x The first vector
     * @param y The second vector
     */
    public static void lte(CUstream idStream,long n, Pointer result, Pointer x, Pointer y)
    {
        call2DGrid("lte", idStream, n, n, result, x, y);
    }

    /**
     * Perform a '==' comparison of the given vectors. 
     * The result will be set to <code>1.0f</code> where the comparison yields 
     * <code>true</code>, and to <code>0.0f</code> where the comparison yields 
     * <code>false</code>.
     *  
     * @param n The size of the vectors
     * @param result The vector that will store the result.
     * @param x The first vector
     * @param y The second vector
     */
    public static void eq(CUstream idStream,long n, Pointer result, Pointer x, Pointer y)
    {
        call2DGrid("eq", idStream, n, n, result, x, y);
    }

    /**
     * Perform a '&gt;=' comparison of the given vectors. 
     * The result will be set to <code>1.0f</code> where the comparison yields 
     * <code>true</code>, and to <code>0.0f</code> where the comparison yields 
     * <code>false</code>.
     *  
     * @param n The size of the vectors
     * @param result The vector that will store the result.
     * @param x The first vector
     * @param y The second vector
     */
    public static void gte(CUstream idStream,long n, Pointer result, Pointer x, Pointer y)
    {
        call2DGrid("gte", idStream, n, n, result, x, y);
    }

    /**
     * Perform a '&gt;' comparison of the given vectors. 
     * The result will be set to <code>1.0f</code> where the comparison yields 
     * <code>true</code>, and to <code>0.0f</code> where the comparison yields 
     * <code>false</code>.
     *  
     * @param n The size of the vectors
     * @param result The vector that will store the result.
     * @param x The first vector
     * @param y The second vector
     */
    public static void gt(CUstream idStream,long n, Pointer result, Pointer x, Pointer y)
    {
        call2DGrid("gt", idStream, n, n, result, x, y);
    }
    
    /**
     * Perform a '!=' comparison of the given vectors. 
     * The result will be set to <code>1.0f</code> where the comparison yields 
     * <code>true</code>, and to <code>0.0f</code> where the comparison yields 
     * <code>false</code>.
     *  
     * @param n The size of the vectors
     * @param result The vector that will store the result.
     * @param x The first vector
     * @param y The second vector
     */
    public static void ne(CUstream idStream,long n, Pointer result, Pointer x, Pointer y)
    {
        call2DGrid("ne", idStream, n, n, result, x, y);
    }
    
    
    //=== Vector-and-scalar comparison =======================================

    /**
     * Perform a '&lt;' comparison of the given vector and scalar. 
     * The result will be set to <code>1.0f</code> where the comparison yields 
     * <code>true</code>, and to <code>0.0f</code> where the comparison yields 
     * <code>false</code>.
     *  
     * @param n The size of the vectors
     * @param result The vector that will store the result.
     * @param x The first vector
     * @param y The scalar
     */
    public static void ltScalar(CUstream idStream,long n, Pointer result, Pointer x, double y)
    {
        call2DGrid("ltScalar", idStream, n, n, result, x, y);
    }

    /**
     * Perform a '&lt;=' comparison of the given vector and scalar. 
     * The result will be set to <code>1.0f</code> where the comparison yields 
     * <code>true</code>, and to <code>0.0f</code> where the comparison yields 
     * <code>false</code>.
     *  
     * @param n The size of the vectors
     * @param result The vector that will store the result.
     * @param x The first vector
     * @param y The scalar
     */
    public static void lteScalar(CUstream idStream,long n, Pointer result, Pointer x, double y)
    {
        call2DGrid("lteScalar", idStream, n, n, result, x, y);
    }

    /**
     * Perform a '==' comparison of the given vector and scalar. 
     * The result will be set to <code>1.0f</code> where the comparison yields 
     * <code>true</code>, and to <code>0.0f</code> where the comparison yields 
     * <code>false</code>.
     *  
     * @param n The size of the vectors
     * @param result The vector that will store the result.
     * @param x The first vector
     * @param y The scalar
     */
    public static void eqScalar(CUstream idStream,long n, Pointer result, Pointer x, double y)
    {
        call2DGrid("eqScalar", idStream, n, n, result, x, y);
    }

    /**
     * Perform a '&gt;=' comparison of the given vector and scalar. 
     * The result will be set to <code>1.0f</code> where the comparison yields 
     * <code>true</code>, and to <code>0.0f</code> where the comparison yields 
     * <code>false</code>.
     *  
     * @param n The size of the vectors
     * @param result The vector that will store the result.
     * @param x The first vector
     * @param y The scalar
     */
    public static void gteScalar(CUstream idStream,long n, Pointer result, Pointer x, double y)
    {
        call2DGrid("gteScalar", idStream, n, n, result, x, y);
    }

    /**
     * Perform a '&gt;' comparison of the given vector and scalar. 
     * The result will be set to <code>1.0f</code> where the comparison yields 
     * <code>true</code>, and to <code>0.0f</code> where the comparison yields 
     * <code>false</code>.
     *  
     * @param n The size of the vectors
     * @param result The vector that will store the result.
     * @param x The first vector
     * @param y The scalar
     */
    public static void gtScalar(CUstream idStream,long n, Pointer result, Pointer x, double y)
    {
        call2DGrid("gtScalar", idStream, n, n, result, x, y);
    }
    
    /**
     * Perform a '!=' comparison of the given vector and scalar. 
     * The result will be set to <code>1.0f</code> where the comparison yields 
     * <code>true</code>, and to <code>0.0f</code> where the comparison yields 
     * <code>false</code>.
     *  
     * @param n The size of the vectors
     * @param result The vector that will store the result.
     * @param x The first vector
     * @param y The scalar
     */
    public static void neScalar(CUstream idStream,long n, Pointer result, Pointer x, double y)
    {
        call2DGrid("neScalar", idStream, n, n, result, x, y);
    }

    
    //=== Vector math (one argument) =========================================
    

    /**
     * Calculate the arc cosine of the input argument.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void acos(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("acos", idStream, n, n, result, x);
    }

    /**
     * Calculate the nonnegative arc hyperbolic cosine of the input argument.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void acosh(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("acosh", idStream, n, n, result, x);
    }

    /**
     * Calculate the arc sine of the input argument.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void asin(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("asin", idStream, n, n, result, x);
    }

    /**
     * Calculate the arc hyperbolic sine of the input argument.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void asinh(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("asinh", idStream, n, n, result, x);
    }

    /**
     * Calculate the arc tangent of the input argument.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void atan(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("atan", idStream, n, n, result, x);
    }

    /**
     * Calculate the arc hyperbolic tangent of the input argument.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void atanh(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("atanh", idStream, n, n, result, x);
    }

    /**
     * Calculate the cube root of the input argument.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void cbrt(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("cbrt", idStream, n, n, result, x);
    }

    /**
     * Calculate ceiling of the input argument.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void ceil(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("ceil", idStream, n, n, result, x);
    }

    /**
     * Calculate the cosine of the input argument.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void cos(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("cos", idStream, n, n, result, x);
    }

    /**
     * Calculate the hyperbolic cosine of the input argument.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void cosh(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("cosh", idStream, n, n, result, x);
    }

    /**
     * Calculate the cosine of the input argument times pi
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void cospi(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("cospi", idStream, n, n, result, x);
    }

    /**
     * Calculate the complementary error function of the input argument.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void erfc(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("erfc", idStream, n, n, result, x);
    }

    /**
     * Calculate the inverse complementary error function of the input argument.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void erfcinv(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("erfcinv", idStream, n, n, result, x);
    }

    /**
     * Calculate the scaled complementary error function of the input argument.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void erfcx(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("erfcx", idStream, n, n, result, x);
    }

    /**
     * Calculate the error function of the input argument.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void erf(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("erf", idStream, n, n, result, x);
    }

    /**
     * Calculate the inverse error function of the input argument.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void erfinv(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("erfinv", idStream, n, n, result, x);
    }

    /**
     * Calculate the base 10 exponential of the input argument.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void exp10(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("exp10", idStream, n, n, result, x);
    }

    /**
     * Calculate the base 2 exponential of the input argument.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void exp2(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("exp2", idStream, n, n, result, x);
    }

    /**
     * Calculate the base e exponential of the input argument.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void exp(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("exp", idStream, n, n, result, x);
    }

    /**
     * Calculate the base e exponential of the input argument, minus 1.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void expm1(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("expm1", idStream, n, n, result, x);
    }

    /**
     * Calculate the absolute value of its argument.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void fabs(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("fabs", idStream, n, n, result, x);
    }

    /**
     * Calculate the largest integer less than or equal to x.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void floor(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("floor", idStream, n, n, result, x);
    }

    /**
     * Calculate the value of the Bessel function of the first kind of 
     * order 0 for the input argument.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void j0(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("j0", idStream, n, n, result, x);
    }

    /**
     * Calculate the value of the Bessel function of the first kind of 
     * order 1 for the input argument.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void j1(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("j1", idStream, n, n, result, x);
    }

    /**
     * Calculate the natural logarithm of the absolute value of the gamma 
     * function of the input argument.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void lgamma(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("lgamma", idStream, n, n, result, x);
    }

    /**
     * Calculate the base 10 logarithm of the input argument.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void log10(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("log10", idStream, n, n, result, x);
    }

    /**
     * Calculate the value of l o g e ( 1 + x ) .
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void log1p(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("log1p", idStream, n, n, result, x);
    }

    /**
     * Calculate the base 2 logarithm of the input argument.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void log2(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("log2", idStream, n, n, result, x);
    }

    /**
     * Calculate the floating point representation of the exponent of the 
     * input argument.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void logb(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("logb", idStream, n, n, result, x);
    }

    /**
     * Calculate the natural logarithm of the input argument.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void log(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("log", idStream, n, n, result, x);
    }

    /**
     * Calculate the standard normal cumulative distribution function.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void normcdf(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("normcdf", idStream, n, n, result, x);
    }

    /**
     * Calculate the inverse of the standard normal cumulative distribution 
     * function.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void normcdfinv(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("normcdfinv", idStream, n, n, result, x);
    }

    /**
     * Calculate reciprocal cube root function.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void rcbrt(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("rcbrt", idStream, n, n, result, x);
    }

    /**
     * Round input to nearest integer value in floating-point.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void rint(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("rint", idStream, n, n, result, x);
    }

    /**
     * Round to nearest integer value in floating-point.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void round(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("round", idStream, n, n, result, x);
    }

    /**
     * Calculate the reciprocal of the square root of the input argument.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void rsqrt(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("rsqrt", idStream, n, n, result, x);
    }

    /**
     * Calculate the sine of the input argument.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void sin(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("sin", idStream, n, n, result, x);
    }

    /**
     * Calculate the hyperbolic sine of the input argument.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void sinh(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("sinh", idStream, n, n, result, x);
    }

    /**
     * Calculate the sine of the input argument times pi
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void sinpi(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("sinpi", idStream, n, n, result, x);
    }

    /**
     * Calculate the square root of the input argument.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void sqrt(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("sqrt", idStream, n, n, result, x);
    }

    /**
     * Calculate the tangent of the input argument.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void tan(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("tan", idStream, n, n, result, x);
    }

    /**
     * Calculate the hyperbolic tangent of the input argument.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void tanh(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("tanh", idStream, n, n, result, x);
    }

    /**
     * Calculate the gamma function of the input argument.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void tgamma(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("tgamma", idStream, n, n, result, x);
    }

    /**
     * Truncate input argument to the integral part.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void trunc(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("trunc", idStream, n, n, result, x);
    }

    /**
     * Calculate the value of the Bessel function of the second kind of 
     * order 0 for the input argument.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void y0(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("y0", idStream, n, n, result, x);
    }

    /**
     * Calculate the value of the Bessel function of the second kind of 
     * order 1 for the input argument.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     */
    public static void y1(CUstream idStream,long n, Pointer result, Pointer x)
    {
        call2DGrid("y1", idStream, n, n, result, x);
    }

    //=== Vector math (two arguments) ========================================
    


    /**
     * Create value with given magnitude, copying sign of second value.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     * @param y The second vector 
     */
    public static void copysign(CUstream idStream,long n, Pointer result, Pointer x, Pointer y)
    {
        call2DGrid("copysign", idStream, n, n, result, x, y);
    }

    /**
     * Compute the positive difference between x and y.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     * @param y The second vector 
     */
    public static void fdim(CUstream idStream,long n, Pointer result, Pointer x, Pointer y)
    {
        call2DGrid("fdim", idStream, n, n, result, x, y);
    }

    /**
     * Divide two floating point values.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     * @param y The second vector 
     */
    public static void fdivide(CUstream idStream,long n, Pointer result, Pointer x, Pointer y)
    {
        call2DGrid("fdivide", idStream, n, n, result, x, y);
    }

    /**
     * Determine the maximum numeric value of the arguments.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     * @param y The second vector 
     */
    public static void fmax(CUstream idStream,long n, Pointer result, Pointer x, Pointer y)
    {
        call2DGrid("fmax", idStream, n, n, result, x, y);
    }

    /**
     * Determine the minimum numeric value of the arguments.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     * @param y The second vector 
     */
    public static void fmin(CUstream idStream,long n, Pointer result, Pointer x, Pointer y)
    {
        call2DGrid("fmin", idStream, n, n, result, x, y);
    }

    /**
     * Calculate the floating-point remainder of x / y.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     * @param y The second vector 
     */
    public static void fmod(CUstream idStream,long n, Pointer result, Pointer x, Pointer y)
    {
        call2DGrid("fmod", idStream, n, n, result, x, y);
    }

    /**
     * Calculate the square root of the sum of squares of two arguments.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     * @param y The second vector 
     */
    public static void hypot(CUstream idStream,long n, Pointer result, Pointer x, Pointer y)
    {
        call2DGrid("hypot", idStream, n, n, result, x, y);
    }

    /**
     * Return next representable single-precision floating-point value 
     * after argument.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     * @param y The second vector 
     */
    public static void nextafter(CUstream idStream,long n, Pointer result, Pointer x, Pointer y)
    {
        call2DGrid("nextafter", idStream, n, n, result, x, y);
    }

    /**
     * Calculate the value of first argument to the power of second argument.
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     * @param y The second vector 
     */
    public static void pow(CUstream idStream,long n, Pointer result, Pointer x, Pointer y)
    {
        call2DGrid("pow", idStream, n, n, result, x, y);
    }

    /**
     * Compute single-precision floating-point remainder.+u*n
     *
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     * @param y The second vector 
     */
    public static void remainder(CUstream idStream,long n, Pointer result, Pointer x, Pointer y)
    {
        call2DGrid("remainder", idStream, n, n, result, x, y);
    }

    
    /**
     * test mul the given vectors.
     * 
     * @param n The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     * @param y The second vector
     */
    public static void testkernel(CUstream idStream,long n, Pointer result, Pointer x, Pointer y)
    {
        call2DGrid("testkernel", idStream, n, n, result, x, y);
    }
    
    
    /**
     * specific function PSF modeling.
     * 
     * @param n The size of the vectors
     * @param realOutput the result real
     * @param imagOutput the result imag
     * @param kx kx
     * @param ky ky
     * @param kz kz
     * @param pupil pupil
     * @param phase phase
     * @param dx x position
     * @param dy y position
     * @param dz z position
     */
    public static void computePSF_phase(CUstream idStream,long n, Pointer realOutput, Pointer imagOutput, Pointer kx, Pointer ky, Pointer kz, Pointer pupil, Pointer phase, double dx, double dy, double dz)
    {
        call2DGrid("computePSF_phase", idStream, n, n, realOutput, imagOutput,  kx, ky, kz, pupil, phase, dx, dy, dz);
    }
    
    
    
    /**
     * specific function PSF modeling. result=(imag/divide)*(imag/divide)+(real/divide)*(real/divide)
     * 
     * @param n The size of the vectors
     * @param result the result real
     * @param real real part
     * @param imag imag part
     * @param divide 
     */
    public static void computePSF_signal(CUstream idStream,long n, Pointer result, Pointer real, Pointer imag,double divide)
    {
        call2DGrid("computePSF_signal", idStream, n, n, result, real, imag,divide);
    }
    
    
    /**
     * specific function PSF modeling. result=sqrt( (imag/divide)*(imag/divide)+(real/divide)*(real/divide) )
     * 
     * @param n The size of the vectors
     * @param result the result real
     * @param real real part
     * @param imag imag part
     * @param divide 
     */
    public static void computePSF_signalsqrt(CUstream idStream,long n, Pointer result, Pointer real, Pointer imag,double divide)
    {
        call2DGrid("computePSF_signalsqrt", idStream, n, n, result, real, imag,divide);
    }
    
    
    /**
     * specific function PSF modeling. result=sqrt( (imag/divide)*(imag/divide)+(real/divide)*(real/divide) )
     * 
     * @param n The size of the vectors
     * @param result the result real
     * @param real real part
     * @param imag imag part
     * @param divide 
     */
    public static void computePSF_signalNsqrt(CUstream idStream,long n, Pointer result, Pointer fft,double divide, Pointer sparseIndexEvenShift2D, Pointer sparseIndexOddShift2D)
    {
        call2DGrid("computePSF_signalNsqrt", idStream, n, n, result, fft,divide, sparseIndexEvenShift2D, sparseIndexOddShift2D);
    }
    
    
    
    
    /**
     * specific function PSF modeling. result=sqrt( (imag/divide)*(imag/divide)+(real/divide)*(real/divide) )
     * 
     * @param n The size of the vectors
     * @param result the result real
     * @param real real part
     * @param imag imag part
     * @param divide 
     */
    public static void computePSF_signalNsqrtMany(CUstream idStream,long n, int sizeSubImage,Pointer result, Pointer fft,double divide, Pointer sparseIndexEvenShift2D, Pointer sparseIndexOddShift2D)
    {
        call2DGrid("computePSF_signalNsqrtMany", idStream, n, n, sizeSubImage,result, fft,divide, sparseIndexEvenShift2D, sparseIndexOddShift2D);
    }
    
    
    public static void computePSF_signalNsqrtMany_f(CUstream idStream,long n, int sizeSubImage,Pointer result, Pointer fft,float divide, Pointer sparseIndexEvenShift2D, Pointer sparseIndexOddShift2D)
    {
        call2DGrid("computePSF_signalNsqrtMany_f", idStream, n, n, sizeSubImage,result, fft,divide, sparseIndexEvenShift2D, sparseIndexOddShift2D);
    }
    
    public static void computePSF_signalNsqrtMany_fcrop(CUstream idStream,long n, int sizeSubImage, int sizeSubImageFull,Pointer result, Pointer fft,float divide, Pointer sparseIndexEvenShift2D, Pointer sparseIndexOddShift2D)
    {
        call2DGrid("computePSF_signalNsqrtMany_fcrop", idStream, n, n, sizeSubImage,sizeSubImageFull,result, fft,divide, sparseIndexEvenShift2D, sparseIndexOddShift2D);
    }
    
    
    /**
     * specific function PSF modeling. result=sqrt( (imag/divide)*(imag/divide)+(real/divide)*(real/divide) )
     * 
     * @param n The size of the vectors
     * @param result the result real
     * @param real real part
     * @param imag imag part
     * @param divide 
     */
//    public static void computePSF_signalNsqrtNormalized(CUstream idStream,long n, Pointer result, Pointer fft,double divide, Pointer sparseIndexEvenShift2D, Pointer sparseIndexOddShift2D,Pointer tmp)
//    {
//        callWithSharedMemory("computePSF_signalNsqrtNormalized", idStream, n, n, result, fft,divide, sparseIndexEvenShift2D, sparseIndexOddShift2D,tmp);
//    }
    
    
    
    /**
     * specific function likelihood. result=M-I*log(M)
     * 
     * @param n The size of the vectors
     * @param result the result real
     * @param real real part
     * @param imag imag part
     * @param divide 
     */
    public static void computePoissonLikelihood(CUstream idStream,long n, Pointer result, Pointer image, Pointer model)
    {
        call2DGrid("computePoissonLikelihood", idStream, n, n, result, image, model);
    }
    
    

    
    /**
     * specific function likelihood. result=M-I*log(M)
     * 
     * @param n The size of the vectors
     * @param result the result real
     * @param real real part
     * @param imag imag part
     * @param divide 
     */
    public static void computePSF_phaseN (CUstream idStream,long n,Pointer kx, Pointer ky, Pointer kz, Pointer pupil, Pointer phase, double dx, double dy, double dz,Pointer device_sparseIndexEvenDisk, Pointer device_sparseIndexOddDisk, Pointer device_fft)
    {
        call2DGrid("computePSF_phaseN", idStream, n, n,kx, ky, kz, pupil, phase, dx, dy, dz,device_sparseIndexEvenDisk,  device_sparseIndexOddDisk,  device_fft); 
    }
    
    
    public static void computePSF_phaseNwithOil (CUstream idStream,long n,Pointer kx, Pointer ky, Pointer kz,Pointer kz_is_imag, Pointer kz_oil,Pointer kz_oil_is_imag, Pointer pupil, Pointer phase, double dx, double dy, double dz, double dz_oil,Pointer device_sparseIndexEvenDisk, Pointer device_sparseIndexOddDisk, Pointer device_fft)
    {
        call2DGrid("computePSF_phaseNwithOil", idStream, n, n,kx, ky, kz,kz_is_imag,kz_oil,kz_oil_is_imag, pupil, phase, dx, dy, dz,dz_oil,device_sparseIndexEvenDisk,  device_sparseIndexOddDisk,  device_fft); 
    }
    
    
    
    public static void computePSF_phaseNMany (CUstream idStream,long n, int sizePart, int sizeTot, Pointer kx, Pointer ky, Pointer kz, Pointer pupil, Pointer phase, Pointer position,Pointer device_sparseIndexEvenDisk, Pointer device_sparseIndexOddDisk, Pointer device_fft,int numberPSF)
                                 
    {
        call2DGrid("computePSF_phaseNMany", idStream, n, n, sizePart,sizeTot, kx, ky, kz, pupil, phase, position,device_sparseIndexEvenDisk,  device_sparseIndexOddDisk,  device_fft,numberPSF); 
    }
    
    
    
    public static void computePSF_phaseNMany_f (CUstream idStream,long n, int sizePart, int sizeTot, Pointer kx, Pointer ky, Pointer kz, Pointer pupil, Pointer phase, Pointer position,Pointer device_sparseIndexEvenDisk, Pointer device_sparseIndexOddDisk, Pointer device_fft,int numberPSF)
                                 
    {
        call2DGrid("computePSF_phaseNMany_f", idStream, n, n, sizePart,sizeTot, kx, ky, kz, pupil, phase, position,device_sparseIndexEvenDisk,  device_sparseIndexOddDisk,  device_fft,numberPSF); 
    }
    
    public static void computePSF_phaseNManywithOil_f (CUstream idStream,long n, int sizePart, int sizeTot, Pointer kx, Pointer ky, Pointer kz,Pointer kz_is_imag, Pointer kz_oil,Pointer kz_oil_is_imag, Pointer pupil, Pointer phase, Pointer position,Pointer device_sparseIndexEvenDisk, Pointer device_sparseIndexOddDisk, Pointer device_fft,int numberPSF)
                                 
    {
        call2DGrid("computePSF_phaseNManywithOil_f", idStream, n, n, sizePart,sizeTot, kx, ky, kz,kz_is_imag,kz_oil,kz_oil_is_imag, pupil, phase, position,device_sparseIndexEvenDisk,  device_sparseIndexOddDisk,  device_fft,numberPSF); 
    }
    
    
    public static void thetest (CUstream idStream,long n, Pointer kz_is_imag)
                                 
    {
        call2DGrid("thetest", idStream, n, n,kz_is_imag); 
    }
    
    
    
    /** 
     * specific function PSF modeling. result=(imag/divide)*(imag/divide)+(real/divide)*(real/divide)
     * 
     * @param n The size of the vectors
     * @param result the result real
     * @param real real part
     * @param imag imag part
     * @param divide 
     */
    public static void computePSF_signalN(CUstream idStream,long n, Pointer result, double divide,Pointer device_sparseIndexEvenShiftOutput, Pointer  device_sparseIndexOddShiftOutput, Pointer  device_fft)
    {
        call2DGrid("computePSF_signalN", idStream, n, n, result, divide,device_sparseIndexEvenShiftOutput,  device_sparseIndexOddShiftOutput,  device_fft);
    }
    
    
    
    
    
    public static void computePSF_signalN2(CUstream idStream,long n,  double divide,Pointer device_sparseIndexEvenShiftOutput, Pointer  device_sparseIndexOddShiftOutput, Pointer  device_fft,Pointer device_sparseIndexEvenShift2DOutputNext,Pointer device_sparseIndexOddShift2DOutputNext,Pointer device_psfFFT)
    {
        call2DGrid("computePSF_signalN2", idStream, n, n,  divide,device_sparseIndexEvenShiftOutput,  device_sparseIndexOddShiftOutput,  device_fft,device_sparseIndexEvenShift2DOutputNext,device_sparseIndexOddShift2DOutputNext,device_psfFFT);
    }
    
    
    
    
    public static void computePSF_signalN2Many(CUstream idStream,long n,int sizePart, int sizeTot,  double divide,Pointer device_sparseIndexEvenShiftOutput, Pointer  device_sparseIndexOddShiftOutput, Pointer  device_fft,Pointer device_sparseIndexEvenShift2DOutputNext,Pointer device_sparseIndexOddShift2DOutputNext,Pointer device_psfFFT)
    {
        call2DGrid("computePSF_signalN2Many", idStream, n, n, sizePart, sizeTot,divide,device_sparseIndexEvenShiftOutput,  device_sparseIndexOddShiftOutput,  device_fft,device_sparseIndexEvenShift2DOutputNext,device_sparseIndexOddShift2DOutputNext,device_psfFFT);
    }
    
    
     public static void computePSF_signalN2Many_f(CUstream idStream,long n,int sizePart, int sizeTot,  float divide,Pointer device_sparseIndexEvenShiftOutput, Pointer  device_sparseIndexOddShiftOutput, Pointer  device_fft,Pointer device_sparseIndexEvenShift2DOutputNext,Pointer device_sparseIndexOddShift2DOutputNext,Pointer device_psfFFT)
    {
        call2DGrid("computePSF_signalN2Many_f", idStream, n, n, sizePart, sizeTot,divide,device_sparseIndexEvenShiftOutput,  device_sparseIndexOddShiftOutput,  device_fft,device_sparseIndexEvenShift2DOutputNext,device_sparseIndexOddShift2DOutputNext,device_psfFFT);
    }
    
    
    
    
//    public static void sum(CUstream idStream,long n, Pointer input,Pointer sum)
//    {
//        callWithSharedMemory("sum", idStream, n, n,input, sum);
//    }
    
    /**
     * Multiply the given vectors.
     * 
     * @param n The size of the vectors
     * @param sizeKernel The size of the vectors
     * @param result The vector that will store the result
     * @param x The first vector
     * @param kernel The second vector
     */
    public static void mulMany(CUstream idStream,long n,int sizeKernel, Pointer result, Pointer x, Pointer kernel)
    {
        call2DGrid("mulMany", idStream, n, n,sizeKernel, result, x, kernel);
    }
    
    
    
    
    
    public static void mulMany_f(CUstream idStream,long n,int sizeKernel, Pointer result, Pointer x, Pointer kernel)
    {
        call2DGrid("mulMany_f", idStream, n, n,sizeKernel, result, x, kernel);
    }
    
    
    public static void computeModelMany1(CUstream idStream,long n,int sizeImage, Pointer result, Pointer x, Pointer amplitude,double background)
    {
        call2DGrid("computeModelMany1", idStream, n, n,sizeImage, result, x, amplitude,background);
    }
    
    public static void computeModelMany2(CUstream idStream,long n,int sizeImage, Pointer result, Pointer x, Pointer amplitude,Pointer background)
    {
        call2DGrid("computeModelMany2", idStream, n, n,sizeImage, result, x, amplitude,background);
    }
    
    public static void computeModelMany3(CUstream idStream,long n,int sizeImage, Pointer result, Pointer x, Pointer amplitude,Pointer background)
    {
        call2DGrid("computeModelMany3", idStream, n, n,sizeImage, result, x, amplitude,background);
    }
    
    
    public static void computeModelMany1_scmos(CUstream idStream,long n,int sizeImage, Pointer result, Pointer x, Pointer amplitude,double background,Pointer scmos)
    {
        call2DGrid("computeModelMany1_scmos", idStream, n, n,sizeImage, result, x, amplitude,background,scmos);
    }
    
    public static void computeModelMany2_scmos(CUstream idStream,long n,int sizeImage, Pointer result, Pointer x, Pointer amplitude,Pointer background,Pointer scmos)
    {
        call2DGrid("computeModelMany2_scmos", idStream, n, n,sizeImage, result, x, amplitude,background,scmos);
    }
    
    public static void computeModelMany3_scmos(CUstream idStream,long n,int sizeImage, Pointer result, Pointer x, Pointer amplitude,Pointer background,Pointer scmos)
    {
        call2DGrid("computeModelMany3_scmos", idStream, n, n,sizeImage, result, x, amplitude,background,scmos);
    }
    
    public static void float2double(CUstream idStream,long n, Pointer output, Pointer input)
    {
        call2DGrid("float2double", idStream, n, n,output, input);
    }
    
    public static void double2float(CUstream idStream,long n, Pointer output, Pointer input)
    {
        call2DGrid("double2float", idStream, n, n,output, input);
    }
    
    
    
    
    public static void divScalarMany(CUstream idStream,long n,int sizeSubImage, Pointer result, Pointer x, Pointer div)
    {
        call2DGrid("divScalarMany", idStream, n, n,sizeSubImage, result, x, div);
    }
     
    
    public static void addPhotonsAndBackgroundMany(CUstream idStream,long n,int sizeSubImage, Pointer output, Pointer input, Pointer photonAndBackground)
    {
        call2DGrid("addPhotonsAndBackgroundMany", idStream, n, n,sizeSubImage, output, input, photonAndBackground);
    }
    
    public static void addPhotonsAndBackgroundMany_scmos(CUstream idStream,long n,int sizeSubImage, Pointer output, Pointer input, Pointer photonAndBackground, Pointer scmos)
    {
        call2DGrid("addPhotonsAndBackgroundMany_scmos", idStream, n, n,sizeSubImage, output, input, photonAndBackground,scmos);
    }
    
    public static void divScalarMany_f(CUstream idStream,long n,int sizeSubImage, Pointer result,Pointer resultDouble, Pointer x, Pointer div)
    {
        call2DGrid("divScalarMany_f", idStream, n, n,sizeSubImage, result,resultDouble, x, div);
    }
     
    
    public static void addPhotonsAndBackgroundMany_f(CUstream idStream,long n,int sizeSubImage, Pointer output, Pointer input, Pointer photonAndBackground)
    {
        call2DGrid("addPhotonsAndBackgroundMany_f", idStream, n, n,sizeSubImage, output, input, photonAndBackground);
    }
    
    
    //this function use FLOAT pointers
    public static void complexeConjugateKernel(CUstream idStream,long n, int sizeInput, Pointer output, Pointer input, Pointer inputKernel)
    {
        call2DGrid("complexeConjugateKernel", idStream, n, n, sizeInput, output, input, inputKernel);
    }
    
    //this function use FLOAT pointers
    public static void makeResultCorrelation(CUstream idStream,long n, int sizeInput,int sizeFullPadded, Pointer output, Pointer input, Pointer sparse)
    {
        call2DGrid("makeResultCorrelation", idStream, n, n, sizeInput,sizeFullPadded, output, input, sparse);
    }
    
    
    //this function use FLOAT pointers
    public static void divCorrelation(CUstream idStream,long n, Pointer correlImage,int sizeFullImage, Pointer varImage,Pointer varPSF)
    {
        call2DGrid("divCorrelation", idStream, n, n, correlImage, sizeFullImage, varImage, varPSF);
    }
    
    
    //this function use FLOAT pointers
    public static void initIndex(CUstream idStream,long n, Pointer index)
    {
        call2DGrid("initIndex", idStream, n, n, index);
    }
    
    
    //this function use FLOAT pointers
    public static void computeLocalMaxima(CUstream idStream,long n, Pointer locMax,int sizeImage,Pointer value,int sizeFilt,int edgeSizeReject)
    {
        call2DGrid("computeLocalMaxima", idStream, n, n, locMax,sizeImage,value,sizeFilt,edgeSizeReject);
    }
    
    
    
    //this function use FLOAT pointers
    public static void eraseNonLocalMaxima(CUstream idStream,long n, Pointer value, Pointer locMax)
    {
        call2DGrid("eraseNonLocalMaxima", idStream, n, n,value, locMax);
    }
    
    
    
    
    //this function use FLOAT pointers
    public static void computeCRLB(CUstream idStream,long n, int sizeMatrix, Pointer output, Pointer input, double h)
    {
        call2DGrid("computeCRLB", idStream, n, n, sizeMatrix,output, input, h);
    }
    
    
    
    //this function use FLOAT pointers
    public static void sortRows(CUstream idStream,Pointer correl,Pointer index,long size)
    {
        
        call2DGrid("sortRows", idStream, size, size, correl,index,size);
    }
    
    
    
    
}
