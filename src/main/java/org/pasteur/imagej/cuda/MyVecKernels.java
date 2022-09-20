/*
 * code with Licence MIT
 */
package org.pasteur.imagej.cuda;


import jcuda.driver.CUstream;
import jcuda.runtime.cudaStream_t;

/**
 * A set of CUDA kernels for vector operations. The kernels are CUfunction 
 * instances that are created from a CUmodule.
 */
interface MyVecKernels
{
    /**
     * Call the kernel that is identified by the given name, with the
     * given arguments. Note that the given name must not necessarily
     * be the name that the kernel has in the CUDA source code.
     * 
     * @param name The name identifying the kernel
     * @param workSize The global work size of the kernel
     * @param arguments The arguments for the kernel
     * @throws CudaException If the kernel could not be called
     */
    void call(String name, long workSize,CUstream idStream, Object ... arguments);
    void callSingle(String name, long workSize,CUstream idStream, Object ... arguments);
    void call2DGrid(String name, long workSize,CUstream idStream, Object ... arguments);
    void callWithSharedMemory(String name, long workSize,CUstream idStream, Object ... arguments);
    /**
     * Perform a shutdown, releasing all resources that have been
     * allocated by this instance.
     */
    void shutdown();
}