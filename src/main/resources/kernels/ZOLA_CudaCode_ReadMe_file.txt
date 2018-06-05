Requirement: NVCC Cuda toolkit is installed, and your gcc/g++ version is 4.8.
To check it, tape in a terminal:
nvcc --version
gcc --version
g++ --version

To compile the source code, tape in a terminal the two following lines.

nvcc -gencode arch=compute_20,code=sm_20 --default-stream per-thread -ptx ZOLA_kernels.cu -o ZOLA_kernels_20.ptx

nvcc -gencode arch=compute_30,code=sm_30 --default-stream per-thread -ptx ZOLA_kernels.cu -o ZOLA_kernels_30.ptx

