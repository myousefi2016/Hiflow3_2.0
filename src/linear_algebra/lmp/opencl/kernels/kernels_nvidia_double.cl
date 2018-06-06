#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void increment(__global double *x, const unsigned int size)
{
	int n = get_global_id(0);
	if(n<size) x[n] = x[n] + 1;
}

__kernel void opencl_setvalues(__global double *x, const unsigned int size, __global int *ind, __global double *val)
{
	int n = get_global_id(0);
	if(n<size) x[ind[n]] = val[n];
}

__kernel void opencl_getvalues(__global double *x, const unsigned int size, __global int *ind, __global double *val)
{
	int n = get_global_id(0);
	if(n<size) val[n] = x[ind[n]];
}

__kernel void opencl_setblockvalues(__global double *x, const unsigned int start_i, const unsigned int start_sub_vec, const unsigned int size, __global double *val)
{
	int n = get_global_id(0);
	if(n<size) x[n+start_i] = val[n+start_sub_vec];
}

__kernel void opencl_getblockvalues(__global double *x, const unsigned int start_i, const unsigned int end_i, __global double *val)
{
	int n = get_global_id(0);
	if(n<end_i-start_i) val[n] = x[n+start_i];
}

__kernel void opencl_zeros(__global double *x, const int size)
{
	int n = get_global_id(0);
	if(n<size) x[n] = 0.0;
}

__kernel void opencl_scaled_add(__global double *x, __global const double *y, const double a, const unsigned int size) {
	// y = x + a*y
	int n = get_global_id(0);
	if (size > n) x[n] = y[n] + a*x[n];	
}

__kernel void opencl_axpy(__global double *x, __global const double *y, const double a, const unsigned int size) {
	// y = a*x + y
	int n = get_global_id(0);
	if (size > n) x[n] += y[n] * a;	
}

__kernel void opencl_scale(__global double *x, const double a, const unsigned int size) {
	// y = a*y
	int n = get_global_id(0);
	if (size > n) x[n] *= a;	
}

__kernel void opencl_multvalues(const unsigned int size, __global const double *values, __global double *buffer) {
	int n = get_global_id(0);
	if (size > n) buffer[n] *= values[n];
}

__kernel void opencl_dot(__global const double *vec_x, __global const double *vec_y, __global double *result, const unsigned int size, __local double* sdata)
{
    // get index into global data array
    unsigned int tid = get_local_id(0);
    unsigned int blockSize = get_local_size(0);
    unsigned int i = get_group_id(0)*(get_local_size(0)*2) + get_local_id(0);

    sdata[tid] = (i<size) ? vec_x[i]*vec_y[i] : 0;

    if ((i + blockSize) < size)
      sdata[tid] += vec_x[i + blockSize]*vec_y[i + blockSize];

    barrier(CLK_LOCAL_MEM_FENCE);

    if (blockSize >= 512) {if (tid < 256) { sdata[tid] += sdata[tid + 256];} barrier(CLK_LOCAL_MEM_FENCE); }
    if (blockSize >= 256) {if (tid < 128) { sdata[tid] += sdata[tid + 128];} barrier(CLK_LOCAL_MEM_FENCE); }
    if (blockSize >= 128) {if (tid <  64) { sdata[tid] += sdata[tid +  64];} barrier(CLK_LOCAL_MEM_FENCE); }

    if (tid < 32) {
     if (blockSize >= 64) { sdata[tid] += sdata[tid +  32]; barrier(CLK_LOCAL_MEM_FENCE); }
     if (blockSize >= 32) { sdata[tid] += sdata[tid +  16]; barrier(CLK_LOCAL_MEM_FENCE); }
     if (blockSize >= 16) { sdata[tid] += sdata[tid +   8]; barrier(CLK_LOCAL_MEM_FENCE); }
     if (blockSize >=  8) { sdata[tid] += sdata[tid +   4]; barrier(CLK_LOCAL_MEM_FENCE); }
     if (blockSize >=  4) { sdata[tid] += sdata[tid +   2]; barrier(CLK_LOCAL_MEM_FENCE); }
     if (blockSize >=  2) { sdata[tid] += sdata[tid +   1]; barrier(CLK_LOCAL_MEM_FENCE); }
    }

    if ( tid==0 ) result[0] = sdata[0];

}

__kernel void opencl_nrm2(__global const double *vec_x, __global double *result, const unsigned int size, __local double* sdata)
{
    // get index into global data array
    unsigned int tid = get_local_id(0);
    unsigned int blockSize = get_local_size(0);
    unsigned int i = get_group_id(0)*(get_local_size(0)*2) + get_local_id(0);

    sdata[tid] = (i<size) ? vec_x[i]*vec_x[i] : 0;

    if ((i + blockSize) < size)
      sdata[tid] += vec_x[i + blockSize]*vec_x[i + blockSize];

    barrier(CLK_LOCAL_MEM_FENCE);

    if (blockSize >= 512) {if (tid < 256) { sdata[tid] += sdata[tid + 256];} barrier(CLK_LOCAL_MEM_FENCE); }
    if (blockSize >= 256) {if (tid < 128) { sdata[tid] += sdata[tid + 128];} barrier(CLK_LOCAL_MEM_FENCE); }
    if (blockSize >= 128) {if (tid <  64) { sdata[tid] += sdata[tid +  64];} barrier(CLK_LOCAL_MEM_FENCE); }

    if (tid < 32) {
     if (blockSize >= 64) { sdata[tid] += sdata[tid +  32]; barrier(CLK_LOCAL_MEM_FENCE); }
     if (blockSize >= 32) { sdata[tid] += sdata[tid +  16]; barrier(CLK_LOCAL_MEM_FENCE); }
     if (blockSize >= 16) { sdata[tid] += sdata[tid +   8]; barrier(CLK_LOCAL_MEM_FENCE); }
     if (blockSize >=  8) { sdata[tid] += sdata[tid +   4]; barrier(CLK_LOCAL_MEM_FENCE); }
     if (blockSize >=  4) { sdata[tid] += sdata[tid +   2]; barrier(CLK_LOCAL_MEM_FENCE); }
     if (blockSize >=  2) { sdata[tid] += sdata[tid +   1]; barrier(CLK_LOCAL_MEM_FENCE); }
    }

    if ( tid==0 ) result[0] = sqrt(sdata[0]);

}

__kernel void opencl_csr_vectormult(__global const double *val, __global const int *col, __global const int *row, __global const double *x, __global double *y, const unsigned int size) {
	int j;
	int tid = get_global_id(0);
	double tmp = 0;

  if (tid < size) {

	y[tid] = 0;
	for(j=row[tid];j<row[tid+1];j++){
		tmp +=  val[j] * x[col[j]];	
	}
	y[tid] = tmp;
  }

}

__kernel void opencl_csr_vectormultadd(__global const double *val, __global const unsigned int *col, __global const unsigned int *row, __global const double *x, __global double *y, const unsigned int size) {
	int j;
	int n = get_global_id(0);
	double tmp = 0;

	if (n < size) {

	for(j=row[n];j<row[n+1];j++){
		tmp += val[j] * x[col[j]];
	}
	y[n] += tmp;
       }
}

// Kernel from NVIDIA SDK
__kernel void DotProduct(__global const double *in_x, __global const double *in_y, __global double *out, __local double *sdata, const unsigned int n, const int blockSize, const int nIsPow2) {
    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    unsigned int tid = get_local_id(0);
    unsigned int i = get_group_id(0)*(get_local_size(0)*2) + get_local_id(0);
    unsigned int gridSize = blockSize*2*get_num_groups(0);
    sdata[tid] = 0;

    // we reduce multiple elements per thread.  The number is determined by the 
    // number of active thread blocks (via gridDim).  More blocks will result
    // in a larger gridSize and therefore fewer elements per thread
    while (i < n)
    {         
        sdata[tid] += in_x[i]*in_y[i];
        // ensure we don't read out of bounds -- this is optimized away for powerOf2 sized arrays
        if (nIsPow2 || i + blockSize < n) 
            sdata[tid] += in_x[i+blockSize]*in_y[i+blockSize];  
        i += gridSize;
    } 

    barrier(CLK_LOCAL_MEM_FENCE);

    // do reduction in shared mem
    if (blockSize >= 512) { if (tid < 256) { sdata[tid] += sdata[tid + 256]; } barrier(CLK_LOCAL_MEM_FENCE); }
    if (blockSize >= 256) { if (tid < 128) { sdata[tid] += sdata[tid + 128]; } barrier(CLK_LOCAL_MEM_FENCE); }
    if (blockSize >= 128) { if (tid <  64) { sdata[tid] += sdata[tid +  64]; } barrier(CLK_LOCAL_MEM_FENCE); }

    if (tid < 32)
    {
        if (blockSize >=  64) { sdata[tid] += sdata[tid + 32]; }
        if (blockSize >=  32) { sdata[tid] += sdata[tid + 16]; }
        if (blockSize >=  16) { sdata[tid] += sdata[tid +  8]; }
        if (blockSize >=   8) { sdata[tid] += sdata[tid +  4]; }
        if (blockSize >=   4) { sdata[tid] += sdata[tid +  2]; }
        if (blockSize >=   2) { sdata[tid] += sdata[tid +  1]; }
    }

    // write result for this block to global mem 
    if (tid == 0) out[get_group_id(0)] = sdata[0];
}

// Kernel from NVIDIA SDK
__kernel void Reduction(__global double *in, __global double *out, __local double *sdata, const unsigned int n, const int blockSize) {
    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    unsigned int tid = get_local_id(0);
    unsigned int i = get_group_id(0)*(get_local_size(0)*2) + get_local_id(0);

    sdata[tid] = (i < n) ? in[i] : 0;
    if (i + blockSize < n) 
        sdata[tid] += in[i+blockSize];  

    barrier(CLK_LOCAL_MEM_FENCE);

    // do reduction in shared mem
    if (blockSize >= 512) { if (tid < 256) { sdata[tid] += sdata[tid + 256]; } barrier(CLK_LOCAL_MEM_FENCE); }
    if (blockSize >= 256) { if (tid < 128) { sdata[tid] += sdata[tid + 128]; } barrier(CLK_LOCAL_MEM_FENCE); }
    if (blockSize >= 128) { if (tid <  64) { sdata[tid] += sdata[tid +  64]; } barrier(CLK_LOCAL_MEM_FENCE); }

    if (tid < 32)
    {
        if (blockSize >=  64) { sdata[tid] += sdata[tid + 32]; }
        if (blockSize >=  32) { sdata[tid] += sdata[tid + 16]; }
        if (blockSize >=  16) { sdata[tid] += sdata[tid +  8]; }
        if (blockSize >=   8) { sdata[tid] += sdata[tid +  4]; }
        if (blockSize >=   4) { sdata[tid] += sdata[tid +  2]; }
        if (blockSize >=   2) { sdata[tid] += sdata[tid +  1]; }
    }

    // write result for this block to global mem 
    if (tid == 0) out[get_group_id(0)] = sdata[0];
}

