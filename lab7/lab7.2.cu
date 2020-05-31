#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda.h>
#include <curand_kernel.h>

const int THREADS = 128; //threads per block
const int trial_number = 1024; //trial numbers per thread
const int BLOCKS = 16; //blocks per grid
const float PI = 3.1415926535;

__global__ void pi_estimation(float *pi, curandState *states)
{
    unsigned int threadID = threadIdx.x + blockDim.x * blockIdx.x;
    int count = 0;
    float x,y,z;

    curand_init(0, threadID, 1, &states[threadID]); //(seed, sequence number, offset, curandState)

    for(int i = 0; i < trial_number; i++)
    {
        x = curand_uniform(&states[threadID]);//return sequence number of pseudorandom uniformly distributed.
        y = curand_uniform(&states[threadID]);
        z = (x*x + y*y);
        if(z <= 1.0f) //if x,y in unit circle
        {
            count += 1;
        }
    }
    pi[threadID] = 4.0f * count/(float)trial_number; //estimate PI value 4*count/number of trial
} 

int main(int argc, char *argv[])
{
    float host[THREADS * BLOCKS];
    float *device;
    curandState *deviceStates;

    cudaMalloc((void **)&device, THREADS * BLOCKS * sizeof(float));
    cudaMalloc((void **)&deviceStates, THREADS * BLOCKS * sizeof(curandState));

    pi_estimation<<<BLOCKS,THREADS>>>(device,deviceStates); //call kernel
    cudaMemcpy(host, device, THREADS * BLOCKS * sizeof(float), cudaMemcpyDeviceToHost); //copy estimated pi value from device to host
    float pi = 0.0;
    for(int i = 0; i < THREADS * BLOCKS; i++)
    {
        pi += host[i];
    }
    pi /= (THREADS * BLOCKS); //get average of each PI value from each thread.
    printf("Monte Carlo PI estimation %d times\n", THREADS * BLOCKS * trial_number);
    printf("PI estimation: %.10f\n",pi);
    printf("Error: %.10\n\n",pi-PI);
    cudaFree(device);
    cudaFree(deviceStates);
    return 0;
}
