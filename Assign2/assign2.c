#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>
#include <mpi.h>

void swap(double* a, double* b) 
{ 
    double t = *a; 
    *a = *b; 
    *b = t; 
} 

int partition (double a[], int low, int high) 
{ 
    double pivot = a[high];
    int i = (low - 1);
  
    for (int j = low; j <= high- 1; j++) 
    { 
        if (a[j] < pivot) 
        { 
            i++; 
            swap(&a[i], &a[j]); 
        } 
    } 
    swap(&a[i + 1], &a[high]); 
    return (i + 1); 
} 

void quicksort(double* a,int low, int high)
{
    int pivot;
    
    if(low < high)
    {
        pivot = partition(a,low,high);
        #pragma omp task
        quicksort(a, low, pivot-1);
        #pragma omp task
        quicksort(a, pivot + 1, high);
    }
}

double * merge(double * a1, int n1, double * a2, int n2)
{
    double * result = (double *)malloc((n1 + n2) * sizeof(double));
    int i = 0, j = 0, k;
    for(k = 0; k < n1+n2; k++)
    {
        if(i >= n1)
        {
            result[k] = a2[j];
            j++;
        }
        else if(j >= n2)
        {
            result[k] = a1[i];
            i++;
        }
        else if(a1[i] < a2[j])
        {
            result[k] = a1[i];
            i++;
        }
        else
        {
            result[k] = a2[j];
            j++;
        }
    }
    return result;
}

int main(int argc,char *argv[])
{
    int rank, size;
    int i = 0;
    int n;
    double start, end;
    double *array;
    double seq_time, par_time, merge_time;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    FILE *inputfile,*outputfile;
    int nThread = atoi(argv[3]);
    omp_set_num_threads(nThread);

    if(size == 1)
    {
        start = omp_get_wtime();
        inputfile = fopen(argv[1],"r");
        fscanf(inputfile,"%d",&n);

        array = (double *)malloc(n * sizeof(double));
        while(!feof(inputfile))
        {
            fscanf(inputfile,"%lf",&array[i]);
            i++;
        }
        fclose(inputfile);
        end = omp_get_wtime();
        seq_time = end - start;

        start = omp_get_wtime();

        #pragma omp parallel
        #pragma omp single
        quicksort(array, 0 , n-1);

        end = omp_get_wtime();
        par_time = end - start;

        start = omp_get_wtime();
        outputfile = fopen(argv[2],"w");
        fprintf(outputfile,"%d\n",n);
        for(i = 0;i < n;i++)
        {
            fprintf(outputfile,"%.4lf\n",array[i]);
        }
        fclose(outputfile);
        end = omp_get_wtime();
        seq_time += end-start;
        printf("read write time: %lf\n",seq_time);
        printf("parallel time: %lf\n",par_time);
        printf("ratio: %lf\n",(seq_time/(seq_time + par_time)));
    }
    else
    {   
        if(rank == 0)
        {
            start = omp_get_wtime();
            inputfile = fopen(argv[1],"r");
            fscanf(inputfile,"%d",&n);
            array = (double *)malloc(n * sizeof(double));
            while(!feof(inputfile))
            {
                fscanf(inputfile,"%lf",&array[i]);
                i++;
            }
            fclose(inputfile);
            end = omp_get_wtime();
            seq_time = end - start;
        }
        
        start = omp_get_wtime();
        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
        int partSize = n/size;

        double *bufArray = (double *)malloc(n * sizeof(double));
        MPI_Scatter(&(*array), partSize, MPI_DOUBLE, &(*bufArray), partSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        #pragma omp parallel
        #pragma omp single
        quicksort(bufArray, 0 , partSize-1);

        MPI_Barrier(MPI_COMM_WORLD);
        
        for(int step = 1; step < size; step *= 2) 
        {
            if(rank%(2*step) == 0)
            {
                if(rank+step < size)
                { 
                    int m;
                    MPI_Recv(&m, 1, MPI_INT, rank+step, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    double *chuck = (double *)malloc(m * sizeof(double));
                    MPI_Recv(&chuck[0], m, MPI_DOUBLE, rank+step, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    double mergestart = omp_get_wtime();
                    bufArray = merge(bufArray, partSize, chuck, m);
                    double mergeend = omp_get_wtime();
                    merge_time += mergeend - mergestart;
                    partSize = partSize + m;
                }
            }
            else
            {
                int near = rank - step;
                MPI_Send(&partSize, 1, MPI_INT, near, 0, MPI_COMM_WORLD);
                MPI_Send(&bufArray[0], partSize, MPI_DOUBLE, near, 1, MPI_COMM_WORLD);
                break;
            }
        }
        end = omp_get_wtime();
        par_time = end - start - merge_time;

        if(rank == 0)
        {
            start = omp_get_wtime();
            outputfile = fopen(argv[2],"w");
            fprintf(outputfile,"%d\n",n);
            for(i = 0;i < n;i++)
            {
                fprintf(outputfile,"%.4lf\n",bufArray[i]);
            }
            fclose(outputfile);
            end = omp_get_wtime();
            seq_time += end-start;
            seq_time += merge_time;
            printf("read write time: %lf\n",seq_time);
            printf("parallel time: %lf\n",par_time);
            printf("ratio: %lf\n",(seq_time/(seq_time + par_time)));
        }
    }

    MPI_Finalize();
    return 0;
}
