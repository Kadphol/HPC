#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>

int main(int argc, char *argv[])
{
    int nt,rank;
    int i = 0, j = 0;
    int rowA, rowB, colA, colB;
    char buf[20];
    FILE *fileA,*fileB,*fileC;
    int nThread = atoi(argv[1]);
    omp_set_num_threads(nThread);

    clock_t start = clock();
    fileA = fopen(argv[2],"r");
    fileB = fopen(argv[3],"r");
    fscanf(fileA,"%d %d",&rowA,&colA);
    fscanf(fileB,"%d %d",&rowB,&colB);
    float* matrixA = (float *)malloc(rowA * colA * sizeof(float));
    float* matrixB = (float *)malloc(rowB * colB * sizeof(float));
    float* matrixC = (float *)malloc(rowA * colA * sizeof(float));

    while(fgets(buf, sizeof(buf), fileA) != NULL)
    {
        fscanf(fileA,"%f",&matrixA[i]);
        i++;
    }
    while(fgets(buf, sizeof(buf), fileB) != NULL)
    {
        fscanf(fileB,"%f",&matrixB[j]);
        j++;
    }

    fclose(fileA);
    fclose(fileB);

    #pragma omp parallel for 
    for(i = 0; i < rowA * colA; i++)
    {
        matrixC[i] = matrixA[i] + matrixB[i];
    }
    
    fileC = fopen(argv[4],"w+");
    fprintf(fileC,"%d %d\n",rowA,colA);
    for(i = 0; i < rowA; i++)
    {
        for(j = 0; j < colA; j++)
        {
            fprintf(fileC,"%f ",matrixC[(i*colA) + j]);
        }
        fprintf(fileC,"\n");
    }
    fclose(fileC);
    clock_t end = clock();
    double time = (double)(end-start)/CLOCKS_PER_SEC;
    printf("Time: %f sec\n",time);
    return 0;
}
