#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#define FILE_NAME_A "matAlarge.txt"
#define FILE_NAME_B "matBlarge.txt"
#define FILE_NAME_C "matClarge.txt"

float **alloc_2d_int(int rows, int cols) 
{
    float *data = (float *)calloc(rows*cols, sizeof(float));
    float **array= (float **)calloc(rows, sizeof(float*));
    for (int i=0; i<rows; i++)
    {
        array[i] = &(data[cols*i]);
    }   

    return array;
}

int main(int argc, char *argv[]) 
{
    
    int rank;
    int size;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    if(rank == 0) 
    {
        double startTime, endTime;
        int i,j;
        int rowA, rowB, colA, colB;
        char buf[1000];
        FILE *fpA, *fpB, *fpC;
        fpA = fopen(FILE_NAME_A,"r");
        fpB = fopen(FILE_NAME_B,"r");
        if(fpA && fpB)
        {   
            fscanf(fpA,"%d %d",&rowA,&colA);
            fscanf(fpB,"%d %d",&rowB,&colB);
            int partRow = (int)rowA/size;
            int remainingRow = rowA - (partRow * size);
            
            float **matrixA = alloc_2d_int(rowA, colA);
            float **matrixB = alloc_2d_int(rowB, colB);
            float **matrixC = alloc_2d_int(rowA, colA);

            for(i = 0; i < rowA; i++) 
            {
                for(j = 0; j < colA; j++) 
                {
                    matrixC[i][j] = 0;
                }
            }

            while(fgets(buf, sizeof(buf), fpA) != NULL) 
            {
                for(i = 0; i < rowA; i++) 
                {
                    for(j = 0; j < colA; j++) 
                    {
                        fscanf(fpA,"%f",&matrixA[i][j]);
                    }
                }
            } 

            while(fgets(buf, sizeof(buf), fpB) != NULL) 
            {
                for(i = 0; i < rowB; i++) 
                {
                    for(j = 0; j < colB; j++) 
                    {
                        fscanf(fpB,"%f",&matrixB[i][j]);
                    }
                }
            } 

            fclose(fpA);
            fclose(fpB);

            startTime = MPI_Wtime();

            for(i = 0;i < partRow; i++) 
            {
                for(j = 0;j < colA; j++) 
                {
                    matrixC[i][j] = matrixA[i][j] + matrixB[i][j];
                }
            }

            if(size > 1) 
            {
                for(int count = 1;count < size;count++)
                {
                    MPI_Send(&partRow,1,MPI_INT,count,1,MPI_COMM_WORLD);
                    MPI_Send(&colA,1,MPI_INT,count,2,MPI_COMM_WORLD);
                }

                if(remainingRow!=0) 
                {
                    for(i = partRow * size;i < rowA;i++) 
                    {
                        for(j = 0;j < colA;j++) 
                        {
                            matrixC[i][j] = matrixA[i][j] + matrixB[i][j];
                        }
                    }
                }

                for(int count = 1;count < size;count++) 
                {
                    MPI_Send(&(matrixA[partRow * count][0]),(partRow * colA),MPI_FLOAT,count,3,MPI_COMM_WORLD);
                    MPI_Send(&(matrixB[partRow * count][0]),(partRow * colA),MPI_FLOAT,count,4,MPI_COMM_WORLD);
                }
                
                for(int count = 1;count < size;count++) 
                {
                    MPI_Recv(&(matrixC[partRow * count][0]),(partRow * colA),MPI_FLOAT,count,5,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                }
            }
            endTime = MPI_Wtime();

            fpC = fopen(FILE_NAME_C,"w+");
            fprintf(fpC,"%d %d",rowA,colA);
            fprintf(fpC, "\n");
            for(i = 0;i < rowA;i++) 
            {
                for(j = 0;j < colA;j++) 
                {
                    fprintf(fpC, "%.1f ",matrixC[i][j]);
                }
                fprintf(fpC, "\n");
            }
            fclose(fpC);
            printf("Time : %f s\n", endTime - startTime);
        }
        else
        {
            printf("File not found!");
        }
    }
    else
    {
        int i,j;
        int partRow,col;
        MPI_Recv(&partRow,1,MPI_INT,0,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MPI_Recv(&col,1,MPI_INT,0,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        float **matrixA = alloc_2d_int(partRow, col);
        float **matrixB = alloc_2d_int(partRow, col);
        float **matrixC = alloc_2d_int(partRow, col);

        MPI_Recv(&(matrixA[0][0]),(partRow * col),MPI_FLOAT,0,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MPI_Recv(&(matrixB[0][0]),(partRow * col),MPI_FLOAT,0,4,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        for(i = 0 ; i < partRow ; i++) 
        {
            for(j = 0; j < col; j++) 
            {
                matrixC[i][j]=matrixA[i][j]+matrixB[i][j];
            }
        }

        MPI_Send(&(matrixC[0][0]),(partRow * col),MPI_FLOAT,0,5,MPI_COMM_WORLD);
    }
    
    MPI_Finalize();
    return 0;
}
