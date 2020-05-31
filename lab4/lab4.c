#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>


void writeMatrix(char* filename,double *matrix,int row,int col)
{
    FILE *file;
    file = fopen(filename,"w+");
    fprintf(file,"%d %d",row,col);
    fprintf(file,"\n");
    for(int i = 0;i < row; i++)
    {
        for(int j = 0;j < col; j++)
        {
            fprintf(file,"%.10lf ",matrix[(i*col)+j]);
        }
        fprintf(file,"\n");
    }
    fclose(file);
}

void readMatrix(FILE *file, double *matrix, int row, int col, int reverse)
{
    int i, j;
    char buf[20];
    while(fgets(buf, sizeof(buf), file) != NULL)
    {
        for(i = 0;i < row; i++)
        {
            for(j = 0;j < col; j++)
            {
                if(reverse == 1)
                {
                    fscanf(file,"%lf",&matrix[(j*row)+i]);
                }
                else
                {
                    fscanf(file,"%lf",&matrix[(i*col)+j]);
                }
            }
        }
    }
    fclose(file);
}

int main(int argc,char *argv[]) 
{   
    int rank;
    int size;
    int temp;
    int i, j, k;
    double sum;
    int rowA, colA, rowB, colB;
    int partRow, remainingRow;
    double *matrixA , *matrixB, *matrixC;
    double *matrixBuffer, *matrixCBuffer;
    FILE *fileA,*fileB,*fileC;
    double start,end;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    char* filename_A = argv[1];
    char* filename_B = argv[2];
    char* filename_C = argv[3];

    if(rank == 0)
    {
        fileA = fopen(filename_A,"r");
        fileB = fopen(filename_B,"r");

        fscanf(fileA,"%d %d",&rowA,&colA);
        fscanf(fileB,"%d %d",&rowB,&colB);

        partRow = (int)rowA / size;
        remainingRow = rowA % size;

        matrixA = (double *)malloc(rowA * colA * sizeof(double));
        matrixB = (double *)malloc(rowB * colB * sizeof(double));
        matrixC = (double *)malloc(rowA * colB * sizeof(double));

        readMatrix(fileA,matrixA,rowA,colA,0);
        readMatrix(fileB,matrixB,rowB,colB,1);

        //swap rowB and colB
        temp = colB;
        colB = rowB;
        rowB = temp;
    }

    start = MPI_Wtime();
    MPI_Bcast(&rowA,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&colA,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&rowB,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&colB,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&partRow,1,MPI_INT,0,MPI_COMM_WORLD);

    if(rank != 0) 
    {
       matrixB = (double *)malloc(rowB * colB * sizeof(double));
    }
    
    matrixBuffer = (double *)malloc(partRow * colA * sizeof(double));
    matrixCBuffer = (double *)malloc(partRow * rowB * sizeof(double));
    
    MPI_Bcast(matrixB,rowB * colB,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Scatter(&(*matrixA),partRow * colA,MPI_DOUBLE,&(*matrixBuffer),partRow * colA,MPI_DOUBLE,0,MPI_COMM_WORLD);

    for (i = 0; i < partRow; i++)
    {
        for(j = 0; j < rowB; j++) 
        {   
            for (k = 0; k < colA; k++)
            {
                sum += matrixBuffer[(i*colA)+k]*matrixB[(j*colA)+k];
            }
            matrixCBuffer[(i*rowB)+j] = sum;
            sum = 0;
        }
    }

    MPI_Gather(&(*matrixCBuffer),partRow * rowB,MPI_DOUBLE,&(*matrixC),partRow * rowB,MPI_DOUBLE,0,MPI_COMM_WORLD);

    if(rank == 0)
    {
        if(remainingRow != 0)
        {
            sum = 0;
            for(i = partRow * size; i< (partRow * size) + remainingRow; i++)
            {
                for(j=0; j < rowB; j++)
                {
                    for(k = 0 ; k < colA; k++)
                    {
                        sum += (matrixA[(i*colA)+k]) * (matrixB[(j*colA)+k]);
                    }
                    (matrixC[(i*rowB)+j]) = sum;
                    sum = 0;
                }
            }
        }
        end = MPI_Wtime();
        writeMatrix(filename_C,matrixC,rowA,rowB);
        printf("Time: %f s\n",end - start);
    }
    
    MPI_Finalize();
    return 0;
}
