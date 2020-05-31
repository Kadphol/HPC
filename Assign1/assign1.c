#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char* argv[])
{
    int rank, size;
    int row, col;
    int i = 0, j, k;
    FILE *inputfile,*outputfile; 

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Request send_req[10],recv_req[20];

    double startTime = MPI_Wtime();
    if(size == 1)
    {
        inputfile = fopen(argv[1],"r");
        fscanf(inputfile,"%d %d",&row,&col);

        float *matrixInput = (float *)calloc(row * col , sizeof(float));
        float *matrixOutput = (float *)calloc(row * col , sizeof(float));

        while(fscanf(inputfile,"%f",&matrixInput[i++])==1);
        fclose(inputfile);
        int setpos = ((row-1)* col)-1;
        int setrow = row-1;
        int setcol = col-1;
        for(int k = atoi(argv[3]); k--;)
        {
            for(i = col + 1; i < setpos; i++)
            {
                matrixOutput[i] = (matrixInput[(i-(col+1))] + matrixInput[(i-(col))] + matrixInput[(i-(col-1))]+
                                    matrixInput[i-1] + matrixInput[i] + matrixInput[i+1]+
                                    matrixInput[(i+(col-1))] + matrixInput[(i+(col))] + matrixInput[(i+(col+1))])/9;
            }
            for(i = 1; i < setrow; i++)
            {
                for(j = 1; j < setcol; j++)
                {
                    matrixInput[(i*col) +j] = matrixOutput[(i*col) + j];
                }
            }
        }

        outputfile = fopen(argv[2],"w+");
        fprintf(outputfile,"%d %d\n",row,col);
        while(k < row*col)
        {
            if((k % col) != (col - 1))
            {
                fprintf(outputfile,"%.f ",matrixInput[k++]);
            }
            else
            {
                fprintf(outputfile,"%.f\n",matrixInput[k++]);
            }
        }
        fclose(outputfile);
        double endTime = MPI_Wtime();
        printf("Time: %f s\n",(endTime - startTime));
    }
    else if(rank == 0)
    {
        inputfile = fopen(argv[1],"r");
        fscanf(inputfile,"%d %d",&row,&col);
        MPI_Bcast(&col, 1, MPI_INT, 0, MPI_COMM_WORLD);
        i = 0;
        int srow = (row-2)/size;
        int partRow[size];
        int remainingRow = (row-2)%size;
        for(i = size;i--;)
        {
            partRow[i] = srow;
            if(remainingRow != 0)
            {
                partRow[i]++;
                remainingRow--;
            }
        }

        float* matrixInput = (float *)malloc(row * col * sizeof(float));
        float* matrixOutput = (float *)malloc(row * col * sizeof(float));
        i = 0;
        while(fscanf(inputfile,"%f",&matrixInput[i++])==1);
        fclose(inputfile);

        int rowNum = partRow[0];
        for(int count = 1; count < size; count++)
        {
            MPI_Isend(&partRow[count], 1, MPI_INT, count, 3, MPI_COMM_WORLD,&send_req[0]);
            MPI_Isend(&matrixInput[(rowNum) * col], (partRow[count] + 2) * col, MPI_FLOAT, count, 4, MPI_COMM_WORLD, &send_req[count]);
            MPI_Wait(&send_req[0], MPI_STATUS_IGNORE);
            MPI_Wait(&send_req[count], MPI_STATUS_IGNORE);
            rowNum += partRow[count];
        }
        int setpos = ((partRow[0]+1)* col)-1;
        int setrow = partRow[0]+1;
        int setcol = col-1;
        for(int k = atoi(argv[3]); k--;)
        {
            for(i = col + 1; i < setpos; i++)
            {
                matrixOutput[i] = (matrixInput[(i-(col+1))] + matrixInput[(i-(col))] + matrixInput[(i-(col-1))]+
                                matrixInput[i-1] + matrixInput[i] + matrixInput[i+1]+
                                matrixInput[(i+(col-1))] + matrixInput[(i+(col))] + matrixInput[(i+(col+1))])/9;
            }
            for(i = 1; i < setrow; i++)
            {
                for(j = 1; j < setcol; j++)
                {
                    matrixInput[(i*col) +j] = matrixOutput[(i*col) + j];
                }
            }
            MPI_Isend(&(matrixInput[partRow[0] * col]), col, MPI_FLOAT, 1, 5, MPI_COMM_WORLD, &send_req[3]);
            MPI_Irecv(&(matrixInput[(partRow[0] + 1) * col]), col, MPI_FLOAT, 1, 6, MPI_COMM_WORLD, &recv_req[4]);
            MPI_Wait(&send_req[3],MPI_STATUS_IGNORE);
            MPI_Wait(&recv_req[4],MPI_STATUS_IGNORE);
        }

        int position = partRow[0] + 1;
        for(int count = 1; count < size; count++)
        {
            MPI_Irecv(&(matrixInput[position * col]), partRow[count] * col, MPI_FLOAT, count, 7, MPI_COMM_WORLD, &recv_req[count]);
            MPI_Wait(&recv_req[count],MPI_STATUS_IGNORE);
            position += partRow[count];
        }
        outputfile = fopen(argv[2],"w+");
        fprintf(outputfile,"%d %d\n",row,col);
        while(k < row*col)
        {
            if((k % col) != (col - 1))
            {
                fprintf(outputfile,"%.f ",matrixInput[k++]);
            }
            else
            {
                fprintf(outputfile,"%.f\n",matrixInput[k++]);
            }
        }
        fclose(outputfile);
        double endTime = MPI_Wtime();
        printf("Time: %f s\n",(endTime - startTime));
    }
    else        
    {
        int col,partRow;
        MPI_Irecv(&partRow, 1, MPI_INT, 0, 3, MPI_COMM_WORLD, &recv_req[rank]);
        MPI_Bcast(&col, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Wait(&recv_req[rank],MPI_STATUS_IGNORE);
        
        float *matrixInputBuf = (float *) malloc ((partRow + 2) * col * sizeof(float));
        float *matrixOutput = (float *) malloc ((partRow + 2)  * col * sizeof(float));
        MPI_Irecv(&matrixInputBuf[0], (partRow + 2) * col, MPI_FLOAT, 0, 4, MPI_COMM_WORLD, &recv_req[rank]);
        MPI_Wait(&recv_req[rank],MPI_STATUS_IGNORE);

        int setpos = ((partRow+1)* col)-1;
        int setrow = partRow+1;
        int setcol = col-1;
        for(k = atoi(argv[3]); k--;)
        {
            for(i = col + 1; i < setpos; i++)
            {
                matrixOutput[i] = (matrixInputBuf[(i-(col+1))] + matrixInputBuf[(i-(col))] + matrixInputBuf[(i-(col-1))]+
                                matrixInputBuf[i-1] + matrixInputBuf[i] + matrixInputBuf[i+1]+
                                matrixInputBuf[(i+(col-1))] + matrixInputBuf[(i+(col))] + matrixInputBuf[(i+(col+1))])/9;
            }
            for(i = 0; i < setrow; i++)
            {
                for(j = 1; j < setcol; j++)
                {
                    matrixInputBuf[(i*col) +j] = matrixOutput[(i*col) + j];
                }
            }
            if(rank == size-1)
            {
                MPI_Isend(&(matrixInputBuf[col]), col, MPI_FLOAT, rank - 1, 6, MPI_COMM_WORLD, &send_req[rank]);
                MPI_Irecv(&(matrixInputBuf[0]), col, MPI_FLOAT, rank - 1, 5, MPI_COMM_WORLD, &recv_req[rank]);
                MPI_Wait(&send_req[rank],MPI_STATUS_IGNORE);
                MPI_Wait(&recv_req[rank],MPI_STATUS_IGNORE);
            }
            else
            {
                MPI_Isend(&(matrixInputBuf[col]), col, MPI_FLOAT, rank - 1, 6, MPI_COMM_WORLD, &send_req[rank]);
                MPI_Isend(&(matrixInputBuf[partRow * col]), col, MPI_FLOAT, rank+1, 5, MPI_COMM_WORLD, &send_req[rank]);
                MPI_Irecv(&(matrixInputBuf[0]), col, MPI_FLOAT, rank - 1, 5, MPI_COMM_WORLD, &recv_req[rank-1]);
                MPI_Irecv(&(matrixInputBuf[(partRow+1) * col]), col, MPI_FLOAT, rank + 1, 6, MPI_COMM_WORLD, &recv_req[rank+1]);
                MPI_Wait(&send_req[rank],MPI_STATUS_IGNORE);
                MPI_Wait(&send_req[rank],MPI_STATUS_IGNORE);
                MPI_Wait(&recv_req[rank-1],MPI_STATUS_IGNORE);
                MPI_Wait(&recv_req[rank+1],MPI_STATUS_IGNORE);
            }
        }
       
        MPI_Send(&(matrixInputBuf[col]), partRow * col, MPI_FLOAT, 0 , 7, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}
