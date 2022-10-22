#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <malloc.h>
#include <string.h>
#include "data.h"
#include "timer.h"

double **alloc2d(unsigned long int n, unsigned long int m)
{
    double *data = malloc(n * m * sizeof(double));
    double **array = malloc(n * sizeof(double *));
    for (unsigned long int i = 0; i < n; i++)
    {
        array[i] = &(data[i * m]);
    }
    return array;
}

void free2d(double **array)
{
    free(array[0]);
    free(array);
}

int main(int argc, char **argv)
{
    if (argc != 3)
    {
        printf("Usage: <num_thread> <vec_a> <vec_b>.%d\n", argc);
        exit(EXIT_FAILURE);
    }
    data_struct *d_1;
    data_struct *d_2;
    data_struct *d_3;

    struct timespec start;

    int ProcessorId, NoOfProcesors;
    unsigned long int RowCountPerProcessor;
    unsigned long int MatrixColumnCount;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &NoOfProcesors);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcessorId);

    if (ProcessorId == 0) // Master Processor
    {
        d_1 = get_data_struct(argv[1]);
        d_2 = get_data_struct(argv[2]);
        d_3 = allocate_data_struct(d_1->rows, d_2->cols);

        if (d_1->cols != d_2->cols || d_1->rows != d_2->rows)
        {
            printf("ERROR: The dimiension of Matrix is not correct.\n");
            exit(EXIT_FAILURE);
        }
        start_timer(&start);
        RowCountPerProcessor = d_1->rows / NoOfProcesors;
        MatrixColumnCount = d_2->cols;

        data_struct *Temp_d_1;
        data_struct *Temp_d_2;
        data_struct *Temp_d_3;
        Temp_d_1 = malloc(sizeof(data_struct));
        Temp_d_2 = malloc(sizeof(data_struct));
        Temp_d_3 = malloc(sizeof(data_struct));
        Temp_d_1->data_point = alloc2d(RowCountPerProcessor, MatrixColumnCount);
        Temp_d_2->data_point = alloc2d(RowCountPerProcessor, MatrixColumnCount);
        Temp_d_3->data_point = alloc2d(RowCountPerProcessor, MatrixColumnCount);

        if (NoOfProcesors > 1)
        {
            for (unsigned long int i = 1; i < (NoOfProcesors - 1); i++)
            {
                unsigned long int ProcessorRowStartndex = i * RowCountPerProcessor;
                for (unsigned long int ii = 0; ii < (RowCountPerProcessor); ii++)
                {
                    for (unsigned long int jj = 0; jj < MatrixColumnCount; jj++)
                    {
                        Temp_d_1->data_point[ii][jj] = d_1->data_point[ProcessorRowStartndex + ii][jj];
                        Temp_d_2->data_point[ii][jj] = d_2->data_point[ProcessorRowStartndex + ii][jj];
                    }
                }
                // for (unsigned long int i = 0; i < (RowCountPerProcessor); i++)
                // printf("AAA processor-%d\t _%lu__%lu_\t%f\t + %f\t = %f\n", ProcessorId, i, RowCountPerProcessor, (Temp_d_1->data_point[i][0]), (Temp_d_2->data_point[i][0]), (Temp_d_3->data_point[i][0]));

                MPI_Send(&RowCountPerProcessor, 1, MPI_UNSIGNED_LONG, i, 0, MPI_COMM_WORLD);
                MPI_Send(&MatrixColumnCount, 1, MPI_UNSIGNED_LONG, i, 10, MPI_COMM_WORLD);
                MPI_Send(&(Temp_d_1->data_point[0][0]), RowCountPerProcessor * MatrixColumnCount, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
                MPI_Send(&(Temp_d_2->data_point[0][0]), RowCountPerProcessor * MatrixColumnCount, MPI_DOUBLE, i, 2, MPI_COMM_WORLD);
            }

            unsigned long int FinalProcessorRowStartndex = (NoOfProcesors - 1) * RowCountPerProcessor;
            unsigned long int RemainingElements = d_1->rows - FinalProcessorRowStartndex;
            // printf("Final Processor %lu  %lu \n", FinalProcessorRowStartndex, RemainingElements);
            for (unsigned long int ii = 0; ii < (RemainingElements); ii++)
            {
                for (unsigned long int jj = 0; jj < MatrixColumnCount; jj++)
                {
                    Temp_d_1->data_point[ii][jj] = d_1->data_point[FinalProcessorRowStartndex + ii][jj];
                    Temp_d_2->data_point[ii][jj] = d_2->data_point[FinalProcessorRowStartndex + ii][jj];
                }
            }
            // for (unsigned long int i = 0; i < (RemainingElements); i++)
            // printf("BBB processor-%d\t _%lu__%lu_\t%f\t + %f\t = %f\n", ProcessorId, i, RemainingElements, (Temp_d_1->data_point[i][0]), (Temp_d_2->data_point[i][0]), (Temp_d_3->data_point[i][0]));

            MPI_Send(&RemainingElements, 1, MPI_UNSIGNED_LONG, NoOfProcesors - 1, 0, MPI_COMM_WORLD);
            MPI_Send(&MatrixColumnCount, 1, MPI_UNSIGNED_LONG, NoOfProcesors - 1, 10, MPI_COMM_WORLD);
            MPI_Send(&(Temp_d_1->data_point[0][0]), RemainingElements * MatrixColumnCount, MPI_DOUBLE, NoOfProcesors - 1, 1, MPI_COMM_WORLD);
            MPI_Send(&(Temp_d_2->data_point[0][0]), RemainingElements * MatrixColumnCount, MPI_DOUBLE, NoOfProcesors - 1, 2, MPI_COMM_WORLD);
        }
        // 1st Processor
        for (unsigned long int ii = 0; ii < RowCountPerProcessor; ii++)
        {
            for (unsigned long int jj = 0; jj < MatrixColumnCount; jj++)
            {
                d_3->data_point[ii][jj] = d_1->data_point[ii][jj] + d_2->data_point[ii][jj];
            }
        }
        if (NoOfProcesors > 1)
        {
            // Intermediate Processors
            for (unsigned long int i = 1; i < (NoOfProcesors - 1); i++)
            {
                MPI_Recv(&(Temp_d_3->data_point[0][0]), RowCountPerProcessor * MatrixColumnCount, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for (unsigned long int ii = 0; ii < RowCountPerProcessor; ii++)
                {
                    for (unsigned long int jj = 0; jj < MatrixColumnCount; jj++)
                    {
                        d_3->data_point[(i * RowCountPerProcessor) + ii][jj] = Temp_d_3->data_point[ii][jj];
                    }
                }
            }

            // Last Processor (to handle last few elements)
            unsigned long int FinalProcessorStartndex = (NoOfProcesors - 1) * RowCountPerProcessor;
            unsigned long int RemainingElements = d_1->rows - FinalProcessorStartndex;
            MPI_Recv(&(Temp_d_3->data_point[0][0]), RemainingElements * MatrixColumnCount, MPI_DOUBLE, NoOfProcesors - 1, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (unsigned long int ii = 0; ii < RowCountPerProcessor; ii++)
            {
                for (unsigned long int jj = 0; jj < MatrixColumnCount; jj++)
                {
                    d_3->data_point[FinalProcessorStartndex + ii][jj] = Temp_d_3->data_point[ii][jj];
                }
            }
        }
        stop_timer(&start);
        fprintf(stderr, " (calculating answer)\n");

        print_data(d_3);

        free2d(Temp_d_1->data_point);
        free(Temp_d_1);
        free2d(Temp_d_2->data_point);
        free(Temp_d_2);
        free2d(Temp_d_3->data_point);
        free(Temp_d_3);

        free_data(d_1);
        free_data(d_2);
        free_data(d_3);
    }
    else // Slave Processors
    {
        unsigned long int RowReceiveCount;
        unsigned long int ColumnReceiveCount;
        MPI_Recv(&RowReceiveCount, 1, MPI_UNSIGNED_LONG, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&ColumnReceiveCount, 1, MPI_UNSIGNED_LONG, 0, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        int PId;
        MPI_Comm_rank(MPI_COMM_WORLD, &PId);
        //printf("Received ProcessorId - %d ; RowCount - %lu ; ColumnCount %lu \n", PId, RowReceiveCount, ColumnReceiveCount);

        d_1 = malloc(sizeof(data_struct));
        d_2 = malloc(sizeof(data_struct));
        d_3 = malloc(sizeof(data_struct));
        d_1->data_point = alloc2d(RowReceiveCount, ColumnReceiveCount);
        d_2->data_point = alloc2d(RowReceiveCount, ColumnReceiveCount);
        d_3->data_point = alloc2d(RowReceiveCount, ColumnReceiveCount);

        MPI_Recv(&(d_1->data_point[0][0]), RowReceiveCount * ColumnReceiveCount, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&(d_2->data_point[0][0]), RowReceiveCount * ColumnReceiveCount, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for (unsigned long int ii = 0; ii < RowReceiveCount; ii++)
        {
            for (unsigned long int jj = 0; jj < ColumnReceiveCount; jj++)
            {
                d_3->data_point[ii][jj] = d_1->data_point[ii][jj] + d_2->data_point[ii][jj];
                //printf("Received PId-%d [%lu][%lu] = %f + %f = %f\n", PId, ii, jj, d_1->data_point[ii][jj], d_2->data_point[ii][jj], d_3->data_point[ii][jj]);
            }
        }
        MPI_Send(&(d_3->data_point[0][0]), RowReceiveCount * ColumnReceiveCount, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
        // free_data(d_1);
        free2d(d_1->data_point);
        free(d_1);
        free2d(d_2->data_point);
        free(d_2);
        free2d(d_3->data_point);
        free(d_3);
    }

    MPI_Finalize();
}