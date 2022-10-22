#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <malloc.h>
#include <string.h>
#include "data.h"
#include "timer.h"
#include <omp.h>

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
    if (argc != 4)
    {
        printf("Usage: <num_thread> <vec_a> <vec_b>.%d\n", argc);
        exit(EXIT_FAILURE);
    }
    data_struct *d_1;
    data_struct *d_2;
    data_struct *d_3;

    struct timespec start;

    int ProcessorId, NoOfProcesors;
    unsigned long int D1RowCountPerProcessor;
    unsigned long int D2ColumnCount;
    unsigned long int D1cD2rCount;
    MPI_Init(NULL, NULL);

    MPI_Comm_size(MPI_COMM_WORLD, &NoOfProcesors);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcessorId);

    if (ProcessorId == 0) // Master Processor
    {
        d_1 = get_data_struct(argv[2]);
        d_2 = get_data_struct(argv[3]);
        d_3 = allocate_data_struct(d_1->rows, d_2->cols);

        if (d_1->cols != d_2->rows)
        {
            printf("ERROR: The dimiension of Matrix is not correct.\n");
            exit(EXIT_FAILURE);
        }
        D1RowCountPerProcessor = d_1->rows / NoOfProcesors;
        D1cD2rCount = d_1->cols;
        D2ColumnCount = d_2->cols;

        data_struct *Temp_d_1;
        data_struct *Temp_d_2;
        data_struct *Temp_d_3;
        Temp_d_1 = malloc(sizeof(data_struct));
        Temp_d_2 = malloc(sizeof(data_struct));
        Temp_d_3 = malloc(sizeof(data_struct));
        Temp_d_1->data_point = alloc2d(D1RowCountPerProcessor, D1cD2rCount);
        Temp_d_2->data_point = alloc2d(D1cD2rCount, D2ColumnCount);
        Temp_d_3->data_point = alloc2d(D1RowCountPerProcessor, D2ColumnCount);
        start_timer(&start);

        if (NoOfProcesors > 1)
        {
            for (unsigned long int i = 1; i < (NoOfProcesors - 1); i++)
            {

#pragma omp parallel for schedule(guided) collapse(2)
                for (unsigned long int ii = 0; ii < (D1RowCountPerProcessor); ii++)
                {
                    for (unsigned long int jj = 0; jj < D1cD2rCount; jj++)
                    {
                        // printf("Temp_D_1[%lu][%lu] = d_1[%lu][%lu]\t \t %lu %lu\n", ii, jj, D1ProcessorRowStartndex + ii, jj, D1RowCountPerProcessor, D1cD2rCount);
                        unsigned long int D1ProcessorRowStartndex = i * D1RowCountPerProcessor;
                        Temp_d_1->data_point[ii][jj] = d_1->data_point[D1ProcessorRowStartndex + ii][jj];
                    }
                }
                // printf("A22A %lu\n", i);

                //  unsigned long int D2ProcessorColumnStartndex = i * D2ColumnCountPerProcessor;
#pragma omp parallel for schedule(guided) collapse(2)
                for (unsigned long int ii = 0; ii < (D1cD2rCount); ii++)
                {
                    for (unsigned long int jj = 0; jj < D2ColumnCount; jj++)
                    {
                        // printf("Temp_D_2[%lu][%lu] = d_2[%lu][%lu]\t \t %lu %lu\n", ii, jj, ii, D2ProcessorColumnStartndex + jj, D1cD2rCount, D2ColumnCountPerProcessor);
                        Temp_d_2->data_point[ii][jj] = d_2->data_point[ii][jj];
                    }
                }

                // for (unsigned long int i = 0; i < (RowCountPerProcessor); i++)
                // printf("AAA processor-%d\t _%lu__%lu_\t%f\t + %f\t = %f\n", ProcessorId, i, RowCountPerProcessor, (Temp_d_1->data_point[i][0]), (Temp_d_2->data_point[i][0]), (Temp_d_3->data_point[i][0]));

                // printf("A33A \n");
                MPI_Send(&D1RowCountPerProcessor, 1, MPI_UNSIGNED_LONG, i, 0, MPI_COMM_WORLD);
                MPI_Send(&D2ColumnCount, 1, MPI_UNSIGNED_LONG, i, 1, MPI_COMM_WORLD);
                MPI_Send(&D1cD2rCount, 1, MPI_UNSIGNED_LONG, i, 2, MPI_COMM_WORLD);
                MPI_Send(&(Temp_d_1->data_point[0][0]), D1RowCountPerProcessor * D1cD2rCount, MPI_DOUBLE, i, 3, MPI_COMM_WORLD);
                MPI_Send(&(Temp_d_2->data_point[0][0]), D1cD2rCount * D2ColumnCount, MPI_DOUBLE, i, 4, MPI_COMM_WORLD);
            }
            // printf("A11A \n");

            unsigned long int FinalProcessorD1RowStartndex = (NoOfProcesors - 1) * D1RowCountPerProcessor;
            unsigned long int D1RemainingRows = d_1->rows - FinalProcessorD1RowStartndex;
            // unsigned long int FinalProcessorD2ColumnStartndex = (NoOfProcesors - 1) * D2ColumnCountPerProcessor;

            // unsigned long int D2RemainingColumns = d_2->cols - FinalProcessorD2ColumnStartndex;
            //  printf("Final Processor %lu  %lu \n", FinalProcessorRowStartndex, RemainingElements);
#pragma omp parallel for schedule(guided) collapse(2)
            for (unsigned long int ii = 0; ii < (D1RemainingRows); ii++)
            {
                for (unsigned long int jj = 0; jj < D1cD2rCount; jj++)
                {
                    Temp_d_1->data_point[ii][jj] = d_1->data_point[FinalProcessorD1RowStartndex + ii][jj];
                }
            }

            /*for (unsigned long int ii = 0; ii < (D1cD2rCount); ii++)
            {
                for (unsigned long int jj = 0; jj < D2RemainingColumns; jj++)
                {
                    Temp_d_2->data_point[ii][jj] = d_2->data_point[ii][FinalProcessorD2ColumnStartndex + jj];
                }
            }*/

            // for (unsigned long int i = 0; i < (RemainingElements); i++)
            // printf("BBB processor-%d\t _%lu__%lu_\t%f\t + %f\t = %f\n", ProcessorId, i, RemainingElements, (Temp_d_1->data_point[i][0]), (Temp_d_2->data_point[i][0]), (Temp_d_3->data_point[i][0]));

            MPI_Send(&D1RemainingRows, 1, MPI_UNSIGNED_LONG, NoOfProcesors - 1, 0, MPI_COMM_WORLD);
            MPI_Send(&D2ColumnCount, 1, MPI_UNSIGNED_LONG, NoOfProcesors - 1, 1, MPI_COMM_WORLD);
            MPI_Send(&D1cD2rCount, 1, MPI_UNSIGNED_LONG, NoOfProcesors - 1, 2, MPI_COMM_WORLD);
            MPI_Send(&(Temp_d_1->data_point[0][0]), D1RemainingRows * D1cD2rCount, MPI_DOUBLE, NoOfProcesors - 1, 3, MPI_COMM_WORLD);
            MPI_Send(&(Temp_d_2->data_point[0][0]), D1cD2rCount * D2ColumnCount, MPI_DOUBLE, NoOfProcesors - 1, 4, MPI_COMM_WORLD);
        }
        // printf("A2A \n");

        // 1st Processor
#pragma omp parallel for schedule(guided) collapse(3)
        for (unsigned long int RowCount = 0; RowCount < D1RowCountPerProcessor; RowCount++)
        {
            for (unsigned long int ColumnCount = 0; ColumnCount < D2ColumnCount; ColumnCount++)
            {
                for (unsigned long int SummationCount = 0; SummationCount < D1cD2rCount; SummationCount++)
                {
#pragma omp atomic
                    d_3->data_point[RowCount][ColumnCount] += (d_1->data_point[RowCount][SummationCount] * d_2->data_point[SummationCount][ColumnCount]);
                    // d_3->data_point[ii][jj] = d_1->data_point[ii][jj] + d_2->data_point[ii][jj];
                }
            }
        }
        if (NoOfProcesors > 1)
        {
            // Intermediate Processors
            for (unsigned long int i = 1; i < (NoOfProcesors - 1); i++)
            {
                MPI_Recv(&(Temp_d_3->data_point[0][0]), D1RowCountPerProcessor * D2ColumnCount, MPI_DOUBLE, i, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#pragma omp parallel for schedule(guided) collapse(2)
                for (unsigned long int ii = 0; ii < D1RowCountPerProcessor; ii++)
                {
                    for (unsigned long int jj = 0; jj < D2ColumnCount; jj++)
                    {
                        d_3->data_point[(i * D1RowCountPerProcessor) + ii][jj] = Temp_d_3->data_point[ii][jj];
                        // printf("D_3[%lu][%lu] = Temp_d_3[%lu][%lu] = %f \t \t %lu %lu %lu\n", (i * D1RowCountPerProcessor) + ii, jj, ii, jj, Temp_d_3->data_point[ii][jj], D1RowCountPerProcessor, D2ColumnCount, i);
                    }
                }
            }

            // Last Processor (to handle last few elements)
            unsigned long int FinalProcessorStartndex = (NoOfProcesors - 1) * D1RowCountPerProcessor;
            unsigned long int RemainingElements = d_1->rows - FinalProcessorStartndex;
            MPI_Recv(&(Temp_d_3->data_point[0][0]), RemainingElements * D2ColumnCount, MPI_DOUBLE, NoOfProcesors - 1, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#pragma omp parallel for schedule(guided) collapse(2)
            for (unsigned long int ii = 0; ii < RemainingElements; ii++)
            {
                for (unsigned long int jj = 0; jj < D2ColumnCount; jj++)
                {
                    d_3->data_point[FinalProcessorStartndex + ii][jj] = Temp_d_3->data_point[ii][jj];
                    // printf("D_3[%lu][%lu] = Temp_d_3[%lu][%lu] = %f \t \t %lu %lu %d Final Processor\n", FinalProcessorStartndex + ii, jj, ii, jj, Temp_d_3->data_point[ii][jj], D1RowCountPerProcessor, D2ColumnCount, NoOfProcesors - 1);
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
        unsigned long int D1RowCount;
        unsigned long int D2ColumnCount;
        unsigned long int D1cD2rCount;
        int PId;
        MPI_Recv(&D1RowCount, 1, MPI_UNSIGNED_LONG, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&D2ColumnCount, 1, MPI_UNSIGNED_LONG, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&D1cD2rCount, 1, MPI_UNSIGNED_LONG, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Comm_rank(MPI_COMM_WORLD, &PId);
        // printf("Received ProcessorId - %d ; D1RowCount - %lu ; D2ColumnCount %lu ; D1cD2R %lu\n", PId, D1RowCount, D2ColumnCount, D1cD2rCount);

        d_1 = malloc(sizeof(data_struct));
        d_2 = malloc(sizeof(data_struct));
        d_3 = malloc(sizeof(data_struct));
        d_1->data_point = alloc2d(D1RowCount, D1cD2rCount);
        d_2->data_point = alloc2d(D1cD2rCount, D2ColumnCount);
        d_3->data_point = alloc2d(D1RowCount, D2ColumnCount);

        MPI_Recv(&(d_1->data_point[0][0]), D1RowCount * D1cD2rCount, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&(d_2->data_point[0][0]), D1cD2rCount * D2ColumnCount, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        /*for (unsigned long int RowCount = 0; RowCount < D1RowCount; RowCount++)
        {
            for (unsigned long int ColumnCount = 0; ColumnCount < D1cD2rCount; ColumnCount++)
            {
                printf("%f \t",d_1->data_point[RowCount][ColumnCount]);
            }
            printf("__%d_d1 \n",ProcessorId);
        }

         for (unsigned long int RowCount = 0; RowCount < D1cD2rCount; RowCount++)
        {
            for (unsigned long int ColumnCount = 0; ColumnCount < D2ColumnCount; ColumnCount++)
            {
                printf("%f \t",d_2->data_point[RowCount][ColumnCount]);
            }
            printf("__%d_d2 \n",ProcessorId);
        }*/

#pragma omp parallel for schedule(guided) collapse(3)
        for (unsigned long int RowCount = 0; RowCount < D1RowCount; RowCount++)
        {
            for (unsigned long int ColumnCount = 0; ColumnCount < D2ColumnCount; ColumnCount++)
            {
                for (unsigned long int SummationCount = 0; SummationCount < D1cD2rCount; SummationCount++)
                {
#pragma omp atomic
                    d_3->data_point[RowCount][ColumnCount] += (d_1->data_point[RowCount][SummationCount] * d_2->data_point[SummationCount][ColumnCount]);
                    // d_3->data_point[ii][jj] = d_1->data_point[ii][jj] + d_2->data_point[ii][jj];
                    // printf("Received PId-%d [%lu][%lu] = %f * %f = %f\n", PId, RowCount, ColumnCount, d_1->data_point[RowCount][ColumnCount], d_2->data_point[RowCount][ColumnCount], d_3->data_point[RowCount][ColumnCount]);
                }
            }
        }
        MPI_Send(&(d_3->data_point[0][0]), D1RowCount * D2ColumnCount, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD);
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
