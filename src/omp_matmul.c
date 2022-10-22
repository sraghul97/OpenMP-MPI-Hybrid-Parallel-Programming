#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include "data.h"
#include <omp.h>
#include "timer.h"

int main(int argc, char **argv)
{
    if (argc != 4)
    {
        printf("Usage: <num_thread> <vec_a> <vec_b>.\n");
        exit(EXIT_FAILURE);
    }

    // int numThreads = omp_get_max_threads();
    // convinient for use in partitioning
    int numThreads = atoi(argv[1]);
    omp_set_num_threads(numThreads);

    struct timespec start;
    start_timer(&start);

    data_struct *d_1 = get_data_struct(argv[2]);
    data_struct *d_2 = get_data_struct(argv[3]);
    data_struct *d_3 = allocate_data_struct(d_1->rows, d_2->cols);

    stop_timer(&start);
    fprintf(stderr, " (reading input)\n");

    if (d_1->cols != d_2->rows)
    {
        printf("ERROR: Matrix dimension mismatch.\n");
        exit(EXIT_FAILURE);
    }

    start_timer(&start);
    /* TODO: Implement Matrix Multiplication  openMP here */
    #pragma omp parallel for schedule(guided) collapse(3)
    for (unsigned long int RowCount = 0; RowCount < d_1->rows; RowCount++)
    {
        for (unsigned long int ColumnCount = 0; ColumnCount < d_2->cols; ColumnCount++)
        {
            for (unsigned long int SummationCount = 0; SummationCount < d_1->cols; SummationCount++)
            {
                #pragma omp atomic
                d_3->data_point[RowCount][ColumnCount] += (d_1->data_point[RowCount][SummationCount] * d_2->data_point[SummationCount][ColumnCount]);
            }
        }
    }
    stop_timer(&start);
    fprintf(stderr, " (calculating answer)\n");

    start_timer(&start);
    /* TODO: Print output */
    print_data(d_3);
    stop_timer(&start);
    fprintf(stderr, " (printing output)\n");
    free_data(d_1);
    free_data(d_2);
    free_data(d_3); 
}
