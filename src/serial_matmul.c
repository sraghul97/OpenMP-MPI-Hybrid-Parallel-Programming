#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include "data.h"
#include "timer.h"

int main(int argc, char **argv)
{
    if (argc != 3)
    {
        printf("ERROR: Please specify only 2 files.\n");
        exit(EXIT_FAILURE);
    }

    struct timespec start;
    start_timer(&start);

    data_struct *d_1 = get_data_struct(argv[1]);
    data_struct *d_2 = get_data_struct(argv[2]);
    data_struct *d_3 = allocate_data_struct(d_1->rows, d_2->cols);

    stop_timer(&start);
    fprintf(stderr, " (reading input)\n");

    if (d_1->cols != d_2->rows)
    {
        printf("ERROR: Matrix dimension mismatch.%d %d\n",d_1->cols,d_2->rows);
        exit(EXIT_FAILURE);
    }

    start_timer(&start);
    /* TODO: Implement serial Matrix Multiplication */
    for (unsigned long int RowCount = 0; RowCount < d_3->rows; RowCount++)
    {
        for (unsigned long int ColumnCount = 0; ColumnCount < d_2->cols; ColumnCount++)
        {
            d_3->data_point[RowCount][ColumnCount] = 0;
            for (unsigned long int SummationCount = 0; SummationCount < d_1->cols; SummationCount++)
            {
                 d_3->data_point[RowCount][ColumnCount] +=  (d_1->data_point[RowCount][SummationCount] * d_2->data_point[SummationCount][ColumnCount]);
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
}
