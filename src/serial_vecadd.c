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

    // Vector is of size (nx1), it will always return coloumn =1
    if (d_1->cols != 1 || d_2->cols != 1)
    {
        printf("ERROR: The dimiension of vector is not correct\n");
        exit(EXIT_FAILURE);
    }

    start_timer(&start);
    /* TODO: Implement serial Vector Addition */
    for (unsigned long int RowCount = 0; RowCount < d_3->rows; RowCount++)
     d_3->data_point[RowCount][0] = d_1->data_point[RowCount][0] + d_2->data_point[RowCount][0];

    stop_timer(&start);
    fprintf(stderr, " (calculating answer)\n");

    start_timer(&start);
    /* TODO: Print output */
    print_data(d_3);
    stop_timer(&start);
    fprintf(stderr, " (printing output)\n");
    free(d_1);
    free(d_2);
    free(d_3);
}
