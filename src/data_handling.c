#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include "data.h"
/* arguments supplied from the test script  */
data_struct *get_data_struct(char data[])
{
    data_struct *d = malloc(sizeof(data_struct));
    d->rows = 0;
    d->cols = 0;
    FILE *myfile = fopen(data, "r");
    /* check the test script for the files */
    if (myfile == NULL)
    {
        printf("Error: The file you entered could not be found.\n");
        exit(EXIT_FAILURE);
    }

    /* Count the number of rows and coloumns */
    /* Matrix nxn */
    /* Vector nx1 */
    int ch = 0;
    do
    {
        ch = fgetc(myfile);

        if (d->rows == 0 && ch == '\t')
            d->cols++;

        if (ch == '\n')
            d->rows++;

    } while (ch != EOF);

    /* Populate the data struct and allocate memory */
    d->cols++;
    d->data_point = calloc(d->rows, sizeof(double *));
    int i;
    for (i = 0; i < d->rows; ++i)
        d->data_point[i] = calloc(d->cols, sizeof(double));

    rewind(myfile);
    int x, y;

    /* The the data and return the data structure */
    for (x = 0; x < d->rows; x++)
    {
        for (y = 0; y < d->cols; y++)
        {
            if (!fscanf(myfile, "%lf", &d->data_point[x][y]))
                break;
        }
    }
    fclose(myfile);
    return d;
}

void print_data(data_struct *data_to_print)
{
    /* Write function to print the matrix/Vector */
    for (int RowCount = 0; RowCount < data_to_print->rows; RowCount++)
    {
        for (int ColumnCount = 0; ColumnCount < data_to_print->cols; ColumnCount++)
        {
            printf("%f ", data_to_print->data_point[RowCount][ColumnCount]);
            // fprintf(stderr, "%f ", data_to_print->data_point[RowCount][ColumnCount]);
        }
        printf("\n");
        // fprintf(stderr, "\n");
    }
}

/* Check how allocated data is free */
void free_data(data_struct *data_to_free)
{

    for (int i = 0; i < data_to_free->rows; i++)
    {
        free(data_to_free->data_point[i]);
    }
    free(data_to_free->data_point);
    free(data_to_free);
}

data_struct *allocate_data_struct(unsigned int RowCount,unsigned int ColumnCount)
{
    data_struct *data_to_allocate = malloc(sizeof(data_struct));
    data_to_allocate->rows = RowCount;
    data_to_allocate->cols = ColumnCount;

    data_to_allocate->data_point = calloc(data_to_allocate->rows, sizeof(double *));

    for (int i = 0; i < data_to_allocate->rows; ++i)
        data_to_allocate->data_point[i] = calloc(data_to_allocate->cols, sizeof(double));

    for (int x = 0; x < data_to_allocate->rows; x++)
    {
        for (int y = 0; y < data_to_allocate->cols; y++)
            data_to_allocate->data_point[x][y] = 0;
    }
    return data_to_allocate;
}
