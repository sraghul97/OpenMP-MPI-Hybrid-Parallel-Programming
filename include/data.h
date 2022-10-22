/* General purpose struct for both Matrix and Vectors  */
/* The dataset has a rectangular matrix (nxn), vector size is (nx1) */

typedef struct {
        unsigned int rows;
        unsigned int cols;
        double **data_point;
} data_struct;

data_struct *get_data_struct(char matrix[]); //reads the data for file and return a nxn matrix
void print_data(data_struct *data_to_print); //prints the data (matrix/vector)
void free_data(data_struct *data_to_free);  //frees up allocated memory
data_struct *allocate_data_struct(unsigned int RowCount,unsigned int ColumnCount);
