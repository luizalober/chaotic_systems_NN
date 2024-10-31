#ifndef AUXILIARY_AND_RECURRENCE_FUNCTIONS_H
#define AUXILIARY_AND_RECURRENCE_FUNCTIONS_H

typedef double (*distance_function)(double[], double[], int);

void print_prog(double percentage);
int alloc_1d_double(double **x, int n);
int alloc_2d_double(double ***x, int n, int m);
int alloc_2d_int(int ***x, int n, int m);
int alloc_1d_int(int **x, int n);
int dealloc_1d_double(double **x);
int dealloc_1d_int(int **x);
int dealloc_2d_double(double ***x, int n);
int dealloc_2d_int(int ***x, int n);
double constrainAngle(double *x);
void create_folder_if_not_exists(const char* folderName);
int create_folders_in_path_if_not_exists(const char* path);
void remove_dot(char *result, size_t resultSize, double number);
double get_time_elapsed_minutes(struct timeval start, struct timeval end);
void embed_time_series(double **embedded_time_series, double *time_series, int n, int dim, int tau);
void manhattan_distance(double **embedded_time_series, int n, int dim, double **D);
void euclidean_distance(double **embedded_time_series, int n, int dim, double **D);
void supremum_distance(double **embedded_time_series, int n, int dim, double **D);
void modulated_distance(double **embedded_time_series, int n, int dim, double **D);
void calculate_distance_matrix(double **embedded_time_series, int n, int dim, double **D, const char *norm);
void compute_recurrence_matrix(double **D, int n, int **R, double threshold);
double compute_recurrence_rate(int **R, int n);
double adjust_threshold_via_recurrence_rate(double **D, int n, int **R, double desired_percentage, int *status);

#endif