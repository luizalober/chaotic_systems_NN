#include "common.h"
#include "auxiliary_and_recurrence_functions.h"

#define PBSTR "=================================================="
#define PBWIDTH 50

#define T 100 //Maximum characteres for file's names.

void print_prog(double percentage)
{
	int val = (int)(percentage * 100);
	int lpad = (int)(percentage * PBWIDTH);
	int rpad = PBWIDTH - lpad;
	printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
	fflush(stdout);
}

int alloc_1d_double(double **x, int n)
{
	*x = (double *)malloc(n * sizeof(double));
	return 0;
}

int alloc_2d_double(double ***x, int n, int m)
{
	*x = (double **)malloc(n * sizeof(double *));
	for (int i = 0; i < n; i++)
	{
		(*x)[i] = (double *)malloc(m * sizeof(double));
	}
	return 0;
}

int alloc_2d_int(int ***x, int n, int m)
{
	*x = (int **)malloc(n * sizeof(int *));
	for (int i = 0; i < n; i++)
	{
		(*x)[i] = (int *)malloc(m * sizeof(int));
	}
	return 0;
}

int alloc_1d_int(int **x, int n)
{
	*x = (int *)malloc(n * sizeof(int));
	return 0;
}

int dealloc_1d_double(double **x)
{
	free(*x);
	return 0;
}

int dealloc_1d_int(int **x)
{
	free(*x);
	return 0;
}

int dealloc_2d_double(double ***x, int n)
{
	for (int i = 0; i < n; i++)
	{
		free((*x)[i]);
	}
	free(*x);
	return 0;
}

int dealloc_2d_int(int ***x, int n)
{
	for (int i = 0; i < n; i++)
	{
		free((*x)[i]);
	}
	free(*x);
	return 0;
}

double constrainAngle(double *x)
{
    *x = fmod(*x, 2 * M_PI);
    if (*x < 0)
        *x += 2 * M_PI;
    return *x;
}

void create_folder_if_not_exists(const char* folder_name) 
{
    // Check if the folder exists
    struct stat st;
    if (stat(folder_name, &st) == -1) {
        // Folder does not exist, create it
        mkdir(folder_name, 0700); // 0700 gives full permissions, you can adjust it accordingly
        printf("Folder '%s' created!\n", folder_name);
    }
}

int create_folders_in_path_if_not_exists(const char* path) 
{
    char temp_path[1024];
    char* p = NULL;
    struct stat st;
    size_t len;
    int status = 0;

    // Copy the path to a temporary variable
    snprintf(temp_path, sizeof(temp_path), "%s", path);
    len = strlen(temp_path);

    // Iterate through the path, creating each part
    for (p = temp_path + 1; p < temp_path + len; p++) {
        if (*p == '/') {
            *p = '\0';  // Temporarily end the string here

            // Create the directory if it doesn't exist
            if (stat(temp_path, &st) == -1) {
                if (mkdir(temp_path, 0700) == 0) {
                    printf("Directory created: %s\n", temp_path);
                } else {
                    perror("Error creating folder");
                    return -1;  // Error creating folder
                }
            }

            *p = '/';  // Restore the original character
        }
    }

    // Create the final directory (if not created in the loop)
    if (stat(temp_path, &st) == -1) {
        if (mkdir(temp_path, 0700) == 0) {
            printf("Full path: %s\n", temp_path);
        } else {
            perror("Error creating folder");
            return -1;  // Error creating folder
        }
    }

    return status;
}

void remove_dot(char *result, size_t resultSize, double number)
{
    char buffer[50];

    // Convert double to string with high precision
    snprintf(buffer, sizeof(buffer), "%.16f", number);

    // Initialize result to empty string
    result[0] = '\0';

    // Copy characters from buffer to result, skipping the dot
    int j = 0;
    for (int i = 0; i < strlen(buffer) && j < resultSize - 1; i++) {
        if (buffer[i] != '.') {
            result[j++] = buffer[i];
        }
    }
    
    // Ensure null termination
    result[j] = '\0';
}

double get_time_elapsed_minutes(struct timeval start, struct timeval end) 
{
    long seconds, useconds;
    double timeElapsed;

    // Calculate the time elapsed in microseconds
    seconds = end.tv_sec - start.tv_sec;
    useconds = end.tv_usec - start.tv_usec;
    timeElapsed = seconds + useconds / 1000000.0;

    // Convert time elapsed to minutes
    return timeElapsed / 60.0;
}

void embed_time_series(double **embedded_time_series, double *time_series, int n, int dim, int tau) 
{
    for (int i = 0; i < n - (dim - 1) * tau; i++) {
        for (int j = 0; j < dim; j++) {
            embedded_time_series[i][j] = time_series[i + j * tau];
        }
    }
}

void manhattan_distance(double **embedded_time_series, int n, int dim, double **D) 
{
    double temp;
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            temp = 0.0;
            for (int k = 0; k < dim; k++) {
                temp += fabs(embedded_time_series[i][k] - embedded_time_series[j][k]);
            }
            D[i][j] = temp;
            D[j][i] = temp;
        }
    }
}

void euclidean_distance(double **embedded_time_series, int n, int dim, double **D) 
{
    double temp;
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            temp = 0.0;
            for (int k = 0; k < dim; k++) {
                temp += pow(embedded_time_series[i][k] - embedded_time_series[j][k], 2);
            }
            D[i][j] = sqrt(temp);
            D[j][i] = sqrt(temp);
        }
    }
}

void supremum_distance(double **embedded_time_series, int n, int dim, double **D) 
{
    double max, temp;
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            max = fabs(embedded_time_series[i][0] - embedded_time_series[j][0]);
            for (int k = 1; k < dim; k++) {
                temp = fabs(embedded_time_series[i][k] - embedded_time_series[j][k]);
                if (temp > max) {
                    max = temp;
                }
            }
            D[i][j] = max;
            D[j][i] = max;
        }
    }
}

void modulated_distance(double **embedded_time_series, int n, int dim, double **D) 
{
    double diff;
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            diff = 0.0;
            for (int k = 0; k < dim; k++) {
                double temp = fmod(fabs(embedded_time_series[i][k] - embedded_time_series[j][k]), 2 * M_PI);
                if (temp > M_PI) {
                    temp = 2 * M_PI - temp;
                }
                diff += temp * temp;
            }
            D[i][j] = sqrt(diff);
            D[j][i] = sqrt(diff);
        }
    }
}

void calculate_distance_matrix(double **embedded_time_series, int n, int dim, double **D, const char *norm) 
{
    if (strcmp(norm, "manhattan") == 0) manhattan_distance(embedded_time_series, n, dim, D);
    else if (strcmp(norm, "euclidean") == 0) euclidean_distance(embedded_time_series, n, dim, D);
    else if (strcmp(norm, "supremum") == 0) supremum_distance(embedded_time_series, n, dim, D);
    else if (strcmp(norm, "modulated") == 0) modulated_distance(embedded_time_series, n, dim, D);
    else
    {
        fprintf(stderr, "Unknown norm type: %s\n", norm);
        exit(1);
    }
}

void compute_recurrence_matrix(double **D, int n, int **R, double threshold) 
{
    for (int i = 0; i < n; i++) 
    {
        for (int j = 0; j < n; j++) 
        {
                if (D[i][j] < threshold) R[i][j] = 1;
                else R[i][j] = 0;
        }
    }
}

double compute_recurrence_rate(int **R, int n) 
{
    int count_recurrent_points = 0;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i != j && R[i][j] == 1) {
                count_recurrent_points++;
            }
        }
    }

    return (double) count_recurrent_points / (double) (n * n);
}

double adjust_threshold_via_recurrence_rate(double **D, int n, int **R, double desired_percentage, int *status) 
{
    if (n <= 0) {
        *status = -1; // Invalid size
        return -1.0;
    }

    double low = 0.0;
    double high = 0.0;

    // Compute the maximum distance in the distance matrix
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (D[i][j] > high) {
                high = D[i][j];
            }
        }
    }

    double mid, recurrence_rate;
    double closest_diff = 1e10; // Start with a large value
    double best_mid = 0.0;

    *status = 1; // Default to closest achievable

    while (high - low > 1e-10) { // Precision threshold
        mid = (low + high) / 2.0;
        compute_recurrence_matrix(D, n, R, mid);
        recurrence_rate = compute_recurrence_rate(R, n) * 100;

        double diff = fabs(recurrence_rate - desired_percentage);
        if (diff < closest_diff) {
            closest_diff = diff;
            best_mid = mid;
        }

        if (diff < 1e-2) { // If within acceptable precision
            *status = 0; // Exact match
            break;
        }

        if (recurrence_rate < desired_percentage) {
            low = mid;
        } else {
            high = mid;
        }
    }

    return *status == 0 ? mid : best_mid;
}