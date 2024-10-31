#include "common.h"
#include "rec_plots_chaotic_systems.h"
#include "auxiliary_and_recurrence_functions.h"

#define GNUPLOT "gnuplot --persist"

#define T 500

//GLOBAL PARAMETERS - DYNAMICS
double parameter;
double initial_conditions[]={0.1, M_PI}; //This is just an example, it will change in the function map_details()
int order = 9; 
int N = 100;
double k_c = 0.971635; // Critical Value for the Standard Map (transition from local to global chaos)

//GLOBAL PARAMETERS - RECURRENCE
double eps = 0.5;
int emb_dim = 2;
int tau = 1;
double percentage = 10; // Percentage of Recurrences (RR)
bool define_threshold_via_percentage = true;

char filename1[T], filename2[T], filename3[T], filename4[T], filename5[T];
char filename6[T], filename7[T], filename8[T], filename9[T], filename10[T];

//Path to Recurrence Plots
const char *results_path = "results/test";

MAP log_map;
MAP std_map;

void logistic_map_eqs(double *x0, double *x1, double *r)
{
    x1[0] = (*r) * x0[0] * (1.0 - x0[0]);
}

void standard_map_eqs(double *x0, double *x1, double *k)
{
    x1[1] = x0[1] + (*k) * sin(x0[0]);
    constrainAngle(&x1[1]);
    x1[0] = x0[0] + x1[1] + M_PI;
    constrainAngle(&x1[0]);
}

void evolve_map(MAP *m, void *orbit, int N, double k)
{
    if (m->dimension == 1)
    {
        double x0[1], x1[1];
        double *orbit_1d = (double *)orbit;

        x0[0] = m->initial_conditions[0];

        for (int i = 0; i < N; i++)
        {
            logistic_map_eqs(x0, x1, &k);
            x0[0] = x1[0];
            orbit_1d[i] = x1[0];
        }
    }
    else if (m->dimension == 2)
    {
        double x0[2], x1[2];
        double **orbit_2d = (double **)orbit;

        x0[0] = m->initial_conditions[0];
        x0[1] = m->initial_conditions[1];

        for (int i = 0; i < N; i++)
        {
            standard_map_eqs(x0, x1, &k);
            x0[0] = x1[0];
            x0[1] = x1[1];
            orbit_2d[i][0] = x1[0];
            orbit_2d[i][1] = x1[1];
        }
    }
}

void give_orbit(MAP *m, void *orbit, int N) 
{
    FILE *out = fopen("orbit.dat", "w");
    if (m->dimension == 1) 
    {
        double *orbit_1d = (double *)orbit;
        for (int i = 0; i < N; i++) 
        {
            fprintf(out, "%1.12f\n", orbit_1d[i]);
        }
    } 
    else if (m->dimension == 2)
    {
        double **orbit_2d = (double **)orbit;
        for (int i = 0; i < N; i++) 
        {
            fprintf(out, "%1.12f %1.12f\n", orbit_2d[i][0], orbit_2d[i][1]);
        }
    }
    fclose(out);
}

int give_recurrence_plot(int **R, int N)
{
    int x_rp[N], y_rp[N];

    sprintf(filename1, "rm.dat");
    FILE *rm = fopen(filename1, "w");
    sprintf(filename2, "rp.dat");
    FILE *rp = fopen(filename2, "w");

    if (rm == NULL || rp == NULL) 
    {
        fprintf(stderr, "Error opening files for writing.\n");
        return 1;
    }

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            fprintf(rm, "%d ", R[i][j]);
            if (R[i][j] == 1)
            {
                x_rp[i]=i;
                y_rp[j]=j;
                fprintf(rp,"%d %d\n", x_rp[i], y_rp[j]);
            }
        }
        fprintf(rm, "\n");
    }

    fclose(rm);
    fclose(rp);
    return 0;
}

int give_bifurcation_diagram(MAP *m, int N)
{
    double parameter_min;
    double parameter_max;
    int number_of_parameters;

    FILE *bifurcation_file;
    sprintf(filename3, "bifurcation_%s.dat", m->name);
    bifurcation_file = fopen(filename3, "w");

    if (bifurcation_file == NULL) 
    {
        fprintf(stderr, "Error opening bifurcation file for writing.\n");
        return 1;
    }

    parameter_min = 3;
    parameter_max = 4; 
    number_of_parameters = 10000;

    double parameter_step = (parameter_max - parameter_min) / number_of_parameters;

    for (double parameter = parameter_min; parameter <= parameter_max; parameter += parameter_step)
    {
        m->parameter = parameter;

        // Initialize the orbit array
        double orbit[N];
        evolve_map(m, orbit, N, m->parameter);

        // Skip the first few iterations (transients) to focus on the attractor
        int transient = N / 2;

        for (int i = transient; i < N; i++) 
        {
            // Save the parameter and the orbit point to the file
            fprintf(bifurcation_file, "%1.16f %1.16f\n", parameter, orbit[i]);
        }
    }

    fclose(bifurcation_file);

    plot_gnuplot_bifurcation_diagram("bifurcation.png", m->name);

    return 0;
}

int ensemble_recurrence_analysis_varying_parameter(MAP *m)
{
    double parameter_min;
    double parameter_max;
    int number_of_parameters;

    double **orbit_2d = NULL;
    double *orbit_1d = NULL;
    double **DM = NULL;
    int **RM = NULL;

    if (strcmp(m->name, "logistic") == 0)
    {
        parameter_min = 3.0;
        parameter_max = 4.0; 
        number_of_parameters = 10000;

        double parameter_step = (parameter_max - parameter_min) / number_of_parameters;

        int p = 0;

        for (double parameter = parameter_min; parameter <= parameter_max; parameter += parameter_step)
        {
            alloc_1d_double(&orbit_1d, N);
            evolve_map(m, orbit_1d, N, parameter);
            double **embedded_time_series = NULL;
            alloc_2d_double(&embedded_time_series, N, emb_dim);
            embed_time_series(embedded_time_series, orbit_1d, N, emb_dim, tau);
            alloc_2d_double(&DM, N, N);
            calculate_distance_matrix(embedded_time_series, N, emb_dim, DM, m->distance_metric);
            dealloc_2d_double(&embedded_time_series, N);
            dealloc_1d_double(&orbit_1d);
            alloc_2d_int(&RM, N, N);
            int status;
            double threshold = define_threshold_via_percentage ? adjust_threshold_via_recurrence_rate(DM, N, RM, percentage, &status) : eps;
            compute_recurrence_matrix(DM, N, RM, threshold);
            give_recurrence_plot(RM, N);

            size_t folder_name_size = strlen(results_path) + 256;
            char *folder_name = (char *)malloc(folder_name_size * sizeof(char));
            snprintf(folder_name, folder_name_size, "%s/%s/x0=%1.2f", results_path, m->name, m->initial_conditions[0]);
            create_folders_in_path_if_not_exists(folder_name);

            char parameter_no_dot[256];
            remove_dot(parameter_no_dot, sizeof(parameter_no_dot), parameter);
            //For making a gif:
            snprintf(filename3, sizeof(filename3), "%s/%d.png", folder_name, p);
            //snprintf(filename3, sizeof(filename3), "%s/%s.png", folder_name, parameter_no_dot);

            plot_gnuplot_bulk_rps(filename3, N);
            
            dealloc_2d_int(&RM, N);
            dealloc_2d_double(&DM, N);

            print_prog((double)(p) / (double)number_of_parameters);
            p++;
        }
    }

    if (strcmp(m->name, "standard") == 0)
    {
        parameter_min = 0.1;
        parameter_max = 2.0; 
        number_of_parameters = 10000;

        double parameter_step = (parameter_max - parameter_min) / number_of_parameters;

        int p = 0;

        for (double parameter = parameter_min; parameter <= parameter_max; parameter += parameter_step)
        {

            alloc_2d_double(&orbit_2d, N, 2);
            evolve_map(m, orbit_2d, N, parameter);
            alloc_2d_double(&DM, N, N);
            calculate_distance_matrix(orbit_2d, N, emb_dim, DM, m->distance_metric);
            dealloc_2d_double(&orbit_2d, N);
            alloc_2d_int(&RM, N, N);
            int status;
            double threshold = define_threshold_via_percentage ? adjust_threshold_via_recurrence_rate(DM, N, RM, percentage, &status) : eps;
            compute_recurrence_matrix(DM, N, RM, threshold);
            give_recurrence_plot(RM, N);
            
            size_t folder_name_size = strlen(results_path) + 256;
            char *folder_name = (char *)malloc(folder_name_size * sizeof(char));
            snprintf(folder_name, folder_name_size, "%s/%s/order=%d/delta=%1.2f", results_path, m->name, order, m->delta);
            create_folders_in_path_if_not_exists(folder_name);

            char parameter_no_dot[256];
            remove_dot(parameter_no_dot, sizeof(parameter_no_dot), parameter);
            //For making a gif:
            snprintf(filename3, sizeof(filename3), "%s/%d.png", folder_name, p);
            //snprintf(filename3, sizeof(filename3), "%s/%s.png", folder_name, parameter_no_dot);

            plot_gnuplot_bulk_rps(filename3, N);
            
            dealloc_2d_int(&RM, N);
            dealloc_2d_double(&DM, N);

            print_prog((double)(p) / (double)number_of_parameters);
            p++;
        }
    }

    return 0;
}

int ensemble_recurrence_analysis_random_parameter(MAP *m)
{
    double parameter_min;
    double parameter_max;
    int number_of_parameters;

    double **orbit_2d = NULL;
    double *orbit_1d = NULL;
    double **DM = NULL;
    int **RM = NULL;

    srand(time(NULL)); // Seed the random number generator

    if (strcmp(m->name, "logistic") == 0)
    {
        parameter_min = 3.56;
        parameter_max = 4.0; 
        number_of_parameters = 1000;

        for (int p = 0; p < number_of_parameters; p++)
        {
            double parameter = parameter_min + ((double)rand() / RAND_MAX) * (parameter_max - parameter_min);

            alloc_1d_double(&orbit_1d, N);
            evolve_map(m, orbit_1d, N, parameter);
            double **embedded_time_series = NULL;
            alloc_2d_double(&embedded_time_series, N, emb_dim);
            embed_time_series(embedded_time_series, orbit_1d, N, emb_dim, tau);
            alloc_2d_double(&DM, N, N);
            calculate_distance_matrix(embedded_time_series, N, emb_dim, DM, m->distance_metric);
            dealloc_2d_double(&embedded_time_series, N);
            dealloc_1d_double(&orbit_1d);
            alloc_2d_int(&RM, N, N);
            int status;
            double threshold = define_threshold_via_percentage ? adjust_threshold_via_recurrence_rate(DM, N, RM, percentage, &status) : eps;
            compute_recurrence_matrix(DM, N, RM, threshold);
            give_recurrence_plot(RM, N);

            size_t folder_name_size = strlen(results_path) + 256;
            char *folder_name = (char *)malloc(folder_name_size * sizeof(char));
            snprintf(folder_name, folder_name_size, "%s/%s/x0=%1.2f", results_path, m->name, m->initial_conditions[0]);
            create_folders_in_path_if_not_exists(folder_name);

            char parameter_no_dot[256];
            remove_dot(parameter_no_dot, sizeof(parameter_no_dot), parameter);
            //For making a gith:
            //snprintf(filename3, sizeof(filename3), "%s/%d.png", folder_name, p);
            snprintf(filename3, sizeof(filename3), "%s/%s.png", folder_name, parameter_no_dot);

            plot_gnuplot_bulk_rps(filename3, N);
            
            dealloc_2d_int(&RM, N);
            dealloc_2d_double(&DM, N);

            print_prog((double)(p) / (double)number_of_parameters);
        }
    }

    if (strcmp(m->name, "standard") == 0)
    {
        parameter_min = 0.1;
        parameter_max = 2.0; 
        number_of_parameters = 1000;

        for (int p = 0; p < number_of_parameters; p++)
        {
            double parameter = parameter_min + ((double)rand() / RAND_MAX) * (parameter_max - parameter_min);

            alloc_2d_double(&orbit_2d, N, 2);
            evolve_map(m, orbit_2d, N, parameter);
            alloc_2d_double(&DM, N, N);
            calculate_distance_matrix(orbit_2d, N, emb_dim, DM, m->distance_metric);
            dealloc_2d_double(&orbit_2d, N);
            alloc_2d_int(&RM, N, N);
            int status;
            double threshold = define_threshold_via_percentage ? adjust_threshold_via_recurrence_rate(DM, N, RM, percentage, &status) : eps;
            compute_recurrence_matrix(DM, N, RM, threshold);
            give_recurrence_plot(RM, N);

            size_t folder_name_size = strlen(results_path) + 256;
            char *folder_name = (char *)malloc(folder_name_size * sizeof(char));
            snprintf(folder_name, folder_name_size, "%s/%s/order=%d/delta=%1.2f", results_path, m->name, order, m->delta);
            create_folders_in_path_if_not_exists(folder_name);

            char parameter_no_dot[256];
            remove_dot(parameter_no_dot, sizeof(parameter_no_dot), parameter);
            //For making a gith:
            //snprintf(filename3, sizeof(filename3), "%s/%d.png", folder_name, p);
            snprintf(filename3, sizeof(filename3), "%s/%s.png", folder_name, parameter_no_dot);

            plot_gnuplot_bulk_rps(filename3, N);
            
            dealloc_2d_int(&RM, N);
            dealloc_2d_double(&DM, N);

            print_prog((double)(p) / (double)number_of_parameters);
        }
    }

    return 0;
}

void ensemble_recurrence_analysis_all_ics(MAP *m)
{
    if (strcmp(m->name, "logistic") == 0)
    {
        double initial_conditions_set[][2] = {{0.1, M_PI}, {0.25, M_PI}, {0.50, M_PI}, {0.75, M_PI}};

        int num_initial_conditions = sizeof(initial_conditions_set) / sizeof(initial_conditions_set[0]);

        for (int i = 0; i < num_initial_conditions; i++)
        {
            double changing_ics[] = {initial_conditions_set[i][0], initial_conditions_set[i][1]};
            m->initial_conditions = changing_ics;
            printf("Running analysis for %d iterations and IC %f\n", N, m->initial_conditions[0]);
            ensemble_recurrence_analysis_varying_parameter(m);
        }

        m->delta = 5;
        printf("Running analysis for %d iterations and IC %f\n", N, initial_conditions[0]);
        ensemble_recurrence_analysis_random_parameter(m);
    }
    
    if (strcmp(m->name, "standard") == 0)
    {
        double delta_set[] = {2, 4, 6, 8};
        int num_deltas = sizeof(delta_set) / sizeof(delta_set[0]);

        for (int i = 0; i < num_deltas; i++)
        {
            m->delta = delta_set[i];
            double changing_ics[] = {m->delta * pow(10, -order), M_PI};
            m->initial_conditions = changing_ics;
            printf("Running analysis for %d iterations, delta = %1.2f, order = %d\n", N, m->delta, order);
            ensemble_recurrence_analysis_varying_parameter(m);
        }

        m->delta = 5;
        printf("Running analysis for %d iterations, delta = %1.2f, order = %d\n", N, m->delta, order);
        ensemble_recurrence_analysis_random_parameter(m);

    }

}

void map_details() 
{
    log_map.dimension = 1;
    log_map.parameter = 3.85;
    log_map.initial_conditions = initial_conditions;
    log_map.name = "logistic";
    log_map.distance_metric = "euclidean";

    std_map.dimension = 2;
    std_map.parameter = 1.80;
    std_map.delta = 1.0;
    double delta_ics[] = {std_map.delta * pow(10, -order), M_PI};
    //printf("%1.10f\n", delta_ics[0]);
    std_map.initial_conditions = delta_ics;
    //std_map.initial_conditions = initial_conditions;
    std_map.name = "standard";
    std_map.distance_metric = "modulated";
}

void plot_gnuplot_bifurcation_diagram(char png_name[], char *map_name)
{	
	FILE *gp;
	gp = popen(GNUPLOT, "w");
	fprintf(gp, "set terminal png enhanced size 2500,800 font \"CMR10\"  50 \n");
	fprintf(gp, "set colors classic \n");
	fprintf(gp, "set autoscale yfix \n");
	fprintf(gp, "set autoscale xfix \n");

    //Remove frame, axis and borders
/* 	fprintf(gp, "unset key \n");
	fprintf(gp, "unset cbtics \n");
	fprintf(gp, "unset xtics \n");
	fprintf(gp, "unset ytics \n");
	fprintf(gp, "unset colorbox \n");
	fprintf(gp, "unset border \n");
	fprintf(gp, "set lmargin 0 \n");
	fprintf(gp, "set rmargin 0 \n");
	fprintf(gp, "set tmargin 0 \n");
	fprintf(gp, "set bmargin 0 \n"); */

	fprintf(gp, "set xtics 3, 0.25, 4\n");

	fprintf(gp, "set output '%s'\n", png_name);

	fprintf(gp, "plot 'bifurcation_%s.dat' w d lc 'black' noti\n", map_name);
	fclose(gp);
}

void plot_gnuplot_bulk_rps(char filename[], int N)
{
    FILE *gp = popen(GNUPLOT, "w");
    if (gp == NULL)
    {
        fprintf(stderr, "Error opening pipe to gnuplot.\n");
        return;
    }
    fprintf(gp, "set terminal png size %d,%d \n", N, N);
    fprintf(gp, "set size square \n");
	fprintf(gp, "unset key \n");
	fprintf(gp, "unset xtics \n");
	fprintf(gp, "unset ytics \n");
	fprintf(gp, "unset cbtics \n");
	fprintf(gp, "unset colorbox \n");
	fprintf(gp, "unset border \n");
	fprintf(gp, "set lmargin 0 \n");
	fprintf(gp, "set rmargin 0 \n");
	fprintf(gp, "set tmargin 0 \n");
	fprintf(gp, "set bmargin 0 \n");
	//fprintf(gp, "set autoscale yfix \n");
	//fprintf(gp, "set autoscale xfix \n");
    fprintf(gp, "set xrange [-0.5:%d.5] \n", N-1);
    fprintf(gp, "set yrange [%d.5:-0.5] \n", N-1);
    fprintf(gp, "set palette maxcolors 2\n");
    fprintf(gp, "set palette defined (0 'white', 1 'black')\n");
    fprintf(gp, "set output '%s'\n", filename);
    fprintf(gp, "plot 'rm.dat' matrix w image\n");
    fclose(gp);
}

int main() 
{
    struct timeval start, end;
    double time_elapsed;
    srand(time(NULL));
    gettimeofday(&start, NULL);

    map_details();

    //Change the desired function for the main:
    //Make sure to set the desired map as well (&log_map) for the Logistic and (&std_map) for the Standard

    //give_bifurcation_diagram(&log_map, N);
    ensemble_recurrence_analysis_varying_parameter(&std_map);
    //ensemble_recurrence_analysis_random_parameter(&std_map);
    //ensemble_recurrence_analysis_all_ics(&log_map);

    gettimeofday(&end, NULL);
    time_elapsed = get_time_elapsed_minutes(start, end);

    printf("Run duration: %1.1f minutes\n", time_elapsed);
    return 0;
}