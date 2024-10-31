#ifndef REC_PLOTS_CHAOTIC_SYSTEMS_H
#define REC_PLOTS_CHAOTIC_SYSTEMS_H

struct map
{
    int dimension;
    double parameter;
    double *initial_conditions;
    double delta;
    char* name;
    const char* distance_metric;

};

typedef struct map MAP;

void standard_map_eqs(double *x0, double *x1, double *k);
void logistic_map_eqs(double *x0, double *x1, double *r);
void evolve_map(MAP *m, void *orbit, int N, double k);
void give_orbit(MAP *m, void *orbit, int N);
int give_recurrence_plot(int **R, int n);
int give_bifurcation_diagram(MAP *m, int N);
int ensemble_recurrence_analysis_varying_parameter(MAP *m);
int ensemble_recurrence_analysis_random_parameter(MAP *m);
void ensemble_recurrence_analysis_all_ics(MAP *m);
void map_details();
void plot_gnuplot_bifurcation_diagram(char png_name[], char *map_name);
void plot_gnuplot_bulk_rps(char filename[], int N);

#endif 