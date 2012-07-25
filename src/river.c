/* river.c*/
/* Simulation of pond or river segment with aquatic trophic chain*/
/* nutrient, phytoplankton, zooplankton, fish */
/* solar radiation and temperature have sinusoidal variation*/
/* flow is triangular */
/* this version does multiple runs*/
/* this version does endpoints*/

#include <stdio.h>
/*#include <stdlib.h>*/
#include <string.h>
#include <math.h>
#define N   11	       /* number of states */
#define DOSES	3     /* number of dosing episodes*/
#define DEBUG_SW  0 /*echo results to scrn when 1; change to 0 for release*/
#define BOUND 1e12 /*avoid overflow*/

/*function prototypes*/
/*rates*/
float per_year (float, float, float, float, float);	/*sinusoidal env factor*/
float triang (float, float, float, float, float); /* triangular flow function*/
float rad_ext (float, float);	 /*light extinction empirical function*/
float smith (float, float, float, float, float); /*primary productivity*/
float steele (float, float, float, float, float); /*primary productivity*/
float temp_mpy (float, float); /*temperature limiting factor*/
float monod (float, float);	/* monod M3 uptake, feeding or grazing*/

/*index*/
float EndMax (float,float);
float EndMin (float,float);
float EndAvg (float,float,float,float);
float EndRMS (float,float,float,float,float);

/*utilities*/
int in_as_round(float,float,float);

/* structures and union for model parameter*/
struct param_name{

/* 15 parameters for env and pond conditions*/

    float mean_temp;	 /* annual mean temp*/
    float amp_temp;	 /* annual temp excursion*/
    float base_temp;	 /* baseline temp */
    float del_temp;	 /* delay for temp wave*/

    float mean_rad;	 /* annual mean rad*/
    float amp_rad;	 /* annual rad excursion*/
    float del_rad;	 /* delay for rad wave*/

    /* make zero for closed pond */
    float base_flow;	 /* baseflow */
    float peak_flow;	 /* water flow peak*/
    float onset_flow;	 /* time for onset of high water flow*/
    float durat_flow;	 /* time for duration of high water flow wave*/

    float volume;	 /* water volume*/
    float depth;	 /* water depth */
    float eupho;	 /* euphotic depth */
    float ext_inorg;	 /* light extinction coefficient inorganics*/

/*3*DOSES + 6 parameters for nutrient discharge and advective input*/

    float nit_init;
    float pho_init;
    float upstream_nit;
    float upstream_pho;
    float discharge_nit;
    float discharge_pho;

    float input_nit[DOSES];	 /* input nitrogen*/
    float input_pho[DOSES];	 /* input phosphorous*/
    float time_nut[DOSES];	 /* time for nutrient input*/

/* 16 parameters for phyto plankton */

    float nit_carb;	 /* nitrogen / carbon ratio for phyto*/
    float pho_carb;	 /* phosphorous / carbon ration for phyto*/
    float max_uptake_nit;	/* maximum uptake rate */
    float max_uptake_pho;	/* maximum uptake rate */
    float half_uptake_nit;	/* half rate for nutrient uptake*/
    float half_uptake_pho;	/* half rate for nutrient uptake*/
    float alpha;	 /* efficiency of photosynthesis rate to light */
    float Lopt; 	 /* opt for Steele max photo rate */
    float Pmax; 	 /* maximum production rate */
    float max_grow_phyto; /* max growth rate for phytoplankton*/
    float carb_chlo;	 /* carbon to chlorophyll ratio*/
    float phyto_temp_mpy; /* coefficient for temp multiplier of pp rate*/
    float phyto_resp;	/* coefficient for respiration dependence on temp*/
    float phyto_death;	/* death rates*/
    float chlo_init; /*initial chlorophyll concentration*/
    float upstream_phyto;

 /*10 parameters for zooplankton*/

    float max_grazing;	/*maximum grazing rate*/
    float half_grazing;  /*half rate for grazing rate*/
    float zoo_eff;	/*efficiency of assimilated food*/
    float zoo_growth_max; /*maximum rate of secondary production*/
    float carb_indiv;	 /* carbon/L per individual */
    float zoo_temp_mpy; /* coeff for growth thermal dependence on temp*/
    float zoo_resp;	 /* coefficient for erspiration dependence on temp*/
    float zoo_death;	/*death rate zooplankton */
    float zoo_init;
    float upstream_zoop;   /* value at upstream segment */

/*9 parameters for fish*/

    float max_feeding;	/*maximum fish feeding rate*/
    float half_feeding;	/*half rate for fish feeding rate*/
    float fish_eff;	 /*efficiency of assimilated food*/
    float fish_growth_max; /*maximum rate of fish production*/
    float carb_length;	   /*ratio of carbon/L to a mm of length */
    float fish_temp_mpy; /* coeff for growth thermal dependence on temp*/
    float fish_resp;	 /* coefficient for respiration dependence on temp*/
    float fish_death;	 /* death rate of fishes */
    float fish_init;	 /* stocking of ponds */

};
struct param_array{
    float par[65];	 /*for control of loop run iteration*/
};
union param_name_array{
    struct param_name par_name;
    struct param_array par_array;
}param;

void river() {
/*variable declarations*/
float dt,time_final,time_zero, time_print_f;	/*simulation parameters*/
float x[11], time;	 /*model variables up to ten states*/
float net_rate[11]={0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00}, non_rate[11];
float par_sens_max,par_sens_min,par_sens_step;
                        /* max, min and step for parameter change*/
float x_max[11], x_min[11], x_avg[11], x_rms[11], x_end[11];
		/* indexes or endpoints */
float ru_ku_aux[11], ru_ku_tim, ru_ku_wt=1.00, ru_ku_avg[11], ru_ku_slope[11];
	/* runge kutta coefficients*/
float ext_coeff, eff_light;

int ru_ku_k;  /* runge kutta index */
unsigned int num_out_data;	/* length of simulation run as determined by sim par*/
unsigned int it_cnt, time_print, prn_cnt, par_cnt=0;
    /*iteration counter, number to print, print loop counter, parameter counter*/
int runs, num_runs;	/* current and max number of runs */
int num_states_print;	/* number of state variables for output file*/
int par_sens_flag=99;	   /*flag to determine parameter for sensitivity*/
int  i,j;   /* used as indexes for loops of states*/
unsigned int dose_cnt=0;
int steele_or_smith, lim_nit, lim_pho;

int error; /* to check fscanf */

char state_name[11][20];
	/* name of output variables for output file header*/
char title[100],sub_title[100], y_axis_title[100], x_axis_title[100];
 /* labels for graphs */
char file_in[30], file_out[30], file_dex[30], file_deb[30]; /* io files*/
char label[20]; /* for reading label of variable or parameter */
char par_sens_name[20];	 /*string for name of sensitivity parameter*/

/* strings for units */

char time_unit[40], x_unit[40];
char temp_unit[40], rad_unit[40], flow_unit[40], dis_unit[40];
char geom_unit[40], nut_unit[40], input_unit[40], coeff_unit[40];
char ratio_unit[40], bio_unit[40], rate_unit[40], prod_unit[40];
char phyto_unit[40], zoo_unit[40], fish_unit[40], eff_unit[40];
char alpha_unit[40], Lopt_unit[40], chlo_unit[40], zoop_unit[40];
char uptake_unit[40], grazing_unit[40], feeding_unit[40];

/* files */

FILE	*ifp,*ofp,*efp,*dfp;

/*get names of input, output and endpoints files */
strcpy(file_in, "river_inp.txt");
strcpy(file_out,"river_out.txt");
strcpy(file_dex,"river_dex.txt");
strcpy(file_deb,"river_deb.txt");

/*open files*/
if (DEBUG_SW == 1) dfp=fopen(file_deb,"w");
ifp=fopen(file_in,"r");
ofp=fopen(file_out,"w");
efp=fopen(file_dex,"w");


/*readin and verify parameters*/
error = fscanf(ifp," %s %s", label, par_sens_name);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, par_sens_name);
error = fscanf(ifp," %s %f", label, &par_sens_max);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, par_sens_max);
error = fscanf(ifp," %s %f", label,&par_sens_min);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, par_sens_min);
error = fscanf(ifp," %s %f", label,&par_sens_step);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, par_sens_step);
error = fscanf(ifp,"%s %[^\n]", label, sub_title);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, sub_title);
error = fscanf(ifp," %s %f",label, &time_zero);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, time_zero);
error = fscanf(ifp," %s %f",label, &time_final);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label,time_final);
error = fscanf(ifp," %s %f",label, &dt);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, dt);
error = fscanf(ifp," %s %d",label, &time_print);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %d\n", label, time_print);
error = fscanf(ifp," %s %s",label, time_unit);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, time_unit);


/*state parameters*/
num_states_print = N;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"num_states_print=%d\n", num_states_print);

/* read model parameters, set flag for sensitivity and verify */

error = fscanf(ifp," %s %f", label, &param.par_name.mean_temp);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.mean_temp);

error = fscanf(ifp," %s %f", label, &param.par_name.amp_temp);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.amp_temp);

error = fscanf(ifp," %s %f", label, &param.par_name.base_temp);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.base_temp);

error = fscanf(ifp," %s %s",label, temp_unit);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, temp_unit);

error = fscanf(ifp," %s %f", label, &param.par_name.del_temp);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.del_temp);

error = fscanf(ifp," %s %s",label, time_unit);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, time_unit);

error = fscanf(ifp," %s %f", label, &param.par_name.mean_rad);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.mean_rad);

error = fscanf(ifp," %s %f", label, &param.par_name.amp_rad);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.amp_rad);

error = fscanf(ifp," %s %s",label, rad_unit);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, rad_unit);

error = fscanf(ifp," %s %f", label, &param.par_name.del_rad);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.del_rad);

error = fscanf(ifp," %s %s",label, time_unit);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, time_unit);

error = fscanf(ifp," %s %f", label, &param.par_name.base_flow);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.base_flow);

error = fscanf(ifp," %s %f", label, &param.par_name.peak_flow);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.peak_flow);

error = fscanf(ifp," %s %s",label, flow_unit);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, flow_unit);

error = fscanf(ifp," %s %f", label, &param.par_name.onset_flow);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.onset_flow);

error = fscanf(ifp," %s %f", label, &param.par_name.durat_flow);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.durat_flow);

error = fscanf(ifp," %s %s",label, time_unit);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, time_unit);

error = fscanf(ifp," %s %f", label, &param.par_name.volume);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.volume);


error = fscanf(ifp," %s %f", label, &param.par_name.depth);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.depth);

error = fscanf(ifp," %s %f", label, &param.par_name.eupho);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.eupho);

error = fscanf(ifp," %s %f", label, &param.par_name.ext_inorg);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.ext_inorg);

error = fscanf(ifp," %s %s",label, geom_unit);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, geom_unit);


error = fscanf(ifp," %s %f", label, &param.par_name.nit_init);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.nit_init);

error = fscanf(ifp," %s %f", label, &param.par_name.pho_init);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.pho_init);

error = fscanf(ifp," %s %f", label, &param.par_name.upstream_nit);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.upstream_nit);

error = fscanf(ifp," %s %f", label, &param.par_name.upstream_pho);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.upstream_pho);


error = fscanf(ifp," %s %s",label, nut_unit);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, nut_unit);


error = fscanf(ifp," %s %f", label, &param.par_name.discharge_nit);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.discharge_nit);

error = fscanf(ifp," %s %f", label, &param.par_name.discharge_pho);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.discharge_pho);

error = fscanf(ifp," %s %s",label, dis_unit);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, dis_unit);

error = fscanf(ifp," %s %f", label, &param.par_name.input_nit[0]);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.input_nit[0]);

error = fscanf(ifp," %s %f", label, &param.par_name.input_nit[1]);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.input_nit[1]);

error = fscanf(ifp," %s %f", label, &param.par_name.input_nit[2]);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.input_nit[2]);

error = fscanf(ifp," %s %f", label, &param.par_name.input_pho[0]);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.input_pho[0]);

error = fscanf(ifp," %s %f", label, &param.par_name.input_pho[1]);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.input_pho[1]);

error = fscanf(ifp," %s %f", label, &param.par_name.input_pho[2]);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.input_pho[2]);

error = fscanf(ifp," %s %s",label, input_unit);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, input_unit);

error = fscanf(ifp," %s %f", label, &param.par_name.time_nut[0]);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.time_nut[0]);

error = fscanf(ifp," %s %f", label, &param.par_name.time_nut[1]);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.time_nut[1]);

error = fscanf(ifp," %s %f", label, &param.par_name.time_nut[2]);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.time_nut[2]);

error = fscanf(ifp," %s %s",label, time_unit);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, time_unit);


error = fscanf(ifp," %s %f", label, &param.par_name.nit_carb);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.nit_carb);

error = fscanf(ifp," %s %f", label, &param.par_name.pho_carb);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.pho_carb);

error = fscanf(ifp," %s %s",label, ratio_unit);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, ratio_unit);

error = fscanf(ifp," %s %s",label, x_unit);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, x_unit);

error = fscanf(ifp," %s %s",label, bio_unit);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, bio_unit);

error = fscanf(ifp," %s %f", label, &param.par_name.max_uptake_nit);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.max_uptake_nit);

error = fscanf(ifp," %s %d",label, &lim_nit);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %d\n", label, lim_nit);

error = fscanf(ifp," %s %f", label, &param.par_name.max_uptake_pho);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.max_uptake_pho);

error = fscanf(ifp," %s %d",label, &lim_pho);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %d\n", label, lim_pho);

error = fscanf(ifp," %s %s",label, uptake_unit);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, uptake_unit);

error = fscanf(ifp," %s %f", label, &param.par_name.half_uptake_nit);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.half_uptake_nit);

error = fscanf(ifp," %s %f", label, &param.par_name.half_uptake_pho);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.half_uptake_pho);

error = fscanf(ifp," %s %s",label, nut_unit);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, nut_unit);

error = fscanf(ifp," %s %f", label, &param.par_name.alpha);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.alpha);

error = fscanf(ifp," %s %s",label, alpha_unit);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, alpha_unit);

error = fscanf(ifp," %s %f", label, &param.par_name.Lopt);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.Lopt);

error = fscanf(ifp," %s %s",label, Lopt_unit);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, Lopt_unit);

error = fscanf(ifp," %s %d",label, &steele_or_smith);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %d\n", label, steele_or_smith);

error = fscanf(ifp," %s %f", label, &param.par_name.Pmax);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.Pmax);

error = fscanf(ifp," %s %s", label, prod_unit);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, prod_unit);

error = fscanf(ifp," %s %f", label, &param.par_name.max_grow_phyto);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.max_grow_phyto);

error = fscanf(ifp," %s %s", label, rate_unit);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, rate_unit);

error = fscanf(ifp," %s %f", label, &param.par_name.carb_chlo);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.carb_chlo);


error = fscanf(ifp," %s %s", label, ratio_unit);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, ratio_unit);

error = fscanf(ifp," %s %f", label, &param.par_name.phyto_temp_mpy);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.phyto_temp_mpy);

error = fscanf(ifp," %s %f", label, &param.par_name.phyto_resp);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.phyto_resp);

error = fscanf(ifp," %s %s",label, coeff_unit);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, coeff_unit);

error = fscanf(ifp," %s %f", label, &param.par_name.phyto_death);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.phyto_death);

error = fscanf(ifp," %s %s",label, rate_unit);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, rate_unit);

error = fscanf(ifp," %s %f", label, &param.par_name.chlo_init);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.chlo_init);

error = fscanf(ifp," %s %f", label, &param.par_name.upstream_phyto);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.upstream_phyto);

error = fscanf(ifp," %s %s",label, chlo_unit);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, chlo_unit);

error = fscanf(ifp," %s %f", label, &param.par_name.max_grazing);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.max_grazing);

error = fscanf(ifp," %s %s",label, grazing_unit);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, grazing_unit);

error = fscanf(ifp," %s %f", label, &param.par_name.half_grazing);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.half_grazing);

error = fscanf(ifp," %s %s",label, phyto_unit);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, phyto_unit);

error = fscanf(ifp," %s %f", label, &param.par_name.zoo_eff);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.zoo_eff);

error = fscanf(ifp," %s %s",label, eff_unit);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, eff_unit);

error = fscanf(ifp," %s %f", label, &param.par_name.zoo_growth_max);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.zoo_growth_max);

error = fscanf(ifp," %s %s",label, prod_unit);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, prod_unit);

error = fscanf(ifp," %s %f", label, &param.par_name.carb_indiv);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.carb_indiv);

error = fscanf(ifp," %s %s",label, ratio_unit);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, ratio_unit);

error = fscanf(ifp," %s %f", label, &param.par_name.zoo_temp_mpy);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.zoo_temp_mpy);

error = fscanf(ifp," %s %f", label, &param.par_name.zoo_resp);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.zoo_resp);

error = fscanf(ifp," %s %s",label, coeff_unit);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, coeff_unit);

error = fscanf(ifp," %s %f", label, &param.par_name.zoo_death);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.zoo_death);

error = fscanf(ifp," %s %s",label, rate_unit);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, rate_unit);

error = fscanf(ifp," %s %f", label, &param.par_name.zoo_init);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.zoo_init);

error = fscanf(ifp," %s %f", label, &param.par_name.upstream_zoop);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.upstream_zoop);

error = fscanf(ifp," %s %s",label, zoo_unit);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, zoo_unit);

error = fscanf(ifp," %s %f", label, &param.par_name.max_feeding);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.max_feeding);

error = fscanf(ifp," %s %s",label, feeding_unit);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, feeding_unit);

error = fscanf(ifp," %s %f", label, &param.par_name.half_feeding);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.half_feeding);

error = fscanf(ifp," %s %s",label, zoop_unit);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, zoop_unit);

error = fscanf(ifp," %s %f", label, &param.par_name.fish_eff);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.fish_eff);

error = fscanf(ifp," %s %s",label, eff_unit);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, eff_unit);

error = fscanf(ifp," %s %f", label, &param.par_name.fish_growth_max);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.fish_growth_max);

error = fscanf(ifp," %s %s",label, prod_unit);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, prod_unit);

error = fscanf(ifp," %s %f", label, &param.par_name.carb_length);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.carb_length);

error = fscanf(ifp," %s %s",label, ratio_unit);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, ratio_unit);

error = fscanf(ifp," %s %f", label, &param.par_name.fish_temp_mpy);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.fish_temp_mpy);

error = fscanf(ifp," %s %f", label, &param.par_name.fish_resp);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.fish_resp);

error = fscanf(ifp," %s %s",label, coeff_unit);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, coeff_unit);

error = fscanf(ifp," %s %f", label, &param.par_name.fish_death);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.fish_death);

error = fscanf(ifp," %s %s",label, rate_unit);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, rate_unit);

error = fscanf(ifp," %s %f", label, &param.par_name.fish_init);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.fish_init);

error = fscanf(ifp," %s %s",label, fish_unit);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, fish_unit);

if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"par_sens_flag=%d\n", par_sens_flag);

/* hard coded information about this simulation */
strcpy(title,"Aquatic Nutrient/Production");
if (DEBUG_SW == 1) fprintf(dfp,"title= %s\n", title);
strcpy(y_axis_title," ");
if (DEBUG_SW == 1) fprintf(dfp,"x_unit= %s\n",x_unit);
strcat(y_axis_title, x_unit);
strcpy(x_axis_title,"Time");
strcat(x_axis_title,time_unit);
if (DEBUG_SW == 1) fprintf(dfp,"y_axis_title=%s\nx_axis_title=%s\n",
   x_axis_title, y_axis_title);

strcpy(state_name[0],"Nitrog");
strcat(state_name[0], nut_unit);
if (DEBUG_SW == 1) fprintf(dfp,"state_name[0]=%s\n",state_name[0]);

strcpy(state_name[1],"Phosph");
strcat(state_name[1], nut_unit);
if (DEBUG_SW == 1) fprintf(dfp,"state_name[1]=%s\n",state_name[1]);

strcpy(state_name[2],"Chlo");
strcat(state_name[2], chlo_unit);
if (DEBUG_SW == 1) fprintf(dfp,"state_name[2]=%s\n",state_name[2]);

strcpy(state_name[3],"Zoop");
strcat(state_name[3], zoo_unit);
if (DEBUG_SW == 1) fprintf(dfp,"state_name[3]=%s\n",state_name[3]);

strcpy(state_name[4],"Fish");
strcat(state_name[4], fish_unit);
if (DEBUG_SW == 1) fprintf(dfp,"state_name[4]=%s\n",state_name[4]);

strcpy(state_name[5],"Phyto_C");
strcat(state_name[5], bio_unit);
if (DEBUG_SW == 1) fprintf(dfp,"state_name[5]=%s\n",state_name[5]);

strcpy(state_name[6],"Zoop_C");
strcat(state_name[6], bio_unit);
if (DEBUG_SW == 1) fprintf(dfp,"state_name[6]=%s\n",state_name[6]);

strcpy(state_name[7],"Fish_C");
strcat(state_name[7], bio_unit);
if (DEBUG_SW == 1) fprintf(dfp,"state_name[7]=%s\n",state_name[7]);

strcpy(state_name[8],"Temp");
strcat(state_name[8], temp_unit);
if (DEBUG_SW == 1) fprintf(dfp,"state_name[8]=%s\n",state_name[8]);

strcpy(state_name[9],"Rad");
strcat(state_name[9], rad_unit);
if (DEBUG_SW == 1) fprintf(dfp,"state_name[9]=%s\n",state_name[9]);

strcpy(state_name[10],"Flow");
strcat(state_name[10], flow_unit);
if (DEBUG_SW == 1) fprintf(dfp,"state_name[10]=%s\n",state_name[10]);



/*calculation of max number of runs */
num_runs=in_as_round(par_sens_max,par_sens_min,par_sens_step);

/* calculation of number of output data*/
time_print_f= ((float) time_print)*dt;
num_out_data= in_as_round(time_final,time_zero,time_print_f);

/* print header in output data file*/
fprintf(ofp," Number_of_runs: %d", num_runs);
fprintf(ofp,"\n Number_of_states: %d", num_states_print);
fprintf(ofp,"\n Number_of_data: %d", num_out_data);
fprintf(ofp,"\n Title: %s", title);
fprintf(ofp,"\n Sub_Title: %s", sub_title);
fprintf(ofp,"\n Y_axis_title_(default): %s", y_axis_title);
fprintf(ofp,"\n %s", x_axis_title);
for(i=0;i<=num_states_print-1;i++){
  fprintf(ofp," %s", state_name[i]);
}
fprintf(ofp,"\n");

/* print header in endpoint data file*/
fprintf(efp," Number_of_runs: %d", num_runs);
fprintf(efp,"\n Number_of_states: %d", num_states_print);
fprintf(efp,"\n Number_of_data: %d", num_out_data);
fprintf(efp,"\n Title: %s", title);
fprintf(efp,"\n Sub_Title: %s", sub_title);
fprintf(efp,"\n Y_axis_title_(default): %s", y_axis_title);
fprintf(efp,"\n Index");
for(i=0;i<=num_states_print-1;i++){
 fprintf(efp," %s", state_name[i]);
}
fprintf(efp,"\n");

/*-------simulation---------*/

/* initialize parameter for multiple runs*/
param.par_array.par[par_sens_flag]=par_sens_min;

/*----loop for multiple runs */
for (runs=0; runs<=num_runs-1; runs++){

    /* header to file_out*/
    fprintf(ofp," Simulation run number=%d", runs+1);
    fprintf(ofp,":   %s=%6.2e\n", par_sens_name, param.par_array.par[par_sens_flag]);

    /* header to endpoint file*/
    fprintf(efp," Simulation run number=%d", runs+1);
    fprintf(efp,":   %s=%6.2e\n", par_sens_name, param.par_array.par[par_sens_flag]);

    /* echo to scrn in debugging mode*/
    if (DEBUG_SW == 1)
      {
      fprintf(dfp," Simulating run number...%d\n", runs+1);	/* msg to scrn*/
      fprintf(dfp,"param.array[flag]=%6.2e\n",param.par_array.par[par_sens_flag]);
      fprintf(dfp," %s\n", x_axis_title);	 /*echo to scrn*/
      for(i=0;i<=N-1;i++){
	 fprintf(dfp," %s", state_name[i]);
	}
	 fprintf(dfp,"\n");
      }	/* end of if debug_sw */


    
    /*initialize simulation clock*/
    time=time_zero;

    /* scale factors */
    param.par_name.volume = param.par_name.volume * 1000.00;
    /*initialization of the dependent variables*/
      x[0]= param.par_name.nit_init;
      x[1]= param.par_name.pho_init;
      x[2]= param.par_name.chlo_init;
      x[3]= param.par_name.zoo_init;
      x[4]= param.par_name.fish_init;

      x[5]= param.par_name.chlo_init * param.par_name.carb_chlo;
      param.par_name.upstream_phyto = param.par_name.upstream_phyto * param.par_name.carb_chlo;
      x[6]= param.par_name.zoo_init * param.par_name.carb_indiv;
      param.par_name.upstream_zoop = param.par_name.upstream_zoop * param.par_name.carb_indiv;
      x[7]= param.par_name.fish_init * param.par_name.carb_length;
      x[8]= per_year (param.par_name.mean_temp,param.par_name.amp_temp,
		    param.par_name.base_temp,param.par_name.del_temp, time);
      x[9]= per_year (param.par_name.mean_rad,param.par_name.amp_rad,
		    0.00, param.par_name.del_rad, time);
      x[10]= triang (param.par_name.base_flow,param.par_name.peak_flow,
		    param.par_name.onset_flow, param.par_name.durat_flow, time);

    /*initialize index or endpoints*/
    for(i=0;i<=N-1;i++){
      x_max[i]=x[i];
      x_min[i]=x[i];
      x_avg[i]=x[i];
      x_rms[i]= (x[i] - x_avg[i]) * (x[i] - x_avg[i]);
      x_end[i]=x[i];
    } /* end of i loop*/

      /*first time to print result*/
      /*write to output file and echo to scrn if debug*/
      it_cnt = 0;
	fprintf(ofp," %6.2f",it_cnt*dt*time_print + time_zero);
if(DEBUG_SW == 1) fprintf(dfp,"%6.2f",it_cnt*dt*time_print + time_zero);
	    for(i=0;i<=N-1;i++){
if(DEBUG_SW == 1) fprintf(dfp," %10.2e",x[i]);
	     fprintf(ofp," %10.2e",x[i]);	  /*write to output file*/
    }  /* end of state loop for printing*/

	     fprintf(ofp,"\n");
	     if(DEBUG_SW == 1) fprintf(dfp,"\n");

    /* ---- loop for one run */
    for (it_cnt=1; it_cnt <= (num_out_data-1); it_cnt++){


      /*calculation loop until time to print*/
      for(prn_cnt=1;prn_cnt <= time_print; prn_cnt++){

      /* numerical integration; runge_kutta order 4 */

      /*initialize weighted average for each state */
      for(i=0;i<=N-1;i++) {ru_ku_avg[i] = 0.000;}

	/*echo to scrn when debugging*/
	if (DEBUG_SW == 1) printf(
	    "     k   ru_ku_tim   wt       aux    net_rate    slope	     avg\n");

      /* loop of 4 changing the order; ru_ku_k is index for ru_ku order */
      for (ru_ku_k=1; ru_ku_k<=4; ru_ku_k++){

	/* assigning weights to each order */
	if(ru_ku_k==1||ru_ku_k==4) ru_ku_wt =1.00;
	if(ru_ku_k==2||ru_ku_k==3) ru_ku_wt =2.00;

	/* updating state and time */
	if(ru_ku_k==1){
	    for(i=0;i<=N-1;i++){ru_ku_aux[i] = x[i];}
	    ru_ku_tim = time;
	}
	else {
	    for(i=0;i<=N-1;i++){
		ru_ku_aux[i] = x[i] + (dt/ru_ku_wt)*net_rate[i];
		if(ru_ku_aux[i] < 0.00) ru_ku_aux[i] = 0.00;
	    }
	    ru_ku_tim = time + dt/ru_ku_wt;
	}

   /*-----------------------------*/
   /*configuring the rate equations*/

	  /*chlorophyll, zooplankton individuals and fish length
	  are observation variables, then the rates are zero.
	  They will be calculated scaling the state variables*/
	  net_rate[2]=0.00;
	  net_rate[3]=0.00;
	  net_rate[4]=0.00;

	  /*calculation of phytoplankton production*/
	  /* calculate extinction coeff */
	  ext_coeff   =	rad_ext (ru_ku_aux[2],param.par_name.ext_inorg);

	  if (steele_or_smith == 0)
	  {
	  /* photosynthesis rate using depth averaged Smith's model*/
	  eff_light   = smith (ru_ku_aux[9],param.par_name.alpha,
			ext_coeff, param.par_name.eupho, ru_ku_aux[2]);
	  net_rate[5] = param.par_name.Pmax * eff_light;
	  }
	  else
	  {
	  /* photosynthesis rate using depth averaged Steele's model*/
	  eff_light   = steele (ru_ku_aux[9],param.par_name.Lopt,
			ext_coeff, param.par_name.eupho, ru_ku_aux[2]);
	  net_rate[5] = param.par_name.Pmax * eff_light;
	  }

	  /* nutrient uptake multiplier */
	  if (lim_nit == 1) net_rate[5] = net_rate[5] * monod(ru_ku_aux[0], param.par_name.half_uptake_nit);
	  if (lim_pho == 1) net_rate[5] = net_rate[5] * monod(ru_ku_aux[1], param.par_name.half_uptake_pho);

	  /* temperature multiplier */
	  net_rate[5] = net_rate[5] * temp_mpy (param.par_name.phyto_temp_mpy,ru_ku_aux[8]);

	  /*calculating losses of nutrients due to uptake by algae*/

	  /*loss of nitrogen due to uptake by algae*/
	  net_rate[0]= - param.par_name.max_uptake_nit * net_rate[5] * ru_ku_aux[5] * param.par_name.nit_carb;
	  /* input of nitrogen due to continuous discharge */
	  net_rate[0] = net_rate[0] + param.par_name.discharge_nit/param.par_name.volume;
	  /* input and loss of nitrogen due to advection by bulk flow*/
	  net_rate[0] = net_rate[0] + 86400.00 * ru_ku_aux[10] *
	     (param.par_name.upstream_nit - ru_ku_aux[0])/param.par_name.volume;

	  /*loss of phosphorous due to uptake by algae*/
	  net_rate[1]= - param.par_name.max_uptake_pho * net_rate[5] * ru_ku_aux[5] * param.par_name.pho_carb;
	  /* input of phosphorous due to continuous discharge */
	  net_rate[1] = net_rate[1] + param.par_name.discharge_pho/param.par_name.volume;
	  /* input and loss of phosphorous due to advection by bulk flow*/
	  net_rate[1] = net_rate[1] + 86400.00 * ru_ku_aux[10]
	    *(param.par_name.upstream_pho - ru_ku_aux[1])/param.par_name.volume;

	  /* subtract phytoplankton respiration losses */
	  net_rate[5] =	net_rate[5] - param.par_name.phyto_resp * ru_ku_aux[8];

	  /* invest net energy in reproduction or population growth*/
	  
	  net_rate[5] = net_rate[5] * param.par_name.max_grow_phyto * ru_ku_aux[5];

	  /*calculation of zooplankton production*/
	  /* food intake */
	  net_rate[6] = param.par_name.max_grazing *
			monod(ru_ku_aux[5], param.par_name.half_grazing);

	  /* temperature multiplier */
	  net_rate[6] = net_rate[6] * temp_mpy (param.par_name.zoo_temp_mpy,
						ru_ku_aux[8]);

	  /* subtract phyto population losses */
	  net_rate[5] = net_rate[5] - param.par_name.phyto_death * ru_ku_aux[5] -
			net_rate[6] * ru_ku_aux[6];

	  /* add net phytoplankton advection by flow */
	  net_rate[5] = net_rate[5] + 86400.00 * ru_ku_aux[10]*
	    (param.par_name.upstream_phyto - ru_ku_aux[5])/param.par_name.volume;

	  /* zooplankton assimilation efficiency*/
	  net_rate[6] = net_rate[6] * param.par_name.zoo_eff;

	  /* subtract zooplankton respiration losses */
	  net_rate[6] =	net_rate[6] - param.par_name.zoo_resp * ru_ku_aux[8];

	  /* zooplankton invest in reproduction, population growth */
	  net_rate[6] = net_rate[6] * param.par_name.zoo_growth_max * ru_ku_aux[6];

	  /*calculation of fish production*/
	  /* food intake */
	  net_rate[7] = param.par_name.max_feeding *
			monod(ru_ku_aux[6], param.par_name.half_feeding);

	  /* temperature multiplier */
	  net_rate[7] = net_rate[7] * temp_mpy (param.par_name.fish_temp_mpy,
							ru_ku_aux[8]);
	  /* subtract zooplankton population losses */
	  net_rate[6] =	net_rate[6] - param.par_name.zoo_death * ru_ku_aux[6];
	  net_rate[6] = net_rate[6] - net_rate[7]* ru_ku_aux[7];

	  /* add net zooplankton advection by flow */
	  net_rate[6] = net_rate[6] + 86400.00 * ru_ku_aux[10]*
	       (param.par_name.upstream_zoop - ru_ku_aux[6])/param.par_name.volume;

	  /* fish assimilation efficiency*/
	  net_rate[7] = net_rate[7] * param.par_name.fish_eff;

	  /* subtract fish respiration losses */
	  net_rate[7] =	net_rate[7] - param.par_name.fish_resp * ru_ku_aux[8];

	  /* fish invest in growth, no reproduction */
	  net_rate[7] = net_rate[7] * param.par_name.fish_growth_max * ru_ku_aux[7];

	  /* subtract fish population losses */
	  net_rate[7] =	net_rate[7] - param.par_name.fish_death * ru_ku_aux[7];

	  /* temperature and radiation, no rates*/
	  net_rate[8]=0.00;
	  net_rate[9]=0.00;
	  net_rate[10]=0.00;

      /*calculating slope k */
      for(i=0;i<=N-1;i++){
      ru_ku_slope[i] = dt * net_rate[i];

      /* updating weighted average*/
      ru_ku_avg[i]+= ru_ku_wt * ru_ku_slope[i];

      /*echo to screen when debugging*/
      if (DEBUG_SW == 1) fprintf(dfp,"%6d  %7.4f %7.4f %8.4e %8.4e %8.4e %8.4e\n",
	ru_ku_k, ru_ku_tim, ru_ku_wt, ru_ku_aux[i], net_rate[i],
	ru_ku_slope[i], ru_ku_avg[i]);

      }	     /* end of i loop */

      }    /* end of ru_ku_k loop; runge kutta order */
	
   /* ----------------------------*/
   /*impulsive disturbance if any*/

   /* impulsive is dosing of nutrients, first two compartments*/

     if (dose_cnt < DOSES){
	  if (time < param.par_name.time_nut[dose_cnt]){
	   non_rate[0]=0.00;
	   non_rate[1]=0.00;
	   }		/* dosing does not occur*/

	  else	   {	/* dosing occurs*/
	   non_rate[0]= param.par_name.input_nit[dose_cnt];
	   non_rate[1]= param.par_name.input_pho[dose_cnt];
	   if (dose_cnt == DOSES) dose_cnt = DOSES;
	   else ++dose_cnt;			 /* counter dosing*/
		   }
     }
     else {
	   non_rate[0]=0.00;
	   non_rate[1]=0.00;
	  }

       /* zero impulsive for the rest*/
	for(i=2;i<=N-1;i++){
	 non_rate[i] = 0.00;
	}

     for(i=0;i<=N-1;i++){
    /*updating state with final average of the four slopes*/
    x[i]+= (1.00/6.00)* ru_ku_avg[i] + non_rate[i];


    if (x[i] < 0.00) x[i]=0.00; /* force to zero if negative*/
    if (x[i] > BOUND) x[i] = BOUND; /* avoid overflow*/

    }	 /*end of state loop*/

      /* calculation of env and observation variables*/
      x[8]=  per_year (param.par_name.mean_temp, param.par_name.amp_temp,
		    param.par_name.base_temp,param.par_name.del_temp, time);
      x[9] = per_year (param.par_name.mean_rad, param.par_name.amp_rad,
		    0.00, param.par_name.del_rad, time);
      x[10]= triang (param.par_name.base_flow,param.par_name.peak_flow,
		    param.par_name.onset_flow, param.par_name.durat_flow, time);

      x[2]=  x[5] / param.par_name.carb_chlo;
      x[3]=  x[6] / param.par_name.carb_indiv;
      x[4]=  x[7] / param.par_name.carb_length;

   /* end of numerical integration runge kutta order 4*/
   /* ----------------------------*/

    if (DEBUG_SW == 1){
	for (j=0;j<=9;j++){
	 fprintf(dfp,"time = %4.2f j= %2d x[j]= %6.2f net_rate[j]= %6.2f non_rate[j]= %6.2f\n",
	 time, j, x[j], net_rate[j], non_rate[j]);
	}
    }

      /*calculate temporary index or endpoints: max, min, mean, rms, end*/
    for(i=0;i<=N-1;i++){
	x_max[i]=EndMax(x[i],x_max[i]);
	x_min[i]=EndMin(x[i],x_min[i]);
	x_avg[i]=EndAvg(x[i],x_avg[i],dt,time-time_zero);
	x_rms[i]=EndRMS(x[i],x_rms[i],x_avg[i],dt,time-time_zero);
	x_end[i]=x[i];
    } /* end of i loop */

	      if(DEBUG_SW == 1) fprintf(dfp," end of i loop");
	      
      time+=dt;	  /*inc time*/

    }  /* end of calculation loop before printing result*/

      /*time to print result*/
      /*write to output file and echo to scrn if debug*/

	fprintf(ofp," %6.2f",it_cnt*dt*time_print + time_zero);
if(DEBUG_SW == 1) fprintf(dfp,"%6.2f",it_cnt*dt*time_print + time_zero);
	    for(i=0;i<=N-1;i++){
if(DEBUG_SW == 1) fprintf(dfp," %10.2e",x[i]);
	     fprintf(ofp," %10.2e",x[i]);	  /*write to output file*/
    }  /* end of state loop for printing*/

	     fprintf(ofp,"\n");
	     if(DEBUG_SW == 1) fprintf(dfp,"\n");

    }  /* ----- end of simulation run loop */

    /*finalize index*/
    for(i=0;i<=N-1;i++){
      x_rms[i]= sqrt(x_rms[i]);
    }

    /* write to index or endpoint file */

	fprintf(efp," max");
if(DEBUG_SW == 1) fprintf(dfp," max");
	    for(i=0;i<=N-1;i++){
if(DEBUG_SW == 1) fprintf(dfp," %10.2e",x_max[i]);
	     fprintf(efp," %10.2e",x_max[i]);	  /*write to endpoint file*/
    }  /* end of state loop for printing*/
	     fprintf(efp,"\n");

	fprintf(efp," min");
if(DEBUG_SW == 1) fprintf(dfp," min");
	    for(i=0;i<=N-1;i++){
if(DEBUG_SW == 1) fprintf(dfp," %10.2e",x_min[i]);
	     fprintf(efp," %10.2e",x_min[i]);	  /*write to endpoint file*/
    }  /* end of state loop for printing*/
	     fprintf(efp,"\n");

	fprintf(efp," avg");
if(DEBUG_SW == 1) fprintf(dfp," avg");
	    for(i=0;i<=N-1;i++){
if(DEBUG_SW == 1) fprintf(dfp," %10.2e",x_avg[i]);
	     fprintf(efp," %10.2e",x_avg[i]);	  /*write to endpoint file*/
    }  /* end of state loop for printing*/
	     fprintf(efp,"\n");

	fprintf(efp," rms");
if(DEBUG_SW == 1) fprintf(dfp," rms");
	    for(i=0;i<=N-1;i++){
if(DEBUG_SW == 1) fprintf(dfp," %10.2e",x_rms[i]);
	     fprintf(efp," %10.2e",x_rms[i]);	  /*write to endpoint file*/
    }  /* end of state loop for printing*/
	     fprintf(efp,"\n");

	fprintf(efp," end");
if(DEBUG_SW == 1) fprintf(dfp," end");
	    for(i=0;i<=N-1;i++){
if(DEBUG_SW == 1) fprintf(dfp," %10.2e",x_end[i]);
	     fprintf(efp," %10.2e",x_end[i]);	  /*write to endpoint file*/
    }  /* end of state loop for printing*/
	     fprintf(efp,"\n");

    /*change value of parameter for next run*/
    param.par_array.par[par_sens_flag]+=par_sens_step;

}    /* -------end of run iteration-------*/

/*close files*/
fclose(ifp);
fclose(ofp);
fclose(efp);
if (DEBUG_SW == 1) fclose(dfp);

}    /* --------end of main program ------ */


/*function to evaluate M3 hyperbolic uptake, grazing and feeding rate*/
float monod (float resource, float half_rate)
{
float rate;
rate=resource/(resource+half_rate);
return rate;
}

/*function to evaluate temperature multiplier*/
float temp_mpy (float temp, float coeff)
{
float rate;
rate = coeff*temp;
return rate;
}

/*function to evaluate periodic annual radiation*/
float per_year(float mean, float amplitude, float base, float delay, float time)
{
float rate;
    rate= mean + amplitude*cos((6.28/365.00)*(time-delay));
		/*sinusoidal annual rate*/
    if (rate <= base) rate = base;
		/* clip if zero */
return rate;
}

/*function to evaluate triangular hydrograph*/
float triang(float base, float peak, float onset, float durat, float time)
{
float rate;
    if ((onset <= time) & (time <= onset + durat)) {
	if ((onset <= time) & (time <= onset + durat/2.00))
	rate= base + ((peak - base)/(durat/2.00)) * (time - onset);
	else
	rate= peak - ((peak - base)/(durat/2.00)) * (time - onset - durat/2.00);
	}
    else rate = base;

return rate;
}


/*function to evaluate algae plus inorganic exticntion coeff*/
float rad_ext (float chlo, float kext)
{
float rate;
    rate =  kext + 0.0088 * chlo + 0.053 * pow(chlo, 0.66);
		/*empirical extinction due to chlorophyll*/
    if (rate < 0.00) rate = 0.00;
		/* clip if zero */
return rate;
}

/*functions to evaluate rad multiplier*/
float smith (float rad, float alpha, float ext, float depth, float chlo)
{
float rate = 0;
float root_num, root_den, rad_down;
float z, y;
      z = ext * depth;
      rad_down = rad * exp ( -z );

      if (alpha == 0.00) y = 0.00;
      else  y = chlo / alpha;

      root_num = sqrt ( y * y + rad * rad );
      root_den = sqrt ( y * y + rad_down * rad_down );

      if ((rad_down + root_den) != 0.00)
       rate = (rad + root_num) / (rad_down + root_den);

      if (rate != 0.00 ) rate = log (rate);

      if(z == 0.00)    rate = 0.00;
      else	       rate= rate * chlo / z ;

return rate;
}

float steele (float rad, float Lopt, float ext, float depth, float chlo)
{
float rate;
float z, ratio;

      z = ext * depth;
      ratio = rad / Lopt;
      rate = exp (-ratio * exp(-z) ) - exp(-ratio);
      rate = rate * exp(1) / z;

return rate;
}


/*function to calculate integer */
int in_as_round(float max_val, float min_val, float step_val)
{
float val;
int n_val;
if (step_val == 0.00) val = 0.00;
  else val= ((max_val - min_val) / step_val);
if (val - floor(val) >= 0.5)
    n_val= ceil(val)+1;
  else
    n_val= floor(val)+1;
return n_val;
}

/*function to evaluate max index or endpoint*/
float EndMax (float x_val, float val_max)
{
float maximum;
if (x_val >= val_max) maximum=x_val;
    else maximum=val_max;
return maximum;
}

/*function to evaluate min index or endpoint */
float EndMin (float x_val, float val_min)
{
float minimum;
if (x_val <= val_min) minimum=x_val;
    else minimum=val_min;
return minimum;
}

/*function to evaluate mean index or endpoint */
float EndAvg (float x_val, float val_mean, float delta_t,
	       float previous_time)
{
float mean;
       mean = (val_mean * previous_time/delta_t + x_val) *
       delta_t/(delta_t + previous_time);
return mean;
}

/*function to evaluate rms index or endpoint */
float EndRMS (float x_val, float val_rms, float val_mean, float delta_t,
	       float previous_time)
{
float rms;
rms = (x_val - val_mean) * (x_val - val_mean);
       rms = (val_rms * previous_time/delta_t + rms) *
       delta_t/(delta_t + previous_time);
return rms;
}
