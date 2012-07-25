/*forsucc.c*/
/*Simulation of forest succession, n species*/
/*Deterministic version of jabowa from Swartzman and Kaluzny page 129*/
/*four species cherry, birch, maple, spruce*/
/*forsucc= forest succession */
/*this version does multiple runs & endpoints*/

#include <stdio.h>
#include <string.h>
#include <math.h>

#define DEBUG_SW 0	/* echo input to scrn when 1; change to 0 for release*/
#define BOUND 1e12	/* to avoid overflow*/
#define MAX 4		/* maximum number of stages*/

/*function prototypes*/

float height_dbh(float, float, float);
float basal(float, float);
float res_light(float, float, float, float);
float res_temp(float, float, float);
float dh_deno(float, float, float);
float mort_rate(float, float, float, float, float);

/* index */
float end_max(float,float);
float end_min(float,float);
float end_avg(float,float,float,float);
float end_rms(float,float,float,float,float);

/* utilities */
int in_as_round1(float,float,float);

/* structures and union for model parameters */
struct param_name{
/*initial 3*MAX*/
    float dbh0[MAX];	    /*initial conditions diameter*/
    float num0[MAX];	    /*initial conditions number of trees*/
    float age0[MAX];	    /*initial condtions age of trees*/
/*species data	10*MAX */
    float dbh_max[MAX];	    /*max dbh */
    float height_max[MAX];  /*max height*/
    float k1_light[MAX];    /*shade tolerance coefficients*/
    float k2_light[MAX];
    float k3_light[MAX];
    float dd_min[MAX];
    float dd_max[MAX];
    float gmax[MAX];
    float la_dbh[MAX];
    float k_mort[MAX];
/*stand and stress data 6 */
    float zlight;
    float kext;
    float temp;
    float basmax;
    float k_slow;
    float thresh_slow;
/*disturbance parameters MAX + 1  */
    float period;	   /*period for disturbance events*/
    float intensity[MAX];   /*intensity for disturbance events*/

    };
struct param_array{
    float par[3*MAX + 10*MAX +	6 + MAX + 1 ];
     /*for control of run loop iteration*/
    };
union param_name_array{
    struct param_name par_name;
    struct param_array par_array;
    }param;

void forsucc() {
/* variable declarations */
float dt,time_final,time_zero, time_print_f;	/*simulation parameters*/
float x[2*MAX], x_prop[2*MAX], time, sum_x;	 /*model variables up to ten states*/
float net_rate[2*MAX], non_rate[2*MAX];
     /*to configure rate eqn for integration*/
/*variable declarations specific to jabowa like models*/
float coeff1[MAX], coeff2[MAX], dh_max[MAX];
     /*allometric coefficients and produtc of d and h max*/
float dbh[MAX], age[MAX], growth[MAX];
    /* dbh, basal area, age, growth rate*/
float height[MAX], mort[MAX], mortage[MAX], deno[MAX], fract[MAX];
    /*height, mortality, denominator of d and h, product of d and h*/
float mult_light[MAX], mult_temp[MAX], mult_space;
/*multipliers*/
float above_leaf[MAX], light[MAX], leaf_area[MAX];
/* leaf area above species i, light received by i, leaf area of i */

float par_sens_max,par_sens_min,par_sens_step;
                        /* max, min and step for parameter change*/
float x_max[2*MAX], x_min[2*MAX], x_avg[2*MAX], x_rms[2*MAX],
      x_end[2*MAX];
/* specific to jabowa */
float num_total, bas_total;
float num_total_max, num_total_min, num_total_avg, num_total_rms, num_total_end;
float bas_total_max, bas_total_min, bas_total_avg, bas_total_rms, bas_total_end;

float x_prop_max[2*MAX], x_prop_min[2*MAX], x_prop_avg[2*MAX],
	 x_prop_rms[2*MAX], x_prop_end[2*MAX];
		/* indexes or endpoints */
float sum_non_rate;  /* accumulator for disturbance to reset succession*/
float dis_cnt=1.0;   /*counter of disturbance events*/

unsigned int num_out_data;/* length of simulation run as determined by sim par*/
unsigned int it_cnt, time_print, prn_cnt, par_cnt=0;
  /*iteration counter, number to print, pnr loop counter and model parameter counter*/
int runs, num_runs;	/* current and max number of runs */
int num_states_print;	/* number of state variables for output file*/
int par_sens_flag=999;	   /*flag to determine parameter for sensitivity*/
int  i,j;   /* used as indexes for loops of states*/
int N;	/* number of compartments */
int display_prop; /* switch to display absolute or proportions*/
int prop_sw;	    /*disturbance type: prop (1), non_prop (0) */
int error; /* to check fscanf */

char state_name[20][20], c_st_name[2];
	/* name of output variables for output file header*/
char title[100],sub_title[100], y_axis_title[100], x_axis_title[100];
 /* labels for graphs */
char file_in[30], file_out[30], file_dex[30], file_deb[30]; /* io files*/
char label[20]; /* for reading label of variable or parameter */
char par_sens_name[20];	 /*string for name of sensitivity parameter*/
char time_unit[40],  x_unit[40],
      x_prop_unit[40];	/*units*/
char period_unit[40],intensity_unit[40];   /* more strings for units*/
char num_unit[40], age_unit[40], k_light_unit[40], dd_unit[40],
     g_unit[40], la_dbh_unit[40], k_mort_unit[40], zlight_unit[40],
     kext_unit[40], temp_unit[40], bas_unit[40], k_slow_unit[40],
     thresh_slow_unit[40];
     /* these are units specific to JABOWA */

FILE	*ifp,*ofp,*efp,*dfp;


/*get names of input and output files */
strcpy(file_in, "forsucc_inp.txt");
strcpy(file_out,"forsucc_out.txt");
strcpy(file_dex,"forsucc_dex.txt");
strcpy(file_deb,"forsucc_deb.txt");

/*open files*/
ifp=fopen(file_in,"r");
ofp=fopen(file_out,"w");
efp=fopen(file_dex,"w");
if(DEBUG_SW == 1) dfp=fopen(file_deb,"w");

/*readin and verify parameters*/
error = fscanf(ifp," %s %s", label, par_sens_name);
if(error == 1 && DEBUG_SW == 1)  fprintf(dfp,"%s %s\n", label, par_sens_name);
error = fscanf(ifp," %s %f", label, &par_sens_max);
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, par_sens_max);
error = fscanf(ifp," %s %f", label,&par_sens_min);
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, par_sens_min);
error = fscanf(ifp," %s %f", label,&par_sens_step);
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, par_sens_step);
error = fscanf(ifp,"%s %[^\n]", label, sub_title);
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, sub_title);
error = fscanf(ifp," %s %f",label, &time_zero);
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, time_zero);
error = fscanf(ifp," %s %f",label, &time_final);
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label,time_final);
error = fscanf(ifp," %s %f",label, &dt);
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, dt);
error = fscanf(ifp," %s %d",label, &time_print);
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %d\n", label, time_print);
error = fscanf(ifp," %s %s", label, time_unit);
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, time_unit);
error = fscanf(ifp," %s %d",label, &N);
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %d\n", label, N);
error = fscanf(ifp," %s %d",label, &display_prop);
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %d\n", label, display_prop);

/*read and verify mode switch parameter*/
error = fscanf(ifp," %s %d", label, &prop_sw);
if (error == 1 && DEBUG_SW == 1) fprintf(dfp,"prop_sw=%d\n", prop_sw);

error = fscanf(ifp," %s %s", label, x_prop_unit);
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, x_prop_unit);

num_states_print= 2*N+ 2;
if(DEBUG_SW == 1) fprintf(dfp,"num_states_print=%d\n", num_states_print);

/* read model parameters, set flag for sensitivity and verify */

for(i=0;i<=N-1;i++){
 error = fscanf(ifp," %s %f", label, &param.par_name.dbh0[i]);
    if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
    par_cnt++;
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.dbh0[i]);
}
par_cnt+=MAX-N;

for(i=0;i<=N-1;i++){
 error = fscanf(ifp," %s %f", label, &param.par_name.dbh_max[i]);
    if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
    par_cnt++;
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.dbh_max[i]);
}
par_cnt+=MAX-N;

for(i=0;i<=N-1;i++){
 error = fscanf(ifp," %s %f", label, &param.par_name.height_max[i]);
    if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
    par_cnt++;
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.height_max[i]);
}
par_cnt+=MAX-N;

error = fscanf(ifp," %s %s", label, x_unit);
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, x_unit);

for(i=0;i<=N-1;i++){
 error = fscanf(ifp," %s %f", label, &param.par_name.num0[i]);
    if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
    par_cnt++;
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.num0[i]);
}
par_cnt+=MAX-N;

error = fscanf(ifp," %s %s", label, num_unit);
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, num_unit);

for(i=0;i<=N-1;i++){
 error = fscanf(ifp," %s %f", label, &param.par_name.age0[i]);
    if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
    par_cnt++;
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.age0[i]);
}
par_cnt+=MAX-N;

error = fscanf(ifp," %s %s", label, age_unit);
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, age_unit);


for(i=0;i<=N-1;i++){
 error = fscanf(ifp," %s %f", label, &param.par_name.k1_light[i]);
    if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
    par_cnt++;
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.k1_light[i]);
}
par_cnt+=MAX-N;

for(i=0;i<=N-1;i++){
 error = fscanf(ifp," %s %f", label, &param.par_name.k2_light[i]);
    if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
    par_cnt++;
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.k2_light[i]);
}
par_cnt+=MAX-N;

for(i=0;i<=N-1;i++){
 error = fscanf(ifp," %s %f", label, &param.par_name.k3_light[i]);
    if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
    par_cnt++;
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.k3_light[i]);
}
par_cnt+=MAX-N;

error = fscanf(ifp," %s %s", label, k_light_unit);
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, k_light_unit);


for(i=0;i<=N-1;i++){
 error = fscanf(ifp," %s %f", label, &param.par_name.dd_min[i]);
    if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
    par_cnt++;
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.dd_min[i]);
}
par_cnt+=MAX-N;

for(i=0;i<=N-1;i++){
 error = fscanf(ifp," %s %f", label, &param.par_name.dd_max[i]);
    if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
    par_cnt++;
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.dd_max[i]);
}
par_cnt+=MAX-N;

error = fscanf(ifp," %s %s", label, dd_unit);
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, dd_unit);


for(i=0;i<=N-1;i++){
 error = fscanf(ifp," %s %f", label, &param.par_name.gmax[i]);
    if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
    par_cnt++;
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.gmax[i]);
}
par_cnt+=MAX-N;

error = fscanf(ifp," %s %s", label, g_unit);
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, g_unit);

for(i=0;i<=N-1;i++){
 error = fscanf(ifp," %s %f", label, &param.par_name.la_dbh[i]);
    if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
    par_cnt++;
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.la_dbh[i]);
}
par_cnt+=MAX-N;

error = fscanf(ifp," %s %s", label, la_dbh_unit);
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, la_dbh_unit);


for(i=0;i<=N-1;i++){
 error = fscanf(ifp," %s %f", label, &param.par_name.k_mort[i]);
    if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
    par_cnt++;
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %8.6f\n", label, param.par_name.k_mort[i]);
}
par_cnt+=MAX-N;

error = fscanf(ifp," %s %s", label, k_mort_unit);
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, k_mort_unit);

error = fscanf(ifp," %s %f", label, &param.par_name.zlight);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.zlight);

error = fscanf(ifp," %s %s", label, zlight_unit);
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, zlight_unit);

error = fscanf(ifp," %s %f", label, &param.par_name.kext);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %8.6f\n", label, param.par_name.kext);

error = fscanf(ifp," %s %s", label, kext_unit);
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, kext_unit);


error = fscanf(ifp," %s %f", label, &param.par_name.temp);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.temp);

error = fscanf(ifp," %s %s", label, temp_unit);
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, temp_unit);

error = fscanf(ifp," %s %f", label, &param.par_name.basmax);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.basmax);

error = fscanf(ifp," %s %s", label, bas_unit);
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, bas_unit);

error = fscanf(ifp," %s %f", label, &param.par_name.k_slow);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.k_slow);

error = fscanf(ifp," %s %s", label, k_slow_unit);
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, k_slow_unit);

error = fscanf(ifp," %s %f", label, &param.par_name.thresh_slow);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.thresh_slow);

error = fscanf(ifp," %s %s", label, thresh_slow_unit);
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, thresh_slow_unit);


error = fscanf(ifp," %s %f", label, &param.par_name.period);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.period);

error = fscanf(ifp," %s %s",label, period_unit);
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, period_unit);

for(i=0;i<=N-1;i++){
error = fscanf(ifp," %s %f", label, &param.par_name.intensity[i]);
if(strcasecmp(label,par_sens_name) == 0) par_sens_flag=par_cnt;
par_cnt++;
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %6.2f\n", label, param.par_name.intensity[i]);
}
par_cnt+= MAX-N;

error = fscanf(ifp," %s %s",label, intensity_unit);
if(error == 1 && DEBUG_SW == 1) fprintf(dfp,"%s %s\n", label, intensity_unit);


/* hard coded information about this simulation */
strcpy(title,"Forest Succession");
if(DEBUG_SW == 1) fprintf(dfp,"title= %s\n", title);
strcpy(y_axis_title,"BA__");

strcat(y_axis_title,x_unit);
strcat(y_axis_title,x_prop_unit);

strcpy(x_axis_title,"Time");
strcat(x_axis_title,time_unit);
if(DEBUG_SW == 1)
fprintf(dfp,"y_axis_title=%s\nx_axis_title=%s\n",x_axis_title, y_axis_title);


/*configuration of state names */
for (i=0;i<=num_states_print - 1;i++){
if (i <	N) {
    if (display_prop==0) strcpy(state_name[i], "BA");
    else strcpy(state_name[i], "RelBA");
    if(i<=9){
	c_st_name[0]=i+49;
	c_st_name[1]='\0';
    }
    else strcpy(c_st_name,"0");
      /* ascii code of the state number assuming base 48*/
    strcat(state_name[i],c_st_name);
    if (display_prop==0) strcat(state_name[i],x_unit);
    else  strcat(state_name[i],x_prop_unit);

}
/* comment out these 4 lines if excluding total*/
     if (i== N){
     strcpy(state_name[i],"BATotal");
     strcat(state_name[i],x_unit);
}

/* repeat for numbers of trees. This applies to Jabowa */
if ((i < (2*N +1)) & (i > N)) {
    if (display_prop==0) strcpy(state_name[i], "Num");
    else strcpy(state_name[i], "RelNum");
    if(i<=9){
	c_st_name[0]= i - N - 1 + 49;
	c_st_name[1]='\0';
    }
    else strcpy(c_st_name,"0");
      /* ascii code of the state number assuming base 48*/

    strcat(state_name[i],c_st_name);
    if (display_prop==0) strcat(state_name[i],num_unit);
    else  strcat(state_name[i],x_prop_unit);

}
/* comment out these 4 lines if excluding total*/
     if (i== 2*N + 1){
     strcpy(state_name[i],"NumTotal");
     strcat(state_name[i],num_unit);
}

if(DEBUG_SW == 1) fprintf(dfp," State=%d state_name=%s\n",i+1,state_name[i]);
}

/*calculation of max number of runs */
num_runs=in_as_round1(par_sens_max,par_sens_min,par_sens_step);
/* calculation of number of output data*/
time_print_f= ((float) time_print)*dt;
num_out_data= in_as_round1(time_final,time_zero,time_print_f);

/* print header in output data file*/
  fprintf(ofp," Number_of_runs: %d", num_runs);
  fprintf(ofp,"\n Number_of Variables (States+1): %d", num_states_print);
  fprintf(ofp,"\n Number of data: %d", num_out_data);
  fprintf(ofp,"\n Title: %s", title);
  fprintf(ofp,"\n Sub_title: %s", sub_title);
  fprintf(ofp,"\n Y_axis_title (depends on absolute or propor: %s", y_axis_title);
  fprintf(ofp,"\n %s", x_axis_title);
for(i=0;i<=num_states_print-1;i++){
fprintf(ofp," %s", state_name[i]);
}
fprintf(ofp,"\n");

  /* print header in index data file*/
  fprintf(efp," Number_of_runs: %d", num_runs);
  fprintf(efp,"\n Number_of Variables (States+1): %d", num_states_print);
  fprintf(efp,"\n Number of data: %d", num_out_data);
  fprintf(efp,"\n Title: %s", title);
  fprintf(efp,"\n Sub_title: %s", sub_title);
  fprintf(efp,"\n Y_axis_title (numbers): %s", y_axis_title);
  fprintf(efp,"\n Index");
for(i=0;i<=num_states_print-1;i++){
fprintf(efp," %s", state_name[i]);
}
fprintf(efp,"\n");

/* simulation*/

/* initialize parameter for multiple runs*/
param.par_array.par[par_sens_flag]=par_sens_min;

/*loop for multiple runs */
for (runs=0; runs<=num_runs-1; runs++){

/* header for file_out */
fprintf(ofp," Simulation run number=%d", runs+1);
fprintf(ofp,":  %s=%6.2e\n",par_sens_name, param.par_array.par[par_sens_flag]);

/* header for file_dex */
fprintf(efp," Simulation run number=%d", runs+1);
fprintf(efp,":  %s=%6.2e\n",par_sens_name, param.par_array.par[par_sens_flag]);

/* echo to scrn for debugging mode */
if(DEBUG_SW == 1){
fprintf(dfp," Simulating run number...%d\n", runs+1);	   /* msg to scrn*/
fprintf(dfp,"param.array[flag]=%6.2e\n",param.par_array.par[par_sens_flag]);
fprintf(dfp," %s", x_axis_title);
for(i=0;i<=N-1;i++){
fprintf(dfp," %s", state_name[i]);
}
fprintf(dfp,"\n");
}      /* end of if debug_sw */

/* calculate height-dbh allometric coeffs and max d and h product*/
    for(i=0;i<=N-1;i++){
	coeff1[i]= 2.00 * (param.par_name.height_max[i] - 137.00)/
			  param.par_name.dbh_max[i];
	coeff2[i]= coeff1[i]/(2.00 * param.par_name.dbh_max[i]);

       /* calculate product of dbh and height max */
       dh_max[i] = param.par_name.dbh_max[i] * param.par_name.height_max[i];

    }	      /*end of i loop*/

    /*initialize external variables*/
    dis_cnt=1.0;

/* initialize simulation clock */
	time=time_zero;

/* initialize dependent variable and derived variables */
    for(i=0;i<=N-1;i++){
	dbh[i]=param.par_name.dbh0[i];
	age[i]= param.par_name.age0[i];
	/* height from dbh using allometric relation*/
	height[i]= height_dbh(coeff1[i], coeff2[i], dbh[i]);
    }	      /*end of i loop*/

/* assigning state */
    for(i=0;i<=N-1;i++){
	/* the second set of N states are number of trees */
	x[i+N]=param.par_name.num0[i];
	/* the first N states are basal areas */
	/* basal area assuming all trees of each species have same dbh*/
	x[i]= basal(dbh[i], x[i+N]);
    }	      /*end of i loop*/

/*calculation of proportions of basal area */
    sum_x = 0.00;
    for (i=0;i<=N-1;i++) sum_x = sum_x + x[i];
    bas_total= sum_x;

    for(i=0;i<=N-1;i++){
	if (bas_total == 0.00) x_prop[i] = 0.00;
	else x_prop[i]= x[i]/bas_total;
    }	      /*end of i loop*/
/* now repeat calc of proportions (second set use numbers) */
    sum_x = 0.00;
    for (i=N;i<=2*N-1;i++) sum_x = sum_x + x[i];
    num_total= sum_x;

    for(i=N;i<=2*N-1;i++){
	if (num_total == 0.00) x_prop[i] = 0.00;
	else x_prop[i]= x[i]/num_total;
    }	      /*end of i loop*/

/*initialize index or endpoints*/
    for(i=0;i<=2*N-1;i++){
      x_max[i]=x[i];
      x_min[i]=x[i];
      x_avg[i]=x[i];
      x_rms[i]= (x[i] - x_avg[i]) * (x[i] - x_avg[i]);
      x_end[i]=x[i];

      x_prop_max[i]=x_prop[i];
      x_prop_min[i]=x_prop[i];
      x_prop_avg[i]=x_prop[i];
      x_prop_rms[i]= (x_prop[i] - x_prop_avg[i]) *
		    (x_prop[i] - x_prop_avg[i]);
      x_prop_end[i]=x_prop[i];

   }   /* end of i loop */
      bas_total_max = bas_total;
      bas_total_min = bas_total;
      bas_total_avg = bas_total;
      bas_total_rms = (bas_total - bas_total_avg)
		    * (bas_total - bas_total_avg);
      bas_total_end = bas_total;

      num_total_max = num_total;
      num_total_min = num_total;
      num_total_avg = num_total;
      num_total_rms = (num_total - num_total_avg)
		    * (num_total - num_total_avg);
      num_total_end = num_total;

      /*first time to print result*/
      it_cnt=0;
      /*write to output file and echo to scrn if debug*/

	fprintf(ofp," %6.2f",it_cnt*dt*time_print + time_zero);

	    for(i=0;i<=N-1;i++){
if (display_prop==0){
	     fprintf(ofp," %10.2e",x[i]);
}
else{
	     fprintf(ofp,"%10.4f",x_prop[i]);	  /*write to output file*/
}
    }  /*-- end of state loop for print */

/*write to output file the total and proportions*/

     /*	 comment out if total not written  */
  fprintf(ofp," %10.2e",bas_total);

for(i=N;i<=2*N-1;i++){
if (display_prop==0){
	     fprintf(ofp," %10.2e",x[i]);
}
else{
	     fprintf(ofp,"%10.4f",x_prop[i]);	  /*write to output file*/
}
    }  /*-- end of state loop for print */

/*write to output file the total and proportions*/

     /*	 comment out if total not written  */
  fprintf(ofp," %10.2e",num_total);

	     fprintf(ofp,"\n");


    /* ---- loop for one run */
    for (it_cnt=1; it_cnt <= (num_out_data-1); it_cnt++){


      /*calculation loop until time to print */
      for(prn_cnt=1;prn_cnt <= time_print; prn_cnt++){

   /*--------------------------------------------*/
      /* numerical integration; euler*/
   /*-----------------------------*/
   /*configuring the rate equation*/
    bas_total= 0.00;
    for(i=0;i<=N-1;i++){
      height[i]= height_dbh(coeff1[i], coeff2[i], dbh[i]);
      x[i]= basal(dbh[i], x[i + N]);
      bas_total+= x[i];
      leaf_area[i] = param.par_name.la_dbh[i] * dbh[i] * dbh[i];
    } /*end of short i loop for first set of variables*/

    for(i=0;i<=N-1;i++){
     above_leaf[i]= 0.00;
     for(j=0;j<=N-1;j++){
      if(height[j] > height[i] ) above_leaf[i]+= leaf_area[j];
     } /*close j loop*/
     /* calculate light attenuation */
     light[i] = param.par_name.zlight;
     light[i] = light[i] * exp( -param.par_name.kext * above_leaf[i] );
     /* light multiplier */
     mult_light[i] = res_light(param.par_name.k1_light[i],
	param.par_name.k2_light[i], param.par_name.k3_light[i], light[i]);
     /* temp multiplier */

     mult_temp[i] = res_temp(param.par_name.temp,
	 param.par_name.dd_min[i], param.par_name.dd_max[i]);

     /* space crowding (lumps moisture + nutrients) */
     if (param.par_name.basmax == 0.00) mult_space=1.00; /* avoid div by zero*/
     else mult_space = 1.00 - bas_total/param.par_name.basmax;
     growth[i] = param.par_name.gmax[i] *
	mult_light[i] * mult_temp[i] * mult_space;

     if (dbh[i] == 0.00) growth[i]=0.00; /* avoid div by zero*/
     else growth[i]= growth[i]*leaf_area[i]/dbh[i];

     if (dh_max[i] == 0.00) fract[i]=0.00; /* avoid div by zero*/
     else fract[i] = dbh[i] * height[i] / dh_max[i];

     growth[i]= growth[i]* (1.00 - fract[i]);
     deno[i] = dh_deno(coeff1[i], coeff2[i], dbh[i]);
     if (deno[i] == 0.00) growth[i]=0.00; /* avoid div by zero*/
     else  growth[i]= growth[i] / deno[i];

     net_rate[i]= growth[i];

    }	     /* end of i loop for diameter integration*/


     for(i=N;i<=2*N-1;i++){	/* loop for second set, the number of trees */

      /*demographics */

      /* updating age*/
      age[i-N]+= dt;

      /* adjusting mortality according to growth rate*/
      mortage[i-N] = mort_rate ( growth[i-N], param.par_name.k_slow,
		param.par_name.thresh_slow, param.par_name.k_mort[i-N],
		age[i-N]);

      /* as in swartzmann book */
      mort[i]= - mortage[i-N] * x[i];

      /*mort[i]= mortage[i-N] * x[i];	as implemented in time zero */

      }	     /* end of i loop for second set, the number of trees*/

   /* ----------------------------*/
    /*impulsive disturbance on the number of trees */
	 /* for diameter is zero */
	 for(i=0;i<=N-1;i++) {non_rate[i]=0.00;}

     if (time < dis_cnt * param.par_name.period){
	   for(i=N;i<=2*N-1;i++) {non_rate[i]=0.00;}
     }
		/* disturbance does not occur*/

     else {					 /* disturbance occurs*/
	   sum_non_rate=0.00;	/* initialize to accumulate cleared area */
	   for(i=N;i<=2*N-1;i++){
	   if (prop_sw == 1) non_rate[i]= param.par_name.intensity[i] * x[i];
			  /*proportional*/
	   else non_rate[i]= param.par_name.intensity[i];
			  /* disturbance fixed */
	   sum_non_rate+= non_rate[i];
	   }
	++dis_cnt;				 /* counter disturbance*/

	 /* harvested numbers go to first pioneer species with small dbh*/
	   x[N]+= -sum_non_rate;
	 /*  dbh[0] = 0.5;	*/	 /* reset pioneer's diameter */

     }

    /*updating state */
    for(i=0;i<=N-1;i++){
    dbh[i]+= dt*net_rate[i] + non_rate[i];
    if (x[i] < 0.00) x[i]=0.00; /* force to zero if negative*/
    if (x[i] > BOUND) x[i] = BOUND; /* avoid overflow*/
   }	 /*end of state loop*/

    for(i=N;i<=2*N-1;i++){

/*    x[i]= mort[i] + non_rate[i];  as in time zero */

    /* adapted from Swartzman but using notion of ODE*/
    x[i]+= dt*mort[i] + non_rate[i];
    if (x[i] < 0.00) x[i]=0.00; /* force to zero if negative*/
    if (x[i] > BOUND) x[i] = BOUND; /* avoid overflow*/
   }	 /*end of state loop*/



   /* end of numerical integration */
   /* ----------------------------*/

	 /* for jabowa update basal area */
    for(i=0;i<=N-1;i++){
    x[i]= basal(dbh[i], x[i+N]);
    }

/*calculation of proportions (for the first set use basal area)*/
    sum_x = 0.00;
    for (i=0;i<=N-1;i++) sum_x = sum_x + x[i];
    bas_total= sum_x;

    for(i=0;i<=N-1;i++){
	if (bas_total == 0.00) x_prop[i] = 0.00;
	else x_prop[i]= x[i]/bas_total;
    }	      /*end of i loop*/
/* now repeat calc of proportions (second set use numbers) */
    sum_x = 0.00;
    for (i=N;i<=2*N-1;i++) sum_x = sum_x + x[i];
    num_total= sum_x;

    for(i=N;i<=2*N-1;i++){
	if (num_total == 0.00) x_prop[i] = 0.00;
	else x_prop[i]= x[i]/num_total;
    }	      /*end of i loop*/


      /*calculate temporary index or endpoints: max, min, mean, rms, end*/
    for(i=0;i<=2*N-1;i++){
	x_max[i]=end_max(x[i],x_max[i]);
	x_min[i]=end_min(x[i],x_min[i]);
	x_avg[i]=end_avg(x[i],x_avg[i],dt,time-time_zero);
	x_rms[i]=end_rms(x[i],x_rms[i],x_avg[i],dt,time-time_zero);
	x_end[i]=x[i];

	x_prop_max[i]=end_max(x_prop[i],x_prop_max[i]);
	x_prop_min[i]=end_min(x_prop[i],x_prop_min[i]);
	x_prop_avg[i]=end_avg(x_prop[i],x_prop_avg[i],dt,time-time_zero);
	x_prop_rms[i]=end_rms(x_prop[i],x_prop_rms[i],x_prop_avg[i],dt,time-time_zero);
	x_prop_end[i]=x_prop[i];

    } /* end of state loop */
	bas_total_max =end_max(bas_total,bas_total_max);
	bas_total_min =end_min(bas_total,bas_total_min);
	bas_total_avg =end_avg(bas_total,bas_total_avg,dt,time-time_zero);
	bas_total_rms =end_rms(bas_total,bas_total_rms,bas_total_avg,dt,time-time_zero);
	bas_total_end =bas_total;

	num_total_max =end_max(num_total,num_total_max);
	num_total_min =end_min(num_total,num_total_min);
	num_total_avg =end_avg(num_total,num_total_avg,dt,time-time_zero);
	num_total_rms =end_rms(num_total,num_total_rms,num_total_avg,dt,time-time_zero);
	num_total_end =num_total;

    time+=dt;		 /*inc time*/

    } /* end of calculation loop before printing results */

      /*time to print result*/
      /*write to output file and echo to scrn if debug*/

	fprintf(ofp," %6.2f",it_cnt*dt*time_print + time_zero);

	    for(i=0;i<=N-1;i++){
if (display_prop==0){
	     fprintf(ofp," %10.2e",x[i]);
}
else{
	     fprintf(ofp,"%10.4f",x_prop[i]);	  /*write to output file*/
}
    }  /*-- end of state loop for print */

/*write to output file the total and proportions*/

     /*	 comment out if total not written  */
  fprintf(ofp," %10.2e",bas_total);

	    for(i=N;i<=2*N-1;i++){
if (display_prop==0){
	     fprintf(ofp," %10.2e",x[i]);
}
else{
	     fprintf(ofp,"%10.4f",x_prop[i]);	  /*write to output file*/
}
    }  /*-- end of state loop for print */

/*write to output file the total and proportions*/

     /*	 comment out if total not written  */
  fprintf(ofp," %10.2e",num_total);

	     fprintf(ofp,"\n");


    }		/* ---------- end of simulation run loop */

    /*finalize index*/
    for(i=0;i<=2*N-1;i++){
      x_rms[i]= sqrt(x_rms[i]);
      x_prop_rms[i] = sqrt(x_prop_rms[i]);
    }
    bas_total_rms = sqrt(bas_total_rms);
    num_total_rms = sqrt(num_total_rms);


    /* write to index or endpoint file */

	fprintf(efp," max");
	    for(i=0;i<=N-1;i++){
if (display_prop==0){
	     fprintf(efp," %10.2e",x_max[i]);	  /*write to endpoint file*/
}
else{
	     fprintf(efp," %10.4f",x_prop_max[i]);	  /*write to endpoint file*/
}
    }  /* end of state loop for printing*/

/* commented out to exclude total*/
	     fprintf(efp," %10.2e",bas_total_max);

	    for(i=N;i<=2*N-1;i++){
if (display_prop==0){
	     fprintf(efp," %10.2e",x_max[i]);	  /*write to endpoint file*/
}
else{
	     fprintf(efp," %10.4f",x_prop_max[i]);	  /*write to endpoint file*/
}
    }  /* end of state loop for printing*/

/* commented out to exclude total*/
	     fprintf(efp," %10.2e",num_total_max);

	     fprintf(efp,"\n");

	fprintf(efp," min");
	    for(i=0;i<=N-1;i++){
if (display_prop ==0){
	     fprintf(efp," %10.2e",x_min[i]);	  /*write to endpoint file*/
}
else{
	     fprintf(efp," %10.4f",x_prop_min[i]);	  /*write to endpoint file*/
}
    }  /* end of state loop for printing*/
/* commented out to exclude total*/
	     fprintf(efp," %10.2e",bas_total_min);
	    for(i=N;i<=2*N-1;i++){
if (display_prop ==0){
	     fprintf(efp," %10.2e",x_min[i]);	  /*write to endpoint file*/
}
else{
	     fprintf(efp," %10.4f",x_prop_min[i]);	  /*write to endpoint file*/
}
    }  /* end of state loop for printing*/
/* commented out to exclude total*/
	     fprintf(efp," %10.2e",num_total_min);

	     fprintf(efp,"\n");

	fprintf(efp," avg");
	    for(i=0;i<=N-1;i++){
if(display_prop ==0){
	     fprintf(efp," %10.2e",x_avg[i]);	  /*write to endpoint file*/
}
else{
	     fprintf(efp," %10.4f",x_prop_avg[i]);	  /*write to endpoint file*/
}
    }  /* end of state loop for printing*/
/* commented out to exclude total*/
	     fprintf(efp," %10.2e",bas_total_avg);

	    for(i=N;i<=2*N-1;i++){
if(display_prop ==0){
	     fprintf(efp," %10.2e",x_avg[i]);	  /*write to endpoint file*/
}
else{
	     fprintf(efp," %10.4f",x_prop_avg[i]);	  /*write to endpoint file*/
}
    }  /* end of state loop for printing*/
/* commented out to exclude total*/
	     fprintf(efp," %10.2e",num_total_avg);

	     fprintf(efp,"\n");

	fprintf(efp," rms");
	    for(i=0;i<=N-1;i++){
if (display_prop==0){
	     fprintf(efp," %10.2e",x_rms[i]);	  /*write to endpoint file*/
}
else{
	     fprintf(efp," %10.4f",x_prop_rms[i]);	  /*write to endpoint file*/
}
    }  /* end of state loop for printing*/
/* commented out to exclude total*/
	     fprintf(efp," %10.2e",bas_total_rms);
	    for(i=N;i<=2*N-1;i++){
if (display_prop==0){
	     fprintf(efp," %10.2e",x_rms[i]);	  /*write to endpoint file*/
}
else{
	     fprintf(efp," %10.4f",x_prop_rms[i]);	  /*write to endpoint file*/
}
    }  /* end of state loop for printing*/
/* commented out to exclude total*/
	     fprintf(efp," %10.2e",num_total_rms);

	     fprintf(efp,"\n");

	fprintf(efp," end");
	    for(i=0;i<=N-1;i++){
if(display_prop==0){
	     fprintf(efp," %10.2e",x_end[i]);	  /*write to endpoint file*/
}
else{
	     fprintf(efp," %10.4f",x_prop_end[i]);	  /*write to endpoint file*/
}
    }  /* end of state loop for printing*/
/* commented out to exclude total*/
	     fprintf(efp," %10.2e",bas_total_end);
	    for(i=N;i<=2*N-1;i++){
if(display_prop==0){
	     fprintf(efp," %10.2e",x_end[i]);	  /*write to endpoint file*/
}
else{
	     fprintf(efp," %10.4f",x_prop_end[i]);	  /*write to endpoint file*/
}
    }  /* end of state loop for printing*/
/* commented out to exclude total*/
	     fprintf(efp," %10.2e",num_total_end);

	     fprintf(efp,"\n");

/*change value of parameter for next run*/
   param.par_array.par[par_sens_flag]+=par_sens_step;
}	       /* ------- end of run loop iteration--------*/

/*close files and exit*/
fclose(ifp);
fclose(ofp);
fclose(efp);
if(DEBUG_SW == 1) fclose(dfp);

}   /* -------------- end of main program -------------*/

/*function to evaluate height-diameter allometry*/
float height_dbh(float coeff1, float coeff2, float dbh)
{
float height;
height= 137 + coeff1*dbh - coeff2*dbh*dbh;
return height;
}

/*function to evaluate basal area*/
float basal(float dbh, float num)
{
float basal_area;
basal_area= (3.1416/4.00) * dbh * dbh * num;
return basal_area;
}

/* function to evaluate shade growth restriction */
float res_light(float k1, float k2, float k3, float light)
{
float mult;
mult = k1 * ( 1.00 - exp( -k2 * (light - k3) ) );
return mult;
}

/* function to evaluate thermal growth restriction */
float res_temp(float temp, float min, float max)
{
float mult;
if ( (max-min) == 0.00) mult = 0.00;
else mult = 4.00 * (temp - min) * (max - temp) / ( (max-min)*(max-min) );
return mult;
}

/* function to eavluate denominator in growth function*/
float dh_deno(float coeff1, float coeff2, float dbh)
{
float deno;
deno = 274.00 + 3.00 * coeff1 * dbh - 4.00 * coeff2 * dbh * dbh;
return deno;
}

/* function to evaluate mortality rates*/
float mort_rate(float growth_rate, float slow, float thresh, float mu, float age)
{
float rate, enhance;
      rate = 1.00 - pow( (1.00 - mu), age );
      enhance = 1.00 - pow( (1.00 - slow), age );
      if (growth_rate < thresh ) rate+= enhance;

/* the following is as implemented in the time0 version */
/*	if (growth_rate > thresh) enhance = 0.00;
	else enhance = slow;
	rate = pow ((1.00 - mu - enhance), age);
*/
return rate;
}

/*function to calculate integer */
int in_as_round1(float max_val, float min_val, float step_val)
{
float val;
int  n_val;
if (step_val == 0.00) val = 0.00;
  else val= ((max_val - min_val) / step_val);
if (val - floor(val) >= 0.5)
    n_val= ceil(val)+1;
  else
    n_val= floor(val)+1;
return n_val;
}

/*function to evaluate max index or endpoint*/
float end_max (float x_val, float val_max)
{
float maximum;
if (x_val >= val_max) maximum=x_val;
    else maximum=val_max;
return maximum;
}

/*function to evaluate min index or endpoint */
float end_min (float x_val, float val_min)
{
float minimum;
if (x_val <= val_min) minimum=x_val;
    else minimum=val_min;
return minimum;
}

/*function to evaluate mean index or endpoint */
float end_avg (float x_val, float val_mean, float delta_t,
	       float previous_time)
{
float mean;
       mean = (val_mean * previous_time/delta_t + x_val) *
       delta_t/(delta_t + previous_time);
return mean;
}

/*function to evaluate rms index or endpoint */
float end_rms (float x_val, float val_rms, float val_mean, float delta_t,
	       float previous_time)
{
float rms;
rms = (x_val - val_mean) * (x_val - val_mean);
       rms = (val_rms * previous_time/delta_t + rms) *
       delta_t/(delta_t + previous_time);
return rms;
}
