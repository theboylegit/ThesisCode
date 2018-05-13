#include<python3.6m/Python.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<unistd.h>
#include <gsl/gsl_spline.h>
#include "imf.h"
double plot_fun(int j, char *Name_of_Run, char *M_star, double low_a, double increment_a, char *Pathto_stellar_mods, int count_R_star){
int i, engulf;
char save_plot[500];
FILE *PLOT_FILE = fopen("./plot_scripts/plot_orbital_evolution.py", "w");
FILE *ENGULFMENTS = fopen("engulfed.txt", "r");

fprintf(PLOT_FILE, "import sys\nimport matplotlib\nmatplotlib.use('Agg')\nimport matplotlib.pyplot as plt\nimport numpy as np\nimport csv\n\nage = []\nradius = []");
fprintf(PLOT_FILE,"\nx=[]\ny=[]\nx1=[]\ny1=[]\nx2=[]\ny2=[]");
for(i=0; i<j; i++){
	fprintf(PLOT_FILE, "\n\nsma%d = []\nRs%d = []\ntime%d = []\n\n", i,i,i);
}
for(i=0; i<j; i++){
	fprintf(PLOT_FILE, "\n\nwith open('./results_and_plots/%s_results/%s%sstar_%lfAU.txt', 'r') as csvfile:\n	plots = csv.reader(csvfile, delimiter=' ')\n	for row in plots:\n		if len(row) > 1:\n			sma%d.append(float(row[0]))\n			Rs%d.append(float(row[1]))\n			time%d.append(float(row[2]))", Name_of_Run, Name_of_Run,M_star, low_a,i,i,i);
	low_a = low_a + increment_a;
}

fprintf(PLOT_FILE, "\n\n\nwith open('./results_and_plots/%s_results/%sRvA.txt', 'r') as csvfile:\n	plots = csv.reader(csvfile, delimiter=' ')\n	for row in plots:\n		if len(row) > 1:\n			age.append(float(row[0]))\n			radius.append(float(row[1]))\n", Name_of_Run, Name_of_Run);
for(i=0; i<j; i++){
	fscanf(ENGULFMENTS,"%d\n", &engulf);
	if (engulf>0 && engulf < 2){
		fprintf(PLOT_FILE, "\nplt.plot(time%d,sma%d, color = 'purple', linewidth=1.0)\n", i, i);
	}
	if(engulf == 2){
		fprintf(PLOT_FILE, "\nplt.plot(time%d,sma%d, color = 'deepskyblue', linewidth=1.0)\n", i, i);
	}
	if(engulf == 0){
		fprintf(PLOT_FILE, "\nplt.plot(time%d,sma%d, color = 'midnightblue', linewidth=1.0)\n", i, i);
	}
}

fprintf(PLOT_FILE,"\nplt.plot(x,y, label = 'RGB engulfment zone', color = 'purple', linewidth=1.0)");
fprintf(PLOT_FILE,"\nplt.plot(x1,y1, label = 'AGB engulfment zone', color = 'deepskyblue', linewidth=1.0)");
fprintf(PLOT_FILE,"\nplt.plot(x2,y2, label = 'Safe semi-major axis', color = 'midnightblue', linewidth=1.0)");
fprintf(PLOT_FILE, "\n\n\nplt.plot(age,radius, label='Stellar radius', color='black', linewidth=2.0)\nplt.grid(False)\nplt.legend(loc='upper left', fontsize=12)\nplt.xlim(age[3],age[%d-10])\nplt.ylim(top=(max(radius)*2))\n\nplt.xlabel('Stellar age(yrs)', fontsize=12)\nplt.ylabel('Semi-major axis, R$_{*}$ ($R_{\\odot}$)', fontsize=12)", count_R_star);

sprintf(save_plot, "\nplt.savefig('./results_and_plots/%s_Plots/%s_Evolution.pdf')", Name_of_Run, Name_of_Run);
fprintf(PLOT_FILE, "%s", save_plot);

fclose(PLOT_FILE);
}

double interpolate(int engulf_point, char *data, char *age_data, char *stellar_mod_path, char *name){

int lines, j, low, high;
char save_new_data[100], save_new_age[100];
char ch;
size_t i;

FILE *LINES = fopen(age_data, "r");
while(!feof(LINES))
{
ch = fgetc(LINES);
	if(ch == '\n')
	{
	lines++;
}															
} 
fclose(LINES);
const size_t N = 20;
double x[lines];
double y[lines];
double data_value[lines];
double age[lines];

sprintf(save_new_data, "%s/interp_%s.txt", stellar_mod_path, name);
sprintf(save_new_age, "%s/interp_age.txt", stellar_mod_path);

FILE *DATA = fopen(data, "r");
FILE *AGEDATA = fopen(age_data, "r");
FILE *INTERP_DATA = fopen(save_new_data, "w");
FILE *INTERPAGE = fopen(save_new_age, "w");


for(i=0;i<(lines);i++){

fscanf(DATA, "%lf", &data_value[i]);
fscanf(AGEDATA, "%lf", &age[i]);
} 


low = engulf_point - 5;
high = engulf_point +15;

for (i=low; i<high; i++)
{
j = i-low;

x[j] = age[i];
y[j] = data_value[i];
}

for (i=0;i<low;i++)
{

fprintf(INTERP_DATA, "%lf\n",data_value[i]);
fprintf(INTERPAGE, "%lf\n", age[i]);
}




  gsl_interp_accel *acc = gsl_interp_accel_alloc();

  gsl_spline *spline_steffen = gsl_spline_alloc(gsl_interp_steffen, N);


  gsl_spline_init(spline_steffen, x, y, N);

  for (i = 0; i < N; ++i)


  for (i = 0; i <= 10000; ++i)
    {
      double xi = (1 - i / 10000.0) * x[0] + (i / 10000.0) * x[N-1];

      double yi_steffen = gsl_spline_eval(spline_steffen, xi, acc);
	 // printf("%ef\n", yi_steffen);
		fprintf(INTERP_DATA, "%lf\n",yi_steffen);
		fprintf(INTERPAGE, "%lf\n", xi);

    }

 for (i=high; i<lines; i++)
{
fprintf(INTERP_DATA, "%lf\n",data_value[i]);
fprintf(INTERPAGE, "%lf\n", age[i]);
}    

gsl_spline_free(spline_steffen);
gsl_interp_accel_free(acc);
fclose(INTERPAGE);
fclose(INTERP_DATA);
fclose(DATA);
fclose(AGEDATA);
  return 0;

}

double RK45_dopri(double dt, double a, double Menv, double Ms, double Mdot, double L, double Rs, double Mpl, double omega_star){
double period, omega_pl, f, Tc;
double a1, a2, a3, a4, a5, a6;
double k1, k2, k3, k4, k5, k6;

Tc = pow((Menv * (Rs * Rs) /(3.0 * L)),1.0/3.0);

//function
//k1
period = sqrt((a/215.0)*(a/215.0)*(a/215.0));
omega_pl = (2.0*3.14) / period;
if(Tc > (period/2.0)){
	f = (period / (2.0*Tc))*(period / (2.0*Tc));
}
else{
	f=1.0;
}
k1 = dt* ( (((-f*a)/Tc) * (Menv/Ms) * (1.0 + (Mpl/Ms)) * (Mpl/Ms) * (pow((Rs/a),8.0)) * (1.0 - (omega_star/omega_pl))) - (a*Mdot / (Ms+Mpl)) );

//k2
a2 = a + ((1.0/5.0)*k1);
period = sqrt((a2/215.0)*(a2/215.0)*(a2/215.0));
omega_pl = (2.0*3.14) / period;
if(Tc > (period/2.0)){
	f = (period / (2.0*Tc))*(period / (2.0*Tc));
}
else{
	f=1.0;
}
k2 = dt* ( (((-f*a2)/Tc) * (Menv/Ms) * (1.0 + (Mpl/Ms)) * (Mpl/Ms) * (pow((Rs/a2),8.0)) * (1.0 - (omega_star/omega_pl))) - (a2*Mdot / (Ms+Mpl)) );
//k3
a3 = a + (((3.0/40.0)*k1) + ((9.0/40.0)*k2));
period = sqrt((a3/215.0)*(a3/215.0)*(a3/215.0));
omega_pl = (2.0*3.14) / period;
if(Tc > (period/2.0)){
	f = (period / (2.0*Tc))*(period / (2.0*Tc));
}
else{
	f=1.0;
}
k3 = dt* ( (((-f*a3)/Tc) * (Menv/Ms) * (1.0 + (Mpl/Ms)) * (Mpl/Ms) * (pow((Rs/a3),8.0)) * (1.0 - (omega_star/omega_pl))) - (a3*Mdot / (Ms+Mpl)) );

//k4
a4 = a + (((44.0/45.0)*k1) - ((56.0/15.0)*k2)+ ((32.0/9.0)*k3));	
period = sqrt((a4/215.0)*(a4/215.0)*(a4/215.0));
omega_pl = (2.0*3.14) / period;
if(Tc > (period/2.0)){
	f = (period / (2.0*Tc))*(period / (2.0*Tc));
}
else{
	f=1.0;
}
k4 = dt* ( (((-f*a4)/Tc) * (Menv/Ms) * (1.0 + (Mpl/Ms)) * (Mpl/Ms) * (pow((Rs/a4),8.0)) * (1.0 - (omega_star/omega_pl))) - (a4*Mdot / (Ms+Mpl)) );

//k5
a5 = a + (((19372.0/6561.0)*k1) - ((25360.0/2187.0)*k2) + ((64448.0/6561.0)*k3) - ((212.0/729.0)*k4));
period = sqrt((a5/215.0)*(a5/215.0)*(a5/215.0));
omega_pl = (2.0*3.14) / period;
if(Tc > (period/2.0)){
	f = (period / (2.0*Tc))*(period / (2.0*Tc));
}
else{
	f=1.0;
}
k5 = dt* ( (((-f*a5)/Tc) * (Menv/Ms) * (1.0 + (Mpl/Ms)) * (Mpl/Ms) * (pow((Rs/a5),8.0)) * (1.0 - (omega_star/omega_pl))) - (a5*Mdot / (Ms+Mpl)) );

//k6
a6 = a + (((9017.0/3168.0)*k1) - ((355.0/33.0)*k2) - ((46732.0/5247.0)*k3) + ((49.0/176.0)*k4) - ((5103.0/18656.0)*k5));
period = sqrt((a6/215.0)*(a6/215.0)*(a6/215.0));
omega_pl = (2.0*3.14) / period;
if(Tc > (period/2.0)){
	f = (period / (2.0*Tc))*(period / (2.0*Tc));
}
else{
	f=1.0;
}
k6 = dt* ( (((-f*a6)/Tc) * (Menv/Ms) * (1.0 + (Mpl/Ms)) * (Mpl/Ms) * (pow((Rs/a6),8.0)) * (1.0 - (omega_star/omega_pl))) - (a6*Mdot / (Ms+Mpl)) );

////calc new value for a///
a = a + (35.0/384.0)*k1 + (500.0/1113.0)*k3 + (125.0/192.0)*k4 - (2187.0/6784.0)*k5 + (11.0/84.0)*k6;

return a;

}

double RK45_fehlberg(double dt, double a, double Menv, double Ms, double Mdot, double L, double Rs, double Mpl, double omega_star){
double period, omega_pl, f, Tc;
double a1, a2, a3, a4, a5, a6;
double k1, k2, k3, k4, k5, k6;

Tc = pow((Menv * (Rs * Rs) /(3.0 * L)),1.0/3.0);

//function
//
//k1//
period = sqrt((a/215.0)*(a/215.0)*(a/215.0));
omega_pl = (2.0*3.14) / period;
		
if(Tc > (period/2.0)){
	f = (period / (2.0*Tc))*(period / (2.0*Tc));
}
else{
	f=1.0;
}
k1 = dt* ( (((-f*a)/Tc) * (Menv/Ms) * (1.0 + (Mpl/Ms)) * (Mpl/Ms) * (pow((Rs/a),8.0)) * (1.0 - (omega_star/omega_pl))) - (a*Mdot / (Ms+Mpl)) );

//k2//
a2 = a + ((1.0/4.0)*k1);
period = sqrt((a2/215.0)*(a2/215.0)*(a2/215.0));
omega_pl = (2.0*3.14) / period;
		
if(Tc > (period/2.0)){
	f = (period / (2.0*Tc))*(period / (2.0*Tc));
}
else{
	f=1.0;
}
k2 = dt* ( (((-f*a2)/Tc) * (Menv/Ms) * (1.0 + (Mpl/Ms)) * (Mpl/Ms) * (pow((Rs/a2),8.0)) * (1.0 - (omega_star/omega_pl))) - (a2*Mdot / (Ms+Mpl)) );

//k3//
a3 = a + (((3.0/32.0)*k1) + ((9.0/32.0)*k2));
period = sqrt((a3/215.0)*(a3/215.0)*(a3/215.0));
omega_pl = (2.0*3.14) / period;
		
if(Tc > (period/2.0)){
	f = (period / (2.0*Tc))*(period / (2.0*Tc));
}
else{
	f=1.0;
}
k3 = dt* ( (((-f*a3)/Tc) * (Menv/Ms) * (1.0 + (Mpl/Ms)) * (Mpl/Ms) * (pow((Rs/a3),8.0)) * (1.0 - (omega_star/omega_pl))) - (a3*Mdot / (Ms+Mpl)) );

//k4//
a4 = a + (((1932.0/2197.0)*k1) - ((7200.0/2197.0)*k2)+ ((7296.0/2197.0)*k3));	
period = sqrt((a4/215.0)*(a4/215.0)*(a4/215.0));
omega_pl = (2.0*3.14) / period;
		
if(Tc > (period/2.0)){
	f = (period / (2.0*Tc))*(period / (2.0*Tc));
}
else{
	f=1.0;
}
k4 = dt* ( (((-f*a4)/Tc) * (Menv/Ms) * (1.0 + (Mpl/Ms)) * (Mpl/Ms) * (pow((Rs/a4),8.0)) * (1.0 - (omega_star/omega_pl))) - (a4*Mdot / (Ms+Mpl)) );

//k5//
a5 = a + (((439.0/216.0)*k1) - (8*k2) + ((3680.0/513.0)*k3) - ((845.0/4104.0)*k4));
period = sqrt((a5/215.0)*(a5/215.0)*(a5/215.0));
omega_pl = (2.0*3.14) / period;
		
if(Tc > (period/2.0)){
	f = (period / (2.0*Tc))*(period / (2.0*Tc));
}
else{
	f=1.0;
}
k5 = dt* ( (((-f*a5)/Tc) * (Menv/Ms) * (1.0 + (Mpl/Ms)) * (Mpl/Ms) * (pow((Rs/a5),8.0)) * (1.0 - (omega_star/omega_pl))) - (a5*Mdot / (Ms+Mpl)) );

//k6//
a6 = a + (((-8.0/27.0)*k1) + (2.0*k2) - ((3544.0/2565.0)*k3) + ((1859.0/4104.0)*k4) - ((11.0/40.0)*k5));
period = sqrt((a6/215.0)*(a6/215.0)*(a6/215.0));
omega_pl = (2.0*3.14) / period;
		
if(Tc > (period/2.0)){
	f = (period / (2.0*Tc))*(period / (2.0*Tc));
}
else{
	f=1.0;
}
k6 = dt* ( (((-f*a6)/Tc) * (Menv/Ms) * (1.0 + (Mpl/Ms)) * (Mpl/Ms) * (pow((Rs/a6),8.0)) * (1.0 - (omega_star/omega_pl))) - (a6*Mdot / (Ms+Mpl)) );

//Calculate new a///
a = a + ((16.0/135.0)*k1) + ((6656.0/12825.0)*k3) + ((28561.0/56430.0)*k4) + ((9.0/50.0)*k5) + ((2.0/55.0)*k6);
return a;
}

double RK4(double dt, double a, double Menv, double Ms, double Mdot, double L, double Rs, double Mpl, double omega_star){
double period, omega_pl, f, Tc;
double a1, a2, a3, a4, a5, a6;
double k1, k2, k3, k4, k5, k6;

Tc = pow((Menv * (Rs * Rs) /(3.0 * L)),1.0/3.0);

//function
//
period = sqrt((a/215.0)*(a/215.0)*(a/215.0));
omega_pl = (2.0*3.14) / period;
		
if(Tc > (period/2.0)){
	f = (period / (2.0*Tc))*(period / (2.0*Tc));
}
else{
	f=1.0;
}
k1 = dt* ( (((-f*a)/Tc) * (Menv/Ms) * (1.0 + (Mpl/Ms)) * (Mpl/Ms) * (pow((Rs/a),8.0)) * (1.0 - (omega_star/omega_pl))) - (a*Mdot / (Ms+Mpl)) );

//k2//
a2 = a + ((1.0/2.0)*k1);
period = sqrt((a2/215.0)*(a2/215.0)*(a2/215.0));
omega_pl = (2.0*3.14) / period;
		
if(Tc > (period/2.0)){
	f = (period / (2.0*Tc))*(period / (2.0*Tc));
}
else{
	f=1.0;
}
k2 = dt* ( (((-f*a2)/Tc) * (Menv/Ms) * (1.0 + (Mpl/Ms)) * (Mpl/Ms) * (pow((Rs/a2),8.0)) * (1.0 - (omega_star/omega_pl))) - (a2*Mdot / (Ms+Mpl)) );

//k3//
a3 = a + (((1.0/2.0)*k2));
period = sqrt((a3/215.0)*(a3/215.0)*(a3/215.0));
omega_pl = (2.0*3.14) / period;
		
if(Tc > (period/2.0)){
	f = (period / (2.0*Tc))*(period / (2.0*Tc));
}
else{
	f=1.0;
}
k3 = dt* ( (((-f*a3)/Tc) * (Menv/Ms) * (1.0 + (Mpl/Ms)) * (Mpl/Ms) * (pow((Rs/a3),8.0)) * (1.0 - (omega_star/omega_pl))) - (a3*Mdot / (Ms+Mpl)) );

//k4//
a4 = a + k3;	
period = sqrt((a4/215.0)*(a4/215.0)*(a4/215.0));
omega_pl = (2.0*3.14) / period;
		
if(Tc > (period/2.0)){
	f = (period / (2.0*Tc))*(period / (2.0*Tc));
}
else{
	f=1.0;
}
k4 = dt* ( (((-f*a4)/Tc) * (Menv/Ms) * (1.0 + (Mpl/Ms)) * (Mpl/Ms) * (pow((Rs/a4),8.0)) * (1.0 - (omega_star/omega_pl))) - (a4*Mdot / (Ms+Mpl)) );

//calc new a
a = a + ((k1 + (2*k2) + (2*k3) + k4))/6.0;
return a;
}

double SMA_EVOL(char *interp, double a0, char *M_star, char *stellar_metallicity, char *Name_of_Run, char *stellar_model, char *Runge_Kutta, char *plotting, char *M_planet, char *Pathto_stellar_mods){

int engulf_phase, file_count, phase_counter, engulf_point, engulf;
char sep_data_dt[1000],sep_data_mass[1000], sep_data_lum[1000], sep_data_radius[1000], sep_data_mdot[1000], sep_data_phase[1000], sep_data_core[1000], sep_data_omega[1000], sep_data_H1[1000], sep_data_HE4[1000];
char result_string[1000], result_summary[1000], python_HR_call[200], clean_data[100];
char Pathto_history[1000], radius_data[1000], dt_data[1000], mass_data[1000], mdot_data[1000], luminosity_data[1000], phases_data[1000], core_data[1000], omega_data[1000], age_path[1000], dt_path[1000], rm_1dt[1000], H1_data[1000], HE4_data[1000];
double Mpl, Ms, Rs, L, Mdot, dt, h, Mc, t, orbit_checker, a, t0, age1, age2, time, omega_star, omega_pl, Rs_AU, a0_AU, Menv, a_final, h1_center, he4_center, RGB_h1, AGB_he4;
double k1, k2, k3, k4, k5, k6, calc_new_a;
float phase, MS_phase, RGB_phase, AGB_phase;
RGB_h1 = pow(10, -3);
AGB_he4 = pow(10, -2);
MS_phase = 0.0;
RGB_phase = 2.0;
AGB_phase = 4.0;
Mpl = strtof(M_planet, NULL);
Mpl = Mpl *9.5e-4;
engulf_point=0;
phase_counter=0;
file_count = 0;	

sprintf(Pathto_history, "%s/history.data", Pathto_stellar_mods);
if(strcmp(interp, "interp_yes")==0){
	sprintf(radius_data,"%s/radius.txt",Pathto_stellar_mods);
	sprintf(mass_data,"%s/mass.txt",Pathto_stellar_mods);
	sprintf(mdot_data,"%s/mdot.txt",Pathto_stellar_mods);											
	sprintf(luminosity_data,"%s/luminosity.txt",Pathto_stellar_mods);											
	sprintf(phases_data,"%s/phases.txt",Pathto_stellar_mods);
	sprintf(age_path, "%s/age.txt", Pathto_stellar_mods);
	sprintf(core_data, "%s/core.txt", Pathto_stellar_mods);
	sprintf(omega_data, "%s/omega.txt", Pathto_stellar_mods);
	sprintf(dt_data,"%s/dt.txt",Pathto_stellar_mods);	
}

if(strcmp(interp, "interp_no")==0){
	sprintf(radius_data,"%s/interp_radius.txt",Pathto_stellar_mods);
	sprintf(mass_data,"%s/interp_mass.txt",Pathto_stellar_mods);
	sprintf(mdot_data,"%s/interp_mdot.txt",Pathto_stellar_mods);											
	sprintf(luminosity_data,"%s/interp_luminosity.txt",Pathto_stellar_mods);											
	sprintf(phases_data,"%s/interp_phases.txt",Pathto_stellar_mods);
	sprintf(age_path, "%s/interp_age.txt", Pathto_stellar_mods);
	sprintf(core_data, "%s/interp_core.txt", Pathto_stellar_mods);
	sprintf(omega_data, "%s/interp_omega.txt", Pathto_stellar_mods);
	sprintf(dt_data,"%s/interp_dt.txt",Pathto_stellar_mods);	
}
									
									
sprintf(result_summary, "./results_and_plots/%s_results/%s_models_summary.txt", Name_of_Run, Name_of_Run);
FILE *RESULT_SUMMARY = fopen(result_summary, "w");
sprintf(result_string, "./results_and_plots/%s_results/%s%sstar_%lfAU.txt", Name_of_Run, Name_of_Run, M_star, a0);
FILE *RESULT = fopen(result_string, "w");	

if(strcmp(interp, "interp_yes")==0){
	sprintf(clean_data, "sed 's/#//' %s/history.data > system_temp.txt", Pathto_stellar_mods);
	system(clean_data);
	system("sed -i -n -e '/star_age/,$p' system_temp.txt");	
	sprintf(sep_data_mass,"awk '{if(NR==1) for(i=1;i<=NF;i++) { if($i~colname) { colnum=i;break} } else print $colnum}' colname=star_mass system_temp.txt > %s/mass.txt", Pathto_stellar_mods);
	system(sep_data_mass);	
	sprintf(sep_data_mdot,"awk '{if(NR==1) for(i=1;i<=NF;i++) { if($i~colname) { colnum=i;break} } else print $colnum}' colname=star_mdot system_temp.txt > %s/mdot.txt", Pathto_stellar_mods);
	system(sep_data_mdot);	
	sprintf(sep_data_lum,"awk '{if(NR==1) for(i=1;i<=NF;i++) { if($i~colname) { colnum=i;break} } else print $colnum}' colname=log_L system_temp.txt > %s/luminosity.txt", Pathto_stellar_mods);
	system(sep_data_lum);	
	sprintf(sep_data_radius,"awk '{if(NR==1) for(i=1;i<=NF;i++) { if($i~colname) { colnum=i;break} } else print $colnum}' colname=log_R system_temp.txt > %s/radius.txt", Pathto_stellar_mods);
	system(sep_data_radius);			
	sprintf(sep_data_core,"awk '{if(NR==1) for(i=1;i<=NF;i++) { if($i~colname) { colnum=i;break} } else print $colnum}' colname=he_core_mass system_temp.txt > %s/core.txt", Pathto_stellar_mods);
	system(sep_data_core);
	sprintf(sep_data_omega,"awk '{if(NR==1) for(i=1;i<=NF;i++) { if($i~colname) { colnum=i;break} } else print $colnum}' colname=surf_avg_omega system_temp.txt > %s/omega.txt", Pathto_stellar_mods);
	system(sep_data_omega);
	sprintf(sep_data_dt,"awk '{if(NR==1) for(i=1;i<=NF;i++) { if($i~colname) { colnum=i;break} } else print $colnum}' colname=star_age system_temp.txt > %s/age.txt", Pathto_stellar_mods);
	system(sep_data_dt);
	if(strcmp(stellar_model, "MIST_NR")==0 || strcmp(stellar_model, "MIST_R")==0){
		sprintf(sep_data_phase,"awk '{if(NR==1) for(i=1;i<=NF;i++) { if($i~colname) { colnum=i;break} } else print $colnum}' colname=phase system_temp.txt > %s/phases.txt", Pathto_stellar_mods);
		system(sep_data_phase);	
	}
	
	if(strcmp(stellar_model, "MESA_R")==0 || strcmp(stellar_model, "CREATE_NEW_MESA")==0 || strcmp(stellar_model, "MESA_NR")==0){	
		sprintf(sep_data_H1,"awk '{if(NR==1) for(i=1;i<=NF;i++) { if($i~colname) { colnum=i;break} } else print $colnum}' colname=center_h1 system_temp.txt > %s/H1.txt", Pathto_stellar_mods);
		system(sep_data_H1);	
		sprintf(sep_data_HE4,"awk '{if(NR==1) for(i=1;i<=NF;i++) { if($i~colname) { colnum=i;break} } else print $colnum}' colname=center_he4 system_temp.txt > %s/HE4.txt", Pathto_stellar_mods);
		system(sep_data_HE4);	
		sprintf(H1_data, "%s/H1.txt", Pathto_stellar_mods);
		sprintf(HE4_data, "%s/HE4.txt", Pathto_stellar_mods);

		FILE *H1_CENTER = fopen(H1_data, "r");
		FILE *HE4_CENTER = fopen(HE4_data, "r");
		FILE *MESA_PHASES = fopen(phases_data, "w");
	
		while(!feof(H1_CENTER)){
			fscanf(H1_CENTER, "%lf", &h1_center);
			fscanf(HE4_CENTER, "%lf", &he4_center);
	
			if(h1_center > RGB_h1){	
				fprintf(MESA_PHASES, "%f\n", MS_phase);	
			}
			if(h1_center <= RGB_h1 && he4_center > AGB_he4){
				fprintf(MESA_PHASES, "%f\n", RGB_phase);
			}	
			if(he4_center <= AGB_he4){
				fprintf(MESA_PHASES, "%f\n", AGB_phase);
			}		
		}
		fclose(H1_CENTER);
		fclose(HE4_CENTER);
		fclose(MESA_PHASES);
	}	
 	FILE *AGE_FILE = fopen(age_path, "r");
	FILE *DT = fopen(dt_data, "w");	
		 dt = 0;
		 age1 = 0;
		 age2 = 0;
	while(!feof(AGE_FILE)){
		age2 = age1;		
		fscanf(AGE_FILE, "%lf", &age1);
		dt = fabs(age2 - age1);
		fprintf(DT,"%lf\n", dt);				
	}		
	fclose(DT);
	fclose(AGE_FILE);
	sprintf(rm_1dt, "sed -i 1,1d %s", dt_data);		
	system(rm_1dt);		
}
if(strcmp(interp, "interp_no")==0){		
	FILE *INTERP_AGE_FILE = fopen(age_path, "r");
	FILE *INTERP_DT = fopen(dt_data, "w");		
	dt = 0;
	age1 = 0;
	age2 = 0;

	while(!feof(INTERP_AGE_FILE)){
		age2 = age1;		
		fscanf(INTERP_AGE_FILE, "%lf\n", &age1);		
		dt = fabs(age2 - age1);
		fprintf(INTERP_DT,"%lf\n", dt);			
	}		
	fclose(INTERP_DT);
	fclose(INTERP_AGE_FILE);
	sprintf(rm_1dt, "sed -i 1,1d %s", dt_data);		
	system(rm_1dt);		
}


FILE *AGE = fopen(age_path, "r");
FILE *input1 = fopen(radius_data, "r");
FILE *input2 = fopen(luminosity_data, "r");
FILE *input3 = fopen(mass_data, "r");
FILE *input4 = fopen(mdot_data, "r");
FILE *input5 = fopen(dt_data, "r");
FILE *input6 = fopen(phases_data, "r");
FILE *input7 = fopen(core_data, "r");
FILE *input8 = fopen(omega_data, "r");
FILE *test = fopen("./age.txt", "r");

phase_counter=0;
a0 = a0*215.0;
orbit_checker = a0;
while(!feof(AGE)){
				
	fscanf(input1, "%lf", &Rs);						
	fscanf(input2, "%lf", &L);
	fscanf(input3, "%lf", &Ms);
	fscanf(input4, "%lf", &Mdot);
	fscanf(input5, "%lf", &dt);
	fscanf(input6, "%f", &phase);
	fscanf(input7, "%f", &Mc);
	fscanf(input8, "%lf", &omega_star);
	fscanf(AGE, "%lf", &time);
	phase = round(phase);							
	Rs = pow(10.0, Rs);
	L = pow(10.0, L);			
	Menv = Ms-Mc;
			
	if(phase < 2.0){
		phase_counter++;											
	}
	//if(a0 < (0.05*215.0)){
		//engulf = 1;							
		//break;							
	//}											
	if(phase >= 0.0 && a0 > Rs){
		if (strcmp(Runge_Kutta,"rk4") == 0){		
			calc_new_a = RK4(dt, a0, Menv, Ms, Mdot, L, Rs, Mpl, omega_star);	
			a0 = calc_new_a;			
		}		
		if (strcmp(Runge_Kutta,"rk45_fehlberg") == 0){		
			calc_new_a = RK45_fehlberg(dt, a0, Menv, Ms, Mdot, L, Rs, Mpl, omega_star);	
			a0 = calc_new_a;			
		}	
		if (strcmp(Runge_Kutta,"rk45_dopri") == 0){		
			calc_new_a = RK45_dopri(dt, a0, Menv, Ms, Mdot, L, Rs, Mpl, omega_star);	
			a0 = calc_new_a;			
		}

		//Keep logging the timesteps until it reaches end of evolution or until it breaks from one of the conditions
		if (phase > 1.0 && a0 > Rs && phase < 6.0){	
			fprintf(RESULT, "%lf %lf %lf\n", a0, Rs, time);
			file_count++;										
		}												
																						
																																																																																																		
			if(strcmp(interp, "interp_yes")==0){										
				if (phase >1.0 && phase < 3.0){																							
					if (a0 < orbit_checker*0.9 || a0 < Rs){														
						engulf_phase = (int)phase;
						printf("\n\n*******This planet was Engulfed on RGB!*********\n\n");																							
						fprintf(RESULT, "%lf %lf %lf\n", a0, Rs, time);												
						printf("\n\n*******Interpolating....********\n\n");
						engulf = 1;	
						engulf_point = file_count+phase_counter-1;															
						interpolate(engulf_point, radius_data, age_path, Pathto_stellar_mods, "radius");															
						interpolate(engulf_point, mass_data, age_path, Pathto_stellar_mods, "mass");															
						interpolate(engulf_point, luminosity_data, age_path, Pathto_stellar_mods, "luminosity");															
						interpolate(engulf_point, mdot_data, age_path, Pathto_stellar_mods, "mdot");															
						interpolate(engulf_point, core_data, age_path, Pathto_stellar_mods, "core");															
						interpolate(engulf_point, omega_data, age_path, Pathto_stellar_mods, "omega");															
						interpolate(engulf_point, phases_data, age_path, Pathto_stellar_mods, "phases");															
						fclose(input1);
						fclose(input2);
						fclose(input3);
						fclose(input4);
						fclose(input5);
						fclose(input6);
						fclose(input7);
						fclose(input8);
						fclose(RESULT);
						fclose(AGE);															
						break;																	
					}
															
				}
												/* 	if (phase >2.0 && phase < 4.0){
														
															
															if(a0 < orbit_checker){
															engulf = 1;
															printf("engulfed on horizontal branch\n");
															
															}
																} */
															
				if (phase > 3.0 && phase < 6.0){
					if (a0 > (orbit_checker+100)){														
						printf("The orbit has increased too far beyond it's initial position.\nIt will not be engulfed\n");
						fprintf(RESULT, "%lf %lf %lf\n", a0, Rs, time);												
						printf("\n\n*******Interpolating....********\n\n");
						engulf = 0;	
						engulf_point = file_count+phase_counter-1;														
						interpolate(engulf_point, radius_data, age_path, Pathto_stellar_mods, "radius");
						interpolate(engulf_point, mass_data, age_path, Pathto_stellar_mods, "mass");
						interpolate(engulf_point, luminosity_data, age_path, Pathto_stellar_mods, "luminosity");
						interpolate(engulf_point, mdot_data, age_path, Pathto_stellar_mods, "mdot");
						interpolate(engulf_point, core_data, age_path, Pathto_stellar_mods, "core");
						interpolate(engulf_point, omega_data, age_path, Pathto_stellar_mods, "omega");
						interpolate(engulf_point, phases_data, age_path, Pathto_stellar_mods, "phases");
						fclose(input1);
						fclose(input2);
						fclose(input3);
						fclose(input4);
						fclose(input5);
						fclose(input6);
						fclose(input7);
						fclose(input8);
						fclose(RESULT);
						fclose(AGE);						
						break;
					}
					if (a0 < orbit_checker || a0 < Rs){
						engulf_phase = (int)phase;
						printf("\n\n*******This planet was Engulfed on AGB!*********\n\n");																							
						fprintf(RESULT, "%lf %lf %lf\n", a0, Rs, time);												
						printf("\n\n*******Interpolating....********\n\n");
						engulf = 1;	
						engulf_point = file_count+phase_counter-1;															
						interpolate(engulf_point, radius_data, age_path, Pathto_stellar_mods, "radius");
						interpolate(engulf_point, mass_data, age_path, Pathto_stellar_mods, "mass");
						interpolate(engulf_point, luminosity_data, age_path, Pathto_stellar_mods, "luminosity");
						interpolate(engulf_point, mdot_data, age_path, Pathto_stellar_mods, "mdot");
						interpolate(engulf_point, core_data, age_path, Pathto_stellar_mods, "core");
						interpolate(engulf_point, omega_data, age_path, Pathto_stellar_mods, "omega");
						interpolate(engulf_point, phases_data, age_path, Pathto_stellar_mods, "phases");															
						fclose(input1);
						fclose(input2);
						fclose(input3);
						fclose(input4);
						fclose(input5);
						fclose(input6);
						fclose(input7);
						fclose(input8);
						fclose(RESULT);
						fclose(AGE);															
						break;	
					}															
				}
														
				if (phase > 5.0){
					printf("The star has reached the end of evolution.\n planet will not be engulfed.\n");
					fprintf(RESULT, "%lf %lf %lf\n", a0, Rs, time);												
					printf("\n\n*******Interpolating....********\n\n");
					engulf = 0;	
					engulf_point = file_count+phase_counter-100;															
					interpolate(engulf_point, radius_data, age_path, Pathto_stellar_mods, "radius");
					interpolate(engulf_point, mass_data, age_path, Pathto_stellar_mods, "mass");
					interpolate(engulf_point, luminosity_data, age_path, Pathto_stellar_mods, "luminosity");
					interpolate(engulf_point, mdot_data, age_path, Pathto_stellar_mods, "mdot");
					interpolate(engulf_point, core_data, age_path, Pathto_stellar_mods, "core");
					interpolate(engulf_point, omega_data, age_path, Pathto_stellar_mods, "omega");
					interpolate(engulf_point, phases_data, age_path, Pathto_stellar_mods, "phases");
					fclose(input1);
					fclose(input2);
					fclose(input3);
					fclose(input4);
					fclose(input5);
					fclose(input6);
					fclose(input7);
					fclose(input8);
					fclose(RESULT);
					fclose(AGE);															
					break;														
				}																									
			}

											
			if(strcmp(interp, "interp_no")==0){
				if (phase > 1.0 && phase < 3.0){
					if ((a0 - Rs) < 50 || a0 < Rs || a0 < 0.3*orbit_checker){
						a_final = Rs;
						printf("\n\n*******This planet was Engulfed on the RGB!*********\n\n");
						fprintf(RESULT, "%lf %lf %lf\n", a_final, Rs, time);																																
						engulf = 1;
						fclose(input1);
						fclose(input2);
						fclose(input3);
						fclose(input4);
						fclose(input5);
						fclose(input6);
						fclose(input7);
						fclose(input8);
						fclose(RESULT);
						fclose(AGE);															
						break;														
					}																														
				}																	
				if (phase >2.0 && phase < 4.0){
					if(a0 < orbit_checker*0.7){
						engulf = 1;																		
						a_final = Rs;
						engulf_phase = (int)phase;
						printf("\n\n*******This planet was Engulfed on HB!*********\n\n");																							
						fprintf(RESULT, "%lf %lf %lf\n", a_final, Rs, time);
						fclose(input1);
						fclose(input2);
						fclose(input3);
						fclose(input4);
						fclose(input5);
						fclose(input6);
						fclose(input7);
						fclose(input8);
						fclose(RESULT);
						fclose(AGE);															
						break;
					}
				}
				if (phase > 3.0 && phase < 6.0){
					if (a0 > (orbit_checker+150)){
						printf("The orbit has increased too far beyond it's initial position.\nIt will not be engulfed\n");
						//fprintf(RESULT, "%lf %lf %lf\n", a0, Rs, time);												
						engulf = 0;	
						fclose(input1);
						fclose(input2);
						fclose(input3);
						fclose(input4);
						fclose(input5);
						fclose(input6);
						fclose(input7);
						fclose(input8);
						fclose(RESULT);
						fclose(AGE);
						break;
					}
					if (a0 < Rs || (a0 - Rs) < 50 || a0 < orbit_checker*0.8){
						a_final = Rs;
						engulf_phase = (int)phase;
						printf("\n\n*******This planet was Engulfed on AGB!*********\n\n");																							
						fprintf(RESULT, "%lf %lf %lf\n", a_final, Rs, time);																								
						engulf = 2;	
						engulf_point = file_count+phase_counter-1;
						fclose(input1);
						fclose(input2);
						fclose(input3);
						fclose(input4);
						fclose(input5);
						fclose(input6);
						fclose(input7);
						fclose(input8);
						fclose(RESULT);
						fclose(AGE);													
						break;
					}
				}
				if (phase > 5.0){
					printf("The star has reached the end of evolution.\n planet will not be engulfed.\n");
					engulf = 0;
					fclose(input1);
					fclose(input2);
					fclose(input3);
					fclose(input4);
					fclose(input5);
					fclose(input6);
					fclose(input7);
					fclose(input8);
					fclose(RESULT);
					fclose(AGE);
					break;
				}
			}
															
		}
}
fclose(input1);
fclose(input2);
fclose(input3);
fclose(input4);
fclose(input5);
fclose(input6);
fclose(input7);
fclose(input8);
fclose(RESULT);
fclose(AGE);

engulf_point = file_count+phase_counter-1;

if(strcmp(plotting, "yes")==0 && strcmp(interp, "interp_yes")==0){

sprintf(python_HR_call, "python3 ./plot_scripts/plot_HR.py %s %s %s %d %d", Pathto_history, Name_of_Run, Name_of_Run, engulf_phase, engulf_point);

system(python_HR_call);
}
fclose(RESULT_SUMMARY);
return engulf;	
}

double calc_exo_probs(){
char text[1000], t1[100], t2[100], t3[100], temp[1000], metallicity[1000], gaidos_py_callA[1000], gaidos_py_callR[1000],planet_pop[1000];
float Z, Prob_Giant, f0,a,Mstar, FEH;
printf("Now calculating probability of stars having a giant planet according to the specified mean Z\n\n");
FILE *INPUT = fopen("./configuration", "r"); 
while(!feof(INPUT))
		{
			fscanf(INPUT,"%[^\n]%*c\n", &text);
				if(text[0] == '#') {
					strcpy(text, "This is comment");}
		
			strcpy(temp,text);
			sscanf(temp, "%s %s %s", &t1, &t2, &t3);
																						
																																																																				
						if (strcmp(t1,"mean_metallicity") == 0){
							strcpy(metallicity,t3);	}	
						if (strcmp(t1,"planet_sample") == 0){
							strcpy(planet_pop,t3);	}	
																								
}

Z = strtof(metallicity, NULL);
printf("%f\n", Z);
if (strcmp(planet_pop,"exoplanets_org")==0){
sprintf(gaidos_py_callA, "python3 ./exo_org_analysis_scripts/gaidosAGB.py %f", Z);
sprintf(gaidos_py_callR, "python3 ./exo_org_analysis_scripts/gaidosRGB.py %f", Z);
system(gaidos_py_callA);
system(gaidos_py_callR);
}

if (strcmp(planet_pop,"IDA_2017")==0){
sprintf(gaidos_py_callA, "python3 ./IDA_data_analysis_scripts/gaidosAGB.py %f", Z);
sprintf(gaidos_py_callR, "python3 ./IDA_data_analysis_scripts/gaidosRGB.py %f", Z);
system(gaidos_py_callA);
system(gaidos_py_callR);
}
//convert Z to FEH
//Loop through star masses
//Calculate and then save to file

//If its a multi epoch simulation then the metallicity may change in each time step. When buliding that part of the code need to make sure this can be implemented in this function
//Gaidos will be calculated for each timestep in the galactic evolution according to metallicity.

}

double Multi_planets(char *sample){
int engulf_point,i, j, engulf, count_R_star, disk_num;
char text[1000], t1[100], t2[100], t3[100], temp[1000], python_ranges_call[100], Pathto_stellar_mods[500], mesa_name[100], new_mesa_path[1000], radius_plot[500], phase_check[500], age_plot[500];
char Name_of_Run[200], radvage[500], stellar_model[20], Runge_Kutta[100], plotting[20], M_planet[20], initial_a[20], M_star[20], planet_pop[50], exoplanets[50], stellar_metallicity[50];
float phase_star;
char results_directory[1000], results_file[1000],metallicity_handle[100], mean_Z[100];
double a0, increment_a, R_star, age_star;
float Z;
count_R_star = 0;
FILE *INPUT = fopen("./configuration", "r"); 
while(!feof(INPUT)){
	fscanf(INPUT,"%[^\n]%*c\n", &text);
		if(text[0] == '#') {
			strcpy(text, "This is comment");}
			strcpy(temp,text);
			sscanf(temp, "%s %s %s", &t1, &t2, &t3);
			if (strcmp(t1,"name_your_experiment") == 0)
				{
				strcpy(Name_of_Run,t3);								
				}																
			if (strcmp(t1,"star_models") == 0)
			{
				strcpy(stellar_model,t3);						
			}
			if (strcmp(t1,"runge_kutta") == 0)
				{
				strcpy(Runge_Kutta,t3);						
				}	
			if (strcmp(t1,"make_plot") == 0){
				strcpy(plotting,t3);	}													
			if (strcmp(t1,"planet_sample") == 0)
			{	strcpy(planet_pop,t3);						
			}
			if(strcmp(t1,"name_your_mesa_dir")==0){
				strcpy(mesa_name,t3);				
			}
			if(strcmp(t1,"metallicity_handle")==0){
			strcpy(metallicity_handle,t3);				
			}
			if(strcmp(t1,"mean_metallicity")==0){
				strcpy(mean_Z,t3);				
			}						
}

sprintf(results_directory, "./results_and_plots/%s_%spop_results", planet_pop, sample);
mkdir(results_directory, 0777);
sprintf(results_file, "%s/%s_%sMULTIS.csv", results_directory, planet_pop, sample);

FILE *pop_engulf = fopen(results_file, "w");

if (strcmp(planet_pop,"exoplanets_org")==0){
		strcpy(exoplanets, "../../Planets_Database/exo_org/Multi.csv");

		//system("python3 exo.py");
		printf("%s\n", exoplanets);
	}
if (strcmp(planet_pop,"IDA_2017")==0){
		printf("You will be calculating engulfments for IDA 2017\n");
		sprintf(exoplanets, "../../Planets_Database/IDA_2017/%s_planets_Multis.dat", sample);
		printf("Use specified file: %s\n", exoplanets);	
	}



	
FILE *PLANETS = fopen(exoplanets, "r"); 
	while(!feof(PLANETS))
		{
		
		count_R_star = 0;
		
		if (strcmp(planet_pop,"exoplanets_org")==0){
		fscanf(PLANETS, "%s %s %lf %f\n", &M_star, &M_planet, &a0, &Z);
		}
		if (strcmp(planet_pop,"IDA_2017")==0){
		fscanf(PLANETS, "%s %s %lf %f %d\n", &M_star, &M_planet, &a0, &Z, &disk_num);
		printf("%d\n", disk_num);
		}
		
		//
		
		
		if(Z < 0.005){
				Z = 0.005;
				strcpy(stellar_metallicity, "0.005");
		}
		if(Z> 0.004 && Z <= 0.0065){
				Z = 0.005;
				strcpy(stellar_metallicity, "0.005");
		}
		if(Z > 0.0065 && Z <= 0.01){
				Z = 0.008;	
				strcpy(stellar_metallicity, "0.008");
		}
		if(Z > 0.01 && Z <= 0.015){
				Z = 0.014;	
				strcpy(stellar_metallicity, "0.014");
		}
		if(Z > 0.015 && Z <= 0.03){
				Z = 0.025;
				strcpy(stellar_metallicity, "0.025");
				}
		if(Z > 0.03 && Z <= 0.05){
				Z = 0.04;
				strcpy(stellar_metallicity, "0.04");
				}

		if(strcmp(metallicity_handle,"mean")==0){
		Z = strtof(mean_Z, NULL);
		strcpy(stellar_metallicity, mean_Z);
		printf("%s and %f\n", stellar_metallicity, Z);
		}	
		
		sprintf(Pathto_stellar_mods,"../../Stellar_Database/%s/%sz_stars/%s_%s", stellar_model, stellar_metallicity, M_star, stellar_metallicity);
		sprintf(radius_plot, "%s/radius.txt", Pathto_stellar_mods);
		sprintf(age_plot, "%s/age.txt", Pathto_stellar_mods);
		sprintf(phase_check, "%s/phases.txt", Pathto_stellar_mods);
		sprintf(radvage, "./results_and_plots/%s_results/%sRvA.txt", Name_of_Run, Name_of_Run);

		j=1;

SMA_EVOL("interp_yes", a0, M_star, stellar_metallicity,Name_of_Run, stellar_model, Runge_Kutta, plotting, M_planet, Pathto_stellar_mods);
engulf=SMA_EVOL("interp_no", a0, M_star, stellar_metallicity,Name_of_Run, stellar_model, Runge_Kutta, plotting, M_planet, Pathto_stellar_mods);

if (strcmp(planet_pop,"exoplanets_org")==0){
fprintf(pop_engulf, "%s %s %f %s %d\n", M_star, M_planet, a0, stellar_metallicity, engulf);
		}
		
if (strcmp(planet_pop,"IDA_2017")==0){
fprintf(pop_engulf, "%s %s %f %s %d %d\n", M_star, M_planet, a0, stellar_metallicity, disk_num, engulf);
		}		
 
FILE *RAD = fopen(radius_plot, "r");
FILE *T = fopen(age_plot, "r");
FILE *PHASE = fopen(phase_check, "r");	
FILE *RVA = fopen(radvage, "w");
	while(!feof(RAD))
		{
		fscanf(RAD, "%lf\n", &R_star);
		fscanf(T, "%lf\n", &age_star);
		fscanf(PHASE, "%f\n", &phase_star);
		R_star = pow(10, R_star);
		if(phase_star > 1.0 && phase_star < 7.0){
		fprintf(RVA, "%lf %lf\n", age_star, R_star);
		count_R_star ++;
		}
		}
		fclose(RVA); 

if(strcmp(plotting, "yes")==0){
increment_a = 0;
plot_fun(j, Name_of_Run, M_star, a0, increment_a, Pathto_stellar_mods, count_R_star);
sprintf(python_ranges_call, "python3 ./plot_scripts/plot_orbital_evolution.py");		
system(python_ranges_call);	
}

fclose(PHASE);
fclose(T);
fclose(RAD);
		} 

fclose(pop_engulf);

}

char * create_new_mesa(){
char text[1000], t1[1000], t2[1000], t3[1000], temp[1000];
char mesa_name[1000];
char new_mesa_dir[1000], new_mesa_path[1000], run_mesa[1000], create_inlist[1000], mesa_controls[1000], add_massloss_controls[1000];
char Z[100], alpha[100], Mstar[100], reimers[10], blocker[10], rotation[10];
FILE *MESA_CONF = fopen("./configuration", "r");

while(!feof(MESA_CONF)){
fscanf(MESA_CONF,"%[^\n]%*c\n", &text);
				if(text[0] == '#') {
					strcpy(text, "This is comment");}
		
			strcpy(temp,text);
			sscanf(temp, "%s %s %s", &t1, &t2, &t3);
			if(strcmp(t1,"name_your_mesa_dir")==0){
				strcpy(mesa_name,t3);				
			}
			if(strcmp(t1,"star_mass")==0){
				strcpy(Mstar,t3);
				printf("%s\n",Mstar);				
			}
			if(strcmp(t1,"star_metallicity")==0){
				strcpy(Z,t3);
				printf("%s\n",Z);				
			}
			if(strcmp(t1,"mixing_length_alpha")==0){
				strcpy(alpha,t3);
				printf("%s\n",alpha);				
			}
			if(strcmp(t1,"Reimers_RGB")==0){
				strcpy(reimers,t3);
				printf("%s\n",reimers);				
			}
			if(strcmp(t1,"Blocker_AGB")==0){
				strcpy(blocker,t3);
				printf("%s\n",blocker);				
			}
			if(strcmp(t1,"rotation")==0){
				strcpy(rotation,t3);
				printf("%s\n",rotation);				
			}
			
			
}
sprintf(mesa_controls, "../../Stellar_database/CREATE_NEW_MESA/%s/mesa_controls.txt", mesa_name);
sprintf(new_mesa_dir, "cp -r ../../Stellar_Database/CREATE_NEW_MESA/Work ../../Stellar_database/CREATE_NEW_MESA/%s", mesa_name);


system(new_mesa_dir);
FILE *MESA_CONTROLS = fopen(mesa_controls, "w"); 


if(strcmp(rotation,"ON")==0){
fprintf(MESA_CONTROLS, " change_initial_rotation_flag = .true. ! change state of rotation flag on run start\n");
fprintf(MESA_CONTROLS, " new_rotation_flag = .true. !new state = rotation ON\n");

fprintf(MESA_CONTROLS, "set_initial_surface_rotation_v = .true. ! set solid-body rotation using surface v \n");
fprintf(MESA_CONTROLS, "new_surface_rotation_v = 2 ! specify surface v. in km/s\nn");

}

fprintf(MESA_CONTROLS, " / !end of star_job namelist\n\n\n&controls\n\n\n! starting specifications\n");
fprintf(MESA_CONTROLS, "\ninitial_mass = %s\n", Mstar);
fprintf(MESA_CONTROLS, "initial_z = %s\n\n", Z);
fprintf(MESA_CONTROLS, "\nmixing_length_alpha = %s\n\n", alpha);


fclose(MESA_CONTROLS);

//EDIT inlist here before running
sprintf(run_mesa, "cd ../../Stellar_database/CREATE_NEW_MESA/%s ; ./mk; ./rn", mesa_name);




///Copy work directory in new mesa diirectory to the name given in config file. This file will now contain everything needed to run the default mesa 
/// Then edit the original inlist for the star that needs to be evolved according to user spec. Then run it.
/// Then copy the history file into a new directory for use by code. Tell simsplash where that file is. 
sprintf(create_inlist, "cat ../../Stellar_database/CREATE_NEW_MESA/%s/make_inlist_starjob.txt ../../Stellar_database/CREATE_NEW_MESA/%s/mesa_controls.txt ../../Stellar_database/CREATE_NEW_MESA/%s/make_inlist_controls.txt > ../../Stellar_database/CREATE_NEW_MESA/%s/inlist_project", mesa_name, mesa_name, mesa_name, mesa_name);

system(create_inlist);


sprintf(add_massloss_controls, "../../Stellar_database/CREATE_NEW_MESA/%s/inlist_project", mesa_name);
FILE *MESA_masslossCONTROLS = fopen(add_massloss_controls, "a+"); 

fprintf(MESA_masslossCONTROLS, "\nReimers_wind_eta = %s\n", reimers);
fprintf(MESA_masslossCONTROLS, "Blocker_wind_eta = %s\n", blocker);
fprintf(MESA_masslossCONTROLS, "/ ! end of controls namelist\n");
//Add a stopping condition in case of convergence issues
fclose(MESA_masslossCONTROLS);
system(run_mesa);
}

int welcomeNote(){

printf("\n!******************* SIMSPLASH (SIMulationS for the PLAnet Shaping Hypothesis) ********************!");
printf("\n!*********************************     Boyle & Redman (2017)     **********************************!\n\n\n");

}

double calculate_turnoff_mass(double pop_age, float mean_Z){

	char text[1000], t1[100], t2[100], t3[100], temp[1000];
    char pop_mass[100], exo_pop[100], Z[100], age[100];	
	double a0, a1, a2, a, b, c;
	double turnoff_mass;



FILE *INPUT = fopen("./configuration", "r"); 
		while(!feof(INPUT))
		{	fscanf(INPUT,"%[^\n]%*c\n", &text);
				if(text[0] == '#') {
					strcpy(text, "This is comment");}		
			strcpy(temp,text);
			sscanf(temp, "%s %s %s", &t1, &t2, &t3);							
					
					
						
							
						}
//printf("The age of the population is: %g\n", pop_age);
//printf("The mean metallicity is: %g\n", mean_Z);

a0 = 10.13 + (0.07547*log10(mean_Z)) - (0.008084*log10(mean_Z)*log10(mean_Z));
a1 = -4.424 - (0.7939*log10(mean_Z)) - (0.1187*log10(mean_Z)*log10(mean_Z));
a2 = 1.262 + (0.3385*log10(mean_Z)) + (0.05417*log10(mean_Z)*log(mean_Z));
c = a0 - log10(pop_age);
b=a1;
a=a2;

turnoff_mass = pow(10,((-1*b) - sqrt((b*b)- (4*a*c))) / (2*a));
return turnoff_mass;
}

double calculate_multi_population(){
IMF imf;
int i,j, count;
char t1[1000],t2[1000],t3[1000],text[1000], temp[1000], name_dir[1000], resultsTotal_dir[1000], resultsVisible_dir[1000], results_dir[1000];
char Call_GaidosAGB[1000], Call_GaidosRGB[1000], timebin_label[1000],mass_loss[100], engulfment_probs_RGB[1000], engulfment_probs_AGB[1000];
double ti, ti_1, mi, mi_1, population_mass, feh;
double N_stars_total, stellar_mass_total;
double N_PN_progs, N_PN_today, N_PN_TODAY;
double N_B, N_S, N_PN_SYSTEMS;
double mean_mass, mean_mass_total;

double mj, mj_1, metalZ;
double binary_f_RGB, binary_f_AGB, PN_prog_mass_low, PN_avg_vis_time;
double AGB_ENGULFMENT_PROBABILITY, RGB_ENGULFMENT_PROBABILITY;
double Nj, Nj_B, Nj_S, Nj_SYSTEMS;
double Nj_PL_AGB, Nj_PL_RGB, Nj_PL_SS;
double Nj_B_AGB, Nj_B_RGB, Nj_B_SS;
double NPN_total, NPN_total_PLRGB, N_PL_AGB_total, N_B_AGB_total, N_SS_total, N_PL_RGB_total, N_B_RGB_total;
double N_PL_AGB_TODAY, N_B_AGB_TODAY, N_SS_TODAY,N_PL_RGB_TODAY ;
double test_pn_total,N_PN_TODAY_PLRGB;
double junk1, junk2;
char test[1000], test2[1000], B_f_RGB[1000], B_f_AGB[1000], prog_mass_low[1000], RGB_pl_prevents[1000], IMF_form[1000], PN_vis[1000];

FILE *INPUT = fopen("./configuration", "r"); 
		while(!feof(INPUT))
		{	fscanf(INPUT,"%[^\n]%*c\n", &text);
				if(text[0] == '#') {
					strcpy(text, "This is comment");}		
			strcpy(temp,text);
			sscanf(temp, "%s %s %s", &t1, &t2, &t3);
			if (strcmp(t1,"name_your_experiment") == 0)
						{sprintf(name_dir,"mkdir ./results_and_plots/%s", t3);
						system(name_dir);
						sprintf(resultsTotal_dir, "./results_and_plots/%s/Timebin_resultsTotal.dat", t3);
						sprintf(resultsVisible_dir, "./results_and_plots/%s/Timebin_resultsVisible.dat", t3);
						sprintf(results_dir, "./results_and_plots/%s/Multi_Epoch_Result.dat", t3);
						}
			if (strcmp(t1,"mass_loss_prescription") == 0)
						{strcpy(mass_loss,t3);
						}
			if (strcmp(t1,"binary_f_RGB") == 0)
						{strcpy(B_f_RGB,t3);
						binary_f_RGB = strtof(B_f_RGB, NULL);
						
						}			
			if (strcmp(t1,"binary_f_AGB") == 0)
						{strcpy(B_f_AGB,t3);
						binary_f_AGB = strtof(B_f_AGB, NULL);
						
						}	
			if (strcmp(t1,"progenitor_mass_low") == 0)
						{strcpy(prog_mass_low,t3);
						PN_prog_mass_low = strtof(prog_mass_low, NULL);
				
						}
			if (strcmp(t1,"RGB_pl_prevents") == 0)
						{strcpy(RGB_pl_prevents,t3);
					
						}	
			if (strcmp(t1,"IMF") == 0)
						{strcpy(IMF_form,t3);
						}	
			if (strcmp(t1,"PN_avg_visibility") == 0)
						{strcpy(PN_vis,t3);
						PN_avg_vis_time = strtof(PN_vis, NULL);
					
						}			
			}
if (strcmp(IMF_form,"kroupa_01") == 0){
imf_init_user(&imf,"0.08(pow:-1.3)0.5(pow:-2.3)120");
}
if (strcmp(IMF_form,"chabrier_03") == 0){
imf_init_chabrier_2003(&imf);
}
			
FILE *POPULATIONS = fopen("../../Star_Formation_History/POPULATIONS.dat", "r");
FILE *POPULATION_TOTAL_RESULTS = fopen(resultsTotal_dir, "w");
FILE *POPULATION_VISIBLE_RESULTS = fopen(resultsVisible_dir, "w");
FILE *MULTI_EPOCH_RESULTS = fopen(results_dir, "w");



fprintf(POPULATION_TOTAL_RESULTS, "ti, ti_1, mi, mi_1, feh, population_mass, N_stars_total, N_PN_SYSTEMS,N_PL_AGB_total, N_B_AGB_total, N_SS_total, N_PL_RGB_total, N_B_RGB_total\n"); 
fprintf(POPULATION_VISIBLE_RESULTS, "ti ti_1 mi mi_1 FEH PopMass NStarsTotal NPN_visible NPL_AGBVisible NPL_RGBVisible NBin_AGBTotal NBin_AGBVisible NSingleTotal\n"); 

fprintf(MULTI_EPOCH_RESULTS, "The following is the resulting PN population currently visible in the Milky Way\n\n");
count=0;
N_PN_TODAY = 0;
mean_mass = 0;
mean_mass_total = 0;
while(!feof(POPULATIONS))
{
count ++;
fscanf(POPULATIONS, "%lf %lf %lf %lf %lf %lf\n", &ti, &ti_1, &mi, &mi_1, &population_mass, &feh);

imf_norm_cl(&imf,population_mass,120);
N_stars_total = imf_int_xi_cl(&imf,0.08,120.0); //save this
stellar_mass_total = imf_int_mxi_cl(&imf,0.08,120.0);



if(mi >= 0.8 && mi < 8.0){ 
N_PN_progs = imf_int_xi_cl(&imf,mi,mi_1);
N_B = 0.67*N_PN_progs; //number of binaries from 50% progenitor fraction
N_S = 0.33*N_PN_progs;
N_PN_SYSTEMS = (N_B/2) + N_S; //Total number of PN progenitor SYSTEMS in this timebin
metalZ = pow(10, feh) * 0.014;
test_pn_total=0;
mj = mi;
mj_1 = mi+0.1;
if(fabs((mi_1*10) - (mi*10)) < 1.0){
mj = mi;
mj_1 = mi_1;
}

for(j=(mi*10); j<((mi_1*10)); j++){


if (strcmp(mass_loss,"reduced") == 0){
if(metalZ < 0.005){
metalZ = 0.005;
}
sprintf(test, "grep -i '%.01f %.03f' ../../Star_Formation_History/all_interpolated_probabilitiesAGBIDA.csv", mi, metalZ);
FILE *fp = popen(test, "r");
fscanf(fp, "%lf %lf %lf", &junk1, &junk2, &AGB_ENGULFMENT_PROBABILITY);
pclose(fp);
sprintf(test2, "grep -i '%.01f %.03f' ../../Star_Formation_History/all_interpolated_probabilitiesRGBIDA.csv", mi, metalZ);
FILE *fp2 = popen(test2, "r");
fscanf(fp2, "%lf %lf %lf", &junk1, &junk2, &RGB_ENGULFMENT_PROBABILITY);
pclose(fp2);
}

if (strcmp(mass_loss,"standard") == 0){
if(metalZ < 0.010){
metalZ = 0.010;
}
sprintf(test, "grep -i '%.01f %.03f' ../../Star_Formation_History/all_interpolated_probabilitiesAGBMESA.csv", mi, metalZ);
FILE *fp = popen(test, "r");
fscanf(fp, "%lf %lf %lf", &junk1, &junk2, &AGB_ENGULFMENT_PROBABILITY);

pclose(fp);
sprintf(test2, "grep -i '%.01f %.03f' ../../Star_Formation_History/all_interpolated_probabilitiesRGBMESA.csv", mi, metalZ);
FILE *fp2 = popen(test2, "r");
fscanf(fp2, "%lf %lf %lf", &junk1, &junk2, &RGB_ENGULFMENT_PROBABILITY);
pclose(fp2);


}


Nj = imf_int_xi_cl(&imf,mj,mj_1);
Nj_B = 0.67*Nj;
Nj_S = 0.33*Nj;
Nj_SYSTEMS = (Nj_B/2) + Nj_S;


Nj_B_AGB = binary_f_AGB*(Nj_B/2);
Nj_B_RGB = binary_f_RGB*(Nj_B/2);
Nj_B_SS = (1- (binary_f_RGB+binary_f_AGB))*(Nj_B/2);

Nj_PL_AGB = AGB_ENGULFMENT_PROBABILITY*Nj_S;
Nj_PL_RGB = RGB_ENGULFMENT_PROBABILITY*Nj_S;
if (strcmp(RGB_pl_prevents,"yes") == 0){
Nj_PL_SS = Nj_S - Nj_PL_AGB - Nj_PL_RGB;
}

if (strcmp(RGB_pl_prevents,"no") == 0){
Nj_PL_SS = Nj_S - Nj_PL_AGB;
}

if(mi<PN_prog_mass_low){
Nj_B_SS = 0;
Nj_PL_SS = 0;
}

//NPN_total = NPN_total+(Nj_B_AGB + Nj_PL_AGB + Nj_B_SS + Nj_PL_SS + Nj_PL_RGB);
NPN_total = NPN_total+(Nj_B_AGB + Nj_PL_AGB + Nj_B_SS + Nj_PL_SS);
mean_mass = mean_mass + (mj * ((Nj_B_AGB + Nj_PL_AGB + Nj_B_SS + Nj_PL_SS + Nj_PL_RGB)*(26000/(fabs(ti_1-ti))))); 
NPN_total_PLRGB = NPN_total_PLRGB+(Nj_B_AGB + Nj_PL_AGB + Nj_B_SS + Nj_PL_SS + Nj_PL_RGB);


N_PL_AGB_total = N_PL_AGB_total + Nj_PL_AGB;

N_B_AGB_total = N_B_AGB_total + Nj_B_AGB;

N_SS_total = N_SS_total + (Nj_PL_SS + Nj_B_SS);

N_PL_RGB_total = N_PL_RGB_total + Nj_PL_RGB;
N_B_RGB_total = N_B_RGB_total + Nj_B_RGB;


mj = mj_1;
mj_1 = mj+0.1;
if(mj_1 > mi_1){
mj_1 = mi_1;
}
}

fprintf(POPULATION_TOTAL_RESULTS,"%g %g %g %g %g %g %g %g %g %g %g %g %g\n", ti, ti_1, mi, mi_1, feh, population_mass, N_stars_total, N_PN_SYSTEMS,N_PL_AGB_total, N_B_AGB_total, N_SS_total, N_PL_RGB_total, N_B_RGB_total);


NPN_total = NPN_total *(PN_avg_vis_time/(fabs(ti_1-ti)));
NPN_total_PLRGB = NPN_total_PLRGB *(PN_avg_vis_time/(fabs(ti_1-ti)));
N_PL_AGB_total = N_PL_AGB_total*(PN_avg_vis_time/(fabs(ti_1-ti)));
N_B_AGB_total = N_B_AGB_total *(PN_avg_vis_time/(fabs(ti_1-ti)));
N_SS_total = N_SS_total*(PN_avg_vis_time/(fabs(ti_1-ti)));
N_PL_RGB_total = N_PL_RGB_total *(PN_avg_vis_time/(fabs(ti_1-ti)));
N_B_RGB_total = N_B_RGB_total*(PN_avg_vis_time/(fabs(ti_1-ti)));
fprintf(POPULATION_VISIBLE_RESULTS,"%g %g %g %g %g %g %g %g %g %g %g %g %g\n", ti, ti_1, mi, mi_1, feh, population_mass, N_stars_total, NPN_total,N_PL_AGB_total, N_B_AGB_total, N_SS_total, N_PL_RGB_total, N_B_RGB_total);

mean_mass_total = mean_mass_total + mean_mass;

N_PN_TODAY = N_PN_TODAY + NPN_total;

N_PN_TODAY_PLRGB = N_PN_TODAY_PLRGB + NPN_total_PLRGB;
N_PL_AGB_TODAY = N_PL_AGB_TODAY + N_PL_AGB_total;
N_PL_RGB_TODAY = N_PL_RGB_TODAY + N_PL_RGB_total;
N_B_AGB_TODAY = N_B_AGB_TODAY + N_B_AGB_total;
N_SS_TODAY = N_SS_TODAY + N_SS_total;

test_pn_total=0;
NPN_total=0;
NPN_total_PLRGB = 0;
N_PL_AGB_total=0;
N_B_AGB_total=0;
N_SS_total=0;
N_PL_RGB_total=0;
N_B_RGB_total=0;
mean_mass=0;

}
}

mean_mass_total = mean_mass_total/N_PN_TODAY;
printf("%lf\n", mean_mass_total);
fprintf(MULTI_EPOCH_RESULTS,"\n\n Visible PNe today: %g\n", N_PN_TODAY);
fprintf(MULTI_EPOCH_RESULTS,"\n\n Visible PNe today (plus RGB engulfment): %g\n", N_PN_TODAY_PLRGB);
fprintf(MULTI_EPOCH_RESULTS,"\n\n Visible binary PNe today: %g\n", N_B_AGB_TODAY);
fprintf(MULTI_EPOCH_RESULTS,"\n\n Visible PLANET PNe today: %g\n", N_PL_AGB_TODAY);
fprintf(MULTI_EPOCH_RESULTS,"\n\n Visible single star.spherical today: %g\n\n", N_SS_TODAY);
printf("%g %g\n", N_PL_AGB_TODAY, N_PL_RGB_TODAY);
fprintf(MULTI_EPOCH_RESULTS,"\n\n Binary Fraction: %g\n", N_B_AGB_TODAY/N_PN_TODAY);
fprintf(MULTI_EPOCH_RESULTS,"\n\n PLANET FRACTION: %g\n", N_PL_AGB_TODAY/N_PN_TODAY);
}

double Multi_epoch(){

double redshift, age_univ, age_population, log_sfr, feh;
double ti, ti_1, timebin;
double SFR, pop_mass;
double mi, mi_1;
double mean_Zi, mean_Zi_1;

printf("You are now in multi-epoch population mode\n");
system("python3 ./SFH_scripts/SFH.py");

FILE *SFH = fopen("../../Star_Formation_History/SFH_DATA.dat", "r");
FILE *POPULATIONS = fopen("../../Star_Formation_History/POPULATIONS.dat", "w");
ti = 0;
ti_1 = 0;

mean_Zi = 0;
mean_Zi_1 = 0;
while(!feof(SFH)){

fscanf(SFH, "%lf %lf %lf %lf %lf\n",&redshift, &age_univ, &age_population, &log_sfr, &feh);

ti_1 = age_population*1e9;
mean_Zi_1 = pow(10, feh) * 0.014;
if(ti == 0){
log_sfr = 0;
}

timebin = fabs(ti_1 - ti);
SFR = pow(10, log_sfr);
pop_mass = SFR * timebin;
//printf("%g %g %g %g\n", ti, ti_1, mean_Zi, mean_Zi_1);
mi = calculate_turnoff_mass(ti,mean_Zi);
mi_1 = calculate_turnoff_mass(ti_1,mean_Zi_1);

if (mi_1 > 8.0){
mi_1 = 8.0;
}
if(mi > 8.0){
mi = 8.0;
}

if(ti !=0 && ti_1 !=0){
fprintf(POPULATIONS, "%g %g %g %g %g %g\n", ti, ti_1, mi, mi_1, pop_mass, feh); 
}


ti = ti_1;
mean_Zi = mean_Zi_1;
}

fclose(POPULATIONS);
calculate_multi_population();
}

double Single_epoch(char *sample){
int engulf_point,i, j, engulf, count_R_star;
char text[1000], t1[100], t2[100], t3[100], temp[1000], python_ranges_call[100], Pathto_stellar_mods[500], mesa_name[100], new_mesa_path[1000], radius_plot[500], phase_check[500], age_plot[500];
char Name_of_Run[200], radvage[500], stellar_model[20], Runge_Kutta[100], plotting[20], M_planet[20], initial_a[20], M_star[20], planet_pop[50], exoplanets[50], stellar_metallicity[50];
char user_specified_planets[1000],results_directory[1000], results_file[1000], metallicity_handle[100], mean_Z[100];
float phase_star;
double a0, increment_a, R_star, age_star;
float Z;
count_R_star = 0;
FILE *engulfed = fopen("./engulfed.txt", "w"); 
FILE *INPUT = fopen("./configuration", "r"); 

while(!feof(INPUT))
	{
	fscanf(INPUT,"%[^\n]%*c\n", &text);
	if(text[0] == '#') {
		strcpy(text, "This is comment");}		
		strcpy(temp,text);
		sscanf(temp, "%s %s %s", &t1, &t2, &t3);
		if (strcmp(t1,"name_your_experiment") == 0)
		{
			strcpy(Name_of_Run,t3);								
		}																
		if (strcmp(t1,"star_models") == 0)
		{
			strcpy(stellar_model,t3);						
		}	
		if (strcmp(t1,"runge_kutta") == 0)
		{
			strcpy(Runge_Kutta,t3);						
		}		
		if (strcmp(t1,"make_plot") == 0){
			strcpy(plotting,t3);	}													
		if (strcmp(t1,"planet_sample") == 0)
			{strcpy(planet_pop,t3);						
		}
		if(strcmp(t1,"name_your_mesa_dir")==0){
			strcpy(mesa_name,t3);				
		}	
		if(strcmp(t1,"metallicity_handle")==0){
			strcpy(metallicity_handle,t3);				
		}
		if(strcmp(t1,"mean_metallicity")==0){
			strcpy(mean_Z,t3);				
		}						
}

sprintf(results_directory, "./results_and_plots/%s_%spop_results", planet_pop, sample);
mkdir(results_directory, 0777);
sprintf(results_file, "%s/%s_%sSingle.csv", results_directory, planet_pop, sample);

FILE *pop_engulf = fopen(results_file, "w");

if (strcmp(planet_pop,"exoplanets_org")==0){
		strcpy(exoplanets, "../../Planets_Database/exo_org/exo_org_cutdata.csv");
		system("python3 ./exo_org_analysis_scripts/exo.py");
		printf("%s\n", exoplanets);
	}
	
if (strcmp(planet_pop,"IDA_2017")==0){
		printf("You will be calculating engulfments for IDA 2017\n");
		sprintf(exoplanets, "../../Planets_Database/IDA_2017/%s_planets_Single.dat", sample);
		printf("Use specified file: %s\n", exoplanets);	
	} 	
	
if (strcmp(planet_pop,"exoplanets_org")!=0 && strcmp(planet_pop,"IDA_2017")!=0){		
	printf("Use specified file: %s\n", planet_pop);	
	sprintf(user_specified_planets, "../../Planets_Database/User_Specified/%s", planet_pop);
	strcpy(exoplanets, user_specified_planets);
}
	
FILE *PLANETS = fopen(exoplanets, "r"); 
	while(!feof(PLANETS))
		{
		count_R_star = 0;
		//scan in and call orbital models
		fscanf(PLANETS, "%s %s %lf %f\n", &M_star, &M_planet, &a0, &Z);
	
		if(Z < 0.005){
				Z = 0.005;
				strcpy(stellar_metallicity, "0.005");
		}
		if(Z> 0.004 && Z <= 0.0065){
				Z = 0.005;
				strcpy(stellar_metallicity, "0.005");
		}
		if(Z > 0.0065 && Z <= 0.01){
				Z = 0.008;	
				strcpy(stellar_metallicity, "0.008");
		}
		if(Z > 0.01 && Z <= 0.015){
				Z = 0.014;	
				strcpy(stellar_metallicity, "0.014");
		}
		if(Z > 0.015 && Z <= 0.03){
				Z = 0.025;
				strcpy(stellar_metallicity, "0.025");
				}
		if(Z > 0.03 && Z <= 0.05){
				Z = 0.04;
				strcpy(stellar_metallicity, "0.04");
				}
		if(strcmp(metallicity_handle,"mean")==0){
		Z = strtof(mean_Z, NULL);
		strcpy(stellar_metallicity, mean_Z);
		printf("%s and %f\n", stellar_metallicity, Z);
		}
		sprintf(Pathto_stellar_mods,"../../Stellar_Database/%s/%sz_stars/%s_%s", stellar_model, stellar_metallicity, M_star, stellar_metallicity);
		sprintf(radius_plot, "%s/radius.txt", Pathto_stellar_mods);
		sprintf(age_plot, "%s/age.txt", Pathto_stellar_mods);
		sprintf(phase_check, "%s/phases.txt", Pathto_stellar_mods);
		sprintf(radvage, "./results_and_plots/%s_results/%sRvA.txt", Name_of_Run, Name_of_Run);
		
		j=1;

SMA_EVOL("interp_yes", a0, M_star, stellar_metallicity,Name_of_Run, stellar_model, Runge_Kutta, plotting, M_planet, Pathto_stellar_mods);
engulf=SMA_EVOL("interp_no", a0, M_star, stellar_metallicity,Name_of_Run, stellar_model, Runge_Kutta, plotting, M_planet, Pathto_stellar_mods);
fprintf(pop_engulf, "%s %s %f %f %d\n", M_star, M_planet, a0, Z, engulf); 
FILE *RAD = fopen(radius_plot, "r");
FILE *T = fopen(age_plot, "r");
FILE *PHASE = fopen(phase_check, "r");	
FILE *RVA = fopen(radvage, "w");
	while(!feof(RAD))
		{
		fscanf(RAD, "%lf\n", &R_star);
		fscanf(T, "%lf\n", &age_star);
		fscanf(PHASE, "%f\n", &phase_star);
		R_star = pow(10, R_star);
		if(phase_star > 1.0 && phase_star < 7.0){
		fprintf(RVA, "%lf %lf\n", age_star, R_star);
		count_R_star ++;
		}
		}
		fclose(RVA); 
fclose(engulfed);
if(strcmp(plotting, "yes")==0){
increment_a = 0;
plot_fun(j, Name_of_Run, M_star, a0, increment_a, Pathto_stellar_mods, count_R_star);
sprintf(python_ranges_call, "python3 ./plot_scripts/plot_orbital_evolution.py");		
system(python_ranges_call);	
}
		} 

fclose(pop_engulf);
}

double Orbital(){
int engulf_point,i, j, engulf, count_R_star;
char text[1000], t1[100], t2[100], t3[100], temp[1000], Pathto_stellar_mods[500], mesa_name[100], new_mesa_path[1000], radius_plot[500], phase_check[500], age_plot[500];
char Name_of_Run[200], radvage[500], stellar_model[20], Runge_Kutta[100], plotting[20], M_planet[20], SMA[50], range_AU[50], initial_a[20], M_star[20], stellar_metallicity[20], planet_pop[50];
char cml_answer[10];
float quantity_of_a, phase_star;
double a0, low_a, high_a, increment_a, start_a, R_star, age_star;
count_R_star = 0;
FILE *engulfed = fopen("./engulfed.txt", "w"); 
FILE *INPUT = fopen("./configuration", "r"); 
while(!feof(INPUT)){
	fscanf(INPUT,"%[^\n]%*c\n", &text);
	if(text[0] == '#') {
		strcpy(text, "This is comment");}
		
	strcpy(temp,text);
	sscanf(temp, "%s %s %s", &t1, &t2, &t3);
	if (strcmp(t1,"name_your_experiment") == 0){
		strcpy(Name_of_Run,t3);								
	}
	if (strcmp(t1,"SMA_type") == 0){
		strcpy(SMA,t3);	
	}
	if (strcmp(t1,"SMA_range") == 0){
		strcpy(range_AU,t3);							
	}	
	if (strcmp(t1,"runge_kutta") == 0){
		strcpy(Runge_Kutta,t3);							
	}
	if (strcmp(t1,"initial_SMA") == 0){
		strcpy(initial_a,t3);
		a0 = strtof(initial_a, NULL);															
	}
	if (strcmp(t1,"star_mass") == 0){
		strcpy(M_star,t3);
	}
	if (strcmp(t1,"star_metallicity") == 0){
		strcpy(stellar_metallicity,t3);
	}																	
	if (strcmp(t1,"star_models") == 0){
		strcpy(stellar_model,t3);						
	}							
	if (strcmp(t1,"planet_mass") == 0){
		strcpy(M_planet,t3);								
	}																																									
	if (strcmp(t1,"make_plot") == 0){
		strcpy(plotting,t3);							
	}	
	if (strcmp(t1,"planet_sample") == 0){
		strcpy(planet_pop,t3);						
	}
	if(strcmp(t1,"name_your_mesa_dir")==0){
		strcpy(mesa_name,t3);
	}
							
}

printf("You have chosen the Orbital mode of SIMSPLASH!\n\n");

if (strcmp(stellar_model,"MIST_R")==0 || strcmp(stellar_model,"MIST_NR")==0){
	if (strcmp(stellar_metallicity,"0.005")==0 || strcmp(stellar_metallicity,"0.008")==0 || strcmp(stellar_metallicity,"0.014")==0 || strcmp(stellar_metallicity,"0.025")==0 || strcmp(stellar_metallicity,"0.04")==0){
	
	}
	else{
		printf("*\n*\n*\n*\n*\nATTENTION:: ERROR: The specified metallicity stellar models do not exist in the original version of SIMSPLASH :\(\nIf you have recently added new stellar models to the database you need to update the source code\n*\n*\n*\n*\n*\nHave you recently added the specified model to the stellar database? (y/n)?\nAns:");
		scanf("%s", &cml_answer);
			if (strcmp(cml_answer, "y") == 0){
			printf("*\n*\n*\n*\n*\nATTENTION: Unfortunately you must update the simsplash_src.c source code file in the Orbital function to include the new models");
			return 0;
			}
			if (strcmp(cml_answer, "n") == 0){
			printf("*\n*\n*\n*\n*\nOk. Would you like to evolve your systems with a similar stellar model? (y/n)\nAns:");
			scanf("%s", &cml_answer);
			}			
			if(strcmp(cml_answer, "y") == 0){
				printf("*\n*\n*\n*\n*\nPlease choose from the following available metallicities in the MIST stellar models in SIMSPLASH:\n 0.005, 0.008, 0.014, 0.025, 0.04\nAns:");
				scanf("%s", &stellar_metallicity);
				if (strcmp(stellar_metallicity,"0.005")==0 || strcmp(stellar_metallicity,"0.008")==0 || strcmp(stellar_metallicity,"0.014")==0 || strcmp(stellar_metallicity,"0.025")==0 || strcmp(stellar_metallicity,"0.04")==0){
							}
				else{
					printf("\n*\n*\n*\n*\n*\nFATAL ERROR: An invalid response was entered.\nPlease check the configuration file is specifying a valid stellar model and start again\n*\n*\n*\n*\n");
					return 0;
				}
				}			
			if(strcmp(cml_answer, "n") == 0){
				printf("\n*\n*\n*\n*\nOk, come back again soon!\n*\n*");
				return 0;				
				}						
		}
	}	
if (strcmp(stellar_model,"MESA_R")==0 || strcmp(stellar_model,"MESA_NR")==0){
	if (strcmp(stellar_metallicity,"0.01")==0 || strcmp(stellar_metallicity,"0.02")==0){	
	}
	else{
		printf("*\n*\n*\n*\n*\nATTENTION: The specified metallicity stellar models do not exist in the original version of SIMSPLASH :\(\nIf you have recently added new stellar models to the database you need to update the source code\n*\n*\n*\n*\n*\nHave you recently added the specified model to the stellar database? (y/n)?\nAns:");
		scanf("%s", &cml_answer);
			if (strcmp(cml_answer, "y") == 0){
			printf("*\n*\n*\n*\n*\nATTENTION: Unfortunately you must update the simsplash_src.c source code file in the Orbital function to include the new models");
			return 0;
			}
			if (strcmp(cml_answer, "n") == 0){
			printf("*\n*\n*\n*\n*\nOk. Would you like to evolve your systems with a similar stellar model? (y/n)\nAns:");
			scanf("%s", &cml_answer);
			}
			
			if(strcmp(cml_answer, "y") == 0){
				printf("*\n*\n*\n*\n*\nPlease choose from the following available metallicities in the MESA(From Boyle) stellar models in SIMSPLASH:\n 0.01 or 0.02\nAns:");
				scanf("%s", &stellar_metallicity);
				if (strcmp(stellar_metallicity,"0.01")==0 || strcmp(stellar_metallicity,"0.02")==0){
						
							}
				else{
				
				printf("\n*\n*\n*\n*\n*\nFATAL ERROR: An invalid response was entered.\nPlease check the configuration file is specifying a valid stellar model and start again\n*\n*\n*\n*\n");
				return 0;
				}
			}
			
			if(strcmp(cml_answer, "n") == 0){
			printf("\n*\n*\n*\n*\nOk, come back again soon!\n*\n*");
			return 0;
				
			}					
		}
	}		
printf("*\n*\n*\n*\n*\nPlease wait a few moments while we evolve your system(s)...\n*\n*\n*\n*\n*\n*\n");
if (strcmp(stellar_model,"CREATE_NEW_MESA")!=0){					
sprintf(Pathto_stellar_mods,"../../Stellar_Database/%s/%sz_stars/%s_%s", stellar_model, stellar_metallicity, M_star, stellar_metallicity);
}
sprintf(radius_plot, "%s/radius.txt", Pathto_stellar_mods);
sprintf(age_plot, "%s/age.txt", Pathto_stellar_mods);
sprintf(phase_check, "%s/phases.txt", Pathto_stellar_mods);
sprintf(radvage, "./results_and_plots/%s_results/%sRvA.txt", Name_of_Run, Name_of_Run);

if (strcmp(SMA,"range_sma") == 0)
{	
	sscanf(range_AU,"%4[^,],%[^,],%[^,]", &t1, &t2, &t3);								
	low_a = strtof(t1, NULL);			
	start_a = low_a;
	high_a = strtof(t2, NULL);
	increment_a = strtof(t3, NULL);		
	quantity_of_a = (high_a - low_a)/increment_a;								
	j = (int)quantity_of_a +1;									
	for(i = 0; i < j; i++){ 
		count_R_star = 0;
		a0 = low_a;
		printf("%lf\n", low_a);
		low_a = low_a + increment_a;
		SMA_EVOL("interp_yes", a0, M_star, stellar_metallicity, Name_of_Run, stellar_model, Runge_Kutta, plotting, M_planet, Pathto_stellar_mods);
		engulf = SMA_EVOL("interp_no", a0, M_star, stellar_metallicity,Name_of_Run, stellar_model, Runge_Kutta, plotting, M_planet, Pathto_stellar_mods);
		fprintf(engulfed, "%d\n", engulf);			
	}
fclose(engulfed);	
FILE *RAD = fopen(radius_plot, "r");
FILE *T = fopen(age_plot, "r");
FILE *PHASE = fopen(phase_check, "r");	
FILE *RVA = fopen(radvage, "w");

while(!feof(RAD)){
	fscanf(RAD, "%lf\n", &R_star);
	fscanf(T, "%lf\n", &age_star);
	fscanf(PHASE, "%f\n", &phase_star);
	R_star = pow(10, R_star);

	if(phase_star > 1.0 && phase_star < 7.0){
		fprintf(RVA, "%lf %lf\n", age_star, R_star);
		count_R_star ++;
	}		
}
	fclose(RVA);
	plot_fun(j, Name_of_Run, M_star, start_a, increment_a, Pathto_stellar_mods, count_R_star);			
}
if (strcmp(SMA,"single_sma") == 0){
	j=1;
	increment_a = 0;

	SMA_EVOL("interp_yes", a0, M_star, stellar_metallicity,Name_of_Run, stellar_model, Runge_Kutta, plotting, M_planet, Pathto_stellar_mods);
	engulf=SMA_EVOL("interp_no", a0, M_star, stellar_metallicity,Name_of_Run, stellar_model, Runge_Kutta, plotting, M_planet, Pathto_stellar_mods);
	fprintf(engulfed, "%d\n", engulf);
	
	FILE *RAD = fopen(radius_plot, "r");
	FILE *T = fopen(age_plot, "r");
	FILE *PHASE = fopen(phase_check, "r");	
	FILE *RVA = fopen(radvage, "w");
	
	while(!feof(RAD)){
		fscanf(RAD, "%lf\n", &R_star);
		fscanf(T, "%lf\n", &age_star);
		fscanf(PHASE, "%f\n", &phase_star);
		R_star = pow(10, R_star);
		if(phase_star > 1.0 && phase_star < 7.0){
			fprintf(RVA, "%lf %lf\n", age_star, R_star);
			count_R_star ++;
		}
	}
	fclose(RVA); 
	fclose(engulfed);
	plot_fun(j, Name_of_Run, M_star, a0, increment_a, Pathto_stellar_mods, count_R_star);	
}
if(strcmp(plotting, "yes") == 0){
printf("\n*\n*\n*\n*\nThe resulting orbital evolution will now be plotted. The plot can be viewed in the results_and_plots directory\n*\n*\n*\n*\n*\n");
system("python3 ./plot_scripts/plot_orbital_evolution.py");	
}
}

double calculate_average_lifetime(double mean_mass_of_star){
	char text[1000], t1[100], t2[100], t3[100], temp[1000];
    char pop_mass[100], exo_pop[100], Z[100], age[100];
	float mean_Z, pop_age;
	double a0, a1, a2;
	double mean_lifetime_of_star;



FILE *INPUT = fopen("./configuration", "r"); 
		while(!feof(INPUT))
		{	fscanf(INPUT,"%[^\n]%*c\n", &text);
				if(text[0] == '#') {
					strcpy(text, "This is comment");}		
			strcpy(temp,text);
			sscanf(temp, "%s %s %s", &t1, &t2, &t3);							
					
					
					if (strcmp(t1,"mean_metallicity") == 0)
						{strcpy(Z,t3);
						mean_Z = strtof(Z, NULL);
						}	
							
						}
												
a0 = 10.13 + (0.07547*log10(mean_Z)) - (0.008084*log10(mean_Z)*log10(mean_Z));
a1 = -4.424 - (0.7939*log10(mean_Z)) - (0.1187*log10(mean_Z)*log10(mean_Z));
a2 = 1.262 + (0.3385*log10(mean_Z)) + (0.05417*log10(mean_Z)*log(mean_Z)); 


mean_lifetime_of_star = pow(10,((a0) + (a1*log10(mean_mass_of_star)) + (a2*(log10(mean_mass_of_star)*log10(mean_mass_of_star)))));

return mean_lifetime_of_star;
}

double calculate_single_population(double turnoff_mass){

IMF imf;
  int i;
  double Total_luminous_mass;
  double N_MW, M_MW;
  double N_B_MW, N_S_MW, N_SYSTEMS_MW;
  double Nt, Mt;
  double N_B, N_S, N_SYSTEMS;
  double Mi, Mi_1, Ni;
  double Ni_B, Ni_S, Ni_SYSTEMS;
  double Ni_B_AGB, Ni_B_RGB, Ni_B_SS;
  double Ni_PL_AGB, Ni_PL_RGB, Ni_PL_SS;
  double Npn_total, N_PL_AGB_total, N_B_AGB_total, N_PL_RGB_total, N_B_RGB_total, N_SS_total, Npn_total_PLRGB, Ni_PL_SSinclRGB;
  double F_planet, F_binary, F_single;
  double F_RGB_planet, F_RGB_binary, F_RGB_single;
	
	
	
  double AGB_engulf_probability, RGB_engulf_probability;
  char text[1000], t1[100], t2[100], t3[100], temp[1000];
  char pop_mass[100], exo_pop[100];
  double PN_lifetime, mean_mass_of_star, mean_lifetime_of_star;
  double N_PN_TODAY, N_PN_TODAY_PLRGB;
PN_lifetime = 20000;  
  

  FILE *INPUT = fopen("./configuration", "r"); 
		while(!feof(INPUT))
		{	fscanf(INPUT,"%[^\n]%*c\n", &text);
				if(text[0] == '#') {
					strcpy(text, "This is comment");}		
			strcpy(temp,text);
			sscanf(temp, "%s %s %s", &t1, &t2, &t3);							
					if (strcmp(t1,"population_mass") == 0)
						{strcpy(pop_mass,t3);
						}	
					if (strcmp(t1,"planet_sample") == 0)
						{strcpy(exo_pop,t3);
						}	
						
						}
imf_init_user(&imf,"0.08(pow:-1.3)0.5(pow:-2.3)120");						
Total_luminous_mass = strtof(pop_mass, NULL);	//Read in total mass				
imf_norm_cl(&imf,Total_luminous_mass,120); //total mass locked up in stars						
N_MW = imf_int_xi_cl(&imf,0.08,120.); //total number of STARS including binary systems 
M_MW = imf_int_mxi_cl(&imf,0.08,120.);//total mass locked up in stars in milky way
N_B_MW = 0.67*N_MW;//Total number of binaries in Milkyway
N_S_MW = 0.33*N_MW;
N_SYSTEMS_MW = (N_B_MW/2) + N_S_MW;
//printf("The total mass in Milky Way: %g\nThe total number of stars is then %g\n", M_MW, N_MW);
//printf("Average mass of star: %g\n", M_MW/N_MW);
//printf("Total number of SYSTEMS in Milky Way: %g\n", N_SYSTEMS_MW);
//printf("Total number of Binary systems in Milky Way: %g\n", (N_B_MW/2));
//printf("Total number of Binary systems in Milky Way: %g\n", N_S_MW);

Mt = imf_int_mxi_cl(&imf,turnoff_mass,8.0);
Nt = imf_int_xi_cl(&imf,turnoff_mass,8.0); //Total # PN progenitor STARS
//printf("The fraction of all stars in PN progenitor range: %g\n", Nt/N_MW);
//printf("Total number of PN progenitor STARS: %g\n", Nt);
//printf("Total stellar mass in PN progenitors: %g\n", Mt);
mean_mass_of_star = Mt/Nt;
//printf("Average mass of star in PN progenitor range: %g\n", Mt/Nt);
mean_lifetime_of_star = calculate_average_lifetime(mean_mass_of_star);
//printf("MEAN_LIFETIMEOF STAR = %g\n", mean_lifetime_of_star);
N_B = 0.67*Nt; //number of binaries from 50% progenitor fraction
N_S = 0.33*Nt;
N_SYSTEMS = (N_B/2) + N_S; //Total number of PN progenitor SYSTEMS 


FILE *EXOPLANET_FRAC_AGB = fopen("./results_and_plots/exoplanets_org_EXO_ORGpop_results/probability_having_and_engulfing_giantAGB.csv", "r"); 
FILE *EXOPLANET_FRAC_RGB = fopen("./results_and_plots/exoplanets_org_EXO_ORGpop_results/probability_having_and_engulfing_giantRGB.csv", "r");


FILE *IDA_FRAC_AGB = fopen("./results_and_plots/IDA_2017_all_mass_results/probability_having_and_engulfing_giantAGB.csv", "r"); 
FILE *IDA_FRAC_RGB = fopen("./results_and_plots/IDA_2017_all_mass_results/probability_having_and_engulfing_giantRGB.csv", "r");
FILE *MASSBINS = fopen("./results_and_plots/exoplanets_org_EXO_ORGpop_results/massbin_results.csv", "w");
fprintf(MASSBINS,"Mi-Mi+1 Ni_SYSTEMS probengulfRGB probengulfAGB Ni_PL_RGB Ni_PL_AGB Ni_B_RGB Ni_B_AGB Ni_SS\n");
for(i=(turnoff_mass*10); i<80; i++){
if(strcmp(exo_pop, "exoplanets_org") == 0){

fscanf(EXOPLANET_FRAC_AGB,"%lf\n", &AGB_engulf_probability);
fscanf(EXOPLANET_FRAC_RGB,"%lf\n", &RGB_engulf_probability);

}

if(strcmp(exo_pop, "IDA_2017") == 0){

fscanf(IDA_FRAC_AGB,"%lf\n", &AGB_engulf_probability);
fscanf(IDA_FRAC_RGB,"%lf\n", &RGB_engulf_probability);
printf("%lf\n", AGB_engulf_probability);

}

Mi = (float)i/10; 
Mi_1 = ((float)i+1)/10;
Ni = imf_int_xi_cl(&imf,Mi,Mi_1);
Ni_B = 0.67*Ni;
Ni_S = 0.33*Ni;
Ni_SYSTEMS = (Ni_B/2) + Ni_S;

Ni_B_AGB = 0.32*(Ni_B/2);
Ni_B_RGB = 0.28*(Ni_B/2);
Ni_B_SS = 0.4*(Ni_B/2);

Ni_PL_AGB = AGB_engulf_probability*Ni_S;
Ni_PL_RGB = RGB_engulf_probability*Ni_S;
Ni_PL_SS = Ni_S - (Ni_PL_AGB + Ni_PL_RGB);
Ni_PL_SSinclRGB = Ni_S - Ni_PL_AGB;

Npn_total = Npn_total+(Ni_B_AGB + Ni_PL_AGB + Ni_B_SS + Ni_PL_SS);
Npn_total_PLRGB = Npn_total_PLRGB+(Ni_B_AGB + Ni_PL_AGB + Ni_B_SS + Ni_PL_SSinclRGB);
N_PL_AGB_total = N_PL_AGB_total + Ni_PL_AGB;
N_B_AGB_total = N_B_AGB_total + Ni_B_AGB;
N_SS_total = N_SS_total + (Ni_PL_SS + Ni_B_SS);

N_PL_RGB_total = N_PL_RGB_total + Ni_PL_RGB;
N_B_RGB_total = N_B_RGB_total + Ni_B_RGB;

fprintf(MASSBINS,"%.01f-%.01f &%.01e &%.03f & %.03f &%.01e &%.01e &%.01e &%.01e &%.01e\\\\ \n",Mi, Mi_1,Ni_SYSTEMS, RGB_engulf_probability, AGB_engulf_probability, Ni_PL_RGB, Ni_PL_AGB, Ni_B_RGB, Ni_B_AGB, (Ni_PL_SS+Ni_B_SS) );
}


printf("total number of Planet rgb interactions: %g\n", N_PL_RGB_total);
printf("total number of Planet agb interactions: %g\n", N_PL_AGB_total);
printf("total number of binary rgb %g\n", N_B_RGB_total);
printf("total number of binary agb %g\n", N_B_AGB_total);
printf("total number of stars evolving as singles %g\n", N_SS_total);
printf("total number of binary rgb %g\n", N_B_RGB_total);
printf("total number of PNe if planet RGB interactions prevent PN: %g\n", Npn_total);
printf("total number of PNe if planet RGB interactions evolve as singles: %g\n", Npn_total_PLRGB);
printf("total number of nonspherical PNe: %g\n\n\n\n", N_B_AGB_total+N_PL_AGB_total);


F_planet = N_PL_AGB_total/Npn_total;
F_binary = N_B_AGB_total/Npn_total;
F_single = N_SS_total/Npn_total;

printf("planet fraction %g\n", F_planet);
printf("binary fraction %g\n", F_binary);
printf("single stars fraction %g\n\n\n\n", F_single);


F_RGB_planet = N_PL_AGB_total/Npn_total_PLRGB;
F_RGB_binary = N_B_AGB_total/Npn_total_PLRGB;
F_RGB_single = (N_SS_total+N_PL_RGB_total)/Npn_total_PLRGB;

printf("planet fraction (without planet rgb) %g\n", F_RGB_planet);
printf("binary fraction %g\n", F_RGB_binary);
printf("single stars fraction %g\n", F_RGB_single);


//printf("The planet fraction, with planet RGB interactions, is: %g\nThe binary fraction is: %g\nThe single/spherical fraction is: %g\n\n\n\n\n", F_planet, F_binary, F_single);

//printf("total number of PNe NOT taking RGB planet interactions into account: %g\n", Npn_total_PLRGB);
//printf("The planet fraction, without planet RGB interactions, is: %g\nThe binary fraction is: %g\nThe single/spherical fraction is: %g\n", F_RGB_planet, F_RGB_binary, F_RGB_single);

N_PN_TODAY = Npn_total*(PN_lifetime/mean_lifetime_of_star);
N_PN_TODAY_PLRGB = Npn_total_PLRGB*(PN_lifetime/mean_lifetime_of_star);

//printf("\n\n\nThe total_number of PNe still visible today: %g\n\n", N_PN_TODAY);
//printf("\n\n\nThe total_number of PNe still visible today (minus RGB planets): %g\n\n", N_PN_TODAY_PLRGB);
fclose(EXOPLANET_FRAC_AGB);
fclose(EXOPLANET_FRAC_RGB);
return 0;
}

int main() {
char text[1000], t1[100], t2[100], t3[100], temp[1000], back_of_envelope[10], age[20], Z[10], user_answer[1000];
char MODE[20], planet_pop[50], Name_of_Run[1000], plotting_string[1000], Results_directory_string[1000], stellar_model[1000];
double pop_age, turnoff_mass, mean_Z;
welcomeNote();
FILE *INPUT = fopen("./configuration", "r"); 
	while(!feof(INPUT)){
		fscanf(INPUT,"%[^\n]%*c\n", &text);
		if(text[0] == '#') {
			strcpy(text, "This is comment");}		
			strcpy(temp,text);
			sscanf(temp, "%s %s %s", &t1, &t2, &t3);							
		if (strcmp(t1,"Mode") == 0)
		{	strcpy(MODE,t3);
		}	
		if (strcmp(t1,"name_your_experiment") == 0)
		{	strcpy(Name_of_Run,t3);						
		}		
		if (strcmp(t1,"planet_sample") == 0)
		{	strcpy(planet_pop,t3);
		}
		if (strcmp(t1,"back_of_envelope") == 0)
		{	strcpy(back_of_envelope,t3);
		}
		if (strcmp(t1,"population_age") == 0)
		{	strcpy(age,t3);
			pop_age = strtof(age, NULL);
		}
		if (strcmp(t1,"mean_metallicity") == 0)
		{	strcpy(Z,t3);
			mean_Z = strtof(Z, NULL);
		}
		if (strcmp(t1,"star_models") == 0)
		{	strcpy(stellar_model,t3);			
		}		
	}
sprintf(Results_directory_string, "./results_and_plots/%s_Results", Name_of_Run);
mkdir(Results_directory_string, 0777);
sprintf(plotting_string, "./results_and_plots/%s_Plots", Name_of_Run);
mkdir(plotting_string, 0777);

if (strcmp(MODE,"Orbital")==0){
	
	if (strcmp(stellar_model,"CREATE_NEW_MESA")==0){		
		printf("!******You have chosen to create a new stellar model with MESA... This may take a while.******\n");
		printf("!****Go grab a coffee, take the dog for a walk, maybe have a bag of cans (not recommended).***\n");
		create_new_mesa();
	}

	Orbital();
	system("rm system_temp.txt");
	system("rm engulfed.txt");
}

if (strcmp(MODE,"Single_epoch_pop")==0){
	if (strcmp(planet_pop,"exoplanets_org")==0 || strcmp(planet_pop,"IDA_2017")==0){	
		if (strcmp(planet_pop,"IDA_2017")==0){
			printf("You have chosen to evolve and analyse synthetic planet populations (of Shigeru Ida) in Single epoch population mode.\n\nThis may take a long time due to large number of systems (1 day).\n\nDo you want to proceed? (yes to proceed, no to use existing evolved data from previous simulations)\n");
			scanf("%s", &user_answer);
			if (strcmp(user_answer,"yes")==0){			
				Single_epoch("08M");//0.8m
				Multi_planets("08M");
				Single_epoch("1M");//1m
				Multi_planets("1M");
				Single_epoch("15M");//1.5
				Multi_planets("15M");
				Single_epoch("2M");//2.0
				Multi_planets("2M");
				Single_epoch("3M");//3.0
				Multi_planets("3M");
				//Single_epoch("5M");//5.0
				//Multi_planets("5Mtest");
				//Single_epoch("8M");//8.0
				//Multi_planets("8M");	
			}
			system("python3 ./IDA_data_analysis_scripts/IDA_after_engulf.py");	
			system("python3 ./IDA_data_analysis_scripts/AGBLogReg.py");
			system("python3 ./IDA_data_analysis_scripts/RGBLogReg.py");
			calc_exo_probs();
			}	
			
		if (strcmp(planet_pop,"exoplanets_org")==0){
			printf("You have chosen to evolve and analyse data from exoplanets.org in Single epoch population mode.\n\nThis will download and evolve a fresh copy of exoplanets.csv.\n\nDo you want to proceed? (yes to proceed, no to use existing evolved data from previous simulations)\n");
			scanf("%s", &user_answer);
			if (strcmp(user_answer,"yes")==0){			
				Single_epoch("EXO_ORG");	
				Multi_planets("EXO_ORG");
			}
			system("python3 ./exo_org_analysis_scripts/exo_after_engulf.py");
			system("python3 ./exo_org_analysis_scripts/AGBLogReg.py");
			system("python3 ./exo_org_analysis_scripts/RGBLogReg.py");
			calc_exo_probs();

			}
			
					
			
}
else{
			Single_epoch("SPECIFIED");


		}
// use turnoff mass is yes
if (strcmp(back_of_envelope,"yes") == 0){
	turnoff_mass = calculate_turnoff_mass(pop_age, mean_Z);
	calculate_single_population(turnoff_mass);
}
//if use turnoff mass is no
if (strcmp(back_of_envelope,"no") == 0){
	calculate_single_population(0.8);
}
system("rm system_temp.txt");
system("rm engulfed.txt");
}

if (strcmp(MODE,"Multi_epoch_pop")==0){
	Multi_epoch();	

}
return 0;
}