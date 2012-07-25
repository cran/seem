        #include <stdio.h>
        #include <string.h>
        #define cerioprint  cerioprint_
        
 	float yd[10];

	void cerioprint(int *irun, char *cname, float *flag, float *t, float *time0,
			int *nruns, int *nstate, int *ndata, char *subtitle, float *yd){
        FILE *ofp;
        char *cp, *cp1, *st[10];
        size_t len, len1;
	int i;

	cp = (cname);    /* convert to C char* */
        len = strlen(cname);

        /* strip trailing blanks */
        while (cp[len-1] == ' ' || cp[len-1] == '\0') --len;

	cp1 = (subtitle);    /* convert to C char* */
        len1 = strlen(subtitle);

        /* strip trailing blanks */
        while (cp1[len1-1] == ' ' || cp1[len1-1] == '\0') --len1;


	st[0]= "Eggs      ";
      	st[1]= "Fem_neon  ";
      	st[2]= "Fem_adult ";
      	st[3]= "Males     ";
      	st[4]= "Algae1    ";
      	st[5]= "Algae2    ";
      	st[6]= "Total     ";
      	st[7]= "Temp      ";
      	st[8]= "Photoper  ";
      	st[9]= "Light     ";


	if(*irun > 1)  {
	 if(*t == *time0){
	 ofp=fopen("cerio_out.txt","a");
	 fprintf(ofp,"Run number=%3d:  %s =%10.4e\n",*irun, cp, *flag);  
         fprintf(ofp,"%10.2f ",*t);
	 for(i=0;i<=*nstate-1;i++){fprintf(ofp, "%10.3e ",yd[i]);}
	 fprintf(ofp,"\n");
	 }
	 else{
	 ofp=fopen("cerio_out.txt","a");
         fprintf(ofp,"%10.2f ",*t);
	 for(i=0;i<=*nstate-1;i++){fprintf(ofp, "%10.3e ",yd[i]);}
	 fprintf(ofp,"\n");
	 }
	} 
	else {
	 if(*t == *time0){	
	 ofp=fopen("cerio_out.txt","w");
         fprintf(ofp, "Number_of_runs: %d \n",*nruns);
 	 fprintf(ofp, "Number_of_states: %d \n",*nstate);
	 fprintf(ofp, "Number_of_data: %d \n", *ndata);
	 fprintf(ofp, "Title: Cladocera_Stage_Structure\n");
	 fprintf(ofp, "Sub_title: %s\n", cp1);
	 fprintf(ofp, "Y_axis_title: Dens_or_Temp_or_Phper_or_Light\n");
	 fprintf(ofp, "Time__[days] ");	 
         for(i=0;i<=*nstate-1;i++){fprintf(ofp, "%s",st[i]);}
	 fprintf(ofp,"\n");
	 fprintf(ofp,"Run number=%3d:  %s =%10.4e\n",*irun, cp, *flag);
         fprintf(ofp,"%10.2f ",*t);
	 for(i=0;i<=*nstate-1;i++){fprintf(ofp, "%10.3e ",yd[i]);}
	 fprintf(ofp,"\n");
	 }
	 else{
	 ofp=fopen("cerio_out.txt","a");
         fprintf(ofp,"%10.2f ",*t);
	 for(i=0;i<=*nstate-1;i++){fprintf(ofp, "%10.3e ",yd[i]);}
	 fprintf(ofp,"\n");
	 }
       }

         fclose(ofp);
        } 
             
        
