        #include <stdio.h>
        #include <string.h>
        #define zprint  zprint_

	float sd[10], dc[10],rsd[10],ba[10],rba[10],riv[10],freq[10];
	float sdc[10][10];
	char cp1[10];
	
	void zprint(int *nspp, char *cname, int *kyr, float *xarea, float *sd, 
		    float (*sdc)[10], float *dc, float *rsd, float *ba, float *rba,
		    float *riv, float *freq, float *tdha, float *tdb, float*tbaha,
		    float *xdbh, float *sdd, float *tbmha, float *xlai, float *xht){
        FILE *ofp;
        char *cp, *cp1[10];
        size_t len, len4;
	int i,j;

	cp = (cname);    /* convert to C char* */
        len = strlen(cname);

        /* strip trailing blanks */
        while (cp[len-1] == ' ' || cp[len-1] == '\0') --len;
	len4 = len/(*nspp);
        /* parse species names */
 	 cp1[0] = strtok(cp,",");	
	 for(j=1;j<=*nspp-1;j++){cp1[j] =strtok(NULL,",");}	


	if(*kyr >0)  ofp=fopen("z_print_out.txt","a");
	else ofp=fopen("z_print_out.txt","w");
        /*     Write stand-level output ... */	 
         fprintf(ofp, "Simulation year: %d \n",*kyr);
 	 fprintf(ofp, "Stand Structure by Species:\n");
	 fprintf(ofp, "Species Dbh Distribution (#/ha, in 10-cm classes):\n");

	  fprintf(ofp, "%s", cp1[0]);
          for(i=0;i<=9;i++){fprintf(ofp, "%10.2f",sdc[0][i]/(*xarea));}
	  fprintf(ofp,"\n");

	 for(j=1;j<=*nspp-1;j++){
          fprintf(ofp, "%s", cp1[j]);
          for(i=0;i<=9;i++){fprintf(ofp, "%10.2f",sdc[j][i]/(*xarea));}
	  fprintf(ofp,"\n");
	 }

	 fprintf(ofp, "All ");
         for(i=0;i<=9;i++){fprintf(ofp, "%10.2f",dc[i]/(*xarea));}
	 fprintf(ofp,"\n");
	 
	 fprintf(ofp,"\n");

	 fprintf(ofp, "Spp.  Dens. Rel.Dens. BA Rel.BA IV200 Freq.:\n");
	  fprintf(ofp, "%s", cp1[0]);
 	 fprintf(ofp, "%6.1f %6.1f %6.1f %6.2f %6.2f %6.2f\n",
	           sd[0]/(*xarea), rsd[0], ba[0]/(*xarea), rba[0], riv[0], freq[0]);	

         for(i=1;i<=*nspp-1;i++){
       		 fprintf(ofp, "%s", cp1[i]);
	    fprintf(ofp, "%6.1f %6.1f %6.1f %6.2f %6.2f %6.2f\n",
	           sd[i]/(*xarea), rsd[i], ba[i]/(*xarea), rba[i], riv[i], freq[i]);	
	 }
	  fprintf(ofp,"\n");

 	  fprintf(ofp,"Stand Aggregates:\n");

	fprintf(ofp,"Density 1/ha: (total) %6.2f and of stems>10cm %6.2f\n", *tdha,*tdb);
	fprintf(ofp,"Basal Area (sq.m/ha) %6.2f\n",*tbaha);
	fprintf(ofp,"Mean Dbh (cm), %6.2f with s.d. %6.2f\n", *xdbh,*sdd);
	fprintf(ofp,"Woody biomass (Mg/ha) %6.2f\n",*tbmha);
	fprintf(ofp,"Leaf-area index %6.2f\n",*xlai);
	fprintf(ofp,"Average canopy height (m) %6.2f\n",*xht);
        fprintf(ofp,"\n");

         fclose(ofp);
        } 
             
        
