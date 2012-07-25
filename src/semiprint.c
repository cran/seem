        #include <stdio.h>
        #include <string.h>
        #define semiprint  semiprint_

	float x[20],xf[20],p[20][20];
	
	void semiprint(int *mst, float *dt, float *tfin, float *xf,
		       float (*p)[20], float *t, float *x){
        FILE *ofp;
	int i,j;
	float totcheck;

	if(*t >0.0) ofp=fopen("semi_out.txt","a");
	else {
	 ofp=fopen("semi_out.txt","w");
         /* print heading*/
	 fprintf(ofp,"Number of states  %d\n",*mst);
	 fprintf(ofp,"Delta t  %6.3f\n",*dt);
	 fprintf(ofp,"Simulation time (decades)  %6.3f\n",*tfin);
 	 fprintf(ofp,"Initial condition  \n");
         for(i=0;i<=*mst-1;i++){fprintf(ofp, "%6.3f", xf[i]);}
         fprintf(ofp,"\n");
 	 fprintf(ofp,"Transition probabilities  \n");
	 for(j=0;j<=*mst-1;j++){
          for(i=0;i<=*mst-1;i++){fprintf(ofp, "%6.3f",p[j][i]);}
	  fprintf(ofp,"\n");
	 }
	 fprintf(ofp,"Time then states from 1 to %d and total\n",*mst);
	}
	totcheck=0.0;
        for(i=0;i<=*mst-1;i++){totcheck=totcheck+x[i];}
 	fprintf(ofp, "%6.3f",*t);	
        for(i=0;i<=*mst-1;i++){fprintf(ofp, "%6.3f", x[i]);}
	fprintf(ofp,"%6.3f\n",totcheck);
        fclose(ofp);
        } 
             
        
