        #include <stdio.h>
        #include <string.h>
        #define ztracer  ztracer_

	float xba[10];
	
	void ztracer(int *nspp, int *kyr, float *xd, float *xb, 
		    float *sdb, float *xbatot, float *xl, float *xh, float *xba){
        FILE *ofp;
	int i;

	if(*kyr >0)  ofp=fopen("z_tracer_out.txt","a");
	else ofp=fopen("z_tracer_out.txt","w");

 	fprintf(ofp, "%3d %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f",
	               *kyr, *xd, *xb, *sdb, *xbatot, *xl, *xh);	

        for(i=0;i<=*nspp-1;i++){fprintf(ofp, "%6.2f", xba[i]);}
	fprintf(ofp,"\n");

         fclose(ofp);
        } 
             
        
