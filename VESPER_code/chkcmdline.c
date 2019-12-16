#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include "struct.h"


int chkcmdline( int argc, char **argv,CMD *cmd){
        char **p;
	int num=0;
	void errmsg();
        if(argc < 2){
		errmsg();
                return(FALSE);
        }
	//default values
        p=argv;
	cmd->map_t=0.00;
	cmd->Nthr=2;
	cmd->dreso=16.00;
	cmd->MergeDist=0.50;
	cmd->Filter=0.10;
	cmd->Mode=0;
	cmd->LocalR=10.0;
	cmd->Dkeep=0.5;
	cmd->Nround=5000;
	cmd->Nnb=30;
	cmd->Ntabu=100;
	cmd->Nsim=10;
	cmd->Allow=1.01;
	cmd->Nbeam=20;
	cmd->ssize=7.0;
	cmd->ang=30.0;
	cmd->TopN=10;
	cmd->ShowGrid=false;

	cmd->th1=0.00;
	cmd->th2=0.00;

	cmd->Mode=1;

        while (--argc){
         p++;
         if(**p == '-'){
          	switch(*(*p+1)){
		 case 'i':
			strcpy(cmd->filename,*(++p));
                	--argc; break;
		 case 'a':
			strcpy(cmd->file1,*(++p));
                	--argc; break;
		 case 'b':
			strcpy(cmd->file2,*(++p));
                	--argc; break;
		 case 'c':
			cmd->Nthr=atoi(*(++p)); 
			--argc; break;
		 case 't':
			cmd->th1=atof(*(++p)); 
			--argc; break;
		 case 'T':
			cmd->th2=atof(*(++p)); 
			--argc; break;
		 case 'g':
			cmd->dreso=atof(*(++p)); 
			--argc; break;
		 case 'A':
                        cmd->ang=atof(*(++p));
                        --argc; break;
		 case 'N':
                        cmd->TopN=atoi(*(++p));
                        --argc; break;
		 case 'R':
                        cmd->LocalR=atof(*(++p));
                        --argc; break;
		 case 'k':
                        cmd->Dkeep=atof(*(++p));
                        --argc; break;
		 case 'r':
                        cmd->Nround=atoi(*(++p));
                        --argc; break;
		 //case 'b':
                 //       cmd->Nnb=atoi(*(++p));
                 //       --argc; break;
		 case 'l':
                        cmd->Ntabu=atoi(*(++p));
                        --argc; break;
		 case 's':
                        cmd->ssize=atof(*(++p));
                        --argc; break;
		 //case 'a':
                 //       cmd->Allow=atoi(*(++p));
                 //       --argc; break;
		 //case 'T':
                 //       cmd->Mode=0;
                 //       break;
		 case 'S':
                        cmd->ShowGrid=true;
                        break;
		 case 'V':
                        cmd->Mode=1;
                        break;
		 case 'L':
                        cmd->Mode=2;
                        break;
		 case 'C':
                        cmd->Mode=3;
                        break;
		 case 'P':
                        cmd->Mode=4;
                        break;
		 default: 
		  	fprintf(stderr,"No such a option: %s\n",*p+1); 
			errmsg(); 
			return(FALSE); 
		 break;
	  }
	 }
        }
	//option check
	printf("#Pairwise matching result from VESPER\n");
	printf("#If this VESPER calculation result is used for publication, please cite:\n");
	printf("#Xusi Han, Genki Terashi, Siyang Chen, & Daisuke Kihara. VESPER: Global and Local Cryo-EM Map Alignment Using Local Density Vectors. Manuscript in preparation. (2019).\n");
	printf("#Number of threads= %d\n",cmd->Nthr);
	printf("#Map1 Threshold= %f\n",cmd->th1);
	printf("#Map2 Threshold= %f\n",cmd->th2);
	printf("#Band Width= %f\n",cmd->dreso);
        return(TRUE);
}

void errmsg(){
	puts("Usage: EMVEC_FIT -a [MAP1.mrc (large)] -b [MAP2.mrc (small)] [(option)]");
	puts("v0.10	Start");
	puts("v0.20	Pure Multi threading Process & refine 5 degree.");
	puts("v0.30	Add Overlap Mode & using single-precision");
	puts("v0.31	Fix Bugs in CC calculation");
	puts("v0.32	Add Pearson CC calculation");
	puts("v0.33	Fix problems in translation vector");

	puts("---Options---");
	printf("-t [float] : Threshold of density map1 def=%.3f\n",0.00);
	printf("-T [float] : Threshold of density map2 def=%.3f\n",0.00);
	printf("-g [float] : Bandwidth of the gaussian filter\n");
        printf("             def=16.0, sigma = 0.5*[float]\n");
	printf("-s [float] : Sampling grid space def=7.0\n");
	printf("-A [float] : Sampling Angle interval def=30.0\n");
	printf("-c [int  ] : Number of cores for threads def=%d\n",2);
	printf("-N [int  ] : Refine Top [int] models def=10\n");
	printf("-S 	   : Show topN models in PDB format def=false\n");
	printf("-V 	   : Vector Products Mode def=true\n");
	printf("-L 	   : Overlap Mode def=false\n");
	printf("-C 	   : Cross Correlation Coefficient Mode def=false\n");
        printf("             Using normalized density value by Gaussian Filter\n");
	printf("-P         : Pearson Correlation Coefficient Mode def=false\n");
        printf("             Using normalized density value by Gaussian Filter and average density\n");
	printf("Thi is Ver %.3f\n",VER);
}
