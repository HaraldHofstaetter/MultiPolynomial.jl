#include <stdio.h>
#include <unistd.h> /* getopt */

#define LIBMODE 2

#include "call_fgb.h"
#include "gmp.h"
#define FGb_MAXI_BASE 100000


int main(int argc,char**argv)
{
    Dpol *input_basis;
    Dpol *output_basis;

    UI32 nb_vars;
    UI32 n_input;
    UI32 n_output;
    UI32 nb_mons;
    UI32 p, m, v;
    UI32 *e;
    mpz_t u;
    Dpol_INT prev;
  
    //char** vars = NULL;
    char* vars[6]={"x1","x2","x3","x4","x5","x6"}; /* name of the variables (can be anything) */
    char buf[1024];

    double t0;

    UI32* Mons;
    mpz_ptr* cfs;    
    int code;

    int n_threads = 1;
    UI32 k2_dlr = 0;
    UI32 n_output_max = FGb_MAXI_BASE;

    SFGB_Options options;
    FGb_set_default_options(&options);
    /* overide some default parameters */
    options._env._force_elim=0; /* if force_elim=1 then return only the result of the elimination 
				    (need to define a monomial ordering DRL(k1,k2) with k2>0 ) */
    options._env._index=1000000; /* This is is the maximal size of the matrices generated by F4 
				     you can increase this value according to your memory */

    {
        int c;
        opterr = 0;
        while ((c = getopt (argc, argv, "t:k:c:o:n:e:m:")) != -1)
            switch (c)
            {
            case 't':
                 sscanf(optarg, "%d", &n_threads);
                 break;
            case 'k':
                 sscanf(optarg, "%u", &k2_dlr);
                 break;
            case 'o':
                 sscanf(optarg, "%u", &n_output_max);
                 break;
            case 'c':
                 {
                     int _compute = 1;
                     sscanf(optarg, "%d", &_compute);
                     options._env._compute = _compute;
                 }
                 break;
            case 'n':
                 sscanf(optarg, "%d", &options._env._nb);
                 break;
            case 'm':
                 sscanf(optarg, "%u", &options._env._index);
                 break;
     /*     case 'v':
                 sscanf(optarg, "%d", &options._verb);
                 break;
      */
            case 'e':
                 sscanf(optarg, "%d", &options._env._force_elim);
                 break;
            default:
            ;
         }
    }
#if 1    
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "-t n_threads %i\n", n_threads);
    fprintf(stderr, "-k k2_dlr %u\n", k2_dlr);
    fprintf(stderr, "-o n_output_max %u\n", n_output_max);
    fprintf(stderr, "-c options._env._compute %i\n", options._env._compute);
    fprintf(stderr, "-n options._env._nb %i\n", options._env._nb);
    fprintf(stderr, "-m options._env._index %u\n", options._env._index);
    fprintf(stderr, "-e options._env._force_elim %i\n", options._env._force_elim);
#endif

    FGB(saveptr)(); /* First thing to do : GMP origmal memory allocators are saved */
    threads_FGb(n_threads);
    init_FGb_Integers(); /* init FGb for integers computation  */

    scanf("%d", &nb_vars);
    FGB(PowerSet)(nb_vars-k2_dlr,k2_dlr,vars); 
                                   /* Define the monomial ordering: DRL(k1,k2) where 
		   		      k1 is the size of the 1st block of variables 
				      k2 is the size of the 2nd block of variables 
				      and vars is the name of the variable
				   */

    scanf("%d", &n_input);

    e = (int*)malloc(sizeof(UI32)*nb_vars);
    input_basis = (Dpol*)malloc(sizeof(Dpol)*n_input);

    mpz_init(u);

    for (p=0; p<n_input; p++) {
        scanf("%d", &nb_mons);
        prev=FGB(creat_poly)(nb_mons);
        input_basis[p]=prev;
        for (m=0; m<nb_mons; m++) {
            for (v=0; v<nb_vars; v++) {
                scanf("%d", &e[v]);
            }    
            FGB(set_expos2)(prev,m,e,nb_vars);
            mpz_inp_str(u, stdin, 10);
            FGB(set_coeff_gmp)(prev,m,u);
        }
        FGB(full_sort_poly2)(prev);/* it is recommended to sort each polynomial */
    }
    free(e);

#if 0
    printf("-- INPUT ----------\n");
    printf("%d %d\n", nb_vars, n_input);
    for (p=0; p<n_input; p++) {
        nb_mons=FGB(nb_terms)(input_basis[p]);
        printf("%d\n", nb_mons);
	Mons=(UI32*)(malloc(sizeof(UI32)*nb_vars*nb_mons));
	cfs=(mpz_ptr*)(malloc(nb_mons*sizeof(mpz_ptr)));
	code=FGB(export_poly_INT_gmp2)(nb_vars,nb_mons,cfs,Mons,input_basis[p]);
        for (m=0; m<nb_mons; m++) {
            e=Mons+m*nb_vars;
            for (v=0; v<nb_vars; v++) {
                printf("%u ", e[v]);
            }    
            mpz_out_str(stdout,10,cfs[m]);
            printf("\n");
        }    
	free(Mons);
	free(cfs);
    }
    printf("-- OUTPUT ----------\n");
#endif    

#if 1
    output_basis = (Dpol*)malloc(sizeof(Dpol)*n_output_max);

    n_output=FGB(fgb)(input_basis,n_input,output_basis,n_output_max,&t0,&options);

    printf("%d %d\n", nb_vars, n_output);
    for (p=0; p<n_output; p++) {
        nb_mons=FGB(nb_terms)(output_basis[p]);
        printf("%d\n", nb_mons);
	Mons=(UI32*)(malloc(sizeof(UI32)*nb_vars*nb_mons));
	cfs=(mpz_ptr*)(malloc(nb_mons*sizeof(mpz_ptr)));
	code=FGB(export_poly_INT_gmp2)(nb_vars,nb_mons,cfs,Mons,output_basis[p]);
        for (m=0; m<nb_mons; m++) {
            e=Mons+m*nb_vars;
            for (v=0; v<nb_vars; v++) {
                printf("%u ", e[v]);
            }    
            mpz_out_str(stdout,10,cfs[m]);
            printf("\n");
        }    
	free(Mons);
	free(cfs);
    }

    free(output_basis);
#endif

    free(input_basis);

    FGB(reset_memory)(); /* to reset Memory */
    FGB(restoreptr)(); /* restore original GMP allocators */
}
