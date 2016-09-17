// g++ -g -o call_giac call_giac.C -lgiac -lgmp

#include <stdio.h>
#include <unistd.h> /* getopt */

#include "gmp.h"
#include "giac.h"

using namespace std;
using namespace giac;


int main(int argc,char**argv)
{
    int nb_vars;
    int n_input;
    int n_output;
    int nb_mons;
    int p, m, v;
    mpz_t u;
    int greduce = 0;

    int c;
    opterr = 0;
    while ((c = getopt (argc, argv, "r:")) != -1)
        switch (c)
        {
        case 'r':
            sscanf(optarg, "%d", &greduce);
            break;
        default:
        ;
    }

    
    scanf("%d", &nb_vars);
    scanf("%d", &n_input);

    index_t e(nb_vars);
    vectpoly input_basis; 
    mpz_init(u);

    for (p=0; p<n_input; p++) {
        scanf("%d", &nb_mons);
        polynome prev(nb_vars);
        for (m=0; m<nb_mons; m++) {
            for (v=0; v<nb_vars; v++) {
                scanf("%hd", &e[v]);
            }    
            mpz_inp_str(u, stdin, 10);
            prev.coord.push_back(monomial<gen>(gen(u) ,e));
        }
        input_basis.push_back(prev);
    }

#if 0    
    printf("-- INPUT ----------\n");
    nb_vars = input_basis.front().dim;
    n_input = input_basis.size();
    printf("%d %d\n", nb_vars, n_input);
    for (p=0; p<n_input; p++) {
        nb_mons=input_basis[p].coord.size();
        printf("%d\n", nb_mons);
        for (m=0; m<nb_mons; m++) {
            index_t::const_iterator it=input_basis[p].coord[m].index.begin();
            index_t::const_iterator itend=input_basis[p].coord[m].index.end();
            for (;it!=itend;++it) {
                 printf("%hd ", *it);
            }     
            mpz_out_str(stdout,10,
                *(input_basis[p].coord[m].value.ref_ZINTptr()));
            printf("\n");
        }    
    }
#endif     

    vecteur variables;
    for (v=0; v<nb_vars; v++) {
#if 0    
        //std::string name = "x" + std::to_string(v);
        std::stringstream name;
        name << "x" << v;
        variables.push_back(gen(identificateur(name.str())));
#else
        variables.push_back(gen(identificateur()));
#endif
    }


    vecteur input_sym;

    vecteur args;
    if (greduce>=1) {
        gen F = r2sym(input_basis[0], variables, 0);
        args.push_back(F);
        p = 1;
    } 
    else {
        p = 0;
    }

    for (; p<n_input; p++) {
        gen s = r2sym(input_basis[p], variables, 0);
        input_sym.push_back(s);
    }    

    args.push_back(gen(input_sym));
    args.push_back(gen(variables));
   
    vecteur *output_vec;
    gen output_sym;
    if (greduce>=1) {
        output_sym = _greduce(args, 0);
        n_output = 1;
    }
    else {
        output_sym = _gbasis(args, 0);
        output_vec = output_sym.ref_VECTptr(); 
        n_output = output_vec->size();
    }    

#if 0
    cout << input_sym << endl;
    cout << variables << endl;
    printf("-- OUTPUT ---------\n");
    cout << output_sym << endl;
#endif
    
    printf("%d %d\n", nb_vars, n_output + (greduce?1:0) );
    for (p=0; p<n_output; p++) {
        fraction ff = greduce?
             sym2r(output_sym, variables, 0) :
             sym2r((*output_vec)[p], variables, 0);
        polynome *pp = ff.num.ref_POLYptr();
        nb_mons=pp->coord.size();
        printf("%d\n", nb_mons);
        for (m=0; m<nb_mons; m++) {
            index_t::const_iterator it=pp->coord[m].index.begin();
            index_t::const_iterator itend=pp->coord[m].index.end();
            for (;it!=itend;++it) {
                 printf("%hd ", *it);
            }    
            switch(pp->coord[m].value.type) {
            case _INT_:                    
                 printf("%i", pp->coord[m].value.val);
                break;
            case _ZINT:                     
                mpz_out_str(stdout,10,
                    *(pp->coord[m].value.ref_ZINTptr()));
                break;
            default:
                printf(" TYPE(%i)", int(pp->coord[m].value.type));
            }
            printf("\n");
        }    
    }
    if (greduce) {  // print denominator
        printf("1\n");
        for (v=0; v<nb_vars; v++) {
            printf("0 ");
        }    
        fraction ff = sym2r(output_sym, variables, 0);
        switch(ff.den.type) {
            case _INT_:                    
                 printf("%i", ff.den.val);
                break;
            case _ZINT:                     
                mpz_out_str(stdout,10,
                    *(ff.den.ref_ZINTptr()));
                break;
            default:
                printf(" TYPE(%i)", int(ff.den.type));
        }
        printf("\n");
    }

}
