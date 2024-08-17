#ifndef EVOLVEI_H
#define EVOLVEI_H

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <string>
#include <gsl/gsl_rng.h>
#include <set>
#include <list>
#include "basics.h"
#include "graphi.h"
#include "fitnessi.h"

using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::ostream;
using std::ifstream;
using std::istream;
using std::string;
using std::string;
using std::set;
using std::list;
using std::ifstream;
using std::ios;

class EvolveI
{
public:
	EvolveI();
  EvolveI(Alea& jacta);
  void start_rng(Alea& jacta);
  void set_params(double muratep, double wolp, int popsizep, int nugenp);
  void freeze(EvolveI &anc, Alea& jacta); //no copia atractores, pero si w
  void print_params(ostream& sal);
  void start_pop(GraphI &founder);
  void start_pop(int n, int e);
  void start_pop(int n, int e, int *ic, int *goal);
  void start_pop(int n, int e, int *ic, int **goal, int per);
  void start_pop(int n, int e, int **maic, int nuic, int **goal);
  void start_pop(int n, int e, int **maic, int nuic, int ***goals, int *vsigo);
  void one_generation();
  void one_generation_mp5();
  void one_generation_mf();
  void one_generation_sex();
  void assign_w(int cual, double cuanto);
  void calc_meanw();
  double return_meanw();
  double return_maxw();
  void optima(set<int> &opt);
  void optima_strict(set<int> &opt);
  int num_optima();
  int num_optima_strict(); 
  void clear();
  void clear_lop();
  double get_murate();
  double get_wol();
  void new_wol(double nwol);
  int get_popsize();
  int get_nugen();
  double get_w(int cual);

  int group_by_adjmat(int *vecid); //new outside
  
  int number_of_genotypes();
  int number_of_genotypes_thr();
  int number_of_genotypes(set<int> &opt);
  int number_of_genotypes_thr(set<int> &opt);
  void genotype_statistics(); //vector must be of size equal to (popsize*(popsize-1))/2
  void genotype_statistics(set<int> &opt, double& megd, double& megdth, double& magd, double& magdth, double *agd, double *agdth);
  int number_of_phenotypes(); //for single-attractor simulations
  int number_of_phenotypes(set<int> &opt, double& mepd, double& mapd, double *apd);
  int number_of_phenotypes(int *nuper); //for single-attractor simulations. size of array =5. Index denotes period. nuper[4] is for periods >= 4.
  int number_of_phenotypes(set<int> &opt, double& mepd, double& mapd, double *apd, int *nuper);
  void build_last_opt_phen(int vagener, int lastgener, int samrat, set<int> &opt);
  double analyze_last_opt_phen(int *hasdo21, double& resnor);

  
  GraphI *population;
  double mean_gen_dist, mean_gen_dist_thr, max_gen_dist, max_gen_dist_thr, *all_gdists, *all_gdists_thr;
  double mean_ph_dist, max_ph_dist, *all_pdists;
  
  int ****last_opt_phen; //
  int **perlop;
  int *nulop;
  
private:

  Alea est;
  Basics basic;
  double murate;
  double wol;
  int popsize;
  int nugen;
  bool yaparams;
  double meanw;
  double maxw;
  double *w;
  GraphI *nepop;
  bool yapop;
  
  bool yalop;//
  int vanop, generini;

};

#endif
