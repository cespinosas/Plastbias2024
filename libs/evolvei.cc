#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <set>
#include <list>
#include <string>
#include "basics.h"
#include "graphi.h"
#include "fitnessi.h"
#include "evolvei.h"

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

EvolveI::EvolveI()
{
  yaparams = false;
  yapop = false;
  yalop = false;
}

EvolveI::EvolveI(Alea& jacta) //creates instance of class EvolveI and assigns rng to jacta
{
  start_rng(jacta);
  yaparams = false;
  yapop = false;
  yalop = false;
}

void EvolveI::start_rng(Alea& jacta) //assigns rng to jacta 
{
  est = jacta;
  basic.start_rng(est);
}

void EvolveI::set_params(double muratep, double wolp, int popsizep, int nugenp) //sets parameters for evolutionary simulations
{
  murate = muratep;
  wol = wolp;
  popsize = popsizep;
  nugen = nugenp;
  yaparams = true;
}

void EvolveI::freeze(EvolveI &anc, Alea& jacta)
{
  if (!yaparams) {
    cout << "[Error]: Parameters had not been set when EvolveI::freeze was called.\n";
    exit(1);
  }
  if (!yapop) {
    cout << "[Error]: Population not already constructed when EvolveI::freeze was called.\n";
    exit(1);
  }
  start_rng(jacta);
  murate = anc.get_murate();
  wol = anc.get_wol();
  popsize = anc.get_popsize();
  nugen = anc.get_nugen();
  int i;
  population = new GraphI[popsize];
  nepop = new GraphI[popsize];
  w = new double[popsize];
  basic.fillv0(w, popsize);
  for (i=0; i<popsize; i++) {
    population[i].start_rng(est);
    nepop[i].start_rng(est);
    population[i].copy(anc.population[i]);
    w[i] = anc.get_w(i);
  }
  meanw = anc.return_meanw();
  maxw = anc.return_maxw();
  yaparams = true;
  yapop = true;
}

void EvolveI::print_params(ostream& sal) //prints parameters
{
  sal << "Integer weights\n";
  sal << "Number of genes: " << nugen << endl;
  sal << "Population size: " << popsize << endl;
  sal << "Mutation rate: " << murate << endl;
  sal << "Expected connectivity: " << wol << endl;
  sal << endl << endl;
}

void EvolveI::start_pop(GraphI &founder)
{
  if (!yaparams) {
    cout << "[Error]: Parameters had not been set when EvolveI::start_pop was called.\n";
    exit(1);
  }
  if (yapop) {
    cout << "[Error]: Population already constructed when EvolveI::start_pop was called.\n";
    exit(1);
  }
  population = new GraphI[popsize];
  nepop = new GraphI[popsize];
  w = new double[popsize];
  basic.fillv0(w, popsize);
  meanw =0;
  maxw = 0;
  int i;
  for (i=0; i<popsize; i++) {
    population[i].start_rng(est);
    nepop[i].start_rng(est);
    population[i].copy(founder);
    population[i].consider_mutation(murate, wol);
  }
  mean_gen_dist=0;
  mean_gen_dist_thr = 0;
  max_gen_dist = 0;
  max_gen_dist_thr = 0;
  int nucomb = (popsize*(popsize-1))/2;
  all_gdists = new double[nucomb];
  basic.fillv0(all_gdists, nucomb);
  all_gdists_thr = new double[nucomb];
  basic.fillv0(all_gdists_thr, nucomb);
  mean_ph_dist = 0;
  max_ph_dist = 0;
  all_pdists = new double[nucomb];
  yapop = true;
}

void EvolveI::start_pop(int n, int e) {
  GraphI adam(n, e, est, true);
  start_pop(adam);
  adam.clear();
}

void EvolveI::start_pop(int n, int e, int *ic, int *goal) {
  int nuic = 1, per = 1;
  int ***goals, *vsigo, **maic;
  basic.create_array(goals, nuic, per, n);
  vsigo = new int[nuic];
  basic.create_array(maic, nuic, n);
  int i, j, k;
  for (i = 0; i < nuic; i++) {
    vsigo[i] = per;
    for (j = 0; j < vsigo[i]; j++) {
      for (k = 0; k < n; k++) {
        goals[i][j][k] = goal[k];
        maic[i][k] = ic[k];
      }
    }
  }
  start_pop(n, e, maic, nuic, goals, vsigo);
  for (i = 0; i < nuic; i++) {
    for (j = 0; j < vsigo[i]; j++) {
      delete [] goals[i][j];
      delete [] maic[j];
    }
    delete [] goals[i];
  }
  delete [] maic;
  delete [] goals;
  delete [] vsigo;
  
}

void EvolveI::start_pop(int n, int e, int *ic, int **goal, int per) {
  int nuic = 1;
  int ***goals, *vsigo, **maic;
  basic.create_array(goals, nuic, per, n);
  vsigo = new int[nuic];
  basic.create_array(maic, nuic, n);
  int i, j, k;
  for (i = 0; i < nuic; i++) {
    vsigo[i] = per;
    for (j = 0; j < vsigo[i]; j++) {
      for (k = 0; k < n; k++) {
        goals[i][j][k] = goal[j][k];
        maic[i][k] = ic[k];
      }
    }
  }
  start_pop(n, e, maic, nuic, goals, vsigo);
  for (i = 0; i < nuic; i++) {
    for (j = 0; j < vsigo[i]; j++) {
      delete [] goals[i][j];
      delete [] maic[j];
    }
    delete [] goals[i];
  }
  delete [] maic;
  delete [] goals;
  delete [] vsigo;
}

void EvolveI::start_pop(int n, int e, int **maic, int nuic, int **goal) {
  int ***goals, *vsigo;
  basic.create_array(goals, nuic, 1, n);
  vsigo = new int[nuic];
  int i, k;
  for (i = 0; i < nuic; i++) {
    vsigo[i] = 1;
    for (k = 0; k < n; k++)
      goals[i][0][k] = goal[i][k];
  }
  start_pop(n, e, maic, nuic, goals, vsigo);
  for (i = 0; i < nuic; i++) {
    delete [] goals[i][0];
    delete [] goals[i];
  }
  delete [] goals;
  delete [] vsigo;
}

void EvolveI::start_pop(int n, int e, int **maic, int nuic, int ***goals, int *vsigo)
{
    if (e< n) {
      cout << "[Error]: Prohibited number of interactions in EvolveI::start_pop.\n";
      exit(1);
    }
    GraphI adam(n, est, true);
    int i,j,k,l;
    double nual;
    for (i=0; i<n; i++) {
      j = est.randint(0, n);
      if (est.toss())
        nual = 1;
      else
        nual = -1;
      adam.force_interaction(j, i, nual);
    }
    for (k = i; k < e; k++) {
      do {
        j = est.randint(0, n);
        l = est.randint(0,n);
      } while (adam.weight(l,j) != 0);
      if (est.toss())
        nual = 1;
      else
        nual = -1;
      adam.force_interaction(l,j, nual);
    }
    FitnessI lafi(est);
    double law = 0;
    int *elatra, *orden;
    bool paiso;
    GraphI vacia(est);
    do {
      if (nuic == 1) {
        adam.set_as_state(maic[0]);
        adam.find_an_attractor();
        law = lafi.strict(adam, goals[0], vsigo[0]);
      }
      else {
        adam.find_attractors(maic, nuic);
        law = lafi.strict_mult(adam, goals, nuic, vsigo);
      }
      if (law<1) {
        paiso = false;
        if ((nuic==1) && (vsigo[0] == 1) && (adam.attractor_size() == 1)) {
          elatra = new int[n];
          orden = new int[n];
          for (i=0; i < n; i++) {
            elatra[i] = adam.attractor_element(0,i);
            orden[i] = i;
          }
          if ((basic.count_in_vector(goals[0][0], n, 1) == basic.count_in_vector(elatra, n, 1)) && (basic.count_in_vector(goals[0][0], n, (-1)) == basic.count_in_vector(elatra, n, (-1)))) {
            for (i=0; i < n; i++)
              if (elatra[i] != goals[0][0][i])
                for (j=i+1; j<n; j++)
                  if ((elatra[j] == goals[0][0][i]) && (elatra[i] == goals[0][0][j])) {
                    elatra[i] = goals[0][0][i];
                    elatra[j] = goals[0][0][j];
                    orden[i] = j;
                    orden[j] = i;
                    break;
                  }
            paiso = true;
          }
        }
        if (paiso) {
          adam.transform_to_isomorph(vacia, orden);
          adam.clear();
          adam.copy(vacia);
          vacia.clear();
          delete [] elatra;
          delete [] orden;
        }
        else {
          adam.clear();
          adam.make_nw(n,true);
          for (i=0; i<n; i++) {
            j = est.randint(0, n);
            if (est.toss())
              nual = 1;
            else
              nual = -1;
            adam.force_interaction(j, i, nual);
          }
          for (k = i; k < e; k++) {
            do {
              j = est.randint(0, n);
              l = est.randint(0,n);
            } while (adam.weight(l,j) != 0);
            if (est.toss())
              nual = 1;
            else
              nual = -1;
            adam.force_interaction(l,j, nual);
          }
        }
      }
    } while(law < 1);
    start_pop(adam);
    adam.clear();

  return;
}

void EvolveI::one_generation()
{
  int i,j;
  double frog, toad;
  for (i=0; i<popsize; i++) {
    frog = est.randreal()*meanw*popsize;
    toad = 0;
    for (j=0; j<popsize; j++) {
      toad = toad + w[j];
      if (toad >= frog)
        break;
    }
    nepop[i].copy(population[j]);
    nepop[i].consider_mutation(murate, wol);
  }
  for (i=0; i < popsize; i++) {
    population[i].clear();
    population[i].copy(nepop[i]);
    nepop[i].clear();
  }
}

void EvolveI::one_generation_mp5()
{
  int i,j;
  double frog, toad;
  for (i=0; i<popsize; i++) {
    frog = est.randreal()*meanw*popsize;
    toad = 0;
    for (j=0; j<popsize; j++) {
      toad = toad + w[j];
      if (toad >= frog)
        break;
    }
    nepop[i].copy(population[j]);
    nepop[i].consider_mutation_mp5(murate);
  }
  for (i=0; i < popsize; i++) {
    population[i].clear();
    population[i].copy(nepop[i]);
    nepop[i].clear();
  }
}

void EvolveI::one_generation_mf()
{
  int i,j;
  double frog, toad;
  for (i=0; i<popsize; i++) {
    frog = est.randreal()*meanw*popsize;
    toad = 0;
    for (j=0; j<popsize; j++) {
      toad = toad + w[j];
      if (toad >= frog)
        break;
    }
    nepop[i].copy(population[j]);
    nepop[i].consider_mutation_mf(murate);
  }
  for (i=0; i < popsize; i++) {
    population[i].clear();
    population[i].copy(nepop[i]);
    nepop[i].clear();
  }
}

void EvolveI::one_generation_sex()
{
  int i,j,k;
  double frog, toad, tadp;
  for (i=0; i<popsize; i++) {
    frog = est.randreal()*meanw*popsize;
    toad = 0;
    for (j=0; j<popsize; j++) {
      toad = toad + w[j];
      if (toad >= frog)
        break;
    }
    frog = est.randreal()*meanw*popsize;
    tadp = 0;
    for (k=0; k< popsize; k++) {
      tadp = tadp + w[k];
      if (tadp >= frog)
        break;
    }
    nepop[i].mate(population[j], population[k]);
    nepop[i].consider_mutation(murate, wol);
  }
  for (i=0; i < popsize; i++) {
    population[i].clear();
    population[i].copy(nepop[i]);
    nepop[i].clear();
  }
}

void EvolveI::assign_w(int cual, double cuanto)
{
  if ((cual >= popsize) || (cual < 0)) {
    cout << "[Error]: Individual does not exist. EvolveI::assign_w.\n";
    exit(1);
  }
  w[cual] = cuanto;
}

void EvolveI::calc_meanw()
{
  meanw = basic.get_mean(w, popsize);
  maxw = basic.find_max(w, popsize);
}

double EvolveI::return_meanw()
{
  return meanw;
}

double EvolveI::return_maxw()
{
  return maxw;
}

void EvolveI::optima(set<int> &opt)
{
  opt.clear();
  int i;
  for (i=0; i < popsize; i++)
    if (w[i]==maxw)
      opt.insert(i);
}

void EvolveI::optima_strict(set<int> &opt)
{
  opt.clear();
  int i;
  for (i=0; i < popsize; i++)
    if (w[i]==1)
      opt.insert(i);
}

int EvolveI::num_optima()
{
  int i, res=0;
  for (i=0; i < popsize; i++)
    if (w[i]==maxw)
      res++;
  return res;
}

int EvolveI::num_optima_strict()
{
  int i, res=0;
  for (i=0; i < popsize; i++)
    if (w[i]==1)
      res++;
  return res;
}

void EvolveI::clear()
{
  int i;
  if (yapop) {
    for (i=0; i<popsize; i++)
      population[i].clear();
    delete [] population;
    delete [] nepop;
    delete [] all_gdists_thr;
    delete [] all_gdists;
    delete [] w;
    delete [] all_pdists;
  }
  yapop = false;
  yaparams = false;
  if (yalop)
    clear_lop();
}

void EvolveI::clear_lop()
{
  int i,j,k;
  for (i=0; i < vanop; i++) {
    for (j=0; j < nulop[i]; j++) {
      for (k=0; k < perlop[i][j]; k++)
        delete [] last_opt_phen[i][j][k];
      delete [] last_opt_phen[i][j];
    }
    delete [] perlop[i];
    delete [] last_opt_phen[i];
  }
  delete [] nulop;
  delete [] perlop;
  delete [] last_opt_phen;
  vanop = 0;
  yalop = false;
}

double EvolveI::get_murate()
{
  return murate;
}

double EvolveI::get_wol()
{
  return wol;
}

void EvolveI::new_wol(double nwol) {
  wol = nwol;
}

int EvolveI::get_popsize()
{
  return popsize;
}

int EvolveI::get_nugen()
{
  return nugen;
}

double EvolveI::get_w(int cual) {
  return w[cual];
}

int EvolveI::group_by_adjmat(int *vecid) {
  basic.fillvm1(vecid, popsize);
  int i, j, voy = 0;
  for (i = 0; i < popsize; i++) {
    if (vecid[i] == -1) {
      vecid[i] = voy;
      for (j = i+1; j < popsize; j++)
        if (population[i].equal_adjmat(population[j]))
          vecid[j] = voy;
      voy++;
    }
  }
  return voy;
}

int EvolveI::number_of_genotypes()
{
  int res = 1;
  int i,j,k,l;
  bool *esnu, igu;
  esnu = new bool[popsize];
  basic.fillv0(esnu, popsize);
  esnu[0] = true;
  for (i=1; i < popsize; i++) {
    esnu[i] = true;
    for (j=0; j<i; j++) {
      if (esnu[j]) {
        igu = true;
        for (k=0; k < nugen; k++) {
          for (l=0; l<nugen; l++) {
            if (population[i].weight(l,k) != population[j].weight(l,k)) {
              igu = false;
              break;
            }
          }
          if (!igu)
            break;
        }
        if (igu) {
          esnu[i] = false;
          break;
        }
      }
    }
    if (esnu[i])
      res++;
  }
  delete [] esnu;
  return res;
}

int EvolveI::number_of_genotypes_thr()
{
  int res = 1;
  int i,j,k,l;
  bool *esnu, igu;
  double prod;
  esnu = new bool[popsize];
  basic.fillv0(esnu, popsize);
  esnu[0] = true;
  for (i=1; i < popsize; i++) {
    esnu[i] = true;
    for (j=0; j<i; j++) {
      if (esnu[j]) {
        igu = true;
        for (k=0; k < nugen; k++) {
          for (l=0; l<nugen; l++) {
            prod = population[i].weight(l,k)*population[j].weight(l,k);
            if ((prod < 0) || ((prod == 0) && (population[i].weight(l,k) != population[j].weight(l,k)))) {
              igu = false;
              break;
            }
          }
          if (!igu)
            break;
        }
        if (igu) {
          esnu[i] = false;
          break;
        }
      }
    }
    if (esnu[i])
      res++;
  }
  delete [] esnu;
  return res;
}

int EvolveI::number_of_genotypes(set<int> &opt)
{
  set<int>::iterator It;
  int tam = opt.size();
  int res,i,j,k,l;
  bool *esnu, igu;
  int *arte;
  if (tam >0) {
    arte = new int[tam];
    i=0;
    for (It = opt.begin(); It != opt.end(); It++) {
      arte[i] = *It;
      i++;
    }
    esnu = new bool[tam];
    basic.fillv0(esnu, tam);
    esnu[0] = true;
    res = 1;
    for (i=1; i<tam; i++) {
      esnu[i] = true;
      for (j=0; j < i; j++) {
        if (esnu[j]) {
          igu = true;
          for (k=0; k < nugen; k++) {
            for (l=0; l<nugen; l++) {
              if (population[arte[i]].weight(l,k) != population[arte[j]].weight(l,k)) {
                igu = false;
                break;
              }
            }
            if (!igu)
              break;
          }
          if (igu) {
            esnu[i] = false;
            break;
          }
        }
      }
      if (esnu[i])
        res++;
    }
    delete [] arte;
    delete [] esnu;
  }
  else
    res = 0;
  return res;
}

int EvolveI::number_of_genotypes_thr(set<int> &opt)
{
  set<int>::iterator It;
  int tam = opt.size();
  int res,i,j,k,l;
  bool *esnu, igu;
  double prod;
  int *arte;
  if (tam >0) {
    arte = new int[tam];
    i=0;
    for (It = opt.begin(); It != opt.end(); It++) {
      arte[i] = *It;
      i++;
    }
    esnu = new bool[tam];
    basic.fillv0(esnu, tam);
    esnu[0] = true;
    res = 1;
    for (i=1; i<tam; i++) {
      esnu[i] = true;
      for (j=0; j < i; j++) {
        if (esnu[j]) {
          igu = true;
          for (k=0; k < nugen; k++) {
            for (l=0; l<nugen; l++) {
              prod = population[arte[i]].weight(l,k)*population[arte[j]].weight(l,k);
              if ((prod < 0) || ((prod == 0) && (population[arte[i]].weight(l,k) != population[arte[j]].weight(l,k)))) {
                igu = false;
                break;
              }
            }
            if (!igu)
              break;
          }
          if (igu) {
            esnu[i] = false;
            break;
          }
        }
      }
      if (esnu[i])
        res++;
    }
    delete [] arte;
  }
  else
    res = 0;
  return res;
}

void EvolveI::genotype_statistics()
{
  int i,j,k,l,van=0;
  double sum = 0, sumth = 0, prod;
  for (i=1; i<popsize; i++) {
    for (j=0; j<i; j++) {
      sum = 0;
      sumth = 0;
      for (k=0; k < nugen; k++)
        for (l=0; l < nugen; l++) {
          sum = sum + abs(population[i].weight(l,k) - population[j].weight(l,k));
          prod = population[i].weight(l,k) * population[j].weight(l,k);
          if (prod < 0)
            sumth = sumth + 2;
          else {
            if ((prod > 0) || (population[i].weight(l,k) == population[j].weight(l,k)))
              sumth = sumth + 0;
            else if (prod == 0)
              sumth = sumth + 1;
          }
        }
      if ((population[i].number_of_edges()+population[j].number_of_edges()+sum)==0)
        all_gdists[van] = 0;
      else
        all_gdists[van] = sum/double(population[i].number_of_edges()+population[j].number_of_edges());
      if ((population[i].number_of_edges()+population[j].number_of_edges()+sumth)==0)
        all_gdists_thr[van] = 0;
      else
        all_gdists_thr[van] = sumth/double(population[i].number_of_edges()+population[j].number_of_edges());
      van++;
    }
  }
  if (van != ((popsize*(popsize-1))/2)) {
    cout << "[Error]: Wrong number of combinations in EvolveI::genotype statistics.\n";
    exit(1);
  }
  mean_gen_dist = basic.get_mean(all_gdists, van);
  mean_gen_dist_thr = basic.get_mean(all_gdists_thr, van);
  max_gen_dist = basic.find_max(all_gdists, van);
  max_gen_dist_thr = basic.find_max(all_gdists_thr, van);
}

void EvolveI::genotype_statistics(set<int> &opt, double& megd, double& megdth, double& magd, double& magdth, double *agd, double *agdth)
{
  int tam = opt.size();
  int nucombi = (tam*(tam-1))/2;
  set<int>::iterator It;
  int i,j,k,l,van=0;
  double sum = 0, sumth = 0, prod;
  int *arte;
  if (tam == 0) {
    cout << "[Error]: No genotypes in set when EvolveI::genotype_statistics was called.\n";
    exit(1);
  }
  arte = new int[tam];
  i=0;
  for (It = opt.begin(); It != opt.end(); It++) {
    arte[i] = *It;
    i++;
  }
  for (i=1; i<tam; i++) {
    for (j=0; j<i; j++) {
      sum = 0;
      sumth = 0;
      for (k=0; k < nugen; k++)
        for (l=0; l < nugen; l++) {
          sum = sum + abs(population[arte[i]].weight(l,k) - population[arte[j]].weight(l,k));
          prod = population[arte[i]].weight(l,k) * population[arte[j]].weight(l,k);
          if (prod < 0)
            sumth = sumth + 2;
          else {
            if ((prod > 0) || (population[arte[i]].weight(l,k) == population[arte[j]].weight(l,k)))
              sumth = sumth + 0;
            else if (prod == 0)
              sumth = sumth + 1;
          }
        }
      if ((population[arte[i]].number_of_edges()+population[arte[j]].number_of_edges()+sum)==0)
        agd[van] = 0;
      else
        agd[van] = sum/double(population[arte[i]].number_of_edges()+population[arte[j]].number_of_edges());
      if ((population[arte[i]].number_of_edges()+population[arte[j]].number_of_edges()+sumth)==0)
        agdth[van]=0;
      else
        agdth[van] = sumth/double(population[arte[i]].number_of_edges()+population[arte[j]].number_of_edges());
      van++;
    }
  }
  if (van != nucombi) {
    cout << "[Error]: Wrong number of combinations in EvolveI::genotype_statistics.\n";
    exit(1);
  }
  delete [] arte;
  megd = basic.get_mean(agd, van);
  megdth = basic.get_mean(agdth, van);
  magd = basic.find_max(agd, van);
  magdth = basic.find_max(agdth, van);
}

int EvolveI::number_of_phenotypes()
{
  if (!population[popsize-1].attractor_exists()) {
    cout << "[Error]: Phenotypes had not been determined when EvolveI::number_of_phenotypes was called.\n";
    exit(1);
  }
  int *tamat;
  tamat = new int[popsize];
  int i,j,k,res =0;
  bool *esnu;
  for (i=0; i < popsize; i++) {
    tamat[i] = population[i].attractor_size();
    if (tamat[i] < 1) {
      cout << "[Error]: Prohibited attractor size in EvolveI::number_of_phenotypes.\n";
      exit(1);
    }
  }
  int ***losat;
  losat = new int**[popsize];
  for (i=0; i < popsize; i++) {
    losat[i] = new int*[tamat[i]];
    for (j=0; j < tamat[i]; j++)
      losat[i][j] =new int[nugen];
  }
  for (i=0; i<popsize; i++)
    for (j=0; j<tamat[i]; j++)
      for (k=0; k < nugen; k++)
        losat[i][j][k] = population[i].attractor_element(j,k);
  esnu = new bool[popsize];
  basic.fillv0(esnu, popsize);
  for (i=0; i<popsize; i++) {
    esnu[i] = true;
    for (j=0; j<i; j++)
      if (esnu[j]) {
        if (basic.eqmatrix_rot(losat[i], tamat[i], nugen, losat[j], tamat[j], nugen)) {
          esnu[i] = false;
          break;
        }
      }
    if (esnu[i])
      res++;
  }
  delete [] esnu;
  int van = 0;
  for (i=1; i < popsize; i++)
    for (j=0; j < i; j++) {
      all_pdists[van] = basic.dist_matrices_rot(losat[i], tamat[i], nugen, losat[j], tamat[j], nugen);
      van++;
    }
  if (van != ((popsize*(popsize-1))/2)) {
    cout << "[Error]: Wrong number of combinations in EvolveI::number_of_phenotypes.\n";
    exit(1);
  }
  mean_ph_dist = basic.get_mean(all_pdists, van);
  max_ph_dist = basic.find_max(all_pdists, van);
  for (i=0; i<popsize; i++) {
    for (j=0; j<tamat[i]; j++)
      delete [] losat[i][j];
    delete [] losat[i];
  }
  delete [] losat;
  delete [] tamat;
  return res;
}

int EvolveI::number_of_phenotypes(set<int> &opt, double& mepd, double& mapd, double *apd)
{
  int tam = opt.size();
  set<int>::iterator It;
  int i,j,k,van=0, res=0;
  bool *esnu;
  int *arte;
  if (tam == 0) {
    cout << "[Error]: No genotypes in set when EvolveI::genotype_statistics was called.\n";
    exit(1);
  }
  arte = new int[tam];
  i=0;
  for (It = opt.begin(); It != opt.end(); It++) {
    arte[i] = *It;
    i++;
  }
  if (!population[arte[tam-1]].attractor_exists()) {
    cout << "[Error]: Phenotypes had not been determined when EvolveI::number_of_phenotypes was called.\n";
    delete [] arte;
    exit(1);
  }
  int *tamat;
  tamat = new int[tam];
  for (i=0; i < tam; i++) {
    tamat[i] = population[arte[i]].attractor_size();
    if (tamat[i] < 1) {
      cout << "[Error]: Prohibited attractor size in EvolveI::number_of_phenotypes.\n";
      exit(1);
    }
  }
  int ***losat;
  losat = new int**[tam];
  for (i=0; i < tam; i++) {
    losat[i] = new int*[tamat[i]];
    for (j=0; j < tamat[i]; j++)
      losat[i][j] =new int[nugen];
  }
  for (i=0; i<tam; i++)
    for (j=0; j<tamat[i]; j++)
      for (k=0; k < nugen; k++)
        losat[i][j][k] = population[arte[i]].attractor_element(j,k);
  esnu = new bool[tam];
  basic.fillv0(esnu, tam);
  for (i=0; i<tam; i++) {
    esnu[i] = true;
    for (j=0; j<i; j++)
      if (esnu[j]) {
        if (basic.eqmatrix_rot(losat[i], tamat[i], nugen, losat[j], tamat[j], nugen)) {
          esnu[i] = false;
          break;
        }
      }
    if (esnu[i])
      res++;
  }
  delete [] esnu;
  for (i=1; i < tam; i++)
    for (j=0; j < i; j++) {
      apd[van] = basic.dist_matrices_rot(losat[i], tamat[i], nugen, losat[j], tamat[j], nugen);
      van++;
    }
  if (van != ((tam*(tam-1))/2)) {
    cout << "[Error]: Wrong number of combinations in EvolveI::number_of_phenotypes.\n";
    exit(1);
  }
  mepd = basic.get_mean(apd, van);
  mapd = basic.find_max(apd, van);
  for (i=0; i<tam; i++) {
    for (j=0; j<tamat[i]; j++)
      delete [] losat[i][j];
    delete [] losat[i];
  }
  delete [] losat;
  delete [] arte;
  delete [] tamat;
  return res;
}

int EvolveI::number_of_phenotypes(int *nuper)
{
  if (!population[popsize-1].attractor_exists()) {
    cout << "[Error]: Phenotypes had not been determined when EvolveI::number_of_phenotypes was called.\n";
    exit(1);
  }
  basic.fillv0(nuper, 5);
  int *tamat;
  tamat = new int[popsize];
  int i,j,k,res =0;
  bool *esnu;
  for (i=0; i < popsize; i++) {
    tamat[i] = population[i].attractor_size();
    if (tamat[i] < 1) {
      cout << "[Error]: Prohibited attractor size in EvolveI::number_of_phenotypes.\n";
      exit(1);
    }
    if (tamat[i] < 4)
      nuper[tamat[i]]++;
    else
      nuper[4]++;
  }
  int ***losat;
  losat = new int**[popsize];
  for (i=0; i < popsize; i++) {
    losat[i] = new int*[tamat[i]];
    for (j=0; j < tamat[i]; j++)
      losat[i][j] =new int[nugen];
  }
  for (i=0; i<popsize; i++)
    for (j=0; j<tamat[i]; j++)
      for (k=0; k < nugen; k++)
        losat[i][j][k] = population[i].attractor_element(j,k);
  esnu = new bool[popsize];
  basic.fillv0(esnu, popsize);
  for (i=0; i<popsize; i++) {
    esnu[i] = true;
    for (j=0; j<i; j++)
      if (esnu[j]) {
        if (basic.eqmatrix_rot(losat[i], tamat[i], nugen, losat[j], tamat[j], nugen)) {
          esnu[i] = false;
          break;
        }
      }
    if (esnu[i])
      res++;
  }
  delete [] esnu;
  int van = 0;
  for (i=1; i < popsize; i++)
    for (j=0; j < i; j++) {
      all_pdists[van] = basic.dist_matrices_rot(losat[i], tamat[i], nugen, losat[j], tamat[j], nugen);
      van++;
    }
  if (van != ((popsize*(popsize-1))/2)) {
    cout << "[Error]: Wrong number of combinations in EvolveI::number_of_phenotypes.\n";
    exit(1);
  }
  mean_ph_dist = basic.get_mean(all_pdists, van);
  max_ph_dist = basic.find_max(all_pdists, van);
  for (i=0; i<popsize; i++) {
    for (j=0; j<tamat[i]; j++)
      delete [] losat[i][j];
    delete [] losat[i];
  }
  delete [] losat;
  delete [] tamat;
  return res;
}

int EvolveI::number_of_phenotypes(set<int> &opt, double& mepd, double& mapd, double *apd, int *nuper)
{
  int tam = opt.size();
  set<int>::iterator It;
  int i,j,k,van=0, res=0;
  int *arte;
  bool *esnu;
  if (tam == 0) {
    cout << "[Error]: No genotypes in set when EvolveI::genotype_statistics was called.\n";
    exit(1);
  }
  arte = new int[tam];
  i=0;
  for (It = opt.begin(); It != opt.end(); It++) {
    arte[i] = *It;
    i++;
  }
  if (!population[arte[tam-1]].attractor_exists()) {
    cout << "[Error]: Phenotypes had not been determined when EvolveI::number_of_phenotypes was called.\n";
    delete [] arte;
    exit(1);
  }
  basic.fillv0(nuper, 5);
  int *tamat;
  tamat = new int[tam];
  for (i=0; i < tam; i++) {
    tamat[i] = population[arte[i]].attractor_size();
    if (tamat[i] < 1) {
      cout << "[Error]: Prohibited attractor size in EvolveI::number_of_phenotypes.\n";
      exit(1);
    }
    if (tamat[i] < 4)
      nuper[tamat[i]]++;
    else
      nuper[4]++;
  }
  int ***losat;
  losat = new int**[tam];
  for (i=0; i < tam; i++) {
    losat[i] = new int*[tamat[i]];
    for (j=0; j < tamat[i]; j++)
      losat[i][j] =new int[nugen];
  }
  for (i=0; i<tam; i++)
    for (j=0; j<tamat[i]; j++)
      for (k=0; k < nugen; k++)
        losat[i][j][k] = population[arte[i]].attractor_element(j,k);
  esnu = new bool[tam];
  basic.fillv0(esnu, tam);
  for (i=0; i<tam; i++) {
    esnu[i] = true;
    for (j=0; j<i; j++)
      if (esnu[j]) {
        if (basic.eqmatrix_rot(losat[i], tamat[i], nugen, losat[j], tamat[j], nugen)) {
          esnu[i] = false;
          break;
        }
      }
    if (esnu[i])
      res++;
  }
  delete [] esnu;
  for (i=1; i < tam; i++)
    for (j=0; j < i; j++) {
      apd[van] = basic.dist_matrices_rot(losat[i], tamat[i], nugen, losat[j], tamat[j], nugen);
      van++;
    }
  if (van != ((tam*(tam-1))/2)) {
    cout << "[Error]: Wrong number of combinations in EvolveI::number_of_phenotypes.\n";
    exit(1);
  }
  mepd = basic.get_mean(apd, van);
  mapd = basic.find_max(apd, van);
  for (i=0; i<tam; i++) {
    for (j=0; j<tamat[i]; j++)
      delete [] losat[i][j];
    delete [] losat[i];
  }
  delete [] losat;
  delete [] arte;
  delete [] tamat;
  return res;
}

void EvolveI::build_last_opt_phen(int vagener, int lastgener, int samrat, set<int> &opt)
{
  if (vagener < (lastgener - (samrat*20)))
    return;
  int nufop = 21;
  if (!yalop) {
    vanop = 0;
    generini = vagener;
    nulop = new int[nufop];
    perlop = new int*[nufop];
    last_opt_phen = new int***[nufop];
    yalop = true;
  }
  if (((vagener-generini)%samrat) == 0) {
    nulop[vanop] = opt.size();
    if (nulop[vanop] >0) {
      int i,j,k;
      int *arte;
      arte = new int[nulop[vanop]];
      set<int>::iterator it;
      perlop[vanop] = new int[nulop[vanop]];
      last_opt_phen[vanop] = new int**[nulop[vanop]];
      i = 0;
      for (it = opt.begin(); it != opt.end(); it++) {
        arte[i] = *it;
        i++;
      }
      for (i=0; i < nulop[vanop]; i++)
        perlop[vanop][i] = population[arte[i]].attractor_size();
      for (i=0; i < nulop[vanop]; i++) {
        last_opt_phen[vanop][i] = new int*[perlop[vanop][i]];
        for (j=0; j < perlop[vanop][i]; j++)
          last_opt_phen[vanop][i][j] = new int[nugen];
        for (j=0; j < perlop[vanop][i]; j++)
          for (k=0; k < nugen; k++)
            last_opt_phen[vanop][i][j][k] = population[arte[i]].attractor_element(j,k);
      }
      delete [] arte;
    }
    vanop++;
  }
  return;
}

double EvolveI::analyze_last_opt_phen(int *hasdo21, double& resnor)
{
  if (!yalop) {
    cout << "[Error]: EvolveI::build_last_opt_phen must be called before EvolveI::analyze_last_opt_phen.\n";
    exit(1);
  }
  int nufop = 21;
  bool **mabo;
  int i,j,k;
  mabo = new bool*[nufop];
  for (i=0; i < nufop; i++)
    mabo[i] = new bool[nulop[i]];
  for (i = 0; i < nufop; i++)
    for (j=0; j < nulop[i]; j++)
      mabo[i][j] = false;
  bool tubo;
  for (i=0; i<nulop[0]; i++) {
    for (j=1; j < nufop; j++) {
      tubo = false;
      for (k=0; k < nulop[j]; k++)
        if (!mabo[j][k]) {
          if (basic.eqmatrix_rot(last_opt_phen[0][i], perlop[0][i], nugen, last_opt_phen[j][k], perlop[j][k], nugen)) {
            mabo[j][k] = true;
            tubo = true;
            break;
          }
        }
      if (!tubo)
        break;
    }
  }
  basic.fillv0(hasdo21, nufop);
  hasdo21[0] = nulop[0];
  double res;
  if (hasdo21[0] > 0) {
    int tonutru = 0;
    int dentopos = 0;
    int denrepos = 0;
    for (i=1; i < nufop; i++) {
      hasdo21[i] = basic.count_in_vector(mabo[i], nulop[i], true);
      tonutru = tonutru + hasdo21[i];
      dentopos = dentopos + nulop[i];
      denrepos = denrepos + basic.find_min(nulop, i+1);
    }
    resnor = tonutru/double(denrepos);
    res = tonutru/double(dentopos);
  }
  else {
    resnor = -1;
    res = -1;
  }
  for (i=0; i < nufop; i++)
    delete [] mabo[i];
  delete [] mabo;
  return res;
}


