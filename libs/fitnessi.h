#ifndef FITNESSI_H
#define FITNESSI_H

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

class FitnessI
{
public:
	FitnessI();
  FitnessI(Alea& jacta);
  void start_rng(Alea& jacta);
  double strict_partial(GraphI &red, int *goal, int selgenes);
  double strict(GraphI &red, int *goal);
  double strict(GraphI &red, int **goal, int lencyc);
  double strict_mult(GraphI &red, int **goals, int nugoals);//only fixed-points
  double strict_mult(GraphI &red, int ***goals, int nugoals, int *gosize);
  double fixed_point(GraphI &red);
  double fixed_points(GraphI &red);
  double period(GraphI &red, int p);
  double periods(GraphI &red, int p);
  double nz(GraphI &red);
  double nz_mult(GraphI &red);
  double fixed_point_nz(GraphI &red);
  double fixed_points_nz(GraphI &red);
  double period_nz(GraphI &red, int p);
  double periods_nz(GraphI &red, int p);
  double linear(GraphI &red, int *goal);
  double linear(GraphI &red, int **goal, int lencyc);
  double linear_nz(GraphI &red, int *goal);
  double linear_nz(GraphI &red, int **goal, int lencyc);
  double pl(GraphI &red, int *goal, double expo);
  double pl(GraphI &red, int **goal, int lencyc, double expo);
  double pl_nz(GraphI &red, int *goal, double expo);
  double pl_nz(GraphI &red, int **goal, int lencyc, double expo);
  double frac1s_strict(GraphI &red, double frad);
  double frac1s_strict_nz(GraphI &red, double frad);
  double frac1s_pl(GraphI &red, double frad, double expo);
  double frac1s_pl_nz(GraphI &red, double frad, double expo); 
  double strict_pathlength(GraphI &red, int opt);
  double greater_pathlength(GraphI &red, int opt, double penalone);
  double greater_pathlength(GraphI &red, int opt);
  double greater_pathlength_nz(GraphI &red, int opt, double penalone);
  double greater_pathlength_nz(GraphI &red, int opt);
  double greater_atsize(GraphI &red, int opt, double penalone);
  double greater_atsize(GraphI &red, int opt);
  double greater_atsize_nz(GraphI &red, int opt, double penalone);
  double greater_atsize_nz(GraphI &red, int opt);

  //require comb
  void linear_mult(GraphI &red, int **goals, int nugoals, double *ws); //only fixed-points
  void linear_mult(GraphI &red, int ***goals, int nugoals, int *gosize, double *ws);
  void linear_mult_nz(GraphI &red, int **goals, int nugoals, double *ws); //only fixed-points
  void linear_mult_nz(GraphI &red, int ***goals, int nugoals, int *gosize, double *ws);
  void pl_mult(GraphI &red, int **goals, int nugoals, double *ws, double expo); //only fixed-points
  void pl_mult(GraphI &red, int ***goals, int nugoals, int *gosize, double *ws, double expo);
  void pl_mult_nz(GraphI &red, int **goals, int nugoals, double *ws, double expo); //only fixed-points
  void pl_mult_nz(GraphI &red, int ***goals, int nugoals, int *gosize, double *ws, double expo);
  void frac1s_strict_mult(GraphI &red, double fad, double *ws);
  void frac1s_strict_mult_nz(GraphI &red, double fad, double *ws);
  void frac1s_pl_mult(GraphI &red, double frad, double *ws, double expo);
  void frac1s_pl_mult_nz(GraphI &red, double frad, double *ws, double expo);
  
  double distance(int *goal, int *atr, int tam);
  double distance(int *goal, int **atr, int numatr, int tam);
  double distance(int **goal, int nugs, int **atr, int numatr, int tam);
  
  double comb_mult(double *ws, int tam);
  double comb_mult(double *ws, int tam, double *weights);
  double comb_add(double *ws, int tam);
  double comb_add(double *ws, int tam, double *weights);

  double mult_multi_goals_under_pert(GraphI& red, int numofgoals, int numcells, int **ci, int ***goal, int *periodgoal, double pertrate, double selcoefperg, double multexp);

  double additive_multi_goals_under_pert(GraphI& red, int numofgoals, int numcells, int **ci, int ***goal, int *periodgoal, double pertrate, double selcoefperg, double multexp);
  double under_pert(GraphI& red, int numcells, int *ci, int **goal, int periodgoal, double pertrate, double selcoefperg, double multexp);
  
private:
  Alea est;
  Basics basic;
  double distance_aux(int **goal, int nugs, int **atr, int numatr, int tam);
  bool zeros_in_attractor(GraphI &red);
  bool zeros_in_attractor(GraphI &red, int wh); //one in multiple
  double get_frac_1s(GraphI &red);
  double get_frac_1s(GraphI &red, int wh);

};

#endif
