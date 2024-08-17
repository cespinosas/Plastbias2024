#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <set>
#include <list>
#include "basics.h"
#include "fitnessi.h"
#include <string>
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

FitnessI::FitnessI()
{
}

FitnessI::FitnessI(Alea& jacta) //creates instance of class FitnessI and assigns rng to jacta
{
  start_rng(jacta);
}

void FitnessI::start_rng(Alea& jacta) //assigns rng to jacta
{
  est = jacta;
  basic.start_rng(est);
}

double FitnessI::strict_partial(GraphI &red, int *goal, int selgenes)
{
  if (!red.attractor_exists()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::strict_partial was called.\n";
    exit(1);
  }
  double w= 1;
  int i, j;
  for (i = 0; i < red.attractor_size(); i++)
    for (j = 0; j < selgenes; j++)
      if (red.attractor_element(i,j) != goal[j]) {
        w=0;
        break;
      }
  return w;
}

double FitnessI::strict(GraphI &red, int *goal) //strict fitness function
{
  if (!red.attractor_exists()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::strict was called.\n";
    exit(1);
  }
  double w= 0;
  if (red.attractor_size() == 1) {
    int i;
    w=1;
    for (i=0; i< red.number_of_nodes(); i++)
      if (goal[i] != red.attractor_element(0,i)) {
        w = 0;
        break;
      }
  }
  return w;
}

double FitnessI::strict(GraphI &red, int **goal, int lencyc)
{
  if (!red.attractor_exists()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::strict was called.\n";
    exit(1);
  }
  double w= 0;
  if (red.attractor_size() == lencyc) {
    int ini,i,j;
    int **elatra;
    elatra = new int*[lencyc];
    for (i=0; i< lencyc; i++)
      elatra[i] = new int[red.number_of_nodes()];
    for (i=0; i< lencyc; i++)
      for (j=0; j< red.number_of_nodes(); j++)
        elatra[i][j] = red.attractor_element(i,j);
    ini = -1;
    for (i=0; i< lencyc; i++)
      if (basic.eqvec(goal[0], red.number_of_nodes(), elatra[i], red.number_of_nodes())) {
        ini = i;
        break;
      }
    if (ini >= 0) {
      w = 1;
      for (i=0; i < lencyc; i++)
        if (!basic.eqvec(goal[i], red.number_of_nodes(), elatra[(ini+i)%lencyc], red.number_of_nodes())) {
          w = 0;
          break;
        }
    }
    for (i=0; i< lencyc; i++)
      delete [] elatra[i];
    delete [] elatra;
  }
  return w;
}

double FitnessI::strict_mult(GraphI &red, int **goals, int nugoals)
{
  if (!red.attractors_exist()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::strict_mult was called.\n";
    exit(1);
  }
  double w= 0;
  int i,j;
  bool salte = false;
  for (i=0; i < red.number_of_attractors(); i++)
    if (red.attractor_size(i) != 1) {
      salte = true;
      break;
    }
  if ((!salte) && (nugoals == red.number_of_attractors())) {
    w = 1;
    for (i = 0; i < nugoals; i++) {
      for (j= 0; j < red.number_of_nodes(); j++)
        if (goals[i][j] != red.attractor_element(i,0,j)) {
          w = 0;
          salte = true;
          break;
        }
      if (salte)
        break;
    }
  }
  return w;
}

double FitnessI::strict_mult(GraphI &red, int ***goals, int nugoals, int *gosize)
{
  if (!red.attractors_exist()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::strict_mult was called.\n";
    exit(1);
  }
  double w= 0;
  int i,j,k,ini;
  bool salte = false;
  if (nugoals == red.number_of_attractors()) {
    for (i=0; i < red.number_of_attractors(); i++)
      if (red.attractor_size(i) != gosize[i]) {
        salte = true;
        break;
      }
    if (!salte) {
      int **elatra;
      for (i=0; i<nugoals; i++) {
        elatra = new int*[gosize[i]];
        for (j=0; j < gosize[i]; j++)
          elatra[j] = new int[red.number_of_nodes()];
        for (j=0; j < gosize[i]; j++)
          for (k= 0; k < red.number_of_nodes(); k++)
            elatra[j][k] = red.attractor_element(i,j,k);
        ini = -1;
        w = 0;
        for (j= 0; j < gosize[i]; j++)
          if (basic.eqvec(goals[i][0], red.number_of_nodes(), elatra[j], red.number_of_nodes())) {
            ini = j;
            break;
          }
        if (ini >=0) {
          w = 1;
          for (j= 0; j < gosize[i]; j++)
            if (!basic.eqvec(goals[i][j], red.number_of_nodes(), elatra[(ini+j)%gosize[i]], red.number_of_nodes())) {
              w=0;
              break;
            }
        }
        for (j= 0; j < gosize[i]; j++)
          delete [] elatra[j];
        delete [] elatra;
        if (w==0)
          break;
      }
    }
  }
  return w;
}

double FitnessI::fixed_point(GraphI &red)
{
  if (!red.attractor_exists()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::fixed_point was called.\n";
    exit(1);
  }
  double w= 0;
  if (red.attractor_size() == 1)
    w = 1;
  return w;
}

double FitnessI::fixed_points(GraphI &red)
{
  if (!red.attractors_exist()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::fixed_points was called.\n";
    exit(1);
  }
  double w = 1;
  int i;
  for (i = 0; i < red.number_of_attractors(); i++)
    if (red.attractor_size(i) != 1) {
      w = 0;
      break;
    }
  return w;
}

double FitnessI::period(GraphI &red, int p)
{
  if (!red.attractor_exists()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::period was called.\n";
    exit(1);
  }
  if (p==1) {
    cout << "[Error]: FitnessI::period not defined for fixed-points.\n";
    exit(1);
  }
  double w= 0;
  if (red.attractor_size() == p)
    w = 1;
  return w;
}

double FitnessI::periods(GraphI &red, int p)
{
  if (!red.attractors_exist()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::periods was called.\n";
    exit(1);
  }
  if (p==1) {
    cout << "[Error]: FitnessI::periods not defined for fixed-points.\n";
    exit(1);
  }
  double w = 1;
  int i;
  for (i = 0; i < red.number_of_attractors(); i++)
    if (red.attractor_size(i) != p) {
      w = 0;
      break;
    }
  return w;
}

//no-zeros
double FitnessI::nz(GraphI &red)
{
  if (!red.attractor_exists()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::nz was called.\n";
    exit(1);
  }
  double w= 1;
  int i,j;
  for (j=0; j < red.attractor_size(); j++) {
    for (i=0; i < red.number_of_nodes(); i++)
      if (red.attractor_element(j,i)==0) {
        w=0;
        break;
      }
    if (w==0)
      break;
  }
  return w;
}


double FitnessI::nz_mult(GraphI &red)
{
  if (!red.attractors_exist()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::nz_mult was called.\n";
    exit(1);
  }
  double w= 1;
  int i,j,k;
  for (k=0; k < red.number_of_attractors(); k++) {
    for (j=0; j < red.attractor_size(k); j++) {
      for (i=0; i < red.number_of_nodes(); i++)
        if (red.attractor_element(k,j,i)==0) {
          w=0;
          break;
        }
      if (w==0)
        break;
    }
    if (w==0)
      break;
  }
  return w;
}

double FitnessI::fixed_point_nz(GraphI &red)
{
  if (!red.attractor_exists()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::fixed_point_nz was called.\n";
    exit(1);
  }
  double w= 0;
  if (red.attractor_size() == 1) {
    w = 1;
    int i;
    for (i=0; i < red.number_of_nodes(); i++)
      if (red.attractor_element(0,i) == 0) {
        w= 0;
        break;
      }
  }
  return w;
}

double FitnessI::fixed_points_nz(GraphI &red)
{
  if (!red.attractors_exist()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::fixed_points_nz was called.\n";
    exit(1);
  }
  double w = 1;
  int i,j;
  for (i = 0; i < red.number_of_attractors(); i++) {
    if (red.attractor_size(i) != 1) {
      w = 0;
      break;
    }
    for (j= 0; j < red.number_of_nodes(); j++)
      if (red.attractor_element(i,0,j) == 0) {
        w = 0;
        break;
      }
    if (w== 0)
      break;
  }
  return w;
}

double FitnessI::period_nz(GraphI &red, int p)
{
  if (!red.attractor_exists()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::period_nz was called.\n";
    exit(1);
  }
  if (p==1) {
    cout << "[Error]: FitnessI::period_nz not defined for fixed-points.\n";
    exit(1);
  }
  double w= 0;
  int j,k;
  if (red.attractor_size() == p) {
    w = 1;
    for (j= 0; j < p; j++) {
      for (k=0; k < red.number_of_nodes(); k++)
        if (red.attractor_element(j,k) == 0) {
          w=0;
          break;
        }
      if (w== 0)
        break;
    }
  }
  return w;
}

double FitnessI::periods_nz(GraphI &red, int p)
{
  if (!red.attractors_exist()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::periods_nz was called.\n";
    exit(1);
  }
  if (p==1) {
    cout << "[Error]: FitnessI::periods_nz not defined for fixed-points.\n";
    exit(1);
  }
  double w = 1;
  int i,j,k;
  for (i = 0; i < red.number_of_attractors(); i++) {
    if (red.attractor_size(i) != p) {
      w = 0;
      break;
    }
    for (j= 0; j < p; j++) {
      for (k= 0; k < red.number_of_nodes(); k++) {
        if (red.attractor_element(i,j,k) == 0) {
          w = 0;
          break;
        }
      }
      if (w== 0)
        break;
    }
    if (w== 0)
      break;
  }
  return w;
}

double FitnessI::linear(GraphI &red, int *goal)
{
  if (!red.attractor_exists()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::linear was called.\n";
    exit(1);
  }
  double w= 0;
  int i,j;
  int nuatr = red.attractor_size();
  if (nuatr==1) {
    int *elatra;
    elatra = new int[red.number_of_nodes()];
    for (i=0; i < red.number_of_nodes(); i++)
      elatra[i] = red.attractor_element(0,i);
    w = 1 - distance(goal, elatra, red.number_of_nodes());
    delete [] elatra;
  }
  else {
    int **elatra;
    elatra = new int*[red.attractor_size()];
    for (i=0; i < red.attractor_size(); i++)
      elatra[i] = new int[red.number_of_nodes()];
    for (i=0; i < red.attractor_size(); i++)
      for (j=0; j < red.number_of_nodes(); j++)
        elatra[i][j] = red.attractor_element(i,j);
    w = 1 - distance(goal, elatra, red.attractor_size(), red.number_of_nodes());
    for (i=0; i < red.attractor_size(); i++)
      delete [] elatra[i];
    delete [] elatra;
  }
  return w;
}

double FitnessI::linear(GraphI &red, int **goal, int lencyc)
{
  if (!red.attractor_exists()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::linear was called.\n";
    exit(1);
  }
  double w= 0;
  int i,j;
  int **elatra;
  elatra = new int*[red.attractor_size()];
  for (i=0; i < red.attractor_size(); i++)
    elatra[i] = new int[red.number_of_nodes()];
  for (i=0; i < red.attractor_size(); i++)
    for (j=0; j < red.number_of_nodes(); j++)
      elatra[i][j] = red.attractor_element(i,j);
  w = 1 - distance(goal, lencyc, elatra, red.attractor_size(), red.number_of_nodes());
  for (i=0; i < red.attractor_size(); i++)
    delete [] elatra[i];
  delete [] elatra;
  return w;
}

double FitnessI::linear_nz(GraphI &red, int *goal)
{
  if (!red.attractor_exists()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::linear was called.\n";
    exit(1);
  }
  double w= 0;
  if (!zeros_in_attractor(red))
    w = linear(red, goal);
  return w;
}

double FitnessI::linear_nz(GraphI &red, int **goal, int lencyc)
{
  if (!red.attractor_exists()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::linear was called.\n";
    exit(1);
  }
  double w= 0;
  if (!zeros_in_attractor(red))
    w = linear(red, goal, lencyc);
  return w;
}
//

double FitnessI::pl(GraphI &red, int *goal, double expo)
{
  if (!red.attractor_exists()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::pl was called.\n";
    exit(1);
  }
  double w= 0;
  int i,j;
  int nuatr = red.attractor_size();
  if (nuatr==1) {
    int *elatra;
    elatra = new int[red.number_of_nodes()];
    for (i=0; i < red.number_of_nodes(); i++)
      elatra[i] = red.attractor_element(0,i);
    w = pow(1 - distance(goal, elatra, red.number_of_nodes()), expo);
    delete [] elatra;
  }
  else {
    int **elatra;
    elatra = new int*[red.attractor_size()];
    for (i=0; i < red.attractor_size(); i++)
      elatra[i] = new int[red.number_of_nodes()];
    for (i=0; i < red.attractor_size(); i++)
      for (j=0; j < red.number_of_nodes(); j++)
        elatra[i][j] = red.attractor_element(i,j);
    w = pow(1 - distance(goal, elatra, red.attractor_size(), red.number_of_nodes()), expo);
    for (i=0; i < red.attractor_size(); i++)
      delete [] elatra[i];
    delete [] elatra;
  }
  return w;
}

double FitnessI::pl(GraphI &red, int **goal, int lencyc, double expo)
{
  if (!red.attractor_exists()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::pl was called.\n";
    exit(1);
  }
  double w= 0;
  int i,j;
  int **elatra;
  elatra = new int*[red.attractor_size()];
  for (i=0; i < red.attractor_size(); i++)
    elatra[i] = new int[red.number_of_nodes()];
  for (i=0; i < red.attractor_size(); i++)
    for (j=0; j < red.number_of_nodes(); j++)
      elatra[i][j] = red.attractor_element(i,j);
  w = pow(1 - distance(goal, lencyc, elatra, red.attractor_size(), red.number_of_nodes()), expo);
  for (i=0; i < red.attractor_size(); i++)
    delete [] elatra[i];
  delete [] elatra;
  return w;
}

double FitnessI::pl_nz(GraphI &red, int *goal, double expo)
{
  if (!red.attractor_exists()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::pl was called.\n";
    exit(1);
  }
  double w= 0;
  if (!zeros_in_attractor(red))
    w = pl(red, goal, expo);
  return w;
}

double FitnessI::pl_nz(GraphI &red, int **goal, int lencyc, double expo)
{
  if (!red.attractor_exists()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::pl was called.\n";
    exit(1);
  }
  double w= 0;
  if (!zeros_in_attractor(red))
    w = pl(red, goal, lencyc, expo);
  return w;
}

double FitnessI::frac1s_strict(GraphI &red, double frad)
{
  if (!red.attractor_exists()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::frac1s_strict was called.\n";
    exit(1);
  }
  if ((frad > 1) || (frad < 0)) {
    cout << "[Error]: Fraction of active genes is out of bounds. FitnessI::frac1s_strict.\n";
    exit(1);
  }
  double w=0, fra;
  fra = get_frac_1s(red);
  if (fra==frad)
    w = 1;
  return w;
}

double FitnessI::frac1s_strict_nz(GraphI &red, double frad)
{
  double w = 0;
  if (!zeros_in_attractor(red))
    w = frac1s_strict(red, frad);
  return w;
}

double FitnessI::frac1s_pl(GraphI &red, double frad, double expo)
{
  if (!red.attractor_exists()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::frac1s_strict was called.\n";
    exit(1);
  }
  if ((frad > 1) || (frad < 0)) {
    cout << "[Error]: Fraction of active genes is out of bounds. FitnessI::frac1s_strict.\n";
    exit(1);
  }
  double w=0, fra,x;
  fra = get_frac_1s(red);
  x = fabs(frad-fra);
  w = pow((1-x), expo);
  return w;
}

double FitnessI::frac1s_pl_nz(GraphI &red, double frad, double expo)
{
  if (!red.attractor_exists()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::frac1s_strict was called.\n";
    exit(1);
  }
  if ((frad > 1) || (frad < 0)) {
    cout << "[Error]: Fraction of active genes is out of bounds. FitnessI::frac1s_strict.\n";
    exit(1);
  }
  double w=0, fra,x;
  if (!zeros_in_attractor(red)) {
    fra = get_frac_1s(red);
    x = fabs(frad-fra);
    w = pow((1-x), expo);
  }
  return w;
}

double FitnessI::strict_pathlength(GraphI &red, int opt)
{
  if (!red.attractor_exists()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::strict_pathlength was called.\n";
    exit(1);
  }
  double w = 0;
  int pal = red.path_length();
  if (pal == opt)
    w = 1;
  return w;
}

double FitnessI::greater_pathlength(GraphI &red, int opt, double penalone)
{
  if (!red.attractor_exists()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::greater_pathlength was called.\n";
    exit(1);
  }
  double w = 0;
  int pal = red.path_length();
  if (pal > opt)
    w = 1;
  else
    w = pow(penalone, double(opt-pal));
  return w;
}

double FitnessI::greater_pathlength(GraphI &red, int opt)
{
  if (!red.attractor_exists()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::greater_pathlength was called.\n";
    exit(1);
  }
  double w = 0;
  int pal = red.path_length();
  if (pal >= opt)
    w = 1;
  return w;
}

double FitnessI::greater_pathlength_nz(GraphI &red, int opt, double penalone)
{
  if (!red.attractor_exists()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::greater_pathlength_nz was called.\n";
    exit(1);
  }
  double w = 0;
  int pal;
  if (!zeros_in_attractor(red)) {
    pal= red.path_length();
    if (pal > opt)
      w = 1;
    else
      w = pow(penalone, double(opt-pal));
  }
  return w;
}

double FitnessI::greater_pathlength_nz(GraphI &red, int opt)
{
  if (!red.attractor_exists()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::greater_pathlength_nz was called.\n";
    exit(1);
  }
  double w = 0;
  int pal;
  if (!zeros_in_attractor(red)) {
    pal= red.path_length();
    if (pal >= opt)
      w = 1;
  }
  return w;
}

double FitnessI::greater_atsize(GraphI &red, int opt, double penalone)
{
  if (!red.attractor_exists()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::greater_atsize was called.\n";
    exit(1);
  }
  double w = 0;
  int pal = red.attractor_size();
  if (pal > opt)
    w = 1;
  else
    w = pow(penalone, double(opt-pal));
  return w;
}

double FitnessI::greater_atsize(GraphI &red, int opt)
{
  if (!red.attractor_exists()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::greater_atsize was called.\n";
    exit(1);
  }
  double w = 0;
  int pal = red.attractor_size();
  if (pal >= opt)
    w = 1;
  return w;
}

double FitnessI::greater_atsize_nz(GraphI &red, int opt, double penalone)
{
  if (!red.attractor_exists()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::greater_atsize_nz was called.\n";
    exit(1);
  }
  double w = 0;
  int pal;
  if (!zeros_in_attractor(red)) {
    pal = red.attractor_size();
    if (pal > opt)
      w = 1;
    else
      w = pow(penalone, double(opt-pal));
  }
  return w;
}

double FitnessI::greater_atsize_nz(GraphI &red, int opt)
{
  if (!red.attractor_exists()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::greater_atsize_nz was called.\n";
    exit(1);
  }
  double w = 0;
  int pal;
  if (!zeros_in_attractor(red)) {
    pal = red.attractor_size();
    if (pal >= opt)
      w = 1;
  }
  return w;
}

//require comb
void FitnessI::linear_mult(GraphI &red, int **goals, int nugoals, double *ws) //only fixed-points as goals
{
  if (!red.attractors_exist()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::linear_mult was called.\n";
    exit(1);
  }
  basic.fillv0(ws, nugoals);
  int i,j,k;
  int *atrtams;
  atrtams = new int[nugoals];
  int **elatra;
  for (i=0; i<nugoals; i++)
    atrtams[i] = red.attractor_size(i);
  for (i=0; i<nugoals; i++) {
    elatra = new int*[atrtams[i]];
    for (j=0; j < atrtams[i]; j++)
      elatra[j] = new int[red.number_of_nodes()];
    for (j=0; j < atrtams[i]; j++)
      for (k=0; k < red.number_of_nodes(); k++)
        elatra[j][k] = red.attractor_element(i,j,k);
    ws[i] = 1 - distance(goals[i], elatra, atrtams[i], red.number_of_nodes());
    for (j=0; j < atrtams[i]; j++)
      delete [] elatra[j];
    delete [] elatra;
  }
  delete [] atrtams;
}

void FitnessI::linear_mult(GraphI &red, int ***goals, int nugoals, int *gosize, double *ws)
{
  if (!red.attractors_exist()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::linear_mult was called.\n";
    exit(1);
  }
  basic.fillv0(ws, nugoals);
  int i,j,k;
  int *atrtams;
  atrtams = new int[nugoals];
  int **elatra;
  for (i=0; i<nugoals; i++)
    atrtams[i] = red.attractor_size(i);
  for (i=0; i<nugoals; i++) {
    elatra = new int*[atrtams[i]];
    for (j=0; j < atrtams[i]; j++)
      elatra[j] = new int[red.number_of_nodes()];
    for (j=0; j < atrtams[i]; j++)
      for (k=0; k < red.number_of_nodes(); k++)
        elatra[j][k] = red.attractor_element(i,j,k);
    ws[i] = 1 - distance(goals[i], gosize[i], elatra, atrtams[i], red.number_of_nodes());
    for (j=0; j < atrtams[i]; j++)
      delete [] elatra[j];
    delete [] elatra;
  }
  delete [] atrtams;
}

//
void FitnessI::linear_mult_nz(GraphI &red, int **goals, int nugoals, double *ws) //only fixed-points as goals
{
  if (!red.attractors_exist()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::linear_mult was called.\n";
    exit(1);
  }
  basic.fillv0(ws, nugoals);
  int i,j,k;
  int *atrtams;
  atrtams = new int[nugoals];
  int **elatra;
  for (i=0; i<nugoals; i++)
    atrtams[i] = red.attractor_size(i);
  for (i=0; i<nugoals; i++) {
    if (!zeros_in_attractor(red, i)) {
      elatra = new int*[atrtams[i]];
      for (j=0; j < atrtams[i]; j++)
        elatra[j] = new int[red.number_of_nodes()];
      for (j=0; j < atrtams[i]; j++)
        for (k=0; k < red.number_of_nodes(); k++)
          elatra[j][k] = red.attractor_element(i,j,k);
      ws[i] = 1 - distance(goals[i], elatra, atrtams[i], red.number_of_nodes());
      for (j=0; j < atrtams[i]; j++)
        delete [] elatra[j];
      delete [] elatra;
    }
    else
      ws[i] = 0;
  }
  delete [] atrtams;
}

void FitnessI::linear_mult_nz(GraphI &red, int ***goals, int nugoals, int *gosize, double *ws)
{
  if (!red.attractors_exist()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::linear_mult was called.\n";
    exit(1);
  }
  basic.fillv0(ws, nugoals);
  int i,j,k;
  int *atrtams;
  atrtams = new int[nugoals];
  int **elatra;
  for (i=0; i<nugoals; i++)
    atrtams[i] = red.attractor_size(i);
  for (i=0; i<nugoals; i++) {
    if (!zeros_in_attractor(red, i)) {
      elatra = new int*[atrtams[i]];
      for (j=0; j < atrtams[i]; j++)
        elatra[j] = new int[red.number_of_nodes()];
      for (j=0; j < atrtams[i]; j++)
        for (k=0; k < red.number_of_nodes(); k++)
          elatra[j][k] = red.attractor_element(i,j,k);
      ws[i] = 1 - distance(goals[i], gosize[i], elatra, atrtams[i], red.number_of_nodes());
      for (j=0; j < atrtams[i]; j++)
        delete [] elatra[j];
      delete [] elatra;
    }
    else
      ws[i] = 0;
  }
  delete [] atrtams;
}

////
void FitnessI::pl_mult(GraphI &red, int **goals, int nugoals, double *ws, double expo) //only fixed-points as goals
{
  if (!red.attractors_exist()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::pl_mult was called.\n";
    exit(1);
  }
  basic.fillv0(ws, nugoals);
  int i,j,k;
  int *atrtams;
  atrtams = new int[nugoals];
  int **elatra;
  for (i=0; i<nugoals; i++)
    atrtams[i] = red.attractor_size(i);
  for (i=0; i<nugoals; i++) {
    elatra = new int*[atrtams[i]];
    for (j=0; j < atrtams[i]; j++)
      elatra[j] = new int[red.number_of_nodes()];
    for (j=0; j < atrtams[i]; j++)
      for (k=0; k < red.number_of_nodes(); k++)
        elatra[j][k] = red.attractor_element(i,j,k);
    ws[i] = pow(1 - distance(goals[i], elatra, atrtams[i], red.number_of_nodes()), expo);
    for (j=0; j < atrtams[i]; j++)
      delete [] elatra[j];
    delete [] elatra;
  }
  delete [] atrtams;
}

void FitnessI::pl_mult(GraphI &red, int ***goals, int nugoals, int *gosize, double *ws, double expo)
{
  if (!red.attractors_exist()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::pl_mult was called.\n";
    exit(1);
  }
  basic.fillv0(ws, nugoals);
  int i,j,k;
  int *atrtams;
  atrtams = new int[nugoals];
  int **elatra;
  for (i=0; i<nugoals; i++)
    atrtams[i] = red.attractor_size(i);
  for (i=0; i<nugoals; i++) {
    elatra = new int*[atrtams[i]];
    for (j=0; j < atrtams[i]; j++)
      elatra[j] = new int[red.number_of_nodes()];
    for (j=0; j < atrtams[i]; j++)
      for (k=0; k < red.number_of_nodes(); k++)
        elatra[j][k] = red.attractor_element(i,j,k);
    ws[i] = pow(1 - distance(goals[i], gosize[i], elatra, atrtams[i], red.number_of_nodes()), expo);
    for (j=0; j < atrtams[i]; j++)
      delete [] elatra[j];
    delete [] elatra;
  }
  delete [] atrtams;
}


void FitnessI::pl_mult_nz(GraphI &red, int **goals, int nugoals, double *ws, double expo) //only fixed-points as goals
{
  if (!red.attractors_exist()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::pl_mult was called.\n";
    exit(1);
  }
  basic.fillv0(ws, nugoals);
  int i,j,k;
  int *atrtams;
  atrtams = new int[nugoals];
  int **elatra;
  for (i=0; i<nugoals; i++)
    atrtams[i] = red.attractor_size(i);
  for (i=0; i<nugoals; i++) {
    if (!zeros_in_attractor(red, i)) {
      elatra = new int*[atrtams[i]];
      for (j=0; j < atrtams[i]; j++)
        elatra[j] = new int[red.number_of_nodes()];
      for (j=0; j < atrtams[i]; j++)
        for (k=0; k < red.number_of_nodes(); k++)
          elatra[j][k] = red.attractor_element(i,j,k);
      ws[i] = pow(1 - distance(goals[i], elatra, atrtams[i], red.number_of_nodes()), expo);
      for (j=0; j < atrtams[i]; j++)
        delete [] elatra[j];
      delete [] elatra;
    }
    else
      ws[i] = 0;
  }
  delete [] atrtams;
}

void FitnessI::pl_mult_nz(GraphI &red, int ***goals, int nugoals, int *gosize, double *ws, double expo)
{
  if (!red.attractors_exist()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::pl_mult was called.\n";
    exit(1);
  }
  basic.fillv0(ws, nugoals);
  int i,j,k;
  int *atrtams;
  atrtams = new int[nugoals];
  int **elatra;
  for (i=0; i<nugoals; i++)
    atrtams[i] = red.attractor_size(i);
  for (i=0; i<nugoals; i++) {
    if (!zeros_in_attractor(red, i)) {
      elatra = new int*[atrtams[i]];
      for (j=0; j < atrtams[i]; j++)
        elatra[j] = new int[red.number_of_nodes()];
      for (j=0; j < atrtams[i]; j++)
        for (k=0; k < red.number_of_nodes(); k++)
          elatra[j][k] = red.attractor_element(i,j,k);
      ws[i] = pow(1 - distance(goals[i], gosize[i], elatra, atrtams[i], red.number_of_nodes()), expo);
      for (j=0; j < atrtams[i]; j++)
        delete [] elatra[j];
      delete [] elatra;
    }
    else
      ws[i] = 0;
  }
  delete [] atrtams;
}

void FitnessI::frac1s_strict_mult(GraphI &red, double fad, double *ws)
{
  if (!red.attractors_exist()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::frac1s_strict_mult was called.\n";
    exit(1);
  }
  basic.fillv0(ws, red.number_of_attractors());
  int i;
  if ((fad > 1) || (fad < 0)) {
    cout << "[Error]: Fraction of active genes is out of bounds. FitnessI::frac1s_strict_mult.\n";
    exit(1);
  }
  double fra;
  for (i=0; i<red.number_of_attractors(); i++) {
    ws[i] = 0;
    fra = get_frac_1s(red, i);
    if (fra == fad)
      ws[i] = 1;
  }
}

void FitnessI::frac1s_strict_mult_nz(GraphI &red, double fad, double *ws)
{
  if (!red.attractors_exist()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::frac1s_strict_mult was called.\n";
    exit(1);
  }
  basic.fillv0(ws, red.number_of_attractors());
  int i;
  if ((fad > 1) || (fad < 0)) {
    cout << "[Error]: Fraction of active genes is out of bounds. FitnessI::frac1s_strict_mult.\n";
    exit(1);
  }
  double fra;
  for (i=0; i<red.number_of_attractors(); i++) {
    ws[i] = 0;
    if (!zeros_in_attractor(red, i)) {
      fra = get_frac_1s(red, i);
      if (fra == fad)
        ws[i] = 1;
    }
  }
}


void FitnessI::frac1s_pl_mult(GraphI &red, double frad, double *ws, double expo)
{
  if (!red.attractors_exist()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::frac1s_strict_mult was called.\n";
    exit(1);
  }
  basic.fillv0(ws, red.number_of_attractors());
  int i;
  if ((frad > 1) || (frad < 0)) {
    cout << "[Error]: Fraction of active genes is out of bounds. FitnessI::frac1s_strict_mult.\n";
    exit(1);
  }
  double fra,x;
  for (i=0; i<red.number_of_attractors(); i++) {
    ws[i] = 0;
    fra = get_frac_1s(red, i);
    x = fabs(frad-fra);
    ws[i] = pow((1-x), expo);
  }
}


void FitnessI::frac1s_pl_mult_nz(GraphI &red, double frad, double *ws, double expo)
{
  if (!red.attractors_exist()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::frac1s_strict_mult was called.\n";
    exit(1);
  }
  basic.fillv0(ws, red.number_of_attractors());
  int i;
  if ((frad > 1) || (frad < 0)) {
    cout << "[Error]: Fraction of active genes is out of bounds. FitnessI::frac1s_strict_mult.\n";
    exit(1);
  }
  double fra,x;
  for (i=0; i<red.number_of_attractors(); i++) {
    ws[i] = 0;
    if (!zeros_in_attractor(red, i)) {
      fra = get_frac_1s(red, i);
      x = fabs(frad-fra);
      ws[i] = pow((1-x), expo);
    }
  }
}

///

double FitnessI::distance(int *goal, int *atr, int tam)
{
  int dif = 0,i;
  for (i=0; i<tam; i++)
    if (goal[i] != atr[i])
      dif++;
  return double(dif)/double(tam);
}

double FitnessI::distance(int *goal, int **atr, int numatr, int tam)
{
  int i,j;
  double dif = 0;
  for (i=0; i<numatr; i++)
    for (j=0; j < tam; j++)
      if (goal[j] != atr[i][j])
        dif = dif + (1.0/double(numatr));
  return dif/double(tam);
}

double FitnessI::distance(int **goal, int nugs, int **atr, int numatr, int tam)
{
  int i,j,k;
  double dis;
  if (nugs != numatr) {
    int mulng, mulna;
    if (nugs < numatr) {
      if ((numatr%nugs)==0) {
        mulna = 1;
        mulng = numatr/nugs;
      }
      else {
        mulng = numatr;
        mulna = nugs;
      }
    }
    else {
      if ((nugs%numatr)==0) {
        mulng = 1;
        mulna = nugs/numatr;
      }
      else {
        mulng = numatr;
        mulna = nugs;
      }
    }
    int **nugoal, **nuatr;
    nugoal = new int*[nugs*mulng];
    for (i= 0; i<(nugs*mulng); i++)
      nugoal[i] = new int[tam];
    nuatr = new int*[numatr*mulna];
    for (i=0; i <(numatr*mulna); i++)
      nuatr[i] = new int[tam];
    for (i=0; i<mulng; i++)
      for (j=0; j<nugs; j++)
        for (k=0; k < tam; k++)
          nugoal[(nugs*i)+j][k] = goal[j][k];
    for (i=0; i < mulna; i++)
      for (j=0; j<numatr; j++)
        for (k=0; k<tam;k++)
          nuatr[(numatr*i)+j][k] = atr[j][k];
    dis = distance_aux(nugoal, (nugs*mulng), nuatr, (numatr*mulna), tam);
    for (i= 0; i<(nugs*mulng); i++)
      delete [] nugoal[i];
    delete [] nugoal;
    for (i=0; i <(numatr*mulna); i++)
      delete [] nuatr[i];
    delete [] nuatr;
  }
  else
    dis = distance_aux(goal, nugs, atr, numatr, tam);
  return dis;
}

double FitnessI::comb_mult(double *ws, int tam)
{
  double w = basic.multinvec(ws, tam);
  return w;
}

double FitnessI::comb_mult(double *ws, int tam, double *weights)
{
  double w = 1;
  int i;
  for (i=0; i < tam; i++)
    w = w*(pow(ws[i], weights[i]));
  return w;
}

double FitnessI::comb_add(double *ws, int tam)
{
  double w = basic.sumatoria(ws, tam);
  w = w/double(tam);
  return w;
}

double FitnessI::comb_add(double *ws, int tam, double *weights)
{
  double w = basic.dotproduct(ws, weights, tam);
  w = w/double(tam);
  return w;
}

double FitnessI::mult_multi_goals_under_pert(GraphI& red, int numofgoals, int numcells, int **ci, int ***goal, int *periodgoal, double pertrate, double selcoefperg, double multexp) {
  double fi = 1;
  int i;
  for (i = 0; i < numofgoals; i++)
    fi *= under_pert(red, numcells, ci[i], goal[i], periodgoal[i], pertrate, selcoefperg, multexp);
  return fi;
}

double FitnessI::additive_multi_goals_under_pert(GraphI& red, int numofgoals, int numcells, int **ci, int ***goal, int *periodgoal, double pertrate, double selcoefperg, double multexp) {
  double fi = 0;
  int i;
  for (i = 0; i < numofgoals; i++)
    fi += under_pert(red, numcells, ci[i], goal[i], periodgoal[i], pertrate, selcoefperg, multexp);
  return fi;
}

double FitnessI::under_pert(GraphI& red, int numcells, int *ci, int **goal, int periodgoal, double pertrate, double selcoefperg, double multexp) {
  int tam = red.number_of_nodes();
  int i, j, k;
  int **cifrec;
  basic.create_array(cifrec, tam+1, tam);
  for (i = 0; i <= tam; i++)
    for (j = 0; j < tam; j++)
      cifrec[i][j] = ci[j];
  for (i = 1; i <= tam; i++)
    cifrec[i][i-1] *= -1;
  
  red.find_attractors(cifrec, tam+1);
  int numat = red.number_of_attractors();
  double *fifrec;
  fifrec = new double[numat];
  basic.fillv0(fifrec, numat);
  int *sizofat, ***losat;
  sizofat = new int[numat];
  for (i = 0; i < numat; i++)
    sizofat[i] = red.attractor_size(i);
  losat = new int**[numat];
  for (i = 0; i < numat; i++) {
    losat[i] = new int*[sizofat[i]];
    for (j = 0; j < sizofat[i]; j++)
      losat[i][j] = new int[tam];
  }
  for (i = 0; i < numat; i++)
    for (j = 0; j < sizofat[i]; j++)
      for (k=0; k < tam; k++)
        losat[i][j][k] = red.attractor_element(i, j, k);
  red.clear_attractors();
  double d;
  for (i = 0; i < numat; i++) {
    d = distance(goal, periodgoal, losat[i], sizofat[i], tam);
    d *= (multexp*tam);
    fifrec[i] = pow((1-selcoefperg), d);
  }
  for (i = 0; i < numat; i++) {
    for (j = 0; j < sizofat[i]; j++)
      delete [] losat[i][j];
    delete [] losat[i];
  }
  delete [] losat;
  delete [] sizofat;
  
  int *cicact, nudi;
  double sumafi = 0;
  cicact = new int[tam];
  int **unat;
  for (i = 0; i < numcells; i++) {
    nudi = 0;
    for (j = 0; j < tam; j++) {
      if (est.randreal() < pertrate) {
        cicact[j] = ci[j]*(-1);
        nudi++;
      }
      else
        cicact[j] = ci[j];
    }
    if (nudi == 0) {
      sumafi += fifrec[0];
    }
    if (nudi == 1) {
      for (j = 0; j < numat; j++)
        if (basic.eqvec(cicact, tam, cifrec[j], tam)) {
          sumafi += fifrec[j];
          break;
        }
    }
    if (nudi > 1) {
      red.set_as_state(cicact);
      red.find_an_attractor();
      basic.create_array(unat, red.attractor_size(), tam);
      for (j = 0; j < red.attractor_size(); j++)
        for (k = 0; k < tam; k++)
          unat[j][k] = red.attractor_element(j, k);
      d = distance(goal, periodgoal, unat, red.attractor_size(), tam);
      d *= (multexp*tam);
      sumafi += pow((1-selcoefperg), d);
      for (j = 0; j < red.attractor_size(); j++)
        delete [] unat[j];
      delete [] unat;
      red.clear_attractor();
    }
  }
  for (i = 0; i <= tam; i++)
    delete [] cifrec[i];
  delete [] cifrec;
  delete [] cicact;
  delete [] fifrec;
    sumafi /= double(numcells);
    return sumafi;
}



//private
double FitnessI::distance_aux(int **goal, int nugs, int **atr, int numatr, int tam)
{
  if (nugs != numatr) {
    cout << "[Error]: Wrong number of matrix rows in FitnessI::distance_aux.\n";
    exit(1);
  }
  int i,j,k;
  double *dists,d;
  dists = new double[nugs];
  basic.fillv0(dists, nugs);
  for (i=0; i<nugs; i++) {
    dists[i] =0;
    for (j=0; j<nugs; j++) {
      for (k=0;k<tam;k++)
        if (goal[j][k] != atr[(j+i)%nugs][k])
          dists[i] = dists[i] + (1.0/double(nugs));
    }
  }
  d = basic.find_min(dists, nugs);
  d = d/double(tam);
  delete [] dists;
  return d;
}

bool FitnessI::zeros_in_attractor(GraphI &red)
{
  if (!red.attractor_exists()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::zeros_in_attractor was called.\n";
    exit(1);
  }
  bool res = false;
  int i,j;
  for (i=0; i<red.attractor_size(); i++)
    for (j=0; j < red.number_of_nodes(); j++)
      if (red.attractor_element(i,j)==0) {
        res = true;
        break;
      }
  return res;
}

bool FitnessI::zeros_in_attractor(GraphI &red, int wh)
{
  if (!red.attractors_exist()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::zeros_in_attractor was called.\n";
    exit(1);
  }
  bool res = false;
  int i,j;
  for (i=0; i<red.attractor_size(wh); i++)
    for (j=0; j < red.number_of_nodes(); j++)
      if (red.attractor_element(wh,i,j)==0) {
        res = true;
        break;
      }
  return res;
}

double FitnessI::get_frac_1s(GraphI &red)
{
  if (!red.attractor_exists()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::get_frac_1s was called.\n";
    exit(1);
  }
  int nuu = 0, i,j;
  for (i=0; i < red.attractor_size(); i++)
    for (j=0; j<red.number_of_nodes(); j++)
      if (red.attractor_element(i,j) == 1)
        nuu++;
  double res = double(nuu)/double(red.attractor_size()*red.number_of_nodes());
  return res;
}

double FitnessI::get_frac_1s(GraphI &red, int wh)
{
  if (!red.attractors_exist()) {
    cout << "[Error]: Attractors had not been searched when FitnessI::get_frac_1s was called.\n";
    exit(1);
  }
  int nuu = 0, i,j;
  for (i=0; i < red.attractor_size(wh); i++)
    for (j=0; j<red.number_of_nodes(); j++)
      if (red.attractor_element(wh,i,j) == 1)
        nuu++;
  double res = double(nuu)/double(red.attractor_size(wh)*red.number_of_nodes());
  return res;
}
