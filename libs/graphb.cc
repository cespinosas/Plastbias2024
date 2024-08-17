#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <set>
#include <list>
#include "basics.h"
#include "graphi.h"
#include "graphb.h"
#include "graphc.h"
#include <string>

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

GraphB::GraphB()
{
  ionary = false;
  yamotdab=false;
  yasgcat=false;
  yaprecosg=false;
}

GraphB::GraphB(Alea& jacta)
{
  start_rng(jacta);
  ionary = false;
  yamotdab=false;
  yasgcat=false;
  yaprecosg=false;
}

void GraphB::start_rng(Alea& jacta)
{
  est = jacta;
  basic.start_rng(est);
}

GraphB::GraphB(int n, Alea& jacta, bool dirnw)
{
  start_rng(jacta);
  make_nw(n, dirnw);
}

void GraphB::make_nw(int n, bool dirnw)
{
  int i;
  size = n;
  nw = new bool*[size];
  for (i = 0; i < size; i++)
    nw[i] = new bool[size];
  basic.fillmat0(nw, size, size);
  ionary = false;
  directed = dirnw;
  set_default_exclusive_vars();
}

void GraphB::copy(GraphB &templ)
{
  int n = templ.number_of_nodes();
  bool dirnw = templ.is_directed();
  make_nw(n, dirnw);
  int i, j;
  for (i = 0; i < size; i++)
    for (j = 0; j < size; j++)
      nw[i][j] = templ.weight(j, i);
}

void GraphB::copy(GraphB &templ, bool withnames)//lo que sigue solo en GraphB
{
  copy(templ);
  int i;
  string *mientras;
  if (withnames && templ.names_exist()) {
    mientras = new string[size];
    for (i=0; i<size; i++)
      mientras[i] = templ.get_name(i);
    assign_names(mientras, size);
    delete [] mientras;
  }
}

void GraphB::make_copy(GraphB& vacia)
{
  vacia.make_nw(size, is_directed());
  int i,j;
  if (is_directed()) {
    for (i=0; i<size; i++)
      for (j=0; j<size; j++)
        vacia.force_interaction(j, i, weight(j, i));
  }
  else {
    for (i=0; i<size; i++)
      for (j=0; j<=i; j++)
        vacia.force_interaction_undir(j, i, weight(j, i));
  }
}

void GraphB::put_subgraph_into(GraphB& nueva, set<int> &cuales)
{
  set<int> cua2 = cuales;
  nueva.make_nw(cuales.size(), is_directed());
  set<int>::iterator it, it2;
  int i,j;
  i=0;
  for (it=cuales.begin(); it!=cuales.end(); it++) {
    j=0;
    for (it2=cua2.begin(); it2!=cua2.end(); it2++) {
      if (directed)
        nueva.force_interaction(j, i, weight(*it2, *it));
      else
        nueva.force_interaction_undir(j, i, weight(*it2, *it));
      j++;
    }
    i++;
  }
}

void GraphB::put_subgraph_into(GraphB& nueva, set<int> &cuales, bool withnames)
{
  put_subgraph_into(nueva, cuales);
  string *mientras;
  if (withnames && names_exist()) {
    mientras = new string[cuales.size()];
    int i=0;
    set<int>::iterator it;
    for (it=cuales.begin(); it!=cuales.end(); it++) {
      mientras[i] = get_name(*it);
      i++;
    }
    nueva.assign_names(mientras, nueva.number_of_nodes());
    delete [] mientras;
  }
}

void GraphB::make_subgraph_of(GraphB &templ, set<int> &cuales)
{
  set<int> cua2 = cuales;
  make_nw(cuales.size(), templ.is_directed());
  set<int>::iterator it, it2;
  int i,j;
  i=0;
  for (it=cuales.begin(); it!=cuales.end(); it++) {
    j=0;
    for (it2=cua2.begin(); it2!=cua2.end(); it2++) {
      nw[i][j] = templ.weight(*it2, *it);
      j++;
    }
    i++;
  }
}

void GraphB::make_subgraph_of(GraphB &templ, set<int> &cuales, bool withnames)
{
  make_subgraph_of(templ, cuales);
  string *mientras;
  if (withnames && templ.names_exist()) {
    mientras = new string[cuales.size()];
    int i=0;
    set<int>::iterator it;
    for (it=cuales.begin(); it!=cuales.end(); it++) {
      mientras[i] = templ.get_name(*it);
      i++;
    }
    assign_names(mientras, number_of_nodes());
    delete [] mientras;
  }
}

void GraphB::copy_wdeg(GraphB &templ)
{
  copy(templ);
  numofe = templ.number_of_edges();
  prepare_for_degrees();
  int i;
  for (i = 0; i < size; i++) {
    indegree[i] = get_indegree(i);
    outdegree[i] = get_outdegree(i);
    degree[i] = get_degree(i);
  }
  yadeg = true;
}

void GraphB::copy_wdeg(GraphB &templ, bool withnames)
{
  copy_wdeg(templ);
  int i;
  string *mientras;
  if (withnames && templ.names_exist()) {
    mientras = new string[size];
    for (i=0; i<size; i++)
      mientras[i] = templ.get_name(i);
    assign_names(mientras, size);
    delete [] mientras;
  }
}

GraphB::GraphB(int n, double c, Alea& jacta, bool dirnw)
{
  start_rng(jacta);
  make_nw(n, c, dirnw);
}

void GraphB::make_nw(int n, double c, bool dirnw)
{
  if (c >= 1) {
    cout << "[Error]: You are trying to construct an overconnected network. This error arose when GraphB::make_nw was called.\n";
    exit(1);
  }
  int i, j;
  double nure;
  make_nw(n, dirnw);
  if (directed) {
    for (i = 0; i < size; i++)
      for(j = 0; j < size; j++) {
        nure = est.randreal();
        if (nure <= c)
          nw[i][j] = true;
      }
  }
  else {
    for (i=0; i<size; i++)
      for (j=0; j<=i; j++) {
        nure = est.randreal();
        if (nure <=c) {
          nw[i][j] = true;
          nw[j][i] = true;
        }
      }
  }
}

GraphB::GraphB(int n, int e, Alea& jacta, bool dirnw)
{
  start_rng(jacta);
  make_nw(n, e, dirnw);
}

void GraphB::make_nw(int n, int e, bool dirnw)
{
  if (e > (n*n)) {
    cout << "[Error]: You are trying to construct an overconnected network. This error arose when GraphB::make_nw was called.\n";
    exit(1);
  }
  int i, j, k;
  make_nw(n, dirnw);
  for (i = 0; i < e; i++) {
    do {
      j = est.randint(0, size);
      k = est.randint(0, size);
    } while (nw[j][k]);
    nw[j][k] = true;
    if (!directed)
      nw[k][j] = true;
  }
}

void GraphB::clear()
{
  int i;
  for (i = 0; i < size; i++)
    delete [] nw[i];
  delete [] nw;
  
  if (ionary)
    clear_dict();
  if (infl)
    clear_ds();
  if (dist)
    clear_dima();
  if (yacomp)
    clear_comp();
  if (cfc)
    clear_scc();
  if (bcya)
    clear_bc();
  if (floops)
    clear_loops();
  if (yadeg)
    clear_deg();
  if (yadegdist)
    clear_degdist();
  if (yae)
    clear_edge_count();
  if (yamotdab)
    clear_db_for_isomorphism();
  if (yasgcat)
    clear_sg_catalogue();
  if (yaprecosg)
    clear_preps_to_count_sg();
  if (yamadya)
    clear_adj();
}

void GraphB::clear_preps_to_count_sg()
{
  delete [] vedosg2;
  delete [] vedosg3;
  delete [] vedosg4;
  yaprecosg = false;
}


void GraphB::clear_sg_catalogue()
{
  int i;
  //2
  for (i=0; i<nupre_2; i++) {
    LM_2[i].clear();
    delete [] nompremot_2[i];
  }
  delete [] LM_2;
  delete [] nompremot_2;
  delete [] socont_2;
  
  //3
  for (i=0; i<nupre_3; i++) {
    LM_3[i].clear();
    delete [] nompremot_3[i];
  }
  delete [] LM_3;
  delete [] nompremot_3;
  delete [] socont_3;
  
  //4
  for (i=0; i<nupre_4; i++) {
    LM_4[i].clear();
    delete [] nompremot_4[i];
  }
  delete [] LM_4;
  delete [] nompremot_4;
  delete [] socont_4;
  
  yasgcat = false;
}

void GraphB::clear_db_for_isomorphism()
{
  int i;
  for (i=0; i<2; i++)
    delete [] comb2[i];
  delete [] comb2;
  for (i=0; i<6; i++)
    delete [] comb3[i];
  delete [] comb3;
  for (i=0; i<24; i++)
    delete [] comb4[i];
  delete [] comb4;
  yamotdab = false;
}

void GraphB::clear_edge_count()
{
  numofe = 0;
  yae = false;
}

void GraphB::clear_loops()
{
  delete [] loops;
  totalloops = 0;
  int i;
  for (i = 0; i < size; i++)
    delete [] guada[i];
  delete [] guada;
  floops = false;
}

void GraphB::clear_comp()
{
  delete [] components;
  yacomp = false;
  numbcomp = 0;
}

void GraphB::clear_bc()
{
  delete [] bc;
  bcya = false;
}

void GraphB::clear_scc()
{
  delete [] scc;
  cfc = false;
  numbscc = 0;
}

void GraphB::clear_adj() {
  int i;
  for (i = 0; i < size; i++)
    delete [] matadya[i];
  delete [] matadya;
  yamadya = false;
}

void GraphB::clear_dima()
{
  int i;
  for (i = 0; i < size; i++)
    delete [] dima[i];
  delete [] dima;
  dist = false;
}

void GraphB::clear_ds()
{
  int i;
  for (i = 0; i < size; i++)
    downstream[i].clear();
  delete [] downstream;
  infl = false;
}

void GraphB::clear_dict()
{
  delete [] dict;
  ionary = false;
}

void GraphB::clear_deg()
{
  delete [] degree;
  delete [] outdegree;
  delete [] indegree;
  yadeg = false;
  if (yadegdist)
    clear_degdist();
}

void GraphB::clear_degdist()
{
  delete [] degdist;
  delete [] odegdist;
  delete [] idegdist;
  yadegdist = false;
}

void GraphB::get_all_degrees()
{
  if (yadeg) {
    cout << "[Error]: Attempting to create degree arrays again. This error arose when GraphB::get_all_degrees was called.\n";
    exit(1);
  }
  int i;
  prepare_for_degrees();
  for (i = 0; i < size; i++) {
    indegree[i] = calc_indegree(i);
    outdegree[i] = calc_outdegree(i);
    degree[i] = calc_degree(i);
  }
  yadeg = true;
}

int GraphB::get_indegree(int no)
{
  if (!yadeg)
    get_all_degrees();
  return indegree[no];
}

int GraphB::get_outdegree(int no)
{
  if (!yadeg)
    get_all_degrees();
  return outdegree[no];
}

int GraphB::get_degree(int no)
{
  if (!yadeg)
    get_all_degrees();
  return degree[no];
}

void GraphB::get_degdist()
{
  if (!yadegdist) {
    prepare_for_degdist();
    priv_get_degdist();
    yadegdist = true;
  }
}

void GraphB::get_degdist(int *empvec, int tam)
{
  if (tam != (size+1)){
    cout << "[Error]: Wrong size of array in GraphB::get_degdist.\n";
    exit(1);
  }
  int i;
  get_degdist();
  for (i=0; i<tam; i++)
    empvec[i] = degdist[i];
}

void GraphB::get_odegdist(int *empvec, int tam)
{
  if (tam != (size+1)){
    cout << "[Error]: Wrong size of array in GraphB::get_degdist.\n";
    exit(1);
  }
  int i;
  get_degdist();
  for (i=0; i<tam; i++)
    empvec[i] = odegdist[i];
}

void GraphB::get_idegdist(int *empvec, int tam)
{
  if (tam != (size+1)){
    cout << "[Error]: Wrong size of array in GraphB::get_degdist.\n";
    exit(1);
  }
  int i;
  get_degdist();
  for (i=0; i<tam; i++)
    empvec[i] = idegdist[i];
}


//access
int GraphB::number_of_edges()
{
  if (!yae) {
    numofe = 0;
    int i, j;
    if (directed) {
      for (i = 0; i < size; i++)
        for (j = 0; j < size; j++)
          if (nw[i][j])
            numofe++;
    }
    else {
      for (i = 0; i < size; i++)
        for (j = 0; j <= i; j++)
          if (nw[i][j])
            numofe++;
    }
    yae = true;
  }
  return numofe;
}

int GraphB::number_of_nodes()
{
  return size;
}

bool GraphB::weight(int source, int target)
{
  return nw[target][source];
}

bool GraphB::names_exist()
{
  return ionary;
}

string GraphB::get_name(int node)
{
  return dict[node];
}

bool GraphB::is_directed()
{
  return directed;
}

bool GraphB::already_deg()
{
  return yadeg;
}

int GraphB::name_to_int(string cual)
{
  int res=-1, i;
  if (!names_exist()) {
    cout << "[Error]: Names do not exist when GraphB::name_to_int was called.\n";
    exit(1);
  }
  for (i=0; i<size; i++)
    if (dict[i] == cual) {
      res = i;
      break;
    }
  return res;
}

int GraphB::name_to_int(string cual, bool checkornot)
{
  int res = name_to_int(cual);
  if ((checkornot) && (res==(-1))) {
    cout << "[Error]: There is no node with the name \'" << cual << "\'. This error arose when GraphB::name_to_int was called.\n";
    exit(1);
  }
  return res;
}

bool GraphB::is_name(string elna)
{
  bool res = false;
  int i;
  if (!names_exist()) {
    cout << "[Error]: Names do not exist when GraphB::is_name was called.\n";
    exit(1);
  }
  for (i=0; i<size; i++) {
    if (dict[i] == elna) {
      res = true;
      break;
    }
  }
  return res;
}

//modification
void GraphB::force_interaction(int source, int target, bool value)
{
  if (directed) {
    if (yae) {
      if ((!value) && nw[target][source])
        numofe--;
      else
        if (value && (!nw[target][source]))
          numofe++;
    }
    if (yadeg) {
      if ((!value) && nw[target][source]) {
        outdegree[source]--;
        indegree[target]--;
        degree[source]--;
        degree[target]--;
      }
      else {
        if (value && (!nw[target][source])) {
          outdegree[source]++;
          indegree[target]++;
          degree[source]++;
          degree[target]++;
        }
      }
    }
    nw[target][source] = value;
  }
  else {
    cout << "[Error]: GraphB::force_interaction does not work for undirected graphs\n";
    exit(1);
  }
  if (infl)
    clear_ds();
  if (dist)
    clear_dima();
  if (yacomp)
    clear_comp();
  if (cfc)
    clear_scc();
  if (bcya)
    clear_bc();
  if (floops)
    clear_loops();
  if (yamadya)
    clear_adj();
  if (yadegdist)
    clear_degdist();
}

void GraphB::force_interaction_undir(int source, int target, bool value)
{
  if (directed) {
    cout << "[Error]: GraphB::force_interaction_undir does not work for directed graphs\n";
    exit(1);
  }
  if (yae) {
    if ((!value) && nw[target][source])
      numofe--;
    else
      if (value && (!nw[target][source]))
        numofe++;
  }
  if (yadeg) {
    if ((!value) && nw[target][source]) {
      outdegree[source]--;
      outdegree[target]--;
      indegree[target]--;
      indegree[source]--;
      degree[source]--;
      degree[target]--;
    }
    else
      if (value && (!nw[target][source])) {
        outdegree[source]++;
        outdegree[target]++;
        indegree[target]++;
        indegree[source]++;
        degree[source]++;
        degree[target]++;
      }
  }
  nw[target][source] = value;
  nw[source][target] = value;
  if (infl)
    clear_ds();
  if (dist)
    clear_dima();
  if (yacomp)
    clear_comp();
  if (cfc)
    clear_scc();
  if (bcya)
    clear_bc();
  if (floops)
    clear_loops();
  if (yamadya)
    clear_adj();
  if (yadegdist)
    clear_degdist();
}

void GraphB::change_interaction(int source, int target, bool value)
{
  if (directed) {
    if (value == nw[target][source]) {
      cout << "[Error]: GraphB::change_interaction changes to stay the same.\n";
      exit(1);
    }
    if (yae) {
      if ((!value) && nw[target][source])
        numofe--;
      else
        if (value && (!nw[target][source]))
          numofe++;
    }
    if (yadeg) {
      if ((!value) && nw[target][source]) {
        outdegree[source]--;
        indegree[target]--;
        degree[source]--;
        degree[target]--;
      }
      else {
        if (value && (!nw[target][source])) {
          outdegree[source]++;
          indegree[target]++;
          degree[source]++;
          degree[target]++;
        }
      }
    }
    nw[target][source] = value;
  }
  else {
    cout << "[Error]: GraphB::change_interaction does not work for undirected graphs\n";
    exit(1);
  }
  if (infl)
    clear_ds();
  if (dist)
    clear_dima();
  if (yacomp)
    clear_comp();
  if (cfc)
    clear_scc();
  if (bcya)
    clear_bc();
  if (floops)
    clear_loops();
  if (yamadya)
    clear_adj();
  if (yadegdist)
    clear_degdist();
}

void GraphB::change_interaction_undir(int source, int target, bool value)
{
  if (directed) {
    cout << "[Error]: GraphB::change_interaction_undir does not work for directed graphs\n";
    exit(1);
  }
  if (value == nw[target][source]) {
    cout << "[Error]: GraphB::change_interaction changes to stay the same.";
    exit(1);
  }
  if (yae) {
    if ((!value) && nw[target][source])
      numofe--;
    else
      if (value && (!nw[target][source]))
        numofe++;
  }
  if (yadeg) {
    if ((!value) && nw[target][source]) {
      outdegree[source]--;
      outdegree[target]--;
      indegree[target]--;
      indegree[source]--;
      degree[source]--;
      degree[target]--;
    }
    else
      if (value && (!nw[target][source])) {
        outdegree[source]++;
        outdegree[target]++;
        indegree[target]++;
        indegree[source]++;
        degree[source]++;
        degree[target]++;
      }
  }
  nw[target][source] = value;
  nw[source][target] = value;
  if (infl)
    clear_ds();
  if (dist)
    clear_dima();
  if (yacomp)
    clear_comp();
  if (cfc)
    clear_scc();
  if (bcya)
    clear_bc();
  if (floops)
    clear_loops();
  if (yamadya)
    clear_adj();
  if (yadegdist)
    clear_degdist();
}

void GraphB::change_name(int node, string nn)
{
  if (!(ionary)) {
    cout << "[Error]: There are no names linked to this network. This network arose when GraphB::change_name was called.\n";
    exit(1);
  }
  dict[node] = nn;
}

void GraphB::set_undirected()
{
  if (!(basic.is_symmetric(nw, size, size)))  {
    cout << "[Error]: Attempt to set a non-symmetric network as undirected. This error arose when GraphB::set_undirected was called.\n";
    exit(1);
  }
  else
    directed = false;
}

void GraphB::set_directed()
{
  directed = true;
}

void GraphB::force_undirected()
{
  int i,j;
  for (i=0; i<size; i++) {
    for (j=0; j< size; j++)
      if (weight(i,j) != weight(j,i)) {
        if (weight(i,j) && (!weight(j,i))) {
          if (directed)
            change_interaction(j,i, weight(i,j));
          else
            change_interaction_undir(j,i, weight(i,j));
        }
      }
  }
  if (yae)
    clear_edge_count();
  if (yadeg)
    clear_deg();
  set_undirected();
}

void GraphB::transform_to_isomorph(GraphB& vacia, int* vec)
{
  vacia.make_nw(size, is_directed());
  int i,j;
  if (is_directed()) {
    for (i=0; i<size; i++)
      for (j=0; j<size; j++)
        vacia.force_interaction(vec[j], vec[i], weight(j, i));
  }
  else {
    for (i=0; i<size; i++)
      for (j=0; j<=i; j++)
        vacia.force_interaction_undir(vec[j], vec[i], weight(j, i));
  }
}

//formatting
void GraphB::get_names_from_file(string arch)
{
  ifstream sal;
  basic.open_ifstream(sal, arch);
  get_names_from_file(sal);
  sal.close();
}

void GraphB::get_names_from_file(istream& en)
{
  int i=0;
  if (ionary) {
    cout << "[Error]: Names have already been assigned. This error arose when GraphB::get_names_from_file was called.\n";
    exit(1);
  }
  dict = new string[size];
  while (en >> dict[i])
    i++;
  while (i<size) {
    dict[i] = "";
    i++;
  }
  ionary = true;
}

void GraphB::assign_names(string *nombres, int ent)
{
  int i;
  if (ionary) {
    cout << "[Error]: Names have already been assigned. This error arose when GraphB::assign_names was called.\n";
    exit(1);
  }
  dict = new string[size];
  if (ent > size) {
    cout << "[Error]: More names than nodes. This error arose when GraphB::assign_names was called.";
    exit(1);
  }
  for (i=0; i<ent; i++)
    dict[i] = nombres[i];
  for (i=ent; i<size; i++) {
    dict[i]= "";
  }
  ionary=true;
}

void GraphB::printnw(ostream& sal)
{
  int i, j;
  for (i = 0; i < size; i++) {
    for (j = 0; j < size; j++) {
      if (!nw[i][j])
        sal << "0 ";
      else
        if (nw[i][j])
          sal << "+ ";
    }
    sal << endl;
  }
}

void GraphB::print_dot(string archi)
{
  int i, j;
  ofstream sal;
  string arch;
  arch = archi+".dot";
  basic.open_ofstream(sal, arch);
  string instr, ke;
  if (directed) {
    instr = "digraph G {\n";
    ke = "->";
  }
  else {
    instr = "graph G {\n";
    ke = "--";
  }
  sal << instr;
  if (directed) {
    for (i=0; i < size; i++)
      for (j=0; j < size; j++)
        if (weight(j,i))
          sal << "\t" << j << "  " << ke << "  " << i << ";\n";
  }
  else
    for (i = 0; i < size; i++)
      for (j = 0; j <= i; j++)
        if (weight(j,i))
          sal << "\t" << j << "  " << ke << "  " << i << ";\n";
  sal << "}\n";
  sal.close();
}

void GraphB::print_dot_wn(string archi)
{
  if (!(ionary)) {
    cout << "[Error]: There are no names linked to this network. This error arose when GraphB::print_dot_wn was called.\n";
    exit(1);
  }
  int i, j;
  ofstream sal;
  string arch;
  arch = archi+".dot";
  basic.open_ofstream(sal, arch);
  string instr, ke;
  if (directed) {
    instr = "digraph G {\n";
    ke = "->";
  }
  else {
    instr = "graph G {\n";
    ke = "--";
  }
  sal << instr;
  if (directed) {
    for (i=0; i < size; i++)
      for (j=0; j < size; j++)
        if (weight(j,i))
          sal << "\t\"" << dict[j] << "\"  " << ke << "  \"" << dict[i] << "\";\n";
  }
  else
    for (i = 0; i < size; i++)
      for (j = 0; j <= i; j++)
        if (weight(j,i))
          sal << "\t\"" << dict[j] << "\"  " << ke << "  \"" << dict[i] << "\";\n";
  sal << "}\n";
  sal.close();
}

void GraphB::print_dot_wo_labels_circ(string archi, double radius)
{
  int i, j;
  ofstream sal;
  string arch;
  arch = archi+".dot";
  basic.open_ofstream(sal, arch);
  string instr, ke,brt,frt;
  if (directed) {
    instr = "digraph G {\n";
    ke = "->";
  }
  else {
    instr = "graph G {\n";
    ke = "--";
  }
  sal << instr;
  sal << "layout = neato;\n";
  sal << "node[shape=circle,label=\"\"];\n";
  double radians[size], x[size], y[size], stera, PI=3.141592, van;
  stera = (2*PI)/(size*1.0);
  radians[0] = 0.75*PI;
  van = radians[0]; //
  if (size != 4) {
    for (i=1; i<size; i++) {
      radians[i] = radians[i-1]+stera;
      if (radians[i] >= (2*PI))
        radians[i] = radians[i] - (2*PI);
    }
  }
  else {
    radians[3] = van + stera;
    van = radians[3];
    radians[1] = van +stera;
    van = radians[1];
    radians[2] = van+stera;
    van = radians[2];
    for (i=0; i<size; i++) {
      if (radians[i] >= (2*PI))
        radians[i] = radians[i] - (2*PI);
    }
  }
  
  for (i=0; i<size; i++) {
    basic.polar_to_cartesian(radians[i], radius, x[i], y[i]);
    sal << i << "[pos = \"" << x[i] << "," << y[i] << "!\"];\n";
  }
  
  if (directed) {
    for (i=0; i < size; i++)
      for (j=0; j < size; j++) {
        if (weight(j,i)) {
          if (i==j) {
            if ((radians[i] >= 0) && (radians[i] <= (PI/2.0))) {
              brt="n";
              frt="e";
            }
            if ((radians[i] > (PI/2.0)) && (radians[i] <= PI)) {
              brt="w";
              frt="n";
            }
            if ((radians[i] > PI) && (radians[i] <= (PI*1.5))) {
              brt="s";
              frt="w";
            }
            if ((radians[i] > (PI*1.5)) && (radians[i] < (2*PI))) {
              brt="e";
              frt="s";
            }
            sal << "\t" << j << "  " << ke << "  " << i << " [headport=" << brt << ",tailport=" << frt << "];\n";
          }
          else
            sal << "\t" << j << "  " << ke << "  " << i << ";\n";
        }
      }
  }
  else {
    for (i = 0; i < size; i++)
      for (j = 0; j <= i; j++) {
        if (weight(j,i)) {
          if (i==j) {
            if ((radians[i] >= 0) && (radians[i] <= (PI/2.0))) {
              brt="n";
              frt="e";
            }
            if ((radians[i] > (PI/2.0)) && (radians[i] <= PI)) {
              brt="w";
              frt="n";
            }
            if ((radians[i] > PI) && (radians[i] <= (PI*1.5))) {
              brt="s";
              frt="w";
            }
            if ((radians[i] > (PI*1.5)) && (radians[i] < (2*PI))) {
              brt="e";
              frt="s";
            }
            sal << "\t" << j << "  " << ke << "  " << i << " [headport=" << brt << ",tailport=" << frt << "];\n";
          }
          else
            sal << "\t" << j << "  " << ke << "  " << i << ";\n";
        }
      }
  }
  sal << "}\n";
  sal.close();
}

void GraphB::print_dot_circ(string archi, double radius)
{
  int i, j;
  ofstream sal;
  string arch;
  arch = archi+".dot";
  basic.open_ofstream(sal, arch);
  string instr, ke,brt,frt;
  if (directed) {
    instr = "digraph G {\n";
    ke = "->";
  }
  else {
    instr = "graph G {\n";
    ke = "--";
  }
  sal << instr;
  sal << "layout = neato;\n";
  sal << "node[shape=circle];\n";
  double radians[size], x[size], y[size], stera, PI=3.141592, van;
  stera = (2*PI)/(size*1.0);
  radians[0] = 0.75*PI;
  van = radians[0]; //
  if (size != 4) {
    for (i=1; i<size; i++) {
      radians[i] = radians[i-1]+stera;
      if (radians[i] >= (2*PI))
        radians[i] = radians[i] - (2*PI);
    }
  }
  else {
    radians[3] = van + stera;
    van = radians[3];
    radians[1] = van +stera;
    van = radians[1];
    radians[2] = van+stera;
    van = radians[2];
    for (i=0; i<size; i++) {
      if (radians[i] >= (2*PI))
        radians[i] = radians[i] - (2*PI);
    }
  }
  
  for (i=0; i<size; i++) {
    basic.polar_to_cartesian(radians[i], radius, x[i], y[i]);
    sal << i << "[pos = \"" << x[i] << "," << y[i] << "!\",label=\"" << i << "\"];\n";
  }
  
  if (directed) {
    for (i=0; i < size; i++)
      for (j=0; j < size; j++) {
        if (weight(j,i)) {
          if (i==j) {
            if ((radians[i] >= 0) && (radians[i] <= (PI/2.0))) {
              brt="n";
              frt="e";
            }
            if ((radians[i] > (PI/2.0)) && (radians[i] <= PI)) {
              brt="w";
              frt="n";
            }
            if ((radians[i] > PI) && (radians[i] <= (PI*1.5))) {
              brt="s";
              frt="w";
            }
            if ((radians[i] > (PI*1.5)) && (radians[i] < (2*PI))) {
              brt="e";
              frt="s";
            }
            sal << "\t" << j << "  " << ke << "  " << i << " [headport=" << brt << ",tailport=" << frt << "];\n";
          }
          else
            sal << "\t" << j << "  " << ke << "  " << i << ";\n";
        }
      }
  }
  else {
    for (i = 0; i < size; i++)
      for (j = 0; j <= i; j++) {
        if (weight(j,i)) {
          if (i==j) {
            if ((radians[i] >= 0) && (radians[i] <= (PI/2.0))) {
              brt="n";
              frt="e";
            }
            if ((radians[i] > (PI/2.0)) && (radians[i] <= PI)) {
              brt="w";
              frt="n";
            }
            if ((radians[i] > PI) && (radians[i] <= (PI*1.5))) {
              brt="s";
              frt="w";
            }
            if ((radians[i] > (PI*1.5)) && (radians[i] < (2*PI))) {
              brt="e";
              frt="s";
            }
            sal << "\t" << j << "  " << ke << "  " << i << " [headport=" << brt << ",tailport=" << frt << "];\n";
          }
          else
            sal << "\t" << j << "  " << ke << "  " << i << ";\n";
        }
      }
  }
  sal << "}\n";
  sal.close();
}

void GraphB::print_cytoscape(string archi)
{
  int i, j;
  ofstream sal;
  string arch;
  if (!basic.string_in_string(".sif", archi))
    arch =archi+".sif";
  else
    arch = archi;
  basic.open_ofstream(sal, arch);
  if (!directed) {
    for (i=0; i<size; i++)
      for (j=0; j<=i; j++)
        if (weight(j,i) || weight(i,j))
          sal << i << "\t" << "n\t" << j << endl;
  }
  else {
    for (i=0; i < size; i++)
      for (j=0; j < size; j++)
        if (weight(j,i))
          sal << j << "\t" << "n\t" << i << endl;
  }
  sal.close();
}

void GraphB::print_cytoscape_wn(string archi)
{
  if (!ionary)
    print_cytoscape(archi);
  else {
    int i, j;
    ofstream sal;
    string arch;
    if (!basic.string_in_string(".sif", archi))
      arch =archi+".sif";
    else
      arch = archi;
    basic.open_ofstream(sal, arch);
    if (!directed) {
      for (i=0; i<size; i++)
        for (j=0; j<=i; j++)
          if (weight(j,i) || weight(i,j))
            sal << dict[i] << "\tn\t" << dict[j] << endl;
    }
    else {
      for (i=0; i < size; i++)
        for (j=0; j < size; j++)
          if (weight(j,i))
            sal << dict[j] << "\tn\t" << dict[i] << endl;
    }
    sal.close();
  }
}

//exclusive
//import
void GraphB::get_dir_nw_from_file(int nn, string arch)
{
  ifstream sal;
  basic.open_ifstream(sal, arch);
  get_dir_nw_from_file(nn, sal);
  sal.close();
}

void GraphB::get_dir_nw_from_file(int nn, istream& en)
{
  make_nw(nn, true);
  int i, j;
  while (en >> j) {
    en >> i;
    if ((i < size) && (j < size) && (i >= 0) && (j >= 0))
      nw[i][j] = true;
    else {
      cout << "[Error]: The network is not specified correctly when GraphB::get_dir_nw_from_file was called.\n";
      exit(1);
    }
  }
  set_default_exclusive_vars();
}

void GraphB::get_undir_nw_from_file(int nn, string arch)
{
  ifstream sal;
  basic.open_ifstream(sal, arch);
  get_undir_nw_from_file(nn, sal);
  sal.close();
}

void GraphB::get_undir_nw_from_file(int nn, istream& en)
{
  make_nw(nn, false);
  int i, j;
  while (en >> j) {
    en >> i;
    if ((i < size) && (j < size) && (i >= 0) && (j >= 0)) {
      nw[i][j] = true;
      nw[j][i] = true;
    }
    else {
      cout << "[Error]: This network is not specified correctly when GraphB::get_undir_nw_from_file was called.\n";
      exit(1);
    }
  }
  set_default_exclusive_vars();
}

//desde aqui
void GraphB::get_dir_nw_from_file_wn(int nn, string arch, string archnam)
{
  ifstream sal, salna;
  basic.open_ifstream(salna, archnam);
  basic.open_ifstream(sal, arch);
  get_dir_nw_from_file_wn(nn, sal, salna);
  sal.close();
  salna.close();
}

void GraphB::get_dir_nw_from_file_wn(int nn, istream& en, istream& enna)
{
  make_nw(nn, true);
  get_names_from_file(enna);
  string i, j;
  int k, l;
  while (en >> j) {
    en >> i;
    for (k=0; k < size; k++)
      if (j==dict[k])
        break;
    for (l=0; l <size; l++)
      if (i==dict[l])
        break;
    if ((k==size) || (l==size)) {
      cout << "[Error]: " << j << " or " << i << " not in the graph's node list. This error arose when calling GraphB::get_dir_nw_from_file_wn.\n";
      exit(1);
    }
    nw[l][k] = true;
  }
  set_default_exclusive_vars();
}

void GraphB::get_dir_nw_from_file_wn(string arch) {
  ifstream fe;
  string basu, *tempo;
  tempo = new string[3000];
  int i,k,l, van = 0;
  basic.open_ifstream(fe, arch);
  while (fe >> basu) {
    for (i = 0; i < van; i++)
      if (basu == tempo[i])
        break;
    if (i == van) {
      tempo[van] = basu;
      van++;
    }
  }
  fe.close();
  make_nw(van, true);
  if (ionary) {
    cout << "[Error]: Names have already been assigned. This error arose when GraphB::get_undir_nw_from_file_wn was called.\n";
    exit(1);
  }
  string sj, si;
  dict = new string[size];
  for (i=0; i<size; i++)
    dict[i] = tempo[i];
  ionary = true;
  delete [] tempo;
  basic.open_ifstream(fe, arch);
  while (fe >> sj) {
    fe >> si;
    for (k=0; k < size; k++)
      if (sj==dict[k])
        break;
    for (l=0; l <size; l++)
      if (si==dict[l])
        break;
    if ((k==size) || (l==size)) {
      cout << "[Error]: " << sj << " or " << si << " not in the graph's node list. This error arose when GraphB::get_undir_nw_from_file_wn was called.\n";
      exit(1);
    }
    nw[l][k] = true;
  }
  set_default_exclusive_vars();
  fe.close();
}

void GraphB::get_undir_nw_from_file_wn(string arch) {
  ifstream fe;
  string basu, *tempo;
  tempo = new string[3000];
  int i,k,l, van = 0;
  basic.open_ifstream(fe, arch);
  while (fe >> basu) {
    for (i = 0; i < van; i++)
      if (basu == tempo[i])
        break;
    if (i == van) {
      tempo[van] = basu;
      van++;
    }
  }
  fe.close();
  make_nw(van, false);
  if (ionary) {
    cout << "[Error]: Names have already been assigned. This error arose when GraphB::get_undir_nw_from_file_wn was called.\n";
    exit(1);
  }
  string si, sj;
  dict = new string[size];
  for (i=0; i<size; i++)
    dict[i] = tempo[i];
  ionary = true;
  delete [] tempo;
  basic.open_ifstream(fe, arch);
  while (fe >> sj) {
    fe >> si;
    for (k=0; k < size; k++)
      if (sj==dict[k])
        break;
    for (l=0; l <size; l++)
      if (si==dict[l])
        break;
    if ((k==size) || (l==size)) {
      cout << "[Error]: " << sj << " or " << si << " not in the graph's node list. This error arose when GraphB::get_undir_nw_from_file_wn was called.\n";
      exit(1);
    }
    nw[l][k] = true;
    nw[k][l] = true;
  }
  set_default_exclusive_vars();
  fe.close();
}



void GraphB::get_undir_nw_from_file_wn(int nn, string arch, string archnam)
{
  ifstream sal, salna;
  basic.open_ifstream(salna, archnam);
  basic.open_ifstream(sal, arch);
  get_undir_nw_from_file_wn(nn, sal, salna);
  sal.close();
  salna.close();
}

void GraphB::get_undir_nw_from_file_wn(int nn, istream& en, istream& enna)
{
  make_nw(nn, false);
  get_names_from_file(enna);
  string i, j;
  int k, l;
  while (en >> j) {
    en >> i;
    for (k=0; k < size; k++)
      if (j==dict[k])
        break;
    for (l=0; l <size; l++)
      if (i==dict[l])
        break;
    if ((k==size) || (l==size)) {
      cout << "[Error]: " << j << " or " << i << " not in the graph's node list. This error arose when GraphB::get_undir_nw_from_file_wn was called.\n";
      exit(1);
    }
    nw[l][k] = true;
    nw[k][l] = true;
  }
  set_default_exclusive_vars();
}


void GraphB::get_from_graphi(GraphI &templ)
{
  int n = templ.number_of_nodes();
  make_nw(n, templ.is_directed());
  int i, j;
  for (i = 0; i < size; i++)
    for (j = 0; j < size; j++)
      if (templ.weight(j, i) != 0)
        nw[i][j] = true;
  if (templ.names_exist()) {
    dict = new string[size];
    for (i = 0; i < size; i++)
      dict[i] = templ.get_name(i);
    ionary = true;
  }
  set_default_exclusive_vars();
}

void GraphB::get_from_graphc(GraphC &templ)
{
  int n = templ.number_of_nodes(); //
  make_nw(n, templ.is_directed()); //
  int i, j;
  for (i = 0; i < size; i++)
    for (j = 0; j < size; j++)
      if (templ.weight(j, i) != 0) //
        nw[i][j] = true;
  if (templ.names_exist()) { //
    dict = new string[size];
    for (i = 0; i < size; i++)
      dict[i] = templ.get_name(i); //
    ionary = true;
  }
  set_default_exclusive_vars();
}

//Export
void GraphB::export_nw(string arch) {
  ofstream fs;
  basic.open_ofstream(fs, arch);
  export_nw(fs);
  fs.close();
}

void GraphB::export_nw(ostream& fs) {
  int i, j;
  if (directed) {
    for (i = 0; i < size; i++)
      for (j=0; j < size; j++)
        if (weight(j,i))
          fs << j << "\t" << i << "\t" << weight(j, i) << endl;
  }
  else {
    for (i = 0; i < size; i++)
      for (j=0; j <= i; j++)
        if (weight(j,i))
          fs << j << "\t" << i << "\t" << weight(j, i) << endl;
  }
}

//status

//analyses
bool GraphB::equal_nw(GraphB &templ)
{
  int i, j;
  bool res = true;
  if ((templ.number_of_nodes() != size) || (templ.is_directed() != is_directed()))
    res = false;
  else {
    for (i = 0; i < size; i++) {
      for (j=0; j<size; j++)
        if (templ.weight(j, i) != weight(j, i)) {
          res = false;
          break;
        }
      if (!res)
        break;
    }
  }
  return res;
}

bool GraphB::isomorphic_pair(GraphB &templ)
{
  bool res = false;
  if (!yamotdab)
    load_database_for_isomorphism();
  if (size > 4) {
    cout << "[Error]: GraphB::isomorphic_pair does not work for graphs with more than 4 nodes.\n";
    exit(1);
  }
  if ((number_of_nodes() == templ.number_of_nodes()) && (templ.is_directed() == is_directed()) && (number_of_edges() == templ.number_of_edges())) {
    if (size==1)
      res = isomorphic_pair1(templ);
    if (size==2)
      res = isomorphic_pair2(templ);
    if (size==3)
      res = isomorphic_pair3(templ);
    if (size==4)
      res = isomorphic_pair4(templ);
  }
  return res;
}

bool GraphB::isomorphic_pair(GraphB &templ, int *ord)
{
  bool res = false;
  if (!yamotdab)
    load_database_for_isomorphism();
  if (size > 4) {
    cout << "[Error]: GraphB::isomorphic_pair does not work for graphs with more than 4 nodes.\n";
    exit(1);
  }
  if ((number_of_nodes() == templ.number_of_nodes()) && (templ.is_directed() == is_directed()) && (number_of_edges() == templ.number_of_edges())) {
    if (size==1)
      res = isomorphic_pair1(templ, ord);
    if (size==2)
      res = isomorphic_pair2(templ, ord);
    if (size==3)
      res = isomorphic_pair3(templ, ord);
    if (size==4)
      res = isomorphic_pair4(templ, ord);
  }
  return res;
}

set<int> GraphB::get_upstream(int no)
{
  int i;
  set<int> res;
  res.clear();
  for (i=0; i<size; i++)
    if (nw[no][i])
      res.insert(i);
  return res;
}

void GraphB::get_downstream()
{
  downstream = new set<int>[size];
  int i, j;
  for (i = 0; i < size; i++) {
    for (j = 0; j < size; j++) {
      if (nw[j][i])
        downstream[i].insert(j);
    }
  }
  infl = true;
}

set<int> GraphB::get_out_component(int no) {
  set<int> res;
  res.clear();
  int i;
  for (i=0; i< size; i++)
    if (there_is_path(no, i))
      res.insert(i);
  return res;
}

set<int> GraphB::get_in_component(int no) {
  set<int> res;
  res.clear();
  int i;
  for (i=0; i< size; i++)
    if (there_is_path(i, no))
      res.insert(i);
  return res;
}

void GraphB::get_adjacency_matrix() {
  if (yamadya) {
    cout << "[Error]: Adjacency matrix already created when GraphB::get_adjacency_matrix was called.\n";
    exit(1);
  }
  basic.create_array(matadya, size, size);
  int i,j;
  for (i = 0; i < size; i++)
    for (j = 0; j < size; j++) {
      if (nw[i][j])
        matadya[i][j] = 1;
      else
        matadya[i][j] = 0;
    }
  if (!directed) {
    for (i = 0; i < size; i++)
      if (matadya[i][i] != 0)
        matadya[i][i]++;
  }
  yamadya = true;
  return;
}

void GraphB::get_distance_matrix_ls()
{
  if (dist) {
    cout << "[Error]: Distance matrix already created. GraphB::get_distance_matrix.\n";
    exit(1);
  }
  int i, j;
  dima = new int*[size];
  for (i = 0; i < size; i++)
    dima[i] = new int[size];
  basic.fillmatm1(dima, size, size);
  if (!infl)
    get_downstream();
  set<int> visited;
  set<int> next;
  set<int> nextnext;
  set<int> toerase;
  set<int>::iterator it;
  int distance;
  for (i = 0; i < size; i++) {
    distance = 0;
    visited.clear();
    next.clear();
    nextnext.clear();
    toerase.clear();
    visited.insert(i);
    dima[i][i] = distance;
    distance++;
    next = downstream[i];
    for (it=next.begin(); it!=next.end(); it++)
      if (visited.count(*it)>0)
        toerase.insert(*it);
    for (it=toerase.begin(); it!=toerase.end(); it++)
      next.erase(*it);
    toerase.clear();
    do {
      for (it=next.begin(); it!=next.end(); it++) {
        nextnext = basic.merge(nextnext, downstream[*it]);
        if (dima[*it][i] != -1) {
          cout << "[Error]: GraphB::get_distance_matrix.\n";
          exit(1);
        }
        dima[*it][i] = distance;
      }
      for (it=next.begin(); it!=next.end(); it++)
        visited.insert(*it);
      next.clear();
      for (it=nextnext.begin(); it!=nextnext.end(); it++)
        if (visited.count(*it)>0)
          toerase.insert(*it);
      for (it=toerase.begin(); it!=toerase.end(); it++)
        nextnext.erase(*it);
      toerase.clear();
      next = nextnext;
      nextnext.clear();
      distance++;
    } while ((!(next.empty())) && (distance < size) && ((signed (visited.size() + next.size())) < size));
    for (it=next.begin(); it!=next.end(); it++) {
      if (dima[*it][i] != -1) {
        cout << "[Error]: GraphB::get_distance_matrix.\n";
        exit(1);
      }
      dima[*it][i] = distance;
    }
  }
  dist = true;
  if (!(directed)) {
    for (i = 0; i < size; i++)
      for (j = 0; j <= i; j++)
        if (dima[i][j] != dima[j][i]) {
          cout << "[Error]: GraphB::get_distance_matrix.\n";
          exit(1);
        }
  }
}

void GraphB::get_distance_matrix(){
  if (dist) {
    cout << "[Error]: Distance matrix already created when GraphB::get_distance_matrix was called.\n";
    exit(1);
  }
  int i,j,k;
  basic.create_array(dima, size, size);
  basic.fillmatm1(dima, size, size);
  int **temain, **temaout;
  if (!yamadya)
    get_adjacency_matrix();
  basic.create_array(temain, size, size);
  basic.fillmat0(temain, size, size);
  
  basic.create_array(temaout, size, size);
  basic.fillmat0(temaout, size, size);
  
  for (i = 0; i < size; i++) {
    for (j = 0; j < size; j++)
      if (nw[i][j]) /////cambiar!!
        dima[i][j] = 1;
    dima[i][i] = 0;
  }
  basic.matxmat(matadya, size, size, matadya, size, size, temaout, size, size);
  k = 2;
  bool to = false;
  for (i = 0; i < size; i++)
    for (j = 0; j < size; j++) {
      temain[i][j] = temaout[i][j];
      if ((temaout[i][j] > 0) && (dima[i][j] < 0)) {
        dima[i][j] = k;
        to = true;
      }
    }
  k = 3;
  while (to) {
    to = false;
    basic.matxmat(temain, size, size, matadya, size, size, temaout, size, size);
    for (i = 0; i < size; i++)
      for (j = 0; j < size; j++) {
        temain[i][j] = temaout[i][j];
        if ((temaout[i][j] > 0) && (dima[i][j] < 0)) {
          dima[i][j] = k;
          to = true;
        }
      }
    k++;
  }
  
  for (i = 0; i < size; i++) {
    delete [] temain[i];
    delete [] temaout[i];
  }
  delete [] temain;
  delete [] temaout;
  dist = true;
}

int GraphB::get_scc() {
  int i, j;
  if (!directed) {
    cout << "[Error]: No SCC's in undirected networks. This error arose during execution of GraphB::get_scc.\n";
    exit(1);
  }
  if (cfc) {
    cout << "[Error]: SCC's already obtained when GraphB::get_scc was called.\n";
    exit(1);
  }
  if (!dist)
    get_distance_matrix();
  bool mu = true;
  numbscc = 0;
  i=0;
  scc = new int[size];
  basic.fillvm1(scc, size);
  while(mu) {
    scc[i] = numbscc;
    for (j = i+1; j < size; j++)
      if ((dima[i][j] > 0) && (dima[j][i] > 0))
        scc[j] = numbscc;
    numbscc++;
    mu = false;
    for (j = i+1; j < size; j++) {
      if (scc[j] < 0) {
        mu = true;
        i = j;
        break;
      }
    }
  }
  cfc = true;
  return numbscc;
}

int GraphB::number_of_sccs() {
  if (!cfc)
    get_scc();
  return numbscc;
}

int GraphB::number_of_components() {
  if (!yacomp)
    get_components();
  return numbcomp;
}

int GraphB::get_components()
{
  int i, j;
  if (!(yacomp)) {
    bool mu = true;
    numbcomp = 0;
    set<int> enco;
    set<int> sig;
    set<int> sig2;
    set<int> lerase;
    set<int>::iterator it;
    enco.clear();
    sig.clear();
    sig2.clear();
    lerase.clear();
    components = new int[size];
    if (!dist)
      get_distance_matrix();
    basic.fillvm1(components, size);
    mu = false;
    for (i = 0; i < size; i++)
      if (components[i] == (-1)) {
        mu = true;
        break;
      }
    while (mu) {
      enco.clear();
      sig.clear();
      sig2.clear();
      lerase.clear();
      for (i = 0; i < size; i++)
        if (components[i] == -1)
          break;
      enco.insert(i);
      for (j = 0; j < size; j++)
        if ((dima[i][j] > (-1)) || (dima[j][i] > (-1)))
          if (j != i)
            sig.insert(j);
      while (!sig.empty()) {
        for (it=sig.begin(); it!=sig.end(); it++)
          for (j = 0; j < size; j++)
            if ((dima[*it][j] > (-1)) || (dima[j][*it] > (-1)))
              if ((enco.count(j)==0) && (sig.count(j)==0) && (components[j]==-1))
                sig2.insert(j);
        for (it=sig.begin(); it!=sig.end(); it++)
          enco.insert(*it);
        sig.clear();
        for (it=sig2.begin(); it!=sig2.end(); it++)
          sig.insert(*it);
        sig2.clear();
      }
      for (it=sig.begin(); it!=sig.end(); it++)
        enco.insert(*it);
      numbcomp++;
      for (it=enco.begin(); it!=enco.end(); it++) {
        if (components[*it] != -1) {
          cout << "[Error]: There was an error during execution of GraphB::get_components.\n";
          exit(1);
        }
        components[*it] = numbcomp-1;
      }
      mu = false;
      for (i = 0; i < size; i++)
        if (components[i] == (-1)) {
          mu = true;
          break;
        }
    }
    yacomp = true;
  }
  return numbcomp;
}


void GraphB::get_bc()
{
  int i, j, dicor, k, l;
  double toshopa = 0;
  int nuste;
  set<int> pagre;
  set<int>::iterator it;
  bc = new double[size];
  if (!dist)
    get_distance_matrix();
  basic.fillv0(bc, size);
  for (i = 0; i < size; i++)
    for (j = 0; j < size; j++) {
      if (dima[i][j] > -1) {
        toshopa = toshopa + 1.0;
        dicor = dima[i][j];
        for (k = 1; k < dicor; k++) {
          nuste = 0;
          pagre.clear();
          for (l = 0; l < size; l++)
            if ((dima[l][j] == k) && (dima[i][l] == (dicor-k)))
              pagre.insert(l);
          nuste = pagre.size();
          for (it=pagre.begin(); it!=pagre.end(); it++)
            bc[*it] = bc[*it] + 1.0/(double(nuste));
        }
      }
    }
  for (i = 0; i < size; i++)
    bc[i] = bc[i]/toshopa;
  bcya = true;
}

double GraphB::get_bc(int nodo)
{
  if (!bcya)
    get_bc();
  return bc[nodo];
}

int GraphB::nupaths(int j, int q, int *vec)
{
  //j is source
  int i, k, l, v, res,guaja;
  res = 0;
  if (!dist)
    get_distance_matrix();
  list<int> Omega;
  list<int> Lambda;
  list<int> *succ;
  list<int>::iterator it;
  bool una;
  bool self = false;
  succ = new list<int>[size];
  for (i=0; i <size; i++)
    succ[i].clear();
  v = j; //step1
  basic.fillv0(vec, size);
  Omega.clear();
  Lambda.clear();
  Omega.push_back(v);
  if (nw[v][v]) {
    Lambda.push_back(v);
    self = true;
  }
  for (k = 0; k < size; k++)
    if (dima[k][v] > 0)
      Lambda.push_back(k);
  
  if (q==v)
    for (k = 0; k < size; k++)
      if (k != q)
        if ((dima[k][q] > 0) && (dima[q][k] > 0) && (!self)) {
          Lambda.push_back(v);
          break;
        }
  if (basic.contains(q, Lambda)) {
    res = 0;
    una = true;
    do { //step2
      if (una) {
        succ[v].clear();
        succ[v] = GetSucc(v, Omega, q);
      }
      una = false; //step 3
      if (!(succ[v].empty())) {
        for (it = succ[v].begin(); it != succ[v].end(); it++) {
          una = false; //step 4,3
          if (*it == q) {
            succ[v].erase(it);
            l = Omega.size() - 1;
            vec[l]++; // here to be modified to identify path sign
            it--; //melate
            break;
          }
          else {
            guaja = *it;
            Omega.push_back(*it);
            succ[v].erase(it);
            it--; //melate
            v = guaja;
            una = true;
            break;
          }
        }
      }
      else { //step 5
        Omega.pop_back();
        if (!(Omega.empty())) {
          it = Omega.end();
          it--; //melate
          v = *it;
        }
      }
    } while (!(Omega.empty()));
    for (l = 0; l < size; l++)
      res = res + vec[l];
  }
  delete [] succ;
  return res;
}

list<int> GraphB::GetSucc(int v, list<int> &omega, int q)
{
  list<int> succ;
  list<int> T;
  list<int> temp;
  list<int>::iterator it;
  list<int>::iterator itd;
  T.clear();
  int u, r;
  //step 2.1
  int c[size];
  int d[size];
  int i;
  basic.fillv0(c, size);
  basic.fillv0(d, size);
  for (it = omega.begin(); it != omega.end(); it++)
    c[*it] = 1;
  T.push_back(q);
  c[q] = 1;
  d[q] = 1;
  //   //step2.4
  while(!(T.empty())) {
    //       step2.3
    for (it = T.begin(); it != T.end(); it++) {
      r = *it;
      temp.clear();
      for (i = 0; i < size; i++)
        if (nw[r][i])
          temp.push_back(i);
      T.erase(it);
      it--; //melate
      for (itd = temp.begin(); itd != temp.end(); itd++) {
        u = *itd;
        if (c[u] == 0) {
          T.push_back(u);
          c[u] = 1;
          d[u] = 1;
        }
      }
    }
  }
  //   //step2.5
  succ.clear();
  temp.clear();
  for (i = 0; i < size; i++)
    if (nw[i][v])
      temp.push_back(i);
  for (itd = temp.begin(); itd != temp.end(); itd++) {
    u = *itd;
    if (d[u] == 1)
      succ.push_back(u);
  }
  return succ;
}

void GraphB::get_loops()
{
  totalloops = 0;
  int i, j;
  if (!floops) {
    loops = new int[size];
    guada = new int*[size];
    for (i = 0; i < size; i++)
      guada[i] = new int[size];
  }
  basic.fillmat0(guada, size, size);
  basic.fillv0(loops, size);
  for (i = 0; i < size; i++)
    nupaths(i, i, guada[i]);
  for (i = 0; i < size; i++)
    for (j = 0; j < size; j++)
      loops[i] = loops[i] + guada[j][i];
  for (i = 0; i < size; i++) {
    if ((loops[i]%(i+1)) != 0) {
      cout << "[Error]: This error arose during execution of GraphB::get_loops.\n";
      exit(1);
    }
    loops[i] = loops[i]/(i+1);
  }
  for (i = 0; i < size; i++)
    totalloops = totalloops + loops[i];
  floops = true;
}

//acces analyses
bool GraphB::there_is_path(int from, int to)
{
  if (!dist)
    get_distance_matrix();
  bool res;
  if (dima[to][from] > -1)
    res = true;
  else
    res = false;
  return res;
}

int GraphB::length_shortest_path(int from, int to)
{
  int res;
  if (!dist)
    get_distance_matrix();
  res = dima[to][from];
  return res;
}

set<int> GraphB::nodes_influenced_by(int n)
{
  if (!infl)
    get_downstream();
  return downstream[n];
}

int GraphB::in_which_component(int n)
{
  if (!yacomp)
    get_components();
  return components[n];
}

int GraphB::in_which_scc(int n)
{
  if (!cfc)
    get_scc();
  return scc[n];
}

int GraphB::number_of_paths(int from, int to, int len, int* vec)
{
  nupaths(from, to, vec);
  return vec[len-1];
}

int GraphB::number_of_loops()
{
  if (!floops)
    get_loops();
  return totalloops;
}
int GraphB::number_of_loops(int len)
{
  if (!floops)
    get_loops();
  return loops[len];
}

void GraphB::switches(int vec)
{
  if (infl)
    clear_ds();
  if (dist)
    clear_dima();
  if (yacomp)
    clear_comp();
  if (cfc)
    clear_scc();
  if (bcya)
    clear_bc();
  if (floops)
    clear_loops();
  if (yamadya)
    clear_adj();
  int i;
  for (i = 0; i < vec; i++)
    aswitch();
}

void GraphB::switches_preserving_sg1(int vec)
{
  if (infl)
    clear_ds();
  if (dist)
    clear_dima();
  if (yacomp)
    clear_comp();
  if (cfc)
    clear_scc();
  if (bcya)
    clear_bc();
  if (floops)
    clear_loops();
  int i;
  for (i = 0; i < vec; i++)
    aswitch_preserving_sg1();
}

void GraphB::switches_preserving_sg2(int vec)
{
  if (infl)
    clear_ds();
  if (dist)
    clear_dima();
  if (yacomp)
    clear_comp();
  if (cfc)
    clear_scc();
  if (bcya)
    clear_bc();
  if (floops)
    clear_loops();
  int i;
  for (i = 0; i < vec; i++)
    aswitch_preserving_sg2();
}

void GraphB::print_list_of_list_of_nodes(list<list<int> > &lalis, ostream& sal)
{
  if (!ionary) {
    cout << "[Error]: Names have not been assigned before calling GraphB::print_list_of_list_of_nodes\n";
    exit(1);
  }
  list<list<int> >::iterator whi;
  list<int> mien;
  list<int>::iterator it;
  sal << "\[";
  for (whi=lalis.begin(); whi!=lalis.end(); whi++) {
    sal << "(";
    mien.clear();
    mien = *whi;
    for (it=mien.begin(); it != mien.end(); it++)
      sal << get_name(*it) << " ";
    sal << ")";
  }
  sal << "]\n";
}

//for evolution
void GraphB::mate(GraphB &mother, GraphB &father)
{
  if (mother.number_of_nodes() != father.number_of_nodes()) {
    cout << "[Error]: Parents with different genome size. GraphB::mate.\n";
    exit(1);
  }
  if (mother.is_directed() != father.is_directed()) {
    cout << "[Error]: One parent is a directed network while the other is an undirected network. GraphB::mate.\n";
    exit(1);
  }
  int i, j;
  bool paoma;
  int n = father.number_of_nodes();
  make_nw(n, father.is_directed());
  if (directed) {
    for (i = 0; i < size; i++) {
      paoma = est.toss();
      if (paoma) {
        for (j = 0; j < size; j++)
          nw[i][j] = father.weight(j, i);
      }
      else {
        for (j = 0; j < size; j++)
          nw[i][j] = mother.weight(j, i);
      }
    }
  }
  else {
    for (i=0; i<size; i++)
      for (j=0; j<=i; j++) {
        paoma = est.toss();
        if (paoma)
          nw[i][j] = father.weight(j,i);
        else
          nw[i][j] = mother.weight(j,i);
      }
    for (i=0; i<size; i++)
      for (j=0; j<i; j++)
        nw[j][i] = nw[i][j];
  }
}

void GraphB::mutate(int gene, double wol)
{
  cout << "[Error]: Mutations are non-sense in the absence of dynamics. GraphB::mutate.\n";
  exit(1);
}

void GraphB::duplicate(int gene)
{
  bool **matemp;
  int i,j;
  int n = size+1;
  matemp = new bool*[n];
  for (i = 0; i < n; i++)
    matemp[i] = new bool[n];
  basic.fillmat0(matemp, n, n);
  for (i = 0; i < size; i++)
    for (j = 0; j < size; j++)
      matemp[i][j] = nw[i][j];
  bool dirnw = is_directed(); //
  clear();
  make_nw(n, dirnw); //
  for (i = 0; i < size; i++)
    for (j = 0; j < size; j++)
      nw[i][j] = matemp[i][j];
  for (j = 0; j < (size-1); j++) {
    nw[size-1][j] = nw[gene][j];
    nw[j][size-1] = nw[j][gene];
  }
  nw[size-1][size-1] = nw[gene][gene];
  for (i = 0; i < size; i++)
    delete [] matemp[i];
  delete [] matemp;
}

void GraphB::duplicate()
{
  int i = est.randint(0, size);
  duplicate(i);
}

//motifs
int GraphB::count_selfinteracting_nodes()
{
  int i,res =0;
  for (i=0; i<size; i++)
    if (weight(i,i)) //
      res++;
  return res;
}

int GraphB::count_non_selfinteracting_nodes()
{
  int i,res =0;
  for (i=0; i<size; i++)
    if (!weight(i,i)) //
      res++;
  return res;
}

//motif exclusive of graphb
bool GraphB::contains_these_interactions(GraphB &subg)
{
  int i,j;
  bool res = true;
  if ((subg.number_of_nodes() != number_of_nodes()) || (subg.is_directed() != is_directed()))
    res= false;
  else {
    for (i=0; i<number_of_nodes(); i++) {
      for (j=0; j<number_of_nodes(); j++)
        if (subg.weight(j, i) && !weight(j,i)) {
          res=false;
          break;
        }
      if (!res)
        break;
    }
  }
  return res;
}

bool GraphB::an_isomorph_contains_these_interactions(GraphB &subg)
{
  bool res=false;
  if (!yamotdab)
    load_database_for_isomorphism();
  if (size > 4) {
    cout << "[Error]: GraphB::an_isomorph_contains_these_interactions does not work for graphs with more than 4 nodes.\n";
    exit(1);
  }
  if ((number_of_nodes() == subg.number_of_nodes()) && (subg.is_directed() == is_directed())) {
    if (size==1)
      res = an_isomorph_contains_these_interactions1(subg);
    if (size==2)
      res = an_isomorph_contains_these_interactions2(subg);
    if (size==3)
      res = an_isomorph_contains_these_interactions3(subg);
    if (size==4)
      res = an_isomorph_contains_these_interactions4(subg);
  }
  return res;
}

void GraphB::load_database_for_isomorphism()
{
  int i,j,fac;
  comb2 = new int*[2];
  for (i=0; i<2; i++)
    comb2[i] = new int[2];
  comb2[0][0] = 0;
  comb2[0][1] = 1;
  comb2[1][0] = 1;
  comb2[1][1] = 0;
  fac = basic.factorial(3);
  comb3 = new int*[fac];
  for (i=0; i<fac; i++)
    comb3[i] = new int[3];
  fac = basic.factorial(4);
  comb4 = new int*[fac];
  for (i=0; i<fac; i++)
    comb4[i] = new int[4];
  int **rings3, **rings4;
  rings3 = new int*[2];
  for (i=0; i<2; i++)
    rings3[i] = new int[3];
  rings4 = new int*[6];
  for (i=0; i<6; i++)
    rings4[i] = new int[4];
  rings3[0][0] = 0;
  rings3[0][1] = 1;
  rings3[0][2] = 2;
  rings3[1][0] = 1;
  rings3[1][1] = 0;
  rings3[1][2] = 2;
  for (i=0; i<2; i++) {
    for (j=0; j<3; j++) {
      comb3[i*3][j] = rings3[i][j];
      comb3[(i*3)+1][(j+1)%3] = rings3[i][j];
      comb3[(i*3)+2][(j+2)%3] = rings3[i][j];
    }
  }
  rings4[0][0] = 0;
  rings4[0][1] = 1;
  rings4[0][2] = 2;
  rings4[0][3] = 3;
  rings4[1][0] = 0;
  rings4[1][1] = 2;
  rings4[1][2] = 1;
  rings4[1][3] = 3;
  rings4[2][0] = 0;
  rings4[2][1] = 1;
  rings4[2][2] = 3;
  rings4[2][3] = 2;
  rings4[3][0] = 0;
  rings4[3][1] = 3;
  rings4[3][2] = 2;
  rings4[3][3] = 1;
  rings4[4][0] = 0;
  rings4[4][1] = 2;
  rings4[4][2] = 3;
  rings4[4][3] = 1;
  rings4[5][0] = 0;
  rings4[5][1] = 3;
  rings4[5][2] = 1;
  rings4[5][3] = 2;
  for (i=0; i<6; i++) {
    for (j=0; j<4; j++) {
      comb4[i*4][j] = rings4[i][j];
      comb4[(i*4)+1][(j+1)%4] = rings4[i][j];
      comb4[(i*4)+2][(j+2)%4] = rings4[i][j];
      comb4[(i*4)+3][(j+3)%4] = rings4[i][j];
    }
  }
  for (i = 0; i < 2; i++)
    delete [] rings3[i];
  delete [] rings3;
  for (i = 0; i < 6; i++)
    delete [] rings4[i];
  delete [] rings4;
  yamotdab = true;
}

void GraphB::load_sg_catalogue()
{
  int i;
  if (directed) {
#include "metacode/Directed/2.cc"
#include "metacode/Directed/3.cc"
  }
  else {
#include "metacode/Undirected/2.cc"
#include "metacode/Undirected/3.cc"
#include "metacode/Undirected/4.cc"
  }
  yasgcat = true;
}

void GraphB::prepare_to_count_subgraphs()
{
  if (!yasgcat)
    load_sg_catalogue();
  vedosg2 = new int[nupre_2];
  basic.fillv0(vedosg2, nupre_2);
  vedosg3 = new int[nupre_3];
  basic.fillv0(vedosg3, nupre_3);
  vedosg4 = new int[nupre_4];
  basic.fillv0(vedosg4, nupre_4);
  yaprecosg = true;
}

int GraphB::loose_count_this_subgraph(GraphB &sg)
{
  int res, nu=sg.number_of_nodes();
  if ( nu> 4) {
    cout << "[Error]: Subgraphs of size greater than 4 can not be counted. GraphB::loose_count_this_subgraph.\n";
    exit(1);
  }
  if (sg.get_components() != 1) {
    cout << "[Error]: Subgraph has more than one component. GraphB::loose_count_this_subgraph.\n";
    exit(1);
  }
  if (nu==2)
    res = loose_count_this_subgraph2(sg);
  else {
    if (nu == 3) {
      res = loose_count_this_subgraph3(sg);
    }
    else {
      res = loose_count_this_subgraph4(sg);
    }
  }
  return res;
}

int GraphB::count_this_subgraph(GraphB &sg)
{
  int res, nu=sg.number_of_nodes();
  if ( nu> 4) {
    cout << "[Error]: Subgraphs of size greater than 4 can not be counted. GraphB::count_this_subgraph.\n";
    exit(1);
  }
  if (sg.get_components() != 1) {
    cout << "[Error]: Subgraph has more than one component. GraphB::count_this_subgraph.\n";
    exit(1);
  }
  if (nu==2)
    res = count_this_subgraph2(sg);
  else {
    if (nu == 3) {
      res = count_this_subgraph3(sg);
    }
    else {
      res = count_this_subgraph4(sg);
    }
  }
  return res;
}

int GraphB::list_subgraph_appearances(GraphB &sg, string arch)
{
  int res;
  ofstream sal;
  basic.open_ofstream(sal, arch);
  res = list_subgraph_appearances(sg, sal);
  sal.close();
  return res;
}


int GraphB::list_subgraph_appearances(GraphB &sg, ostream& os)
{
  int res, nu=sg.number_of_nodes();
  if ( nu> 4) {
    cout << "[Error]: Subgraphs of size greater than 4 can not be counted. GraphB::list_subgraph_appearances.\n";
    exit(1);
  }
  if (sg.get_components() != 1) {
    cout << "[Error]: Subgraph has more than one component. GraphB::list_subgraph_appearances.\n";
    exit(1);
  }
  if (nu==2)
    res = list_subgraph_appearances2(sg, os);
  else {
    if (nu == 3) {
      res = list_subgraph_appearances3(sg, os);
    }
    else {
      if (nu==4)
        res = list_subgraph_appearances4(sg, os);
      else {
        cout << "[Error]: GraphB::list_subgraph_appearances does not work for subgraphs of size greater than 4.\n";
        exit(1);
      }
    }
  }
  return res;
}

int GraphB::list_subgraph_appearances(GraphB &sg, string arch, bool **me, int ge, int tis)
{
  int res;
  ofstream sal;
  basic.open_ofstream(sal, arch);
  res = list_subgraph_appearances(sg, sal, me, ge, tis);
  sal.close();
  return res;
}


int GraphB::list_subgraph_appearances(GraphB &sg, ostream& os, bool **me, int ge, int tis)
{
  int res, nu=sg.number_of_nodes();
  if ( nu> 4) {
    cout << "[Error]: Subgraphs of size greater than 4 can not be counted. GraphB::list_subgraph_appearances.\n";
    exit(1);
  }
  if (sg.get_components() != 1) {
    cout << "[Error]: Subgraph has more than one component. GraphB::list_subgraph_appearances.\n";
    exit(1);
  }
  if (nu==2)
    res = list_subgraph_appearances2(sg, os, me, ge, tis);
  else {
    if (nu == 3) {
      res = list_subgraph_appearances3(sg, os, me, ge, tis);
    }
    else {
      if (nu==4)
        res = list_subgraph_appearances4(sg, os, me, ge, tis);
      else {
        cout << "[Error]: GraphB::list_subgraph_appearances does not work for subgraphs of size greater than 4.\n";
        exit(1);
      }
    }
  }
  return res;
}

int GraphB::list_subgraph_appearances_using_names(GraphB &sg, string arch)
{
  int res;
  ofstream sal;
  basic.open_ofstream(sal, arch);
  res = list_subgraph_appearances_using_names(sg, sal);
  sal.close();
  return res;
}

int GraphB::list_subgraph_appearances_using_names(GraphB &sg, ostream& os)
{
  int res, nu=sg.number_of_nodes();
  if ( nu> 4) {
    cout << "[Error]: Subgraphs of size greater than 4 can not be counted. GraphB::list_subgraph_appearances_using_names\n";
    exit(1);
  }
  if (sg.get_components() != 1) {
    cout << "[Error]: Subgraph has more than one component. GraphB::list_subgraph_appearances_using_names.\n";
    exit(1);
  }
  if (nu==2)
    res = list_subgraph_appearances_using_names2(sg, os);
  else {
    if (nu == 3) {
      res = list_subgraph_appearances_using_names3(sg, os);
    }
    else {
      res = list_subgraph_appearances_using_names4(sg, os);
    }
  }
  return res;
}
//

int GraphB::list_subgraph_appearances_using_names(GraphB &sg, string arch, bool **me, int ge, int tis)
{
  int res;
  ofstream sal;
  basic.open_ofstream(sal, arch);
  res = list_subgraph_appearances_using_names(sg, sal, me, ge, tis);
  sal.close();
  return res;
}

int GraphB::list_subgraph_appearances_using_names(GraphB &sg, ostream& os, bool **me, int ge, int tis)
{
  int res, nu=sg.number_of_nodes();
  if ( nu> 4) {
    cout << "[Error]: Subgraphs of size greater than 4 can not be counted. GraphB::list_subgraph_appearances_using_names\n";
    exit(1);
  }
  if (sg.get_components() != 1) {
    cout << "[Error]: Subgraph has more than one component. GraphB::list_subgraph_appearances_using_names.\n";
    exit(1);
  }
  if (nu==2)
    res = list_subgraph_appearances_using_names2(sg, os, me, ge, tis);
  else {
    if (nu == 3) {
      res = list_subgraph_appearances_using_names3(sg, os, me, ge, tis);
    }
    else {
      res = list_subgraph_appearances_using_names4(sg, os, me, ge, tis);
    }
  }
  return res;
}



void GraphB::count_subgraphs_2()
{
  if (!yaprecosg)
    prepare_to_count_subgraphs();
  int i, j,k;
  basic.fillv0(vedosg2, nupre_2);
  set<int> losque;
  GraphB redtemp(est);
  
  for (i=0; i<(number_of_nodes()-1); i++) {
    for (j=(i+1); j < (number_of_nodes()); j++) {
      losque.clear();
      losque.insert(i);
      losque.insert(j);
      put_subgraph_into(redtemp, losque);
      if (redtemp.get_components() == 1) {
        for (k=0; k<nupre_2; k++)
          if (LM_2[k].isomorphic_pair(redtemp)){
            vedosg2[k]++;
            break;
          }
      }
      redtemp.clear();
    }
  }
}

void GraphB::count_subgraphs_3()
{
  if (!yaprecosg)
    prepare_to_count_subgraphs();
  int i, j,k,l;
  basic.fillv0(vedosg3, nupre_3);
  set<int> losque;
  GraphB redtemp(est);
  for (i=0; i<(number_of_nodes()-2); i++) {
    for (j=(i+1); j < (number_of_nodes()-1); j++) {
      for (l=(j+1); l<number_of_nodes(); l++) {
        losque.clear();
        losque.insert(i);
        losque.insert(j);
        losque.insert(l);
        put_subgraph_into(redtemp, losque);
        if (redtemp.get_components() == 1) {
          for (k=0; k<nupre_3; k++)
            if (LM_3[k].isomorphic_pair(redtemp)){
              vedosg3[k]++;
              break;
            }
        }
        redtemp.clear();
      }
    }
  }
}

void GraphB::count_subgraphs_4()
{
  if (!yaprecosg)
    prepare_to_count_subgraphs();
  int i, j,k,l,m;
  if (directed) {
    cout << "[Error]: Subgraphs of size 4 for directed networks have not been generated. GraphB::count_subgraphs_4.\n";
    exit(1);
  }
  basic.fillv0(vedosg4, nupre_4);
  set<int> losque;
  GraphB redtemp(est);
  for (i=0; i<(number_of_nodes()-3); i++) {
    for (j=(i+1); j < (number_of_nodes()-2); j++) {
      for (l=(j+1); l<(number_of_nodes()-1); l++) {
        for (m=(l+1); m < number_of_nodes(); m++) {
          losque.clear();
          losque.insert(i);
          losque.insert(j);
          losque.insert(l);
          losque.insert(m);
          put_subgraph_into(redtemp, losque);
          if (redtemp.get_components() == 1) {
            for (k=0; k<nupre_4; k++)
              if (LM_4[k].isomorphic_pair(redtemp)){
                vedosg4[k]++;
                break;
              }
          }
          redtemp.clear();
        }
      }
    }
  }
}

void GraphB::count_subgraphs_2(bool **me, int ge, int tis)
{
  if (!yaprecosg)
    prepare_to_count_subgraphs();
  if (ge != number_of_nodes()) { //
    cout << "[Error]: Number of genes does not match the number of rows in expression matrix.\n";//
    exit(1);//
  }//
  bool shati; //
  int i, j,k, tt;//
  basic.fillv0(vedosg2, nupre_2);
  set<int> losque;
  GraphB redtemp(est);
  
  for (i=0; i<(number_of_nodes()-1); i++) {
    for (j=(i+1); j < (number_of_nodes()); j++) {
      losque.clear();
      losque.insert(i);
      losque.insert(j);
      shati = false; //
      for (tt=0; tt<tis; tt++) {//
        if (me[i][tt] && me[j][tt]) {//
          shati = true;//
          break;//
        }//
      }//
      if (shati) {//
        put_subgraph_into(redtemp, losque);
        if (redtemp.get_components() == 1) {
          for (k=0; k<nupre_2; k++)
            if (LM_2[k].isomorphic_pair(redtemp)){
              vedosg2[k]++;
              break;
            }
        }
        redtemp.clear();
      }		//
      
    }
  }
}

void GraphB::count_subgraphs_3(bool **me, int ge, int tis)
{
  if (!yaprecosg)
    prepare_to_count_subgraphs();
  if (ge != number_of_nodes()) { //
    cout << "[Error]: Number of genes does not match the number of rows in expression matrix.\n";//
    exit(1);//
  }//
  bool shati; //
  int i, j,k,l,tt;//
  basic.fillv0(vedosg3, nupre_3);
  set<int> losque;
  GraphB redtemp(est);
  for (i=0; i<(number_of_nodes()-2); i++) {
    for (j=(i+1); j < (number_of_nodes()-1); j++) {
      for (l=(j+1); l<number_of_nodes(); l++) {
        losque.clear();
        losque.insert(i);
        losque.insert(j);
        losque.insert(l);
        shati = false; //
        for (tt=0; tt<tis; tt++) {//
          if (me[i][tt] && me[j][tt] && me[l][tt]) {//
            shati = true;//
            break;//
          }//
        }//
        if (shati) { //
          put_subgraph_into(redtemp, losque);
          if (redtemp.get_components() == 1) {
            for (k=0; k<nupre_3; k++)
              if (LM_3[k].isomorphic_pair(redtemp)){
                vedosg3[k]++;
                break;
              }
          }
          redtemp.clear();
        }//
      }
    }
  }
}

void GraphB::count_subgraphs_4(bool **me, int ge, int tis)
{
  if (!yaprecosg)
    prepare_to_count_subgraphs();
  
  if (ge != number_of_nodes()) { //
    cout << "[Error]: Number of genes does not match the number of rows in expression matrix.\n";//
    exit(1);//
  }//
  bool shati; //
  int i, j,k,l,m,tt;//
  if (directed) {
    cout << "[Error]: Subgraphs of size 4 for directed networks have not been generated. GraphB::count_subgraphs_4.\n";
    exit(1);
  }
  basic.fillv0(vedosg4, nupre_4);
  set<int> losque;
  GraphB redtemp(est);
  for (i=0; i<(number_of_nodes()-3); i++) {
    for (j=(i+1); j < (number_of_nodes()-2); j++) {
      for (l=(j+1); l<(number_of_nodes()-1); l++) {
        for (m=(l+1); m < number_of_nodes(); m++) {
          losque.clear();
          losque.insert(i);
          losque.insert(j);
          losque.insert(l);
          losque.insert(m);
          shati = false; //
          for (tt=0; tt<tis; tt++) {//
            if (me[i][tt] && me[j][tt] && me[l][tt] && me[m][tt]) {//////aqui habia un error. habia l en lugar de m
              shati = true;//
              break;//
            }//
          }//
          if (shati) { //
            put_subgraph_into(redtemp, losque);
            if (redtemp.get_components() == 1) {
              for (k=0; k<nupre_4; k++)
                if (LM_4[k].isomorphic_pair(redtemp)){
                  vedosg4[k]++;
                  break;
                }
            }
            redtemp.clear();
          }//
        }
      }
    }
  }
}

int GraphB::how_many_subgraphs(int tam, int which)
{
  bool ya=false;
  int res=0, i;
  if (tam==2)
    res = vedosg2[which];
  else {
    if (tam==3)
      res = vedosg3[which];
    else {
      if (tam==4) {
        if (is_directed()) {
          cout << "[Error]: GraphB::how_many_subgraphs is not defined for directed subgraphs of size 4.\n";
          exit(1);
        }
        res = vedosg4[which];
      }
      else {
        cout << "[Error]: GraphB::how_many_subgraphs is not defined for subgraphs of size greater than 4.\n";
        exit(1);
      }
    }
  }
  if (res==0) {
    if (tam==2) {
      for (i=0; i<nupre_2; i++)
        if (vedosg2[i] > 0) {
          ya = true;
          break;
        }
    }
    if (tam==3) {
      for (i=0; i<nupre_3; i++)
        if (vedosg3[i] > 0) {
          ya = true;
          break;
        }
    }
    if (tam==4) {
      for (i=0; i<nupre_4; i++)
        if (vedosg4[i] > 0) {
          ya = true;
          break;
        }
    }
    if (!ya) {
      cout << "[Error]: Subgraphs must be counted first. GraphB::how_many_subgraphs.\n";
      exit(1);
    }
  }
  return res;
}

int GraphB::count_this_subgraph(int tam, int which)
{
  if (!yaprecosg)
    prepare_to_count_subgraphs();
  int res;
  if (tam==2)
    res = count_this_subgraph(LM_2[which]);
  else {
    if (tam==3)
      res = count_this_subgraph(LM_3[which]);
    else {
      if (tam==4)
        res = count_this_subgraph(LM_4[which]);
      else {
        cout << "[Error]: GraphB::count_this_subgraph is not defined for subgraphs of size greater than 4.\n";
        exit(1);
      }
    }
  }
  return res;
}

int GraphB::loose_count_this_subgraph(int tam, int which)
{
  if (!yaprecosg)
    prepare_to_count_subgraphs();
  int res;
  if (tam==2)
    res = loose_count_this_subgraph(LM_2[which]);
  else {
    if (tam==3)
      res = loose_count_this_subgraph(LM_3[which]);
    else {
      if (tam==4)
        res = loose_count_this_subgraph(LM_4[which]);
      else {
        cout << "[Error]: GraphB::loose_count_this_subgraph is not defined for subgraphs of size greater than 4.\n";
        exit(1);
      }
    }
  }
  return res;
}

set<int> GraphB::loose_list(GraphB &sg)
{
  if (is_directed() != sg.is_directed()) {
    cout << "[Error]: either the focal graph is directed and the subgraph not, or vice versa. GraphB::loose_list.\n";
    exit(1);
  }
  set<int> res;
  res.clear();
  if (!yaprecosg)
    prepare_to_count_subgraphs();
  int tam = sg.number_of_nodes();
  int k;
  if (tam==2) {
    for (k=0; k<nupre_2; k++)
      if (LM_2[k].an_isomorph_contains_these_interactions(sg))
        res.insert(k);
  }
  else {
    if (tam==3) {
      for (k=0; k<nupre_3; k++)
        if (LM_3[k].an_isomorph_contains_these_interactions(sg))
          res.insert(k);
    }
    else {
      if (tam==4) {
        for (k=0; k<nupre_4; k++)
          if (LM_4[k].an_isomorph_contains_these_interactions(sg))
            res.insert(k);
      }
      else {
        cout << "[Error]: GraphB::loose_list does not work for subgraphs of size greater than 4.\n";
        exit(1);
      }
    }
  }
  return res;
}

int GraphB::number_of_sgs(int nnodes, int cuasg)
{
  int res = 0;
  if (!yaprecosg) {
    cout << "[Error]: subgraph count must be performed before calling GraphB::number_of_sgs.\n";
    exit(1);
  }
  if ((nnodes!=2)&&(nnodes!=3)&&(nnodes!=4)) {
    cout << "[Error]: subgraphs of size other than 2, 3 or 4 can not be counted. GraphB::number_of_sgs.\n";
    exit(1);
  }
  if (cuasg<0) {
    cout << "[Error]: subgraph id number must be equal or greater than zero. GraphB::number_of_sgs.\n";
    exit(1);
  }
  
  if (nnodes==2) {
    if (cuasg >= nupre_2) {
      cout << "[Error]: No subgraph of size " << nnodes << " with id " << cuasg << ". GraphB::number_of_sgs.\n";
      exit(1);
    }
    else
      res = vedosg2[cuasg];
  }
  if (nnodes==3) {
    if (cuasg >= nupre_3) {
      cout << "[Error]: No subgraph of size " << nnodes << " with id " << cuasg << ". GraphB::number_of_sgs.\n";
      exit(1);
    }
    else
      res = vedosg3[cuasg];
  }
  if (nnodes==4) {
    if (cuasg >= nupre_4) {
      cout << "[Error]: No subgraph of size " << nnodes << " with id " << cuasg << ". GraphB::number_of_sgs.\n";
      exit(1);
    }
    else
      res = vedosg4[cuasg];
  }
  return res;
}

//modularity

void GraphB::build_moma_u(double **momaundi) {
  if (!yamadya)
    get_adjacency_matrix();
  if (!directed) {
    if (!yadeg)
      get_all_degrees();
    build_moma_u_aux(matadya, degree, number_of_edges(), momaundi);
  }
  else {
    int **madjundi, *degundi, edchundi = 0, i,j;
    basic.create_array(madjundi, size, size);
    degundi = new int[size];
    basic.fillv0(degundi, size);
    basic.fillmat0(madjundi, size, size);
    for (i = 0; i < size; i++)
      for (j = 0; j < size; j++)
        if (matadya[i][j] != 0) {
          madjundi[i][j]++;
          madjundi[j][i]++;
          degundi[i]++;
          degundi[j]++;
          edchundi++;
        }
    //Control!
    if (edchundi != number_of_edges()) {
      cout << "[Error]: Mira, no sale igual!\n";
      exit(1);
    }
    build_moma_u_aux(madjundi, degundi, edchundi, momaundi);
    for (i = 0; i < size; i++)
      delete [] madjundi[i];
    delete [] madjundi;
    delete [] degundi;
  }
}

void GraphB::build_moma_d(double **momadi) {
  //Momadi corresponds to the final modularity matrix (B + B^T according to Leicht and Newman 2008)
  if (!directed) {
    cout << "[Error]: Attempt to build directed modularity matrix from undirected adjacency matrix in GraphB::build_moma_d\n";
    exit(1);
  }
  if (!yamadya)
    get_adjacency_matrix();
  if (!yadeg)
    get_all_degrees();
  int i, j;
  for (i = 0; i < size; i++)
    for (j= 0; j < size; j++)
      momadi[i][j] = matadya[i][j] + matadya[j][i] - ((indegree[i]*outdegree[j]) + (indegree[j]*outdegree[i]))/(double(number_of_edges()));
}

void GraphB::build_moma_out(double **momaout) {
  if (!directed) {
    cout << "[Error]: Attempt to build directed modularity matrix from undirected adjacency matrix in GraphB::build_moma_out\n";
    exit(1);
  }
  if (!yamadya)
    get_adjacency_matrix();
  if (!yadeg)
    get_all_degrees();
  double deno1= 0, deno2 = 0;
  int **c;
  int i, j, k;
  basic.create_array(c, size, size);
  for (i = 0; i < size; i++) {
    for (j = 0; j < size; j++) {
      if (i == j)
        c[i][j] = outdegree[i];
      else {
        c[i][j] = 0;
        for (k = 0; k < size; k++)
          if ((matadya[k][i] > 0) && matadya[k][j] > 0)
            c[i][j] += 1.0;
      }
    }
  }
  deno2 = double(number_of_edges()*number_of_edges());
  deno1 = 0;
  for (i = 0; i < size; i++)
    if (indegree[i] > 0)
      deno1 += (indegree[i]*(indegree[i] - 1.0));
  for (i = 0; i < size; i++)
    for (j = 0; j < size; j++)
      momaout[i][j] = (c[i][j]/deno1) - ((outdegree[i]*outdegree[j])/deno2);
  for (i = 0; i < size; i++)
    delete [] c[i];
  delete [] c;
}

void GraphB::build_moma_in(double **momain) {
  if (!directed) {
    cout << "[Error]: Attempt to build directed modularity matrix from undirected adjacency matrix in GraphB::build_moma_in\n";
    exit(1);
  }
  if (!yamadya)
    get_adjacency_matrix();
  if (!yadeg)
    get_all_degrees();
  double deno1= 0, deno2 = 0;
  int **c;
  int i, j, k;
  basic.create_array(c, size, size);
  for (i = 0; i < size; i++) {
    for (j = 0; j < size; j++) {
      if (i == j)
        c[i][j] = indegree[i];
      else {
        c[i][j] = 0;
        for (k = 0; k < size; k++)
          if ((matadya[i][k] > 0) && matadya[j][k] > 0)//
            c[i][j] += 1.0;
      }
    }
  }
  deno2 = double(number_of_edges()*number_of_edges());
  deno1 = 0;
  for (i = 0; i < size; i++)
    if (outdegree[i] > 0)
      deno1 += (outdegree[i]*(outdegree[i] - 1.0));
  
  for (i = 0; i < size; i++)
    for (j = 0; j < size; j++)
      momain[i][j] = (c[i][j]/deno1) - ((indegree[i]*indegree[j])/deno2);
  for (i = 0; i < size; i++)
    delete [] c[i];
  delete [] c;
}

int GraphB::components_gt1() {
  if (!yacomp)
    get_components();
  int i, j = 0;
  for (i = 0; i < number_of_components(); i++)
    if (basic.count_in_vector(components, size, i) == 1)
      j++;
  return number_of_components() - j;
}

int GraphB::sccs_gt1() {
  if (!cfc)
    get_scc();
  int i, j = 0;
  for (i = 0; i < number_of_sccs(); i++)
    if (basic.count_in_vector(scc, size, i) == 1)
      j++;
  return number_of_sccs() - j;
}

void GraphB::check_props_of_rdset(int samplesiz, int minsw, int maxsw, double **meansdminmax) {
  //  hilera 0: comp
  //  hilera 1: comp gt1
  //  2: scc
  //  3: scc gt1
  //  4: modund
  //  5: moddir
  //  6: moddirmax
  //  7:modin
  //  8:modout
  double **todo;
  int e = number_of_edges();
  basic.create_array(todo, 9, samplesiz);
  int i, j;
  GraphB atiza(est);
  make_copy(atiza);
  
  double **matmod;
  basic.create_array(matmod, size, size);
  set<set<int> > teams;
  for (i = 0; i < samplesiz; i++) {
    j = est.randint(minsw, maxsw+1);
    j *= e;
    atiza.switches(j);
    todo[0][i] = atiza.number_of_components();
    todo[1][i] = atiza.components_gt1();
    todo[2][i] = atiza.number_of_sccs();
    todo[3][i] = atiza.sccs_gt1();
    teams.clear();
    atiza.build_moma_u(matmod);
    atiza.iterative_newman06(matmod, teams);
    todo[4][i] = atiza.eval_mod(matmod, teams);
    teams.clear();
    atiza.build_moma_d(matmod);
    atiza.iterative_newman06(matmod, teams);
    todo[5][i] = atiza.eval_mod(matmod, teams);
    todo[6][i] = atiza.maxmod_d_per_part_and_dd(teams);
    teams.clear();
    atiza.build_moma_in(matmod);
    atiza.iterative_newman06(matmod, teams);
    todo[7][i] = atiza.eval_mod(matmod, teams);
    teams.clear();
    atiza.build_moma_out(matmod);
    atiza.iterative_newman06(matmod, teams);
    todo[8][i] = atiza.eval_mod(matmod, teams);
  }
  atiza.clear();
  for (i = 0; i < size; i++)
    delete [] matmod[i];
  delete [] matmod;
  teams.clear();
  for (i = 0; i < 9; i++) {
    meansdminmax[i][0] = basic.get_mean(todo[i], samplesiz);
    meansdminmax[i][1] = basic.get_sample_stddev(todo[i], samplesiz);
  }
  for (i = 0; i < 9; i++)
    delete [] todo[i];
  delete [] todo;
}

void GraphB::check_props_of_rdset(int samplesiz, int minsw, int maxsw, double **meansdminmax, set<set<int> > &predpar) {
  //  hilera 0: comp
  //  hilera 1: comp gt1
  //  2: scc
  //  3: scc gt1
  //  4: modund
  //  5: moddir
  //  6: moddirmax
  //  7:modin
  //  8:modout
  double **todo;
  int e = number_of_edges();
  basic.create_array(todo, 10, samplesiz);
  int i, j;
  GraphB atiza(est);
  make_copy(atiza);
  
  double **matmod;
  basic.create_array(matmod, size, size);
  set<set<int> > teams;
  for (i = 0; i < samplesiz; i++) {
    j = est.randint(minsw, maxsw+1);
    j *= e;
    atiza.switches(j);
    todo[0][i] = atiza.number_of_components();
    todo[1][i] = atiza.components_gt1();
    todo[2][i] = atiza.number_of_sccs();
    todo[3][i] = atiza.sccs_gt1();
    teams.clear();
    atiza.build_moma_u(matmod);
    atiza.iterative_newman06(matmod, teams);
    todo[4][i] = atiza.eval_mod(matmod, teams);
    teams.clear();
    atiza.build_moma_d(matmod);
    atiza.iterative_newman06(matmod, teams);
    todo[5][i] = atiza.eval_mod(matmod, teams);
    todo[6][i] = atiza.maxmod_d_per_part_and_dd(teams);
    todo[9][i] = atiza.eval_mod(matmod, predpar);
    teams.clear();
    atiza.build_moma_in(matmod);
    atiza.iterative_newman06(matmod, teams);
    todo[7][i] = atiza.eval_mod(matmod, teams);
    teams.clear();
    atiza.build_moma_out(matmod);
    atiza.iterative_newman06(matmod, teams);
    todo[8][i] = atiza.eval_mod(matmod, teams);
  }
  atiza.clear();
  for (i = 0; i < size; i++)
    delete [] matmod[i];
  delete [] matmod;
  teams.clear();
  for (i = 0; i < 10; i++) {
    meansdminmax[i][0] = basic.get_mean(todo[i], samplesiz);
    meansdminmax[i][1] = basic.get_sample_stddev(todo[i], samplesiz);
  }
  for (i = 0; i < 10; i++)
    delete [] todo[i];
  delete [] todo;
}

void GraphB::check_Qd_in_rdset(int samplesiz, int minsw, int maxsw, double *meansdminmax) {
  //  hilera 0: moddir
  double *todo;
  int e = number_of_edges();
  todo = new double[samplesiz];
  int i, j;
  GraphB atiza(est);
  make_copy(atiza);
  
  double **matmod;
  basic.create_array(matmod, size, size);
  set<set<int> > teams;
  for (i = 0; i < samplesiz; i++) {
    j = est.randint(minsw, maxsw+1);
    j *= e;
    atiza.switches(j);
    teams.clear();
    atiza.build_moma_d(matmod);
    atiza.iterative_newman06(matmod, teams);
    todo[i] = atiza.eval_mod(matmod, teams);
  }
  atiza.clear();
  for (i = 0; i < size; i++)
    delete [] matmod[i];
  delete [] matmod;
  teams.clear();
  meansdminmax[0] = basic.get_mean(todo, samplesiz);
  meansdminmax[1] = basic.get_sample_stddev(todo, samplesiz);
  delete [] todo;
}

void GraphB::check_Qd_in_rdset(int samplesiz, int minsw, int maxsw, double **meansdminmax, set<set<int> > &predpar) {
  //  hilera 0: moddir
  //hilera1: moddirpred
  double **todo;
  int e = number_of_edges();
  basic.create_array(todo, 2, samplesiz);
  int i, j;
  GraphB atiza(est);
  make_copy(atiza);
  
  double **matmod;
  basic.create_array(matmod, size, size);
  set<set<int> > teams;
  for (i = 0; i < samplesiz; i++) {
    j = est.randint(minsw, maxsw+1);
    j *= e;
    atiza.switches(j);
    teams.clear();
    atiza.build_moma_d(matmod);
    atiza.iterative_newman06(matmod, teams);
    todo[0][i] = atiza.eval_mod(matmod, teams);
    todo[1][i] = atiza.eval_mod(matmod, predpar);
  }
  atiza.clear();
  for (i = 0; i < size; i++)
    delete [] matmod[i];
  delete [] matmod;
  teams.clear();
  for (i = 0; i < 2; i++) {
    meansdminmax[i][0] = basic.get_mean(todo[i], samplesiz);
    meansdminmax[i][1] = basic.get_sample_stddev(todo[i], samplesiz);
  }
  for (i = 0; i < 2; i++)
    delete [] todo[i];
  delete [] todo;
}

void GraphB::build_moma_u_aux(int **adjundi, int *degundi, int edgundim, double **momaundi) {
  int i,j;
  for (i=0; i < size; i++)
    for (j = 0; j < size; j++)
      momaundi[i][j] = adjundi[i][j] - ((degundi[i]*degundi[j])/(2.0*edgundim));
}

void GraphB::mod_partitions_from_vector(double *vecs, int tam, int *diccion, set<int> &parta, set<int> &partb) {
  
  parta.clear();
  partb.clear();
  int i;
  for (i=0; i < tam; i++) {
    if (vecs[i] == 1)
      parta.insert(diccion[i]);
    else if (vecs[i] == -1)
      partb.insert(diccion[i]);
    else {
      cout << "[Error]: Vector contains numbers other than 1 or -1 in GraphB::mod_partitions_from_vector.\n";
      exit(1);
    }
  }
}

double GraphB::spectral_method(double **moduloc, int tam, double *vecs) {
  if (!basic.is_symmetric(moduloc, tam, tam)) {
    cout << "[Error]: Non-symmetric matrix in GrahB::spectral_method.\n";
    exit(1);
  }
  double eval = basic.get_leading_evector_symmat(moduloc, tam, vecs);
  return eval;
}

double GraphB::mod_after_spectral(double **matloc, int tam, double *vecs) {
  int i;
  for (i = 0; i < tam; i++) {
    if (vecs[i] > 0)
      vecs[i] = 1;
    if (vecs[i] < 0)
      vecs[i] = -1;
    if (vecs[i] == 0) {
      if (est.toss())
        vecs[i] = 1;
      else
        vecs[i] = -1;
    }
  }
  double mienQ = eval_mod_subm(matloc, vecs, tam);
  return mienQ;
}

void GraphB::shake_kl(double mienQ, double **matloc, int tam, double *vecs) {
  double mamol, miemol;
  int i, besti, cont = 0;
  mamol = mienQ;
  do {
    besti = -1;
    for (i = 0; i < tam; i++) {
      vecs[i] *= -1;
      miemol = eval_mod_subm(matloc, vecs, tam);
      if (miemol > mamol) {
        besti = i;
        mamol = miemol;
      }
      vecs[i] *= -1;
    }
    if (besti >= 0)
      vecs[besti] *= -1;
    cont++;
    //    if (cont >= (10000*tam)) {
    //      cout << tam << " muchas vueltas\n";
    //      exit(1);
    //    }
  } while (besti >= 0);
  i = cont;
}

double GraphB::eval_mod_subm(double **mamod, double *vecs, int tam) {
  double *vt;
  vt = new double[tam];
  basic.matxvec(vecs, tam, tam , tam, mamod, vt);
  double Q = basic.dotproduct(vecs, vt, tam);
  delete [] vt;
  Q /= (4.0*number_of_edges());
  return Q;
}

void GraphB::partition_to_vector(const set<set<int> > &equipos, int *vop) {
  set<set<int> >::iterator it;
  set<int>::iterator ite;
  set<int> whi;
  basic.fillvm1(vop, size);
  int conc = 0;
  for (it = equipos.begin(); it != equipos.end(); it++) {
    whi = *it;
    for (ite = whi.begin(); ite != whi.end(); ite++)
      vop[*ite] = conc;
    conc++;
  }
}

double GraphB::eval_mod(double **mamod, int *vop) {
  int i, j;
  double Q = 0;
  for (i = 0; i < size; i++)
    for (j = 0; j < size; j++)
      if (vop[i] == vop[j])
        Q += mamod[i][j];
  Q /= (2.0*number_of_edges());
  return Q;
}

double GraphB::eval_mod(double **mamod, set<set<int> > &equipos) {
  int *arr;
  arr = new int[size];
  partition_to_vector(equipos, arr);
  double Q = eval_mod(mamod, arr);
  delete [] arr;
  return Q;
}

double GraphB::maxmod_u_per_part_and_dd(int *vop) {
  double Qmax = 2.0*number_of_edges();
  if (!yadeg)
    get_all_degrees();
  int i, j;
  for (i = 0; i < size; i++)
    for (j = 0; j < size; j++)
      if ((vop[i] >= 0) &&(vop[i] == vop[j]))
        Qmax -= (degree[i]*degree[j])/(2.0*number_of_edges());
  Qmax /= 2.0*number_of_edges();
  return Qmax;
}

double GraphB::maxmod_u_per_part_and_dd(set<set<int> > &equipos) {
  int *arr;
  arr = new int[size];
  partition_to_vector(equipos, arr);
  double Qmax = maxmod_u_per_part_and_dd(arr);
  delete [] arr;
  return Qmax;
}
//
double GraphB::maxmod_d_per_part_and_dd(int *vop) {
  double Qmax = 2.0*number_of_edges();
  if (!yadeg)
    get_all_degrees();
  int i, j;
  for (i = 0; i < size; i++)
    for (j = 0; j < size; j++)
      if ((vop[i] >= 0) &&(vop[i] == vop[j]))
        Qmax -= (indegree[i]*outdegree[j])/(2.0*number_of_edges());
  Qmax /= 2.0*number_of_edges();
  return Qmax;
}

double GraphB::maxmod_d_per_part_and_dd(set<set<int> > &equipos) {
  int *arr;
  arr = new int[size];
  partition_to_vector(equipos, arr);
  double Qmax = maxmod_d_per_part_and_dd(arr);
  delete [] arr;
  return Qmax;
}

void GraphB::adjust_modmat(double **mator, set<int> &parta, int tam, int *diccion, double **matloc) {
  set<int>::iterator It;
  bool *arte;
  arte = new bool[size];
  basic.fillv0(arte, size);
  for (It = parta.begin(); It != parta.end(); It++)
    arte[*It] = true;
  int i, j=0;
  for (i = 0; i < size; i++)
    if (arte[i]) {
      diccion[j] = i;
      j++;
    }
  double aver;
  delete [] arte;
  for (i = 0; i < tam; i++)
    for (j = 0; j < tam; j++)
      matloc[i][j] = mator[diccion[i]][diccion[j]];
  for (i = 0; i < tam; i++) {
    aver = matloc[i][i];
    for (j = 0; j < tam; j++)
      aver -= matloc[i][j];
    matloc[i][i] = aver;
  }
}

void GraphB::iterative_newman06(double **mator, set<set<int> > &res) {
  int i;
  int cuantu;
  set<int> todos;
  for (i = 0; i < size; i++)
    todos.insert(i);
  double mQ, *vecs, evalue;
  bool otra = false;
  vecs  = new double[size];
  evalue = spectral_method(mator, size, vecs);
  if (evalue > 0) {
    mQ = mod_after_spectral(mator, size, vecs);
    if (mQ > 0) {
      cuantu = basic.count_in_vector(vecs, size, 1);
      if ((cuantu != 0) && (cuantu != size)) {
        shake_kl(mQ, mator, size, vecs);
        cuantu = basic.count_in_vector(vecs, size, 1);
        if ((cuantu != 0) && (cuantu != size))
          otra = true;
      }
    }
  }
  if (otra) {
    todos.clear();
    set<int> parta1;
    set<int> parta2;
    parta1.clear();
    parta2.clear();
    for (i = 0; i < size; i++) {
      if (vecs[i] == 1)
        parta1.insert(i);
      else if (vecs[i] == -1)
        parta2.insert(i);
      else {
        cout << "[Error]: Vector contains numbers other than 1 or -1 in GraphB::mod_partitions_from_vector.\n";
        exit(1);
      }
    }
    set<set<int> > res2;
    delete [] vecs;
    iterative_newman06(mator, parta1, res);
    iterative_newman06(mator, parta2, res2);
    basic.append_sets_of_sets(res2, res);
    res2.clear();
  }
  else {
    delete [] vecs;
    res.insert(todos);
    todos.clear();
  }
}

void GraphB::iterative_newman06(double **mator, set<int> &parta, set<set<int> > &res) {
  int i,*diccion, tam = parta.size();
  bool otra = false;
  int cuantu;
  
  double mQ, **matloc, *vecs, evalue;
  diccion = new int[tam];
  basic.create_array(matloc, tam, tam);
  vecs = new double[tam];
  adjust_modmat(mator, parta, tam, diccion, matloc);
  evalue = spectral_method(matloc, tam, vecs);
  if (evalue > 0) {
    mQ = mod_after_spectral(matloc, tam, vecs);
    if (mQ > 0) {
      cuantu = basic.count_in_vector(vecs, tam, 1);
      if ((cuantu != 0) && (cuantu != tam)) {
        shake_kl(mQ, matloc, tam, vecs);
        cuantu = basic.count_in_vector(vecs, tam, 1);
        if ((cuantu != 0) && (cuantu != tam))
          otra = true;
      }
    }
  }
  for (i = 0; i < tam; i++)
    delete [] matloc[i];
  delete [] matloc;
  
  if (otra) {
    parta.clear();
    set<int> parta1;
    set<int> parta2;
    mod_partitions_from_vector(vecs, tam, diccion, parta1, parta2);
    set<set<int> > res2;
    delete [] vecs;
    delete [] diccion;
    iterative_newman06(mator, parta1, res);
    iterative_newman06(mator, parta2, res2);
    basic.append_sets_of_sets(res2, res);
    res2.clear();
    
  }
  else {
    delete [] vecs;
    delete [] diccion;
    res.insert(parta);
    parta.clear();
  }
}

void GraphB::split_in_2mods(double **mator, set<set<int> > &res) {
  int i;
  set<int> todos;
  for (i = 0; i < size; i++)
    todos.insert(i);
  double mQ, *vecs, evalue;
  bool otra = false;
  vecs = new double[size];
  evalue = spectral_method(mator, size, vecs);
  if (evalue > 0) {
    mQ = mod_after_spectral(mator, size, vecs);
    if (mQ > 0) {
      otra = true;
      shake_kl(mQ, mator, size, vecs);
    }
  }
  if (otra) {
    todos.clear();
    set<int> parta1;
    set<int> parta2;
    parta1.clear();
    parta2.clear();
    for (i = 0; i < size; i++) {
      if (vecs[i] == 1)
        parta1.insert(i);
      else if (vecs[i] == -1)
        parta2.insert(i);
      else {
        cout << "[Error]: Vector contains numbers other than 1 or -1 in GraphB::mod_partitions_from_vector.\n";
        exit(1);
      }
    }
    res.insert(parta1);
    res.insert(parta2);
    parta1.clear();
    parta2.clear();
    delete [] vecs;
  }
  else {
    delete [] vecs;
    res.insert(todos);
    todos.clear();
  }
}

//private

void GraphB::set_default_exclusive_vars()
{
  infl = false;
  dist = false;
  yacomp = false;
  cfc = false;
  bcya = false;
  floops = false;
  yae = false;
  yadeg = false;
  yamotdab=false;
  yasgcat=false;
  yaprecosg=false;
  yadegdist=false;
  yamadya = false;
}


int GraphB::calc_indegree(int no)
{
  int id = 0;
  int i;
  for (i = 0; i < size; i++)
    if (nw[no][i])
      id++;
  return id;
}

int GraphB::calc_outdegree(int no)
{
  int od = 0;
  int i;
  for (i = 0; i < size; i++)
    if (nw[i][no])
      od++;
  return od;
}

int GraphB::calc_degree(int no)
{
  int d;
  if (directed)
    d = indegree[no] + outdegree[no];
  else {
    if (indegree[no] != outdegree[no]) {
      cout << "[Error]: Different indegree and outdegree in an undirected network. GraphB::calc_degree.\n";
      exit(1);
    }
    d = indegree[no];
    if (weight(no,no))
      d++;
  }
  return d;
}

void GraphB::prepare_for_degrees()
{
  degree = new int[size];
  outdegree = new int[size];
  indegree = new int[size];
}

void GraphB::prepare_for_degdist()
{
  degdist = new int[size+1];
  basic.fillv0(degdist, size+1);
  
  odegdist = new int[size+1];
  basic.fillv0(odegdist, size+1);
  
  idegdist = new int[size+1];
  basic.fillv0(idegdist, size+1);
}

void GraphB::priv_get_degdist()
{
  if (!yadeg)
    get_all_degrees();
  int i;
  for (i=0; i<size; i++) {
    degdist[degree[i]]++;
    odegdist[outdegree[i]]++;
    idegdist[indegree[i]]++;
  }
}

void GraphB::aswitch()
{
  int i, j, k, l, cont=0;
  bool v1, v2, v3, v4; //
  int nesq = number_of_edges()*number_of_edges();
  bool swdone = false;
  do {
    do {
      i = est.randint(0, size);
      j = est.randint(0, size);
    }while (!nw[i][j]);
    do {
      k = est.randint(0, size);
      l = est.randint(0, size);
    }while (!nw[k][l]);
    if (directed) { //
      if (!((nw[i][l]) || (nw[k][j]) || (i==k) || (j==l) )) {
        v1 = nw[i][j];
        v2 = nw[k][l];
        nw[i][j] = false;
        nw[k][l] = false;
        nw[i][l] = v1;
        nw[k][j] = v2;
        swdone = true; //
      }
    } //////
    else {
      if (!((nw[i][l]) || (nw[k][j]) || (i==k) || (j==l) || ((i==j) && (k==l)) || ((i==l) && (k==j)))) {
        v1 = nw[i][j];
        v3 = nw[j][i];
        v2 = nw[k][l];
        v4 = nw[l][k];
        nw[i][j] = false;
        nw[k][l] = false;
        nw[j][i] = false;
        nw[l][k] = false;
        nw[i][l] = v1;
        nw[k][j] = v2;
        nw[l][i] = v3;
        nw[j][k] = v4;
        swdone = true;
      }
    }
    cont++;
  }while ((!swdone) && (cont < (2*nesq)));
  if (!swdone) {
    cout << "[Error]: It is hard to find switchable edges. GraphB::aswitch.\n";
    exit(1);
  }
}

void GraphB::aswitch_preserving_sg1()
{
  int i, j, k, l, cont=0;
  bool v1, v2, v3, v4; //
  int nesq = number_of_edges()*number_of_edges();
  bool swdone = false;
  do {
    do {
      i = est.randint(0, size);
      j = est.randint(0, size);
    }while (!nw[i][j]);
    do {
      k = est.randint(0, size);
      l = est.randint(0, size);
    }while (!nw[k][l]);
    if (directed) { //																	//
      if (!((nw[i][l]) || (nw[k][j]) || (i==k) || (j==l)  || (i==j) || (k==l) || (i==l) || (k==j)   )) {
        v1 = nw[i][j];
        v2 = nw[k][l];
        nw[i][j] = false;
        nw[k][l] = false;
        nw[i][l] = v1;
        nw[k][j] = v2;
        swdone = true; //
      }
    } //////
    else {																																													 //
      if (!((nw[i][l]) || (nw[k][j]) || (i==k) || (j==l) || ((i==j) && (k==l)) || ((i==l) && (k==j))   || (i==j) || (k==l) || (i==l) || (k==j)  )) {
        v1 = nw[i][j];
        v3 = nw[j][i];
        v2 = nw[k][l];
        v4 = nw[l][k];
        nw[i][j] = false;
        nw[k][l] = false;
        nw[j][i] = false;
        nw[l][k] = false;
        nw[i][l] = v1;
        nw[k][j] = v2;
        nw[l][i] = v3;
        nw[j][k] = v4;
        swdone = true;
      }
    }
    cont++;
  }while ((!swdone) && (cont < (2*nesq)));
  if (!swdone) {
    cout << "[Error]: It is hard to find switchable edges. GraphB::aswitch_preserving_sg1.\n";
    exit(1);
  }
}

void GraphB::aswitch_preserving_sg2()
{
  int i, j, k, l, cont=0;
  bool v1, v2, v3, v4; //
  int nesq = number_of_edges()*number_of_edges();
  bool swdone = false;
  do {
    do {
      i = est.randint(0, size);
      j = est.randint(0, size);
    }while (!nw[i][j]);
    do {
      k = est.randint(0, size);
      l = est.randint(0, size);
    }while (!nw[k][l]);
    if (directed) { //																	//																			 //
      if (!((nw[i][l]) || (nw[k][j]) || (i==k) || (j==l) || (i==j) || (k==l) || (i==l) || (k==j) || (nw[i][i]!=nw[k][k]) || (nw[j][j]!=nw[l][l]) )) {
        v1 = nw[i][j];
        v2 = nw[k][l];
        nw[i][j] = false;
        nw[k][l] = false;
        nw[i][l] = v1;
        nw[k][j] = v2;
        swdone = true; //
      }
    } //////
    else {																																													 //																				//
      if (!((nw[i][l]) || (nw[k][j]) || (i==k) || (j==l) || ((i==j) && (k==l)) || ((i==l) && (k==j))   || (i==j) || (k==l) || (i==l) || (k==j) || (nw[i][i]!=nw[k][k]) || (nw[j][j]!=nw[l][l]) )) {
        v1 = nw[i][j];
        v3 = nw[j][i];
        v2 = nw[k][l];
        v4 = nw[l][k];
        nw[i][j] = false;
        nw[k][l] = false;
        nw[j][i] = false;
        nw[l][k] = false;
        nw[i][l] = v1;
        nw[k][j] = v2;
        nw[l][i] = v3;
        nw[j][k] = v4;
        swdone = true;
      }
    }
    cont++;
  }while ((!swdone) && (cont < (2*nesq)));
  if (!swdone) {
    cout << "[Error]: It is hard to find switchable edges. GraphB::aswitch_preserving_sg2.\n";
    exit(1);
  }
}

//exclusive of graphb
int GraphB::count_this_subgraph2(GraphB &sg)
{
  int i, j, res=0;
  set<int> losque;
  GraphB redtemp(est);
  if (sg.is_directed() != is_directed()) {
    cout << "[Error]: Either the chosen subgraph or the focal graph is directed while the other is undirected. GraphB::count_this_subgraph2.\n";
    exit(1);
  }
  for (i=0; i<(number_of_nodes()-1); i++) {
    for (j=(i+1); j < (number_of_nodes()); j++) {
      losque.clear();
      losque.insert(i);
      losque.insert(j);
      put_subgraph_into(redtemp, losque);
      if (redtemp.get_components() == 1)
        if (redtemp.isomorphic_pair(sg))
          res++;
      redtemp.clear();
    }
  }
  return res;
}

int GraphB::count_this_subgraph3(GraphB &sg)
{
  if (sg.is_directed() != is_directed()) {
    cout << "[Error]: Either the chosen subgraph or the focal graph is directed while the other is undirected. GraphB::count_this_subgraph3.\n";
    exit(1);
  }
  int i, j,l, res=0;
  set<int> losque;
  GraphB redtemp(est);
  for (i=0; i<(number_of_nodes()-2); i++) {
    for (j=(i+1); j < (number_of_nodes()-1); j++) {
      for (l=(j+1); l<number_of_nodes(); l++) {
        losque.clear();
        losque.insert(i);
        losque.insert(j);
        losque.insert(l);
        put_subgraph_into(redtemp, losque);
        if (redtemp.get_components() == 1)
          if (redtemp.isomorphic_pair(sg))
            res++;
        redtemp.clear();
      }
    }
  }
  return res;
}

int GraphB::count_this_subgraph4(GraphB &sg)
{
  if (sg.is_directed() != is_directed()) {
    cout << "[Error]: Either the chosen subgraph or the focal graph is directed while the other is undirected. GraphB::count_this_subgraph4.\n";
    exit(1);
  }
  int i, j,l,m,res=0;
  set<int> losque;
  GraphB redtemp(est);
  for (i=0; i<(number_of_nodes()-3); i++) {
    for (j=(i+1); j < (number_of_nodes()-2); j++) {
      for (l=(j+1); l<(number_of_nodes()-1); l++) {
        for (m=(l+1); m < number_of_nodes(); m++) {
          losque.clear();
          losque.insert(i);
          losque.insert(j);
          losque.insert(l);
          losque.insert(m);
          put_subgraph_into(redtemp, losque);
          if (redtemp.get_components() == 1)
            if (redtemp.isomorphic_pair(sg))
              res++;
          redtemp.clear();
        }
      }
    }
  }
  return res;
}

int GraphB::loose_count_this_subgraph2(GraphB &sg)
{
  int i, j, res=0;
  set<int> losque;
  GraphB redtemp(est);
  if (sg.is_directed() != is_directed()) {
    cout << "[Error]: Either the chosen subgraph or the focal graph is directed while the other is undirected. GraphB::loose_count_this_subgraph2.\n";
    exit(1);
  }
  for (i=0; i<(number_of_nodes()-1); i++) {
    for (j=(i+1); j < (number_of_nodes()); j++) {
      losque.clear();
      losque.insert(i);
      losque.insert(j);
      put_subgraph_into(redtemp, losque);
      if (redtemp.get_components() == 1)
        if (redtemp.an_isomorph_contains_these_interactions(sg)) ///?????
          res++;
      redtemp.clear();
    }
  }
  return res;
}

int GraphB::loose_count_this_subgraph3(GraphB &sg)
{
  if (sg.is_directed() != is_directed()) {
    cout << "[Error]: Either the chosen subgraph or the focal graph is directed while the other is undirected. GraphB::loose_count_this_subgraph3.\n";
    exit(1);
  }
  int i, j,l, res=0;
  set<int> losque;
  GraphB redtemp(est);
  for (i=0; i<(number_of_nodes()-2); i++) {
    for (j=(i+1); j < (number_of_nodes()-1); j++) {
      for (l=(j+1); l<number_of_nodes(); l++) {
        losque.clear();
        losque.insert(i);
        losque.insert(j);
        losque.insert(l);
        put_subgraph_into(redtemp, losque);
        if (redtemp.get_components() == 1)
          if (redtemp.an_isomorph_contains_these_interactions(sg))
            res++;
        redtemp.clear();
      }
    }
  }
  return res;
}

int GraphB::loose_count_this_subgraph4(GraphB &sg)
{
  if (sg.is_directed() != is_directed()) {
    cout << "[Error]: Either the chosen subgraph or the focal graph is directed while the other is undirected. GraphB::loose_count_this_subgraph4.\n";
    exit(1);
  }
  int i, j,l,m,res=0;
  set<int> losque;
  GraphB redtemp(est);
  for (i=0; i<(number_of_nodes()-3); i++) {
    for (j=(i+1); j < (number_of_nodes()-2); j++) {
      for (l=(j+1); l<(number_of_nodes()-1); l++) {
        for (m=(l+1); m < number_of_nodes(); m++) {
          losque.clear();
          losque.insert(i);
          losque.insert(j);
          losque.insert(l);
          losque.insert(m);
          put_subgraph_into(redtemp, losque);
          if (redtemp.get_components() == 1)
            if (redtemp.an_isomorph_contains_these_interactions(sg))
              res++;
          redtemp.clear();
        }
      }
    }
  }
  return res;
}

int GraphB::list_subgraph_appearances2(GraphB &sg, string arch)
{
  int res;
  ofstream sal;
  basic.open_ofstream(sal, arch);
  res = list_subgraph_appearances2(sg, sal);
  sal.close();
  return res;
}

int GraphB::list_subgraph_appearances3(GraphB &sg, string arch)
{
  int res;
  ofstream sal;
  basic.open_ofstream(sal, arch);
  res = list_subgraph_appearances3(sg, sal);
  sal.close();
  return res;
}

int GraphB::list_subgraph_appearances4(GraphB &sg, string arch)
{
  int res;
  ofstream sal;
  basic.open_ofstream(sal, arch);
  res = list_subgraph_appearances4(sg, sal);
  sal.close();
  return res;
}

int GraphB::list_subgraph_appearances2(GraphB &sg, ostream& os)
{
  int *ord;
  ord = new int[2];
  int i, j, res=0,ai,ci;
  set<int> losque;
  set<int>::iterator It;
  GraphB redtemp(est);
  if (sg.is_directed() != is_directed()) {
    cout << "[Error]: Either the chosen subgraph or the focal graph is directed while the other is undirected. GraphB::list_subgraph_appearances2.\n";
    exit(1);
  }
  for (i=0; i<(number_of_nodes()-1); i++) {
    for (j=(i+1); j < (number_of_nodes()); j++) {
      losque.clear();
      losque.insert(i);
      losque.insert(j);
      put_subgraph_into(redtemp, losque);
      if (redtemp.get_components() == 1)
        if (redtemp.isomorphic_pair(sg, ord)) {
          for (ai=0; ai<2; ai++) {
            ci=0;
            It = losque.begin();
            for (ci=0; ci<ord[ai]; ci++)
              It++;
            os << *It << "\t";
          }
          os << endl;
          res++;
        }
      redtemp.clear();
    }
  }
  delete [] ord;
  return res;
}

int GraphB::list_subgraph_appearances3(GraphB &sg, ostream& os)
{
  if (sg.is_directed() != is_directed()) {
    cout << "[Error]: Either the chosen subgraph or the focal graph is directed while the other is undirected. GraphB::list_subgraph_appearances3.\n";
    exit(1);
  }
  int *ord;
  ord = new int[3];
  basic.fillv0(ord, 3);
  int i, j,l, res=0,ai,ci;
  set<int> losque;
  set<int>::iterator It;
  GraphB redtemp(est);
  for (i=0; i<(number_of_nodes()-2); i++) {
    for (j=(i+1); j < (number_of_nodes()-1); j++) {
      for (l=(j+1); l<number_of_nodes(); l++) {
        losque.clear();
        losque.insert(i);
        losque.insert(j);
        losque.insert(l);
        put_subgraph_into(redtemp, losque);
        if (redtemp.get_components() == 1)
          if (redtemp.isomorphic_pair(sg, ord)){
            for (ai=0; ai<3; ai++) {
              ci=0;
              It = losque.begin();
              for (ci=0; ci<ord[ai]; ci++)
                It++;
              os << *It << "\t";
            }
            os << endl;
            res++;
          }
        redtemp.clear();
      }
    }
  }
  delete [] ord;
  return res;
}

int GraphB::list_subgraph_appearances4(GraphB &sg, ostream& os)
{
  if (sg.is_directed() != is_directed()) {
    cout << "[Error]: Either the chosen subgraph or the focal graph is directed while the other is undirected. GraphB::list_subgraph_appearances4.\n";
    exit(1);
  }
  int *ord;
  ord = new int[4];
  int i, j,l,m,res=0,ai,ci;
  set<int> losque;
  set<int>::iterator It;
  GraphB redtemp(est);
  for (i=0; i<(number_of_nodes()-3); i++) {
    for (j=(i+1); j < (number_of_nodes()-2); j++) {
      for (l=(j+1); l<(number_of_nodes()-1); l++) {
        for (m=(l+1); m < number_of_nodes(); m++) {
          losque.clear();
          losque.insert(i);
          losque.insert(j);
          losque.insert(l);
          losque.insert(m);
          put_subgraph_into(redtemp, losque);
          if (redtemp.get_components() == 1)
            if (redtemp.isomorphic_pair(sg, ord)){
              for (ai=0; ai<4; ai++) {
                ci=0;
                It = losque.begin();
                for (ci=0; ci<ord[ai]; ci++)
                  It++;
                os << *It << "\t";
              }
              os << endl;
              res++;
            }
          redtemp.clear();
        }
      }
    }
  }
  delete [] ord;
  return res;
}

//>

int GraphB::list_subgraph_appearances2(GraphB &sg, string arch, bool **me, int ge, int tis)
{
  int res;
  ofstream sal;
  basic.open_ofstream(sal, arch);
  res = list_subgraph_appearances2(sg, sal, me, ge, tis);
  sal.close();
  return res;
}

int GraphB::list_subgraph_appearances3(GraphB &sg, string arch, bool **me, int ge, int tis)
{
  int res;
  ofstream sal;
  basic.open_ofstream(sal, arch);
  res = list_subgraph_appearances3(sg, sal, me, ge, tis);
  sal.close();
  return res;
}

int GraphB::list_subgraph_appearances4(GraphB &sg, string arch, bool **me, int ge, int tis)
{
  int res;
  ofstream sal;
  basic.open_ofstream(sal, arch);
  res = list_subgraph_appearances4(sg, sal, me, ge, tis);
  sal.close();
  return res;
}

int GraphB::list_subgraph_appearances2(GraphB &sg, ostream& os, bool **me, int ge, int tis)
{
  if (ge != number_of_nodes()) {
    cout << "[Error]: Number of genes does not match the number of rows in expression matrix.\n";//
    exit(1);
  }
  bool siex;
  int *ord;
  ord = new int[2];
  int i, j, res=0,ai,ci,tt;
  set<int> losque;
  set<int>::iterator It;
  GraphB redtemp(est);
  if (sg.is_directed() != is_directed()) {
    cout << "[Error]: Either the chosen subgraph or the focal graph is directed while the other is undirected. GraphB::list_subgraph_appearances2.\n";
    exit(1);
  }
  for (i=0; i<(number_of_nodes()-1); i++) {
    for (j=(i+1); j < (number_of_nodes()); j++) {
      losque.clear();
      losque.insert(i);
      losque.insert(j);
      put_subgraph_into(redtemp, losque);
      if (redtemp.get_components() == 1) {
        siex=false;
        for (tt=0; tt<tis; tt++)
          if (me[i][tt] && me[j][tt]) {
            siex = true;
            break;
          }
        if ((redtemp.isomorphic_pair(sg, ord)) && siex) { //esta ya quedó
          for (ai=0; ai<2; ai++) {
            ci=0;
            It = losque.begin();
            for (ci=0; ci<ord[ai]; ci++)
              It++;
            os << *It << "\t";
          }
          os << endl;
          res++;
        }
      }
      redtemp.clear();
    }
  }
  delete [] ord;
  return res;
}

int GraphB::list_subgraph_appearances3(GraphB &sg, ostream& os, bool **me, int ge, int tis)
{
  if (ge != number_of_nodes()) {
    cout << "[Error]: Number of genes does not match the number of rows in expression matrix.\n";//
    exit(1);
  }
  bool siex;
  if (sg.is_directed() != is_directed()) {
    cout << "[Error]: Either the chosen subgraph or the focal graph is directed while the other is undirected. GraphB::list_subgraph_appearances3.\n";
    exit(1);
  }
  int *ord;
  ord = new int[3];
  basic.fillv0(ord, 3);
  int i, j,l, res=0,ai,ci, tt;
  set<int> losque;
  set<int>::iterator It;
  GraphB redtemp(est);
  for (i=0; i<(number_of_nodes()-2); i++) {
    for (j=(i+1); j < (number_of_nodes()-1); j++) {
      for (l=(j+1); l<number_of_nodes(); l++) {
        losque.clear();
        losque.insert(i);
        losque.insert(j);
        losque.insert(l);
        put_subgraph_into(redtemp, losque);
        if (redtemp.get_components() == 1) {
          siex = false;
          for (tt=0; tt<tis; tt++)
            if (me[i][tt] && me[j][tt] && me[l][tt]) {
              siex = true;
              break;
            }
          if ((redtemp.isomorphic_pair(sg, ord)) && siex) {
            for (ai=0; ai<3; ai++) {
              ci=0;
              It = losque.begin();
              for (ci=0; ci<ord[ai]; ci++)
                It++;
              os << *It << "\t";
            }
            os << endl;
            res++;
          }
        }
        redtemp.clear();
      }
    }
  }
  delete [] ord;
  return res;
}

int GraphB::list_subgraph_appearances4(GraphB &sg, ostream& os, bool **me, int ge, int tis)
{
  if (ge != number_of_nodes()) {
    cout << "[Error]: Number of genes does not match the number of rows in expression matrix.\n";//
    exit(1);
  }
  bool siex;
  if (sg.is_directed() != is_directed()) {
    cout << "[Error]: Either the chosen subgraph or the focal graph is directed while the other is undirected. GraphB::list_subgraph_appearances4.\n";
    exit(1);
  }
  int *ord;
  ord = new int[4];
  int i, j,l,m,res=0,ai,ci, tt;
  set<int> losque;
  set<int>::iterator It;
  GraphB redtemp(est);
  for (i=0; i<(number_of_nodes()-3); i++) {
    for (j=(i+1); j < (number_of_nodes()-2); j++) {
      for (l=(j+1); l<(number_of_nodes()-1); l++) {
        for (m=(l+1); m < number_of_nodes(); m++) {
          losque.clear();
          losque.insert(i);
          losque.insert(j);
          losque.insert(l);
          losque.insert(m);
          put_subgraph_into(redtemp, losque);
          if (redtemp.get_components() == 1) {
            siex = false;
            for (tt=0; tt<tis; tt++)
              if (me[i][tt] && me[j][tt] && me[l][tt] && me[m][tt]) {
                siex = true;
                break;
              }
            if ((redtemp.isomorphic_pair(sg, ord)) && siex) {
              for (ai=0; ai<4; ai++) {
                ci=0;
                It = losque.begin();
                for (ci=0; ci<ord[ai]; ci++)
                  It++;
                os << *It << "\t";
              }
              os << endl;
              res++;
            }
          }
          redtemp.clear();
        }
      }
    }
  }
  delete [] ord;
  return res;
}
//<


int GraphB::list_subgraph_appearances_using_names2(GraphB &sg, string arch)
{
  int res;
  ofstream sal;
  basic.open_ofstream(sal, arch);
  res = list_subgraph_appearances_using_names2(sg, sal);
  sal.close();
  return res;
}

int GraphB::list_subgraph_appearances_using_names3(GraphB &sg, string arch)
{
  int res;
  ofstream sal;
  basic.open_ofstream(sal, arch);
  res = list_subgraph_appearances_using_names3(sg, sal);
  sal.close();
  return res;
}

int GraphB::list_subgraph_appearances_using_names4(GraphB &sg, string arch)
{
  int res;
  ofstream sal;
  basic.open_ofstream(sal, arch);
  res = list_subgraph_appearances_using_names4(sg, sal);
  sal.close();
  return res;
}

int GraphB::list_subgraph_appearances_using_names2(GraphB &sg, ostream& os)
{
  int *ord;
  ord = new int[2];
  int i, j, res=0,ai,ci;
  set<int> losque;
  set<int>::iterator It;
  GraphB redtemp(est);
  if (sg.is_directed() != is_directed()) {
    cout << "[Error]: Either the chosen subgraph or the focal graph is directed while the other is undirected. GraphB::list_subgraph_appearances_using_names2.\n";
    exit(1);
  }
  for (i=0; i<(number_of_nodes()-1); i++) {
    for (j=(i+1); j < (number_of_nodes()); j++) {
      losque.clear();
      losque.insert(i);
      losque.insert(j);
      put_subgraph_into(redtemp, losque);
      if (redtemp.get_components() == 1)
        if (redtemp.isomorphic_pair(sg, ord)) {
          for (ai=0; ai<2; ai++) {
            ci=0;
            It = losque.begin();
            for (ci=0; ci<ord[ai]; ci++)
              It++;
            os << get_name(*It) << "\t";
          }
          os << endl;
          res++;
        }
      redtemp.clear();
    }
  }
  delete [] ord;
  return res;
}

int GraphB::list_subgraph_appearances_using_names3(GraphB &sg, ostream& os)
{
  if (sg.is_directed() != is_directed()) {
    cout << "[Error]: Either the chosen subgraph or the focal graph is directed while the other is undirected. GraphB::list_subgraph_appearances_using_names3.\n";
    exit(1);
  }
  int *ord;
  ord = new int[3];
  int i, j,l, res=0,ai,ci;
  set<int> losque;
  set<int>::iterator It;
  GraphB redtemp(est);
  for (i=0; i<(number_of_nodes()-2); i++) {
    for (j=(i+1); j < (number_of_nodes()-1); j++) {
      for (l=(j+1); l<number_of_nodes(); l++) {
        losque.clear();
        losque.insert(i);
        losque.insert(j);
        losque.insert(l);
        put_subgraph_into(redtemp, losque);
        if (redtemp.get_components() == 1)
          if (redtemp.isomorphic_pair(sg, ord)){
            for (ai=0; ai<3; ai++) {
              ci=0;
              It = losque.begin();
              for (ci=0; ci<ord[ai]; ci++)
                It++;
              os << get_name(*It) << "\t";
            }
            os << endl;
            res++;
          }
        redtemp.clear();
      }
    }
  }
  delete [] ord;
  return res;
}

int GraphB::list_subgraph_appearances_using_names4(GraphB &sg, ostream& os)
{
  if (sg.is_directed() != is_directed()) {
    cout << "[Error]: Either the chosen subgraph or the focal graph is directed while the other is undirected. GraphB::list_subgraph_appearances_using_names4.\n";
    exit(1);
  }
  int *ord;
  ord = new int[4];
  int i, j,l,m,res=0,ai,ci;
  set<int>::iterator It;
  set<int> losque;
  GraphB redtemp(est);
  for (i=0; i<(number_of_nodes()-3); i++) {
    for (j=(i+1); j < (number_of_nodes()-2); j++) {
      for (l=(j+1); l<(number_of_nodes()-1); l++) {
        for (m=(l+1); m < number_of_nodes(); m++) {
          losque.clear();
          losque.insert(i);
          losque.insert(j);
          losque.insert(l);
          losque.insert(m);
          put_subgraph_into(redtemp, losque);
          if (redtemp.get_components() == 1)
            if (redtemp.isomorphic_pair(sg, ord)){
              for (ai=0; ai<4; ai++) {
                ci=0;
                It = losque.begin();
                for (ci=0; ci<ord[ai]; ci++)
                  It++;
                os << get_name(*It) << "\t";
              }
              os << endl;
              res++;
            }
          redtemp.clear();
        }
      }
    }
  }
  delete [] ord;
  return res;
}

//##

int GraphB::list_subgraph_appearances_using_names2(GraphB &sg, string arch, bool **me, int ge, int tis)
{
  int res;
  ofstream sal;
  basic.open_ofstream(sal, arch);
  res = list_subgraph_appearances_using_names2(sg, sal, me, ge, tis);
  sal.close();
  return res;
}

int GraphB::list_subgraph_appearances_using_names3(GraphB &sg, string arch, bool **me, int ge, int tis)
{
  int res;
  ofstream sal;
  basic.open_ofstream(sal, arch);
  res = list_subgraph_appearances_using_names3(sg, sal, me, ge, tis);
  sal.close();
  return res;
}

int GraphB::list_subgraph_appearances_using_names4(GraphB &sg, string arch, bool **me, int ge, int tis)
{
  int res;
  ofstream sal;
  basic.open_ofstream(sal, arch);
  res = list_subgraph_appearances_using_names4(sg, sal, me, ge, tis);
  sal.close();
  return res;
}

int GraphB::list_subgraph_appearances_using_names2(GraphB &sg, ostream& os, bool **me, int ge, int tis)
{
  if (ge != number_of_nodes()) {
    cout << "[Error]: Number of genes does not match the number of rows in expression matrix.\n";//
    exit(1);
  }
  bool siex;
  int *ord;
  ord = new int[2];
  int i, j, res=0,ai,ci,tt;
  set<int> losque;
  set<int>::iterator It;
  GraphB redtemp(est);
  if (sg.is_directed() != is_directed()) {
    cout << "[Error]: Either the chosen subgraph or the focal graph is directed while the other is undirected. GraphB::list_subgraph_appearances_using_names2.\n";
    exit(1);
  }
  for (i=0; i<(number_of_nodes()-1); i++) {
    for (j=(i+1); j < (number_of_nodes()); j++) {
      losque.clear();
      losque.insert(i);
      losque.insert(j);
      put_subgraph_into(redtemp, losque);
      if (redtemp.get_components() == 1) {
        siex=false;
        for (tt=0; tt<tis; tt++)
          if (me[i][tt] && me[j][tt]) {
            siex = true;
            break;
          }
        if ((redtemp.isomorphic_pair(sg, ord)) && siex) {
          for (ai=0; ai<2; ai++) {
            ci=0;
            It = losque.begin();
            for (ci=0; ci<ord[ai]; ci++)
              It++;
            os << get_name(*It) << "\t";
          }
          os << endl;
          res++;
        }
      }
      redtemp.clear();
    }
  }
  delete [] ord;
  return res;
}

int GraphB::list_subgraph_appearances_using_names3(GraphB &sg, ostream& os, bool **me, int ge, int tis)
{
  if (ge != number_of_nodes()) {
    cout << "[Error]: Number of genes does not match the number of rows in expression matrix.\n";//
    exit(1);
  }
  bool siex;
  if (sg.is_directed() != is_directed()) {
    cout << "[Error]: Either the chosen subgraph or the focal graph is directed while the other is undirected. GraphB::list_subgraph_appearances_using_names3.\n";
    exit(1);
  }
  int *ord;
  ord = new int[3];
  int i, j,l, res=0,ai,ci,tt;
  set<int> losque;
  set<int>::iterator It;
  GraphB redtemp(est);
  for (i=0; i<(number_of_nodes()-2); i++) {
    for (j=(i+1); j < (number_of_nodes()-1); j++) {
      for (l=(j+1); l<number_of_nodes(); l++) {
        losque.clear();
        losque.insert(i);
        losque.insert(j);
        losque.insert(l);
        put_subgraph_into(redtemp, losque);
        if (redtemp.get_components() == 1) {
          siex=false;
          for (tt=0; tt<tis; tt++)
            if (me[i][tt] && me[j][tt] && me[l][tt]) {
              siex = true;
              break;
            }
          if ((redtemp.isomorphic_pair(sg, ord)) && siex) {
            for (ai=0; ai<3; ai++) {
              ci=0;
              It = losque.begin();
              for (ci=0; ci<ord[ai]; ci++)
                It++;
              os << get_name(*It) << "\t";
            }
            os << endl;
            res++;
          }
        }
        redtemp.clear();
      }
    }
  }
  delete [] ord;
  return res;
}

int GraphB::list_subgraph_appearances_using_names4(GraphB &sg, ostream& os, bool **me, int ge, int tis)
{
  if (ge != number_of_nodes()) {
    cout << "[Error]: Number of genes does not match the number of rows in expression matrix.\n";//
    exit(1);
  }
  bool siex;
  if (sg.is_directed() != is_directed()) {
    cout << "[Error]: Either the chosen subgraph or the focal graph is directed while the other is undirected. GraphB::list_subgraph_appearances_using_names4.\n";
    exit(1);
  }
  int *ord;
  ord = new int[4];
  int i, j,l,m,res=0,ai,ci,tt;
  set<int>::iterator It;
  set<int> losque;
  GraphB redtemp(est);
  for (i=0; i<(number_of_nodes()-3); i++) {
    for (j=(i+1); j < (number_of_nodes()-2); j++) {
      for (l=(j+1); l<(number_of_nodes()-1); l++) {
        for (m=(l+1); m < number_of_nodes(); m++) {
          losque.clear();
          losque.insert(i);
          losque.insert(j);
          losque.insert(l);
          losque.insert(m);
          put_subgraph_into(redtemp, losque);
          if (redtemp.get_components() == 1) {
            siex=false;
            for (tt=0; tt<tis; tt++)
              if (me[i][tt] && me[j][tt] && me[l][tt] && me[m][tt]) {
                siex = true;
                break;
              }
            if ((redtemp.isomorphic_pair(sg, ord)) && siex) {
              for (ai=0; ai<4; ai++) {
                ci=0;
                It = losque.begin();
                for (ci=0; ci<ord[ai]; ci++)
                  It++;
                os << get_name(*It) << "\t";
              }
              os << endl;
              res++;
            }
          }
          redtemp.clear();
        }
      }
    }
  }
  delete [] ord;
  return res;
}
//##



bool GraphB::an_isomorph_contains_these_interactions1(GraphB &subg)
{
  bool res;
  res = contains_these_interactions(subg);
  return res;
}

bool GraphB::an_isomorph_contains_these_interactions2(GraphB &subg)
{
  int tam = size;
  int rows = basic.factorial(tam);
  bool res = false;
  GraphB whig(est);
  int i;
  for (i=0; i<rows; i++) {
    subg.transform_to_isomorph(whig, comb2[i]);
    if (contains_these_interactions(whig)) {
      res = true;
      whig.clear();
      break;
    }
    whig.clear();
  }
  return res;
}

bool GraphB::an_isomorph_contains_these_interactions3(GraphB &subg)
{
  int tam = size;
  int rows = basic.factorial(tam);
  bool res = false;
  GraphB whig(est);
  int i;
  for (i=0; i<rows; i++) {
    subg.transform_to_isomorph(whig, comb3[i]);
    if (contains_these_interactions(whig)) {
      res = true;
      whig.clear();
      break;
    }
    whig.clear();
  }
  return res;
}

bool GraphB::an_isomorph_contains_these_interactions4(GraphB &subg)
{
  int tam = size;
  int rows = basic.factorial(tam);
  bool res = false;
  GraphB whig(est);
  int i;
  for (i=0; i<rows; i++) {
    subg.transform_to_isomorph(whig, comb4[i]);
    if (contains_these_interactions(whig)) {
      res = true;
      whig.clear();
      break;
    }
    whig.clear();
  }
  return res;
}


bool GraphB::isomorphic_pair1(GraphB &templ, int *ord)
{
  int i,j;
  i = count_selfinteracting_nodes();
  j = templ.count_selfinteracting_nodes();
  bool res = false;
  if (i==j)
    res = true;
  ord[0]=0;
  return res;
}

bool GraphB::isomorphic_pair2(GraphB &templ, int *ord)
{
  int tam = size;
  int rows = basic.factorial(tam);
  bool res = false;
  GraphB whig(est);
  int i,se,sd,j;
  se = count_selfinteracting_nodes();
  sd = templ.count_selfinteracting_nodes();
  if (se == sd) {
    for (i=0; i<rows; i++) {
      templ.transform_to_isomorph(whig, comb2[i]);
      if (equal_nw(whig)) {
        res = true;
        whig.clear();
        for (j=0; j<2; j++)
          ord[j] = comb2[i][j];
        break;
      }
      whig.clear();
    }
  }
  return res;
}

bool GraphB::isomorphic_pair3(GraphB &templ, int *ord)
{
  int tam = size;
  int rows = basic.factorial(tam);
  bool res = false;
  GraphB whig(est);
  int i,j,se,sd;
  se = count_selfinteracting_nodes();
  sd = templ.count_selfinteracting_nodes();
  if (se == sd) {
    for (i=0; i<rows; i++) {
      templ.transform_to_isomorph(whig, comb3[i]);
      if (equal_nw(whig)) {
        res = true;
        whig.clear();
        for (j=0; j<3; j++)
          ord[j] = comb3[i][j];
        break;
      }
      whig.clear();
    }
  }
  return res;
}

bool GraphB::isomorphic_pair4(GraphB &templ, int *ord)
{
  int tam = size;
  int rows = basic.factorial(tam);
  bool res = false;
  GraphB whig(est);
  int i,j,se,sd;
  se = count_selfinteracting_nodes();
  sd = templ.count_selfinteracting_nodes();
  if (se == sd) {
    for (i=0; i<rows; i++) {
      templ.transform_to_isomorph(whig, comb4[i]);
      if (equal_nw(whig)) {
        res = true;
        whig.clear();
        for (j=0; j<4; j++)
          ord[j] = comb4[i][j];
        break;
      }
      whig.clear();
    }
  }
  return res;
}

bool GraphB::isomorphic_pair1(GraphB &templ)
{
  int i,j;
  i = count_selfinteracting_nodes();
  j = templ.count_selfinteracting_nodes();
  bool res = false;
  if (i==j)
    res = true;
  return res;
}

bool GraphB::isomorphic_pair2(GraphB &templ)
{
  int tam = size;
  int rows = basic.factorial(tam);
  bool res = false;
  GraphB whig(est);
  int i,se,sd;
  se = count_selfinteracting_nodes();
  sd = templ.count_selfinteracting_nodes();
  if (se == sd) {
    for (i=0; i<rows; i++) {
      templ.transform_to_isomorph(whig, comb2[i]);
      if (equal_nw(whig)) {
        res = true;
        whig.clear();
        break;
      }
      whig.clear();
    }
  }
  return res;
}

bool GraphB::isomorphic_pair3(GraphB &templ)
{
  int tam = size;
  int rows = basic.factorial(tam);
  bool res = false;
  GraphB whig(est);
  int i,se,sd;
  se = count_selfinteracting_nodes();
  sd = templ.count_selfinteracting_nodes();
  if (se == sd) {
    for (i=0; i<rows; i++) {
      templ.transform_to_isomorph(whig, comb3[i]);
      if (equal_nw(whig)) {
        res = true;
        whig.clear();
        break;
      }
      whig.clear();
    }
  }
  return res;
}

bool GraphB::isomorphic_pair4(GraphB &templ)
{
  int tam = size;
  int rows = basic.factorial(tam);
  bool res = false;
  GraphB whig(est);
  int i,se,sd;
  se = count_selfinteracting_nodes();
  sd = templ.count_selfinteracting_nodes();
  if (se == sd) {
    for (i=0; i<rows; i++) {
      templ.transform_to_isomorph(whig, comb4[i]);
      if (equal_nw(whig)) {
        res = true;
        whig.clear();
        break;
      }
      whig.clear();
    }
  }
  return res;
}
