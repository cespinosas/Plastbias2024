#ifndef GRAPHI_H
#define GRAPHI_H

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

class GraphI
{
public:
  GraphI();
  GraphI(Alea& jacta);
  void start_rng(Alea& jacta);
  GraphI(int n, Alea& jacta, bool dirnw);
  void make_nw(int n, bool dirnw);
  void make_robust_nw(int n, int e, int *cia, int *pa, int *cib, int *pb);
  void copy(GraphI &templ);
  void make_copy(GraphI& vacia);
	void put_subgraph_into(GraphI& nueva, set<int> &cuales, bool withnames);
	void put_subgraph_into(GraphI& nueva, set<int> &cuales);
	void make_subgraph_of(GraphI &templ, set<int> &cuales, bool withnames);
	void make_subgraph_of(GraphI &templ, set<int> &cuales);
  void copy_wdeg(GraphI &templ);
  GraphI(int n, double c, Alea& jacta, bool dirnw);
  void make_nw(int n, double c, bool dirnw);
  GraphI(int n, int e, Alea& jacta, bool dirnw);
  void make_nw(int n, int e, bool dirnw);
  
  
  void rwalk_1traj(int n, int e, int *fp, int *ic, int cuan, string arch, int mame, int numut, bool estri, int numarch, int archdesde);
  void rwalk_2fp1t(int n, int e, int *fp1, int *fp2, int *ic1, int cuan, string arch, int mame, int numut, bool estri, int numarch, int archdesde);
  void rwalk_2fp2t(int n, int e, int *fp1, int *fp2, int *ic1, int *ic2, int cuan, string arch, int mame, int numut, bool estri, int numarch, int archdesde);
  void rwalk_2fp1tn1(int n, int e, int *fp1, int *fp2, int *ic1, int cuan, string arch, int mame, int numut, bool estri, int numarch, int archdesde);
  void rwalk_2fp1tn2(int n, int e, int *fp1, int *fp2, int *ic1, int cuan, string arch, int mame, int numut, bool estri, int numarch, int archdesde);
  
  void clear();
	void clear_edge_count();
  void clear_loops();
  void clear_comp();
  void clear_adj();
  void clear_scc();
  void clear_bc();
  void clear_dima();
  void clear_attractor();
  void clear_attractors();
  void clear_ds();
  void clear_dict();
  void clear_deg();
	void clear_degdist(); //f
  void clear_allatts(); //2024
  int get_indegree(int no);
  int get_outdegree(int no);
  int get_degree(int no);
  void get_all_degrees();
	void get_degdist(); //f
	void get_degdist(int *empvec, int tam); //f
	void get_odegdist(int *empvec, int tam); //f
	void get_idegdist(int *empvec, int tam); //f
  
  //access
  int number_of_edges();
  int number_of_nodes();
  int weight(int source, int target);
  bool names_exist();
  string get_name(int node);
  bool is_directed();
  bool already_deg();
	int name_to_int(string cual);
	int name_to_int(string cual, bool checkornot);
	bool is_name(string elna);
  bool attractor_exists();
  bool attractors_exist();
  
  //modification
  void force_interaction(int source, int target, int value);
  void force_interaction_undir(int source, int target, int value);
  void change_interaction(int source, int target, int value);
  void change_interaction_undir(int source, int target, int value);
  void change_name(int node, string nn);
  void set_undirected(); //only changes label
  void set_directed();
	void force_undirected(); //makes undirected and changes label
	void transform_to_isomorph(GraphI& vacia, int* vec);
  
  //formatting
  void get_names_from_file(string arch);
	void get_names_from_file(istream& en);
	void assign_names(string *nombres, int ent);
  void printnw(ostream& sal);
  void print_dot(string arch);
  void print_dot_wn(string arch);
	void print_dot_wo_labels_circ(string arch, double radius);
  void print_dot_circ(string arch, double radius);
  void print_cytoscape(string arch);
  void print_cytoscape_wn(string arch);
  void print_state(ostream& sal); //nuevisimo

  
	//exclusive
  //import
  void get_dir_nw_from_file(int nn, istream& en, int e); //nuevisimo
  void get_dir_nw_from_file(int nn, string arch);
  void get_dir_nw_from_file(int nn, istream& en);
  void get_undir_nw_from_file(int nn, string arch);
  void get_undir_nw_from_file(int nn, istream& en);
	void get_dir_nw_from_file_wn(int nn, string arch, string archnam);
  void get_dir_nw_from_file_wn(int nn, istream& en, istream& enna);
  void get_undir_nw_from_file_wn(int nn, string arch, string archnam);
  void get_undir_nw_from_file_wn(int nn, istream& en, istream& enna);
	
  //export
  void export_nw(string arch);
  void export_nw(ostream& fs);

  //status
	
  //analyses
  int distance_from_nw(GraphI& otra);
  bool equal_adjmat(GraphI& otra);
  bool equal_nw(GraphI &templ);
	set<int> get_upstream(int no);
  void get_downstream();
  
  set<int> get_out_component(int no);
  set<int> get_in_component(int no);
  void get_adjacency_matrix();
  void get_distance_matrix_ls(); //for large & sparse networks
  void get_distance_matrix();
  int get_scc();
  int number_of_sccs();
  int number_of_components();
  int get_components();
  void get_bc();
  double get_bc(int nodo);
  int nupaths(int u, int v, int *emptyvec);
  list<int> GetSucc(int v, list<int> &omega, int q);
  void get_loops();
	
  //access analyses
  bool there_is_path(int from, int to);
  int length_shortest_path(int from, int to);
  set<int> nodes_influenced_by(int n); //directly downstream
  int in_which_component(int n);
  int in_which_scc(int n);
  int number_of_paths(int from, int to, int len, int* vec);
  int number_of_loops();
  int number_of_loops(int len);
  void switches(int vec);
  void synchrony();
  void set_as_state(int *vec);
  void find_an_attractor();
  void find_attractors(int **ic, int nic);
  int attractor_size();
  int path_length();
  int attractor_element(int row, int node);
  void print_attractor(ostream& sal);
  int attractor_size(int wh);
  int attractor_element(int wh, int row, int node);
  void print_attractor(int wh, ostream& sal);
  int number_of_attractors();
  double robustness(int *ic, int n, double wol);
  double robustness(int *ic, int n, double& roper, double& roperz, double wol);  //roperz excludes zeros
  double robustness(int *ic, int n, double& roper, double& roperz, double& totnnf, int& nufedi, double& totnnfz, int& nufediz, double& ronoz, double& rof1, double& rof1z, double wol); //roperz excludes zeros
  double robustness(int *ic, int n, double& roper, double& roperz, double& totnnf, int& nufedi, double& totnnfz, int& nufediz, double& ronoz, double& rof1, double& rof1z, double& ropl, double& roplz, double& roplper, double& roplperz, double wol); 
  void robustness_mult(int **ic, double *robs, int nic, int n, double wol);
  int penetrance_ap_n1(int *cior, int *altphen);
  void penetrance_n1_2p(int *cior, int *pA, int *pB, int& countA, int& countB);
  int penetrance_ap_n2_ex(int *cior, int *altphen);
  void penetrance_n2_ex_2p(int *cior, int *pA, int *pB, int& countA, int& countB);
  int mutational_access(int *ci, int *altphen);
  void mutational_access_2p(int *ci, int *pA, int *pB, int& countA, int& countB);
  void mutational_access_2p_eqmu(int *ci, int *pA, int *pB, int& countA, int& countB);

  //for evolution
  void build_offspring(GraphI &mother, GraphI &father, bool *frommom);
  void mate(GraphI &mother, GraphI &father);
  void mutate_old(double con); //old; one mutation per genotype
  void mutate(double con); //201707-\infty; It is already decided that mutation takes place. a single mutation somewhere. (mutation may be without any effect).
  void mutate_gene(int gene, double con);
  void mutate_gene_mp5(int gene);
  void mutate_gene_mf(int gene);
  void consider_mutation(double muratepg, double con); //probably multiple mutations per genotype
  void consider_mutation_mp5(double muratepg);
  void consider_mutation_mf(double muratepg);
  void duplicate(int gene);
  void duplicate();
	
	
	//motifs
	int count_selfinteracting_nodes();
	int count_non_selfinteracting_nodes();
	
  //modularity
  void build_moma_u(double **momaundi);
  void build_moma_d(double **momadi);
  void build_moma_out(double **momadi);
  void build_moma_in(double **momadi);
  int components_gt1();
  int sccs_gt1();
  void check_props_of_rdset(int samplesiz, int minsw, int maxsw, double **meansdminma);
  void check_props_of_rdset(int samplesiz, int minsw, int maxsw, double **meansdminma, set<set<int> > &predpar);
  void check_Qd_in_rdset(int samplesiz, int minsw, int maxsw, double *meansdminma); //1row
  void check_Qd_in_rdset(int samplesiz, int minsw, int maxsw, double **meansdminma, set<set<int> > &predpar); //2rows
  void build_moma_u_aux(int **adjundi, int *degundi, int edgundi, double **momaundi); //undirected networks
  void mod_partitions_from_vector(double *vecs, int tam, int *dicc, set<int> &parta, set<int> &partb); //it creates sets according to partition in vecs with names (int numbers) defined by dicc
  double spectral_method(double **moduloc, int tam, double *vecs);
  double mod_after_spectral(double **matloc, int tam, double *vecs);
  void shake_kl(double mienQ, double **matloc, int tam, double *vecs);
  double eval_mod_subm(double **mamod, double *vecs, int tam); //Ojo, con la dirigida de Guimera tal vez falta un factor. Al construirla, multiplicar por 4 edges
  void partition_to_vector(const set<set<int> > &equipos, int *vop); //requires full partition. vop of size size
  double eval_mod(double **mamod, int *vop);
  double eval_mod(double **mamod, set<set<int> > &equipos);
  double maxmod_u_per_part_and_dd(int *vop); //para otra forma de normalizar. Cuando se hace siempre la misma partición... (habría que defender por que aleatorias con la misma dd)
  double maxmod_u_per_part_and_dd(set<set<int> > &equipos);
  double maxmod_d_per_part_and_dd(int *vop);
  double maxmod_d_per_part_and_dd(set<set<int> > &equipos);
  void adjust_modmat(double **mator, set<int> &parta, int tam, int *dicc, double **matloc);
  void iterative_newman06(double **mator, set<set<int> > &equipos);
  void iterative_newman06(double **mator, set<int> &parta, set<set<int> > &res);
  void split_in_2mods(double **mator, set<set<int> > &equipos);

  void get_all_attractors(); //2024
  int data_of_all_attractors(int** &tosat, int* &tobasinsizes, int* &toperiods);//2024 //regresa numero de atractores
  void all_basins(int* &cuencas);//2024
  int phens_from_neighbor_cis(int *ci, int** &neigh_atts, int* &neigh_per, int* &nuneighci) ; //2024 returns #of diff attractors 
  
private:
  bool allattsya; //2024
  int **allatts; //2024
  int *basinsizes; //2024 de tamaño efectivo como numero de atractores numtotat
  int *allatt_periods; //2024
  int *basins; //2024 de tamaño 2 a la N
 // int numaxat; //2024
  int numtotat; //2024
  int numallci; //2024
  
  int size;
  int **nw;
  int **dima;
  int **matadya;
  bool yamadya;
  Alea est;
  string *dict;
  bool ionary;
  Basics basic;
  bool directed;
  //exclusive
  set<int> *downstream;
  bool infl;
  bool dist;
  bool yacomp;
  int numbcomp;
  int *components;
  int *scc;
  int numbscc;
  bool cfc;
  double *bc;
  bool bcya;
  int *loops;
  int **guada;
  int totalloops;
  bool floops;
  int numofe;
  bool yae;
  int *outdegree;
  int *indegree;
  int *degree;
  bool yadeg;
	int *degdist; // falta poner en graphi y graphc
	int *odegdist;
	int *idegdist;
	bool yadegdist; // falta poner en graphi y graphc... recordar que falta tambien en inicializadas
  int *edo;
  int *ima;
  int **attractor;
  int atsize;
  int palen;
  bool yatra;
  int ***attractors;
  int *atsizes;
  int numatrs;
  bool yatras;
  
  //auxiliary functions
	/*   void get_distance_matrix_sym(); */
	/*   void get_distance_sym_(); */
	void set_default_exclusive_vars();
  int calc_indegree(int no);
  int calc_outdegree(int no);
  int calc_degree(int no);
  void prepare_for_degrees();
	void prepare_for_degdist();
	void priv_get_degdist(); // falta poner en graphi y graphc// esta y la anterior siempre se deben llamar juntas
  void aswitch();
  double get_maxcomb(int n, int e); //nuevisimo subopt
  double get_maxcomb2(int n, int e, int maxmo); //nuevisimo subopt
  int formaxmo(int n, int *fp1, int *fp2);
};

//Notas de cosas que falta por hacer:
//-función para obtener todos los atractores
//-función para muestrear atractores
//-Incluir distintas formas de actualización de estado. (tal vez estado no booleano?) (otra lib u opción aquí??)

#endif

