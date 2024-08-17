#ifndef GRAPHB_H
#define GRAPHB_H

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
#include "graphc.h"

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

//desconfío de get_bc en todas


class GraphB
{
public:
	GraphB();
	GraphB(Alea& jacta);
	void start_rng(Alea& jacta);
	GraphB(int n, Alea& jacta, bool dirnw);
	void make_nw(int n, bool dirnw);
	void copy(GraphB &templ);
	void copy(GraphB &templ, bool withnames);
  void make_copy(GraphB& nueva);
	void put_subgraph_into(GraphB& nueva, set<int> &cuales, bool withnames);
	void put_subgraph_into(GraphB& nueva, set<int> &cuales);
	void make_subgraph_of(GraphB &templ, set<int> &cuales, bool withnames);
	void make_subgraph_of(GraphB &templ, set<int> &cuales);
	void copy_wdeg(GraphB &templ);
	void copy_wdeg(GraphB &templ, bool withnames);
	GraphB(int n, double c, Alea& jacta, bool dirnw);
	void make_nw(int n, double c, bool dirnw);
	GraphB(int n, int e, Alea& jacta, bool dirnw);
	void make_nw(int n, int e, bool dirnw);
	void clear();
	void clear_preps_to_count_sg(); //exclusive of graphB
	void clear_sg_catalogue(); //exclusive of graphB
	void clear_db_for_isomorphism(); //exclusive of graphB
	void clear_edge_count();
	void clear_loops();
	void clear_comp();
	void clear_scc();
  void clear_adj();
	void clear_bc();
	void clear_dima();
	void clear_ds();
	void clear_dict();
	void clear_deg();
	void clear_degdist();
	int get_indegree(int no);
	int get_outdegree(int no);
	int get_degree(int no);
	void get_all_degrees();
	void get_degdist();
	void get_degdist(int *empvec, int tam);
	void get_odegdist(int *empvec, int tam);
	void get_idegdist(int *empvec, int tam);
	
	//access
	int number_of_edges();
	int number_of_nodes();
	bool weight(int source, int target);
	bool names_exist();
	string get_name(int node);
	bool is_directed();
	bool already_deg();
	int name_to_int(string cual);
	int name_to_int(string cual, bool checkornot);
	
	bool is_name(string elna);
	
	//modification
	void force_interaction(int source, int target, bool value);
	void force_interaction_undir(int source, int target, bool value);
	void change_interaction(int source, int target, bool value);
	void change_interaction_undir(int source, int target, bool value);
	void change_name(int node, string nn);
	void set_undirected(); //only changes label
	void set_directed();
	void force_undirected(); //makes undirected and changes label
	void transform_to_isomorph(GraphB& vacia, int* vec);
	
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
	
	//exclusive
	//import
	void get_dir_nw_from_file(int nn, string arch);
	void get_dir_nw_from_file(int nn, istream& en);
  void get_dir_nw_from_file_wn(string arch);
  void get_undir_nw_from_file_wn(string arch);
	void get_undir_nw_from_file(int nn, string arch);
	void get_undir_nw_from_file(int nn, istream& en);
	
	void get_dir_nw_from_file_wn(int nn, string arch, string archnam);
	void get_dir_nw_from_file_wn(int nn, istream& en, istream& enna);
	void get_undir_nw_from_file_wn(int nn, string arch, string archnam);
	void get_undir_nw_from_file_wn(int nn, istream& en, istream& enna);
	void get_from_graphi(GraphI &templ);
	void get_from_graphc(GraphC &templ);
  
  //export
  void export_nw(string arch);
  void export_nw(ostream& fs);
	
  //status
	
	//analyses
	bool equal_nw(GraphB &templ);
	bool isomorphic_pair(GraphB &templ);
	bool isomorphic_pair(GraphB &templ, int *ord);
	
	
	set<int> get_upstream(int no); //directly upstream
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
	void get_bc(); //no estoy seguro de que esté bien
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
	void switches_preserving_sg1(int vec);
	void switches_preserving_sg2(int vec);
	void print_list_of_list_of_nodes(list<list<int> > &lalis, ostream& sal);
	
	//for evolution
	void mate(GraphB &mother, GraphB &father);
	void mutate(int gene, double wol);
	void duplicate(int gene);
	void duplicate();
	
	//motifs
	int count_selfinteracting_nodes();
	int count_non_selfinteracting_nodes();
	//motifexclusive of graphb
	bool contains_these_interactions(GraphB &subg);
	bool an_isomorph_contains_these_interactions(GraphB &subg);
	
	
	void load_database_for_isomorphism(); // remember to check that subgraphcount is set to zero every time the nw is rearranged
	void load_sg_catalogue();
	void prepare_to_count_subgraphs();
	int loose_count_this_subgraph(GraphB &sg);
	int count_this_subgraph(GraphB &sg);
	int list_subgraph_appearances(GraphB &sg, string arch);
	int list_subgraph_appearances(GraphB &sg, ostream& os);
	int list_subgraph_appearances(GraphB &sg, string arch, bool **me, int ge, int tis);
	int list_subgraph_appearances(GraphB &sg, ostream& os, bool **me, int ge, int tis);
  
	int list_subgraph_appearances_using_names(GraphB &sg, string arch);
	int list_subgraph_appearances_using_names(GraphB &sg, ostream& os);
  
	int list_subgraph_appearances_using_names(GraphB &sg, string arch, bool **me, int ge, int tis);
	int list_subgraph_appearances_using_names(GraphB &sg, ostream& os, bool **me, int ge, int tis);
  
	void count_subgraphs_2();
	void count_subgraphs_3();
	void count_subgraphs_4();
	
	void count_subgraphs_2(bool **me, int ge, int tis);
	void count_subgraphs_3(bool **me, int ge, int tis);
	void count_subgraphs_4(bool **me, int ge, int tis);
	
	
	int how_many_subgraphs(int tam, int which); //access to previously counted subgraphs
	int count_this_subgraph(int tam, int which);
	int loose_count_this_subgraph(int tam, int which);
	set<int> loose_list(GraphB &sg); //constructs list of subgraphs that include subgraph sg
	int number_of_sgs(int nnodes, int cuasg);
	
	
	int nupre_2, nupre_3, nupre_4;
	GraphB *LM_2, *LM_3, *LM_4;
  
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
  
  
private:
	int size;
	bool **nw;
	int **dima;
  int **matadya;
  bool yamadya;
  
	Alea est;
	string *dict;
	bool ionary;
	Basics basic;
	bool directed;
	//exclusive
	bool yamotdab; //only gb. Initialization of gb includes setting this to false
	int **comb2, **comb3, **comb4;
	///////
	int *vedosg2, *vedosg3, *vedosg4;
	bool yaprecosg;
	bool yasgcat;
	
	int **nompremot_2, **nompremot_3, **nompremot_4;
	int numnompremot_2, numnompremot_3, numnompremot_4;
	int *socont_2, *socont_3, *socont_4;
	
	//
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
	int *degdist;
	int *odegdist;
	int *idegdist;
	bool yadegdist;
  
  
	//functions
	void set_default_exclusive_vars();
	int calc_indegree(int no);
	int calc_outdegree(int no);
	int calc_degree(int no);
	void prepare_for_degrees();
	void prepare_for_degdist();
	void priv_get_degdist(); // esta y la anterior siempre se deben llamar juntas

  void aswitch();
  void aswitch_preserving_sg1();
  void aswitch_preserving_sg2();

  
  //modularity


	//exclusive of graphb:
	int count_this_subgraph2(GraphB &sg);
	int count_this_subgraph3(GraphB &sg);
	int count_this_subgraph4(GraphB &sg);
	int loose_count_this_subgraph2(GraphB &sg);
	int loose_count_this_subgraph3(GraphB &sg);
	int loose_count_this_subgraph4(GraphB &sg);
	int list_subgraph_appearances2(GraphB &sg, string arch);
	int list_subgraph_appearances3(GraphB &sg, string arch);
	int list_subgraph_appearances4(GraphB &sg, string arch);
	int list_subgraph_appearances2(GraphB &sg, ostream& os);
	int list_subgraph_appearances3(GraphB &sg, ostream& os);
	int list_subgraph_appearances4(GraphB &sg, ostream& os);
  
  
	int list_subgraph_appearances2(GraphB &sg, string arch, bool **me, int ge, int tis);
	int list_subgraph_appearances3(GraphB &sg, string arch, bool **me, int ge, int tis);
	int list_subgraph_appearances4(GraphB &sg, string arch, bool **me, int ge, int tis);
	int list_subgraph_appearances2(GraphB &sg, ostream& os, bool **me, int ge, int tis);
	int list_subgraph_appearances3(GraphB &sg, ostream& os, bool **me, int ge, int tis);
	int list_subgraph_appearances4(GraphB &sg, ostream& os, bool **me, int ge, int tis);
  
	int list_subgraph_appearances_using_names2(GraphB &sg, string arch);
	int list_subgraph_appearances_using_names3(GraphB &sg, string arch);
	int list_subgraph_appearances_using_names4(GraphB &sg, string arch);
	int list_subgraph_appearances_using_names2(GraphB &sg, ostream& os);
	int list_subgraph_appearances_using_names3(GraphB &sg, ostream& os);
	int list_subgraph_appearances_using_names4(GraphB &sg, ostream& os);
  
	int list_subgraph_appearances_using_names2(GraphB &sg, string arch, bool **me, int ge, int tis);
	int list_subgraph_appearances_using_names3(GraphB &sg, string arch, bool **me, int ge, int tis);
	int list_subgraph_appearances_using_names4(GraphB &sg, string arch, bool **me, int ge, int tis);
	int list_subgraph_appearances_using_names2(GraphB &sg, ostream& os, bool **me, int ge, int tis);
	int list_subgraph_appearances_using_names3(GraphB &sg, ostream& os, bool **me, int ge, int tis);
	int list_subgraph_appearances_using_names4(GraphB &sg, ostream& os, bool **me, int ge, int tis);
  
	bool an_isomorph_contains_these_interactions1(GraphB &subg);
	bool an_isomorph_contains_these_interactions2(GraphB &subg);
	bool an_isomorph_contains_these_interactions3(GraphB &subg);
	bool an_isomorph_contains_these_interactions4(GraphB &subg);
	bool isomorphic_pair1(GraphB &templ, int *ord);
	bool isomorphic_pair2(GraphB &templ, int *ord);
	bool isomorphic_pair3(GraphB &templ, int *ord);
	bool isomorphic_pair4(GraphB &templ, int *ord);
	bool isomorphic_pair1(GraphB &templ);
	bool isomorphic_pair2(GraphB &templ);
	bool isomorphic_pair3(GraphB &templ);
	bool isomorphic_pair4(GraphB &templ);
  
  
  
};


#endif
