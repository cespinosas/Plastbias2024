#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <set>
#include <list>
#include "basics.h"
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

GraphC::GraphC()
{
  ionary = false;
}

GraphC::GraphC(Alea& jacta)
{
  start_rng(jacta);
  ionary = false;
}

void GraphC::start_rng(Alea& jacta)
{
  est = jacta;
  basic.start_rng(est);
}

GraphC::GraphC(int n, Alea& jacta, bool dirnw)
{
  start_rng(jacta);
  make_nw(n, dirnw);
}

void GraphC::make_nw(int n, bool dirnw)
{
  int i;
  size = n;
  nw = new double*[size];
  for (i = 0; i < size; i++)
		nw[i] = new double[size];
  basic.fillmat0(nw, size, size);
  ionary = false;
  directed = dirnw;
  set_default_exclusive_vars();
  edo = new int[size];
  ima = new int[size];
}

void GraphC::make_robust_nw(int n, int e, int *cia, int *pa, int *cib, int *pb) {
  int nurcom=0,nurdif=0,comd=0,difd=0; //counter for genes in each class
  int i,j,l;
  for (i=0; i < n; i++) {
    if (((cia[i]*cib[i])>0) && ((pa[i]*pb[i])>0) && ((cia[i]*pa[i])>0))
      nurcom++;
    if (((cia[i]*cib[i])<0) && ((pa[i]*pb[i])<0) && ((cia[i]*pa[i])>0))
      nurdif++;
    if (pa[i]==pb[i])
      comd++;
    else
      difd++;
  }
  if ((nurcom==0) || ((nurdif==0)&&(difd>0))) {
    cout << "[Error]: Impossible network for GraphC::make_robust_nw.\n";
    exit(1);
  }
  int *vrcom, *vrdif, mc=0, md=0, *igudown, *difdown, cid=0, cdd = 0;
  vrcom = new int[nurcom];
  vrdif = new int[nurdif];
  igudown = new int[comd];
  difdown = new int[difd];
  for (i=0; i < n; i++) {
    if (((cia[i]*cib[i])>0) && ((pa[i]*pb[i])>0) && ((cia[i]*pa[i])>0)) {
      vrcom[mc] = i;
      mc++;
    }
    if (((cia[i]*cib[i])<0) && ((pa[i]*pb[i])<0) && ((cia[i]*pa[i])>0)) {
      vrdif[md] = i;
      md++;
    }
    if ((pa[i]*pb[i])>0) {
      igudown[cid] = i;
      cid++;
    }
    if ((pa[i]*pb[i])<0) {
      difdown[cdd] = i;
      cdd++;
    }
  }
  int maxd = nurdif*difd;
  int maxc = nurcom*comd;
  int k;//,minmax, maxmax;
  bool usobc;
  if (maxd > maxc)
    usobc = false;
  else
    usobc = true;
  
  if ((comd != n) && (difd != n)) {
    if (nurdif < nurcom)
      k = nurdif;
    else
      k = nurcom;
  } else {
    if (nurdif==0)
      k = nurcom;
    else
      k = nurdif;
  }
  while (e < (k*n)) {
    k--;
  }
  if (k < 1) {
    cout << "[Error]: Impossible network for GraphC::make_robust_nw.\n";
    exit(1);
  }
  int val, maxinbue, vanin, conta=0;
  bool yast = false;
  do {
    make_nw(n, true);
    for (l = 0; l < k; l++) {
      for (i=0; i < n; i++) {
        if ((pa[i]*pb[i])>0) {
          do {
            j = est.randint(0,nurcom);
          } while (weight(vrcom[j],i) != 0);
          if ((pa[vrcom[j]]*pa[i])>0)
            val = 1;
          else
            val = -1;
          force_interaction(vrcom[j],i,val);
        } else {
          do {
            j = est.randint(0,nurdif);
          } while (weight(vrdif[j],i) != 0);
          if ((pa[vrdif[j]]*pa[i])>0)
            val = 1;
          else
            val = -1;
          force_interaction(vrdif[j],i,val);
        }
      }
    }
    vanin = k*n;
    if (e < (maxd+maxc))
      maxinbue = e;
    else
      maxinbue = maxc+maxd;
    if (usobc) {
      for (i = vanin; i < maxinbue; i++) {
        do {
          j = est.randint(0,nurcom);
          l = est.randint(0,comd);
        } while (weight(vrcom[j],igudown[l]) != 0) ;
        if ((pa[vrcom[j]]*pa[igudown[l]]) > 0)
          val = 1;
        else
          val = -1;
        force_interaction(vrcom[j],igudown[l],val);
      }
    } else {
      for (i = vanin; i < maxinbue; i++) {
        do {
          j = est.randint(0,nurdif);
          l = est.randint(0,difd);
        } while (weight(vrdif[j],difdown[l]) != 0) ;
        if ((pa[vrdif[j]]*pa[difdown[l]]) > 0)
          val = 1;
        else
          val = -1;
        force_interaction(vrdif[j],difdown[l],val);
      }
    }
    vanin = maxinbue;
    for (i = vanin; i < e; i++) {
      do {
        j = est.randint(0, n);
        l = est.randint(0, n);
      } while (weight(j,l) != 0);
      if (est.toss())
        val = 1;
      else
        val = -1;
      force_interaction(j,l,val);
    }
    yast = true;
    set_as_state(cia);
    find_an_attractor();
    if (attractor_size() != 1)
      yast = false;
    else {
      if (!basic.eqvec(pa, n, attractor[0], n))
        yast = false;
      else {
        clear_attractor();
        set_as_state(cib);
        find_an_attractor();
        if (attractor_size() != 1)
          yast = false;
        else
          if (!basic.eqvec(pb, n, attractor[0], n))
            yast = false;
        clear_attractor();
      }
    }
    if (!yast)
      clear();
    cout << ++conta << endl;
    if (conta > 200) {
      cout << "[Error]: GraphC::make_robust_nw surpassed limit number of tries.\n";
      exit(1);
    }
  }while(!yast);
  delete [] vrcom;
  delete [] vrdif;
  delete [] igudown;
  delete [] difdown;
}

void GraphC::copy(GraphC &templ)
{
  int n = templ.number_of_nodes();
	bool dirnw = templ.is_directed();
  make_nw(n, dirnw);
  int i, j;
  for (i = 0; i < size; i++)
    for (j = 0; j < size; j++)
      nw[i][j] = templ.weight(j, i);
}

void GraphC::make_copy(GraphC& vacia)
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

void GraphC::put_subgraph_into(GraphC& nueva, set<int> &cuales)
{
	set<int> cua2 = cuales;
	nueva.make_nw(cuales.size(), is_directed());
	set<int>::iterator it, it2;
	int i,j;
	i=0;
	for (it=cuales.begin(); it!=cuales.end(); it++) {
		j=0;
		for (it2=cua2.begin(); it2!=cua2.end(); it2++) {
			nueva.force_interaction(j, i, weight(*it2, *it));
			j++;
		}
		i++;
	}
}

void GraphC::put_subgraph_into(GraphC& nueva, set<int> &cuales, bool withnames)
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

void GraphC::make_subgraph_of(GraphC &templ, set<int> &cuales)
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

void GraphC::make_subgraph_of(GraphC &templ, set<int> &cuales, bool withnames)
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

void GraphC::copy_wdeg(GraphC &templ)
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

GraphC::GraphC(int n, double c, Alea& jacta, bool dirnw)
{
  start_rng(jacta);
  make_nw(n, c, dirnw);
}

void GraphC::make_nw(int n, double c, bool dirnw)
{
  if (c >= 1) {
		cout << "[Error]: You are trying to construct an overconnected network. GraphC::make_nw.\n";
		exit(1);
	}
  int i, j;
  double nure,newval;
  make_nw(n, dirnw);
	if (directed) {
		for (i = 0; i < size; i++)
			for(j = 0; j < size; j++) {
				nure = est.randreal();
				if (nure <= c) {
          newval = est.randgauss();
          change_interaction(j,i,newval);
				}
			}
	}
	else {
		for (i=0; i<size; i++)
			for (j=0; j<=i; j++) {
				nure = est.randreal();
				if (nure <= c) {
          newval = est.randgauss();
          change_interaction_undir(j,i,newval);
				}
			}
	}
}

GraphC::GraphC(int n, int e, Alea& jacta, bool dirnw)
{
  start_rng(jacta);
  make_nw(n, e, dirnw);
}

void GraphC::make_nw(int n, int e, bool dirnw)
{
  if (e > (n*n)) {
		cout << "[Error]: You are trying to construct an overconnected network. GraphC::make_nw.\n";
		exit(1);
	}
  int i, j, k;
  double newval;
  make_nw(n, dirnw);
  for (i = 0; i < e; i++) {
		do {
			j = est.randint(0, size);
			k = est.randint(0, size);
		} while (nw[j][k] != 0);
    newval = est.randgauss();
    if (directed)
      change_interaction(k,j,newval);
		else
      change_interaction_undir(j,k,newval);
	}
}

void GraphC::rwalk_1traj(int n, int e, int *fp, int *ic, int cuan, string nomarch, int mame, int numut, bool estri, int numarch, int archdesde) {
  make_robust_nw(n, e, ic, fp, ic, fp); //nuevo

  int h,i,j,k,l,sig, vval, numutr,mr, enarch;
  ofstream fs;
  if ((cuan%numarch)!=0) {
    cout << "[Error]:  Files of different sizes in GraphC::rwalk_1fp.\n";
    exit(1);
  }
  enarch = cuan/numarch;
  bool sube, estasi;
  for (h = 0; h < numarch; h++) {
    basic.open_ofstream(fs, nomarch+basic.inttostring(h+archdesde)+".txt");
    for (i=0; i<enarch; i++) {
      do {
        mr = n*n*numut;
        numutr = est.randint(int(mr*8/10), int(mr*12/10));
        for (j=0; j<numutr; j++) {
          if (number_of_edges() == e+mame)
            sube = false;
          else {
            if (number_of_edges() == e-mame)
              sube = true;
            else {
              sube = false;
              if (est.toss())
                sube = true;
            }
          }
          if (!sube) {
            do {
              k = est.randint(0,size);
              l = est.randint(0, size);
            }while(weight(k,l)==0);
            vval = weight(k,l);
            force_interaction(k,l,0);
            set_as_state(ic);
            find_an_attractor();
            if ((attractor_size() != 1) || (!basic.eqvec(attractor[0], size, fp, size)))
              force_interaction(k,l,vval);
          } else {
            do {
              k = est.randint(0,size);
              l = est.randint(0, size);
            }while(weight(k,l)!=0);
            vval = weight(k,l);
            sig = est.randgauss();
            force_interaction(k,l,sig);
            set_as_state(ic);
            find_an_attractor();
            if ((attractor_size() != 1) || (!basic.eqvec(attractor[0], size, fp, size)))
              force_interaction(k,l,vval);
          }
        }
        estasi = true;
        if (estri)
          if (number_of_edges() != e)
            estasi = false;
      } while(!estasi);
      export_nw(fs);
      fs << "*\n";
    }
    fs.close();
  }
}

void GraphC::rwalk_2fp1t(int n, int e, int *fp1, int *fp2, int *ic1, int cuan, string nomarch, int mame, int numut, bool estri, int numarch, int archdesde) {
  make_robust_nw(n, e, ic1, fp1, fp2, fp2); //nuevo
  int h,i,j,k,l,sig, vval, numutr,mr, enarch;
  ofstream fs;
  if ((cuan%numarch)!=0) {
    cout << "[Error]:  Files of different sizes in GraphC::rwalk_2fp1t.\n";
    exit(1);
  }
  enarch = cuan/numarch;
  bool sube, estasi;
  for (h = 0; h < numarch; h++) {
    basic.open_ofstream(fs, nomarch+basic.inttostring(h+archdesde)+".txt");
    for (i=0; i<enarch; i++) {
      do {
        mr = n*n*numut;
        numutr = est.randint(int(mr*8/10), int(mr*12/10));
        for (j=0; j<numutr; j++) {
          if (number_of_edges() == e+mame)
            sube = false;
          else {
            if (number_of_edges() == e-mame)
              sube = true;
            else {
              sube = false;
              if (est.toss())
                sube = true;
            }
          }
          if (!sube) {
            do {
              k = est.randint(0,size);
              l = est.randint(0, size);
            }while(weight(k,l)==0);
            vval = weight(k,l);
            force_interaction(k,l,0);
            set_as_state(ic1);
            find_an_attractor();
            if ((attractor_size() != 1) || (!basic.eqvec(attractor[0], size, fp1, size)))
              force_interaction(k,l,vval);
            else {
              clear_attractor();
              set_as_state(fp2);
              find_an_attractor();
              if ((attractor_size() != 1) || (!basic.eqvec(attractor[0], size, fp2, size)))
                force_interaction(k,l,vval);
              else
                clear_attractor();
            }
          } else {
            do {
              k = est.randint(0,size);
              l = est.randint(0, size);
            }while(weight(k,l)!=0);
            vval = weight(k,l);
            sig = est.randgauss();
            force_interaction(k,l,sig);
            set_as_state(ic1);
            find_an_attractor();
            if ((attractor_size() != 1) || (!basic.eqvec(attractor[0], size, fp1, size)))
              force_interaction(k,l,vval);
            else {
              clear_attractor();
              set_as_state(fp2);
              find_an_attractor();
              if ((attractor_size() != 1) || (!basic.eqvec(attractor[0], size, fp2, size)))
                force_interaction(k,l,vval);
              else
                clear_attractor();
            }
          }
        }
        estasi = true;
        if (estri)
          if (number_of_edges() != e)
            estasi = false;
      } while(!estasi);
      export_nw(fs);
      fs << "*\n";
    }
    fs.close();
  }
}

void GraphC::rwalk_2fp1tn1(int n, int e, int *fp1, int *fp2, int *ic1, int cuan, string nomarch, int mame, int numut, bool estri, int numarch, int archdesde) {
  int h,i,j,k,l,sig, vval, numutr,mr, enarch,ii;
  int *ic2;
  ic2 = new int[n];
  for (ii=0; ii<n; ii++) //nuevo
    ic2[ii] = ic1[ii];
  for (i=0; i < n; i++) {
    if ((ic1[i]==fp1[i]) && (fp1[i] != fp2[i])) {
      ic2[i] *= (-1);
      break;
    }
  }
  if (i==n) {
    cout << "[Error]: Impossible initial network for GraphC::rwalk_2fp1tn1.\n";
    exit(1);
  }
  make_robust_nw(n,e, ic1, fp1, ic2, fp2);//hastaaqui
  ofstream fs;
  bool ucon;
  if ((cuan%numarch)!=0) {
    cout << "[Error]:  Files of different sizes in GraphC::rwalk_2fp1tn1.\n";
    exit(1);
  }
  enarch = cuan/numarch;
  bool sube, estasi;
  for (h = 0; h < numarch; h++) {
    basic.open_ofstream(fs, nomarch+basic.inttostring(h+archdesde)+".txt");
    for (i=0; i<enarch; i++) {
      do {
        mr = n*n*numut;
        numutr = est.randint(int(mr*8/10), int(mr*12/10));
        for (j=0; j<numutr; j++) {
          if (number_of_edges() == e+mame)
            sube = false;
          else {
            if (number_of_edges() == e-mame)
              sube = true;
            else {
              sube = false;
              if (est.toss())
                sube = true;
            }
          }
          if (!sube) {
            do {
              k = est.randint(0,size);
              l = est.randint(0, size);
            }while(weight(k,l)==0);
            vval = weight(k,l);
            force_interaction(k,l,0);
            set_as_state(ic1);
            find_an_attractor();
            if ((attractor_size() != 1) || (!basic.eqvec(attractor[0], size, fp1, size)))
              force_interaction(k,l,vval);
            else {
              clear_attractor();
              for (ii=0; ii<size; ii++)
                ic2[ii] = ic1[ii];
              ucon = false;
              for (ii=0; ii<size; ii++) {
                ic2[ii]*=(-1);
                set_as_state(ic2);
                find_an_attractor();
                ic2[ii]*=(-1);
                if ((attractor_size() == 1) && (basic.eqvec(attractor[0], size, fp2, size))) {
                  clear_attractor();
                  ucon = true;
                  break;
                }
                clear_attractor();
              }
              if (!ucon)
                force_interaction(k,l,vval);
            }
          } else {
            do {
              k = est.randint(0,size);
              l = est.randint(0, size);
            }while(weight(k,l)!=0);
            vval = weight(k,l);
            sig = est.randgauss();
            force_interaction(k,l,sig);
            set_as_state(ic1);
            find_an_attractor();
            if ((attractor_size() != 1) || (!basic.eqvec(attractor[0], size, fp1, size)))
              force_interaction(k,l,vval);
            else {
              clear_attractor();
              for (ii=0; ii<size; ii++)
                ic2[ii] = ic1[ii];
              ucon = false;
              for (ii=0; ii<size; ii++) {
                ic2[ii]*=(-1);
                set_as_state(ic2);
                find_an_attractor();
                ic2[ii]*=(-1);
                if ((attractor_size() == 1) && (basic.eqvec(attractor[0], size, fp2, size))) {
                  clear_attractor();
                  ucon = true;
                  break;
                }
                clear_attractor();
              }
              if (!ucon)
                force_interaction(k,l,vval);
            }
          }
        }
        estasi = true;
        if (estri)
          if (number_of_edges() != e)
            estasi = false;
      } while(!estasi);
      export_nw(fs);
      fs << "*\n";
    }
    fs.close();
  }
  delete [] ic2;
}

void GraphC::rwalk_2fp1tn2(int n, int e, int *fp1, int *fp2, int *ic1, int cuan, string nomarch, int mame, int numut, bool estri, int numarch, int archdesde) {
  int h,i,j,k,l,sig, vval, numutr,mr, enarch,ii,jj;
  int *ic2;
  ic2 = new int[n];
  int uy = 0;//nuevo
  for (ii=0; ii<n; ii++)
    ic2[ii] = ic1[ii];
  for (i=0; i < n; i++) {
    if ((ic1[i]==fp1[i]) && (fp1[i] != fp2[i])) {
      ic2[i] *= (-1);
      if (uy == 0)
        uy++;
      if (uy > 0)
        break;
    }
  }
  if (i==n) {
    cout << "[Error]: Impossible initial network for GraphI::rwalk:2fp1tn2.\n";
    exit(1);
  }
  make_robust_nw(n,e, ic1, fp1, ic2, fp2);//hastaaqui
  ofstream fs;
  bool ucon;
  if ((cuan%numarch)!=0) {
    cout << "[Error]:  Files of different sizes in GraphC::rwalk_2fp1tn2.\n";
    exit(1);
  }
  enarch = cuan/numarch;
  bool sube, estasi;
  for (h = 0; h < numarch; h++) {
    basic.open_ofstream(fs, nomarch+basic.inttostring(h+archdesde)+".txt");
    for (i=0; i<enarch; i++) {
      do {
        mr = n*n*numut;
        numutr = est.randint(int(mr*8/10), int(mr*12/10));
        for (j=0; j<numutr; j++) {
          if (number_of_edges() == e+mame)
            sube = false;
          else {
            if (number_of_edges() == e-mame)
              sube = true;
            else {
              sube = false;
              if (est.toss())
                sube = true;
            }
          }
          if (!sube) {
            do {
              k = est.randint(0,size);
              l = est.randint(0, size);
            }while(weight(k,l)==0);
            vval = weight(k,l);
            force_interaction(k,l,0);
            set_as_state(ic1);
            find_an_attractor();
            if ((attractor_size() != 1) || (!basic.eqvec(attractor[0], size, fp1, size)))
              force_interaction(k,l,vval);
            else {
              clear_attractor();
              for (ii=0; ii<size; ii++)
                ic2[ii] = ic1[ii];
              ucon = false;
              for (ii=0; ii<size; ii++) {
                ic2[ii]*=(-1);
                for (jj = ii; jj < size; jj++) {
                  if (jj != ii)
                    ic2[jj] *= (-1);
                  set_as_state(ic2);
                  find_an_attractor();
                  if (jj != ii)
                    ic2[jj] *= (-1);
                  if ((attractor_size() == 1) && (basic.eqvec(attractor[0], size, fp2, size))) {
                    clear_attractor();
                    ucon = true;
                    break;
                  }
                  clear_attractor();
                }
                ic2[ii]*=(-1);
                if (ucon)
                  break;
              }
              if (!ucon)
                force_interaction(k,l,vval);
            }
          } else {
            do {
              k = est.randint(0,size);
              l = est.randint(0, size);
            }while(weight(k,l)!=0);
            vval = weight(k,l);
            sig = est.randgauss();
            force_interaction(k,l,sig);
            set_as_state(ic1);
            find_an_attractor();
            if ((attractor_size() != 1) || (!basic.eqvec(attractor[0], size, fp1, size)))
              force_interaction(k,l,vval);
            else {
              clear_attractor();
              for (ii=0; ii<size; ii++)
                ic2[ii] = ic1[ii];
              ucon = false;
              for (ii=0; ii<size; ii++) {
                ic2[ii]*=(-1);
                for (jj = ii; jj < size; jj++) {
                  if (jj != ii)
                    ic2[jj] *= (-1);
                  set_as_state(ic2);
                  find_an_attractor();
                  if (jj != ii)
                    ic2[jj] *= (-1);
                  if ((attractor_size() == 1) && (basic.eqvec(attractor[0], size, fp2, size))) {
                    clear_attractor();
                    ucon = true;
                    break;
                  }
                  clear_attractor();
                }
                ic2[ii]*=(-1);
                if (ucon)
                  break;
              }
              if (!ucon)
                force_interaction(k,l,vval);
            }
          }
        }
        estasi = true;
        if (estri)
          if (number_of_edges() != e)
            estasi = false;
      } while(!estasi);
      export_nw(fs);
      fs << "*\n";
    }
    fs.close();
  }
  delete [] ic2;
}

void GraphC::clear()
{
  int i;
  for (i = 0; i < size; i++)
    delete [] nw[i];
  delete [] nw;
  delete [] ima;
  delete [] edo;
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
  if (yatra)
    clear_attractor();
	if (yae)
		clear_edge_count();
  if (yatras)
    clear_attractors();
  if (yamadya)
    clear_adj();
}

void GraphC::clear_edge_count()
{
	numofe = 0;
	yae = false;
}

void GraphC::clear_loops()
{
  delete [] loops;
  totalloops = 0;
  int i;
  for (i = 0; i < size; i++)
    delete [] guada[i];
  delete [] guada;
  floops = false;
}

void GraphC::clear_comp()
{
  delete [] components;
  yacomp = false;
}

void GraphC::clear_bc()
{
  delete [] bc;
  bcya = false;
}

void GraphC::clear_scc()
{
  delete [] scc;
  cfc = false;
  numbscc = 0;
}

void GraphC::clear_adj() {
  int i;
  for (i = 0; i < size; i++)
    delete [] matadya[i];
  delete [] matadya;
  yamadya = false;
}

void GraphC::clear_attractor()
{
  int i;
  for (i = 0; i < atsize; i++)
    delete [] attractor[i];
  delete [] attractor;
  atsize = 0;
  palen = 0;
  yatra = false;
}

void GraphC::clear_attractors()
{
  int i,j;
  for (i=0; i < numatrs; i++) {
    for (j=0; j < atsizes[i]; j++)
      delete [] attractors[i][j];
    delete [] attractors[i];
  }
  delete [] attractors;
  delete [] atsizes;
  numatrs = 0;
  yatras = false;
}

void GraphC::clear_dima()
{
  int i;
  for (i = 0; i < size; i++)
    delete [] dima[i];
  delete [] dima;
  dist = false;
}

void GraphC::clear_ds()
{
  int i;
  for (i = 0; i < size; i++)
    downstream[i].clear();
  delete [] downstream;
  infl = false;
}

void GraphC::clear_dict()
{
  delete [] dict;
  ionary = false;
}

void GraphC::clear_deg()
{
  delete [] degree;
  delete [] outdegree;
  delete [] indegree;
  yadeg = false;
	if (yadegdist)
		clear_degdist();
}

void GraphC::clear_degdist()
{
	delete [] degdist;
	delete [] odegdist;
	delete [] idegdist;
	yadegdist = false;
}

void GraphC::get_all_degrees()
{
  if (yadeg) {
		cout << "[Error]: Attempting to create degree arrays again. GraphC::get_all_degrees.\n";
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

int GraphC::get_indegree(int no)
{
  if (!yadeg)
    get_all_degrees();
  return indegree[no];
}

int GraphC::get_outdegree(int no)
{
  if (!yadeg)
    get_all_degrees();
  return outdegree[no];
}

int GraphC::get_degree(int no)
{
  if (!yadeg)
    get_all_degrees();
  return degree[no];
}

void GraphC::get_degdist()
{
	if (!yadegdist) {
		prepare_for_degdist();
		priv_get_degdist();
		yadegdist = true;
	}
}

void GraphC::get_degdist(int *empvec, int tam)
{
	if (tam != (size+1)){
		cout << "[Error]: Wrong size of array in GraphC::get_degdist.\n";
		exit(1);
	}
	int i;
	get_degdist();
	for (i=0; i<tam; i++)
		empvec[i] = degdist[i];
}

void GraphC::get_odegdist(int *empvec, int tam)
{
	if (tam != (size+1)){
		cout << "[Error]: Wrong size of array in GraphC::get_degdist.\n";
		exit(1);
	}
	int i;
	get_degdist();
	for (i=0; i<tam; i++)
		empvec[i] = odegdist[i];
}

void GraphC::get_idegdist(int *empvec, int tam)
{
	if (tam != (size+1)){
		cout << "[Error]: Wrong size of array in GraphC::get_degdist.\n";
		exit(1);
	}
	int i;
	get_degdist();
	for (i=0; i<tam; i++)
		empvec[i] = idegdist[i];
}

//access
int GraphC::number_of_edges()
{
  if (!yae) {
		numofe = 0;
		int i, j;
		if (directed) {
			for (i = 0; i < size; i++)
				for (j = 0; j < size; j++)
					if (nw[i][j] != 0)
						numofe++;
		}
		else {
			for (i = 0; i < size; i++)
				for (j = 0; j <= i; j++)
					if (nw[i][j] != 0)
						numofe++;
		}
		yae = true;
	}
  return numofe;
}

int GraphC::number_of_nodes()
{
  return size;
}

double GraphC::weight(int source, int target)
{
  return nw[target][source];
}

bool GraphC::names_exist()
{
  return ionary;
}

string GraphC::get_name(int node)
{
  return dict[node];
}

bool GraphC::is_directed()
{
  return directed;
}

bool GraphC::already_deg()
{
  return yadeg;
}

int GraphC::name_to_int(string cual)
{
	int res=-1, i;
	if (!names_exist()) {
		cout << "[Error]: Names do not exist. GraphC::name_to_int.\n";
		exit(1);
	}
	for (i=0; i<size; i++)
		if (dict[i] == cual) {
			res = i;
			break;
		}
	return res;
}

int GraphC::name_to_int(string cual, bool checkornot)
{
	int res = name_to_int(cual);
	if ((checkornot) && (res==-1)) {
		cout << "[Error]: There is no node with the name \'" << cual << "\'. GraphC::name_to_int.\n";
		exit(1);
	}
	return res;
}

bool GraphC::is_name(string elna)
{
	bool res = false;
	int i;
	if (!names_exist()) {
		cout << "[Error]: Names do not exist when GraphC::is_name was called.\n";
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

bool GraphC::attractor_exists()
{
  return yatra;
}

bool GraphC::attractors_exist()
{
  return yatras;
}

//modification
void GraphC::force_interaction(int source, int target, double value)
{
	if (directed) {
		if (yae) {
      if ((value == 0) && (nw[target][source] != 0))
				numofe--;
      else
				if ((value != 0) && (nw[target][source] == 0))
					numofe++;
		}
		if (yadeg) {
      if ((value == 0) && (nw[target][source] != 0)) {
				outdegree[source]--;
				indegree[target]--;
				degree[source]--;
				degree[target]--;
			}
      else {
				if ((value != 0) && (nw[target][source] == 0)) {
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
		cout << "[Error]: force_interaction does not work for undirected graphs. GraphC::force_interaction.\n";
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
  if (yatra)
    clear_attractor();
  if (yatras)
    clear_attractors();
  if (yadegdist)
      clear_degdist();
}

void GraphC::force_interaction_undir(int source, int target, double value)
{
	if (directed) {
		cout << "[Error]: force_interaction_undir does not work for directed graphs. GraphC::force_interaction_undir.\n";
		exit(1);
	}
	if (yae) {
		if ((value==0) && (nw[target][source] != 0))
			numofe--;
		else
			if ((value != 0) && (nw[target][source] == 0))
				numofe++;
	}
	if (yadeg) {
		if ((value==0) && (nw[target][source] != 0)) {
			outdegree[source]--;
			outdegree[target]--;
			indegree[target]--;
			indegree[source]--;
			degree[source]--;
			degree[target]--;
		}
		else
			if ((value!=0) && (nw[target][source] == 0)) {
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
  if (yatra)
    clear_attractor();
  if (yatras)
    clear_attractors();
  if (yadegdist)
    clear_degdist();
}

void GraphC::change_interaction(int source, int target, double value)
{
	if (directed) {
		if (value == nw[target][source]) {
			cout << "[Error]: GraphC::change_interaction changes to stay the same.\n";
			exit(1);
		}
		if (yae) {
      if ((value == 0) && (nw[target][source] != 0))
				numofe--;
      else
				if ((value != 0) && (nw[target][source] == 0))
					numofe++;
		}
		if (yadeg) {
      if ((value == 0) && (nw[target][source] != 0)) {
				outdegree[source]--;
				indegree[target]--;
				degree[source]--;
				degree[target]--;
			}
      else {
				if ((value != 0) && (nw[target][source] == 0)) {
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
		cout << "[Error]: GraphC::change_interaction does not work for undirected graphs.\n";
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
  if (yatra)
    clear_attractor();
  if (yatras)
    clear_attractors();
  if (yadegdist)
    clear_degdist();
}

void GraphC::change_interaction_undir(int source, int target, double value)
{
	if (directed) {
		cout << "[Error]: GraphC::change_interaction_undir does not work for directed graphs\n";
		exit(1);
	}
	if (value == nw[target][source]) {
		cout << "[Error]: GraphC::change_interaction changes to stay the same.";
		exit(1);
	}
	if (yae) {
		if ((value==0) && (nw[target][source] != 0))
			numofe--;
		else
			if ((value != 0) && (nw[target][source] == 0))
				numofe++;
	}
	if (yadeg) {
		if ((value==0) && (nw[target][source] != 0)) {
			outdegree[source]--;
			outdegree[target]--;
			indegree[target]--;
			indegree[source]--;
			degree[source]--;
			degree[target]--;
		}
		else
			if ((value!=0) && (nw[target][source] == 0)) {
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
  if (yatra)
    clear_attractor();
  if (yatras)
    clear_attractors();
  if (yadegdist)
    clear_degdist();
}

void GraphC::change_name(int node, string nn)
{
  if (!(ionary)) {
		cout << "[Error]: There are no names linked to this network. GraphC::change_name.\n";
		exit(1);
	}
  dict[node] = nn;
}

void GraphC::set_undirected()
{
  if (!(basic.is_symmetric(nw, size, size))) {
		cout << "[Error]: Attempt to set a non-symmetric network as undirected. GraphC::set_undirected.\n";
		exit(1);
	}
  else
    directed = false;
}

void GraphC::set_directed()
{
  directed = true;
}

void GraphC::force_undirected()
{
	int i,j;
	for (i=0; i<size; i++) {
		for (j=0; j< size; j++)
			if (weight(i,j) != weight(j,i)) {
        if (directed)
          change_interaction(j,i, weight(i,j));
        else
          change_interaction_undir(j,i, weight(i,j));
			}
	}
	if (yae)
		clear_edge_count();
	if (yadeg)
		clear_deg();
	set_undirected();
}
                                                                                                                                        
void GraphC::transform_to_isomorph(GraphC& vacia, int* vec)
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
void GraphC::get_names_from_file(string arch)
{
	ifstream sal;
	basic.open_ifstream(sal, arch);
	get_names_from_file(sal);
	sal.close();
}

void GraphC::get_names_from_file(istream& en)
{
  int i=0;
	if (ionary) { //
		cout << "[Error]: Names have already been assigned. GraphC::get_names_from_file.\n"; //
		exit(1);//
	}//
  dict = new string[size];
	while (en >> dict[i])
		i++;
	while (i<size) {
		dict[i] = "";
		i++;
	}
  ionary = true;
}

void GraphC::assign_names(string *nombres, int ent)
{
	int i;
	if (ionary) { //
		cout << "[Error]: Names have already been assigned. GraphC::assign_names.\n"; //
		exit(1);//
	}//
  dict = new string[size];
	if (ent > size) {
		cout << "[Error]: More names than nodes. GraphC::assign_names.";
		exit(1);
	}
	for (i=0; i<ent; i++)
		dict[i] = nombres[i];
	for (i=ent; i<size; i++) {
		dict[i]= "";
	}
	ionary=true;
}

void GraphC::printnw(ostream& sal)
{
  int i, j;
  for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++)
      sal << nw[i][j] << "\t";
		sal << endl;
	}
}

void GraphC::printnw_bin(ostream& sal)
{
  int i, j;
  for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++) {
			if (nw[i][j] == 0)
				sal << 0 << " ";
			else
				if (nw[i][j] >0 )
					sal << "+ ";
				else
					if (nw[i][j] < 0)
						sal << "- ";
					else {
						cout << "[Error]: error while printing network. GraphC::printnw.\n";
						exit(1);
					}
		}
		sal << endl;
	}
}

void GraphC::print_dot(string archi)
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
				if (weight(j,i) !=0)
					sal << "\t" << j << "  " << ke << "  " << i << ";\n";
	}
	else
		for (i = 0; i < size; i++)
			for (j = 0; j <= i; j++)
				if (weight(j,i) !=0)
					sal << "\t" << j << "  " << ke << "  " << i << ";\n";
	sal << "}\n";
	sal.close();
}

void GraphC::print_dot_wn(string archi)
{
  if (!(ionary)) {
		cout << "[Error]: There are no names linked to this network. GraphC::print_dot_wn.\n";
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
				if (weight(j,i) !=0)
					sal << "\t\"" << dict[j] << "\"  " << ke << "  \"" << dict[i] << "\";\n";
	}
	else
		for (i = 0; i < size; i++)
			for (j = 0; j <= i; j++)
				if (weight(j,i) !=0)
					sal << "\t\"" << dict[j] << "\"  " << ke << "  \"" << dict[i] << "\";\n";
	sal << "}\n";
	sal.close();
}

void GraphC::print_dot_wo_labels_circ(string archi, double radius)
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
	double radians[size], x[size], y[size], stera, PI=3.141592;
	stera = (2*PI)/(size*1.0);
	radians[0] = 0.75*PI;
	for (i=1; i<size; i++) {
		radians[i] = radians[i-1]+stera;
		if (radians[i] >= (2*PI))
			radians[i] = radians[i] - (2*PI);
	}
	for (i=0; i<size; i++) {
		basic.polar_to_cartesian(radians[i], radius, x[i], y[i]);
		sal << i << "[pos = \"" << x[i] << "," << y[i] << "!\"];\n";
	}
	if (directed) {
		for (i=0; i < size; i++)
			for (j=0; j < size; j++) {
				if ((weight(j,i)) != 0) {
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
				if ((weight(j,i)) != 0) {
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

void GraphC::print_dot_circ(string archi, double radius)
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
  double radians[size], x[size], y[size], stera, PI=3.141592;
  stera = (2*PI)/(size*1.0);
  radians[0] = 0.75*PI;
  for (i=1; i<size; i++) {
    radians[i] = radians[i-1]+stera;
    if (radians[i] >= (2*PI))
      radians[i] = radians[i] - (2*PI);
  }
  for (i=0; i<size; i++) {
    basic.polar_to_cartesian(radians[i], radius, x[i], y[i]);
    sal << i << "[pos = \"" << x[i] << "," << y[i] << "!\",label=\"" << i << "\"];\n";
  }
  
  if (directed) {
    for (i=0; i < size; i++)
      for (j=0; j < size; j++) {
        if ((weight(j,i)) != 0) {
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
        if ((weight(j,i)) != 0) {
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

void GraphC::print_cytoscape(string archi)
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
        if ((weight(j,i) != 0) || (weight(i,j)!=0)) {
          sal << i << "\t" ;
          if (weight(j,i) > 0)
            sal << "a\t";
          else
            sal << "r\t";
          sal << j << "\t" << weight(j,i) << endl;
        }
  }
  else {
    for (i=0; i < size; i++)
      for (j=0; j < size; j++)
        if (weight(j,i)) {
          sal << j << "\t";
          if (weight(j,i) > 0)
            sal << "a\t";
          else
            sal << "r\t";
          sal << i << "\t" << weight(j,i) << endl;
        }
  }
  sal.close();
}

void GraphC::print_cytoscape_wn(string archi)
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
          if ((weight(j,i) != 0) || (weight(i,j)!=0)) {
            sal << dict[i] << "\t" ;
            if (weight(j,i) > 0)
              sal << "a\t";
            else
              sal << "r\t";
            sal << dict[j] << "\t" << weight(j,i) << endl;
          }
    }
    else {
      for (i=0; i < size; i++)
        for (j=0; j < size; j++)
          if (weight(j,i)) {
            sal << dict[j] << "\t";
            if (weight(j,i) > 0)
              sal << "a\t";
            else
              sal << "r\t";
            sal << dict[i] << "\t" << weight(j,i) << endl;
          }
    }
    sal.close();
  }
}

//exclusive
//import
void GraphC::get_dir_nw_from_file(int nn, istream& en, int e)
{
  make_nw(nn, true);
  int i, j, l;
  double k;
  for (l = 0; l < e; l++) {
    en >> j;
    en >> i;
    en >> k;
    if (((k>0) || (k<0)) && (i < size) && (j < size) && (i >= 0) && (j >= 0))
      nw[i][j] = k;
    else {
      cout << "[Error]: This is not an int network. GraphC::get_dir_nw_from_file.\n";
      exit(1);
    }
  }
  set_default_exclusive_vars();
}

void GraphC::get_dir_nw_from_file(int nn, string arch)
{
	ifstream sal;
	basic.open_ifstream(sal, arch);
	get_dir_nw_from_file(nn, sal);
	sal.close();
}

void GraphC::get_dir_nw_from_file(int nn, istream& en)
{
  make_nw(nn, true);
  int i, j;
  while (en >> j) {
		en >> i;
		if ((i < size) && (j < size) && (i >= 0) && (j >= 0))
			en >> nw[i][j];
		else {
			cout << "[Error]: This is not a double network. GraphC::get_dir_nw_from_file.\n";
			exit(1);
		}
	}
  set_default_exclusive_vars();
}

void GraphC::get_undir_nw_from_file(int nn, string arch)
{
	ifstream sal;
	sal.open(arch.c_str());
	get_undir_nw_from_file(nn, sal);
	sal.close();
}

void GraphC::get_undir_nw_from_file(int nn, istream& en)
{
  make_nw(nn, false);
  int i, j;
  while (en >> j) {
		en >> i;
		if ((i < size) && (j < size) && (i >= 0) && (j >= 0)) {
			en >> nw[i][j];
			nw[j][i] = nw[i][j];
		}
		else {
			cout << "[Error]: This is not a double network. GraphC::get_undir_nw_from_file.\n";
			exit(1);
		}
	}
  set_default_exclusive_vars();
}

void GraphC::get_dir_nw_from_file_wn(int nn, string arch, string archnam)
{
	ifstream sal, salna;
	basic.open_ifstream(salna, archnam);
	basic.open_ifstream(sal, arch);
	get_dir_nw_from_file_wn(nn, sal, salna);
	sal.close();
	salna.close();
}

void GraphC::get_dir_nw_from_file_wn(int nn, istream& en, istream& enna)
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
			cout << "[Error]: " << j << " or " << i << " not in the graph's node list. GraphC::get_dir_nw_from_file_wn.\n";
			exit(1);
		}
		en >> nw[l][k];
	}
  set_default_exclusive_vars();
}

void GraphC::get_undir_nw_from_file_wn(int nn, string arch, string archnam)
{
	ifstream sal, salna;
	basic.open_ifstream(salna, archnam);
	basic.open_ifstream(sal, arch);
	get_undir_nw_from_file_wn(nn, sal, salna);
	sal.close();
	salna.close();
}

void GraphC::get_undir_nw_from_file_wn(int nn, istream& en, istream& enna)
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
			cout << "[Error]: " << j << " or " << i << " not in the graph's node list. GraphC::get_undir_nw_from_file_wn.\n";
			exit(1);
		}
		en >> nw[l][k];
		nw[k][l] = nw[l][k];
	}
  set_default_exclusive_vars();
}

//Export
void GraphC::export_nw(string arch) {
  ofstream fs;
  basic.open_ofstream(fs, arch);
  export_nw(fs);
  fs.close();
}

void GraphC::export_nw(ostream& fs) {
  int i, j;
  if (directed) {
    for (i = 0; i < size; i++)
      for (j=0; j < size; j++)
        if (weight(j,i) != 0)
          fs << j << "\t" << i << "\t" << weight(j, i) << endl;
  }
  else {
    for (i = 0; i < size; i++)
      for (j=0; j <= i; j++)
        if (weight(j,i) != 0)
          fs << j << "\t" << i << "\t" << weight(j, i) << endl;
  }
}

//status

//analyses
double GraphC::distance_from_nw(GraphC& otra) {
  double res = 0;
  int i, j;
  for (i=0; i < size; i++)
    for (j = 0; j < size; j++)
      res += fabs(weight(i,j) - otra.weight(i,j));
  return res;
}

bool GraphC::equal_adjmat(GraphC& otra) {
  bool res = true;
  int i, j;
  if ((otra.number_of_nodes() != size) || (otra.number_of_edges() != number_of_edges()) || (otra.is_directed() != is_directed()))
    res = false;
  else {
    for (i = 0; i < size; i++) {
      if (!res)
        break;
      for (j = 0; j < size; j++) {
        if (((weight(i,j) == 0) && (otra.weight(i, j) != 0)) || ((weight(i,j) != 0) && (otra.weight(i, j) == 0))) {
          res = false;
          break;
        }
      }
    }
  }
  return res;
}

bool GraphC::equal_nw(GraphC &templ)
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

set<int> GraphC::get_upstream(int no)
{
	int i;
	set<int> res;
	res.clear();
	for (i=0; i<size; i++)
		if (nw[no][i] != 0)
			res.insert(i);
	return res;
}

void GraphC::get_downstream()
{
  downstream = new set<int>[size];
  int i, j;
  for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++) {
			if (nw[j][i] != 0)
				downstream[i].insert(j);
		}
	}
  infl = true;
}


set<int> GraphC::get_out_component(int no) {
  set<int> res;
  res.clear();
  int i;
  for (i=0; i< size; i++)
    if (there_is_path(no, i))
      res.insert(i);
  return res;
}

set<int> GraphC::get_in_component(int no) {
  set<int> res;
  res.clear();
  int i;
  for (i=0; i< size; i++)
    if (there_is_path(i, no))
      res.insert(i);
  return res;
}

void GraphC::get_adjacency_matrix() {
  if (yamadya) {
    cout << "[Error]: Adjacency matrix already created when GraphC::get_adjacency_matrix was called.\n";
    exit(1);
  }
  basic.create_array(matadya, size, size);
  int i,j;
  for (i = 0; i < size; i++)
    for (j = 0; j < size; j++) {
      if (nw[i][j] != 0)
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

void GraphC::get_distance_matrix_ls()
{
  if (dist) {
    cout << "[Error]: Distance matrix already created. GraphC::get_distance_matrix.\n";
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
          cout << "[Error]: GraphC::get_distance_matrix.\n";
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
        cout << "[Error]: GraphC::get_distance_matrix.\n";
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
          cout << "[Error]: GraphC::get_distance_matrix.\n";
          exit(1);
        }
  }
}

void GraphC::get_distance_matrix() {
  if (dist) {
    cout << "[Error]: Distance matrix already created when GraphC::get_distance_matrix was called.\n";
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
      if (nw[i][j] != 0)
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

int GraphC::get_scc() {
  int i, j;
  if (!directed) {
    cout << "[Error]: No SCC's in undirected networks. This error arose during execution of GraphC::get_scc.\n";
    exit(1);
  }
  if (cfc) {
    cout << "[Error]: SCC's already obtained when GraphC::get_scc was called.\n";
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

int GraphC::number_of_sccs() {
  if (!cfc)
    get_scc();
  return numbscc;
}

int GraphC::number_of_components() {
  if (!yacomp)
    get_components();
  return numbcomp;
}

int GraphC::get_components()
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
					cout << "[Error]: GraphC::get_components.\n";
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

void GraphC::get_bc()
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

double GraphC::get_bc(int nodo)
{
  if (!bcya)
    get_bc();
  return bc[nodo];
}

int GraphC::nupaths(int j, int q, int *vec)
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
  if (nw[v][v] != 0) {
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

list<int> GraphC::GetSucc(int v, list<int> &omega, int q)
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
        if (nw[r][i] != 0)
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
    if (nw[i][v] != 0)
      temp.push_back(i);
  for (itd = temp.begin(); itd != temp.end(); itd++) {
    u = *itd;
    if (d[u] == 1)
      succ.push_back(u);
  }
  return succ;
}

void GraphC::get_loops()
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
			cout << "[Error]: GraphC::get_loops.\n";
			exit(1);
		}
		loops[i] = loops[i]/(i+1);
	}
  for (i = 0; i < size; i++)
    totalloops = totalloops + loops[i];
  floops = true;
}

//acces analyses
bool GraphC::there_is_path(int from, int to)
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

int GraphC::length_shortest_path(int from, int to)
{
  int res;
  if (!dist)
    get_distance_matrix();
  res = dima[to][from];
  return res;
}

set<int> GraphC::nodes_influenced_by(int n)
{
  if (!infl)
    get_downstream();
  return downstream[n];
}

int GraphC::in_which_component(int n)
{
  if (!yacomp)
    get_components();
  return components[n];
}

int GraphC::in_which_scc(int n)
{
  if (!cfc)
    get_scc();
  return scc[n];
}

int GraphC::number_of_paths(int from, int to, int len, int* vec)
{
  nupaths(from, to, vec);
  return vec[len-1];
}

int GraphC::number_of_loops()
{
  if (!floops)
    get_loops();
  return totalloops;
}
int GraphC::number_of_loops(int len)
{
  if (!floops)
    get_loops();
  return loops[len];
}

void GraphC::switches(int vec)
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
  if (yatra)
    clear_attractor();
  if (yatras)
    clear_attractors();
  int i;
  for (i = 0; i < vec; i++)
    aswitch();
}

void GraphC::synchrony()
{
  int i, j;
  double sum;
  for (i = 0; i < size; i++) {
		sum = 0;
		for (j = 0; j < size; j++) {
			if (weight(j,i) != 0)
        sum = sum + (edo[j]*weight(j,i));
		}
		if (sum > 0)
			ima[i] = 1;
		if (sum < 0)
			ima[i] = -1;
		if (sum == 0)
      ima[i] = edo[i]; //it was 0 before 2017
	}
  for (i = 0; i < size; i++)
    edo[i] = ima[i];
}

void GraphC::set_as_state(int *vec)
{
  int i;
  for (i=0; i< size; i++)
    edo[i] = vec[i];
}

void GraphC::find_an_attractor()
{
  if (yatra) {
    cout << "[Error]: Attractor already found when GraphC::find_an_attractor was called.\n";
    exit(1);
  }
  int i, j, k, l;
  int maxsta = int(pow(2.0, size));
  int **trajectory;
  trajectory = new int*[maxsta];
  for (i = 0; i < maxsta; i++)
		trajectory[i] = new int[size];
  k = 0;
  do {
    for (i = 0; i < size; i++)
      trajectory[k][i] = edo[i];
    synchrony();
    k++;
    if (basic.eqvec(edo, size, trajectory[k-1], size)) {
			l = k-1;
			break;
		}
    l = basic.vecinmat(trajectory, k, size, edo, size);
  } while (l==(-1));
  attractor = new int*[k-l];
  for (i= 0; i<(k-l); i++)
		attractor[i] = new int[size];
  for (i = 0; i < (k-l); i++)
		for (j = 0; j < size; j++)
			attractor[i][j] = trajectory[i+l][j];
  atsize = k-l;
  palen = l; 
  for (i = 0; i < maxsta; i++)
    delete [] trajectory[i];
  delete [] trajectory;
  yatra = true;
}

void GraphC::find_attractors(int **ic, int nic)
{
  if (yatras) {
    cout << "[Error]: Attractors already found when GraphC::find_attractors was called.\n";
    exit(1);
  }
  numatrs = nic;
  attractors = new int**[numatrs];
  atsizes = new int[numatrs];
  int i,j,k;
  for (i=0; i<nic; i++) {
    set_as_state(ic[i]);
    find_an_attractor();
    atsizes[i] = atsize;
    attractors[i] = new int*[atsizes[i]];
    for (j=0; j<atsizes[i]; j++)
      attractors[i][j] = new int[size];
    for (j=0; j<atsizes[i]; j++)
      for (k=0; k< size; k++)
        attractors[i][j][k] = attractor_element(j,k);
    clear_attractor();
  }
  yatras=true;
}

int GraphC::attractor_size()
{
  if (!yatra) {
    cout << "[Error]: Attractors have not been defined when GraphC::attractor_size was called.\n";
    exit(1);
  }
  return atsize;
}

int GraphC::path_length()
{
  if (!yatra) {
    cout << "[Error]: Attractors have not been defined when GraphC::path_length was called.\n";
    exit(1);
  }
  return palen;
}

int GraphC::attractor_element(int row, int node)
{
  if (!yatra) {
    cout << "[Error]: Attractors have not been defined when GraphC::attractor_element was called.\n";
    exit(1);
  }
  if ((row >=atsize)||(row < 0)) {
    cout << "[Error]: Attractor has no " << row << "th row. GraphC::attractor_element.\n";
    exit(1);
  }
  if ((node >= size) || (node < 0)) {
    cout << "[Error]: Node " << node << " does not exist. GraphC::attractor_element.\n";
  }
  return attractor[row][node];
}

void GraphC::print_attractor(ostream& sal)
{
  if (!yatra) {
    cout << "[Error]: Attractors have not been defined when GraphC::print_attractor was called.\n";
    exit(1);
  }
  int i,j;
  for (i=0; i < atsize; i++) {
    for (j=0; j<size; j++)
      sal << attractor_element(i,j) << " ";
    sal << endl;
  }
}

int GraphC::attractor_size(int wh)
{
  if (!yatras) {
    cout << "[Error]: Attractors have not been defined when GraphC::attractor_size was called.\n";
    exit(1);
  }
  return atsizes[wh];
}

int GraphC::attractor_element(int wh, int row, int node)
{
  if (!yatras) {
    cout << "[Error]: Attractors have not been defined when GraphC::attractor_element was called.\n";
    exit(1);
  }
  if ((row >=atsizes[wh])||(row < 0)) {
    cout << "[Error]: Attractor has no " << row << "th row. GraphC::attractor_element.\n";
    exit(1);
  }
  if ((node >= size) || (node < 0)) {
    cout << "[Error]: Node " << node << " does not exist. GraphC::attractor_element.\n";
  }
  return attractors[wh][row][node];
}

void GraphC::print_attractor(int wh, ostream& sal)
{
  if (!yatras) {
    cout << "[Error]: Attractors have not been defined when GraphC::print_attractor was called.\n";
    exit(1);
  }
  int i,j;
  for (i=0; i < attractor_size(wh); i++) {
    for (j=0; j<size; j++)
      sal << attractor_element(wh,i,j) << " ";
    sal << endl;
  }
}

int GraphC::number_of_attractors()
{
  if (!yatras) {
    cout << "[Error]: Attractors have not been defined when GraphC::number_of_attractors was called.\n";
    exit(1);
  }
  return numatrs;
}

double GraphC::robustness(int *ic, int n, double wol)
{
  double rob;
  int igual=0,i,j,k,ini;
  bool sigual;
  GraphC vac(est);
  int **elatror, **elatrnu;
  elatror = new int*[attractor_size()];
  for (i=0; i < attractor_size(); i++)
    elatror[i] = new int[size];
  for (i=0; i < attractor_size(); i++)
    for (j= 0; j < size; j++)
      elatror[i][j] = attractor_element(i,j);
  for (i=0; i<n;i++) {
    make_copy(vac);
    vac.mutate(wol);
    vac.set_as_state(ic);
    vac.find_an_attractor();
    sigual = false;
    if (attractor_size() == vac.attractor_size()) {
      elatrnu = new int*[attractor_size()];
      for (j=0; j< attractor_size(); j++)
        elatrnu[j] = new int[size];
      for (j=0; j< attractor_size(); j++)
        for (k= 0; k < size; k++)
          elatrnu[j][k] = vac.attractor_element(j,k);
      ini=-1;
      for (j=0; j< attractor_size(); j++)
        if (basic.eqvec(elatror[0], size, elatrnu[j], size)) {
          ini = j;
          break;
        }
      if (ini >=0) {
        sigual = true;
        for (j=0; j< attractor_size(); j++)
          if (!basic.eqvec(elatror[j], size, elatrnu[(ini+j)%attractor_size()], size)) {
            sigual = false;
            break;
          }
      }
      for (j=0; j< attractor_size(); j++)
        delete [] elatrnu[j];
      delete [] elatrnu;
    }
    vac.clear();
    if (sigual)
      igual++;
  }
  for (i=0; i < attractor_size(); i++)
    delete [] elatror[i];
  delete [] elatror;
  rob = double(igual)/double(n);
  return rob;
}

double GraphC::robustness(int *ic, int n, double& roper, double& roperz, double wol)
{
  double rob;
  int igual=0,i,j,k,ini, coper=0, coperz=0;
  bool sigual,noz;
  GraphC vac(est);
  int **elatror, **elatrnu;
  elatror = new int*[attractor_size()];
  for (i=0; i < attractor_size(); i++)
    elatror[i] = new int[size];
  for (i=0; i < attractor_size(); i++)
    for (j= 0; j < size; j++)
      elatror[i][j] = attractor_element(i,j);
  for (i=0; i<n;i++) {
    make_copy(vac);
    vac.mutate(wol);
    vac.set_as_state(ic);
    vac.find_an_attractor();
    sigual = false;
    if (attractor_size() == vac.attractor_size()) {
      coper++;
      elatrnu = new int*[attractor_size()];
      for (j=0; j< attractor_size(); j++)
        elatrnu[j] = new int[size];
      for (j=0; j< attractor_size(); j++)
        for (k= 0; k < size; k++)
          elatrnu[j][k] = vac.attractor_element(j,k);
      noz = true;
      for (j=0; j< attractor_size(); j++)
        if (basic.vector_contains(elatrnu[j], 0, size)) {
          noz = false;
          break;
        }
      if (noz)
        coperz++;
      ini=-1;
      for (j=0; j< attractor_size(); j++)
        if (basic.eqvec(elatror[0], size, elatrnu[j], size)) {
          ini = j;
          break;
        }
      if (ini >=0) {
        sigual = true;
        for (j=0; j< attractor_size(); j++)
          if (!basic.eqvec(elatror[j], size, elatrnu[(ini+j)%attractor_size()], size)) {
            sigual = false;
            break;
          }
      }
      for (j=0; j< attractor_size(); j++)
        delete [] elatrnu[j];
      delete [] elatrnu;
    }
    vac.clear();
    if (sigual)
      igual++;
  }
  for (i=0; i < attractor_size(); i++)
    delete [] elatror[i];
  delete [] elatror;
  rob = double(igual)/double(n);
  roper = coper/double(n);
  roperz = coperz/double(n);
  return rob;
}

double GraphC::robustness(int *ic, int n, double& roper, double& roperz, double& dtotnnf, int& nufedi, double& dtotnnfz, int& nufediz, double& ronoz, double& rof1, double& rof1z, double wol)
{
  /*
   roper: fraction of single-mutations that preserve attractor size
   roperz: fraction of single-mutations that preserve attractor size and that do not produce 0s
   dtotnnf: fraction of single-mutations that preserve attractor size but that produce a different phenotype
   dtotnnfz: fraction of single-mutations that preserve attractor size but that produce a different phenotype without any 0s
   nufedi: number of different phenotypes with the same attractor size produced by single mutations
   nufediz: number of different phenotypes (lacking 0s) with the same attractor size produced by single mutations
   ronoz: fraction of single-mutations that do not produce 0s
   rof1: fraction of single-mutations that preserve fraction of 1s
   rof1z: fraction of single-mutations that preserve fraction of 1s and that do not produce 0s
  */
  double rob, fra1or, fra1n;
  int igual=0,i,j,k,l,ini, coper=0, coperz=0, totnnf, totnnfz, conoz, cof1, cof1z;
  int *nnf, *nnfz;
  nnf = new int[n];
  nnfz = new int[n];
  basic.fillv0(nnf, n);
  nufedi = 0;
  basic.fillv0(nnfz, n);
  nufediz = 0;
  totnnf = 0;
  totnnfz = 0;
  conoz = 0;
  cof1 = 0;
  cof1z = 0;
  
  int ***nf, ***nfz;
  nf = new int**[n];
  nfz = new int**[n];
  for (i=0; i < n; i++) {
    nf[i] = new int*[attractor_size()];
    nfz[i] = new int*[attractor_size()];
    for (j= 0; j < attractor_size(); j++) {
      nf[i][j] = new int[size];
      nfz[i][j] = new int[size];
    }
  }
  
  bool sigual,noz,sista;
  GraphC vac(est);
  int **elatror, **elatrnu;
  elatror = new int*[attractor_size()];
  for (i=0; i < attractor_size(); i++)
    elatror[i] = new int[size];
  l=0;
  for (i=0; i < attractor_size(); i++)
    for (j= 0; j < size; j++) {
      elatror[i][j] = attractor_element(i,j);
      if (elatror[i][j] == 1)
        l++;
    }
  fra1or = double(l)/double(size*attractor_size());
  for (i=0; i<n;i++) {
    make_copy(vac);
    vac.mutate(wol);
    vac.set_as_state(ic);
    vac.find_an_attractor();
    sigual = false;
    elatrnu = new int*[vac.attractor_size()];
    for (j=0; j< vac.attractor_size(); j++)
      elatrnu[j] = new int[size];
    l = 0;
    for (j=0; j< vac.attractor_size(); j++)
      for (k= 0; k < size; k++) {
        elatrnu[j][k] = vac.attractor_element(j,k);
        if (elatrnu[j][k] == 1)
          l++;
      }
    fra1n = double(l)/double(size*vac.attractor_size());
    if (fra1n==fra1or)
      cof1++;
    noz = true;
    for (j=0; j< vac.attractor_size(); j++)
      if (basic.vector_contains(elatrnu[j], 0, size)) {
        noz = false;
        break;
      }
    if (noz) {
      conoz++;
      if (fra1n==fra1or)
        cof1z++;
    }
    if (attractor_size() == vac.attractor_size()) {
      coper++;
      if (noz)
        coperz++;
      ini=-1;
      for (j=0; j< attractor_size(); j++)
        if (basic.eqvec(elatror[0], size, elatrnu[j], size)) {
          ini = j;
          break;
        }
      if (ini >=0) {
        sigual = true;
        for (j=0; j< attractor_size(); j++)
          if (!basic.eqvec(elatror[j], size, elatrnu[(ini+j)%attractor_size()], size)) {
            sigual = false;
            break;
          }
      }
      if (!sigual) {
        sista = false;
        totnnf++;
        for (j=0; j < nufedi; j++) {
          if (basic.eqmatrix_rot(elatrnu, attractor_size(), size, nf[j], attractor_size(), size)) {
            sista = true;
            nnf[j]++;
            break;
          }
        }
        if (!sista) {
          for (j=0; j< attractor_size(); j++)
            for (k= 0; k < size; k++)
              nf[nufedi][j][k] = elatrnu[j][k];
          nnf[nufedi] = 1;
          nufedi++;
        }
        if (noz) {
          sista = false;
          totnnfz++;
          for (j=0; j < nufediz; j++) {
            if (basic.eqmatrix_rot(elatrnu, attractor_size(), size, nfz[j], attractor_size(), size)) {
              sista = true;
              nnfz[j]++;
              break;
            }
          }
          if (!sista) {
            for (j=0; j< attractor_size(); j++)
              for (k= 0; k < size; k++)
                nfz[nufediz][j][k] = elatrnu[j][k];
            nnfz[nufediz] = 1;
            nufediz++;
          }
        }
      }
    }
    for (j=0; j< vac.attractor_size(); j++)
      delete [] elatrnu[j];
    delete [] elatrnu;
    vac.clear();
    if (sigual)
      igual++;
  }
  for (i=0; i < attractor_size(); i++)
    delete [] elatror[i];
  delete [] elatror;
  rob = double(igual)/double(n);
  roper = coper/double(n);
  roperz = coperz/double(n);
  dtotnnf = totnnf/double(n);
  dtotnnfz = totnnfz/double(n);
  ronoz = conoz/double(n);
  rof1 = cof1/double(n);
  rof1z = cof1z/double(n);
  for (i=0; i < n; i++) {
    for (j=0; j < attractor_size(); j++) {
      delete [] nf[i][j];
      delete [] nfz[i][j];
    }
    delete [] nf[i];
    delete [] nfz[i];
  }
  delete [] nf;
  delete [] nfz;
  delete [] nnf;
  delete [] nnfz;
  return rob;
}

double GraphC::robustness(int *ic, int n, double& roper, double& roperz, double& dtotnnf, int& nufedi, double& dtotnnfz, int& nufediz, double& ronoz, double& rof1, double& rof1z, double& ropl, double& roplz, double& roplper, double& roplperz, double wol)
{
  /*
   roper: fraction of single-mutations that preserve attractor size
   roperz: fraction of single-mutations that preserve attractor size and that do not produce 0s
   dtotnnf: fraction of single-mutations that preserve attractor size but that produce a different phenotype
   dtotnnfz: fraction of single-mutations that preserve attractor size but that produce a different phenotype without any 0s
   nufedi: number of different phenotypes with the same attractor size produced by single mutations
   nufediz: number of different phenotypes (lacking 0s) with the same attractor size produced by single mutations
   ronoz: fraction of single-mutations that do not produce 0s
   rof1: fraction of single-mutations that preserve fraction of 1s
   rof1z: fraction of single-mutations that preserve fraction of 1s and that do not produce 0s
   ropl: fraction of single-mutations that preserve path-length (number of steps from initial condition to first state in attractor)
   roplz: fraction of single-mutations that do not produce 0s and that preserve path-length (number of steps from initial condition to first state in attractor)
   roplper: fraction of single-mutations that preserve the sum of path-length and attractor size
   roplperz: fraction of single-mutations that preserve the sum of path-length and attractor size and that do not produce 0s
   */
  double rob, fra1or, fra1n, nupalen, nuplmat, orpalen, orplmat;
  int igual=0,i,j,k,l,ini, coper=0, coperz=0, totnnf, totnnfz, conoz, cof1, cof1z, copl=0, coplz =0, coplper=0,coplperz=0;
  int *nnf, *nnfz;
  nnf = new int[n];
  nnfz = new int[n];
  basic.fillv0(nnf, n);
  nufedi = 0;
  basic.fillv0(nnfz, n);
  nufediz = 0;
  totnnf = 0;
  totnnfz = 0;
  conoz = 0;
  cof1 = 0;
  cof1z = 0;
  
  orpalen = path_length();
  orplmat = orpalen+attractor_size();
  int ***nf, ***nfz;
  nf = new int**[n];
  nfz = new int**[n];
  for (i=0; i < n; i++) {
    nf[i] = new int*[attractor_size()];
    nfz[i] = new int*[attractor_size()];
    for (j= 0; j < attractor_size(); j++) {
      nf[i][j] = new int[size];
      nfz[i][j] = new int[size];
    }
  }
  
  bool sigual,noz,sista;
  GraphC vac(est);
  int **elatror, **elatrnu;
  elatror = new int*[attractor_size()];
  for (i=0; i < attractor_size(); i++)
    elatror[i] = new int[size];
  l=0;
  for (i=0; i < attractor_size(); i++)
    for (j= 0; j < size; j++) {
      elatror[i][j] = attractor_element(i,j);
      if (elatror[i][j] == 1)
        l++;
    }
  fra1or = double(l)/double(size*attractor_size());
  for (i=0; i<n;i++) {
    make_copy(vac);
    vac.mutate(wol);
    vac.set_as_state(ic);
    vac.find_an_attractor();
    sigual = false;
    elatrnu = new int*[vac.attractor_size()];
    for (j=0; j< vac.attractor_size(); j++)
      elatrnu[j] = new int[size];
    l = 0;
    for (j=0; j< vac.attractor_size(); j++)
      for (k= 0; k < size; k++) {
        elatrnu[j][k] = vac.attractor_element(j,k);
        if (elatrnu[j][k] == 1)
          l++;
      }
    fra1n = double(l)/double(size*vac.attractor_size());
    if (fra1n==fra1or)
      cof1++;
    noz = true;
    for (j=0; j< vac.attractor_size(); j++)
      if (basic.vector_contains(elatrnu[j], 0, size)) {
        noz = false;
        break;
      }
    if (noz) {
      conoz++;
      if (fra1n==fra1or)
        cof1z++;
    }
    nupalen = vac.path_length();
    nuplmat = nupalen+vac.attractor_size();
    if (nupalen==orpalen) {
      copl++;
      if (noz)
        coplz++;
    }
    if (nuplmat== orplmat) {
      coplper++;
      if (noz)
        coplperz++;
    }
    if (attractor_size() == vac.attractor_size()) {
      coper++;
      if (noz)
        coperz++;
      ini=-1;
      for (j=0; j< attractor_size(); j++)
        if (basic.eqvec(elatror[0], size, elatrnu[j], size)) {
          ini = j;
          break;
        }
      if (ini >=0) {
        sigual = true;
        for (j=0; j< attractor_size(); j++)
          if (!basic.eqvec(elatror[j], size, elatrnu[(ini+j)%attractor_size()], size)) {
            sigual = false;
            break;
          }
      }
      if (!sigual) {
        sista = false;
        totnnf++;
        for (j=0; j < nufedi; j++) {
          if (basic.eqmatrix_rot(elatrnu, attractor_size(), size, nf[j], attractor_size(), size)) {
            sista = true;
            nnf[j]++;
            break;
          }
        }
        if (!sista) {
          for (j=0; j< attractor_size(); j++)
            for (k= 0; k < size; k++)
              nf[nufedi][j][k] = elatrnu[j][k];
          nnf[nufedi] = 1;
          nufedi++;
        }
        if (noz) {
          sista = false;
          totnnfz++;
          for (j=0; j < nufediz; j++) {
            if (basic.eqmatrix_rot(elatrnu, attractor_size(), size, nfz[j], attractor_size(), size)) {
              sista = true;
              nnfz[j]++;
              break;
            }
          }
          if (!sista) {
            for (j=0; j< attractor_size(); j++)
              for (k= 0; k < size; k++)
                nfz[nufediz][j][k] = elatrnu[j][k];
            nnfz[nufediz] = 1;
            nufediz++;
          }
        }
      }
    }
    for (j=0; j< vac.attractor_size(); j++)
      delete [] elatrnu[j];
    delete [] elatrnu;
    vac.clear();
    if (sigual)
      igual++;
  }
  for (i=0; i < attractor_size(); i++)
    delete [] elatror[i];
  delete [] elatror;
  rob = double(igual)/double(n);
  roper = coper/double(n);
  roperz = coperz/double(n);
  dtotnnf = totnnf/double(n);
  dtotnnfz = totnnfz/double(n);
  ronoz = conoz/double(n);
  rof1 = cof1/double(n);
  rof1z = cof1z/double(n);
  ropl = double(copl)/double(n);
  roplz = double(coplz)/double(n);
  roplper = double(coplper)/double(n);
  roplperz = double(coplperz)/double(n);
  for (i=0; i < n; i++) {
    for (j=0; j < attractor_size(); j++) {
      delete [] nf[i][j];
      delete [] nfz[i][j];
    }
    delete [] nf[i];
    delete [] nfz[i];
  }
  delete [] nf;
  delete [] nfz;
  delete [] nnf;
  delete [] nnfz;
  return rob;
}

void GraphC::robustness_mult(int **ic, double *robs, int nic, int n, double wol)
{
  int i,j,k,l,ini;
  int *atrorsizes;
  atrorsizes = new int[nic];
  for (i=0; i< nic; i++)
    atrorsizes[i] = attractor_size(i);
  int ***losatror;
  int *vanigu;
  vanigu = new int[nic];
  basic.fillv0(vanigu, nic);
  losatror = new int**[nic];
  for (i=0; i< nic; i++) {
    losatror[i] = new int*[atrorsizes[i]];
    for (j=0; j < atrorsizes[i]; j++)
      losatror[i][j] = new int[size];
  }
  for (i=0; i< nic; i++)
    for (j= 0; j < atrorsizes[i]; j++)
      for (k= 0; k < size; k++)
        losatror[i][j][k] = attractor_element(i,j,k);
  GraphC vac(est);
  bool sigual;
  int ***nueatrs;
  nueatrs = new int**[nic];
  int *nuesizes;
  nuesizes = new int[nic];
  
  for (i=0; i<n; i++) {
    make_copy(vac);
    vac.mutate(wol);
    vac.find_attractors(ic, nic);
    for (j=0; j < nic; j++) {
      nuesizes[j] = vac.attractor_size(j);
      nueatrs[j] = new int*[nuesizes[j]];
      for (k=0; k < nuesizes[j]; k++) {
        nueatrs[j][k] = new int[size];
        for (l=0; l < size; l++)
          nueatrs[j][k][l] = vac.attractor_element(j,k,l);
      }
    }
    for (j=0; j < nic; j++) {
      sigual = false;
      ini = -1;
      if (attractor_size(j) == vac.attractor_size(j)) {
        for (k=0; k < nuesizes[j]; k++)
          if (basic.eqvec(losatror[j][0], size, nueatrs[j][k], size)) {
            ini = k;
            break;
          }
        if (ini >=0) {
          sigual = true;
          for (k=0; k < atrorsizes[j]; k++)
            if (!basic.eqvec(losatror[j][k], size, nueatrs[j][(ini+k)%attractor_size()], size))  {
              sigual = false;
              break;
            }
        }
      }
      if (sigual)
        vanigu[j]++;
    }
    for (j=0; j<nic; j++) {
      for (k=0; k < nuesizes[j]; k++)
        delete [] nueatrs[j][k];
      delete [] nueatrs[j];
    }
    vac.clear();
  }
  for (i=0; i < nic; i++) {
    for (j=0; j < atrorsizes[i]; j++)
      delete [] losatror[i][j];
    delete [] losatror[i];
  }
  for (i=0; i < nic; i++)
    robs[i] = double(vanigu[i])/double(n);
  
  delete [] losatror;
  delete [] atrorsizes;
  delete [] vanigu;
  delete [] nuesizes;
  delete [] nueatrs;
}
//desdeaqui


int GraphC::penetrance_ap_n1(int *cior, int *altphen) {
  int cuen = 0;
  int i;
  int *vt;
  vt = new int[size];
  for (i=0; i < size; i++)
    vt[i] = cior[i];
  for (i=0; i < size; i++) {
    vt[i]*=(-1);
    set_as_state(vt);
    find_an_attractor();
    if (attractor_size() == 1)
      if (basic.eqvec(attractor[0], size, altphen, size))
        cuen++;
    clear_attractor();
    vt[i]*=(-1);
  }
  return cuen;
}

void GraphC::penetrance_n1_2p(int *cior, int *pA, int *pB, int& countA, int& countB) {
  int i;
  countA = 0;
  countB = 0;
  int *vt;
  vt = new int[size];
  for (i=0; i < size; i++)
    vt[i] = cior[i];
  for (i=0; i < size; i++) {
    vt[i]*=(-1);
    set_as_state(vt);
    find_an_attractor();
    if (attractor_size() == 1) {
      if (basic.eqvec(attractor[0], size, pA, size))
        countA++;
      else if (basic.eqvec(attractor[0], size, pB, size))
        countB++;
    }
    clear_attractor();
    vt[i]*=(-1);
  }
  return;
}

int GraphC::penetrance_ap_n2_ex(int *cior, int *altphen) {
  int cuen = 0;
  int i,j;
  int *vt;
  vt = new int[size];
  for (i=0; i < size; i++)
    vt[i] = cior[i];
  for (i=0; i < (size-1); i++) {
    vt[i]*=(-1);
    for (j = i+1; j < size; j++) {
      vt[j] *=(-1);
      set_as_state(vt);
      find_an_attractor();
      if (attractor_size() == 1)
        if (basic.eqvec(attractor[0], size, altphen, size))
          cuen++;
      clear_attractor();
      vt[j] *=(-1);
    }
    vt[i]*=(-1);
  }
  return cuen;
}

void GraphC::penetrance_n2_ex_2p(int *cior, int *pA, int *pB, int& countA, int& countB) {
  countA = 0;
  countB = 0;
  int i,j;
  int *vt;
  vt = new int[size];
  for (i=0; i < size; i++)
    vt[i] = cior[i];
  for (i=0; i < (size-1); i++) {
    vt[i]*=(-1);
    for (j = i+1; j < size; j++) {
      vt[j] *=(-1);
      set_as_state(vt);
      find_an_attractor();
      if (attractor_size() == 1) {
        if (basic.eqvec(attractor[0], size, pA, size))
          countA++;
        else if (basic.eqvec(attractor[0], size, pB, size))
          countB++;
      }
      clear_attractor();
      vt[j] *=(-1);
    }
    vt[i]*=(-1);
  }
  return;
}

int GraphC::mutational_access(int *ci, int *altphen) {
  int cuen = 0;
  int i,j,k, nume, numut;
  double vval, nval;
  nume = number_of_edges();
  numut = (2*size*size)-nume;
  for (k = 0; k < numut; k++) {
    i = est.randint(0, size);
    j = est.randint(0, size);
    vval = weight(i,j);
    if (vval!=0) {
      if (est.toss()) {
        force_interaction(i,j,0);
        set_as_state(ci);
        find_an_attractor();
        if (attractor_size() == 1)
          if (basic.eqvec(attractor[0], size, altphen, size))
            cuen++;
      } else {
        nval = est.randgauss();
        if ((nval*vval)<0)
          nval*=(-1);
        force_interaction(i,j,nval);
        set_as_state(ci);
        find_an_attractor();
        if (attractor_size() == 1)
          if (basic.eqvec(attractor[0], size, altphen, size))
            cuen++;
      }
      clear_attractor();
    } else {
      nval = est.randgauss();
      force_interaction(i,j,nval);
      set_as_state(ci);
      find_an_attractor();
      if (attractor_size() == 1)
        if (basic.eqvec(attractor[0], size, altphen, size))
          cuen++;
      clear_attractor();
    }
    force_interaction(i,j,vval);
  }
  return cuen;
}

void GraphC::mutational_access_2p(int *ci, int *pA, int *pB, int& countA, int& countB) {
  countA = 0;
  countB = 0;
  int i,j,k, nume, numut;
  double vval, nval;
  nume = number_of_edges();
  numut = (2*size*size)-nume;
  for (k = 0; k < numut; k++) { //aqui
    i = est.randint(0, size);
    j = est.randint(0, size);
    vval = weight(i,j);
    if (vval!=0) {
      if (est.toss()) {
        force_interaction(i,j,0);
        set_as_state(ci);
        find_an_attractor();
        if (attractor_size() == 1) {
          if (basic.eqvec(attractor[0], size, pA, size))
            countA++;
          else if (basic.eqvec(attractor[0], size, pB, size))
            countB++;
        }
      } else {
        nval = est.randgauss();
        if ((nval*vval)<0)
          nval*=(-1);
        force_interaction(i,j,nval);
        set_as_state(ci);
        find_an_attractor();
        if (attractor_size() == 1) {
          if (basic.eqvec(attractor[0], size, pA, size))
            countA++;
          else if (basic.eqvec(attractor[0], size, pB, size))
            countB++;
        }
      }
      clear_attractor();
    } else {
      nval = est.randgauss();
      force_interaction(i,j,nval);
      set_as_state(ci);
      find_an_attractor();
      if (attractor_size() == 1) {
        if (basic.eqvec(attractor[0], size, pA, size))
          countA++;
        else if (basic.eqvec(attractor[0], size, pB, size))
          countB++;
      }
      clear_attractor();
    }
    force_interaction(i,j,vval);
  }
  return;
}

void GraphC::mutational_access_2p_eqmu(int *ci, int *pA, int *pB, int& countA, int& countB) {
  countA = 0;
  countB = 0;
  int i,j,k, nume, numut;
  double vval, nval;
  nume = number_of_edges();
  numut = (2*size*size)-nume;
  for (k = 0; k < numut; k++) { //aqui
    i = est.randint(0, size);
    j = est.randint(0, size);
    vval = weight(i,j);
    if (vval!=0) {
      if (est.toss()) {
        force_interaction(i,j,0);
        set_as_state(ci);
        find_an_attractor();
        if (attractor_size() == 1) {
          if (basic.eqvec(attractor[0], size, pA, size)) {
            countA++;
            countA++;
          }
          else if (basic.eqvec(attractor[0], size, pB, size)) {
            countB++;
            countB++;
          }
        }
      } else {
        nval = est.randgauss();
        if ((nval*vval)<0)
          nval*=(-1);
        force_interaction(i,j,nval);
        set_as_state(ci);
        find_an_attractor();
        if (attractor_size() == 1) {
          if (basic.eqvec(attractor[0], size, pA, size))
            countA++;
          else if (basic.eqvec(attractor[0], size, pB, size))
            countB++;
        }
      }
      clear_attractor();
    } else {
      nval = est.randgauss();
      force_interaction(i,j,nval);
      set_as_state(ci);
      find_an_attractor();
      if (attractor_size() == 1) {
        if (basic.eqvec(attractor[0], size, pA, size))
          countA++;
        else if (basic.eqvec(attractor[0], size, pB, size))
          countB++;
      }
      clear_attractor();
    }
    force_interaction(i,j,vval);
  }
  return;
}

void GraphC::build_offspring(GraphC &mother, GraphC &father, bool *frommom) {
  if (mother.number_of_nodes() != father.number_of_nodes()) {
    cout << "[Error]: Parents with different genome size. GraphC::build_offspring.\n";
    exit(1);
  }
  if ((!mother.is_directed()) || (!father.is_directed())) {
    cout << "[Error]: One parent is not a directed network. GraphC::build_offspring.\n";
    exit(1);
  }
  int i, j, n;
  n = father.number_of_nodes();
  make_nw(n, true);
  for (i=0; i < size; i++) {
    if (frommom[i]) {
      for (j= 0; j < size; j++)
        force_interaction(j,i,mother.weight(j,i));
    } else {
      for (j= 0; j < size; j++)
        force_interaction(j,i,father.weight(j,i));
    }
  }
  return;
}

//hastaaqui
void GraphC::mate(GraphC mother, GraphC father)
{
  if (mother.number_of_nodes() != father.number_of_nodes()) {
		cout << "[Error]: Parents with different genome size. GraphC::mate.\n";
		exit(1);
	}
	if (mother.is_directed() != father.is_directed()) {
		cout << "[Error]: One parent is a directed network while the other is an undirected network. GraphC::mate.\n";
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
          force_interaction(j,i,father.weight(j,i));
			}
			else {
				for (j = 0; j < size; j++)
          force_interaction(j,i,mother.weight(j,i));
			}
		}
	}
	else {
		for (i=0; i<size; i++)
			for (j=0; j<=i; j++) {
				paoma = est.toss();
				if (paoma)
          force_interaction(j,i,father.weight(j,i));
				else
          force_interaction(j,i,mother.weight(j,i));
			}
	}
}

void GraphC::mutate_old(double con)
{
  if (!directed) {
    cout << "[Error]: GraphC::mutate is not defined for undirected networks.\n";
    exit(1);
  }
  int i,j;
  double newval;
  double expon, conreal;
  if (est.toss()) {
    if (number_of_edges()!=0) {
      do {
        i = est.randint(0,size);
        j = est.randint(0,size);
      } while (weight(j,i)==0);
      newval = fabs(est.randgauss());
      if (weight(j,i) < 0)
        newval = newval*(-1);
      change_interaction(j,i,newval);
    }
  }
  else {
    conreal = double(number_of_edges())/double(size*size);
    expon = log(0.5)/log(con);
    if (est.randreal() < pow(conreal, expon)) {
      do {
        i = est.randint(0,size);
        j = est.randint(0,size);
      } while (weight(j,i)==0);
      change_interaction(j,i,0);
    }
    else {
      do {
        i = est.randint(0,size);
        j = est.randint(0,size);
      } while (weight(j,i)!=0);
      newval = est.randgauss();
      change_interaction(j,i,newval);
    }
  }
}

void GraphC::mutate(double con) {
  int i = est.randint(0, size);
  mutate_gene(i, con);
}

void GraphC::mutate_gene(int gene, double con) {
  if (!directed) {
    cout << "[Error]: GraphC::mutate_gene is not defined for undirected networks.\n";
    exit(1);
  }
  bool sube = false;
  int j;
  j = est.randint(0, size);
  double newval, oldval = weight(j,gene);
  if (est.randreal() < con)
    sube = true;
  if (sube) {
    newval = est.randgauss();
    if (oldval < 0)
      newval = fabs(newval)*(-1);
    else {
      if (oldval > 0)
        newval = fabs(newval);
    }
  }
  else
    newval = 0;
  force_interaction(j, gene, newval);
}

void GraphC::mutate_gene_mp5(int gene) {
  if (!directed) {
    cout << "[Error]: GraphC::mutate_gene is not defined for undirected networks.\n";
    exit(1);
  }
  int j;
  j = est.randint(0, size);
  double newval, oldval = weight(j,gene);
  if (oldval == 0)
    newval = est.randgauss();
  else
    newval = 0;
  force_interaction(j, gene, newval);
}

void GraphC::mutate_gene_mf(int gene) {
  if (!directed) {
    cout << "[Error]: GraphC::mutate_gene is not defined for undirected networks.\n";
    exit(1);
  }
  int j, numz;
  double newval = 0;
  if (est.toss())
    newval = est.randgauss();
  numz = basic.count_in_vector(nw[gene], size, 0);
  if (((newval==0) && (numz == size)) || ((newval!=0) && (numz == 0)))
    return;
  if (newval!=0) {
    do {
      j = est.randint(0, size);
    } while (weight(j, gene)!=0);
    force_interaction(j, gene, newval);
  }
  else {
    do {
      j = est.randint(0, size);
    } while (weight(j, gene)==0);
    force_interaction(j, gene, newval);
  }
}

void GraphC::consider_mutation(double muratepg, double con) {
  int i;
  for (i = 0; i < size; i++)
    if (est.randreal() < muratepg)
      mutate_gene(i, con);
}

void GraphC::consider_mutation_mp5(double muratepg) {
  int i;
  for (i = 0; i < size; i++)
    if (est.randreal() < muratepg)
      mutate_gene_mp5(i);
}

void GraphC::consider_mutation_mf(double muratepg) {
  int i;
  for (i = 0; i < size; i++)
    if (est.randreal() < muratepg)
      mutate_gene_mf(i);
}

void GraphC::duplicate(int gene)
{
  double **matemp;
  int i,j;
  int n = size+1;
  matemp = new double*[n];
  for (i = 0; i < n; i++)
		matemp[i] = new double[n];
  basic.fillmat0(matemp, n, n);
  for (i = 0; i < size; i++)
    for (j = 0; j < size; j++)
      matemp[i][j] = nw[i][j];
	bool dirnw = is_directed(); //
  clear();
  make_nw(n, dirnw);
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

void GraphC::duplicate()
{
  int i = est.randint(0, size);
  duplicate(i);
}

//motifs
int GraphC::count_selfinteracting_nodes()
{
	int i,res =0;
	for (i=0; i<size; i++)
		if (weight(i,i) != 0) //
			res++;
	return res;
}

int GraphC::count_non_selfinteracting_nodes()
{
	int i,res =0;
	for (i=0; i<size; i++)
		if (weight(i,i)==0) //
			res++;
	return res;
}

//modularity

void GraphC::build_moma_u(double **momaundi) {
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

void GraphC::build_moma_d(double **momadi) {
  //Momadi corresponds to the final modularity matrix (B + B^T according to Leicht and Newman 2008)
  if (!directed) {
    cout << "[Error]: Attempt to build directed modularity matrix from undirected adjacency matrix in GraphC::build_moma_d\n";
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

void GraphC::build_moma_out(double **momaout) {
  if (!directed) {
    cout << "[Error]: Attempt to build directed modularity matrix from undirected adjacency matrix in GraphC::build_moma_out\n";
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

void GraphC::build_moma_in(double **momain) {
  if (!directed) {
    cout << "[Error]: Attempt to build directed modularity matrix from undirected adjacency matrix in GraphC::build_moma_in\n";
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

int GraphC::components_gt1() {
  if (!yacomp)
    get_components();
  int i, j = 0;
  for (i = 0; i < number_of_components(); i++)
    if (basic.count_in_vector(components, size, i) == 1)
      j++;
  return number_of_components() - j;
}

int GraphC::sccs_gt1() {
  if (!cfc)
    get_scc();
  int i, j = 0;
  for (i = 0; i < number_of_sccs(); i++)
    if (basic.count_in_vector(scc, size, i) == 1)
      j++;
  return number_of_sccs() - j;
}

void GraphC::check_props_of_rdset(int samplesiz, int minsw, int maxsw, double **meansdminmax) {
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
  GraphC atiza(est);
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

void GraphC::check_props_of_rdset(int samplesiz, int minsw, int maxsw, double **meansdminmax, set<set<int> > &predpar) {
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
  GraphC atiza(est);
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

void GraphC::check_Qd_in_rdset(int samplesiz, int minsw, int maxsw, double *meansdminmax) {
  //  hilera 0: moddir
  double *todo;
  int e = number_of_edges();
  todo = new double[samplesiz];
  int i, j;
  GraphC atiza(est);
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

void GraphC::check_Qd_in_rdset(int samplesiz, int minsw, int maxsw, double **meansdminmax, set<set<int> > &predpar) {
  //  hilera 0: moddir
  //hilera1: moddirpred
  double **todo;
  int e = number_of_edges();
  basic.create_array(todo, 2, samplesiz);
  int i, j;
  GraphC atiza(est);
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

void GraphC::build_moma_u_aux(int **adjundi, int *degundi, int edgundim, double **momaundi) {
  int i,j;
  for (i=0; i < size; i++)
    for (j = 0; j < size; j++)
      momaundi[i][j] = adjundi[i][j] - ((degundi[i]*degundi[j])/(2.0*edgundim));
}

void GraphC::mod_partitions_from_vector(double *vecs, int tam, int *diccion, set<int> &parta, set<int> &partb) {
  
  parta.clear();
  partb.clear();
  int i;
  for (i=0; i < tam; i++) {
    if (vecs[i] == 1)
      parta.insert(diccion[i]);
    else if (vecs[i] == -1)
      partb.insert(diccion[i]);
    else {
      cout << "[Error]: Vector contains numbers other than 1 or -1 in GraphC::mod_partitions_from_vector.\n";
      exit(1);
    }
  }
}

double GraphC::spectral_method(double **moduloc, int tam, double *vecs) {
  if (!basic.is_symmetric(moduloc, tam, tam)) {
    cout << "[Error]: Non-symmetric matrix in GrahB::spectral_method.\n";
    exit(1);
  }
  double eval = basic.get_leading_evector_symmat(moduloc, tam, vecs);
  return eval;
}

double GraphC::mod_after_spectral(double **matloc, int tam, double *vecs) {
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

void GraphC::shake_kl(double mienQ, double **matloc, int tam, double *vecs) {
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

double GraphC::eval_mod_subm(double **mamod, double *vecs, int tam) {
  double *vt;
  vt = new double[tam];
  basic.matxvec(vecs, tam, tam , tam, mamod, vt);
  double Q = basic.dotproduct(vecs, vt, tam);
  delete [] vt;
  Q /= (4.0*number_of_edges());
  return Q;
}

void GraphC::partition_to_vector(const set<set<int> > &equipos, int *vop) {
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

double GraphC::eval_mod(double **mamod, int *vop) {
  int i, j;
  double Q = 0;
  for (i = 0; i < size; i++)
    for (j = 0; j < size; j++)
      if (vop[i] == vop[j])
        Q += mamod[i][j];
  Q /= (2.0*number_of_edges());
  return Q;
}

double GraphC::eval_mod(double **mamod, set<set<int> > &equipos) {
  int *arr;
  arr = new int[size];
  partition_to_vector(equipos, arr);
  double Q = eval_mod(mamod, arr);
  delete [] arr;
  return Q;
}

double GraphC::maxmod_u_per_part_and_dd(int *vop) {
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

double GraphC::maxmod_u_per_part_and_dd(set<set<int> > &equipos) {
  int *arr;
  arr = new int[size];
  partition_to_vector(equipos, arr);
  double Qmax = maxmod_u_per_part_and_dd(arr);
  delete [] arr;
  return Qmax;
}
//
double GraphC::maxmod_d_per_part_and_dd(int *vop) {
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

double GraphC::maxmod_d_per_part_and_dd(set<set<int> > &equipos) {
  int *arr;
  arr = new int[size];
  partition_to_vector(equipos, arr);
  double Qmax = maxmod_d_per_part_and_dd(arr);
  delete [] arr;
  return Qmax;
}

void GraphC::adjust_modmat(double **mator, set<int> &parta, int tam, int *diccion, double **matloc) {
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
  //prueba
  //  double suma = 0;
  //  for (i = 0; i < tam; i++)
  //    for (j = 0; j < tam; j++)
  //      suma += matloc[i][j];
  //  if (fabs(suma) > 0.0001) {
  //    cout << "No suma cero al hacer particion " << suma;
  //    basic.printset(parta, cout);
  //    cout << endl;
  //    exit(1);
  //  }
}

void GraphC::iterative_newman06(double **mator, set<set<int> > &res) {
  int i;
  //  //prue
  //  if (res.size() != 0) {
  //    cout << "ahjijo en iterative...\n";
  //    exit(1);
  //  }
  int cuantu;
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
        cout << "[Error]: Vector contains numbers other than 1 or -1 in GraphC::mod_partitions_from_vector.\n";
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

void GraphC::iterative_newman06(double **mator, set<int> &parta, set<set<int> > &res) {
  int i,*diccion, tam = parta.size();
  bool otra = false;
  int cuantu;
  //  //prue
  //  if (res.size() != 0) {
  //    cout << "ahjijo en iterative...\n";
  //    exit(1);
  //  }
  
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

void GraphC::split_in_2mods(double **mator, set<set<int> > &res) {
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
        cout << "[Error]: Vector contains numbers other than 1 or -1 in GraphC::mod_partitions_from_vector.\n";
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
void GraphC::set_default_exclusive_vars()
{
  infl = false;
  dist = false;
  yacomp = false;
  cfc = false;
  bcya = false;
  floops = false;
  yae = false;
  yadeg = false;
  yatra = false;
	yadegdist=false;
  yatras= false;
  yamadya = false;
}

int GraphC::calc_indegree(int no)
{
  int id = 0;
  int i;
  for (i = 0; i < size; i++)
    if (nw[no][i] != 0)
      id++;
  return id;
}

int GraphC::calc_outdegree(int no)
{
  int od = 0;
  int i;
  for (i = 0; i < size; i++)
    if (nw[i][no] != 0)
      od++;
  return od;
}

int GraphC::calc_degree(int no)
{
	int d;
	if (directed)
		d = indegree[no] + outdegree[no];
	else {
		if (indegree[no] != outdegree[no]) {
			cout << "[Error]: Different indegree and outdegree in an undirected network. GraphC::calc_degree.\n";
			exit(1);
		}
		d = indegree[no];
		if (weight(no, no) != 0)
			d++;
	}
	return d;
}

void GraphC::prepare_for_degrees()
{
	degree = new int[size];
	outdegree = new int[size];
	indegree = new int[size];
}

void GraphC::prepare_for_degdist()
{
	degdist = new int[size+1];
	basic.fillv0(degdist, size+1);
	
	odegdist = new int[size+1];
	basic.fillv0(odegdist, size+1);
	
	idegdist = new int[size+1];
	basic.fillv0(idegdist, size+1);
}

void GraphC::priv_get_degdist()
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

void GraphC::aswitch()
{
  int i, j, k, l, cont=0;
  double v1, v2,v3,v4;
  int nesq = number_of_edges()*number_of_edges();
  bool swdone = false;
  do {
    do {
      i = est.randint(0, size);
      j = est.randint(0, size);
    }while (nw[i][j] == 0);
    do {
      k = est.randint(0, size);
      l = est.randint(0, size);
    }while (nw[k][l] == 0);
    if (directed) { //
      if (!((nw[i][l] != 0) || (nw[k][j] != 0) || (i==k) || (j==l) )) {
        v1 = nw[i][j];
        v2 = nw[k][l];
        nw[i][j] = 0;
        nw[k][l] = 0;
        nw[i][l] = v1;
        nw[k][j] = v2;
        swdone = true; //
      }
    }
    else {
      if (!((nw[i][l] != 0) || (nw[k][j] != 0) || (i==k) || (j==l) || ((i==j) && (k==l)) || ((i==l) && (k==j)))) {
        v1 = nw[i][j];
        v2 = nw[k][l];
        v3 = nw[j][i];
        v4 = nw[l][k];
        nw[i][j] = 0;
        nw[k][l] = 0;
        nw[j][i] = 0;
        nw[l][k] = 0;
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
    cout << "[Error]: It is hard to find switchable edges. GraphC::aswitch.\n";
    exit(1);
  }
}

