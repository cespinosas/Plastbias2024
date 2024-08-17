#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <string>
#include "alea.h"
#include "basics.h"
#include "graphb.h"
#include "graphi.h"
#include "graphc.h"
#include "fitnessi.h"
#include "evolvei.h"


using std::cout;
using std::cin;
using std::getline;
using std::string;
using std::endl;
using std::ofstream;
using std::ifstream;
using namespace std;

int main(int argc, char **argv) {
  if ((argc != 10) && (argc != 16)) {
    cout << "[Error]; Wrong number of arguments in program. It should be ./asu seed nodes edges popsize murate generations pertrate selcoefpg gfromintsctn_one ...\n";
    exit(1);
  }
  ofstream fs;
  int seed = atoi(argv[1]);
  int nodes = atoi(argv[2]);
  int edges = atoi(argv[3]);
  int popsize = atoi(argv[4]);
  double murate = atof(argv[5]);
  int generations = atoi(argv[6]);
  double wol = double(edges)/double(nodes*nodes);
  double pertrate = atof(argv[7]);
  double selcoefpg = atof(argv[8]);
  //  int numofcells = 1;//atoi(argv[9]); UNA SOLA CELULA
  int numint = atoi(argv[9]); //numint
  bool genfromint = false;
  if (numint == 1)
    genfromint = true;
  
  ifstream fe;
  int dcia, dcib, dab, dcic, dbc, dac;
  if (argc == 10) {
    dcia = 3;
    dcib = 3;
    dcic = 3;
    dab = 6;
    dbc = 6;
    dac = 6;
  } else {
    dcia = atoi(argv[10]);
    dcib = atoi(argv[11]);
    dcic = atoi(argv[12]);
    dab = atoi(argv[13]);
    dbc = atoi(argv[14]);
    dac = atoi(argv[15]);
  }
  
  Alea jacta(seed);
  Basics bas(jacta);
  EvolveI chw(jacta);
  FitnessI fifu(jacta);
  string orden;
  
  int i, j,k,l;
  double *defaultA, *defaultB, *defaultC, *defaultO, *wmax, *meanedges;
  
  defaultA = new double[generations];
  defaultB = new double[generations];
  defaultC = new double[generations];
  defaultO = new double[generations];
  wmax = new double[generations];
  meanedges = new double[generations];
  
  bas.fillvm1(defaultA, generations);
  bas.fillvm1(defaultB, generations);
  bas.fillvm1(defaultC, generations);
  bas.fillvm1(defaultO, generations);
  bas.fillvm1(wmax, generations);
  bas.fillvm1(meanedges, generations);
  int tpC = -1, tpB = -1, t50C = -1, t50B = -1, tdC = -1, tdB = -1;
  
  //nuevo
  int t5C = -1, t80C = -1, t5B = -1, t80B = -1;
  int numredaccAfromB, numredaccAfromC, numperaccAfromB, numperaccAfromC, contemp;
  
  bool *cualesdefB, *cualesdefC;
  cualesdefB = new bool[popsize];
  cualesdefC = new bool[popsize];
  double meanaccpredAfromB[3] = {-1, -1, -1}; //mod 5, 50, 80
  double meanaccpredAfromC[3] = {-1, -1, -1}; //mod 5, 50, 80
  double meanaccpperAfromB[3] = {-1, -1, -1}; //mod 5, 50, 80
  double meanaccpperAfromC[3] = {-1, -1, -1}; //mod 5, 50, 80
  //hasta aqui
  
  int *ci, *phenA, *phenB, *phenC;
  ci = new int[nodes];
  phenA = new int[nodes];
  phenB = new int[nodes];
  phenC = new int[nodes];
  for (i=0; i < nodes; i++) {
    if ((i%2)==0)
      ci[i] = 1;
    else
      ci[i] = -1;
    if (i < (nodes-dcia))
      phenA[i] = ci[i];
    else
      phenA[i] = ci[i]*(-1);
    if (i < dcic)
      phenC[i] = ci[i]*(-1);
    else
      phenC[i]  = ci[i];
    phenB[i] = ci[i];
  }
  if ((dcia+dcib)==dab) {
    for (i = (nodes-dcia)-1; i >= (nodes-dab); i--)
    phenB[i]*=(-1);
  } else {
    if ((dcia+dcib)>dab) {
      for (i = (nodes - ((dcia+dcib+dab)/2)); i < (nodes + ((dcib-dab-dcia)/2)); i++)
      phenB[i]*=(-1);
    } else {
      cout << "[Error]: Impossible phenotypes.\n";
      exit(1);
    }
  }
  if ((bas.difsinvecs(ci, nodes, phenA, nodes) != dcia) || (bas.difsinvecs(ci, nodes, phenB, nodes) != dcib) ||
      (bas.difsinvecs(phenB, nodes, phenA, nodes) != dab) || (bas.difsinvecs(ci, nodes, phenC, nodes) != dcic) || (bas.difsinvecs(phenC, nodes, phenB, nodes) != dbc) || (bas.difsinvecs(phenC, nodes, phenA, nodes) != dac)) {
    cout << "[Error]: Distances do not match.\n";
    exit(1);
  }
  
  int *ci2;
  ci2 = new int[nodes];
  for (i=0; i < nodes; i++) {
    ci2[i] = ci[i];
  }
  for (i=0; i < nodes; i++) {
    if ((ci[i] == phenA[i]) && (phenA[i] == phenB[i]) && (phenC[i] == phenB[i])) {
      ci2[i] *= (-1);
      break;
    }
  }
  
  
  string cuerda;
  
  if (genfromint)
    cuerda = "Results/pn1/";
  else
    cuerda = "Results/sinplas/";
  
  chw.set_params(murate, wol, popsize, nodes);
  bas.open_ofstream(fs, cuerda+"pars/"+bas.inttostring(seed)+"pars.txt");
  chw.print_params(fs);
  fs.close();
  
  GraphI adam(jacta);
  if (genfromint)
    cuerda = "networks/plas/"+bas.inttostring(seed)+".txt";
  else
    cuerda = "networks/sinp/"+bas.inttostring(seed)+".txt";
  
  bas.open_ifstream(fe, cuerda);
  adam.get_dir_nw_from_file(nodes, fe, edges);
  fe.close();
  
  int **matatr;
  
  chw.start_pop(adam);//nodes, edges, ics[0], goals[0][0]);
  double w, wb, wc, prope, disb, disc;
  int defa, defb, defc, sume;
  
  int *oci;
  oci = new int[nodes];
  for (i=0; i< nodes; i++)
  oci[i] = ci[i];
  
  //term//
  for (i=0; i < 20; i++) {
    for (j=0; j < popsize; j++) {
      w = 0;
      chw.population[j].set_as_state(ci);
      chw.population[j].find_an_attractor();
      if (fifu.strict(chw.population[j], phenA) > 0)
        w=1;
      chw.population[j].clear_attractor();
      if (genfromint) {
        if (w==1) {
          chw.population[j].set_as_state(ci2);
          chw.population[j].find_an_attractor();
          if (!(fifu.strict(chw.population[j], phenB) > 0))
            w=0;
          chw.population[j].clear_attractor();
        }
      }
      chw.assign_w(j, w);
    }
    chw.calc_meanw();
    chw.one_generation();
  }
  
  
  
  
  for (i=0; i<generations; i++) {
    //    if ((i%100)==0)
    //      cout << i << endl;
    defa = 0;
    defb = 0;
    defc = 0;
    sume = 0;
    bas.fillv0(cualesdefB, popsize); //nuevo
    bas.fillv0(cualesdefC, popsize); //nuevo
    for (j=0; j < popsize; j++) {
      sume += chw.population[j].number_of_edges();
      //para contar
      for (k=0; k < nodes; k++) {
        oci[k] = ci[k];
      }
      chw.population[j].set_as_state(oci);
      chw.population[j].find_an_attractor();
      if (fifu.strict(chw.population[j], phenA) > 0) {
        defa++;
      }
      else {
        if (fifu.strict(chw.population[j], phenB) > 0) {
          defb++;
          cualesdefB[j] = true; //nuevo
        }
        else {
          if (fifu.strict(chw.population[j], phenC) > 0) {
            defc++;
            cualesdefC[j] = true; //nuevo
          }
        }
      }
      chw.population[j].clear_attractor();
      //paraevolu:
      for (k=0; k < nodes; k++) {
        prope = jacta.randreal();
        if (prope < pertrate)
          oci[k] = ci2[k]*(-1); //mod
        else
          oci[k] = ci2[k]; //mod
      }
      chw.population[j].set_as_state(oci);
      chw.population[j].find_an_attractor();
      bas.create_array(matatr, chw.population[j].attractor_size(), nodes);
      for (k=0; k < chw.population[j].attractor_size(); k++) {
        for (l=0; l< nodes; l++) {
          matatr[k][l] = chw.population[j].attractor_element(k,l);
        }
      }
      disb = nodes*fifu.distance(phenB, matatr, chw.population[j].attractor_size(), nodes);
      disc = nodes*fifu.distance(phenC, matatr, chw.population[j].attractor_size(), nodes);
      for (k=0; k < chw.population[j].attractor_size(); k++) {
        delete [] matatr[k];
      }
      delete [] matatr; ////
      wb = pow((1 - selcoefpg), disb);
      wc = pow((1 - selcoefpg), disc);
      if (wb > wc)
        w = wb;
      else
        w = wc;
      chw.assign_w(j, w);
      chw.population[j].clear_attractor();
    }
    chw.calc_meanw();
    defaultC[i] = defc/(1.0*popsize);
    defaultB[i] = defb/(1.0*popsize);
    defaultA[i] = defa/(1.0*popsize);
    defaultO[i] = (popsize - defa - defb - defc)/(1.0*popsize);
    wmax[i] = chw.return_maxw();
    meanedges[i] = sume/(1.0*popsize);
    if ((defc > 0) && (tpC < 0)) {
      tpC = i;
      tdC = i;
    }
    if ((defc == 0) && (tpC >= 0))
      tdC = -1;
    if ((defc > 0) && (tpC >= 0))
      tdC = i;
    
    if ((defb > 0) && (tpB < 0)) {
      tpB = i;
      tdB = i;
    }
    if ((defb == 0) && (tpB >= 0))
      tdB = -1;
    if ((defb > 0) && (tpB >= 0) && (tdB < 0))
      tdB = i;
    ///
    if ((defaultC[i] >= 0.05) && (t5C < 0)) {//mod
      t5C = i;
      for (k=0; k < nodes; k++) {
        oci[k] = ci2[k];
      }
      numperaccAfromC = 0;
      numredaccAfromC =0;
      for (j=0; j < popsize; j++) {
        if (cualesdefC[j]) {
          contemp = 0;
          for (k=0; k < nodes; k++) {
            oci[k] *= (-1);
            chw.population[j].set_as_state(oci);
            chw.population[j].find_an_attractor();
            if (fifu.strict(chw.population[j], phenA) > 0) {
              contemp++;
            }
            oci[k] *= (-1);
            chw.population[j].clear_attractor();
          }
          if (contemp > 0) {
            numperaccAfromC += contemp;
            numredaccAfromC++;
          }
        }
      }
      meanaccpperAfromC[0] = (numperaccAfromC*1.0)/(defc*1.0);
      meanaccpredAfromC[0] = (numredaccAfromC*1.0)/(defc*1.0);
    }
    if ((defaultB[i] >= 0.05) && (t5B < 0)) {//mod
      t5B = i;
      for (k=0; k < nodes; k++) {
        oci[k] = ci2[k];
      }
      numperaccAfromB = 0;
      numredaccAfromB =0;
      for (j=0; j < popsize; j++) {
        if (cualesdefB[j]) {
          contemp = 0;
          for (k=0; k < nodes; k++) {
            oci[k] *= (-1);
            chw.population[j].set_as_state(oci);
            chw.population[j].find_an_attractor();
            if (fifu.strict(chw.population[j], phenA) > 0) {
              contemp++;
            }
            oci[k] *= (-1);
            chw.population[j].clear_attractor();
          }
          if (contemp > 0) {
            numperaccAfromB += contemp;
            numredaccAfromB++;
          }
        }
      }
      meanaccpperAfromB[0] = (numperaccAfromB*1.0)/(defb*1.0);
      meanaccpredAfromB[0] = (numredaccAfromB*1.0)/(defb*1.0);
    }
    //
    
    
    ///
    if ((defaultC[i] >= 0.50) && (t50C < 0)) {//mod
      t50C = i;
      for (k=0; k < nodes; k++) {
        oci[k] = ci2[k];
      }
      numperaccAfromC = 0;
      numredaccAfromC =0;
      for (j=0; j < popsize; j++) {
        if (cualesdefC[j]) {
          contemp = 0;
          for (k=0; k < nodes; k++) {
            oci[k] *= (-1);
            chw.population[j].set_as_state(oci);
            chw.population[j].find_an_attractor();
            if (fifu.strict(chw.population[j], phenA) > 0) {
              contemp++;
            }
            oci[k] *= (-1);
            chw.population[j].clear_attractor();
          }
          if (contemp > 0) {
            numperaccAfromC += contemp;
            numredaccAfromC++;
          }
        }
      }
      meanaccpperAfromC[1] = (numperaccAfromC*1.0)/(defc*1.0);
      meanaccpredAfromC[1] = (numredaccAfromC*1.0)/(defc*1.0);
    }
    if ((defaultB[i] >= 0.50) && (t50B < 0)) {//mod
      t50B = i;
      for (k=0; k < nodes; k++) {
        oci[k] = ci2[k];
      }
      numperaccAfromB = 0;
      numredaccAfromB =0;
      for (j=0; j < popsize; j++) {
        if (cualesdefB[j]) {
          contemp = 0;
          for (k=0; k < nodes; k++) {
            oci[k] *= (-1);
            chw.population[j].set_as_state(oci);
            chw.population[j].find_an_attractor();
            if (fifu.strict(chw.population[j], phenA) > 0) {
              contemp++;
            }
            oci[k] *= (-1);
            chw.population[j].clear_attractor();
          }
          if (contemp > 0) {
            numperaccAfromB += contemp;
            numredaccAfromB++;
          }
        }
      }
      meanaccpperAfromB[1] = (numperaccAfromB*1.0)/(defb*1.0);
      meanaccpredAfromB[1] = (numredaccAfromB*1.0)/(defb*1.0);
    }
    //
    ///
    if ((defaultC[i] >= 0.80) && (t80C < 0)) {//mod
      t80C = i;
      for (k=0; k < nodes; k++) {
        oci[k] = ci2[k];
      }
      numperaccAfromC = 0;
      numredaccAfromC =0;
      for (j=0; j < popsize; j++) {
        if (cualesdefC[j]) {
          contemp = 0;
          for (k=0; k < nodes; k++) {
            oci[k] *= (-1);
            chw.population[j].set_as_state(oci);
            chw.population[j].find_an_attractor();
            if (fifu.strict(chw.population[j], phenA) > 0) {
              contemp++;
            }
            oci[k] *= (-1);
            chw.population[j].clear_attractor();
          }
          if (contemp > 0) {
            numperaccAfromC += contemp;
            numredaccAfromC++;
          }
        }
      }
      meanaccpperAfromC[2] = (numperaccAfromC*1.0)/(defc*1.0);
      meanaccpredAfromC[2] = (numredaccAfromC*1.0)/(defc*1.0);
    }
    if ((defaultB[i] >= 0.80) && (t80B < 0)) {//mod
      t80B = i;
      for (k=0; k < nodes; k++) {
        oci[k] = ci2[k];
      }
      numperaccAfromB = 0;
      numredaccAfromB =0;
      for (j=0; j < popsize; j++) {
        if (cualesdefB[j]) {
          contemp = 0;
          for (k=0; k < nodes; k++) {
            oci[k] *= (-1);
            chw.population[j].set_as_state(oci);
            chw.population[j].find_an_attractor();
            if (fifu.strict(chw.population[j], phenA) > 0) {
              contemp++;
            }
            oci[k] *= (-1);
            chw.population[j].clear_attractor();
          }
          if (contemp > 0) {
            numperaccAfromB += contemp;
            numredaccAfromB++;
          }
        }
      }
      meanaccpperAfromB[2] = (numperaccAfromB*1.0)/(defb*1.0);
      meanaccpredAfromB[2] = (numredaccAfromB*1.0)/(defb*1.0);
    }
    //
    if (genfromint) {
      if ((t80B > 0) || (t80C > 0)) //mod
        break;
    } else {
      if ((t50B > 0) || (t50C > 0)) //mod
        break;
    }
    
    if ((i < (generations-1)) && (t80C < 0) && (t80B < 0)) //mod
      chw.one_generation_sex();
  }
  //checar
  
  int hasta = i;
  
  
  if (genfromint)
    cuerda = "Results/pn1/";
  else
    cuerda = "Results/sinplas/";
  
  
  
  if (hasta == generations) {
    bas.open_ofstream(fs, cuerda+"errors/"+bas.inttostring(seed)+"neverended.txt");
    fs << "never ended.\n";
    fs.close();
  }
  hasta++;
  
  bas.open_ofstream(fs, cuerda+"defaultA/"+bas.inttostring(seed)+"defaultA.txt");
  for (i=0; i< hasta; i++)
  fs << i << "\t" << defaultA[i] << endl;
  fs.close();
  
  bas.open_ofstream(fs, cuerda+"defaultB/"+bas.inttostring(seed)+"defaultB.txt");
  for (i=0; i< hasta; i++)
  fs << i << "\t" << defaultB[i] << endl;
  fs.close();
  
  bas.open_ofstream(fs, cuerda+"defaultC/"+bas.inttostring(seed)+"defaultC.txt");
  for (i=0; i< hasta; i++)
  fs << i << "\t" << defaultC[i] << endl;
  fs.close();
  
  bas.open_ofstream(fs, cuerda+"defaultO/"+bas.inttostring(seed)+"defaultO.txt");
  for (i=0; i< hasta; i++)
  fs << i << "\t" << defaultO[i] << endl;
  fs.close();
  
  
  
  bas.open_ofstream(fs, cuerda+"wmax/"+bas.inttostring(seed)+"wmax.txt");
  for (i=0; i< hasta; i++)
  fs << i << "\t" << wmax[i] << endl;
  fs.close();
  
  bas.open_ofstream(fs, cuerda+"meanedges/"+bas.inttostring(seed)+"meanedges.txt");
  for (i=0; i< hasta; i++)
  fs << i << "\t" << meanedges[i] << endl;
  fs.close();
  
  bas.open_ofstream(fs, cuerda+"tpC/"+bas.inttostring(seed)+"tpC.txt");
  fs << tpC << endl;
  fs.close();
  bas.open_ofstream(fs, cuerda+"tdC/"+bas.inttostring(seed)+"tdC.txt");
  fs << tdC << endl;
  fs.close();
  bas.open_ofstream(fs, cuerda+"t50C/"+bas.inttostring(seed)+"t50C.txt");
  fs << t50C << endl;
  fs.close();
  
  bas.open_ofstream(fs, cuerda+"tpB/"+bas.inttostring(seed)+"tpB.txt");
  fs << tpB << endl;
  fs.close();
  bas.open_ofstream(fs, cuerda+"tdB/"+bas.inttostring(seed)+"tdB.txt");
  fs << tdB << endl;
  fs.close();
  bas.open_ofstream(fs, cuerda+"t50B/"+bas.inttostring(seed)+"t50B.txt");
  fs << t50B << endl;
  fs.close();
  
  
  bas.open_ofstream(fs, cuerda+"t5C/"+bas.inttostring(seed)+"t5C.txt");
  fs << t5C << endl;
  fs.close();
  bas.open_ofstream(fs, cuerda+"t5B/"+bas.inttostring(seed)+"t5B.txt");
  fs << t5B << endl;
  fs.close();
  
  bas.open_ofstream(fs, cuerda+"t80C/"+bas.inttostring(seed)+"t80C.txt");
  fs << t80C << endl;
  fs.close();
  bas.open_ofstream(fs, cuerda+"t80B/"+bas.inttostring(seed)+"t80B.txt");
  fs << t80B << endl;
  fs.close();
  
  bas.open_ofstream(fs, cuerda+"meanaccpperAfromB/"+bas.inttostring(seed)+".txt");
  fs << 5 << "\t" << meanaccpperAfromB[0] << endl;
  fs << 50 << "\t" << meanaccpperAfromB[1] << endl;
  fs << 80 << "\t" << meanaccpperAfromB[2] << endl;
  fs.close();
  
  bas.open_ofstream(fs, cuerda+"meanaccpredAfromB/"+bas.inttostring(seed)+".txt");
  fs << 5 << "\t" << meanaccpredAfromB[0] << endl;
  fs << 50 << "\t" << meanaccpredAfromB[1] << endl;
  fs << 80 << "\t" << meanaccpredAfromB[2] << endl;
  fs.close();
  
  bas.open_ofstream(fs, cuerda+"meanaccpperAfromC/"+bas.inttostring(seed)+".txt");
  fs << 5 << "\t" << meanaccpperAfromC[0] << endl;
  fs << 50 << "\t" << meanaccpperAfromC[1] << endl;
  fs << 80 << "\t" << meanaccpperAfromC[2] << endl;
  fs.close();
  
  bas.open_ofstream(fs, cuerda+"meanaccpredAfromC/"+bas.inttostring(seed)+".txt");
  fs << 5 << "\t" << meanaccpredAfromC[0] << endl;
  fs << 50 << "\t" << meanaccpredAfromC[1] << endl;
  fs << 80 << "\t" << meanaccpredAfromC[2] << endl;
  fs.close();
  
  delete [] phenA;
  delete [] phenB;
  delete [] phenC;
  delete [] ci;
  delete [] oci;
  delete [] defaultA;
  delete [] defaultB;
  delete [] defaultC;
  delete [] defaultO;
  delete [] wmax;
  delete [] meanedges;
  
  jacta.close_rng();
  cout << 1 << endl;
  return 0;
}


