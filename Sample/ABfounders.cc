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


using std::cout;
using std::cin;
using std::getline;
using std::string;
using std::endl;
using std::ofstream;
using std::ifstream;
using namespace std;

int main(int argc, char **argv) {
  ifstream fe;
  ofstream fs;

  int seed, nuredes, nodes, edges, redesxarch, numarch, mame, mult, dcia, dcib, dab;
  seed = atoi(argv[1]);
  nuredes = atoi(argv[2]); //por job 20
  nodes = atoi(argv[3]);
  edges = atoi(argv[4]);
  redesxarch = atoi(argv[5]); //20
  mame = atoi(argv[6]); //2
  mult = atoi(argv[7]); //20
  dcia = atoi(argv[8]); //3
  dcib = atoi(argv[9]); //3
  dab = atoi(argv[10]); //6
  
    int dcic, dbc, dac;
  dcic = dcib;//3;
  dbc = dab;//6;
  dac = dab;//6;
  
  if ((nuredes%redesxarch)!=0) {
    cout << "[Error]: The number of networks is not a multiple integer of the number of networks per file.\n";
    exit(1);
  }
  numarch = nuredes/redesxarch;

  Alea jacta(seed+34);
  Basics bas(jacta);
  int i;
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

  GraphI co(jacta);
  co.rwalk_2fp2t(nodes, edges, phenA, phenB, ci, ci2, nuredes, "networks/", mame, mult, true, numarch, ((seed-1)*numarch)+1);
  
  delete [] ci;
  delete [] ci2;
  delete [] phenB;
  delete [] phenA;
  jacta.close_rng();
  return 0;
}
