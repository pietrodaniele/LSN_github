#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "random.h"
#include "funzioni.h"
#include "MC_NVT.h"

using namespace std;

// defining objects
Random rnd;
Mc_nvt nvt;

int main(int argc, char const *argv[]){
  if(argc!=2){
    cerr << "./file.exe state_of_the_system(1-->gas;2-->liquid;3-->solid)" << endl;
    exit(1);
  }
  int state_of_the_system = atoi(argv[1]);
  nvt.Input(rnd,state_of_the_system); // carica i file con i parametri e ed esegue il ciclo di equilibrazione
  for(int i=0; i<nvt.Get_nblk(); i++){
    for(int j=0; j<nvt.Get_nstep(); j++){
      nvt.Move(rnd); // muove il sistema con accettazione di metropolis
      nvt.Measure(); // misura U ista e pres ista e accumula sui vari hist(block,finale)

    }
    cout << "----- Block number " << i+1 << " -----" << endl;
    cout << "Accettazione = " << nvt.Get_acc()/(double)nvt.Get_att() << endl << endl;
    nvt.Print_hist_block(); // stampa hist per ogni singolo blocco in un unico file, su python verranno differenziati
    nvt.Clean_hist_block(); // pulisce hist_block alla fine di ogni blocco
  }
  //nvt.Print_hist_block();
  cout << "Calculating final average value of 𝑔(𝑟) with uncert. ..." << endl;
  nvt.Final_g_err(); // stampa hist finale creato accumulando tutti i vari blocchi e le incertezze
return 0;
}
