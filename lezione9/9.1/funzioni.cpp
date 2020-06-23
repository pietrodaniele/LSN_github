#include "funzioni.h"
#include "random.h"
#include <string>

void createrandom(Random& rnd){
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;
   rnd.SaveSeed();
};

void Cities_config(Random& rnd, int config, vector<double>& x, vector<double>& y){
  ofstream Write;
  if(config==0){
    Write.open("data/Cities_config_Circ.0");
    for(int i=0; i<32; i++){
      double teta = rnd.Rannyu(0,2*M_PI);
      x.push_back(cos(teta));
      y.push_back(sin(teta));
      Write << i+1 << " " << x[i] << " " << y[i] << endl;
    }
  }else if(config==1){
    Write.open("data/Cities_config_Square.0");
    for(int i=0; i<32; i++){
      x.push_back(rnd.Rannyu());
      y.push_back(rnd.Rannyu());
      Write << i+1 << " " << x[i] << " " << y[i] << endl;
    }
  }
  Write.close();
};

/*void New_Generation(Population* population, Population* newgen){
  Random rnd;
  for(int i=0; i<50; i++){
    //population->Set_an_Individual(i,newgen->Get_an_Individual(i));
    Individual ind = newgen->Get_an_Individual(i);
    population->Set_an_Individual(i,ind);
    cout<<i<<endl;
  }
  population->In_Fitness_order();
};*/

bool Is_New_Son(Individual Prova, vector<Individual>& newgen){
  bool Son=true;
  int diversity=0;
  for(int i=0; i<newgen.size(); i++){
    for(int j=0; j<32; j++){
      if(Prova.Get_gene(j)==newgen[i].Get_gene(j)){
        diversity++;
      }
    }
    if(diversity==32){
      Son=false;
      return Son;
    }
  }
  return Son;
}

void Print(Population* pop, vector<double>& x, vector<double>& y){
  for(int i=0; i<pop->Get_Number_ind(); i++){
    Individual ind = pop->Get_an_Individual(i);
    cout << ind.Calc_Fitness2_oncirc(x,y) << ' ';
    for(int j=0; j<ind.Get_Number_gene(); j++){
      cout << ind.Get_gene(j) << ' ';
    }
    cout << endl;
  }


};

Individual Mutation(Population* pop, Random& rnd, vector<double>& x, vector<double>& y, double prop_mut){
  Individual Son1(rnd);
  Individual Son2(rnd);
  double r = rnd.Rannyu();
  if(r<prop_mut){
    pop->Mutation1(rnd,Son1);
  }else if (prop_mut<=r && r<2*prop_mut){
    pop->Mutation2(rnd,Son1);
  }else if (2*prop_mut<=r && r<3*prop_mut){
    pop->Mutation3(rnd,Son1);
  }else if (3*prop_mut<=r && r<4*prop_mut){
    pop->Mutation3(rnd,Son1);
  }else if (r>=4*prop_mut){
    pop->Crossover(rnd,Son1,Son2);
    if(Son1.Calc_Fitness2_oncirc(x,y)>Son2.Calc_Fitness2_oncirc(x,y)){
      Son1=Son2;
    }
  }
  return Son1;
};
/*Individual Mutation(Population* pop, Random& rnd, vector<double>& x, vector<double>& y, double prop_mut){
  Individual Son1(rnd);
  Individual Son2(rnd);
  double r = rnd.Rannyu();
  pop->Crossover(rnd,Son1,Son2);
  return Son1;
};*/
