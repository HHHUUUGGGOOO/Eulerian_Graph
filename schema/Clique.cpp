#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <cstring>
#include "sat.h"
#include <ctime>
#include "Graph.h"
//why the running time is proportional to k value?
//i.e when k large , it should be easy to be unsat , but the running time still increase.
using namespace std;
//(a1*...*an)->f can be written as (~a1->~f)(~a2->~f)...(~an->f)(a1...an->f)
//create clauses (a1+~f)(a2+~f)...(an+~f)(~a1+~a2+...+~an+f)

class Gate{
public:
   Gate(unsigned i = 0): _gid(i) {}
   ~Gate() {}

   Var getVar() const { return _var; }
   void setVar(const Var& v) { _var = v; }

private:
   unsigned   _gid;  // for debugging purpose...
   Var        _var;
};

vector<Gate *> gates;
int numGates;

void initCircuit(int k , Graph* g){ //send by const referrences(read only but no copy)
   // Init gates
   int V = g->numofV();
   int E = g->numofE();
   numGates = 2 * V * k + 2*(k*(k - 1)/2)*(pow(V , 2) - (E * 2)) - 1;
   for(int i = 0 ; i < numGates ; i++) gates.push_back(new Gate(i));
}

void genProofModel(SatSolver& s , int k , Graph* g){
   // Allocate and record variables; No Var ID for POs
   int V = g->numofV();
   int E = g->numofE();
   vector<vector<edge>> adjL = g->findAdjL();
   for (size_t i = 0, n = gates.size(); i < n; ++i) {
      Var v = s.newVar();
      gates[i]->setVar(v);
   }

   int temp = V * k - 1;
   for(int i = 0 ; i < k ; i++){ //can modify addAigCNF to acommodate multiple input AND graph (send vectors)
      temp += 1;
      s.addAigCNF(gates[temp]->getVar() , gates[V*i]->getVar() , true , gates[V*i + 1]->getVar() , true);
      for(int j = 2 ; j <= (V-1) ; j++){
         temp += 1;
         s.addAigCNF(gates[temp]->getVar() , gates[temp - 1]->getVar() , false , gates[V*i + j]->getVar() , true);
      }
   }
   int stage1 = temp;
   for(int i = 0 ; i < V ; i++){
      for(int j = 0 ; j < k ; j++){
         for(int m = j + 1 ; m < k ; m++){
            temp += 1;
            s.addAigCNF(gates[temp]->getVar() , gates[V*j + i]->getVar() , false , gates[V*m + i]->getVar() , false);
         }
      }
   }
   for(int i = 0 ; i < adjL.size() ; i++){
      for(int j = 0 ; j < adjL[i].size() ; j++){
         int st = adjL[i][j].start->index;
         int en = adjL[i][j].end->index;
         for(int m = 0 ; m < k ; m++){
            for(int n = m + 1 ; n < k ; n++){
               temp += 1;
               s.addAigCNF(gates[temp]->getVar() , gates[V*m + st]->getVar() , false , gates[V*n + en]->getVar() , false);
               temp += 1;
               s.addAigCNF(gates[temp]->getVar() , gates[V*m + en]->getVar() , false , gates[V*n + st]->getVar() , false);
            }
         }
      }
   }
   temp += 1;
   s.addAigCNF(gates[temp]->getVar() , gates[stage1 - (k - 1)*(V - 1)]->getVar() , true , gates[stage1 - (k - 2)*(V - 1)]->getVar() , true);
   for(int i = 3 ; i <= k ; i++){
      temp += 1;
      s.addAigCNF(gates[temp]->getVar() , gates[temp - 1]->getVar() , false , gates[stage1 - (k - i)*(V - 1)]->getVar() , true);
   }
   for(int i = 1 ; i < (k*(k - 1)/2)*(pow(V , 2) - (E * 2)) + 1 ; i++){
      temp += 1;
      s.addAigCNF(gates[temp]->getVar() , gates[temp-1]->getVar() , false , gates[stage1 + i]->getVar() , true);
   }
   //cout << "finaltemp = " << temp << endl;
   cout << "model complete!" << endl;
}

void reportResult(const SatSolver& solver, bool result , int V , int k){
   solver.printStats();
   cout << (result? "SAT" : "UNSAT") << endl;
   if(result){
      int counter = 0;
      for(size_t i = 0 ; i < V*k; ++i){
         counter++;
         cout << solver.getValue(gates[i]->getVar()) << " ";
         if(counter%V == 0) cout << endl;
      }
   }
   cout << endl;
}

int main(int argc , char* argv[]){
   clock_t start , end;
   start = clock();
   fstream fin(argv[1]);
   int n , m;
   fin >> n >> m;
   Graph g;
   g.initialize(n , m);
	int a , b , c;
   while (fin >> a >> b){
      g.buildNonEdge(a , b);
   }

   int x;
   cin >> x;
   int V = g.numofV();
   int E = g.numofE();

   initCircuit(x , &g);

   SatSolver solver;
   solver.initialize();

   //
   genProofModel(solver , x , &g);

   bool result;
   solver.assumeProperty(gates[2 * V * x + 2*(x*(x - 1)/2)*(pow(V , 2) - (E * 2)) - 2]->getVar() , true);
   result = solver.assumpSolve();
   reportResult(solver, result , V , x);
   end = clock();
   double time = end - start;
   cout << "time = " << time/CLOCKS_PER_SEC << endl;
}