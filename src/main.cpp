/*****************************************************************************************
  FileName     [ main.cpp ]
  PackageName  [ Hamiltonian-cycle-with-SAT ]
  Synopsis     [ Main function ]
  Author       [ Hugo, Chen ]
  Professor    [ Ric Huang ]
  Copyright    [ Copyleft(c) 2021, NTUEE, Taiwan ]
*****************************************************************************************/
//----------------------------------------------------------------------------------------
//    Include 
//----------------------------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <vector>
#include "sat.h"
#include "ham_cycle.h"
// #include "sat/sat.h"

using namespace std;

//----------------------------------------------------------------------------------------
//    Global Variable / Function
//----------------------------------------------------------------------------------------
ham_cycle* hamcycle = new ham_cycle();

//----------------------------------------------------------------------------------------
//    Main Function
//----------------------------------------------------------------------------------------
int main(int argc, char** argv)
{
    // Calculate operation time
    clock_t start , end;
    start = clock();
    // Read in file (Read input file here, argc = 3 --> bin/input/output)
    if (argc == 3) 
    {  
      if (!hamcycle->read_input(argv[1])) 
      {
        cerr << "Error: cannot open input file: " << argv[1] << endl;
        return 0;
      }
    }
    else if (argc != 3) 
    { 
      cerr << "Error: illegal number of argument (" << argc << ") " << endl;
      return 0;
    }
    // Transform input graph into CNF (HC_to_CNF())
    hamcycle->initCircuit();
    SatSolver solver;
    solver.initialize();
    // Do SAT solver (include sat.h and create a function Do_SAT())
    hamcycle->HC_to_CNF(solver);
    // Print whether satisfied and output the corresponding file (Print_HC_Result())
    bool result;
    // Var newV = solver.newVar();  // create a new variable
    // solver.addXorCNF(newV, hamcycle->gates[hamcycle->gates.size()-1]->getVar(), false, hamcycle->gates[hamcycle->gates.size()-2]->getVar(), false); // true means "inverted", SAT gate
    // solver.assumeRelease();  // Clear assumptions
    // solver.assumeProperty(newV , true);
    solver.assumeProperty(hamcycle->gates[hamcycle->gates.size()-1]->getVar() , true); // true: 結果為1叫做SAT
    result = solver.assumpSolve();
    hamcycle->report_result(solver, result); // 最後必須是1
    // Calculate end time
    end = clock();
    double time = end - start;
    cout << "time = " << time/CLOCKS_PER_SEC << endl;
    return 0;
}
