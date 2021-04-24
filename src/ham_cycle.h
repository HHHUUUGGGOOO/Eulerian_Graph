/*****************************************************************************************
  FileName     [ ham_cycle.h ]
  PackageName  [ Hamiltonian-cycle-with-SAT ]
  Synopsis     [ header file ]
  Author       [ Hugo, Chen ]
  Professor    [ Ric Huang ]
  Copyright    [ Copyleft(c) 2021, NTUEE, Taiwan ]
*****************************************************************************************/
//----------------------------------------------------------------------------------------
//    Include 
//----------------------------------------------------------------------------------------
#ifndef HAM_CYCLE_H
#define HAM_CYCLE_H

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <utility>
#include <vector>
#include <map>
#include "sat.h"

using namespace std;

//----------------------------------------------------------------------------------------
//    Global Variable
//----------------------------------------------------------------------------------------
class ham_cycle;
extern ham_cycle* hamcycle;

//----------------------------------------------------------------------------------------
//    Main Function
//----------------------------------------------------------------------------------------
// Edge //
struct Edge
{
    // u: begin node, v: end node //
    int u, v;
};

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

class ham_cycle
{
  public:
    friend class Gate;
    ham_cycle();
    virtual ~ham_cycle() {}

    // read input file --> call in main.cpp //
    bool                  read_input(const string&);
    void                  initCircuit();
    void                  HC_to_CNF(SatSolver& s);
    void                  report_result(const SatSolver&, bool result);

    // variables //
    int                   V_num;
    int                   E_num;
    map<int, int>         count_node;
    vector<Edge>          edges; 
    vector<Gate *>        gates;

};

#endif