/*****************************************************************************************
  FileName     [ ham_cycle.cpp ]
  PackageName  [ Hamiltonian-cycle-with-SAT ]
  Synopsis     [ Define function ]
  Author       [ Hugo, Chen ]
  Professor    [ Ric Huang ]
  Copyright    [ Copyleft(c) 2021, NTUEE, Taiwan ]
*****************************************************************************************/
//----------------------------------------------------------------------------------------
//    Include 
//----------------------------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <algorithm>
#include <utility>
#include <vector>
#include <map>
#include "sat.h"
#include "ham_cycle.h"
// #include "../sat/sat.h"

using namespace std;

//----------------------------------------------------------------------------------------
//    Global Variable
//----------------------------------------------------------------------------------------
// constructor //
ham_cycle::ham_cycle() {}
int numGates;

//----------------------------------------------------------------------------------------
//    Main Function
//----------------------------------------------------------------------------------------
// read file //
bool
ham_cycle::read_input(const string& infile)
{
  // For MobaXTerm, I need to use 'c_str()' to make it successfully compile
  fstream fin(infile.c_str()); 
  // store number of nodes/edges //
  fin >> V_num;
  fin >> E_num;
  // store edges //
  edges.assign(E_num, Edge());
  for (int i = 0; i < E_num; i++) 
  { 
    fin >> edges[i].u >> edges[i].v; 
    // Calculate number of each node
    if (count_node.count(edges[i].u) == 0) { count_node[edges[i].u] = 1; }
    else if (count_node.count(edges[i].u) == 1) { count_node[edges[i].u] = count_node[edges[i].u] + 1; }
    if (count_node.count(edges[i].v) == 0) { count_node[edges[i].v] = 1; }
    else if (count_node.count(edges[i].v) == 1) { count_node[edges[i].v] = count_node[edges[i].v] + 1; }
  }
  fin.close();
  return true;
}

void 
ham_cycle::initCircuit()
{  
   // Init gates, send by const referrences(read only but no copy)
   // V_num: 每個node input要設置一個gate variable
   // 2*E_num: 每個邊的兩個端點要在cycle上 (edge = a <-> b = XNOR = XOR + INV = (ab')'*(a'b)')
   // 2*E_num-1: 每個clause的bar互相AND, 兩兩一組, 這邊作法太暴力, 耗memory
   // 3: 前面預處理每個node一定要出現偶數次, 這個gate直接為0, 若到最後都是偶數, 這個gate設為1, 再跟前面結果輸出
   numGates = V_num + 2*E_num + (2*E_num-1) + 3;
   for(int i = 0 ; i < numGates ; i++) { gates.push_back(new Gate(i)); }
}

void
ham_cycle::HC_to_CNF(SatSolver& s)
{
  // Allocate and record variables; No Var ID for POs
  for(size_t i = 0, n = gates.size(); i < n; ++i) 
  {
    Var v = s.newVar();
    gates[i]->setVar(v);
  }
  // Constraint: 每個邊的兩個端點要在cycle上
  // 2*E_num: 每個邊的兩個端點要在cycle上 (edge = a <-> b = XNOR = XOR + INV = (ab')'*(a'b)') 
  int stage_1 = V_num;
  for(int i = 0 ; i < E_num ; i++)
  {
    s.addAigCNF(gates[2*i+stage_1]->getVar(), gates[2*i]->getVar(), true, gates[2*i+1]->getVar(), false);
    s.addAigCNF(gates[2*i+1+stage_1]->getVar(), gates[2*i]->getVar(), false, gates[2*i+1]->getVar(), true);
  }
  // 2*E_num-1: 兩兩邊一組互相 AND (a' & b' & c' = ((a' & b') & c'))
  int stage_2 = V_num + 2*E_num;
  for(int i = 0 ; i < 2*E_num-1 ; i++)
  {
    s.addAigCNF(gates[i+stage_2]->getVar(), gates[2*i+stage_1]->getVar(), true, gates[2*i+1+stage_1]->getVar(), true);
  }
  // Constraint: 每個node只能出現偶數次
  // 3: 前面預處理每個node一定要出現偶數次, 這個gate直接為0, 若到最後都是偶數, 這個gate設為1, 再跟前面結果輸出
  int mark = 0, stage_2_parallel = V_num + 2*E_num + (2*E_num-1);
  map<int, int>::iterator it;
	for(it = count_node.begin() ; it != count_node.end() ; it++)
  {
    // cout << "node: " << it->first << " --> " << it->second << endl;
    if (it->second % 2 == 1)  // gates[it->first+stage_2_parallel] = const 0
    { 
      cout << "Failed.........(odd node)" << endl;
      s.addAigCNF(gates[stage_2_parallel+1]->getVar(), gates[stage_2_parallel]->getVar(), true, gates[stage_2_parallel]->getVar(), false);
      s.addAigCNF(gates[stage_2_parallel+2]->getVar(), gates[stage_2_parallel+1]->getVar(), false, gates[stage_2_parallel-1]->getVar(), false);
      mark = 1;
      break; 
    } 
  }
  if (mark == 0) 
  { 
    cout << "Success!!!!!!!!!! (even node)" << endl; 
    s.addAigCNF(gates[stage_2_parallel+1]->getVar(), gates[stage_2_parallel]->getVar(), true, gates[stage_2_parallel]->getVar(), false);
    s.addAigCNF(gates[stage_2_parallel+2]->getVar(), gates[stage_2_parallel+1]->getVar(), true, gates[stage_2_parallel-1]->getVar(), false);
  }
  return;
}

void 
ham_cycle::report_result(const SatSolver& solver, bool result)
{
  solver.printStats();
  cout << (result? "SAT" : "UNSAT") << endl;
  if (result) {
    for (size_t i = 0, n = gates.size(); i < n; ++i)
        cout << solver.getValue(gates[i]->getVar()) << endl;
  }
}