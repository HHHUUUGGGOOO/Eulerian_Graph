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
#include <cmath>
#include <utility>
#include <vector>
#include <algorithm>
#include <map>
#include "sat.h"

using namespace std;

//----------------------------------------------------------------------------------------
//    Class
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
    void                  initCircuit();
    void                  HC_to_CNF(SatSolver& s, vector<int>* start, vector<int>* end);
    void                  report_result(const SatSolver&, bool result);

    // variables //
    int                   V_num;
    int                   E_num;
    int                   gate_use_size;
    vector<int>           path_index;
    vector<int>           edges_u;
    vector<int>           edges_v;
    vector<int>           selected_e; 
    vector<int>           multi_cycle_node;
    vector<Edge>          edges; 
    vector<Gate *>        gates;
    map<int, int>         count_node;

};

//----------------------------------------------------------------------------------------
//    Global Variable / Function
//----------------------------------------------------------------------------------------
ham_cycle::ham_cycle() {}
ham_cycle* hamcycle = new ham_cycle();
int numGates;

//----------------------------------------------------------------------------------------
//    Declare Function
//----------------------------------------------------------------------------------------

void 
ham_cycle::initCircuit()
{  
   // Initialize new gates
   numGates = 2*E_num + E_num + (E_num-1) + 2*(pow(V_num, 3) - pow(V_num, 2)) + 1 + 2*(pow(V_num, 2) - V_num);
   for(int i = 0 ; i < numGates ; i++) { gates.push_back(new Gate(i)); }
   // debug
   cout << "# of gates = " << numGates << endl;
   cout << endl;
}

void
ham_cycle::HC_to_CNF(SatSolver& s, vector<int>* start, vector<int>* end)
{
// Allocate and record variables; No Var ID for POs
  for(size_t i = 0, n = gates.size(); i < n; ++i) 
  {
    Var v = s.newVar();
    gates[i]->setVar(v);
  }
// input: 2*E_num ( original & inverse edge )
  int stage_xor = 2*E_num;
  for(int i = 0 ; i < E_num ; i++)
  {
    s.addXorCNF(gates[i+stage_xor]->getVar(), gates[2*i]->getVar(), false, gates[2*i+1]->getVar(), false);
  }
// xor -> and all the xor ( result = gates[stage_xor_and + E_num - 2] )
  int stage_xor_and = stage_xor + E_num;
  s.addAigCNF(gates[stage_xor_and]->getVar(), gates[stage_xor]->getVar(), false, gates[stage_xor+1]->getVar(), false);
  for(int i = 2 ; i < E_num ; i++)
  {
    s.addAigCNF(gates[stage_xor_and+i-1]->getVar(), gates[stage_xor_and+i-2]->getVar(), false, gates[stage_xor+i]->getVar(), false);
  }
// start edges: <a, b, c> --> (a+b+c)(a'+b')(a'+c')(b'+c') --> (a'+b') = (ab)'
  int count_in = 0;
  int stage_start = stage_xor_and + E_num - 1;
  for (int i = 0 ; i < V_num ; i++)
  {
    for (int j = 0 ; j < start[i].size() ; j++)
    {
      if (start[i].size() > 1)
      {
        for (int k = j+1 ; k < start[i].size() ; k++)
        {
          s.addAigCNF(gates[stage_start+count_in]->getVar(), gates[start[i][j]]->getVar(), false, gates[start[i][k]]->getVar(), false);
          ++count_in;
        }
      }
      else if (start[i].size() == 1)
      {
        // const_1 = (a*a')'
        s.addAigCNF(gates[stage_start+count_in]->getVar(), gates[start[i][j]]->getVar(), true, gates[start[i][j]]->getVar(), false);
        ++count_in;
      }
    }
  }
// end edges: <a, b, c> --> (a+b+c)(a'+b')(a'+c')(b'+c') --> (a'+b') = (ab)'
  int count_out = 0;
  int stage_end = stage_start + count_in;
  for (int i = 0 ; i < V_num ; i++)
  {
    for (int j = 0 ; j < end[i].size() ; j++)
    {
      if (end[i].size() > 1)
      {
        for (int k = j+1 ; k < end[i].size() ; k++)
        {
          s.addAigCNF(gates[stage_end+count_out]->getVar(), gates[end[i][j]]->getVar(), false, gates[end[i][k]]->getVar(), false);
          ++count_out;
        }
      }
      else if (end[i].size() == 1)
      {
        s.addAigCNF(gates[stage_end+count_out]->getVar(), gates[end[i][j]]->getVar(), true, gates[end[i][j]]->getVar(), false);
        ++count_out;
      }
    }
  }
  
// start edges && end edges ( result = gates[stage_s_e_and + count_in + count_out - 3] )
  int stage_s_e_and = stage_end + count_out;
  s.addAigCNF(gates[stage_s_e_and]->getVar(), gates[stage_start]->getVar(), true, gates[stage_start+1]->getVar(), true);
  for(int i = 2 ; i < count_in + count_out - 1 ; i++)
  {
    s.addAigCNF(gates[stage_s_e_and+i-1]->getVar(), gates[stage_start+i]->getVar(), true, gates[stage_s_e_and+i-2]->getVar(), false);
  }

// stage_xor_and && stage_s_e_and ( result = gates[stage_xor_s_e_and] )
  int stage_xor_s_e_and = stage_s_e_and + count_in + count_out - 2;
  s.addAigCNF(gates[stage_xor_s_e_and]->getVar(), gates[stage_s_e_and+count_in+count_out-3]->getVar(), false, gates[stage_xor_and+E_num-2]->getVar(), false);

// start edges "OR all variable" term 
  int count_or_start = 0;
  int stage_or_start = stage_xor_s_e_and + 1;
  // store every "OR all" term's result gate index
  vector<int> result_gate;
  for (int i = 0 ; i < V_num ; i++)
  {
    if (start[i].size() > 1)
    {
      s.addAigCNF(gates[stage_or_start+count_or_start]->getVar(), gates[start[i][0]]->getVar(), true, gates[start[i][1]]->getVar(), true);
      ++count_or_start;
      if (start[i].size() > 2)
      {
        for (int j = 2 ; j < start[i].size() ; j++)
        {
          s.addAigCNF(gates[stage_or_start+count_or_start]->getVar(), gates[start[i][j]]->getVar(), true, gates[stage_or_start+count_or_start-1]->getVar(), false);
          ++count_or_start;
        }
      }
    }
  }
  result_gate.push_back(stage_or_start + count_or_start - 1); // 要inverse
// end edges "OR all variable" term 
  int count_or_end = 0;
  int stage_or_end = stage_or_start + count_or_start;
  for (int i = 0 ; i < V_num ; i++)
  {
    if (end[i].size() > 1)
    {
      s.addAigCNF(gates[stage_or_end+count_or_end]->getVar(), gates[end[i][0]]->getVar(), true, gates[end[i][1]]->getVar(), true);
      ++count_or_end;
      if (end[i].size() > 2)
      {
        for (int j = 2 ; j < end[i].size() ; j++)
        {
          s.addAigCNF(gates[stage_or_end+count_or_end]->getVar(), gates[end[i][j]]->getVar(), true, gates[stage_or_end+count_or_end-1]->getVar(), false);
          ++count_or_end;
        }
      }
    }
  }
  result_gate.push_back(stage_or_end + count_or_end - 1);
// AND all result
  int stage_and_all = stage_or_end + count_or_end;
  // AND previous result with the first OR_all gate
  s.addAigCNF(gates[stage_and_all]->getVar(), gates[stage_xor_s_e_and]->getVar(), false, gates[result_gate[0]]->getVar(), true);
  s.addAigCNF(gates[stage_and_all+1]->getVar(), gates[stage_and_all]->getVar(), false, gates[result_gate[1]]->getVar(), true);
// Must have "even" node
  int mark = 0, stage_even_node = stage_and_all + 2;
  map<int, int>::iterator it;
	for(it = count_node.begin() ; it != count_node.end() ; it++)
  {
    if (it->second % 2 == 1)  
    { 
      cout << "Failed.........(odd node)" << endl;
      s.addAigCNF(gates[stage_even_node+1]->getVar(), gates[stage_even_node]->getVar(), true, gates[stage_even_node]->getVar(), false);
      s.addAigCNF(gates[stage_even_node+2]->getVar(), gates[stage_even_node+1]->getVar(), false, gates[stage_and_all+1]->getVar(), false);
      mark = 1;
      break; 
    } 
  }
  if (mark == 0) 
  { 
    cout << "Success!!!!!!!!!! (even node)" << endl; 
    s.addAigCNF(gates[stage_even_node+1]->getVar(), gates[stage_even_node]->getVar(), true, gates[stage_even_node]->getVar(), false);
    s.addAigCNF(gates[stage_even_node+2]->getVar(), gates[stage_even_node+1]->getVar(), true, gates[stage_and_all+1]->getVar(), false);
  }

  gate_use_size = stage_even_node + 2;
  return;
}

void 
ham_cycle::report_result(const SatSolver& solver, bool result)
{
  solver.printStats();
  bool multi_cycle = false;
  if (!result) 
  { 
    cout << "======" << endl;
    cout << "Result" << endl;
    cout << "======" << endl;
    cout << "UNSAT (odd node)" << endl; 
    cout << "Failed......It's not a Hamiltonian cycle !" << endl;
    cout << endl;
    return; 
  }
  else if (result) {
  // debug for only edges
    cout << "===========" << endl;
    cout << "Edge result" << endl;
    cout << "===========" << endl;
    for (size_t i = 0, n = 2*E_num; i < n; ++i)
        cout << "Edge[" << i << "] -> " << solver.getValue(gates[i]->getVar()) << endl;
    cout << endl;
  // debug for all result
    /*
    cout << "==========" << endl;
    cout << "All result" << endl;
    cout << "==========" << endl;
    for (size_t i = 0, n = gate_use_size; i < n; ++i)
        cout << "gate[" << i << "] -> " << solver.getValue(gates[i]->getVar()) << endl;
    cout << endl;
    */
  // Print the path, start at the first edge
    for (int i = 0 ; i < 2*E_num ; i++)
    {
      if (solver.getValue(gates[i]->getVar()) == 1) 
      { 
        selected_e.push_back(i); 
        edges_u.push_back(edges[i].u);
        edges_v.push_back(edges[i].v);
      }
    }
    // Choose the first edge to be the start point
    path_index.push_back(selected_e[0]);
    int start_point = edges[selected_e[0]].v;
    for (int i = 0 ; i < selected_e.size()-1 ; i++)
    {
      vector<int>::iterator it = std::find(edges_u.begin(), edges_u.end(), start_point); 
      if (it != edges_u.end())
      {
        // if forms a loop, break
        if (edges[selected_e[std::distance(edges_u.begin(), it)]].v == edges[path_index[0]].u)
        {
          path_index.push_back(selected_e[std::distance(edges_u.begin(), it)]);
          start_point = edges[selected_e[std::distance(edges_u.begin(), it)]].v;
          break;
        }
        else
        {
          path_index.push_back(selected_e[std::distance(edges_u.begin(), it)]);
          start_point = edges[selected_e[std::distance(edges_u.begin(), it)]].v;
        }
      }
      else { break; }
    }
    cout << "====" << endl;
    cout << "Path" << endl;
    cout << "====" << endl;
    for (int i = 0 ; i < path_index.size() ; i++)
    {
      cout << edges[path_index[i]].u << " -> " << edges[path_index[i]].v << endl;
      multi_cycle_node.push_back(edges[path_index[i]].u);
    }
    cout << endl;
  // Checking if there is multi-cycle
    if (multi_cycle_node.size() < V_num) { multi_cycle = true; }
    else if (multi_cycle_node.size() == V_num) { multi_cycle = false; }
  // Print final result
    cout << "======" << endl;
    cout << "Result" << endl;
    cout << "======" << endl;
    if (multi_cycle) { cout << "UNSAT (multi-cycle)" << endl; cout << "Failed......It's not a Hamiltonian cycle !" << endl; }
    else if (!multi_cycle) { cout << "SAT" << endl; cout << "Success!!! It's a Hamiltonian cycle !" << endl; }
    cout << endl;
    return;
  }
}

//////////////////////////////////////////////////////////////////////////////////////////

//----------------------------------------------------------------------------------------
//    Main Function
//----------------------------------------------------------------------------------------
int main(int argc, char** argv)
{
  // Calculate operation time
    clock_t start , end;
    start = clock();
  // Read in file (Read input file here, argc = 3 --> bin/input/output)
    fstream fin(argv[1]); 
    // store number of nodes/edges //
    fin >> hamcycle->V_num;
    fin >> hamcycle->E_num;
    // s[0] = <0, 2> --> node_0 為起點的邊為 edge_0, edge_2
    vector<int> s[hamcycle->V_num];
    vector<int> e[hamcycle->V_num];
    int u, v;
    // default count node
    for (int i = 0 ; i < hamcycle->V_num ; i++) { hamcycle->count_node[i] = 0; }
    // store edges //
    hamcycle->edges.assign(2*hamcycle->E_num, Edge());
    for (int i = 0; i < hamcycle->E_num; i++) 
    { 
      fin >> u >> v; 
      hamcycle->edges[2*i].u = u;
      hamcycle->edges[2*i].v = v;
      hamcycle->edges[2*i+1].u = v;
      hamcycle->edges[2*i+1].v = u;
      // original edge
      s[u].push_back(2*i);
      e[v].push_back(2*i);
      // inverse edge
      s[v].push_back(2*i+1);
      e[u].push_back(2*i+1);
      // count node
      hamcycle->count_node[u] = hamcycle->count_node[u] + 1;
      hamcycle->count_node[v] = hamcycle->count_node[v] + 1;
    }
    // debug
    cout << "==============================================" << endl;
    cout << "Node -> Edges start from it -> Edges end to it" << endl;
    cout << "==============================================" << endl;
    for (int i = 0 ; i < hamcycle->V_num ; i++)
    {
      cout << "node " << i << " -> s = < ";
      for (int j = 0 ; j < s[i].size() ; j++) { cout << s[i][j] << " "; }
      cout << "> -> e = < ";
      for (int j = 0 ; j < e[i].size() ; j++) { cout << e[i][j] << " "; }
      cout << ">" << endl; 
    }
    cout << endl;
    fin.close();
  // New gates with numGates
    hamcycle->initCircuit();
  // Initialize satsolver
    SatSolver solver;
    solver.initialize();
  // Transform input graph into CNF (HC_to_CNF())
    hamcycle->HC_to_CNF(solver, s, e);
  // Print whether satisfied and output the corresponding file (Print_HC_Result())
    bool result;
    solver.assumeProperty(hamcycle->gates[hamcycle->gate_use_size]->getVar() , true); // true: 結果為1叫做SAT
    result = solver.assumpSolve();
    hamcycle->report_result(solver, result); // 最後必須是1
  // Calculate end time
    end = clock();
    double time = end - start;
    cout << "====" << endl;
    cout << "Time" << endl;
    cout << "====" << endl;
    cout << "time = " << time/CLOCKS_PER_SEC << endl;
    cout << endl;
    return 0;
}
