//
// Created by bacon on 6/25/18.
//

#include <gurobi_c++.h>
using namespace std;
/* Copyright 2018, Gurobi Optimization, LLC */

/* This example formulates and solves the following simple MIP model:

     maximize    x +   y + 2 z
     subject to  x + 2 y + 3 z <= 4
                 x +   y       >= 1
     x, y, z binary
*/

int main()
{
    try {

        int num_of_curves = 100;
        int num_of_nodes = num_of_curves + 1; // Add V0 as Source
        int num_of_edges = num_of_nodes * num_of_nodes;  // We have to consider the common source V0 when calculate all edges

        // initiate all vars
        vector<GRBVar> variables;
        GRBEnv env2 = GRBEnv();
        GRBModel mTSPmodel = GRBModel(env2);

        variables.resize(num_of_edges + num_of_curves + 1); // besides all edge pairs, we have take u_i and k into consideration except u_0=0
        for (int i = 0; i < num_of_edges; i++)
        {
            variables[i] = mTSPmodel.addVar(0.0, 1.0, 0.0, GRB_BINARY);     // x_ij is binary
        }
        for (int i = num_of_edges; i < num_of_edges + num_of_curves + 1; i++)
        {
            variables[i] = mTSPmodel.addVar(0.0, num_of_curves, 0.0, GRB_INTEGER);   // u_i and k are integer
        }

        // Set objectives
        vector<double> cost;
        vector<GRBVar> vars;
        cost.resize(num_of_edges);
        for (int i = 0; i < num_of_nodes; i++)    // Set w_ij
        {
            for (int j = 0; j < num_of_nodes; j++)
            {
                if (i == j or i == 0 or j == 0)
                {
                    cost[i*(num_of_nodes) + j] = 0;
                }
                else
                {
                    cost[i*(num_of_nodes) + j] = 1;
                }
            }
        }
        cost.emplace_back(0.1);     // Set Ita
        vars.insert(vars.begin(), variables.begin(), variables.begin()+num_of_edges);       // add x_ij
        vars.emplace_back(variables.back());        // add k

        GRBLinExpr mTSP_obj;

        mTSP_obj.addTerms(&cost[0], &vars[0], num_of_edges + 1);

        mTSPmodel.setObjective(mTSP_obj, GRB_MINIMIZE);

        // Add constraint
        vector<double> coeffs;
        vector<GRBVar> terms1;
        vector<GRBVar> terms2;

        // Constraint 1: exactly one of the incoming and outgoing edges of a node needs to be selected in the solution
        terms1.resize(num_of_curves);
        terms2.resize(num_of_curves);
        coeffs.resize(num_of_curves);
        fill(coeffs.begin(), coeffs.end(), 1);
        for (int i = 1; i < num_of_nodes; i++)
        {
            for (int j = 0; j < num_of_nodes; j++)
            {
                if (i == j)
                {
                    continue;
                }
                int idx = (j<i ? j : j-1);
                terms1[idx] = variables[j*num_of_nodes + i];
                terms2[idx] = variables[i*num_of_nodes + j];
            }
            GRBLinExpr constraint1, constraint2;
            constraint1.addTerms(&coeffs[0], &terms1[0], num_of_curves);
            mTSPmodel.addConstr(constraint1, GRB_EQUAL, 1.0);
            constraint2.addTerms(&coeffs[0], &terms2[0], num_of_curves);
            mTSPmodel.addConstr(constraint2, GRB_EQUAL, 1.0);
        }

        // Constraint 2: each of the k paths is required to start and end at the start node V0
        {
            for (int i = 1; i < num_of_nodes; i++)
            {
                terms1[i-1] = variables[i];
                terms2[i-1] = variables[i*num_of_nodes];
            }
            coeffs.emplace_back(-1);
            terms1.emplace_back(variables.back());
            terms2.emplace_back(variables.back());
            GRBLinExpr constraint1, constraint2;
            constraint1.addTerms(&coeffs[0], &terms1[0], num_of_curves + 1);
            constraint2.addTerms(&coeffs[0], &terms2[0], num_of_curves + 1);
            mTSPmodel.addConstr(constraint1, GRB_EQUAL, 0.0);
            mTSPmodel.addConstr(constraint2, GRB_EQUAL, 0.0);
        }

        // Constraint 3: Subtour elimination constraints
        // Important: remember that there is no u_0 as u_0 is set as 0
        // -k*x_0i + (b-1)x_0i -x_i0 + u_i + k <= b
        for (int i = 1; i < num_of_nodes; i++)
        {
            GRBQuadExpr constraint;
            constraint += - variables.back() * variables[i] + (num_of_curves-1) * variables[i]
                          - variables[i*num_of_nodes] + variables[num_of_edges+i-1] + variables.back();
            mTSPmodel.addQConstr(constraint, GRB_LESS_EQUAL, num_of_curves);
        }

        // Constraint 4: Subtour elimination constraints
        // u_i + x_0i >= 2
        for (int i = 1; i < num_of_nodes; i++)
        {
            GRBLinExpr constraint;
            constraint += variables[i] + variables[num_of_edges+i-1];
            mTSPmodel.addConstr(constraint, GRB_GREATER_EQUAL, 2);
        }

        // Constraint 5: Subtour elimination constraints
        // u_i - u_j + (b-k-1)x_ij + (b-1-k)x_ji + k <= b
        for (int i = 1; i < num_of_nodes; i++)
        {
            for (int j = 1; j < num_of_nodes; j++)
            {
                if (i==j)
                {
                    continue;
                }
                GRBQuadExpr constraint;
                constraint += variables[num_of_edges+i-1] - variables[num_of_edges+j-1]
                              + (num_of_curves + 1 - variables.back()) * variables[i*num_of_nodes+j]
                              + (num_of_curves - 1 - variables.back()) * variables[j*num_of_nodes+i]
                              + variables.back();
                mTSPmodel.addQConstr(constraint, GRB_LESS_EQUAL, num_of_curves);
            }
        }

        // Optimize
        mTSPmodel.optimize();

    } catch(GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch(...) {
        cout << "Exception during optimization" << endl;
    }

    return 0;
}