#ifndef FEM_H
#define FEM_H

#include <vector>
#include <string>

using namespace std;

struct Node {
    int id;
    double x;
    double y;
};

struct Element {
    int id, n1, n2;
    double E, A, I,q;
    string load_type;

    // computed
    double L;
    double c;
    double s;
};

// Geometry
void computeGeometry(Element &e, const vector<Node> &nodes);

// Element routines
void localStiffness(const Element &e, double k[6][6]);
void transformationMatrix(const Element &e, double T[6][6]);
void globalStiffness(const Element &e, double kg[6][6]);
void equivalentLoad(const Element &e, double fe[6]);

// Assembly
void assembleGlobal(vector<vector<double>>& K,double kg[6][6],int n1,int n2);
void assembleLoad(vector<double>& F,double fe[6],int n1,int n2);

// DOF handling
void partitionDOF(int totalDOF,const vector<int>& constrained,
                  vector<int>& freeDOF,vector<int>& fixedDOF);

void buildReducedSystem(const vector<vector<double>>& K,
                        const vector<double>& F,
                        const vector<int>& freeDOF,
                        vector<vector<double>>& Kff,
                        vector<double>& Ff);

void expandDisplacements(int totalDOF,
                         const vector<int>& freeDOF,
                         const vector<double>& Uf,
                         vector<double>& U);

void computeReactions(const vector<vector<double>>& K,
                      const vector<double>& F,
                      const vector<double>& U,
                      const vector<int>& fixedDOF,
                      vector<double>& R);

// Linear solver
void solveSystem(const vector<vector<double>>& K,
                 const vector<double>& F,
                 vector<double>& U);

void elementInternalForces(
    const Element& e,
    const vector<double>& U,
    double fl[6]);

#endif

