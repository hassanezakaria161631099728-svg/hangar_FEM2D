#include "fem.h"
#include "io.h"
#include <iostream>
#include <fstream>
int main()
{
vector<Node> nodes = {
    {1,0,0},
    {2,0,7},
    {3,9,9},
    {4,18,7},
    {5,18,0}
};

double E=2.1e8;
double A1=131.4e-4;
double I1=19270e-8;
double A2=84.5e-4;
double I2=23130e-8;
double q=-195;

vector<Element> elements;

elements.push_back({1,1,2,E,A1,I1,0,"local"});
elements.push_back({2,2,3,E,A2,I2,q,"global"});
elements.push_back({3,3,4,E,A2,I2,q,"global"});
elements.push_back({4,4,5,E,A1,I1,0,"local"});


    int totalDOF = nodes.size()*3;

    vector<vector<double>> K(totalDOF,
        vector<double>(totalDOF,0.0));

    vector<double> F(totalDOF,0.0);

    // Geometry + Assembly
    for(auto &e : elements){
        computeGeometry(e,nodes);

        double kg[6][6];
        globalStiffness(e,kg);

        assembleGlobal(K,kg,e.n1,e.n2);

        double fe[6];
        equivalentLoad(e,fe);
        assembleLoad(F,fe,e.n1,e.n2);
    }

    // Boundary conditions
    vector<int> bc = {1,2,3,13,14,15};
// nodal loads
F[4]+= -1387;
F[10]+= -1387;

    vector<int> freeDOF,fixedDOF;
    partitionDOF(totalDOF,bc,freeDOF,fixedDOF);

    vector<vector<double>> Kff;
    vector<double> Ff;

    buildReducedSystem(K,F,freeDOF,Kff,Ff);

    vector<double> Uf(freeDOF.size());
    solveSystem(Kff,Ff,Uf);

    vector<double> U;
    expandDisplacements(totalDOF,freeDOF,Uf,U);

    vector<double> R;
    computeReactions(K,F,U,fixedDOF,R);

   ofstream out("results.txt");

    writeNodalDisplacements(U,nodes.size(),out);
    writeReactions(R,fixedDOF,nodes.size(),out);
    writeElementForces(elements,U,out);

out.close();
    return 0;
}

