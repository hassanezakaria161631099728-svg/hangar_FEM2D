#ifndef IO_H
#define IO_H
#include <fstream>
#include <vector>
#include "fem.h"
using namespace std;

void writeNodalDisplacements(const vector<double>& U, int nNodes,ofstream& out);
void writeReactions(const vector<double>& R,const vector<int>& fixedDOF,int nNodes,ofstream& out);
void writeElementForces(const vector<Element>& elements,const vector<double>& U,ofstream& out);

#endif

