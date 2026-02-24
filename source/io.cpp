#include "io.h"
#include "fem.h"
#include <fstream>
#include <iomanip>

void writeNodalDisplacements(const vector<double>& U, int nNodes, ofstream& out)
{

    out << "=========================================\n";
    out << "        NODAL DISPLACEMENTS [cm]\n";
    out << "=========================================\n\n";

    out << setw(10) << "Node"
        << setw(18) << "ux"
        << setw(18) << "uy"
        << setw(18) << "rotz\n";

    out << "--------------------------------------------------------------\n";

    for(int i=0;i<nNodes;i++)
    {
        out << setw(10) << i+1
            << setw(18) << U[3*i]
            << setw(18) << U[3*i+1]
            << setw(18) << U[3*i+2]
            << "\n";
    }

    out << "\n";
}

void writeReactions(const vector<double>& R,
                    const vector<int>& fixedDOF,
                    int nNodes,
                    ofstream& out)
{
    out << "=========================================\n";
    out << "        SUPPORT REACTIONS [kg & m]\n";
    out << "=========================================\n\n";

    out << setw(10) << "Node"
        << setw(18) << "Rx"
        << setw(18) << "Ry"
        << setw(18) << "Mz\n";

    out << "--------------------------------------------------------------\n";

    vector<bool> isFixed(3*nNodes,false);
    for(int d : fixedDOF)
        isFixed[d] = true;

    for(int i=0;i<nNodes;i++)
    {
        if(isFixed[3*i] || isFixed[3*i+1] || isFixed[3*i+2])
        {
            out << setw(10) << i+1
                << setw(18) << R[3*i]
                << setw(18) << R[3*i+1]
                << setw(18) << R[3*i+2]
                << "\n";
        }
    }

    out << "\n";
}

void writeElementForces(
    const vector<Element>& elements,
    const vector<double>& U,
    ofstream& out)
{
    out << "ELEMENT FORCES (LOCAL AXIS) [kg & m]\n";
    out << "--------------------------------------------------------------\n";
    out << setw(10) << "Element"
        << setw(10) << "Node"
        << setw(15) << "N"
        << setw(15) << "V"
        << setw(15) << "M\n";
    out << "--------------------------------------------------------------\n";

    for(size_t i=0;i<elements.size();i++)
    {
        double fl[6];
        elementInternalForces(elements[i],U,fl);

        // Node 1 end
        out << setw(10) << i+1
            << setw(10) << elements[i].n1
            << setw(15) << fl[0]
            << setw(15) << fl[1]
            << setw(15) << fl[2]
            << "\n";

        // Node 2 end
        out << setw(10) << i+1
            << setw(10) << elements[i].n2
            << setw(15) << fl[3]
            << setw(15) << fl[4]
            << setw(15) << fl[5]
            << "\n";
    }
}
