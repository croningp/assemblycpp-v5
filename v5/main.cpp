#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <stack>
#include <sys/stat.h>
#include <iomanip>
#include <bitset>
#include <csignal>
#include <vector>
#ifdef _WIN32
    #include <windows.h>
#else
    #include <unistd.h>
#endif

using namespace std;
typedef vector<int> vi;
typedef vector<bool> vb;
typedef pair<int, int> pii;
#define BITSET_LENGTH 512
#define MAX_INT 2147483647
#define HASH_DEPTH_MAX 7
using standardBitset = bitset<BITSET_LENGTH>;

#include "globalPrimitives.h"
#include "ufds.h"
#include "molGraph.h"
#include "vf2.h"
#include "molfileParser.h"
#include "graphio.h"
#include "treeCanon.h"
#include "assemblyState.h"
#include "graphHashes.h"
#include "dagEnumeration.h"
#include "duplicateMatching.h"
#include "fragmentation.h"
#include "pathwayGenerator.h"
#include "improvedBnB.h"
#include "signalHandler.h"
#include "ioflag.h"
#include "help.h"

/**
 * @brief Function to write out intermediate MAs before the calculation has terminated
 * 
 * @param filename output filename
 */
void writeoutIntermediateMAs(string &filename)
{
    ofstream ofs(filename);
    int compensation = 1;
    if (disjointCompensation) compensation = disjointFragments;
    for (size_t i = 0; i < intermediateMAs.size(); i++)
    {
        ofs << intermediateMAs[i].first << ' ' << intermediateMAs[i].second  - disjointFragments + 1 << '\n';
    }
}

/**
 * @brief Function that takes a molfile fileName.mol and finds its assembly index
 *
 * This function will open the file, parse it, and then calls improvedBnb
 *
 * @param fileName The molfile in question
 */
void assemblyCalculator(string &fileName)
{
    cout << "opening " << fileName << '\n';
    string in = fileName + ".mol";
    string out = fileName + "Out";
    string intermediates = fileName + "IntermediateMAs";
    ifstream molfile(in.c_str());
    if (molfile.is_open())
    {
        ofstream outputFile(out.c_str());
        moleculeName = fileName + "Pathway";
        molGraph mol_graph;
        vector<double> coords;
        molfileParser(molfile, mol_graph);
        clock_t t1 = clock();
        outputFile << fileName << " has assembly index: ";
        improvedBnB(mol_graph, outputFile);
        clock_t t2 = clock();
        outputFile << "time to completion: " << t2 - t1 << '\n';
        if (writeIntermediateMAs) writeoutIntermediateMAs(intermediates);
    }
    else
    {
        ifstream graphFile(fileName.c_str());
        if (graphFile.is_open())
        {
            ofstream outputFile(out.c_str());
            molGraph mol_graph;
            graphio(graphFile, mol_graph);
            moleculeName = fileName + "Pathway";
            clock_t startTime = clock();
            outputFile << fileName << " has assembly index: ";
            improvedBnB(mol_graph, outputFile);
            outputFile << "time to completion: " << clock() - startTime << '\n';
            if (writeIntermediateMAs) writeoutIntermediateMAs(intermediates);
        }
        else cout << "No file found\n";
    }
}

int main(int argc, char** argv)
{
    #ifdef _WIN32
        if (SetConsoleCtrlHandler((PHANDLER_ROUTINE)CtrlHandler, TRUE)) {}
    #else
        signal(SIGINT, signalHandler);
    #endif
    if (argc > 1)
    {
        cout << "argv1: " << argv[1] << "\n";
        string s = argv[1];

        if (argc > 2) flagParser(argc, argv);

        if (s == "--help") help();
        else assemblyCalculator(s);
    }
    else cout << "no file selected\n";
}