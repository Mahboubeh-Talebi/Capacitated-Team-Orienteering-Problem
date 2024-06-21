#include <algorithm>
#include <chrono>
#include <unistd.h>
#include <ios>
#include <iostream>
#include <fstream>
#include<string> 
#include <cassert>
#include <sstream>
#include <map>
#include <cmath>
#include <getopt.h>
#include <iomanip>
#include <random>
#include <filesystem>
//using namespace std::chrono;
//#include "Graph.cpp"
#include "Subset.cpp"
#include "Subroutine.cpp"
#include "PD.cpp"
#include "SqrBApproxAlg.cpp"
#include "graph_generator.cpp"
#include "greedy.cpp"
#include "greedy_random.cpp"


class OutputResults{

          int numVehicles;
          double budget;
          double capacity;
          int numThresholds;
          int numThreads;
          double subsetSizeFactor;
          double guessedUpperOpt;
          
          //Declare 2-aaprox parameters
          std::chrono::high_resolution_clock::time_point startOf2Approx;// = high_resolution_clock::now();
          std::chrono::high_resolution_clock::time_point stopOf2Approx;// = high_resolution_clock::now();
          std::vector<std::pair<double, double>> verticesCoordinationsOf2Approx;
          std::vector<double> verticesPrizesOf2Approx;
          std::vector<double> verticesSizesOf2Approx;
          
          //Declare sqrB-approx parameters
          std::chrono::high_resolution_clock::time_point startOfSqrBApprox;// = high_resolution_clock::now();
          std::chrono::high_resolution_clock::time_point stopOfSqrBApprox;// =  high_resolution_clock::now();
          std::vector<std::pair<double, double>> verticesCoordinationsOfSqrBApprox;
          std::vector<double> verticesPrizesOfSqrBApprox;
          std::vector<double> verticesSizesOfSqrBApprox;
          
          //Declare Greedy parameters
          std::chrono::high_resolution_clock::time_point startOfGreedy;// = high_resolution_clock::now();
          std::chrono::high_resolution_clock::time_point stopOfGreedy;// = high_resolution_clock::now();
          std::vector<std::pair<double, double>> verticesCoordinationsOfGreedy;
          std::vector<double> verticesPrizesOfGreedy;
          std::vector<double> verticesSizesOfGreedy;
          
          //Decclare GreedyRandom algorithms
          std::chrono::high_resolution_clock::time_point startOfGreedyRandom;// = high_resolution_clock::now();
          std::chrono::high_resolution_clock::time_point stopOfGreedyRandom;// = high_resolution_clock::now();
          std::vector<std::pair<double, double>> verticesCoordinationsOfGreedyRandom;
          std::vector<double> verticesPrizesOfGreedyRandom;
          std::vector<double> verticesSizesOfGreedyRandom;
          
          //Declare compute functions for algorithms 
          void compute2Approx(Graph G, std::vector<std::pair<double, double>> verticesCoordinationsOfAlg, std::vector<double> verticesPrizesOfAlg, std::vector<double> verticesSizesOfAlg, std::vector<int>& verticesOf2ApproxTemp, double& costOf2ApproxTemp, double& prizeOf2ApproxTemp, double& sizeOf2ApproxTemp);
          
          void computeSqrBApprox(std::vector<std::pair<double, double>> verticesCoordinationsOfAlg, std::vector<double> verticesPrizesOfAlg, std::vector<double> verticesSizesOfAlg, std::vector<int>& verticesOfAlg, double& costOfAlg, double& prizeOfAlg, double& sizeOfAlg);
          
          void computeGreedy(std::vector<std::pair<double, double>> verticesCoordinationsOfAlg, std::vector<double> verticesPrizesOfAlg, std::vector<double> verticesSizesOfAlg, std::vector<int>& verticesOfAlg, double& costOfAlg, double& prizeOfAlg, double& sizeOfAlg);
          
          void computeGreedyRandom(std::vector<std::pair<double, double>> verticesCoordinationsOfAlg, std::vector<double> verticesPrizesOfAlg, std::vector<double> verticesSizesOfAlg, std::vector<int>& verticesOfAlg, double& costOfAlg, double& prizeOfAlg, double& sizeOfAlg);
                        
          public:
               //constructor
               OutputResults(std::string graphT, int numVertices, int valNumVehicles, std::vector<std::pair<double, double>> valVerticesCoordinations, std::vector<double> valVerticesPrizes, std::vector<double> valVerticesSizes, double valBudget, double valCapacity, int valNumThresholds, int valNumThreads, int solve2Approx, double valSubsetSizeFactor, double factorViolation, std::string outputDirectory);
               
               
               void removeCoveredVertices(std::vector<int> coveredVertices, std::vector<std::pair<double, double>>& verticesCoordinationOfAlg, std::vector<double>& verticesPrizesOfAlg, std::vector<double>& verticesSizeOfAlg);
               
               
               void binarySearch(std::vector<std::pair<double, double>> verticesCoordinationsOfAlg, std::vector<double> verticesPrizesOfAlg, std::vector<double> verticesSizesOfAlg, std::string nameOfAlg, std::vector<int>& verticesOfAlg, double& costOfAlg, double& prizeOfAlg, double& sizeOfAlg);
               
               
               double sizeTree(std::list<std::shared_ptr<Edge>>& tree, std::vector<double> verticesSizesOfAlg, std::vector<int>& verticesTree);            
               
               void findPrevOfTree(std::list<std::shared_ptr<Edge> > & edges, std::vector<bool>& visited, std::vector<int>& prev, int v);
};

