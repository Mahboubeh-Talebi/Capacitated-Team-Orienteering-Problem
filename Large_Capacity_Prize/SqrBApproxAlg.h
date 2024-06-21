#ifndef ALG_H_
#define ALG_H_
#include <chrono>
#include <iostream>
#include <fstream>
#include<string> 
#include <cassert>
#include <sstream>
#include <map> 
#include <queue>
#include <algorithm>
#include <chrono>
#include <stdio.h>
#include <cstdlib>
#include "NodeInfo.h"
class SqrBApproxAlg{
     private:
         double sqrBPrize = -1;
         double sqrBCost = 0;
         double sqrBSize = 0;
         std::vector<int> verticesSqrB;
     public:
         //find the tree started at s including goalVertices. The tree represents the shortest paths from s to the goal vertices. prev is obtained from the Dijkestra algorithm in which prev[u] shows the previous vertex on the shortest path from s to u
//         void addPaths(int s, const Graph &G, std::vector<int> prev, std::vector<int> goalVertices, std::vector<std::vector<double>> edgesWeights, std::vector<int>& pathVertices, double& weightOfTree, double& prizeOfTree);
         
         //Return the cost from the vertex current to vertex u
         std::vector<bool> addCost(int u, std::vector<int> prev, int current, std::vector<std::vector<double>> edgesWeights, std::vector<bool>& visitedV, double& addedWeights);
         
         
         //Constructor
         SqrBApproxAlg(double capacity, double budget, std::vector<std::pair<double, double>> verticesCoordinations, std::vector<double> verticesPrizes, std::vector<double> verticesSizes, std::vector<double> verticesServicesTimeOfAlg, int numThresholds, double subsetSizeFactor, int numThreads);
         
//         void addMoreVertices(std::vector<ratioPrizeCost> vecRatioPrizeCost, std::vector<int>& verticesRespectingCapacity, double& size_x, double& prize_x, double& cost_x, double budget, double capacity, std::vector<int> vectTemp, std::vector<std::vector<double>> edgesWeights, std::vector<double> verticesPrizes, std::vector<double> verticesSizes);
//         double getSqrBPrize();
//         double getSqrBCost();
//         double getSqrBSize();
         std::vector<int> getVerticesSqrB();
//         std::vector<std::vector<int>> getShortestPathTree();
 };
 
 #endif /* ALG_H_ */
