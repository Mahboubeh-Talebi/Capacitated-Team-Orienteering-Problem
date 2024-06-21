#ifndef GreedyRandomAlg_H_
#define GreedyRandomAlg_H_
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
class GreedyRandomAlg{
     private:
         double greedyRandomPrize;
         double greedyRandomCost;
         std::vector<int> verticesGreedyRandom;
         std::vector<std::vector<int>> prevUGreedyRandom;
     public:
         //find the tree started at s including goalVertices. The tree represents the shortest paths from s to the goal vertices. prev is obtained from the Dijkestra algorithm in which prev[u] shows the previous vertex on the shortest path from s to u
         //void Add_Paths(int s, std::vector<int> prev, std::vector<int> goalVertices, std::vector<std::vector<double>> EdgesWeights, std::vector<int>& PathVertices, double& sumWeights);
         
         //Find the path from v to the root vertex
         void addVertex(const Graph G, int v, std::vector<int> prev, std::vector<std::vector<double>> edgesWeights, std::vector<bool> visitedV, std::vector<bool>& tempVisitedV, double& addedWeights, double& addedPrize);
         void addPaths(int s, const Graph G, std::vector<int> prev, std::vector<int> goalVertices, std::vector<std::vector<double>> edgesWeights, std::vector<int>& shortTree, std::vector<bool>& visitedV, double& weightOfTree, double& prizeOfTree);
         
         //Constructor
         GreedyRandomAlg(double capacity, double Budget, std::vector<std::pair<double, double>> verticesCoordinations, std::vector<double> verticesPrizes, std::vector<double> verticesSizes, std::vector<double> verticesServicesTime);
         
         double getGreedyRandomPrize();
         double getGreedyRandomCost();
         
         std::vector<int> getVerticesGreedyRandom();
//         std::vector<std::vector<int>> getShortestPathTreeGreedyRandom();
 };
 
 #endif /* GreedyAlg_H_ */
