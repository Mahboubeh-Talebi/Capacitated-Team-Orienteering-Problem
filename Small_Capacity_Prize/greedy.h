#ifndef GreedyAlg_H_
#define GreedyAlg_H_
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
class GreedyAlg{
     private:
         double greedyPrize;
         double greedyCost;
         std::vector<int> verticesGreedy;
     public:
         //find the tree started at s including goalVertices. The tree represents the shortest paths from s to the goal vertices. prev is obtained from the Dijkestra algorithm in which prev[u] shows the previous vertex on the shortest path from s to u
         //void Add_Paths(int s, std::vector<int> prev, std::vector<int> goalVertices, std::vector<std::vector<double>> EdgesWeights, std::vector<int>& PathVertices, double& sumWeights);
         
         //Find the path from v to the root vertex
//         void addVertex(const Graph G, int v, std::vector<int> prev, std::vector<std::vector<double>> edgesWeights, std::vector<bool> visitedV, std::vector<bool>& tempVisitedV, double& addedWeights, double& addedPrize);

         //Constructor
         GreedyAlg(double capacity, double budget, std::vector<std::pair<double, double>> verticesCoordinationsOfAlg, std::vector<double> verticesPrizes, std::vector<double> verticesSizes, int numThreads);
         
         double getGreedyPrize();
         double getGreedyCost();
         
         std::vector<int> getVerticesGreedy();
 };
 
 #endif /* GreedyAlg_H_ */
