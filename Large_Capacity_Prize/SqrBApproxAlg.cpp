#include "SqrBApproxAlg.h"
#define INF 0x3f3f3f3f

/*
This is a modifcation of the algorithm by D'Angel et al. budgeted directed Steiner problems with submodular prize.

In this algorithm, for each vertex u, we select a subset of vertices SU such that for given parameters x and y, u is in SU, |SU| \le y*B^{1-x} and the distance from u to any vertex v \in SU is at most B^{x}.
Then we select a vertex z, such that S_z maximizes the prize, a minimum spanning tree spanning _z costs and most the budget

*/

SqrBApproxAlg::SqrBApproxAlg(double capacity, double budget, std::vector<std::pair<double, double>> verticesCoordinations, std::vector<double> verticesPrizes, std::vector<double> verticesSizes, std::vector<double> verticesServicesTime, int numThresholds, double subsetSizeFactor, int numThreads) {

        //Compute the edge cost based on vertices coordinations
        std::vector<std::vector<double>> edgesWeights = coordinationsToEdges(verticesCoordinations);
//        std::vector<int> gVertices;
//        for (int i = 0; i < verticesSizes.size(); i++) 
//            gVertices.push_back(i);
        //Build the graph by verties prizes and edges weights
//        Graph G = buildGraph(verticesPrizes, edgesWeights, gVertices);
//         //Find the shortest path from each vertex u to all vertices in g
//         std::vector<std::vector<int>> prev_u(gVertices.size());//prev<> represents the previous vertex for each vertex on the shortest path 
//         std::vector<std::vector<double>> dist_u(gVertices.size());//dist<> represents the distance of each vertex from the root vertex
//         std::cout<<"\n Before all_to_all_shortest_paths \n";
//         all_to_all_shortest_paths(G, gVertices, verticesServicesTime, prev_u, dist_u);         
        double temp_max_prize = -1;
        #pragma omp parallel for default(none)  shared(capacity, budget, sqrBCost, sqrBSize, sqrBPrize, verticesSizes, verticesPrizes, verticesServicesTime, edgesWeights, verticesCoordinations, subsetSizeFactor, numThresholds, std::cout) /*private()*/ num_threads(numThreads)
        for (int u = 0; u < verticesPrizes.size(); u++) {//Guess which vertex u is in the soloution
             std::vector<std::vector<ratioPrizeCost>> vecRatioPrizeCost(2);
             for (int j = 1; j < verticesPrizes.size(); j++) {
                if (j != u) {
                   double  tempDistSize = verticesSizes[j]*edgesWeights[u][j]*verticesServicesTime[j];
                   double tempDivisor = tempDistSize != 0 ? tempDistSize : 1/INT_MAX;
                   struct ratioPrizeCost a = {j, verticesPrizes[j]/tempDistSize};
                   vecRatioPrizeCost[0].push_back(a);
                   a = {j, verticesPrizes[j]};
                   vecRatioPrizeCost[1].push_back(a);
                }
             }
             std::sort(vecRatioPrizeCost[0].begin(), vecRatioPrizeCost[0].end(), [](const auto& i, const auto& j) { return i.ratioPrizeCost > j.ratioPrizeCost; } );
             std::sort(vecRatioPrizeCost[1].begin(), vecRatioPrizeCost[1].end(), [](const auto& i, const auto& j) { return i.ratioPrizeCost > j.ratioPrizeCost; } );
         
             std::vector<int> SU;
             double prizeU = -1;
             double costU = 0;
             double sizeU = 0;
             for (int i = 0; i < vecRatioPrizeCost.size(); i++) {
                  std::vector<int> SI;
                  double prizeI = -1;
                  double costI = 0;
                  double sizeI = 0;
                  for (int x = 0; x <= numThresholds; ++x) {
                     for (int y = subsetSizeFactor; y >= 1; y--) {//We set the size vertices we can pick by y variable
                      std::vector<int> vectorY;            
                      std::vector<ratioPrizeCost> vecRPCTemp;//Saving those vertices not in vectorY
                      vectorY.push_back(u);//Add vertex u
                      vectorY.push_back(0);//Add the depot
                      double prizeY = verticesPrizes[u]; //G.getVertex(u)->getPrize();
                      double costY = 0;//verticesServicesTime[u];
                      double sizeY = verticesSizes[u];
                      for (auto v : vecRatioPrizeCost[i]) {
                          //Save at most pow(budget, 1.0-(x/numThresholds))*y*0.5 vertices with larger prizes that are at distance no more than pow(budget, x/numThresholds) from u
                          double costV = edgesWeights[u][v.id] + verticesServicesTime[v.id] + verticesServicesTime[u];
                          if (costV <= std::pow(budget, x/numThresholds) && vectorY.size() <= std::pow(budget, 1.0-(x/numThresholds))*y) {
                             vectorY.push_back(v.id);
                             prizeY += verticesPrizes[v.id];
                             sizeY += verticesSizes[v.id];
                           }
                           else vecRPCTemp.push_back(v);
                      }//End_for_v
                      //Set SU
                      if(prizeI < prizeY) {
                         
                         //First check if the total size of vertices in vectorY is more than capacity, then run dynamic programming
                         std::vector<int> verticesRespectingCapacity;
                         if (sizeY > capacity)
                            dpForSize(vectorY, capacity, verticesCoordinations, verticesPrizes, verticesSizes, verticesServicesTime, sizeY, prizeY, costY, verticesRespectingCapacity);
                         else 
                            {
                               costY = findTSP(vectorY, verticesCoordinations, verticesServicesTime);
                               verticesRespectingCapacity = vectorY;
                            }
                            
                         if (sizeY < capacity && costY < budget)
                             addMoreVertices(vecRPCTemp, verticesRespectingCapacity, sizeY, prizeY, costY, budget, capacity, vectorY, verticesCoordinations, verticesPrizes, verticesSizes, verticesServicesTime, "SqrB");
     
                         if (costY <= budget && prizeI < prizeY){
                             prizeI = prizeY;
                             costI = costY;
                             SI = verticesRespectingCapacity;
                             sizeI = sizeY;
                             break;
                         }
                      }
                    }//End_For_y
                    
                   }//End_for_x
                   if (prizeU < prizeI){
                       prizeU = prizeI;
                       costU = costI;
                       SU = SI;
                       sizeU = sizeI;
                   }
             }//End_For_i
             if (sqrBPrize < prizeU) {
                  sqrBCost = costU;
                  sqrBPrize = prizeU;
                  sqrBSize = sizeU;
                  verticesSqrB = SU;
            }
        }//End_for_u
        //If the budget is not enought to select tree, the vertex with the highest prize is selected
        if(sqrBPrize == -1) {
          //std::cout<<"Vertex "<<verticesPrizes[0].id <<" with the highest prize is selected as output\n";
          sqrBPrize = verticesPrizes[0];
          verticesSqrB.push_back(0);
          sqrBCost = verticesServicesTime[0];
        }
}

//Return the obtained shortest paths from the root vertex s to the goal vertices by prev
//void SqrBApproxAlg::addPaths(int s, const Graph &G, std::vector<int> prev, std::vector<int> goalVertices, std::vector<std::vector<double>> edgesWeights, std::vector<int>& pathVertices, double& weightOfTree, double& prizeOfTree){
//     std::vector<bool> visitV(prev.size(), false);
//     //std::cout<<"Vertices on the path from "<<s<<" to goal verrtices:\n";
//     for(int i=0; i < goalVertices.size(); i++){
//         int current = goalVertices[i];
//         while(!visitV[current] && current!=s) 
//            {  
//               //std::cout<<current<<" ";
//               visitV[current] = true;//the current vetex is visited
//               //add_edge(prev[current], current, Temp_Out_Tree);
//               weightOfTree += edgesWeights[prev[current]][current];
//               prizeOfTree += G.getVertex(current)->getPrize();
//               pathVertices.push_back(current);
//               current = prev[current];
//            }
//            //std::cout<<"\n";
//     }
//}

//double SqrBApproxAlg::getSqrBPrize(){return sqrBPrize;}
//double SqrBApproxAlg::getSqrBCost(){return sqrBCost;}
//double SqrBApproxAlg::getSqrBSize(){return sqrBSize;}
std::vector<int> SqrBApproxAlg::getVerticesSqrB(){return verticesSqrB;}
//std::vector<std::vector<int>> SqrBApproxAlg::getShortestPathTree(){return prevU;}
