#include "SqrBApproxAlg.h"
#define INF 0x3f3f3f3f

/*
This is a modifcation of the algorithm by D'Angel et al. budgeted directed Steiner problems with submodular prize.

In this algorithm, for each vertex u, we select a subset of vertices SU such that for given parameters x and y, u is in SU, |SU| \le y*B^{1-x} and the distance from u to any vertex v \in SU is at most B^{x}.
Then we select a vertex z, such that S_z maximizes the prize, a minimum spanning tree spanning _z costs and most the budget

*/

SqrBApproxAlg::SqrBApproxAlg(double capacity, double budget, std::vector<std::pair<double, double>> verticesCoordinations, std::vector<double> verticesPrizes, std::vector<double> verticesSizes, int numThresholds, double subsetSizeFactor, int numThreads) {
         double sumWeights=0;
         
         std::vector<std::vector<double>> edgesWeights;
         //Compute edges costs based on the coordinations of vertices
         edgesWeights = coordinationsToEdges(verticesCoordinations);
//         //Find the shortest path from each vertex u to all vertices in g
//         std::vector<std::vector<int>> prev_u(G_Vertices.size());//prev<> represents the previous vertex for each vertex on the shortest path 
//         std::vector<std::vector<double>> dist_u(G_Vertices.size());//dist<> represents the distance of each vertex from the root vertex
////         std::cout<<"\n Before all_to_all_shortest_paths \n";
//         all_to_all_shortest_paths(G, G_Vertices, prev_u, dist_u);
         
//         std::cout<<"\n Before vecRatioPrizeCost \n";
     
        double temp_max_prize=-1;
        #pragma omp parallel for default(none)  shared(capacity, budget, sumWeights, verticesSizes, verticesPrizes, edgesWeights, verticesCoordinations, subsetSizeFactor, temp_max_prize, numThresholds, std::cout) /*private()*/ num_threads(numThreads)
         
        for (int u = 0; u < verticesPrizes.size(); u++) {//Guess which vertex u is in the soloution
           std::vector<ratioPrizeCost> vecRatioPrizeCost;
           for(int v = 1; v < verticesPrizes.size(); v++) {
               if (v != u) {
                  double  tempDistSize = verticesSizes[v]*edgesWeights[u][v];
                  double tempDivisor = tempDistSize != 0 ? tempDistSize : 1/INT_MAX;
                  struct ratioPrizeCost a = {v, verticesPrizes[v]/tempDistSize};
                  vecRatioPrizeCost.push_back(a);
               }
            }
            std::sort(vecRatioPrizeCost.begin(), vecRatioPrizeCost.end(), [](const auto& i, const auto& j) { return i.ratioPrizeCost > j.ratioPrizeCost; } );
              std::vector<int> SU;
           double prizeU = -1;
           double costU = 0;
           for (int x = 0; x <= numThresholds; ++x) {
              for (int y = subsetSizeFactor; y >= 1; --y) {//Fist we set the size vertices we can pick by y variable
                 std::vector<int> vectorY;
                 vectorY.push_back(u);
                 vectorY.push_back(0);
                 double costY = 0;
                 double sizeY = verticesSizes[u];
                 double prizeY = verticesPrizes[u];
                 std::vector<ratioPrizeCost> vecRPCost;
                 for (auto v : vecRatioPrizeCost) {
                     //Save at most pow(budget, 1.0-(x/numThresholds))*y*0.5 vertices with larger prizes that are at distance no more than pow(budget, x/numThresholds) from u
                     if(edgesWeights[u][v.id] <= std::pow(budget, x/numThresholds) && vectorY.size() <= std::pow(budget, 1.0-(x/numThresholds))*y){
                        vectorY.push_back(v.id);
                        prizeY += verticesPrizes[v.id];
                        sizeY += verticesSizes[v.id];
                      }
                     else vecRPCost.push_back(v);
                 }//End_for_v
                 //check if u is not in vectorY
//                 if(!(find(vectorY.begin(), vectorY.end(), u)!=vectorY.end())){
//                     vectorY.push_back(u);
//                     prizeY += G.getVertex(u)->getPrize();
//                 }
                 //Set SU
                 if(prizeU < prizeY) {
                    std::vector<int> verticesRespectingCapacity;
                    if (sizeY > capacity)
                          dpForSize(vectorY, capacity, verticesCoordinations, verticesPrizes, verticesSizes, sizeY, prizeY, costY, verticesRespectingCapacity);
                    else {
                           costY = findTSPApprox(vectorY, verticesCoordinations);
                           verticesRespectingCapacity = vectorY;
                    }
                    if (costY < budget && sizeY < capacity) 
                         addMoreVertices(vecRPCost, verticesRespectingCapacity, sizeY, prizeY, costY, budget, capacity, vectorY, verticesCoordinations, verticesPrizes, verticesSizes, "SqrB");
                    if(costY <= budget && prizeU < prizeY){
                        prizeU = prizeY;
                        costU = costY;
                        SU = verticesRespectingCapacity;
                        break;
                    }
                 }
               }//End_For_y
               if(temp_max_prize < prizeU){
                   sumWeights = costU;
                   temp_max_prize = prizeU;
                   verticesSqrB = SU;
               }
              }//End_for_x
        }//End_For_u
        //If the budget is not enought to select tree, the vertex with the highest prize is selected
        if(temp_max_prize == -1) {
          //std::cout<<"Vertex "<<verticesPrizes[0].id <<" with the highest prize is selected as output\n";
          sqrB_prize = verticesPrizes[0];
          verticesSqrB.push_back(0);
          sumWeights = 0;
        }
        else sqrB_prize = temp_max_prize;
    sqrB_cost = sumWeights;
    //std::cout<<"\n-sqrB_alg's cost is: "<<sqrB_cost;
    
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

double SqrBApproxAlg::getSqrBPrize(){return sqrB_prize;}
double SqrBApproxAlg::getSqrBCost(){return sqrB_cost;}
std::vector<int> SqrBApproxAlg::getVerticesSqrB(){return verticesSqrB;}
//std::vector<std::vector<int>> SqrBApproxAlg::getShortestPathTree(){return prevU;}
