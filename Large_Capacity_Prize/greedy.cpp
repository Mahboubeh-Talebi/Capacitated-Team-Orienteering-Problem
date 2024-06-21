#include "greedy.h"
#define INF 0x3f3f3f3f

/*
Here we implement the greedy algorithm in which for each vertex v, first we suppose v is in the solution and add other vertices to the solution
to maximize the total prize until the minimum spanning tree spanning the vertices deosn't exceed the budget constraint
*/
GreedyAlg::GreedyAlg(double capacity, double budget, std::vector<std::pair<double, double>> verticesCoordinations, std::vector<double> verticesPrizes, std::vector<double> verticesSizes, std::vector<double> verticesServicesTime, int numThreads){
     
     

      //Compute the edge cost based on vertices coordinations
      std::vector<std::vector<double>> edgesWeights = coordinationsToEdges(verticesCoordinations);
//      std::vector<int> gVertices;
//      //Build the graph by verties prizes and edges weights
////      Graph G = buildGraph(verticesPrizes, edgesWeights, gVertices);

      double sumWeights = 0;
       
     //Find the shortest path from each vertex u to all vertices in g
//     std::vector<std::vector<int>> prev_u(gVertices.size());//prev<> represents the previous vertex for each vertex on the shortest path 
//     std::vector<std::vector<double>> dist_u(gVertices.size());//dist<> represents the distance of each vertex from the root vertex
//     all_to_all_shortest_paths(G, gVertices, verticesServicesTime, prev_u, dist_u);
     //Set the paprameter
     
     double temp_max_prize=-1;
     #pragma omp parallel for default(none)  shared(budget, capacity, verticesPrizes, verticesSizes, verticesServicesTime, sumWeights, edgesWeights, verticesCoordinations, temp_max_prize) /*private()*/ num_threads(numThreads)
     for (int u = 0; u < verticesPrizes.size(); u++) {
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
         double weightOfTree = 0;//Save the weight of the tree in each vertex u
         double prizeOfTree = 0;//Save the prize of the tree in each vertex u
         std::vector<int> solU;
         for (int i = 0; i < vecRatioPrizeCost.size(); i++) {
            double weightOfTreeI = 0;//Save the weight of the tree in each vertex u
            double prizeOfTreeI = 0;//G.getVertex(u)->getPrize();//Save the prize of the tree in each vertex u
            double sizeOfTreeI = 0;
            std::vector<int> solUI;
            std::vector<int> tempSolU;
            tempSolU.push_back(u);//u should be in the solution
            //Check vertices one by one whether it is possible to add them to the solution
            for (auto v : vecRatioPrizeCost[i]) {                                            
               double tempWeightOfTree = 0;
               double tempPrizeOfTree = 0;
               double tempSizeOfTree = 0;
               //Add vertex v.id to the temp solution if it's not added yet
    //           if(!(find(tempSolU.begin(), tempSolU.end(), v.id)!=tempSolU.end()))
               tempSolU.push_back(v.id);
    //           for ( int vert : tempSolU) tempPrizeOfTree += verticesPrizes[vert];
               //check if the minimum spanning tree on the vertices of tempSolU is feasible and has a better prize
               std::vector<int> verticesRespectingCapacity;
               //Find a subset of vertices in tempSolU with total size at most capacity and maximum prize
               dpForSize(tempSolU, capacity, verticesCoordinations, verticesPrizes, verticesSizes, verticesServicesTime, tempSizeOfTree, tempPrizeOfTree, tempWeightOfTree, verticesRespectingCapacity);
               if (tempWeightOfTree <= budget && prizeOfTreeI < tempPrizeOfTree) {
                   prizeOfTreeI = tempPrizeOfTree;
                   weightOfTreeI = tempWeightOfTree;
                   solUI = verticesRespectingCapacity;
                   sizeOfTreeI = tempSizeOfTree;
               }
               
               tempSolU = solUI;
               
            }//End_for_v
            if(prizeOfTree < prizeOfTreeI){
                 weightOfTree = weightOfTreeI;
                 prizeOfTree = prizeOfTreeI;
                 solU = solUI;   
            }
        }//End_For_i
        if(temp_max_prize < prizeOfTree){
             sumWeights = weightOfTree;
             temp_max_prize = prizeOfTree;
             verticesGreedy = solU;   
        }
     }//End_For_u
    greedyPrize=temp_max_prize;
//    std::cout<<"Greedy's prize is: "<<greedyPrize<<"\n";
    greedyCost = sumWeights;
//    std::cout<<"\n-Greedy's cost is: "<<greedyCost;
    
}


//Return the cost from the vertex v to the root vetex u
//void GreedyAlg::addVertex(const Graph G, int v, std::vector<int> prev, std::vector<std::vector<double>> edgesWeights, std::vector<bool> visitedV, std::vector<bool>& tempVisitedV, double& addedWeights, double& addedPrize){
//       std::vector<bool> newVisitedV(tempVisitedV.size(), false);
//       int current = v;
//       while(!visitedV[current]) 
//          {
//             newVisitedV[current] = true;
//             //add_edge(prev[current], current, Temp_Out_Tree);
//             addedWeights+=edgesWeights[prev[current]][current];
//             addedPrize += G.getVertex(current)->getPrize();
//             //std::cout<<"prev[current] is "<<prev[current]<<" ";
//             current = prev[current];
//          }
//          //std::cout<<"\n";
//        tempVisitedV = newVisitedV;
//}


double GreedyAlg::getGreedyPrize(){return greedyPrize;}
double GreedyAlg::getGreedyCost(){return greedyCost;}


std::vector<int> GreedyAlg::getVerticesGreedy(){return verticesGreedy;}
