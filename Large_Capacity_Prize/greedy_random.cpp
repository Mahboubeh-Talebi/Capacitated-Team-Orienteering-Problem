#include "greedy_random.h"
#define INF 0x3f3f3f3f


/*
Here we implement the randomized algorithm by Arora and Scherer, "Randomized Algorithm for Informative Path Planning with Budget Constraints"

In this algorithm, we first select a subset of vertices S \subseteq V uniformly randomly, then if the tree spanning S costs at most the budget we consider S in the solution bestTree.
The the following process is repeated 3*|V|. A vertex v_new is selected randomly, if v is present in bestTree, it is removed from bestTree. Otherwise, it is added to bestTree if 
the minimum spanning tree covering bestTree and v_new costs at most the budget

*/
GreedyRandomAlg::GreedyRandomAlg(double capacity, double budget, std::vector<std::pair<double, double>> verticesCoordinations, std::vector<double> verticesPrizes, std::vector<double> verticesSizes, std::vector<double> verticesServicesTime){
     
     //Compute the edge cost based on vertices coordinations
      std::vector<std::vector<double>> edgesWeights = coordinationsToEdges(verticesCoordinations);
//      std::vector<int> gVertices;
      //Build the graph by verties prizes and edges weights
//      Graph G = buildGraph(verticesPrizes, edgesWeights, gVertices);
      
      
     double sumWeights = 0;
     
     //Uniformly sampling a set S
     srand(time(0));
     std::vector<int> S;
     double weightOfTree = 0;
     double prizeOfTree = 0;
     double sizeOfTree = 0;
//     std::vector<int> verticesOfSRespectingCapacity;
//     for(int i = 0; i < gVertices.size(); i++){
//         double ranNumber= ((double)rand()) / RAND_MAX;
////                std::cout<<" ranNumber " <<ranNumber<<" ";
//         if(ranNumber >= 0.5) {
//            S.push_back(i);
//            prizeOfTree += G.getVertex(i)->getPrize();
//         }
//     }
//     
//     dpForSize(S, capacity, edgesWeights, verticesPrizes, verticesSizes, sizeOfTree, prizeOfTree, weightOfTree, verticesOfSRespectingCapacity);
     S.push_back(0);
     //Save the shorTree if it costs at most the budegt 
     std::vector<int> bestTree;
     double bestPrize = -1;
     double bestSize = 0;
     //int bestCost = 0;
     if(weightOfTree <= 0.5*budget) {
        bestPrize = prizeOfTree;
        bestTree = S;//SverticesOfSRespectingCapacity;  
        sumWeights = weightOfTree;
        verticesGreedyRandom = bestTree;
     }
     
     std::random_device rdw; // obtain a random number from hardware
     std::mt19937 genVertex(rdw()); // seed the generator
     std::uniform_int_distribution<> distrV(1, verticesPrizes.size() - 1); // define the range
     
////     std::cout<<"\nList of vertices: ";
////     for(int v : gVertices) std::cout<< v <<" ";
     std::vector<ratioPrizeCost> vecRatioPrizeCost(verticesPrizes.size());
    for (int u = 1; u < verticesPrizes.size(); u++){
       double tempDivisor = verticesSizes[u] != 0 ? verticesSizes[u] : 1/INT_MAX;
       struct ratioPrizeCost a = {u, verticesPrizes[u]/*tempDivisor*/};
       vecRatioPrizeCost.push_back(a);
       std::sort(vecRatioPrizeCost.begin(), vecRatioPrizeCost.end(), [](const auto& i, const auto& j) { return i.ratioPrizeCost > j.ratioPrizeCost; } );
    }
     for(int i = 0; i <= 3 * verticesPrizes.size(); i++) {
        //Select a vertex randomly
        int newVertex = distrV(genVertex);
        //check if newVertex is in S
        std::vector<int>::iterator position = std::find(S.begin(), S.end(), newVertex);
        if (position != S.end()) { // == S.end() means the element was not found
            S.erase(position);
        }
        else{
           S.push_back(newVertex);
        }
        weightOfTree = 0;
        prizeOfTree = 0;
        sizeOfTree = 0;
        for ( int vert : S) {
            prizeOfTree += verticesPrizes[vert];
            sizeOfTree += verticesSizes[vert];
        }
//        weightOfTree = primeMST(S, edgesWeights);
//        //Find the minimum spanning tree on the graph induced by tempSolU by prime's alg
//        weightOfTree = primeMST(S, edgesWeights);
        if (bestPrize < prizeOfTree) {
           std::vector<int> verticesRespectingCapacity;
           //Find a subset of vertices in tempSolU with total size at most capacity and maximum prize
           if (sizeOfTree > capacity) 
             dpForSize(S, capacity, verticesCoordinations, verticesPrizes, verticesSizes, verticesServicesTime, sizeOfTree, prizeOfTree, weightOfTree, verticesRespectingCapacity);
           else {
                  weightOfTree = findTSP(S, verticesCoordinations, verticesServicesTime);
                  verticesRespectingCapacity = S;
           }
           
//           if (sizeOfTree < capacity && weightOfTree < budget)
//              addMoreVertices(vecRatioPrizeCost, verticesRespectingCapacity, sizeOfTree, prizeOfTree, weightOfTree, budget, capacity, S, verticesCoordinations, verticesPrizes, verticesSizes, verticesServicesTime, "GreeyRandom");
           if (weightOfTree <= budget && bestPrize < prizeOfTree) {
               bestPrize = prizeOfTree;
               bestSize = sizeOfTree;
               bestTree = verticesRespectingCapacity;  
               sumWeights = weightOfTree;
               verticesGreedyRandom = bestTree;
            }
        }
     }//End_For_i
    //Set the paprameter
    greedyRandomPrize = bestPrize;
////           std::cout<<"\n";
//    std::cout<<"\n-GreedyRandom's prize is: "<<greedyRandomPrize;
    greedyRandomCost = sumWeights;
//    std::cout<<"\n-GreedyRandom's cost is: "<<greedyRandomCost;
    
}


//Return the cost from the vertex v to the root vetex u
void GreedyRandomAlg::addVertex(const Graph G, int v, std::vector<int> prev, std::vector<std::vector<double>> edgesWeights, std::vector<bool> visitedV, std::vector<bool>& tempVisitedV, double& addedWeights, double& addedPrize){
       std::vector<bool> newVisitedV(tempVisitedV.size(), false);
       int current = v;
       while(!visitedV[current]) 
          {
             newVisitedV[current] = true;
             //add_edge(prev[current], current, Temp_Out_Tree);
             addedWeights+=edgesWeights[prev[current]][current];
             addedPrize += G.getVertex(current)->getPrize();
             //std::cout<<"prev[current] is "<<prev[current]<<" ";
             current = prev[current];
          }
          //std::cout<<"\n";
        tempVisitedV = newVisitedV;
}

//Return the obtained shortest paths from the root vertex s to the goal vertices by prev
//void GreedyRandomAlg::addPaths(int s, const Graph G, std::vector<int> prev, std::vector<int> goalVertices, std::vector<std::vector<double>> edgesWeights, std::vector<int>& shortTree, std::vector<bool>& visitedV, double& weightOfTree, double& prizeOfTree){
//     //std::cout<<"Vertices on the path from "<<s<<" to goal verrtices:\n";
////     std::vector<bool> newVisitedV(prev.size(), false);
//     for(int i=0; i<goalVertices.size(); ++i){
//         int current = goalVertices[i];
//         while(!visitedV[current] && current!=s) 
//            {  
//               //std::cout<<current<<" ";
//               visitedV[current]=true;//the current vetex is visited
//               //add_edge(prev[current], current, Temp_Out_Tree);
//               weightOfTree+=edgesWeights[prev[current]][current];
//               prizeOfTree+=G.getVertex(current)->getPrize();
//               shortTree.push_back(current);
//               current = prev[current];
//            }
//            //std::cout<<"\n";
//     }
//}


double GreedyRandomAlg::getGreedyRandomPrize(){return greedyRandomPrize;}
double GreedyRandomAlg::getGreedyRandomCost(){return greedyRandomCost;}

//std::vector<std::vector<int>> GreedyRandomAlg::getShortestPathTreeGreedyRandom(){return prevUGreedyRandom;}
std::vector<int> GreedyRandomAlg::getVerticesGreedyRandom(){return verticesGreedyRandom;}
