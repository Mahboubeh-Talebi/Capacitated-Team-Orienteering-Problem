#include "output.h"
#include <stdio.h>
#define INF 0x3f3f3f3f
OutputResults::OutputResults(std::string graphT, int numVertices, int valNumVehicles, std::vector<std::pair<double, double>> valVerticesCoordinations, std::vector<double> valVerticesPrizes, std::vector<double> valVerticesSizes, std::vector<double> valVerticesServicesTime, double valBudget, double valCapacity, int valNumThresholds, int valNumThreads, int solve2Approx, double valSubsetSizeFactor, double factorViolation, std::string outputDirectory){
    //Set the parameters
    numVehicles = valNumVehicles;
    budget = (1 + factorViolation/2)*valBudget;
    capacity = valCapacity;
    numThresholds = valNumThresholds;
    numThreads = valNumThreads;
    //Set the number of vertices and edges
    std::pair<int, int> numVerticesEdges = std::make_pair(numVertices, std::pow(numVertices, 2)/2);
    
    
    //Set 2-approx parameters
    verticesCoordinationsOf2Approx = valVerticesCoordinations;
    verticesPrizesOf2Approx = valVerticesPrizes;
    verticesSizesOf2Approx = valVerticesSizes;
    verticesServicesTimeOf2Approx = valVerticesServicesTime;
    std::vector<int> verticesOf2Approx;
    
    
    //Set sqrB-approx parameters
    subsetSizeFactor = valSubsetSizeFactor;
    verticesCoordinationsOfSqrBApprox = valVerticesCoordinations;
    verticesPrizesOfSqrBApprox = valVerticesPrizes;
    verticesSizesOfSqrBApprox = valVerticesSizes;
    verticesServicesTimeOfSqrBApprox = valVerticesServicesTime;
    std::vector<int> verticesOfSqrBApprox;
    
    //Set Greedy parameters
    verticesCoordinationsOfGreedy = valVerticesCoordinations;
    verticesPrizesOfGreedy = valVerticesPrizes;
    verticesSizesOfGreedy = valVerticesSizes;
    verticesServicesTimeOfGreedy = valVerticesServicesTime;
    std::vector<int> verticesOfGreedy;
    
    //Set GreedyRandom algorithms
    verticesCoordinationsOfGreedyRandom = valVerticesCoordinations;
    verticesPrizesOfGreedyRandom = valVerticesPrizes;
    verticesSizesOfGreedyRandom = valVerticesSizes;
    verticesServicesTimeOfGreedyRandom = valVerticesServicesTime;
    std::vector<int> verticesOfGreedyRandom;
     
     
    std::vector<double> prizeOf2Approx(numVehicles, 0);
    std::vector<double> prizeOfSqrBApprox(numVehicles, 0);
    std::vector<double> prizeOfGreedy(numVehicles, 0);
    std::vector<double> prizeOfGreedyRandom(numVehicles, 0);
    
    std::vector<double> costOf2Approx(numVehicles, 0);
    std::vector<double> costOfSqrBApprox(numVehicles, 0);
    std::vector<double> costOfGreedy(numVehicles, 0);
    std::vector<double> costOfGreedyRandom(numVehicles, 0);
    
    std::vector<double> sizeOf2Approx(numVehicles, 0);
    std::vector<double> sizeOfSqrBApprox(numVehicles, 0);
    std::vector<double> sizeOfGreedy(numVehicles, 0);
    std::vector<double> sizeOfGreedyRandom(numVehicles, 0);
    
    std::vector<double> runningTimeOf2Approx(numVehicles, 0);
    std::vector<double> runningTimeOfSqrBApprox(numVehicles, 0);
    std::vector<double> runningTimeOfGreedy(numVehicles, 0);
    std::vector<double> runningTimeOfGreedyRandom(numVehicles, 0);
    
    //Find paths for vehicles by 2-approx
    if (solve2Approx == 1)
      for (int k = 0; k < numVehicles; k++) {
          if (verticesCoordinationsOf2Approx.size() > 1) {

             std::cout << "\n-----------Starting the 2-approx for vehicle "<< k + 1 <<" ----------- \n\n";
             startOf2Approx = high_resolution_clock::now();
             //First remove the covered vertices by previuse vehicle
             removeCoveredVertices(verticesOf2Approx, verticesCoordinationsOf2Approx, verticesPrizesOf2Approx, verticesSizesOf2Approx, verticesServicesTimeOf2Approx);
             verticesOf2Approx.clear();
             binarySearch(verticesCoordinationsOf2Approx, verticesPrizesOf2Approx, verticesSizesOf2Approx, verticesServicesTimeOf2Approx, "2Approx", verticesOf2Approx, costOf2Approx[k], prizeOf2Approx[k], sizeOf2Approx[k]);
             stopOf2Approx = high_resolution_clock::now();
             runningTimeOf2Approx[k] = duration_cast<seconds>(stopOf2Approx - startOf2Approx).count();
             std::cout << "\n-The total prize of items in 2-approx collected by vehicle " << k + 1 << " is: " << prizeOf2Approx[k];
             std::cout << "\n-The total size of items in 2-approx collected by vehicle " << k + 1 << " is: " << sizeOf2Approx[k];
             std::cout << "\n-The total cost of path in 2-approx taken by vehicle " << k + 1 << " is: " << costOf2Approx[k];
             std::cout << "\n-The running time of 2-approx alg for computing a path for vehicle "<< k + 1 << " is: " << runningTimeOf2Approx[k] <<" seconds.";
             std::cout << "\n\n-Ending the 2-approx for vehicle "<< k + 1 <<" ----------- \n\n\n"; 
          }
     }
   
   //Find paths for vehicles by sqrB-approx
   for (int k = 0; k < numVehicles; k++) {
      if (verticesCoordinationsOfSqrBApprox.size() > 1) {
           
           std::cout << "\n-----------Starting the sqrB-approx for vehicle "<< k + 1 <<" ----------- \n\n";
           startOfSqrBApprox = high_resolution_clock::now();
           //First remove the covered vertices by previuse vehicle
           removeCoveredVertices(verticesOfSqrBApprox, verticesCoordinationsOfSqrBApprox, verticesPrizesOfSqrBApprox, verticesSizesOfSqrBApprox, verticesServicesTimeOfSqrBApprox);
           verticesOfSqrBApprox.clear();
           binarySearch(verticesCoordinationsOfSqrBApprox, verticesPrizesOfSqrBApprox, verticesSizesOfSqrBApprox, verticesServicesTimeOfSqrBApprox, "SqrBApprox", verticesOfSqrBApprox, costOfSqrBApprox[k], prizeOfSqrBApprox[k], sizeOfSqrBApprox[k]);
//           std::cout<<"Size of vectors " <<verticesCoordinationsOfSqrBApprox.size() <<" "<< verticesPrizesOfSqrBApprox.size()<<" "<< verticesSizesOfSqrBApprox.size()<<"\n";
           stopOfSqrBApprox =  high_resolution_clock::now();
           runningTimeOfSqrBApprox[k] = duration_cast<seconds>(stopOfSqrBApprox - startOfSqrBApprox).count();
           std::cout << "\n-The total prize of items in the sqrB-approx alg collected by vehicle " << k + 1 << " is: " << prizeOfSqrBApprox[k];
           std::cout << "\n-The total size of items in the sqrB-approx alg collected by vehicle " << k + 1  << " is: " << sizeOfSqrBApprox[k];
           std::cout << "\n-The total cost of path in the sqrB-approx alg taken by vehicle " << k + 1  << " is: " << costOfSqrBApprox[k];
           std::cout << "\n-The running time of the sqrB-approx alg for computing a path for vehicle "<< k + 1 << " is: " << runningTimeOfSqrBApprox[k]<<" seconds.";
           std::cout << "\n\n-----------Ending the sqrB-approx for vehicle "<< k + 1 <<" ----------- \n\n\n";
           
        }
   }
   
   //Find paths for vehicles by Greedy
   for (int k = 0; k < numVehicles; k++) {
         if (verticesCoordinationsOfGreedy.size() > 1) {
          
           std::cout << "\n-----------Starting the Greedy for vehicle "<< k + 1 <<" ----------- \n\n";
           startOfGreedy = high_resolution_clock::now();
           //First remove the covered vertices by previuse vehicle
           removeCoveredVertices(verticesOfGreedy, verticesCoordinationsOfGreedy, verticesPrizesOfGreedy, verticesSizesOfGreedy, verticesServicesTimeOfGreedy);
           verticesOfGreedy.clear();
           binarySearch(verticesCoordinationsOfGreedy, verticesPrizesOfGreedy, verticesSizesOfGreedy, verticesServicesTimeOfGreedy, "Greedy", verticesOfGreedy, costOfGreedy[k], prizeOfGreedy[k], sizeOfGreedy[k]);
           stopOfGreedy = high_resolution_clock::now();
           runningTimeOfGreedy[k] = duration_cast<seconds>(stopOfGreedy - startOfGreedy).count();
           std::cout << "\n-The total prize of items in the Greedy alg collected by vehicle " << k + 1 << " is: " << prizeOfGreedy[k];
           std::cout << "\n-The total size of items in the Greedy alg collected by vehicle " << k + 1 << " is: " << sizeOfGreedy[k];
           std::cout << "\n-The total cost of path in the Greedy alg taken by vehicle " << k + 1 << " is: " << costOfGreedy[k];
           std::cout << "\n-The running time of the Greedy alg for computing a path for vehicle "<< k + 1 << " is: " << runningTimeOfGreedy[k]<<" seconds.";
           std::cout << "\n\n-----------Ending the Greedy for vehicle "<< k + 1 <<" ----------- \n\n\n";
           
        }
   }
   
   //Find paths for vehicles by GreedyRandom
   for (int k = 0; k < numVehicles; k++){
    
       if (verticesCoordinationsOfGreedyRandom.size() > 1) {
        
           std::cout << "\n-----------Starting the GreedyRandom for vehicle "<< k + 1 <<" ----------- \n\n";
           startOfGreedyRandom = high_resolution_clock::now();
           //First remove the covered vertices by previous vehicle
           removeCoveredVertices(verticesOfGreedyRandom, verticesCoordinationsOfGreedyRandom, verticesPrizesOfGreedyRandom, verticesSizesOfGreedyRandom, verticesServicesTimeOfGreedyRandom);
           verticesOfGreedyRandom.clear();
           //BinarySearch for GeedyRandom
           binarySearch(verticesCoordinationsOfGreedyRandom, verticesPrizesOfGreedyRandom, verticesSizesOfGreedyRandom, verticesServicesTimeOfGreedyRandom, "GreedyRandom", verticesOfGreedyRandom, costOfGreedyRandom[k], prizeOfGreedyRandom[k], sizeOfGreedyRandom[k]);
           stopOfGreedyRandom = high_resolution_clock::now();
           runningTimeOfGreedyRandom[k] = duration_cast<seconds>(stopOfGreedyRandom - startOfGreedyRandom).count();
           std::cout << "\n-The total prize of items in the GreedyRandom alg collected by vehicle " << k + 1 << " is: " << prizeOfGreedyRandom[k];
           std::cout << "\n-The total size of items in the GreedyRandom alg collected by vehicle " << k + 1 << " is: " << sizeOfGreedyRandom[k];
           std::cout << "\n-The total cost of path in the GreedyRandom alg taken by vehicle " << k + 1 << " is: " << costOfGreedyRandom[k];
           std::cout << "\n-The running time of the GreedyRandom alg for computing a path for vehicle "<< k + 1 << " is: " << runningTimeOfGreedyRandom[k]<<" seconds.";
           std::cout << "\n\n-----------Ending the GreedyRandom for vehicle "<< k + 1 <<" ----------- \n\n\n";
        }
     }//End-Of-For-numVehicles

     std::ofstream ofs;
     //Open the file in the append mode, meaning that in case the file exists we add data
     ofs.open(outputDirectory, std::ios::app);
     if (std::filesystem::is_empty(outputDirectory)) {
         ofs<<"g"<<std::setw(25)<<"n_nodes"<<std::setw(25)<< "numVehicles"<<std::setw(25)<<"capacity"<<std::setw(25)<<"budget"<<std::setw(25)<<"alg"<<std::setw(25)<<"score"<<std::setw(25)<<"cost"<<std::setw(25)<<"size"<<std::setw(25)<<"time(s)"<<std::setw(25)<<"n_thre"<<std::setw(25)<<"factorViolation"<<std::setw(25)<<"subsetSizeFactor"<<std::endl;
     }
     std::vector<std::string> vecAlgs = {"2Approx", "up_bound", "GreedyRandom", "SqrBApprox", "Greedy"};
     for(int i = 0; i < vecAlgs.size(); ++i){
          double tempRunningTime = 0; 
          double tempCost = 0;
          double tempPrize  =  0;
          double tempSize  =  0;
          std::string tempAlgName;
          if(vecAlgs[i] == "2Approx"){
            if(solve2Approx==1){
                tempAlgName = "2-approx";
                for(double c : costOf2Approx) tempCost += c;
                tempCost = 2*tempCost;
                for(double p : prizeOf2Approx) tempPrize += p;
                for(double s : sizeOf2Approx) tempSize += s;
                for(double r: runningTimeOf2Approx) tempRunningTime += r;
            }
            else continue;
          }
         else if(vecAlgs[i] == "up_bound"){
             if(solve2Approx == 1){
               tempAlgName="up_bound";
               for(double c : costOf2Approx) tempCost += c;
               tempCost = 2*tempCost;
//               for(double p : prizeOf2Approx) tempPrize += p;
               tempPrize = (1+0.5)*capacity*numVehicles;
//               4*tempPrize;//This is based on the fact that 2-approx looses a factor of 2 when capacity is added to prize 
               for(double s : sizeOf2Approx) tempSize += s;
               for(double r: runningTimeOf2Approx) tempRunningTime += r;
             }
             else continue;
         }
         else if(vecAlgs[i] == "Greedy"){
               tempAlgName = "greedy";
               for(double c : costOfGreedy) tempCost += c;
               tempCost = tempCost;
               for(double p : prizeOfGreedy) tempPrize += p;
               for(double s : sizeOfGreedy) tempSize += s;
               for(double r: runningTimeOfGreedy) tempRunningTime += r;
          }
         else if(vecAlgs[i] == "GreedyRandom"){
               tempAlgName = "greedyrand";
               for(double c : costOfGreedyRandom) tempCost += c;
               tempCost = tempCost;
               for(double p : prizeOfGreedyRandom) tempPrize += p;
               for(double s : sizeOfGreedyRandom) tempSize += s;
               for(double r: runningTimeOfGreedyRandom) tempRunningTime += r;
         }
         else {
            if(vecAlgs[i] == "SqrBApprox"){
               tempAlgName = "sqrB-approx";
               tempCost = 0;
               for(double c : costOfSqrBApprox) tempCost += c;
               tempCost = tempCost;
               for(double p : prizeOfSqrBApprox) tempPrize += p;
               for(double s : sizeOfSqrBApprox) tempSize += s;
               for(double r: runningTimeOfSqrBApprox) tempRunningTime += r;
            }
         }
         ofs<<graphT<<std::setw(25)<<numVerticesEdges.first<<std::setw(25)<< numVehicles <<std::setw(25)<<capacity<<std::setw(25)<<valBudget<<std::setw(25)<<tempAlgName<<std::setw(25)<<tempPrize<<std::setw(25)<<tempCost<<std::setw(25)<<tempSize<<std::setw(25)<<tempRunningTime<<std::setw(25)<<numThresholds<<std::setw(25)<<factorViolation<<std::setw(25)<<subsetSizeFactor<<std::endl;
         
     }
     //}//Running in parallel algorithms
     ofs.close();
}


void OutputResults::removeCoveredVertices(std::vector<int> coveredVertices, std::vector<std::pair<double, double>>& verticesCoordinationOfAlg, std::vector<double>& verticesPrizesOfAlg, std::vector<double>& verticesSizeOfAlg, std::vector<double>& verticesServicesTimeOfAlg){
         //First Sort coveredVertices in non decreasing order 
         sort(coveredVertices.begin(), coveredVertices.end()); 
//         for (int v : coveredVertices) std::cout<<" "<<v <<" ";
         //Remov covered vertices   
         //Note that we don't remove the deppot 0
         for(int i = coveredVertices.size() - 1; i >= 1; i--) {
                
                verticesCoordinationOfAlg.erase(verticesCoordinationOfAlg.begin() + coveredVertices[i]);
                verticesPrizesOfAlg.erase(verticesPrizesOfAlg.begin() + coveredVertices[i]);
                verticesSizeOfAlg.erase(verticesSizeOfAlg.begin() + coveredVertices[i]);
                verticesServicesTimeOfAlg.erase(verticesServicesTimeOfAlg.begin() + coveredVertices[i]);
        }
//        std::cout<<"\n";
}


//In binarySearch, we guess the opt for COTP and add the sizes to prizes and run each algorithm for the Orienteeiing problem
void OutputResults::binarySearch(std::vector<std::pair<double, double>> verticesCoordinationsOfAlg, std::vector<double> verticesPrizesOfAlg, std::vector<double> verticesSizesOfAlg, std::vector<double> verticesServicesTimeOfAlg, std::string nameOfAlg, std::vector<int>& verticesOfAlg, double& costOfAlg, double& prizeOfAlg, double& sizeOfAlg){

        
        double MaxRepeated = 0;
        if(nameOfAlg == "2Approx"){
//            for(double p : verticesPrizesOfAlg) prizeOfG += p;
            MaxRepeated = 2.0;//(1+0.5)*capacity;//static_cast<int>(std::log2((1+0.5)*capacity));
        }
        else MaxRepeated = 0.0;
        std::vector<int> verticesOfG;
        double sizeOfG = 0;
        double prizeOfG = 0;
        for (int i = 0; i < verticesPrizesOfAlg.size(); i++) {
            sizeOfG += verticesSizesOfAlg[i];
            prizeOfG += verticesPrizesOfAlg[i];
            verticesOfG.push_back(i);
        }

        for(double countBinary = 0.0; countBinary <= MaxRepeated; countBinary += 0.2) { 
           guessedUpperOpt = countBinary;//std::pow(2, countBinary);
//           //////////////////////////////////////////////////////////////////////////////////////////
           double coeficient =  countBinary;//guessedUpperOpt/(2*capacity);
//           double mst_w = findTSP2Approx(verticesOfG, edgesWeightsOfAlg, verticesCoordinationsOfAlg, verticesServicesTimeOfAlg);
           double tsp_w = findTSP(verticesOfG, verticesCoordinationsOfAlg, verticesServicesTimeOfAlg);
           if (tsp_w <= budget) { 
                  std::cout << "\n-MST weight is at most the budget and is equal to: " << tsp_w <<"\n";
                  
                  double costOfMST = 0;
                  double prizeOfMST = 0;
                  double sizeOfMST = 0;
//                  double sizeOfMST = sizeTree(mst, verticesSizesOfAlg, tempVerticesOfMST);
                  std::vector<int> verticesOfMSTRespectingCapacity;
                  if(sizeOfG <= capacity){
                       sizeOfMST = sizeOfG;
                       costOfMST = tsp_w;
                       prizeOfMST = prizeOfG;
                       verticesOfMSTRespectingCapacity = verticesOfG;              
                       
                  }
                  else{
//                       std::vector<std::vector<double>> edgesWeights = coordinationsToEdges(verticesCoordinationsOfAlg);
                       dpForSize(verticesOfG, capacity, verticesCoordinationsOfAlg, verticesPrizesOfAlg, verticesSizesOfAlg, verticesServicesTimeOfAlg, sizeOfMST, prizeOfMST, costOfMST, verticesOfMSTRespectingCapacity);
                  }
                  verticesOfAlg = verticesOfMSTRespectingCapacity;
                  costOfAlg = costOfMST;
                  prizeOfAlg = prizeOfMST;
                  sizeOfAlg = sizeOfMST;
                  break;
                  
           }
           else {
                 if(nameOfAlg == "2Approx") {
                    
                    //Compute the edge cost based on vertices coordinations
                    std::vector<std::vector<double>> edgesWeights = coordinationsToEdges(verticesCoordinationsOfAlg);
                    Graph G;
                    G.addVertex(0, 4*prizeOfG);//For 2-approx, We set the prize of depot to high to include it in the solution
                    for(int i = 1; i < verticesPrizesOfAlg.size(); i++) {
                           double modifyPriz = 100*(verticesPrizesOfAlg[i] - coeficient*verticesSizesOfAlg[i]);
                           if(modifyPriz < 1) 
                              G.addVertex(i, 1);
                           else  G.addVertex(i, modifyPriz);
                    }
              
                    //Add edges to the graph
                    for(int head = 0; head < edgesWeights.size(); head++)
                       for(int tail = head + 1; tail < edgesWeights[head].size(); tail++)
                           G.addEdge(head, tail, edgesWeights[head][tail]);
                           
                    double sizeOf2ApproxTemp = 0;
                    double prizeOf2ApproxTemp = 0;
                    double costOf2ApproxTemp = 0;
                    std::vector<int> verticesOf2ApproxTemp;
                    compute2Approx(G, edgesWeights, verticesCoordinationsOfAlg, verticesPrizesOfAlg, verticesSizesOfAlg, verticesServicesTimeOfAlg, verticesOf2ApproxTemp, costOf2ApproxTemp, prizeOf2ApproxTemp, sizeOf2ApproxTemp);
                    if(costOf2ApproxTemp <= 0.5*budget && prizeOfAlg < prizeOf2ApproxTemp){
                        verticesOfAlg = verticesOf2ApproxTemp; 
                        costOfAlg = costOf2ApproxTemp;; 
                        prizeOfAlg = prizeOf2ApproxTemp; 
                        sizeOfAlg = sizeOf2ApproxTemp;
                    }
                 }
                 else if(nameOfAlg == "SqrBApprox"){
                    double sizeOfSqrBApproxTemp = 0;
                    double prizeOfSqrBApproxTemp = 0;
                    double costOfSqrBApproxTemp = 0;
                    std::vector<int> verticesOfSqrBApproxTemp;
                    computeSqrBApprox(verticesCoordinationsOfAlg, verticesPrizesOfAlg, verticesSizesOfAlg, verticesServicesTimeOfAlg, verticesOfSqrBApproxTemp, costOfSqrBApproxTemp, prizeOfSqrBApproxTemp, sizeOfSqrBApproxTemp);
                    if(prizeOfAlg < prizeOfSqrBApproxTemp){
                        verticesOfAlg = verticesOfSqrBApproxTemp; 
                        costOfAlg = costOfSqrBApproxTemp; 
                        prizeOfAlg = prizeOfSqrBApproxTemp; 
                        sizeOfAlg = sizeOfSqrBApproxTemp;
                    }
                 }
                 
                 else if(nameOfAlg == "Greedy"){
                    double sizeOfGreedyTemp = 0;
                    double prizeOfGreedyTemp = 0;
                    double costOfGreedyTemp = 0;
                    std::vector<int> verticesOfGreedyTemp;
                    computeGreedy(verticesCoordinationsOfAlg, verticesPrizesOfAlg, verticesSizesOfAlg, verticesServicesTimeOfAlg, verticesOfGreedyTemp, costOfGreedyTemp, prizeOfGreedyTemp, sizeOfGreedyTemp);
                    if(prizeOfAlg < prizeOfGreedyTemp){
                        verticesOfAlg = verticesOfGreedyTemp; 
                        costOfAlg = costOfGreedyTemp; 
                        prizeOfAlg = prizeOfGreedyTemp; 
                        sizeOfAlg = sizeOfGreedyTemp;
                    }
                 }
                 else
                    if(nameOfAlg == "GreedyRandom"){ 
                    for (int i = 0; i < 10; i++) {
                       double sizeOfGreedyRandomTemp = 0;
                       double prizeOfGreedyRandomTemp = 0;
                       double costOfGreedyRandomTemp = 0;
                       std::vector<int> verticesOfGreedyRandomTemp;    
                       computeGreedyRandom(verticesCoordinationsOfAlg, verticesPrizesOfAlg, verticesSizesOfAlg, verticesServicesTimeOfAlg, verticesOfGreedyRandomTemp, costOfGreedyRandomTemp, prizeOfGreedyRandomTemp, sizeOfGreedyRandomTemp);
                       if(prizeOfAlg < prizeOfGreedyRandomTemp){
                           verticesOfAlg = verticesOfGreedyRandomTemp; 
                           costOfAlg = costOfGreedyRandomTemp; 
                           prizeOfAlg = prizeOfGreedyRandomTemp; 
                           sizeOfAlg = sizeOfGreedyRandomTemp;
                       }
                    }
                 }
           }
     }//End-Of-Binary-For
}//End-Of-BinaryFunction

     
void OutputResults::compute2Approx(Graph G, std::vector<std::vector<double>> edgesWeightsOfAlg, std::vector<std::pair<double, double>> verticesCoordinationsOfAlg, std::vector<double> verticesPrizesOfAlg, std::vector<double> verticesSizesOfAlg, std::vector<double> verticesServicesTimeOfAlg, std::vector<int>& verticesOf2ApproxTemp, double& costOf2ApproxTemp, double& prizeOf2ApproxTemp, double& sizeOf2ApproxTemp) {
     
     
                  
     //set the parameters
     double l = 0.0020055;
     double r = 0.0020057;
     std::list<std::shared_ptr<Subset> > subsetsL = growSubsets(G,l);
     std::list<std::shared_ptr<Subset> > subsetsR = growSubsets(G,r);
     double wplus = reverseDelete(subsetsL,false);
     double wminus = reverseDelete(subsetsR,false);
     std::list<std::shared_ptr<Edge> > edges;
     double upper = 0;
     int recursions = 0;
     double lambda;
     bool found;
     std::pair<double, double> approxPrizeCost =  PD(G, budget, edges, upper, recursions,lambda, found, true);
//     prizeOf2ApproxTemp = approxPrizeCost.first;
//     costOf2ApproxTemp = approxPrizeCost.second;
     //Save the vertices returned by 2-approx with a budget B before applying the capaccity constraint
     std::vector<int> selectedVerticesOf2Approx;
     sizeOf2ApproxTemp = sizeTree(edges, verticesSizesOfAlg, selectedVerticesOf2Approx);
     
     //Save the vertices in selectedVerticesOf2Approx with the maximum prize that respects the capacity constraint
     dpForSize(selectedVerticesOf2Approx, capacity, verticesCoordinationsOfAlg, verticesPrizesOfAlg, verticesSizesOfAlg, verticesServicesTimeOfAlg, sizeOf2ApproxTemp, prizeOf2ApproxTemp, costOf2ApproxTemp, verticesOf2ApproxTemp);
}


void OutputResults::computeSqrBApprox(std::vector<std::pair<double, double>> verticesCoordinationsOfAlg, std::vector<double> verticesPrizesOfAlg, std::vector<double> verticesSizesOfAlg, std::vector<double> verticesServicesTimeOfAlg, std::vector<int>& verticesOfSqrBApproxTemp, double& costOfSqrBApproxTemp, double& prizeOfSqrBApproxTemp, double& sizeOfSqrBApproxTemp){

    SqrBApproxAlg sqrB(capacity, budget, verticesCoordinationsOfAlg, verticesPrizesOfAlg, verticesSizesOfAlg, verticesServicesTimeOfAlg, numThresholds, subsetSizeFactor, numThreads);
    verticesOfSqrBApproxTemp = sqrB.getVerticesSqrB();
    sizeOfSqrBApproxTemp = 0;
    prizeOfSqrBApproxTemp = 0;
    for (int v : verticesOfSqrBApproxTemp){
        sizeOfSqrBApproxTemp += verticesSizesOfAlg[v];
        prizeOfSqrBApproxTemp += verticesPrizesOfAlg[v];
    }
//    std::cout<< "\n"<<0 <<" "<<verticesSizesOfAlg[0] <<"\n";
//    costOfSqrBApproxTemp = findTSP2Approx(verticesOfSqrBApproxTemp, edgesWeightsOfAlg, verticesCoordinationsOfAlg, verticesServicesTimeOfAlg);
    costOfSqrBApproxTemp = findTSP(verticesOfSqrBApproxTemp, verticesCoordinationsOfAlg, verticesServicesTimeOfAlg);
}


void OutputResults::computeGreedy(std::vector<std::pair<double, double>> verticesCoordinationsOfAlg, std::vector<double> verticesPrizesOfAlg, std::vector<double> verticesSizesOfAlg, std::vector<double> verticesServicesTimeOfAlg, std::vector<int>& verticesOfGreedyTemp, double& costOfGreedyTemp, double& prizeOfGreedyTemp, double& sizeOfGreedyTemp) {

    GreedyAlg greedyAlg(capacity, budget, verticesCoordinationsOfAlg, verticesPrizesOfAlg, verticesSizesOfAlg, verticesServicesTimeOfAlg, numThreads);
    //Get the selected vertices in SqrBApproxAlg before appplying the capacity constraint
    verticesOfGreedyTemp = greedyAlg.getVerticesGreedy();
    
    sizeOfGreedyTemp = 0;
    prizeOfGreedyTemp = 0;
    for (int v : verticesOfGreedyTemp){
        sizeOfGreedyTemp += verticesSizesOfAlg[v];
        prizeOfGreedyTemp += verticesPrizesOfAlg[v];
    }
//    costOfGreedyTemp = findTSP2Approx(verticesOfGreedyTemp, edgesWeightsOfAlg, verticesCoordinationsOfAlg, verticesServicesTimeOfAlg);
    costOfGreedyTemp = findTSP(verticesOfGreedyTemp, verticesCoordinationsOfAlg, verticesServicesTimeOfAlg);

}

void OutputResults::computeGreedyRandom(std::vector<std::pair<double, double>> verticesCoordinationsOfAlg, std::vector<double> verticesPrizesOfAlg, std::vector<double> verticesSizesOfAlg, std::vector<double> verticesServicesTimeOfAlg, std::vector<int>& verticesOfGreedyRandomTemp, double& costOfGreedyRandomTemp, double& prizeOfGreedyRandomTemp, double& sizeOfGreedyRandomTemp){
//    std::cout<<"\nBefore GreedyRandom\n";
    GreedyRandomAlg greedyRandomAlg(capacity, budget, verticesCoordinationsOfAlg, verticesPrizesOfAlg, verticesSizesOfAlg, verticesServicesTimeOfAlg);
//    std::cout<<"\nAfter GreedyRandom\n";
    //Get the selected vertices in SqrBApproxAlg before appplying the capacity constraint
    verticesOfGreedyRandomTemp = greedyRandomAlg.getVerticesGreedyRandom();
    
    sizeOfGreedyRandomTemp = 0;
    prizeOfGreedyRandomTemp = 0;
    for (int v : verticesOfGreedyRandomTemp){
        sizeOfGreedyRandomTemp += verticesSizesOfAlg[v];
        prizeOfGreedyRandomTemp += verticesPrizesOfAlg[v];
    }
//    costOfGreedyRandomTemp = findTSP2Approx(verticesOfGreedyRandomTemp, edgesWeightsOfAlg, verticesCoordinationsOfAlg, verticesServicesTimeOfAlg);
    costOfGreedyRandomTemp = findTSP(verticesOfGreedyRandomTemp, verticesCoordinationsOfAlg, verticesServicesTimeOfAlg);
    
//    //Save vertices in selectedVerticesOfGreedyRandom with the maximum prize respecting the capacity constraints
//    totalSizeVertices(G, edgesWeightsOfAlg, verticesPrizesOfAlg, verticesSizesOfAlg, selectedVerticesOfGreedyRandom, sizeOfGreedyRandomTemp, prizeOfGreedyRandomTemp, costOfGreedyRandomTemp, verticesOfGreedyRandomTemp);
}




// Calculate the total size of vertices in tree
double OutputResults::sizeTree(std::list<std::shared_ptr<Edge>>& tree, std::vector<double> verticesSizesOfAlg, std::vector<int>& verticesTree){
    double c = 0;//total capacity
    std::vector<bool> visited(verticesSizesOfAlg.size(), false); 
    double sumWeight = 0; 
    std::vector<int> headV;
    std::vector<int> tailV;
    for (auto e:tree){
        int i = e->getHead(), j = e->getTail();
//        std::cout<< "("<<i <<", "<< j<<")  ";
        if(visited[i] == false){
            verticesTree.push_back(i);
            c += verticesSizesOfAlg[i];
            visited[i] = true;
            headV.push_back(i);
        }
        
        if(visited[j] == false){
            verticesTree.push_back(j);
            c += verticesSizesOfAlg[j];
            visited[j] = true;
            tailV.push_back(j);
        }
    }
    return c;
}

////prev
//void  OutputResults::findPrevOfTree(std::list<std::shared_ptr<Edge> > & edges, std::vector<bool>& visited, std::vector<int>& prev, int v){
//    // Mark v as visited and add to tour
//    visited[v] = true;
//    
//    // Iterate through edges
//    for (auto e:edges)
//    {
//        if ((e->getTail() == v) || (e->getHead() == v)){
//            int u = e->getOther(v);
//            if (visited[u] == false){
//                prev[u] = v;
//                findPrevOfTree(edges, visited, prev, u);
//            }
//        }
//    }
//    
//}