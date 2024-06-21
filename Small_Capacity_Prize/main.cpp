//
//  main.cpp

//  Created by Esmaeil Delfaraz on 3/14/17.



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
#include<float.h>
#include <filesystem>

using namespace std::chrono;
#include "Graph.cpp"
#include "ReadFile.cpp"
#include "output.cpp"
  
//Read_OutPut reads the input file and provides the graph to the OutPut class in order to run the algorithms
void Read_OutPut(std::string inputFileName, std::string graphT, int numVehicles, double budget, double capacity, int numThresholds, int numThreads, int solve2Approx, double subsetSizeFactor, double factorViolation, std::string outputDirectory){
         
         std::vector<double> verticesPrizes;
         std::vector<double> verticesSizes;
         std::vector<std::pair<double, double>> verticesCoordinations;
         int numVertices = graphFromFile(inputFileName, verticesCoordinations, verticesSizes, verticesPrizes);
         std::vector<std::vector<double>> edgesWeights(numVertices, std::vector<double> (numVertices, 0));
         
         /* 
           In the rooted Capacitated Team Orienteering Problem (CTOP), we only consider those vertices at a distance no more than budget/2 from the root vertex.
           As a vertex with a distance more than budget/2 cannot be in a feasible solution to CTOP
         */
         double dist;
//         double totalPrize = 0;
//         for(double p : verticesPrizes) totalPrize += p;
//         std::cout<<"\nTotal prize is : "<<totalPrize <<"\n";
         for(int i = verticesPrizes.size() - 1; i >= 1; i--){
                dist = std::sqrt(std::pow((verticesCoordinations[i].first - verticesCoordinations[0].first), 2) + std::pow((verticesCoordinations[i].second - verticesCoordinations[0].second), 2));
                if(dist > 0.5*budget){
                   verticesCoordinations.erase(verticesCoordinations.begin() + i);
                   verticesPrizes.erase(verticesPrizes.begin() + i);
                   verticesSizes.erase(verticesSizes.begin() + i);
                }
         }    
         size_t lastIndexPoint = inputFileName.find_last_of("."); 
         size_t lastIndexSlash = inputFileName.find_last_of("/");
         graphT  = inputFileName.substr(lastIndexSlash + 1, lastIndexPoint);
         //If number of vertices is more than zero, run the algorithm
         if(!verticesSizes.empty() && !verticesPrizes.empty() && !verticesCoordinations.empty()){
               OutputResults ORClass(graphT, numVertices, numVehicles, verticesCoordinations, verticesPrizes, verticesSizes, budget, capacity, numThresholds, numThreads, solve2Approx, subsetSizeFactor, factorViolation, outputDirectory);
         }   
         else{
                  std::cout<< "\nThe input or the command line or both are not valid. Please try again...\n";
                  //return 1;
        }
//        }//Running in parallel algorithms
}

int main(int argn, char *argv[]){
    
    if(argn < 2){
     		std::cout << "usage " << argv[0] << "\n"
     				<< " -d the directory to save the generated graph[gen] or the directory of collection of files for which algorithms will be run on[algs]" << "\n"
     				<< " -i input file name [when d is not set]" << "\n"
     				<< " -K number of vehicles to use [algs]" << "\n"
     				<< " -B budget[algs]"  << "\n"
         << " -C capacity[algs]"  << "\n"
         << " -f directory/outputfilename in which to save the results of algorithms[algs]" << "\n"
         << " -l number of of thresholds to run the sqrB_approx algo (default=1)[algs]" << "\n"
         << " -v the factor by which we set the size of subset of vertices we select in sqrB-approx (default=1)[algs]" << "\n"
         << " -e the factor by which we violate the budget in algorithms (default=0)[algs]" << "\n"
         << " -s  s=1 means that we solve 2-approx, s=0 means otherwise (default=0) [algs]" << "\n"
         << " -z number of threads to compute approx_alg (default=1)[algs]" << "\n";
     		return 0;
 }
  std::string generatedFileDirectory, inputFileName, outputDirectory, graphT, runT; 
  int n, k, solve2Approx=0, numVehicles, numThreads, numThresholds, minPrize, maxPrize;
  double p=1, a, b, minWeight, maxWeight, budget, capacity, subsetSizeFactor=1, factorViolation = 0;
  int c=1;
      	while ((c = getopt(argn, argv, "r:g:K:n:k:p:a:b:m:M:w:W:i:d:B:C:f:l:v:e:s:z:")) != -1) {
    		switch(c) {
       case 'i':
           inputFileName = strdup(optarg);
           break;
       case 'd':
           generatedFileDirectory = strdup(optarg);
           break;
       case 'K':
         numVehicles = atoi(optarg);
         std::cout<<"Number of vehicles K = "<< numVehicles <<"\n";
         break;
       case 'B':
         budget = atof(optarg);
         std::cout<<"Budget B= "<<budget<<"\n";
         break;
       case 'C':
         capacity = atof(optarg);
         std::cout<<"Capacity C= "<<budget<<"\n";
         break;
       case 'e':
         factorViolation = atof(optarg);
         std::cout<<"The factor by which we violate the budget = "<<budget<<"\n";
         break;
       case 'l':
         numThresholds = atoi(optarg);
         std::cout<<"Number of thresholds to run the sqrB_approx alg is "<< numThresholds<<" \n";
         if(numThresholds<1)
         {
           numThresholds = 1;
           std::cout << "As default l is set to 1" << std::endl;
         }
         break;
       case 'v':
         subsetSizeFactor = atof(optarg);
         std::cout<<"The factor by which we set the size S_u in sqrB_approx is "<< subsetSizeFactor<<" \n";
         if(subsetSizeFactor<1)
         {
           subsetSizeFactor = 1;
           std::cout << "As default l is set to 1" << std::endl;
         }
         break;
       case 's':
         solve2Approx = atoi(optarg);
         std::cout<<"Solve 2-approx alg s= "<< solve2Approx<<" \n";
         if(solve2Approx != 1 && solve2Approx !=0)
         {
           solve2Approx = 0;
           std::cout << "As default s is set to 0" << std::endl;
         }
         break;
       case 'z':
         numThreads = atoi(optarg);
         std::cout<<"Number of threads z to compute approx_alg "<< numThreads<<" \n";
         if(numThreads<1)
         {
           numThreads = 1;
           std::cout << "As default z is set to 1" << std::endl;
         }
         break;
       case 'f':
         outputDirectory = strdup(optarg);
         std::cout<<"output file name f= "<<outputDirectory<<"\n";
         break;
       default:
    			  std::cout<<"error"<<std::endl;
  		  }
  	}
  	    
  std::map<std::string,int> RunType;
  const int GEN = 1, ALG=2;
  const std::string GENstr = "gen", ALGstr = "algs";
  RunType[GENstr] = GEN;
  RunType[ALGstr] = ALG;
//  switch(RunType[runT]) {
//      case GEN:{
//           Graph_Generator CGGGraph_Generator(graphT, n, k, a, b, p, minPrize, maxPrize, minWeight, maxWeight, generatedFileDirectory);
//           break;
//      }  
//      case ALG:{
 const std::filesystem::path path=generatedFileDirectory;
 if(!generatedFileDirectory.empty())
   for (auto const& dir_entry : std::filesystem::directory_iterator{path}){
       inputFileName=dir_entry.path().string();
       //Read the input file from the folder and give the graph to OutPut class in oder to run the algorithms
       Read_OutPut(inputFileName, graphT, numVehicles, budget, capacity, numThresholds, numThreads, solve2Approx, subsetSizeFactor, factorViolation, outputDirectory);
       }//END_FOR
 else{//when we set the inputfile directly
       //Read the input file which is set in the command line and give the graph to OutPut class in oder to run the algorithms
       Read_OutPut(inputFileName, graphT, numVehicles, budget, capacity, numThresholds, numThreads, solve2Approx, subsetSizeFactor, factorViolation, outputDirectory);
 }
//         break;
//      }
//      default:
//        std::cout << "\n" << "Error: Select " << GENstr << " or " << ALGstr <<" as run type" << "\n";
//        break;
//  }
   return 0;
}
