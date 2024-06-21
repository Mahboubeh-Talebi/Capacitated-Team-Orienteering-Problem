//
//  ReadFile.cpp
//  
//
//  Created by Alice Paul on 4/10/17.
//
//  MIT License
//  Copyright (c) 2020 alicepaul
//
//
//
#include<iostream>  
#include<vector> // for vector  
#include<algorithm> // for copy() and assign()  
#include<iterator> // for back_inserter
#include "ReadFile.hpp"
#define INF 0x3f3f3f3f

// reads the graphs from a file
int graphFromFile(std::string fileName, std::vector<std::pair<double, double>>& verticesCoordinations, std::vector<double>& verticesSizes, std::vector<double>& verticesPrizes, std::vector<double>& verticesServicesTime, int& numVehicles, double& budget, double& capacity) {
    std::string line;
    std::ifstream infile(fileName);
    std::ifstream infile_Second(fileName);
    
    int nodesAdded = 0;
    double x;//x coordination
    double y;//y coordination
    getline(infile,line);
//    std::cout<<line <<"\n";
    while(!(line.find("CUSTOMERDATA") != std::string::npos) && infile.peek() != EOF) {
      
     if (line.find("MAXVEHICLES") != std::string::npos) {
         std::istringstream iss(line);
         std::string str;
         iss>> str >> numVehicles;
     }
     else if (line.find("MAXCAPACITY") != std::string::npos) {
         std::istringstream iss(line);
         std::string str;
         iss>> str >> capacity;
     }
     else if (line.find("MAXDURATION") != std::string::npos) {
         std::istringstream iss(line);
         std::string str;
         iss>> str >> budget;
     }
     else 
       if (line.find("DEPOT") != std::string::npos) {
         std::istringstream iss(line);
         std::string str;
         iss>> str >> x >> y;
         verticesCoordinations.push_back(std::make_pair(x, y));
         verticesSizes.push_back(0);
         verticesPrizes.push_back(0);
         verticesServicesTime.push_back(0);
         ++nodesAdded;
       }
      getline(infile,line);
//      std::cout<<line <<"\n";
    }

    double p;//profit
    double s;//size
    double t;//service time
    while (getline(infile, line)){
          std::istringstream iss(line);
          iss>> x >> y >> s >> t >> p;
          verticesCoordinations.push_back(std::make_pair(x, y));
          verticesSizes.push_back(s);
          verticesPrizes.push_back(p);
          verticesServicesTime.push_back(t);
          ++nodesAdded;
//          std::cout<< x <<" "<< y<< " " << s <<" " << t <<" "<< p <<	"\n";
    }//END-While
   if(verticesSizes.empty() && verticesPrizes.empty() && verticesCoordinations.empty()){
        //int nodesAdded = 0;
        getline(infile_Second,line);
    //    std::cout<<line <<"\n";
//        int count = 0;
        while(!(line.find("NODE_COORD_SECTION") != std::string::npos) && infile_Second.peek() != EOF) {
          getline(infile_Second,line);
        }
        getline(infile_Second,line);
        std::cout<<line <<"\n";
        while(!(line.find("DEMAND_SECTION") != std::string::npos) && infile_Second.peek() != EOF){
              //std::string str;
              int e; //edge e
              double x;
              double y;
              std::istringstream iss(line);
              
              iss>> e >> x >> y;
              verticesCoordinations.push_back(std::make_pair(x, y));
    //          std::cout<< x <<" "<< y<<"\n";
              getline(infile_Second,line);
    //          std::cout<<line <<"\n";
              //++nodesAdded;
    
        }//END-While
        
        getline(infile_Second,line);
        while(!(line.find("PRIZE_SECTION") != std::string::npos) && infile_Second.peek() != EOF){
              //std::string str;
              int v; //vertex v
              double c;//capacity c
              std::istringstream iss(line);
              
              iss>> v >> c;
              verticesSizes.push_back(c);
    //          std::cout<< v <<" " << c <<"\n";
              getline(infile_Second,line);
    //          std::cout<<line <<"\n";
    //          ++nodesAdded;
    
        }//END-While
        getline(infile_Second,line);
        while(!(line.find("DEPOT_SECTION") != std::string::npos) && infile_Second.peek() != EOF){
              //std::string str;
              int v; //vertex v
              double p;//prize p
              std::istringstream iss(line);
              
              iss>> v >> p;
              verticesPrizes.push_back(p);
    //          std::cout<< v <<" " << c <<"\n";
              getline(infile_Second,line);
    //          std::cout<<line <<"\n";
              ++nodesAdded;
    
        }//END-While    
}
    std::cout<<"nodes added:"<<nodesAdded<<"\n";

    return nodesAdded;
    
}

//std::pair<int, int> graphFromFile(const std::string &fileName, double& meanEdgeWeight, int& numNodes, std::vector<std::vector<double>>& edgesWeights, std::vector<double>& verticesPrizes, std::vector<double>& verticesCapacities){
//    std::string line;
//    std::ifstream infile(fileName);
//    double totalEdgeWeight = 0;
//    int nodesAdded = 0;
//    int edgesAdded = 0;
//    int dimension;
//    int edgeWeight;
//    int head;
//    int tail;
//    double distance;
//    
//    //std::string line;
//    // for (int i=0; i<dimension; ++i){
//    // for (int j=0; j<dimension; ++j){
//    while(std::getline(infile,line)){
//       if(line[0]=='N' && line[1]=='o'){
//           int num_Nodes;//number of nodes
//           std::string str;
//           std::istringstream iss(line);
//           iss>> str >> num_Nodes;//dimension;
//           std::vector<std::vector<double>> vec(num_Nodes, std::vector<double> (num_Nodes, INF));
//           //copy(vec.begin(), vec.end(), back_inserter(edgesWeights));
//           edgesWeights=vec;
//       }
//       else if(line[0]=='N' && line[1]==' '){
//          int u;
//          double p;//prize
//          std::string str;
//          std::istringstream iss(line);
//          iss>> str >> u >> p;//dimension;
//          verticesPrizes.push_back(p);
////          G.addVertex(u, p);
//          nodesAdded++;
//       }
//       else if(line[0]=='S' && line[1]==' '){
//          int u;
//          double c;//capacity
//          std::string str;
//          std::istringstream iss(line);
//          iss>> str >> u >> c;//dimension;
//          verticesCapacities.push_back(c);
//       }
//       else 
//        if(line[0]=='E' && line[1]==' '){
//         std::istringstream iss(line);
//         std::string str;
//         iss>> str >> head >> tail >> distance;
//         //std::cout<<"adding edge:"<<head<<","<<tail<<","<<distance<<std::endl;
//         //G.addEdge(head,tail,distance);
//         edgesWeights[head][tail]=distance;
//         edgesWeights[tail][head]=distance;
//         //edgesWeights.push_back(std::make_pair());
//         totalEdgeWeight += distance;
//         edgesAdded++;
//       }
//    }//END-While
//    meanEdgeWeight = totalEdgeWeight/edgesAdded;
//    std::cout<<"meanEdgeWeight: "<<meanEdgeWeight<<std::endl;
//    
//    //std::cout<<"graph:"<<graphName<<std::endl;
//    numNodes = nodesAdded;
//    std::cout<<"nodes added:"<<nodesAdded<<"\n"
//    <<"edges added:"<<edgesAdded<<std::endl;
//    return std::make_pair(numNodes, edgesAdded);
    
//}

