//
//  ReadFile.hpp
//  
//
//  Created by Alice Paul on 4/10/17.
//
//  MIT License
//  Copyright (c) 2020 alicepaul
//
#include <stdio.h>
#include <map>
#include <list>
#include <iostream>
#include <fstream>
#include <memory>
#include <set>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>

#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>


#include "Graph.hpp"
#include "Subset.hpp"
#include "Subroutine.hpp"
#include "PD.hpp"
#ifndef ReadFile_hpp
#define ReadFile_hpp

#include <stdio.h>

int graphFromFile(std::string fileName, std::vector<std::pair<double, double>>& coordination, std::vector<double>& verticesSizes, std::vector<double>& verticesPrizes);
//std::pair<int, int> graphFromFile(const std::string &fileName, double& meanEdgeWeight, int& numNodes, std::vector<std::vector<double>>& edgesWeights, std::vector<double>& verticesPrizes, std::vector<double>& verticesCapacities);


#endif /* ReadFile_hpp */
