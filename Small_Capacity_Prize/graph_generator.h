#ifndef Graph_Generator_H_
#define Graph_Generator_H_

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/graph/small_world_generator.hpp>
#include <boost/graph/plod_generator.hpp>
#include <boost/random/linear_congruential.hpp>
#include <iostream>
#include <random>
#include <ctime>
#include <stdexcept>
#include <sstream>
#include <algorithm>
#include<vector>

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, boost::no_property, boost::property< boost::edge_weight_t, double> > UndirectedGraph;

class Graph_Generator{
     private:
         int n, k, minPrize, maxPrize;
         std::ofstream outputfile;//Generate weight and edge files
         double a, b, p, minWeight, maxWeight;
         std::string filename_t;//Will set the file name based on the machine time
         std::vector<double> val_EdgesWeights;
         UndirectedGraph val_g;
         std::vector<int> val_VerticesPrizes;
         //Generate Random Terminals
         void genRandomPrizes();
         //Generate random weight
         void genRandomWeight();
         //set edges and terminals
         void setInputFile();
         std::vector<bool> Visit_BFS(UndirectedGraph g, int src);
     public:
         //Generating Graph
         Graph_Generator(std::string graphT, int val_n, int val_k, double val_a, double val_b, double val_p, int val_minPrize, int val_maxPrize, double val_minWeight, double val_maxWeight, std::string val_outputfiledirectory);
 };
 
 #endif /* ALG_H_ */
