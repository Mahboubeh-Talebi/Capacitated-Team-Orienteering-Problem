//
//  Graph.hpp
//  Prize Collecting TSP
//
//  Created by Alice Paul on 3/14/17.
//
//  MIT License
//  Copyright (c) 2020 alicepaul
//

#ifndef Graph_hpp
#define Graph_hpp

#include <stdio.h>
#include <map>
#include <list>
#include <iostream>
#include <memory>
#include <set>
#include <climits>
#include <vector>
#include <queue>
#include <algorithm>
#include "NodeInfo.h"
/* -------------------------EDGE--------------------------*/


class Edge
{
private:
    int head;                               // head vertex id
    int tail;                               // tail vertex id
    double weight;                          // weight of edge
    
public:
    // Constructors and Destructors
    Edge(int h, int t, double w);
    ~Edge();
    
    // Get Functions
    int getHead() const {return head;}
    int getTail() const {return tail;}
    double getWeight() const {return weight;}
    void setWeight(double w) {weight = w;}
    int getOther(int id);
    
    // Print Function
    friend std::ostream& operator<< (std::ostream &out, const Edge &edge);
    
    // Comparison Function
    friend bool operator==(const Edge &e1, const Edge &e2);
};


/* -------------------------VERTEX--------------------------*/


class Vertex
{
private:
    int id;                                     // id
    double prize;                                  // prize value
    std::list<std::shared_ptr<Edge> > neighbors; // list of incident edge pointers
    double degree;                              // sum of weights for incident edges
    
public:
    // Constructors and Destructors
    Vertex(int i);
    Vertex(int i, double p);
    ~Vertex();
    
    // Get Functions
    int getId() const {return id;}
    double getPrize() const {return prize;}
    std::list<std::shared_ptr<Edge> > const & getIncEdges() const {return neighbors;}
    double getDegree() const {return degree;}
    
    // Add and Remove Edges
    void addEdge(std::shared_ptr<Edge> e);
    void removeEdge(std::shared_ptr<Edge> e);
};


/* -------------------------GRAPH--------------------------*/


class Graph
{
private:
    std::map<int, std::shared_ptr<Vertex> > vertexmap;   // map of vertex ids to vertex data structures
    std::list<std::shared_ptr<Edge> > edges;             // list of edge pointers
    std::list<int> vertices;                            // list of vertex ids
    double W;                                           // total weight of edges
    double P;                                              // sum of prizes of vertices
    
    
    // Used internally for removing vertices
    void removeVertexLists(int id);
    
public:
    // Constructors and Destructors
    Graph();
    ~Graph();
    Graph(const Graph &G);                              // copy constructor
    Graph(const Graph &G, const std::list<int> &S);            // subgraph G(S)
    
    // Get Functions
    double getWeight() const {return W;}
    double getPrize() const {return P;}
    std::list<int> const & getVertices() const {return vertices;}
    std::list<std::shared_ptr<Edge> > const & getEdges() const {return edges;}
    std::shared_ptr<Vertex> const & getVertex(int i) const;
    double getVertexDegree(int it) const;
    
    // Add and Remove Functions
    void addVertex(int id);
    void addVertex(int id, double prize);
    void addEdge(int id1, int id2, double weight);
    void deleteVertex(int id);
    
    // Print Function
    friend std::ostream& operator<< (std::ostream &out, const Graph &graph);
    
    // Minimum spanning tree
    // Returns weight and tree is saved to edges
    double MST(std::list<std::shared_ptr<Edge> > & edges) const;
    
};

//Create edges from the coordinations of vertices
std::vector<std::vector<double>> coordinationsToEdges(std::vector<std::pair<double, double>> verticesCoordinations);

//build the graph by vertices prizes and edges weights
Graph buildGraph(std::vector<double> verticesPrizes, std::vector<std::vector<double>> edgesWeights, std::vector<int>& verticesOfG);

// Calculate tour length
double getTourLength(const Graph & G, const std::vector<int> & tour);

// DFS
void DFS(std::list<std::shared_ptr<Edge> > & edges, std::map<int,bool> & visited,
         std::list<int> & tour, int v);

// Find ordered list of tour given a tree of edges
std::vector<int> tourList(std::list<std::shared_ptr<Edge> > & edges);

// Find length of shortest path from i to j in G
double shortestPath(const Graph & G, int i, int j);

double findMinRoute(std::vector<std::vector<double>> tsp);
double tspDynamic(const std::vector<std::vector<double>>& cities, int pos, int visited, std::vector<std::vector<double>>& state);

double findTSP2Approx(std::vector<int> S_u, std::vector<std::pair<double, double>> verticesCoordinations);

double findTSP(std::vector<int> S_u, std::vector<std::pair<double, double>> verticesCoordinations, std::vector<double> verticesServicesTime);

//Find the shortest path from each vertex to all other vertices
void all_to_all_shortest_paths(const Graph G, std::vector<int> G_Vertices, std::vector<double> verticesServicesTime, std::vector<std::vector<int>>& prev_u, std::vector<std::vector<double>>& dist_u);
void shortestPath_to_DestNodes(const Graph & G, int src, std::vector<int> destVertices, std::vector<double> verticesServicesTime, std::vector<int>& prev, std::vector<double>& dist);

//The dynamic programmign for finding a subset of vertices in verticesOfAlg with size at most the capacity and maximum prize
void dpForSize(std::vector<int> verticesOfAlg, double capacity, std::vector<std::pair<double, double>> verticesCoordinationsOfAlg, std::vector<double> verticesPrizesOfAlg, std::vector<double> verticesSizesOfAlg, std::vector<double> verticesServicesTimeOfAlg, double& totalSize, double& totalPrize, double& totalCost, std::vector<int>& verticesOfAlgRespectingCapacity);

//Add more vertices to the solution if it is possible
void addMoreVertices(std::vector<ratioPrizeCost> vecRatioPrizeCost, std::vector<int>& verticesRespectingCapacity, double& size_x, double& prize_x, double& cost_x, double budget, double capacity, std::vector<int> vectTemp, std::vector<std::pair<double, double>> verticesCoordinationsOfAlg, std::vector<double> verticesPrizes, std::vector<double> verticesSizes, std::vector<double> verticesServicesTime, std::string algName);

double subGraphMST(std::vector<int> vecTempVertices, std::vector<std::vector<double>> edgesWeightsOfAlg);



#endif /* Graph_hpp */
