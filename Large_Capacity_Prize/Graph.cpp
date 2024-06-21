//
//  Graph.cpp
//  Prize Collecting TSP
//
//  Created by Alice Paul on 3/14/17.
//
//  MIT License
//  Copyright (c) 2020 alicepaul
//
#define INF 0x3f3f3f3f
#include "Graph.hpp"
#include "tsp.cpp"
/* -------------------------EDGE--------------------------*/

// Create an edge
Edge::Edge(int h, int t, double w)
{
    //std::cout << "Creating Edge (" << h << "," << t << ") \n";
    head = h;
    tail = t;
    weight = w;
}

//Delete an edge
Edge::~Edge()
{
    //std::cout << "Deleting Edge (" << head << "," << tail << ") \n";
}

// Find the other end of an edge
int Edge::getOther(int id)
{
    if (head == id){
        return tail;
    }
    else{
        return head;
    }
}


// Comparison of two edges
bool operator== (const Edge &e1, const Edge &e2)
{
    return (((e1.getHead()==e2.getHead()) && (e1.getTail() == e2.getTail())) && (e1.getWeight() == e2.getWeight()));
}

// How to print an edge
std::ostream& operator<< (std::ostream &out, const Edge &edge)
{
    out << " ( " << edge.head << " , " << edge.tail << " , " << edge.weight << ") ";
    return out;
}

/* -------------------------VERTEX--------------------------*/

// Create a vertex (edge list is empty)
Vertex::Vertex(int i)
{
    //std::cout << "Making vertex " << i << "\n";
    degree = 0;
    id = i;
    prize = 1; // Default prize
}

// Create a vertex (edge list is empty)
Vertex::Vertex(int i, double p)
{
    //std::cout << "Making vertex " << i << "\n";
    degree = 0;
    id = i;
    prize = p;
}

// Delete a vertex
Vertex::~Vertex()
{
    //std::cout << "Destroying vertex " << id << "\n";
}

// Add an edge to the list of neighbors
void Vertex::addEdge(std::shared_ptr<Edge> e)
{
    neighbors.push_back(e);
    if (e->getWeight() < 0){
        degree -= e->getWeight();
    }
    else{
        degree += e->getWeight();
    }
}

// Remove an edge from the list of neighbors
void Vertex::removeEdge(std::shared_ptr<Edge> e)
{
    neighbors.remove(e);
    if (e->getWeight() < 0){
        degree += e->getWeight();
    }
    else{
        degree -= e->getWeight();
    }
}



/* -------------------------GRAPH--------------------------*/

Graph::Graph(){
    W = 0.0; // Nothing else to do
    P = 0.0;
}

Graph::~Graph(){
    // Nothing to do
}

// copy constructor
Graph::Graph(const Graph &G){
    W = G.W;
    P = G.P;
    std::list<int>::const_iterator it;
    for(it = G.vertices.begin(); it != G.vertices.end(); it++){
        vertices.push_back(*it);
        int i = *it;
        double p = G.getVertex(i)->getPrize();
        std::shared_ptr<Vertex> v = std::make_shared<Vertex>(i,p);
        vertexmap[*it] = v;
    }
    
    std::list<std::shared_ptr<Edge> >::const_iterator it2;
    for(it2 = G.edges.begin(); it2 != G.edges.end(); it2++){
        int id1 = (**it2).getHead();
        int id2 = (**it2).getTail();
        double w = (**it2).getWeight();
        std::shared_ptr<Edge> e = std::make_shared<Edge>(id1,id2,w);
        vertexmap[id1]->addEdge(e);
        vertexmap[id2]->addEdge(e);
        edges.push_back(e);
    }
}

// subgraph constructor
Graph::Graph(const Graph &G, const std::list<int> &S){
    W = 0;
    P = 0;
    // add vertices in S
    for (auto x : S){
        double p = G.getVertex(x)->getPrize();
        addVertex(x,p);
    }
    
    // add edges with both endpts in S
    for(auto e : G.getEdges()){
        int id1 = e->getHead(), id2 = e->getTail();
        double w = e->getWeight();
        
        if ((vertexmap.count(id1) > 0) && (vertexmap.count(id2) > 0)){
            addEdge(id1, id2, w);
        }
    }
}


// Add a vertex to a graph with a prize
void Graph::addVertex(int id, double p)
{
    vertices.push_back(id);
    std::shared_ptr<Vertex> v = std::make_shared<Vertex>(id,p);
    vertexmap[id] = v;
    P += p;
}

// Add a vertex to a graph
void Graph::addVertex(int id)
{
    addVertex(id,1);
}


// Add an edge to a graph
void Graph::addEdge(int id1, int id2, double weight)
{
    if ((vertexmap[id1]==NULL) || (vertexmap[id2]==NULL)){
        return;
    }
    std::shared_ptr<Edge> e = std::make_shared<Edge>(id1, id2, weight);
    vertexmap[id1]->addEdge(e);
    vertexmap[id2]->addEdge(e);
    edges.push_back(e);
    if (weight < 0){
        W -= weight;
    }
    else{
        W += weight;
    }
}



// Return ptr to vertex
std::shared_ptr<Vertex> const & Graph::getVertex(int id) const
{
    if (vertexmap.count(id) == 0){
        throw std::invalid_argument("Vertex does not exist");
    }
    return vertexmap.at(id);
}

// Return vertex degree
double Graph::getVertexDegree(int id) const
{
    if (vertexmap.count(id) == 0){
        throw std::invalid_argument("Vertex does not exist");
    }
    return vertexmap.at(id)->getDegree();
}


// Remove vertex from list of vertices and vertex map
void Graph::removeVertexLists(int id){
    vertices.remove(id);
    vertexmap.erase(id);
}

// Delete a vertex from a graph
void Graph::deleteVertex(int id)
{
    std::shared_ptr<Vertex> p_v = vertexmap[id];
    if (p_v == NULL){
        return;
    }
    
    // Update prize
    P -= p_v->getPrize();
    
    // If no edges to go through just remove references and delete vertex
    if (p_v->getIncEdges().size() == 0){
        removeVertexLists(id);
        return;
    }
    
    // Iterate through edges
    std::list<std::shared_ptr<Edge> >::iterator it;
    for (it = edges.begin(); it != edges.end();)
    {
        if (((*it)->getHead() == id) || ((*it)->getTail() == id)){
            // Delete other vertice's pointer to the edge and the lists's pointer
            if ((*it)->getWeight() < 0){
                W += (*it)->getWeight();
            }
            else{
                W -= (*it)->getWeight();
            }
            int j = (*it)->getOther(id);
            vertexmap[j]->removeEdge(*it);
            edges.erase(it++);
        }
        else{
            ++it;
        }
    }
    
    // Delete references in the map/list
    removeVertexLists(id);
}


// Print a graph
std::ostream& operator<< (std::ostream &out, const Graph &graph)
{
    // out << "Vertices: " << "\n";
    std::list<int>::const_iterator it;
    for (it = graph.vertices.begin(); it != graph.vertices.end(); ++it){
        out << *it << ", ";
    }
    out << "\n";
    //out << "Edges: " << "\n";
    //std::list<std::shared_ptr<Edge> >::const_iterator it2;
    //for (it2 = graph.edges.begin(); it2 != graph.edges.end(); ++it2){
    //    out << **it2;
    //}
    //out << "\n";
    return out;
}


// Minimum spanning tree
double Graph::MST(std::list<std::shared_ptr<Edge> > & edges) const
{
    std::map<int,bool> inTree;
    std::map<int,double> key;
    std::map<int,std::shared_ptr<Edge> > edge_keys;
    int numInTree = 0;
    double weightTree = 0;
    
    // initialize all vertices
    for (auto i : vertices){
        key[i] = INT_MAX, inTree[i] = false, edge_keys[i] = NULL;
    }
    
    // add first vertex
    key[vertices.front()] = 0;
    
    // while not all vertices are in the tree
    while (numInTree < vertices.size()){
        // find vertex with min key
        double min = INT_MAX;
        int min_v;
        
        for (auto v : vertices){
            if ((inTree[v] == false) && (key[v] < min)){
                min = key[v], min_v = v;
            }
        }
        
        // add min vertex to tree and edge
        inTree[min_v] = true;
        if (edge_keys[min_v] != NULL){
            edges.push_back(edge_keys[min_v]);
            weightTree += edge_keys[min_v]->getWeight();
        }
        
        // update key values
        for (auto e: getVertex(min_v)->getIncEdges()){
            int u = e->getOther(min_v);
            if (e->getWeight() < key[u]){
                key[u] = e->getWeight();
                edge_keys[u] = e;
            }
        }
        numInTree += 1;
        
    }
    
    return weightTree;
}


// Calculate tour length
double getTourLength(const Graph & G, const std::vector<int> & tour)
{
    int n = tour.size();
    double l = 0;
    for (int i = 0; i < tour.size()-1; i++){
        double w = shortestPath(G, tour[i], tour[i+1]);
        l += w;
        //std::cout << tour[i] << "," << tour[i+1] << "," << w << "\n";
    }
    l += shortestPath(G, tour[n-1], tour[0]);
    return l;
}

// DFS
void DFS(std::list<std::shared_ptr<Edge> > & edges, std::map<int,bool> & visited,
                   std::vector<int> & tour, int v)
{
    // Mark v as visited and add to tour
    visited[v] = true;
    tour.push_back(v);
    
    // Iterate through edges
    for (auto e:edges)
    {
        if ((e->getTail() == v) || (e->getHead() == v)){
            int u = e->getOther(v);
            if (visited[u] == false){
                DFS(edges, visited, tour, u);
            }
        }
    }
    
}

// Find ordered list of tour given a tree of edges
std::vector<int> tourList(std::list<std::shared_ptr<Edge> > & edges)
{
    // Find starting vertex and create bool map
    int v = edges.front()->getHead();
    std::map<int,bool> visited;
    for (auto e:edges){
        visited[e->getHead()] = false, visited[e->getTail()] = false;
    }
    std::vector<int> tour;
    
//    std::cout << "Constructing tour of length " << edges.size() +1 << "\n";
//    
//    for (auto e:edges){
//        std::cout << *e << ",";
//    }
//    std::cout << "\n";
    
    DFS(edges, visited, tour, v);
    return tour;
}

// Find length of shortest path from i to j in G
double shortestPath(const Graph & G, int i, int j)
{
    // Set up distance maps
    std::map<int,double> distances;
    std::map<int,bool> setDist;
    for (auto v:G.getVertices()){
        setDist[v] = false;
        distances[v] = INT_MAX;
    }
    
    // Start at i
    distances[i] = 0;
    
    // While j is not set
    while (setDist[j] == false){
        // Find min unmarked
        int minIndex = INT_MAX, u;
        for (auto v:G.getVertices()){
            if ((setDist[v] == false) && (distances[v] < minIndex)){
                u = v;
            }
        }
        
        // Set distances of u
        setDist[u] = true;
        
        // Iterate through incident edges to u
        if (distances[u] != INT_MAX){
            for (auto e : G.getVertex(u)->getIncEdges()){
                int v = e->getOther(u);
                double w = e->getWeight();
                if ((setDist[v] == false) && (distances[v] > distances[u]+w)){
                    distances[v] = distances[u]+w;
                }
            }
        }
    }
    return distances[j];
    
}

//Find the shortest path from each vertex u to all vertices in g
void all_to_all_shortest_paths(const Graph G, std::vector<int> G_Vertices, std::vector<double> verticesServicesTime, std::vector<std::vector<int>>& prev_u, std::vector<std::vector<double>>& dist_u) {
       for (auto u:G.getVertices()){
          std::vector<int> prev;//prev<> represents the previous vertex for each vertex on the shortest path 
          std::vector<double> dist;//dist<> represents the distance of each vertex from the root vertex
          shortestPath_to_DestNodes(G, u, G_Vertices, verticesServicesTime, prev, dist);
          prev_u[u] = prev;
          dist_u[u] = dist;  
       }
}


typedef std::pair<double, int> iPair;

//Using the Dijkestra algorithm to find shortest paths from the source vertex src to the destination vertices destVertices
void shortestPath_to_DestNodes(const Graph & G, int src, std::vector<int> destVertices, std::vector<double> verticesServicesTime, std::vector<int>& prev, std::vector<double>& dist) {
    // Create a priority queue to store vertices that
    // are being preprocessed. This is weird syntax in C++.
    // Refer below link for details of this syntax
    // https://www.geeksforgeeks.org/implement-min-heap-using-stl/
    std::priority_queue<iPair, std::vector<iPair>, std::greater<iPair> >
        pq;
   
    // Create a vector for distances and initialize all
    // distances as infinite (INF)
    std::vector<double> distTemp;
    //A vector to save parents of vertices in the obtained shortest paths
    std::vector<int> prevTemp;
    for (auto v:G.getVertices()){
        prevTemp.push_back(-1);
        distTemp.push_back(INT_MAX);
    }
    // Insert source itself in priority queue and initialize
    // its distance as 0.
    pq.push(std::make_pair(verticesServicesTime[src], src));
    distTemp[src] = verticesServicesTime[src];
    std::vector<int> tempdestVertices = destVertices;   
     /* Looping till priority queue becomes empty (or all
    distances are not finalized) or destVertices are visited*/
    while (!pq.empty() && !tempdestVertices.empty()) {
        // The first vertex in pair is the minimum distance
        // vertex, extract it from priority queue.
        // vertex label is stored in second of pair (it
        // has to be done this way to keep the vertices
        // sorted distance (distance must be first item
        // in pair)
        int u = pq.top().second;
        pq.pop();
        //Remove u from tempdestVertices
//        tempdestVertices.erase(remove(tempdestVertices.begin(), tempdestVertices.end(), u), tempdestVertices.end());
        std::vector<int>::iterator position = std::find(tempdestVertices.begin(), tempdestVertices.end(), u);
        if (position != tempdestVertices.end()) { // == S.end() means the element was not found
            tempdestVertices.erase(position);
        }
        // Get all adjacent of u.
						  for (auto e : G.getVertex(u)->getIncEdges()){
				  	   // Get vertex label and weight of current
          // adjacent of u.
						    int v = e->getOther(u);
						    double weight = e->getWeight() + verticesServicesTime[v];
						    //std::cout<<"Weight "<< *ei <<" is "<<EdgeWeightMap[*ei]<<"\n";
            // If there is shorted path to v through u.
            if (distTemp[v] > distTemp[u] + weight) {
                // Updating distance of v
                distTemp[v] = distTemp[u] + weight;
                pq.push(std::make_pair(distTemp[v], v));
                prevTemp[v] = u;
            }
        }
    }
 
    // Print shortest distances stored in dist[]
    //printf("Vertex Distance from Source\n");
    //dist=distTemp;
    for (int i = 0; i < distTemp.size(); ++i){
        dist.push_back(distTemp[i]);
        //prev.push_back(prevTemp[i]);
    }
    
    //prev=prevTemp;
    for (int i = 0; i < prevTemp.size(); ++i){
        prev.push_back(prevTemp[i]);
    }
}


double findTSP(std::vector<int> S_u, std::vector<std::pair<double, double>> verticesCoordinations, std::vector<double> verticesServicesTime) {
         //If the number of vertices is at most one, the cost of minimum spanning tree is 0
         if(S_u.size() == 0) return 0;
         if(S_u.size() == 1) return verticesServicesTime[S_u[0]];
         std::vector<point_2d> pointsOfTSP;
     //    std::cout<<"\nBefore setting the points\n";

         for (int i = 0; i < S_u.size(); i++) {
                  double xp = verticesCoordinations[S_u[i]].first;
                  double yp = verticesCoordinations[S_u[i]].second;
                  point_2d a = {xp, yp, i};
                  pointsOfTSP.push_back(a);
         }
//    std::cout<<"\nAfter setting the points\n";
         TSP tsp(pointsOfTSP);
         tsp.solve();
         double lengthOfTour = tsp.get_length();
//         
         double totalServieTime = 0;
         for (int v : S_u) totalServieTime += verticesServicesTime[v];
         lengthOfTour += totalServieTime;
         
         
         //2-approximation algorithm for tsp
//         double lengthOfTour = findTSP2Approx(S_u, verticesCoordinations);
//         lengthOfTour+= totalServieTime;
//         Check which one among 3/2 and 2 aprox for tsp costs less
         return lengthOfTour;// <= tempLength ? lengthOfTour : tempLength;        
}
typedef std::pair<double, int> pii;
// Function to find sum of weights of edges of the Minimum Spanning Tree.
double findTSP2Approx(std::vector<int> S_u, std::vector<std::pair<double, double>> verticesCoordinations) {
    //If the number of vertices is at most one, the cost of minimum spanning tree is 0
    if(S_u.size() <= 1) return 0;
//    if(S_u.size() == 1) return verticesServicesTime[S_u[0]];   
    
    std::vector<std::vector<double>> edgesWeights(S_u.size(), std::vector<double>(S_u.size(), 0));
    std::vector<pii> adj[S_u.size()];
    for (int i = 0; i < S_u.size(); i++)
      for (int j = i + 1; j < S_u.size(); j++){
        double dist = std::sqrt(std::pow((verticesCoordinations[S_u[j]].first - verticesCoordinations[S_u[i]].first), 2) + std::pow((verticesCoordinations[S_u[j]].second - verticesCoordinations[S_u[i]].second), 2));
        edgesWeights[i][j] = dist;
        edgesWeights[j][i] = dist;
        adj[i].push_back(std::make_pair(dist, j));
        adj[j].push_back(std::make_pair(dist, i));
      }
    
    // Create a map of keys to save minimum-cost of an edge conneted to each vertex in a mst
    std::map<int, double> key;
    // Create a map of paren to save the parent of each vertex in a mst
    std::map<int, int> parent;

    for (int u = 0; u < S_u.size(); u++) {
        key[u] = INF;
    }

        
    // Create a priority queue to store edges with their weights
    std::priority_queue<pii, std::vector<pii>, std::greater<pii>> pq;
     
    // Create a visited array to keep track of visited vertices
    std::vector<bool> visited(S_u.size(), false);
     
    // Variable to store the result (sum of edge weights)
    double res = 0;
     
    // Start with vertex 0
    pq.push(std::make_pair(0, 0));
    parent[0] = 0;

    // Perform Prim's algorithm to find the Minimum Spanning Tree
    while(!pq.empty()){
        auto p = pq.top();
        pq.pop();
         
        double wt = p.first;  // Weight of the edge
        int u = p.second;  // Vertex connected to the edge
         
        if(visited[u] == true) {
            continue;  // Skip if the vertex is already visited
        }
         
        res += wt;  // Add the edge weight to the result
        visited[u] = true;  // Mark the vertex as visited
         
        // Explore the adjacent vertices
        for(auto v : adj[u]) {
            // v[0] represents the vertex and v[1] represents the edge weight
            if(visited[v.second] == false) {
                pq.push(std::make_pair(v.first, v.second));  // Add the adjacent edge to the priority queue
                //Set the parent and key of v.second
                if (key[v.second] > edgesWeights[u][v.second]) {
                    key[v.second] = edgesWeights[u][v.second];
                    parent[v.second] = u;
                }
            }
        }
    }
    
    
//    std::cout<<"\S_u: ";
    std::list<std::shared_ptr<Edge> > mstEdges;
    for (int u = 0; u < S_u.size(); u++) {
//        std::cout<<" "<<u;
//        std::cout<<"("<<parent[u]<<", "<<u<<", "<<edgesWeights[parent[u]][u]<<") ";
        std::shared_ptr<Edge> e = std::make_shared<Edge>(parent[u], u, edgesWeights[parent[u]][u]);
        mstEdges.push_back(e);
    }
    std::vector<int> mstTour = tourList(mstEdges);
    
    int n = mstTour.size();
    double lengthOfTour = 0;
    for (int i = 0; i < mstTour.size()-1; i++){
//        std::cout<<"("<<mstTour[i]<<", "<<mstTour[i+1]<<", "<<edgesWeights[mstTour[i]][mstTour[i+1]]<<") ";
        double w = edgesWeights[mstTour[i]][mstTour[i+1]]; //shortestPath(G, tour[i], tour[i+1]);
        lengthOfTour += w;
        //std::cout << tour[i] << "," << tour[i+1] << "," << w << "\n";
    }
//    std::cout<<"("<<mstTour[n-1]<<", "<<mstTour[0]<<", "<<edgesWeights[mstTour[n-1]][mstTour[0]]<<")\n";
    lengthOfTour += edgesWeights[mstTour[n-1]][mstTour[0]];//shortestPath(G, tour[n-1], tour[0]);
    //Adding the cost of nodes to the mst
//    for (int v : S_u) lengthOfTour = lengthOfTour + verticesServicesTime[v];
    return lengthOfTour;  // Return the sum of edge weights of the Minimum Spanning Tree



//    return res;
//    Graph subG;
//          
//    for(int i = 0; i < S_u.size(); i++){
//        subG.addVertex(i);
//    }
//    
//   for(int head = 0; head < S_u.size(); head++)
//      for(int tail = head + 1; tail < S_u.size(); tail++)
//          subG.addEdge(head, tail, edgesWeights[S_u[head]][S_u[tail]]);
//
//   //first check if the cost of minimum spanning tree is at most half of the budget    
//   std::list<std::shared_ptr<Edge>> mst;
//   double mst_w = subG.MST(mst); 
//   std::vector<int> mstTour = tourList(mst);
////   double l = getTourLength(subG, mstTour);
//   
//   int n = mstTour.size();
//   double l = 0;
//   for (int i = 0; i < mstTour.size()-1; i++){
////        std::cout<<"("<<mstTour[i]<<", "<<mstTour[i+1]<<", "<<edgesWeights[mstTour[i]][mstTour[i+1]]<<") ";
//       double w = edgesWeights[S_u[mstTour[i]]][S_u[mstTour[i+1]]]; //shortestPath(G, tour[i], tour[i+1]);
//       l += w;
//       //std::cout << tour[i] << "," << tour[i+1] << "," << w << "\n";
//   }
////    std::cout<<"("<<mstTour[n-1]<<", "<<mstTour[0]<<", "<<edgesWeights[mstTour[n-1]][mstTour[0]]<<")\n";
//   l += edgesWeights[S_u[mstTour[n-1]]][S_u[mstTour[0]]];//shortestPath(G, tour[n-1], tour[0]);    
}



void dpForSize(std::vector<int> verticesOfAlg, double capacity, std::vector<std::pair<double, double>> verticesCoordinationsOfAlg, std::vector<double> verticesPrizesOfAlg, std::vector<double> verticesSizesOfAlg, std::vector<double> verticesServicesTimeOfAlg, double& totalSize, double& totalPrize, double& totalCost, std::vector<int>& verticesOfAlgRespectingCapacity) {

    totalSize = 0;
    totalPrize = 0;
    totalCost = 0;
//    std::cout<<"\nVerticesOfAlg: ";
//    for (int v : verticesOfAlg) std::cout<< v <<" ";
//    std::cout<<"\n";
    //Select the items with the maximum prize respecting the capacity constraint by dynamic programming as the same technique for the Knapsack problem
    std::vector<std::vector<double>> DP(verticesOfAlg.size() + 1, std::vector<double>(capacity + 1));
        
    for (int i = 0; i <= verticesOfAlg.size(); i++) {
        for (double w = 0; w <= capacity; w++) {
            if (i == 0 || w == 0)
                DP[i][w] = 0;
            else if (verticesSizesOfAlg[verticesOfAlg[i - 1]] <= w)
                DP[i][w] = std::max(verticesPrizesOfAlg[verticesOfAlg[i - 1]] + DP[i - 1][w - verticesSizesOfAlg[verticesOfAlg[i - 1]]], DP[i - 1][w]);
            else
                DP[i][w] = DP[i - 1][w];
        }
    }
    
    totalPrize =  DP[verticesOfAlg.size()][capacity];
    double res = totalPrize;
    double w = capacity;
    std::vector<int> vecTempVertices;
    double tempSize = 0;
    for (int i = verticesOfAlg.size(); i > 0 && res > 0; i--) {
        if (res == DP[i - 1][w])
            continue;    
        else {
 
            // This item verticesOfAlg[i - 1] is included.
            vecTempVertices.push_back(verticesOfAlg[i - 1]);
            // Since this weight is included its
            // value is deducted
            res = res - verticesPrizesOfAlg[verticesOfAlg[i - 1]];
            w = w - verticesSizesOfAlg[verticesOfAlg[i - 1]];
            tempSize += verticesSizesOfAlg[verticesOfAlg[i - 1]];
        }
    }
    
    totalSize = tempSize;
    
    //Add the depot to vecTempVertices if it's not selected
    if(!(std::find(vecTempVertices.begin(), vecTempVertices.end(), 0) != vecTempVertices.end()))
         vecTempVertices.insert(vecTempVertices.begin(), 0);
//    totalCost = findTSP2Approx(vecTempVertices, edgesWeightsOfAlg, verticesCoordinationsOfAlg, verticesServicesTimeOfAlg);
    totalCost = findTSP(vecTempVertices, verticesCoordinationsOfAlg, verticesServicesTimeOfAlg);
    verticesOfAlgRespectingCapacity = vecTempVertices;
}



//Add more vertices to the solution if it is possible
void addMoreVertices(std::vector<ratioPrizeCost> vecRatioPrizeCost, std::vector<int>& verticesRespectingCapacity, double& size_x, double& prize_x, double& cost_x, double budget, double capacity, std::vector<int> vectTemp, std::vector<std::pair<double, double>> verticesCoordinations, std::vector<double> verticesPrizes, std::vector<double> verticesSizes, std::vector<double> verticesServicesTime, std::string algName) {
   double tempCost_x = 0;
   double tempSize_x = size_x;
   std::vector<int> tempVerticesRespectingCapacity = verticesRespectingCapacity;
   //Add the depot to vecTempVertices if it's not selected
//   if(!(std::find(vectTemp.begin(), vectTemp.end(), 0) != vectTemp.end()))
//         vectTemp.insert(vectTemp.begin(), 0);
   if (cost_x < budget && size_x < capacity) {
     for (auto vrpc : vecRatioPrizeCost) {
         if (cost_x < budget && size_x < capacity) {
            if(cost_x + verticesServicesTime[vrpc.id] < budget && size_x + verticesSizes[vrpc.id] <= capacity) {
                 bool flag;
                 if(algName == "SqrB") flag = true;
                 else flag = !(std::find(vectTemp.begin(), vectTemp.end(), vrpc.id) != vectTemp.end());
                 if (flag) {
                    tempVerticesRespectingCapacity.push_back(vrpc.id);
                    tempCost_x = findTSP(tempVerticesRespectingCapacity, verticesCoordinations, verticesServicesTime);//findTSP2Approx(tempVerticesRespectingCapacity, edgesWeights, verticesCoordinations, verticesServicesTime);
                    tempSize_x += verticesSizes[vrpc.id];
                    if (tempCost_x <= budget && tempSize_x <= capacity) {
                        prize_x += verticesPrizes[vrpc.id];
                        size_x = tempSize_x;
                        cost_x = tempCost_x;
                        verticesRespectingCapacity.push_back(vrpc.id);
                    }
                    tempSize_x = size_x;
                    tempVerticesRespectingCapacity = verticesRespectingCapacity;
               }
            }
         }
         else break;
     }
   }
}

std::vector<std::vector<double>> coordinationsToEdges(std::vector<std::pair<double, double>> verticesCoordinations) {
     //Compute the edge cost based on vertices coordinations
        std::vector<std::vector<double>> edgesWeights(verticesCoordinations.size(), std::vector<double>(verticesCoordinations.size()));
        //Compute edge costs based on the the vertices coordinations
        for(int i = 0; i < verticesCoordinations.size(); i++){
          for(int j = i + 1; j < verticesCoordinations.size(); j++){
              double dist = std::sqrt(std::pow((verticesCoordinations[j].first - verticesCoordinations[i].first), 2) + std::pow((verticesCoordinations[j].second - verticesCoordinations[i].second), 2));
              edgesWeights[i][j] = dist;
              edgesWeights[j][i] = dist;
          }
        }
        return edgesWeights;
}

//build the graph by vertices prizes and edges weights
Graph buildGraph(std::vector<double> verticesPrizes, std::vector<std::vector<double>> edgesWeights, std::vector<int>& verticesOfG) {
           Graph G;
           for (int i = 0; i < verticesPrizes.size(); i++) {
               G.addVertex(i, verticesPrizes[i]);
               verticesOfG.push_back(i);
           } 
           
           //Add edges to the graph
           for(int head = 0; head < edgesWeights.size(); head++)
              for(int tail = head + 1; tail < edgesWeights[head].size(); tail++)
                  G.addEdge(head, tail, edgesWeights[head][tail]);
                     
           return G; 
}
double subGraphMST(std::vector<int> vecTempVertices, std::vector<std::vector<double>> edgesWeightsOfAlg){
   //Find the minimum spanning tree covering verticesOfAlgRespectingCapacity
   //First creat a subgraph on the vertices in verticesOfAlgRespectingCapacity
   Graph subG;
   for(int i = 0; i < vecTempVertices.size(); i++){
      subG.addVertex(i);
   }
   for(int head = 0; head < vecTempVertices.size(); head++)
      for(int tail = head + 1; tail < vecTempVertices.size(); tail++)
          subG.addEdge(head, tail, edgesWeightsOfAlg[vecTempVertices[head]][vecTempVertices[tail]]);
   //Now we find the minimum spanning tree of subG
   std::list<std::shared_ptr<Edge>> mstSubG;
   return subG.MST(mstSubG);
}

