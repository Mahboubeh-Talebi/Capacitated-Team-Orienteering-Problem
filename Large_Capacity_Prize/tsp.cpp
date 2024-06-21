#include "tsp.h"

#include <utility>


//Constructor
TSP::TSP(vector<point_2d> aPointList) : points(std::move(aPointList)) {
    n = points.size();

    graph = new double *[n];
    for (int i = 0; i < n; i++) {
        graph[i] = new double[n];
        for (int j = 0; j < n; j++) {
            graph[i][j] = 0.;
        }
    }

    adj_list = new vector<int>[n];
}


//Destructor
TSP::~TSP() {
    for (int i = 0; i < n; i++) {
        delete[] graph[i];
    }

    delete[] graph;
    delete[] adj_list;
}


void TSP::solve() {
    fill_matrix();
    find_MST();

    perfect_matching();

    // Loop through each index and find the shortest path
    double best = numeric_limits<double>::max();
    int bestIndex;
    for (long t = 0; t < n; t++) {
        double result = find_best_path(t);
        if (result < best) {
            bestIndex = t;
            best = result;
        }
    }

    //Create path for best tour
    euler_tour(bestIndex, circuit);
    make_hamiltonian(circuit, path_length);
}

double TSP::get_length() {
    return path_length;
}

double TSP::get_distance(struct point_2d c1, struct point_2d c2) {
    double dx = pow(c1.x - c2.x, 2);
    double dy = pow(c1.y - c2.y, 2);
    double dz = 0; // = pow(c1.id - c2.id, 2);

    return sqrt(dx + dy + dz);
}

void TSP::fill_matrix() {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            graph[i][j] = graph[j][i] = get_distance(points[i], points[j]);
        }
    }
}

/******************************************************************************
  This function uses Prim's algorithm to determine a minimum spanning tree on
    the graph
******************************************************************************/
typedef pair<double, int> pii;
void TSP::find_MST() {
    auto *key = new double[n];
    bool *included = new bool[n];
    int *parent = new int[n];

    for (int i = 0; i < n; i++) {

        // set each key to infinity
        key[i] = numeric_limits<double>::max();

        // node node yet included in MST
        included[i] = false;

    }

    // root of MST has distance of 0 and no parent
    key[0] = 0.;
    parent[0] = -1;
    
    // Create an adjacency list representation of the graph
    vector<pii> adj[n];
     
    // Fill the adjacency list with edges and their weights
    for (int i = 0; i < n; i++) 
       for (int j = 0; j < n; j++) {
        adj[i].push_back(std::make_pair(graph[i][j], j));
        adj[j].push_back(std::make_pair(graph[i][j], i));
    }
     
    // Create a priority queue to store edges with their weights
    priority_queue<pii, vector<pii>, greater<pii>> pq;
     
    // Create a visited array to keep track of visited vertices
    vector<bool> visited(n, false);
     
    // Start with vertex 0
    pq.push(std::make_pair(0, 0));
     
    // Perform Prim's algorithm to find the Minimum Spanning Tree
    while(!pq.empty()){
        auto p = pq.top();
        pq.pop();
         
        int wt = p.first;  // Weight of the edge
        int u = p.second;  // Vertex connected to the edge
         
        if(visited[u] == true){
            continue;  // Skip if the vertex is already visited
        }
         
//        res += wt;  // Add the edge weight to the result
        visited[u] = true;  // Mark the vertex as visited
         
        // Explore the adjacent vertices
        for(auto v : adj[u]){
            // v[0] represents the vertex and v[1] represents the edge weight
            if(visited[v.second] == false){
                pq.push(std::make_pair(v.first, v.second));  // Add the adjacent edge to the priority queue
                if (key[v.second] > graph[u][v.second]) {
                    key[v.second] = graph[u][v.second];
                    parent[v.second] = u;
                }
            }
        }
    }
    
//    for (int i = 0; i < n - 1; i++) {
//
//        // find closes vertex not already in tree
//        int k = get_min_index(key, included);
//
//        // set included to true for this vertex
//        included[k] = true;
//
//        // examine each unexamined vertex adjacent to most recently added
//        for (int j = 0; j < n; j++) {
//
//            // node exists, is unexamined, and graph[k][j] less than previous
//            // key for u
//            if (graph[k][j] && !included[j] && graph[k][j] < key[j]) {
//
//                // update parent
//                parent[j] = k;
//
//                // update key
//                key[j] = graph[k][j];
//            }
//        }
//    }

    delete[] key;
    delete[] included;

    // construct a tree by forming adjacency matrices
    for (int i = 0; i < n; i++) {
        int j = parent[i];
        if (j != -1) {
            adj_list[i].push_back(j);
            adj_list[j].push_back(i);
        }
    }

    delete[] parent;
}

/******************************************************************************
  find the index of the closest unexamined node
******************************************************************************/

int TSP::get_min_index(const double key[], const bool mst[]) {
    // initialize min and min_index
    double min = numeric_limits<double>::max();
    int min_index;

    // iterate through each vertex
    for (int i = 0; i < n; i++) {

        // if vertex hasn't been visited and has a smaller key than min
        if (!mst[i] && key[i] < min) {

            // reassign min and min_index to the values from this node
            min = key[i];
            min_index = i;
        }
    }
    return min_index;
}

/******************************************************************************
  find all vertices of odd degree in the MST. Store them in an subgraph O
******************************************************************************/

void TSP::find_odds() {
    for (int i = 0; i < n; i++) {

        // if degree of vertex i is odd
        if ((adj_list[i].size() % 2) != 0) {

            // push vertex to odds, which is a representation of subgraph O
            odds.push_back(i);
        }
    }
}

void TSP::perfect_matching() {
    /************************************************************************************
     find a perfect matching M in the subgraph O using greedy algorithm but not minimum
    *************************************************************************************/
    int closest;
    double length; //int d;
    vector<int>::iterator tmp, first;

    // Find nodes with odd degrees in double to get subgraph O
    find_odds();

    // for each odd node
    while (!odds.empty()) {
        first = odds.begin();
        auto it = odds.begin() + 1;
        auto end = odds.end();
        length = numeric_limits<int>::max();
        for (; it != end; ++it) {
            // if this node is closer than the current closest, update closest and length
            if (graph[*first][*it] < length) {
                length = graph[*first][*it];
                closest = *it;
                tmp = it;
            }
        } // two nodes are matched, end of list reached
        adj_list[*first].push_back(closest);
        adj_list[closest].push_back(*first);
        odds.erase(tmp);
        odds.erase(first);
    }
}

//find an euler circuit
void TSP::euler_tour(int start, vector<int> &path) {
    //Create copy of adj. list
    vector<vector<int> > tempList;

    for (int i = 0; i < n; i++) {
        tempList.push_back(adj_list[i]);
    }

    stack<int> stack;
    int pos = start;
    path.push_back(start);
    while (!stack.empty() || !tempList[pos].empty()) {
        //Current node has no neighbors
        if (tempList[pos].empty()) {
            //add to circuit
            path.push_back(pos);
            //remove last vertex from stack and set it to current
            pos = stack.top();
            stack.pop();
        }
            //If current node has neighbors
        else {
            //Add vertex to stack
            stack.push(pos);
            //Take a neighbor
            int neighbor = tempList[pos].back();
            //Remove edge between neighbor and current vertex
            tempList[pos].pop_back();
            for (unsigned int i = 0; i < tempList[neighbor].size(); i++) {
                if (tempList[neighbor][i] == pos) {
                    tempList[neighbor].erase(tempList[neighbor].begin() + i);
                }
            }
            //Set neighbor as current vertex
            pos = neighbor;
        }
    }
    path.push_back(pos);
}

//Make euler tour Hamiltonian
void TSP::make_hamiltonian(vector<int> &path, double &pathCost) {
    //remove visited nodes from Euler tour
    bool *visited = new bool[n];

    for (int i = 0; i < n; i++) {
        visited[i] = false;
    }

    pathCost = 0.;

    int root = path.front();
    auto cur = path.begin();
    auto iter = path.begin() + 1;
    visited[root] = true;

    //iterate through circuit
    while (iter != path.end()) {
        if (!visited[*iter]) {
            pathCost += graph[*cur][*iter];
            cur = iter;
            visited[*cur] = true;
            iter = cur + 1;
        } else {
            iter = path.erase(iter);
        }
    }

    delete[] visited;

    //Add distance to root
    pathCost += graph[*cur][*iter];
}

double TSP::find_best_path(int start) {
    vector<int> path;
    euler_tour(start, path);

    double length;
    make_hamiltonian(path, length);

    return length;
}

void TSP::print_result() {
    for (int &it: circuit) {
        cout << it << endl;
    }
}

void TSP::print_path() {
    cout << endl;
    for (auto it = circuit.begin(); it != circuit.end() - 1; ++it) {
        cout << *it << "->" << *(it + 1) << " ";
        cout << graph[*it][*(it + 1)] << endl;
    }
    // Print the last edge separately
    auto lastVertex = *(circuit.end() - 1);
    auto firstVertex = circuit.front();
    cout << lastVertex << "->" << firstVertex << " ";
    cout << graph[lastVertex][firstVertex] << endl;

    cout << "\nLength: " << path_length << endl << endl;
}

vector<int> TSP::get_path_id() {
    vector<int> tmp;

    for (auto p: circuit) {
        tmp.push_back(p);
    }
    // Also the first one
    tmp.push_back(circuit[0]);

    return tmp;
}

vector<point_2d> TSP::get_path() {
    vector<point_2d> path;

    for (auto i: circuit) {
        path.push_back(points[i]);
    }
    // Also the first one
    path.push_back(points[circuit[0]]);

    return path;
}

void TSP::print_adj_list() {
    for (int i = 0; i < n; i++) {
        cout << i << ": "; //print which vertex's edge list follows
        for (int j: adj_list[i]) {
            cout << j << " "; //print each item in edge list
        }
        cout << endl;
    }
}
