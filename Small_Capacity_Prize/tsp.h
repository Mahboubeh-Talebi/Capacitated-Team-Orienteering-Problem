#ifndef TSP_H
#define TSP_H

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <queue>
#include <stack>
#include <string>
#include <cstdio>
#include <vector>
#include <limits>

using namespace std;

struct point_2d {
    double x;
    double y;
    int id;
};

//---------------------------------------------------------------------------
class TSP {
public:

protected:
    // List of odd nodes
    vector<int> odds;

    //Adjacency list
    vector<int> *adj_list;

    // Number of points
    int n;

    //Shortest path length
    double path_length{};

    //euler circuit
    vector<int> circuit;

    // n x n, pairwise distances between points
    double **graph;

    // Point list
    vector<point_2d> points;

    // -
    void find_odds();

    //Find perfect matching
    void perfect_matching();

    //Find Euler tour
    void euler_tour(int start, vector<int> &path);

    //Find Hamiltonian path
    void make_hamiltonian(vector<int> &path, double &pathCost);

    double get_distance(struct point_2d c1, struct point_2d c2);

    //T get_distance(struct point_2d c1, struct point_2d c2);

    void find_MST();

    int get_min_index(const double key[], const bool mst[]);

    void fill_matrix();

    double find_best_path(int start);

public:
    TSP(vector<point_2d> aPointList);

    ~TSP();

    void solve();

    vector<int> get_path_id();

    vector<point_2d> get_path();

    double get_length();

    void print_result();

    void print_path();

    void print_adj_list();
};

#endif
