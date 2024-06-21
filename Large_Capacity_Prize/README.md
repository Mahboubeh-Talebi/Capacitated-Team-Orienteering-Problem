CTOPAlgoSimulation
=================

C++ code to build and run evaluation of the Capacitated Team Orienteering Problem

The algorithms were implemented and evaluated in order to write the paper which title is
"Improved Algorithms for Capacitated Budgeted Prize-Collecting Vehicle Routing". 

# Run this project

Note that you need to run this project on large instances given by Tarantilis et al. “The capacitated team orienteering problem: A bi-level filter-and-fan method.”

If you want to run this project, first run the make file by running make in the command line
then run the following command line
    ./main  -i (or -d )input file  -f output file  -s bool(bool=0 or 1) -C capacity -B budget -l set_x_parameter(>=1) -v number (>= 1) -e violation_factor(>=0) -z number of threads(>=0)
    , where
        
    -i specifies the input file (i.e., the graph), which should be selected from set of small instances. Instead of -i we can use -d to read all the files of a directory
    
    -f is the output file in which we save the results
      
    -l enables us to set parameter x in SBAA algorithm in the paper and it should be at least 1
    
    -v enables us to set parameter y in SBAA algorithm in the paper and it should be at least 1
    
     -s specifies whether we want to solve the 4-AA algorithm in the paper, -s 0 says no, -s 1 says yes
     
     -e (default=0) is the factor by which we violate the budget in the algorithms and it should be at least 0
     
     -z specifies the number of threads by which we solve SBAA and GA in the paper
     
 *Note that running ./main will give you a help to run the algorithms 

# Evaluate the algorithms

This project contains the implementation of 3 heuristic algorithms 

# CTOP Heuristic Algorithms
The 3 approximation algorithms are
- SBAA, a modification of the algorithm from “Budgeted prize-collecting traveling salesman and minimum spanning tree problems” by Paul et al.
- GA, a greedy algorithm
- GRA, modification of the algorithm from “Randomized algorithm for informative path planning with budget constraints” by Arora and A. Scherer.


