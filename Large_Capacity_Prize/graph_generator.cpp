#include "graph_generator.h"

typedef boost::small_world_iterator<boost::minstd_rand, UndirectedGraph> SWGen;
typedef boost::sorted_erdos_renyi_iterator<boost::minstd_rand, UndirectedGraph> ERGen;
typedef boost::plod_iterator<boost::minstd_rand, UndirectedGraph> SFGen;
//Constructor
Graph_Generator::Graph_Generator(std::string graphT, int val_n, int val_k, double val_a, double val_b, double val_p, int val_minPrize, int val_maxPrize, double val_minWeight, double val_maxWeight, std::string val_outputfiledirectory){
    //Set variables
    n=val_n;
    k=val_k;
    p=val_p;
    a=val_a;
    b=val_b;
    minWeight=val_minWeight;
    maxWeight=val_maxWeight;
    minPrize=val_minPrize;
    maxPrize=val_maxPrize;
    std::map<std::string, int> graphType;
    const int ER = 1, SW=2, SF= 3;
    const std::string SFstr = "sf", ERstr = "er", SWstr="sw";
    graphType[SFstr] = SF;
    graphType[ERstr] = ER;
    graphType[SWstr] = SW;
    //A std::vector showing that which vertices are reachable from r
    std::vector<bool> RVisited(n, false);
    int iter=0;//A threshold indicating the number of times to generate an Erdos-Reny graph in which n-1 vertices that are reachable from r
    //Generating a graph gTemp with n vertices reachable from r
    UndirectedGraph gTemp;
    while(count(RVisited.begin(), RVisited.end(), true)<n && iter<pow(n, 1.5)){
         srand(time(0));
         boost::minstd_rand gen(rand());
         switch(graphType[graphT]){
               case SF:{
                    // Create graph with n nodes
                    UndirectedGraph gSF(SFGen(gen, n, a, b), SFGen(), n);
                    gTemp=gSF;
                    break;
               }
               case SW:{
                    // Create graph with n nodes, and edges with probability p
                    UndirectedGraph gSW(SWGen(gen, n, k, p), SWGen(), n);
                    gTemp=gSW;
                    break;
               }
               case ER:{
                    // Create graph with n nodes, and edges with probability p
                    UndirectedGraph gER(ERGen(gen, n, p), ERGen(), n);
                    gTemp=gER;
                    break;
               }
               default:
                     std::cout << "\n" << "Error: Select " << SFstr << " or " << ERstr << " or "<<SWstr<<" as graph type" << "\n";
         }
         //Computing a std::vector showing that which vertices are reachable from the root in gTemp.
         RVisited= Visit_BFS(gTemp, 0);
         ++iter;
         std::cout<<"iteration "<<iter<<"\n";
         //We make sure that the last vertex in gTemp is reachable from the root. This is due to the graph representation
   }
   //Check either the number of reachable vertices from r is n or vertex n is reachable from the root, we check the latter condition due to the representation of graph in boost library
   if(count(RVisited.begin(), RVisited.end(), true)<n && !RVisited[n-1]){
     std::cout<<"Either the number of reachable vertices from the root is less than n  or it's not possible to run the algorithm on the generated graph due to reachability of vertices.\n";
     std::cout<<"Please set the parameters for generating the graph in a different way and try again.\n";
   }
   else{//The number of reachable vertices from the root is reasonable
   std::cout<<"Numbe of reachable vertices from the root "<<count(RVisited.begin(), RVisited.end(), true)<<"\n";
   boost::graph_traits<UndirectedGraph>::vertex_iterator vi, vi_end;
   for (boost::tie(vi, vi_end) = vertices(gTemp); vi != vi_end; ++vi){
       if(!RVisited[*vi]){
          remove_vertex(*vi, gTemp);
          RVisited.erase(RVisited.begin()+*vi);
          --vi;
       }
   }
   //Create a grach containing reachable vertices from vertex 0 in gTemp
   typedef boost::graph_traits<UndirectedGraph>::edge_descriptor edge_t;
   std::vector <edge_t> vecEdgeG;
   //Adding reachable vertices from root and the corresponding edges  from gTemp to val_g
   typename boost::graph_traits<UndirectedGraph>::edge_iterator e, e_end;
   for (boost::tie(e, e_end) = edges(gTemp); e!= e_end; ++e){
      if(RVisited[source(*e, gTemp)] && RVisited[target(*e, gTemp)] && source(*e, gTemp)!=target(*e, gTemp) && target(*e, gTemp)!=0 && !(find(vecEdgeG.begin(), vecEdgeG.end(), *e)!=vecEdgeG.end())){
          add_edge(source(*e, gTemp), target(*e, gTemp), val_g);
          vecEdgeG.push_back(*e);//Indicating that the edge (source(*e, g), target(*e, g)) is visited
      }
   }
  std::cout<<"The number of edges in the generated graph is: "<<boost::num_edges(val_g)<<"\n";
  std::cout<<"The number of vertices in the generated graph is: "<<boost::num_vertices(val_g)<<"\n";
  //Set the terminals

  //Select t reachable random vertices from root as terminals
  genRandomPrizes();
  //Generate random weight (set the weight file)
  genRandomWeight();
  //Name the output file based on the machine time
  time_t t = time(NULL);
  tm *ptm = localtime(&t);
  std::stringstream ss;
  ss << ptm->tm_year+1900<<"-"<<ptm->tm_mon+1<<"-"<< ptm->tm_mday <<"-"<< ptm->tm_hour+1 <<"-"<< ptm->tm_min+1 <<"-"<< ptm->tm_sec+1;
  std::cout << "Generating " << std::flush;
  switch(graphType[graphT]) {
      case SF:{
        std::cout << "Scale-Free\n";
        filename_t="SF_"+ss.str();
        outputfile.open((val_outputfiledirectory+filename_t+".txt").c_str(), std::ios::app);
        break;
        }
      case SW:{
        std::cout << "Small-World\n";
        filename_t="SW_"+ss.str();
        outputfile.open((val_outputfiledirectory+filename_t+".txt").c_str(), std::ios::app);
        break;}
      case ER:{
        std::cout << "Erdos-Renyi\n";
        filename_t="ER_"+ss.str();
        outputfile.open((val_outputfiledirectory+filename_t+".txt").c_str(), std::ios::app);
        break;
        }
      default:
        std::cout << "\n" << "Error: Select " << SFstr << " or " << ERstr << " as graph type" << "\n";
    }
    //Set the edges and terminals in the outputfile
   setInputFile();
   outputfile.close();
  }//End-Else
}

//Select t reachable vertices fom the root vetex as terminals randomly in case we have at least t reachable vertice in g
void Graph_Generator::genRandomPrizes(){
    std::random_device rdw; // obtain a random number from hardware
    std::mt19937 genPrize(rdw()); // seed the generator
    std::uniform_int_distribution<> distrW(minPrize, maxPrize); // define the range 
    for(int i = 0; i < boost::num_vertices(val_g); i++) 
        val_VerticesPrizes.push_back(distrW(genPrize));
}

//Generate random weights for vertices
void Graph_Generator::genRandomWeight(){
       std::random_device rdw; // obtain a random number from hardware
       std::mt19937 genWeight(rdw()); // seed the generator
       std::uniform_int_distribution<> distrW(minWeight, maxWeight); // define the range
       for(int i=0; i<boost::num_edges(val_g); ++i)
            val_EdgesWeights.push_back(distrW(genWeight));
}

//Write the graph (its edges, weights of vertices and terminals) in the file
void Graph_Generator::setInputFile(){
      outputfile<<"\n";
      outputfile<<"SECTION Comment\n Name    "<<filename_t<<"\n"<<"END\n\n";
      outputfile<<"SECTION Graph"<<"\n"<<"Nodes: "<< boost::num_vertices(val_g)<<"\n";
      //outputfile<<"SECTION Prizes\n";
     for(int i=0; i<boost::num_vertices(val_g); ++i){
         outputfile<<"N "<< i <<" "<<val_VerticesPrizes[i]<<"\n";
     }
     outputfile<<"\n\n";
     outputfile<<"Edges "<<boost::num_edges(val_g)<<"\n";
     typename boost::graph_traits<UndirectedGraph>::edge_iterator e, e_end;
     int count_edges=0;
     for (boost::tie(e, e_end) = boost::edges(val_g); e!= e_end; ++e){
        outputfile<<"E "<< source(*e, val_g)<<" "<< target(*e, val_g)<<" "<<val_EdgesWeights[count_edges]<<"\n";
        ++count_edges;
     }
     outputfile<<"\n\n";
     //outputfile<<"EOF";
}

 // Return reachable vertices from a source node in agraph by BFS algorithm
std::vector<bool> Graph_Generator::Visit_BFS(UndirectedGraph g, int src) {
  //A vector to save the reachability of vertices from src in g
  std::vector<bool> vecVisited(num_vertices(g), false);
  vecVisited[src] = true;
  //A queue to set the order of visitting the vertices
  std::list<int> queue;
  queue.push_back(src);

  std::list<int>::iterator i;

  while (!queue.empty()) {
    int currVertex = queue.front();
    //cout << "Visited " << currVertex << " ";
    queue.pop_front();

    typename boost::graph_traits <UndirectedGraph >::out_edge_iterator ei, ei_end;
    //Check on out_edges of  vertex
				for (tie(ei, ei_end) = out_edges(currVertex, g); ei != ei_end; ++ei) {
      int adjVertex = target (*ei, g );
      if (!vecVisited[adjVertex]) {
        vecVisited[adjVertex] = true;
        queue.push_back(adjVertex);
      }
    }
  }
  
  return vecVisited;
}

