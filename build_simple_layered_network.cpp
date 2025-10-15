#include "myfunctions.h"
#include <vector>
#include <fstream>
#include <cmath>
#include <string>
#include <random>
#include <iostream>
#include <cstdlib>
#include "Rand.h"
#include "graph_theory.h"

struct Params
{
    int width = 10;
    int length = 100;
};



Params read_params(int argc, char** argv);
Graph simple_layered_network(std::size_t width, std::size_t length);


int main(int argc, char** argv){

    
    Params p = read_params(argc, argv);

    std::cout << "width = " << p.width << std::endl << "length = " << p.length << std::endl;
    std::cout << "njuncs = " << p.width*p.length << std::endl;
    
    
    //We create the network:
    Graph G = simple_layered_network(p.width, p.length);
    G.save_links("links_simple_layered.txt");
    


    return 0;

}

Params read_params(int argc, char** argv)
{
   Params p;
   
   if(argc > 1) p.width = std::atoi(argv[1]);
   if(argc > 2) p.length = std::atoi(argv[2]);
   return p;
}

Graph simple_layered_network(std::size_t width, std::size_t length)
{
    //This function creates a simple layered network in which each junction gives its flow to the 
    //closest two junctions of the next layer.

    Graph G;
    std::size_t N = length*width;
    
    G.nodes.resize(N);
   
    for(std::size_t i = 0; i < (N - width); ++i)
    {
        //"bulk" cases:
        G.add_link(i, i + width);
        G.add_link(i, i + width + 1);
    }
    
    //we correct edge cases and impose PBC: 
    for(std::size_t i = (width-1); i < (N-width); i += width)
    {
        G.erase_link(i, i + width + 1);
        G.add_link(i, i + 1);
    }
    
    if (G.nodes.size() == N + 1) G.erase_node(N);


    return G;
             
}


