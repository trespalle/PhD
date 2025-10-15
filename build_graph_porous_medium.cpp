#include <vector>
#include <fstream>
#include <cmath>
#include <string>
#include <random>
#include <iostream>
#include <cstdlib>
#include "Rand.h"
#include "graph_theory.h"
#include <algorithm>


RNG Rand;

bool read_flows_and_junctions(const std::string& filename, std::vector<double>& flows, std::vector<std::vector<std::size_t>>& juncs, std::vector<std::size_t>& junc_type);
void build_graph(Graph& G, std::vector<double>& flows, std::vector<std::vector<std::size_t>>& juncs, std::vector<std::size_t>& junc_type);


int main(int argc, char** argv){


    Graph G;
    std::vector<double> flows;
    std::vector<std::vector<std::size_t>> juncs;
    std::vector<std::size_t> junc_types;

    if (read_flows_and_junctions("flows_and_junctions.txt", flows, juncs, junc_types )) 
    {
        std::cout << "Flows read: " << flows.size() << "\n";
        if (!flows.empty())
            std::cout << "First flow: " << flows[0] << "\n";

        std::cout << "Triplets read: " << juncs.size() << "\n";

        if (!juncs.empty()) 
            std::cout << "First triplet: " << juncs[0][0] << " " << juncs[0][1] << " " << juncs[0][2] << "\n"; 
    }

    

    build_graph(G, flows, juncs, junc_types);

    double n_type1 = 0.0;
    for(auto& link : G.links) if((1.0 - link.weight) < 1.0e-5) n_type1 += 1.0;
    n_type1 /= G.links.size();
    std::cout<<n_type1<<std::endl;

    G.save_links("links_real.txt");

    return 0;

}

bool read_flows_and_junctions(const std::string& filename, std::vector<double>& flows, std::vector<std::vector<std::size_t>>& juncs, std::vector<std::size_t>& junc_type) 
{
    std::ifstream fin(filename);
    if (!fin) 
    {
        std::cerr << "Error: could not open file " << filename << "\n";
        return false;
    }

    int nflows;
    if (!(fin >> nflows)) 
    {
        std::cerr << "Error: could not read number of flows\n";
        return false;
    }

    flows.resize(nflows);
    for (int i = 0; i < nflows; i++) fin >> flows[i];


    juncs.clear();
    junc_type.clear();
    std::size_t type, a, b, c;
    while (fin >> type >> a >> b >> c) 
    {
        junc_type.push_back(type);
        juncs.push_back({a, b, c});  // each triplet as a vector<int>
    }

    return true;
}

void build_graph(Graph& G, std::vector<double>& flows, std::vector<std::vector<std::size_t>>& juncs, std::vector<std::size_t>& junc_type)
{
    std::size_t N = juncs.size();
    G.nodes.resize(N);
    
    for(std::size_t i = 0; i < juncs.size(); ++i)
    {
        std::size_t jmax = 0;
        double qmax = flows[juncs[i][0]];
        for(std::size_t j = 1; j < juncs[i].size(); ++j) 
        {
            if(std::abs(flows[juncs[i][j]]) > std::abs(qmax))
            {
                qmax = flows[juncs[i][j]];
                jmax = j;
            }
        }

        std::swap(juncs[i][jmax], juncs[i].back());

    
        if(junc_type[i] == 2) // type 2 junction
        {
            double weight = Rand();
            G.add_link(juncs[i][2], juncs[i][0], weight);
            G.add_link(juncs[i][2], juncs[i][1], 1.0 - weight);
        }

        else //type 1 junction
        {
            G.add_link(juncs[i][0], juncs[i][2], 1.0);
            G.add_link(juncs[i][1], juncs[i][2], 1.0);
        }

    }
}







