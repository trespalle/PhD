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
#include <filesystem>
#include <map>
#include <functional>

struct Params
{
    std::string filename;
    int nreals = 1;
    int nbins = 50;
    int randomize_flag = 0;
};

RNG Rand;

void fix_disorder(Graph& G);
Params read_params(int argc, char** argv);
void initialize_flows(Graph& G);

// Recursive helper function for Depth-First Search (DFS)
void dfs_sum_products( std::size_t current_node_idx, double current_product, const Graph& G_original, const std::vector<bool>& is_big_tube, std::map<size_t, double>& target_weights, std::vector<bool>& visited) 
{
    // Mark the current node as visited for this branch of the search
    visited[current_node_idx] = true;

    // Explore all links leaving the current node
    for (std::size_t out_link_idx : G_original.nodes[current_node_idx].outlnks) 
    {
        const auto& out_link = G_original.links[out_link_idx];
        
        // Calculate the product of weights for the next step
        double next_product = current_product * out_link.weight;

        if (is_big_tube[out_link_idx]) 
        {
            // Success! We found a path to another "big link".
            // Add the product of this path to the total weight for this destination.
            target_weights[out_link_idx] += next_product;
            
            // DO NOT continue the search down this branch.
        } 
        else 
        {
            // It's a "small" link, so we continue the search if the next node hasn't been visited.
            size_t next_node_idx = out_link.out;
            if (!visited[next_node_idx]) 
                dfs_sum_products(next_node_idx, next_product, G_original, is_big_tube, target_weights, visited);
        }
    }

    // Unmark the current node when backtracking.
    // This is VITAL to allow other paths to use this node.
    visited[current_node_idx] = false;
}



Graph create_big_tubes_graph(const Graph& G_original)
{
    Graph G_big;
    if (G_original.links.empty()) 
    {
        std::cout << "Warning: Original graph has no links. Returning empty big_graph." << std::endl;
        return G_big;
    }

    // Identify nodes and create the map  
    std::vector<bool> is_big_tube(G_original.links.size(), false);
    for(std::size_t i = 0; i < G_original.links.size(); ++i) 
        if(G_original.nodes[G_original.links[i].in].innbrs.size() == 2)  is_big_tube[i] = true;
    
    // Links connected to original inlets are now "big"
    for(std::size_t inlet_node_idx : G_original.inlet)
        for(std::size_t link_idx : G_original.nodes[inlet_node_idx].outlnks) is_big_tube[link_idx] = true;

    for (std::size_t outlet_node_idx : G_original.outlet) 
        for (std::size_t link_idx : G_original.nodes[outlet_node_idx].inlnks) is_big_tube[link_idx] = true;

    std::map<size_t, size_t> original_to_big_map;
    std::vector<size_t> big_tube_indices;
    for(std::size_t i = 0; i < G_original.links.size(); ++i) 
    {
        if (is_big_tube[i]) 
        {
            // Create a new node and get a reference to it
            Node& new_big_node = G_big.nodes.emplace_back();
            
            // Get the original link that this new node represents
            const auto& original_link = G_original.links[i];
            
            // Set the new node's mass 
            new_big_node.mass = G_original.nodes[original_link.in].mass*original_link.weight;

            // Map the original link index to the new node index
            std::size_t new_node_index = G_big.nodes.size() - 1;
            original_to_big_map[i] = new_node_index;
            big_tube_indices.push_back(i);
        }
    }
    
    //Find links and calculate weights with DFS
    for(std::size_t start_link_idx : big_tube_indices) 
    {
        std::size_t start_node_in_big_graph = original_to_big_map.at(start_link_idx);
        std::size_t dfs_start_node = G_original.links[start_link_idx].out;

        // Map to store the results: <destination_link_index> -> <accumulated_weight>
        std::map<std::size_t, double> target_weights;
        std::vector<bool> visited(G_original.nodes.size(), false);

        // Start the recursive search
        dfs_sum_products(dfs_start_node, 1.0, G_original, is_big_tube, target_weights, visited);

        // Once the DFS has explored all paths, we add the found links to G_big
        for(const auto& pair : target_weights) 
        {
            std::size_t end_link_idx = pair.first;
            double total_weight = pair.second;
            std::size_t end_node_in_big_graph = original_to_big_map.at(end_link_idx);
            
            G_big.add_link(start_node_in_big_graph, end_node_in_big_graph, total_weight);
        }
    }

    G_big.update_boundaries();

    return G_big;
}

void create_output_directory(const std::string& dir_path)
{
    try 
    {
        // Attempt to create the directory (and any parent directories needed)
        std::filesystem::create_directories(dir_path);
    } 
    catch (const std::exception& e) 
    {
        // If it fails, print an error message
        std::cerr << "Error creating directory '" << dir_path << "': " << e.what() << std::endl;
        // Exit the program, as the output directory is essential
        std::exit(1); 
    }
}

int main(int argc, char** argv){

    // Read parameters (number of histogram bins & number of realizations of the network:
    Params p = read_params(argc, argv);
    std::cout << "nbins = " << p.nbins << std::endl << "nreals = " << p.nreals << std::endl;

    // Create output directory:
    const std::string output_dir = "output/";
    create_output_directory(output_dir);

    // We read the graph from its list of links:
    Graph G(p.filename);
    
    // We compute the stationary flows for p.nreals realizations of the disorder (weights):
    std::vector<double> mass_history;
    std::vector<double> flow_history, big_flow_history, small_flow_history; 
    std::vector<double> histogram(p.nbins);
    double minflow, maxflow;
    
    for(std::size_t e = 0; e < p.nreals; ++e)
    {
        
        initialize_flows(G); 
        if(p.randomize_flag == 1) fix_disorder(G);
        G.compute_stationary_flows();
        
        //We collect flow-rate data at pores (junctions):
        for(auto& pore : G.nodes) if(!pore.is_boundary) mass_history.push_back(pore.mass);
        //We collect flow-rate data at tubes:
        for(auto& tube : G.links) 
        {
            if((!G.nodes[tube.in].is_boundary)&&(!G.nodes[tube.out].is_boundary))
            {
                double flow = G.nodes[tube.in].mass*tube.weight;
                flow_history.push_back(flow);
                if(G.nodes[tube.in].innbrs.size() == 2) big_flow_history.push_back(flow);
                else small_flow_history.push_back(flow);
                
            }
                
        }         
    }

    // Compute correspoinding histograms
    histogram = myfun::buildHistogram(mass_history, p.nbins, minflow, maxflow);
    myfun::guardaHistograma(histogram, minflow, maxflow, p.nbins, output_dir + "pore_flow_dist.dat");
    histogram = myfun::buildHistogram(flow_history, p.nbins, minflow, maxflow);
    myfun::guardaHistograma(histogram, minflow, maxflow, p.nbins, output_dir + "tube_flow_dist.dat");
    histogram = myfun::buildHistogram(big_flow_history, p.nbins, minflow, maxflow);
    myfun::guardaHistograma(histogram, minflow, maxflow, p.nbins, output_dir + "big_tube_flow_dist.dat");
    histogram = myfun::buildHistogram(small_flow_history, p.nbins, minflow, maxflow);
    myfun::guardaHistograma(histogram, minflow, maxflow, p.nbins, output_dir + "small_tube_flow_dist.dat");

    // Checking mass conservation: we display total inlet mass vs. total outlet mass
    double inlet_mass, outlet_mass;
    inlet_mass = outlet_mass = 0.0;
    for(auto& v : G.inlet) inlet_mass += G.nodes[v].mass;
    for(auto& v : G.outlet) outlet_mass += G.nodes[v].mass;
    std::cout << "Total inlet flow = " << inlet_mass << ", Total outlet flow = " << outlet_mass << std::endl;
    std::cout << "Graph has " << G.inlet.size() << " inlets." << std::endl;

    
    // BIG TUBE ANALYISIS:
    std::cout << "\nBIG TUBES GRAPH ANALYSIS" << std::endl;

    //Create the big-tube graph
    Graph G_big = create_big_tubes_graph(G);
    std::cout << "Big graph has " << G_big.nodes.size() << " nodes and " << G_big.links.size() << " links." << std::endl;
    std::cout << "Big graph has " << G_big.inlet.size() << " inlets." << std::endl;
    // Save the structure of the new graph
    G_big.save_links(output_dir + "big_graph_links.dat");




    /* // Again, (just to check consistency) we compute the stationary flows for p.nreals realizations of the disorder (weights):
    mass_history.clear();
    flow_history.clear();
    big_flow_history.clear(); 
    small_flow_history.clear(); 
    
    for(std::size_t e = 0; e < p.nreals; ++e)
    {
        
        //initialize_flows(G_big); 
        G_big.compute_stationary_flows();
        
        //We collect flow-rate data at tubes:
        for(auto& big_tube : G_big.nodes) if(!big_tube.is_boundary) 
        {
            flow_history.push_back(big_tube.mass);
            big_flow_history.push_back(big_tube.mass);
        }
        for(auto& small_tube : G_big.links) 
        {
            double flow = G_big.nodes[small_tube.in].mass*small_tube.weight;
            flow_history.push_back(flow);
            small_flow_history.push_back(flow);
        }

        
    }

    // Compute correspoinding histograms
    histogram = myfun::buildHistogram(flow_history, p.nbins, minflow, maxflow);
    myfun::guardaHistograma(histogram, minflow, maxflow, p.nbins, output_dir + "tube_flow_dist_2.dat");
    histogram = myfun::buildHistogram(big_flow_history, p.nbins, minflow, maxflow);
    myfun::guardaHistograma(histogram, minflow, maxflow, p.nbins, output_dir + "big_tube_flow_dist_2.dat");
    histogram = myfun::buildHistogram(small_flow_history, p.nbins, minflow, maxflow);
    myfun::guardaHistograma(histogram, minflow, maxflow, p.nbins, output_dir + "small_tube_flow_dist_2.dat");

    // Checking mass conservation: we display total inlet mass vs. total outlet mass
    inlet_mass = outlet_mass = 0.0;
    for(auto& v : G.inlet) inlet_mass += G.nodes[v].mass;
    for(auto& v : G.outlet) outlet_mass += G.nodes[v].mass;
    std::cout << "Total inlet flow = " << inlet_mass << ", Total outlet flow = " << outlet_mass << std::endl;
 */
    
    std::cout << "BIG TUBES GRAPH ANALYSIS FINISHED" << std::endl;


    std::cout << "\nAll output files have been saved to the '" << output_dir << "' directory." << std::endl;

    return 0;

}


void fix_disorder(Graph& G)
{
    std::size_t N = G.nodes.size();
    
    //out-weights:
    for(std::size_t i = 0; i < N; i++)
    {
        //if(G.nodes[i].outnbrs.size() > 2) std::cout << "Hola_mundo" << std::endl;
        if(G.nodes[i].outnbrs.size() == 2)
        {
            G.nodes[i].outwgs[0] = Rand();
            G.nodes[i].outwgs[1] = 1.0 - G.nodes[i].outwgs[0];
            G.links[G.nodes[i].outlnks[0]].weight =  G.nodes[i].outwgs[0];
            G.links[G.nodes[i].outlnks[1]].weight =  G.nodes[i].outwgs[1];
        }
    }

    //in-weights:
    for (std::size_t i = 0; i < N; ++i) for (std::size_t j = 0; j < G.nodes[i].outnbrs.size(); ++j)
    {
        auto outnbr  = G.nodes[i].outnbrs[j];
        auto outwg   = G.nodes[i].outwgs[j];
        auto& innbrs = G.nodes[outnbr].innbrs;
        auto& inwgs = G.nodes[outnbr].inwgs;
        // we find the position of i in innbrs
        auto pos = std::find(innbrs.begin(), innbrs.end(), i);
        if (pos != innbrs.end()) 
        {
            std::size_t k = std::size_t(pos - innbrs.begin());
            if (k < inwgs.size()) inwgs[k] = outwg;           
        } 
    
}

}

Params read_params(int argc, char** argv)
{
   Params p;

   if(argc < 2)
   {
    std::cerr << "Error: you must provide a filename as the first argument." << std::endl;
    std::cerr << "Usage: " << argv[0] << " <filename> [n_realizations] [n_bins] [0 for quenched weights, 1 for annealed]" << std::endl;
    std::exit(1);
   }

   p.filename = argv[1];
   if(argc > 2) p.nreals = std::atoi(argv[2]);
   if(argc > 3) p.nbins = std::atoi(argv[3]);
   if(argc > 4) p.randomize_flag = std::atoi(argv[4]);

   return p;
}

void initialize_flows( Graph& G)
{
    for(auto& v : G.nodes) v.mass = 0.0;
    double sum = 0.0;
    for(auto& i : G.inlet) 
    {
        double mass = Rand();
        G.nodes[i].mass = mass;
        sum += mass;
    }
    
    for(auto& i : G.inlet) G.nodes[i].mass /= sum;
}

