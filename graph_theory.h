#ifndef GRAPH_THEORY_H
#define GRAPH_THEORY_H

#include <vector>
#include <cstddef>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>


struct Node
{
    double mass = 0.0;
    int state = 0;
    bool is_boundary = false;
    std::vector<double> coordinates;
    std::vector<std::size_t> outnbrs, innbrs, outlnks, inlnks; ///list of outnbrs and innbrs labels
    std::vector<double> outwgs, inwgs;
};

struct Link
{
    double weight;
    std::size_t in, out;
};

struct Graph
{
    Graph() = default; /// if the constructor I will define below is not used, then the default constructor is used 
    std::vector<Node> nodes;
    std::vector<Link> links;
    std::vector<std::size_t> inlet, outlet;

private:
    std::size_t position(std::size_t guy, const std::vector<std::size_t>& people)
    {
        std::size_t pos = 0;
        for(pos = 0; pos < people.size(); ++pos) if(people[pos] == guy) return pos;
        std::cout << "We are sorry, but the element was not found in the vector." << std::endl;
        std::exit(1);
    }

    void erase_duplicates(std::vector<std::size_t>& my_vector)
    {
        std::sort(my_vector.begin(), my_vector.end());
        auto last = std::unique(my_vector.begin(), my_vector.end()); // "unique" puts the unique elements of the vector first and then the rest
        // last is an iterator (the C++ evolution of the pointer, a pointer with gadgets) that point to the begining of that "rest" in the new vector
        my_vector.erase(last, my_vector.end()); ///daraa
    }

    


public:
    void add_link(std::size_t i, std::size_t j, double weight=1.0)
    {
        if(i >= nodes.size()) nodes.resize(i+1);
        if(j >= nodes.size()) nodes.resize(j+1);

        nodes[i].outnbrs.push_back(j);
        nodes[i].outwgs.push_back(weight);
        nodes[i].outlnks.push_back(links.size());
        nodes[j].innbrs.push_back(i);
        nodes[j].inwgs.push_back(weight);
        nodes[j].inlnks.push_back(links.size());
        
        Link newlink;
        newlink.weight = weight;
        newlink.in = i;
        newlink.out = j;
        links.push_back(newlink);
    }

    void erase_link(std::size_t in, std::size_t out)
    {
        ///First, we delete "out" from the list of outnbrs of "in"
        std::size_t pos = position(out, nodes[in].outnbrs);
        nodes[in].outnbrs.erase(nodes[in].outnbrs.begin() + pos);
        nodes[in].outwgs.erase(nodes[in].outwgs.begin() + pos);
        std::size_t link_label = nodes[in].outlnks[pos];
        nodes[in].outlnks.erase(nodes[in].outlnks.begin() + pos);

        ///Then, we delete "in" from the list of innbrs of "out":
        pos = position(link_label, nodes[out].inlnks);
        nodes[out].innbrs.erase(nodes[out].innbrs.begin() + pos);
        nodes[out].inwgs.erase(nodes[out].inwgs.begin() + pos);
        nodes[out].inlnks.erase(nodes[out].inlnks.begin() + pos);

        ///Finally, we delete the link from the list of links, and we shift the labels of all the other links accordingly:
        links.erase(links.begin() + link_label);
        for(std::size_t i = 0; i < nodes.size(); ++i) 
        {
            for(std::size_t j = 0; j < nodes[i].outlnks.size(); ++j) 
                if(nodes[i].outlnks[j] > link_label) -- nodes[i].outlnks[j];
            for(std::size_t j = 0; j < nodes[i].inlnks.size(); ++j) 
                if(nodes[i].inlnks[j] > link_label) -- nodes[i].inlnks[j];
        }

    }

    void erase_link(std::size_t link_label)
    {
        ///First, we delete "out" from the list of outnbrs of "in"
        if(link_label >= links.size()) 
        {
            std::cout << "You almost kill us all." << std::endl;
            std::exit(2);
        }
        std::size_t in, out, pos;
        in = links[link_label].in;
        out = links[link_label].out;
        pos = position(link_label, nodes[in].outlnks);
        nodes[in].outnbrs.erase(nodes[in].outnbrs.begin() + pos);
        nodes[in].outwgs.erase(nodes[in].outwgs.begin() + pos);
        nodes[in].outlnks.erase(nodes[in].outlnks.begin() + pos);

        ///Then, we delete "in" from the list of innbrs of "out":
        pos = position(link_label, nodes[out].inlnks);
        nodes[out].innbrs.erase(nodes[out].innbrs.begin() + pos);
        nodes[out].inwgs.erase(nodes[out].inwgs.begin() + pos);
        nodes[out].inlnks.erase(nodes[out].inlnks.begin() + pos);

        ///Finally, we delete the link from the list of links, and we shift the labels of all the other links accordingly:
        links.erase(links.begin() + link_label);
        for(std::size_t i = 0; i < nodes.size(); ++i) 
        {
            for(std::size_t j = 0; j < nodes[i].outlnks.size(); ++j) 
                if(nodes[i].outlnks[j] > link_label) -- nodes[i].outlnks[j];
            for(std::size_t j = 0; j < nodes[i].inlnks.size(); ++j) 
                if(nodes[i].inlnks[j] > link_label) -- nodes[i].inlnks[j];
        }

    }

    void erase_node(std::size_t label)
    {
        //First, we erase all the links in which it participates
        auto outlnks = nodes[label].outlnks;
        std::sort(outlnks.rbegin(), outlnks.rend()); //we order the vector in decreasing order, so that sequential elimination of links is not problematic
        for(std::size_t i = 0; i < outlnks.size(); ++i) erase_link(outlnks[i]);
        auto inlnks = nodes[label].inlnks;
        std::sort(inlnks.rbegin(), inlnks.rend());
        for(std::size_t i = 0; i < inlnks.size(); ++i) erase_link(inlnks[i]);

        // Then, we erase the node from the list of nodes
        nodes.erase(nodes.begin() + label);

        //Finally, we shift the labels of all the "bigger nodes" one unit down:
        for(std::size_t i = 0; i < nodes.size(); ++i)
        {
            for(std::size_t j = 0; j < nodes[i].outnbrs.size(); ++j) 
                if(nodes[i].outnbrs[j] > label) -- nodes[i].outnbrs[j];

            for(std::size_t j = 0; j < nodes[i].innbrs.size(); ++j) 
                if(nodes[i].innbrs[j] > label) -- nodes[i].innbrs[j];   
        }

        for(std::size_t i = 0; i < links.size(); ++i)
        {
            if(links[i].in > label) -- links[i].in;
            if(links[i].out > label) -- links[i].out;
            
        }
        
    }

    void check_and_renormalize_weights()
    {
        std::cout << "Starting K1L check and renormalization..." << std::endl;
        bool any_node_failed = false;
        std::size_t N = nodes.size();

        //Loop 1: check, renormalize outwgs, and update links ---
        for(std::size_t i = 0; i < N; ++i)
        {
            // We only care about nodes with outgoing links
            if(nodes[i].outwgs.empty()) continue; 

            // Calculate the sum of outgoing weights
            double sum = 0.0;
            for(double w : nodes[i].outwgs) sum += w;

            
            // Case 1: Sum is 0 (problematic)
            if (std::abs(sum) < 1e-10) 
            {
                if(!any_node_failed) 
                { 
                    // Print header only the first time
                    std::cerr << "--- K1L WARNINGS (Critical Errors) ---" << std::endl;
                    any_node_failed = true;
                }
                std::cerr << "  ERROR Node " << i << ": Has " << nodes[i].outwgs.size() << " outgoing links, but their weight sum is 0." << " Cannot renormalize! Flows will be incorrect." << std::endl;
                continue; // We cannot divide by zero, skip this node
            }
            
            // Case 2: Sum is not 1.0 (standard violation)
            if(std::abs(sum - 1.0) > 1e-9)
            {
                if(!any_node_failed) 
                { 
                    // Print header only the first time
                    std::cerr << "--- K1L WARNINGS (Corrections) ---" << std::endl;
                    any_node_failed = true;
                }
                std::cerr << "  Node " << i << ": Weight sum = " << sum << ". Forcing renormalization." << std::endl;
            }

            // Force renormalization (always, if sum != 0)
            // This corrects violations from Case 2 and also ensures
            // nodes that already summed to 1 are unaffected.
            for(std::size_t j = 0; j < nodes[i].outwgs.size(); ++j)
            {
                // 1. Renormalize the outwg of node i
                nodes[i].outwgs[j] /= sum;
                
                // 2. Update the weight in the global links list
                std::size_t link_idx = nodes[i].outlnks[j];
                links[link_idx].weight = nodes[i].outwgs[j];
            }
        }

        // Loop 2: update in-weights (copied from your logic in fix_disorder) 
        // This loop ensures that the 'inwgs' of neighbors reflect
        // the renormalized weights we just calculated.
        for (std::size_t i = 0; i < N; ++i) 
        {
            for (std::size_t j = 0; j < nodes[i].outnbrs.size(); ++j)
            {
                auto outnbr  = nodes[i].outnbrs[j];
                auto outwg   = nodes[i].outwgs[j]; // The already renormalized weight
                auto& innbrs = nodes[outnbr].innbrs;
                auto& inwgs = nodes[outnbr].inwgs;
                
                // Find 'i' in the 'innbrs' list of 'outnbr'
                auto pos = std::find(innbrs.begin(), innbrs.end(), i);
                if (pos != innbrs.end()) 
                {
                    std::size_t k = std::size_t(pos - innbrs.begin());
                    if (k < inwgs.size()) inwgs[k] = outwg;           
                } 
                // No 'else' is needed; if not found, it's a graph
                // consistency error that 'add_link' should prevent.
            }
        }

        // Final Message 
        if (any_node_failed) 
        {
            std::cerr << "------------------------------------------" << std::endl;
            std::cout << "K1L check finished. Violations were found and corrected (or reported)." << std::endl;
        } 
        else std::cout << "K1L check passed. All nodes comply." << std::endl;
    
    }


    void update_boundaries()
    {
        // Reset the entire boundary state first
        inlet.clear();
        outlet.clear();
        for (auto& node : nodes) node.is_boundary = false;

        // Recalculate inlets and outlets in a single pass
        for(std::size_t i = 0; i < nodes.size(); ++i) 
        {
            // Check for inlet condition
            if (nodes[i].innbrs.empty()) 
            {
                inlet.push_back(i);
                nodes[i].is_boundary = true;
            }
            // Check for outlet condition
            if (nodes[i].outnbrs.empty()) 
            {
                outlet.push_back(i);
                nodes[i].is_boundary = true;
            }
        }
    }

    void read_coordinates(const std::string& filename)
    {
        std::ifstream fin(filename);
        if (!fin) throw std::runtime_error("Error: could not open coordinate file: " + filename);
    
        std::string line;
        std::size_t node_id;
        double coord_val;
        int line_num = 0;

        while (std::getline(fin, line)) 
        {
            line_num++;
            std::stringstream ss(line);

            // Read Node ID
            if (!(ss >> node_id)) throw std::runtime_error("Error reading node ID from line " + std::to_string(line_num) + " in " + filename);
            if (node_id >= nodes.size()) throw std::runtime_error("Error: Node ID " + std::to_string(node_id) + " from coordinate file is out of bounds (max is " + std::to_string(nodes.size() - 1) + ").");
        
            // Clear previous coordinates and read new ones
            nodes[node_id].coordinates.clear();
            while (ss >> coord_val) nodes[node_id].coordinates.push_back(coord_val);
            
            if (nodes[node_id].coordinates.empty()) throw std::runtime_error("Error: No coordinates found for node ID " + std::to_string(node_id) + " on line " + std::to_string(line_num) + " in " + filename);
            
        }
        fin.close();
        std::cout << "Successfully read coordinates for " << nodes.size() << " nodes from " << filename << std::endl;
    }
    void check_mass_conservation()
    {
        std::cout << "Checking mass conservation for non-outlet nodes..." << std::endl;
        bool conservation_failed = false;

        for (std::size_t i = 0; i < nodes.size(); ++i)
        {
            // Solo nos importan los nodos que no son 'outlets'
            // (ya que los outlets por definición no tienen salida)
            if (nodes[i].outnbrs.empty()) continue; 

            double weight_sum = 0.0;
            for (double w : nodes[i].outwgs)
            {
                weight_sum += w;
            }

            // Comprobamos si la suma se desvía significativamente de 1.0
            if (std::abs(weight_sum - 1.0) > 1e-9)
            {
                std::cerr << "WARNING: Mass not conserved at node " << i << ". "
                        << "Sum of output weights is: " << weight_sum 
                        << " (Difference: " << (weight_sum - 1.0) << ")" << std::endl;
                conservation_failed = true;
                std::cout << "Outnbrs: " << std::endl;
                for( auto& v : nodes[i].outnbrs) std::cout << v << " ";
                std::cout << std::endl;

            }
        }

        if (!conservation_failed)
        {
            std::cout << "Mass conservation check passed for all nodes." << std::endl;
        }
    }


    void compute_stationary_flows()
    {
        std::vector<std::size_t> active_nodes;
        // First, we define some functions:

        auto initialize_active_nodes = [&]()
        {
            if(inlet.empty()) (*this).update_boundaries();
            if(inlet.empty()) 
            {
                std::cerr << "The graph has no inlet nodes. This is problematic." << std::endl;
                std::exit(1);
            }
            
            for(std::size_t i = 0; i < inlet.size(); ++i) for(std::size_t j = 0; j < nodes[inlet[i]].outnbrs.size(); ++j) 
                active_nodes.push_back(nodes[inlet[i]].outnbrs[j]);

            erase_duplicates(active_nodes);
        };

        auto update_active_nodes = [&](const std::vector<double>& new_flows)
        {
            std::vector<std::size_t> new_active_nodes;

            for(std::size_t i = 0; i < active_nodes.size(); ++i) 
            {
                if(std::abs(new_flows[i] - nodes[active_nodes[i]].mass) > 1e-10)
                {
                    for(std::size_t j = 0; j < nodes[active_nodes[i]].outnbrs.size(); ++j) 
                        new_active_nodes.push_back(nodes[active_nodes[i]].outnbrs[j]);
                }
            }

            erase_duplicates(new_active_nodes);
            for(std::size_t i = 0; i < active_nodes.size(); ++i) nodes[active_nodes[i]].mass = new_flows[i];
            active_nodes = new_active_nodes;
        };


        auto sweep = [&]()
        {
            std::vector<double> new_flows(active_nodes.size(), 0.0);
            for(std::size_t i = 0; i < active_nodes.size(); ++i)
            {
                std::size_t label = active_nodes[i];

                for(std::size_t j = 0; j < nodes[label].innbrs.size(); ++j)
                {
                    std::size_t innbr = nodes[label].innbrs[j];
                    new_flows[i] += nodes[label].inwgs[j]*nodes[innbr].mass;
                }    
            }
            update_active_nodes(new_flows);
        };

        // Now, the core of the function itself:

        initialize_active_nodes();
        while(!active_nodes.empty()) sweep();
        
    }

    Graph(const std::string& file_name, const std::string& coords_filename = "") ///this is the constructor
    {
        //We read the list of links from a file:
        std::ifstream fin(file_name);
        if(!fin)
        {
            std::cerr << "Error: could not open file.";
            std::exit(5);
        }
        
        std::size_t i, j;
        double weight;

        while(fin >> i >> j >> weight) add_link(i, j, weight);
        fin.close();

        check_and_renormalize_weights();

        // Read node coordinates (if provided)
        if (!coords_filename.empty()) 
        {
            try 
            {
                read_coordinates(coords_filename);
            } 
            catch (const std::runtime_error& e) 
            {
                 // Re-throw or handle more gracefully if needed
                 throw std::runtime_error("Error during coordinate reading: " + std::string(e.what()));
            }
        }

        update_boundaries();
        compute_stationary_flows();
    }

    void print_links()
    {
        for(std::size_t i = 0; i < links.size(); ++i) 
            std::cout << links[i].in << " " << links[i].out << " " << links[i].weight << std::endl;
        std::cout << std::endl;
    }

    void save_links(const std::string& name_file)
    {
        std::ofstream fout(name_file);
        if(!fout)
        {
            std::cerr << "Unable to open output file.";
            std::exit(1);
        }
        for(std::size_t i = 0; i < links.size(); ++i) 
            fout << links[i].in << " " << links[i].out << " " << links[i].weight << std::endl;

    }

    void print_outnbrs()
    {
        for(std::size_t i = 0; i < nodes.size(); ++i) 
        {
            std::cout << i << " :";
            for(std::size_t j = 0; j < nodes[i].outnbrs.size(); ++j)
                std::cout << " " << nodes[i].outnbrs[j];
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    void save_outnbrs(const std::string& name_file)
    {
        std::ofstream fout(name_file);
        if(!fout)
        {
            std::cerr << "Unable to open output file.";
            std::exit(1);
        }
        for(std::size_t i = 0; i < nodes.size(); ++i) 
        {
            fout << i << " ";
            for(std::size_t j = 0; j < nodes[i].outnbrs.size(); ++j)
                fout << " " << nodes[i].outnbrs[j];
            fout << std::endl;
        }
    }

    void print_innbrs()
    {
        for(std::size_t i = 0; i < nodes.size(); ++i) 
        {
            std::cout << i << " :";
            for(std::size_t j = 0; j < nodes[i].innbrs.size(); ++j)
                std::cout << " " << nodes[i].innbrs[j];
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    void save_innbrs(const std::string& name_file)
    {
        std::ofstream fout(name_file);
        if(!fout)
        {
            std::cerr << "Unable to open output file.";
            std::exit(1);
        }
        for(std::size_t i = 0; i < nodes.size(); ++i) 
        {
            fout << i << " ";
            for(std::size_t j = 0; j < nodes[i].innbrs.size(); ++j)
                fout << " " << nodes[i].innbrs[j];
            fout << std::endl;
        }
    }

};

#endif