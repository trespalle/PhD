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


struct Node
{
    double mass = 0.0;
    int state = 0;
    bool is_boundary = false;
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
                if(std::abs(new_flows[i] - nodes[active_nodes[i]].mass) > 1e-05)
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

    Graph(const std::string& file_name) ///this is the constructor
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