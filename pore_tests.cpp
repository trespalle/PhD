#include <iostream>
#include <vector>
#include <string>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <filesystem> 
#include "pore_networks.h"
#include "myfunctions.h" 

namespace fs = std::filesystem;

int main(int argc, char** argv) 
{
    std::cout << "==========================================" << std::endl;
    std::cout << "      PORE NETWORK CONSISTENCY CHECK      " << std::endl;
    std::cout << "==========================================" << std::endl;

    // 0. PREPARE OUTPUT
    std::string output_dir = "check_output/";
    try {
        if (!fs::exists(output_dir)) {
            fs::create_directories(output_dir);
        }
    } catch (const std::exception& e) {
        std::cerr << "Error creating directory: " << e.what() << std::endl;
        return 1;
    }

    // 1. DEFINE FILE SCHEMAS & LOAD
    // We read Weight, Width and FlowRate from links (columns 3, 4, 5)
    std::vector<LinkData> link_fmt = { 
        LinkData::Weight, 
        LinkData::HalfWidth, 
        LinkData::FlowRate 
    };

    // We read Pressure from nodes (column 3, after ID, X, Y)
    std::vector<NodeData> node_fmt = { 
        NodeData::Pressure 
    };

    std::string link_file = "junction_links.txt";
    std::string coord_file = "junction_coordinates.txt";

    std::cout << "Loading network from " << link_file << " and " << coord_file << "..." << std::endl;
    
    // Load topology and experimental data (Q and P)
    // '2' indicates 2D coordinates (skip X, Y before reading Pressure)
    PoreNetwork PN(link_file, coord_file, link_fmt, node_fmt, 2); 

    std::cout << "\n[Topology Loaded]" << std::endl;
    std::cout << "  Nodes:   " << PN.nodes.size() << std::endl;
    std::cout << "  Links:   " << PN.links.size() << std::endl;
    std::cout << "  Inlets:  " << PN.inlet.size() << std::endl;
    std::cout << "  Outlets: " << PN.outlet.size() << std::endl;


    // 2. INITIAL CHECKS (Pressure & Mass Balance)
    // -------------------------------------------------------
    
    // 2a. Check Pressure Uniformity at Boundaries
    auto check_pressure_stats = [&](const std::vector<std::size_t>& nodes, std::string label) {
        if(nodes.empty()) return;
        double min_p = 1e200, max_p = -1e200, sum_p = 0.0;
        for(auto id : nodes) {
            double p = PN.pressures[id];
            if(p < min_p) min_p = p;
            if(p > max_p) max_p = p;
            sum_p += p;
        }
        double avg = sum_p / nodes.size();
        std::cout << "  " << label << " Pressure: " 
                  << "Min=" << min_p << ", Max=" << max_p << ", Avg=" << avg << " Pa" << std::endl;
    };

    std::cout << "\n[Boundary Pressure Check (Experimental Data)]" << std::endl;
    check_pressure_stats(PN.inlet, "Inlet");
    check_pressure_stats(PN.outlet, "Outlet");

    // 2b. Check Global Flow Balance
    // Calculate total flow entering/leaving the network boundaries
    double total_inlet_Q = 0.0;
    for (std::size_t id : PN.inlet) {
        for (std::size_t link_idx : PN.nodes[id].outlnks) {
            // We sum absolute values of flow on boundary links
            total_inlet_Q += std::abs(PN.flows[link_idx]);
        }
    }

    double total_outlet_Q = 0.0;
    for (std::size_t id : PN.outlet) {
        for (std::size_t link_idx : PN.nodes[id].inlnks) {
            total_outlet_Q += std::abs(PN.flows[link_idx]);
        }
    }
    
    // We will use the average of In/Out as the target for the re-simulation
    double avg_Q_target = (total_inlet_Q + total_outlet_Q) / 2.0;
    
    std::cout << "\n[Global Mass Balance Check (Experimental Data)]" << std::endl;
    std::cout << "  Total Inlet Flow:  " << total_inlet_Q << " m^3/s" << std::endl;
    std::cout << "  Total Outlet Flow: " << total_outlet_Q << " m^3/s" << std::endl;
    std::cout << "  -> Target Q for Simulation: " << avg_Q_target << " m^3/s" << std::endl;


    // 3. INVERSE PROBLEM: Calculate Effective Conductances
    // -------------------------------------------------------
    // Use loaded Q and P to calibrate g_geom = (Q*mu)/dP for every tube.
    // This populates PN.geometric_conductances based on the experimental reality.
    std::cout << "\n[Inverse Problem: Calibrating Conductances]" << std::endl;
    PN.init_physics_from_experimental_data();


    // 4. SAVE INITIAL HISTOGRAMS (Experimental Data)
    // -------------------------------------------------------
    std::cout << "\n[Saving Experimental Histograms]" << std::endl;
    
    std::size_t nbins = 50;
    double min_val, max_val, delta; 

    // Helper lambda for histogram saving using 'myfun' namespace
    auto save_histogram = [&](const std::vector<double>& data, std::string name) {
        if(data.empty()) return;
        // Convert to absolute values for log-plotting consistency
        std::vector<double> abs_data;
        abs_data.reserve(data.size());
        for(double v : data) abs_data.push_back(std::abs(v));

        // Build histogram (Note: buildHistogram sorts the input vector, so we pass the copy 'abs_data')
        auto h = myfun::buildHistogram(abs_data, nbins, min_val, max_val, delta, 0.99); 
        myfun::guardaHistograma(h, min_val, delta, output_dir + name);
        std::cout << "  Saved " << name << std::endl;
    };

    save_histogram(PN.geometric_conductances, "hist_conductance.dat");
    save_histogram(PN.pressure_drops, "hist_pressure_drop_exp.dat");
    save_histogram(PN.flows, "hist_flow_rate_exp.dat");


    // 5. FORWARD PROBLEM: Re-simulate Flow
    // -------------------------------------------------------
    // Now we forget the read pressures and flows.
    // We calculate a new equilibrium field based ONLY on the conductances we just calibrated.
    // We impose the 'avg_Q_target' we calculated in step 2b.
    
    std::cout << "\n[Forward Problem: Re-simulating Flow Field]" << std::endl;
    
    // Use FixedFlowRate (Neumann at inlet) to match the experimental condition
    PN.solve_flow_field(BoundaryCondition::FixedFlowRate, avg_Q_target);


    // 6. SAVE FINAL HISTOGRAMS (Simulated Data)
    // -------------------------------------------------------
    std::cout << "\n[Saving Simulated Histograms]" << std::endl;

    // Conductances haven't changed, so we only save the new dynamic variables
    save_histogram(PN.pressure_drops, "hist_pressure_drop_sim.dat");
    save_histogram(PN.flows, "hist_flow_rate_sim.dat");

    std::cout << "\nCheck finished. Outputs saved to '" << output_dir << "'" << std::endl;
    return 0;
}