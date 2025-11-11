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
#include <numeric>

struct Params
{
    std::string filename;
    int nreals = 1;
    int nbins = 50;
    int randomize_flag = 0;
    std::string coords_file = "";
    int pbc_flag = 0;
};

RNG Rand;

void fix_disorder(Graph& G);
Params read_params(int argc, char** argv);
void initialize_flows(Graph& G);
void dfs_sum_products( std::size_t current_node_idx, double current_product, const Graph& G_original, const std::vector<bool>& is_big_tube, std::map<size_t, double>& target_weights, std::vector<bool>& visited);
Graph create_big_tubes_graph(const Graph& G_original);
void create_output_directory(const std::string& dir_path);


// Creates a Gnuplot script to plot simulation vs. analytical data.
void create_all_plots(const std::string& output_dir, const std::string& base_name_tubes, const std::string& base_name_pores, const std::string& analytical_base_name, double mu, double sigma2) 
{
    // Derive the other histogram filenames automatically
    std::string big_in_tubes_hist_file = "big_in_" + base_name_tubes;
    std::string big_out_tubes_hist_file = "big_out_" + base_name_tubes;
    std::string small_in_tubes_hist_file = "small_in_" + base_name_tubes;
    std::string small_out_tubes_hist_file = "small_out_" + base_name_tubes;

    // SCRIPT 1: The main plot with P(q) and Alim fit 
    {
        std::string script_path = output_dir + "tube_flow_dist_plot.gp";
        std::ofstream gpout(script_path);
        
        gpout << "# --- Main Plot: Simulation vs. Full Analytical Fit and Gamma (Alim) Fit ---" << std::endl;
        gpout << "reset" << std::endl;
        
        // Alim Fit parameters (Standard Gamma from moments)
        gpout << "k_alim = (" << mu << "**2) / " << sigma2 << std::endl;
        gpout << "theta_alim = " << sigma2 << " / " << mu << std::endl;
        gpout << "GammaPDF(x) = (1.0/(gamma(k_alim)*theta_alim**k_alim)) * x**(k_alim-1) * exp(-x/theta_alim)" << std::endl;
        
        gpout << "set terminal wxt size 900,600 enhanced" << std::endl;
        gpout << "set title 'Full Distribution vs. Fits'" << std::endl;
        gpout << "set xlabel 'Flow Rate (q)'" << std::endl;
        gpout << "set ylabel 'Probability Density'" << std::endl;
        
        gpout << "plot '" << base_name_tubes << "' u 2:3 with linespoints title 'All Tubes Data', \\" << std::endl;
        gpout << "     '" << analytical_base_name << "_P_fit.dat' with lines lw 2 title 'P(q) Analytical Fit', \\" << std::endl;
        gpout << "     GammaPDF(x) with lines lw 2 dashtype 2 title 'Alim Fit (Gamma PDF)'" << std::endl;
        gpout.close();
        std::cout << "Gnuplot script saved to " << script_path << std::endl;
    }

    // SCRIPT 2: Big tubes vs. f(q)
    {
        std::string script_path = output_dir + "big_tube_flow_dist_plot.gp";
        std::ofstream gpout(script_path);

        gpout << "# --- Big Tubes Plot: Simulation vs. f(q) component ---" << std::endl;
        gpout << "set xlabel 'Flow rate (q)'" << std::endl;
        gpout << "set ylabel 'Probability density'" << std::endl;

        gpout << "plot '" << big_in_tubes_hist_file << "' u 2:3 with linespoints title 'Big-in tubes data', \\" << std::endl;
        gpout << "     '" << big_out_tubes_hist_file << "' u 2:3 with linespoints title 'Big-out tubes data', \\" << std::endl;
        gpout << "     '" << base_name_pores << "' u 2:3 with linespoints title 'Pore data', \\" << std::endl;
        gpout << "     '" << analytical_base_name << "_g_fit.dat' with lines lw 2 title 'Model fit'" << std::endl;
        gpout.close();
        std::cout << "Gnuplot script saved to " << script_path << std::endl;
    }

    // SCRIPT 3: Small tubes vs. g(q) 
    {
        std::string script_path = output_dir + "small_tube_flow_dist_plot.gp";
        std::ofstream gpout(script_path);

        gpout << "# Small Tubes Plot: Simulation vs. g(q) component" << std::endl;
        gpout << "set xlabel 'Flow rate (q)'" << std::endl;
        gpout << "set ylabel 'Probability density'" << std::endl;

        gpout << "plot '" << small_in_tubes_hist_file << "' u 2:3 with linespoints title 'Small-in tubes data', \\" << std::endl;
        gpout << "     '" << small_out_tubes_hist_file << "' u 2:3 with linespoints title 'Small-out tubes data', \\" << std::endl;
        gpout << "     '" << analytical_base_name << "_f_fit.dat' with lines lw 2 title 'Model prediction'" << std::endl;
        gpout.close();
        std::cout << "Gnuplot script saved to " << script_path << std::endl;
    }
}

void save_inlet_outlet_coordinates(const Graph& G, const std::string& output_dir)
{
    // Save Inlet Coordinates
    std::string inlet_filename = output_dir + "inlet_coordinates.dat";
    std::ofstream fout_in(inlet_filename);
    if(!fout_in) std::cerr << "Error: Unable to open file: " << inlet_filename << std::endl;
    else 
    {
        for(std::size_t node_id : G.inlet) 
        {
            fout_in << node_id;
            for (double coord : G.nodes[node_id].coordinates) fout_in << " " << coord;
            fout_in << std::endl;
        }
        fout_in.close();
        std::cout << "Inlet coordinates saved to " << inlet_filename << std::endl;
    }

    // Save Outlet Coordinates
    std::string outlet_filename = output_dir + "outlet_coordinates.dat";
    std::ofstream fout_out(outlet_filename);
    if(!fout_out) std::cerr << "Error: Unable to open file: " << outlet_filename << std::endl;
    else 
    {
        for(std::size_t node_id : G.outlet) 
        {
            fout_out << node_id;
            for (double coord : G.nodes[node_id].coordinates) fout_out << " " << coord;
            fout_out << std::endl;
        }
        fout_out.close();
        std::cout << "Outlet coordinates saved to " << outlet_filename << std::endl;
    }
}

double Pearson_correlation(const std::vector<double>& pair_products, const std::vector<double>& series_1, const std::vector<double>& series_2)
{
    double pair_products_mean, trash, mean_1, mean_2, var_1, var_2, covariance;
    myfun::med_var(series_1, mean_1, var_1);
    myfun::med_var(series_2, mean_2, var_2);
    myfun::med_var(pair_products, pair_products_mean, trash);

    covariance = pair_products_mean - mean_1*mean_2;
    return covariance/std::sqrt(var_1)/std::sqrt(var_2);
}

void save_k_data(const std::string& filename, const std::vector<std::size_t>& k_values, const std::vector<double>& y_values, const std::vector<double>& y_errors)
{
    std::ofstream fout(filename);
    if(!fout) 
    {
        std::cerr << "Error: Unable to open file: " << filename << std::endl;
        return;
    }

    fout << "# k\tvalue\tstd_error_of_mean\n";
    for(std::size_t i = 0; i < k_values.size(); ++i) 
        fout << k_values[i] << "\t" << y_values[i] << "\t" << y_errors[i] << "\n";
    
    fout.close();
    std::cout << "k-dependent data saved to " << filename << std::endl;
}

int main(int argc, char** argv){

    // Read parameters (number of histogram bins & number of realizations of the network:
    Params p = read_params(argc, argv);
    std::cout << "nbins = " << p.nbins << std::endl << "nreals = " << p.nreals << std::endl;

    // Create output directory:
    const std::string output_dir = "output/";
    create_output_directory(output_dir);

    // We read the graph from its list of links:
    Graph G(p.filename, p.coords_file);

    // Check if the graph is a DAG
    bool is_dag = G.is_DAG();
    std::cout << "Graph is a DAG: " << (is_dag ? "true" : "false") << std::endl;

    

    // Save boundary coordinates if coordinates were loaded 
    if (!p.coords_file.empty() && !G.nodes.empty()) save_inlet_outlet_coordinates(G, output_dir);
    

    // We print some interesting information about the pore network
    std::size_t n_bulk = G.nodes.size() - G.inlet.size() - G.outlet.size();
    std::cout << "Number of inlets = " << G.inlet.size() << ", Number of outlets = " << G.outlet.size() << std::endl;
    std::cout << "Number of bulk nodes = " << n_bulk << std::endl;
    std::cout << "Among the bulk nodes there are: " << std::endl;
    

    std::size_t n_type_1, n_type_2;
    n_type_1 = n_type_2 = 0;
    for(auto& junction : G.nodes)
    {
        if((junction.outnbrs.size() == 2)&&(junction.innbrs.size() == 1)) 
        {
            ++ n_type_2;
            junction.state = 2;
        } 
        
        if((junction.outnbrs.size() == 1)&&(junction.innbrs.size() == 2))
        {
            ++ n_type_1;
            junction.state = 1;
        }
         
    }
    std::cout << n_type_1 << " type 1 junctions(merge) " << std::endl;
    std::cout << n_type_2 << " type 2 junctions(split) " << std::endl;
    std::cout << n_bulk - n_type_1 - n_type_2 << " others" << std::endl;

    // We compute the (k=1) degree distribution:
    std::vector<double> in_degrees, out_degrees;
    double min_in_degree, max_in_degree, min_out_degree, max_out_degree, delta;
    for(const auto& node : G.nodes)
    {
        in_degrees.push_back(static_cast<double>(node.innbrs.size()));
        out_degrees.push_back(static_cast<double>(node.outnbrs.size()));
    }
    
    auto degree_histogram = myfun::buildHistogram(in_degrees, p.nbins, min_in_degree, max_in_degree, delta, 1.0);
    myfun::guardaHistograma(degree_histogram, min_in_degree, delta, output_dir + "in_degree_dist.dat");
    std::cout << "Max in-degree = " << max_in_degree << std::endl;
    std::cout << "Min in-degree = " << min_in_degree << std::endl;

    degree_histogram = myfun::buildHistogram(out_degrees, p.nbins, min_out_degree, max_out_degree, delta, 1.0);
    myfun::guardaHistograma(degree_histogram, min_out_degree, delta, output_dir + "out_degree_dist.dat");
    std::cout << "Max out-degree = " << max_out_degree << std::endl;
    std::cout << "Min out-degree = " << min_out_degree << std::endl;


    // We compute some statistics about order-k neighbours:
    std::cout << "\n--- STARTING K-ORDER ANALYSIS LOOP ---" << std::endl;

    // Define k values
    std::vector<std::size_t> k_values;
    for(std::size_t i = 1; i <= 150; ++i) k_values.push_back(i);
    std::size_t max_k_for_analysis = *std::max_element(k_values.begin(), k_values.end());
    
    // Define summary vectors
    std::vector<double> k_avg_in_degree, k_avg_out_degree, k_avg_in_degree_err, k_avg_out_degree_err;
    std::vector<double> k_corr_evenness, k_corr_mass;
    std::vector<double> k_avg_cum_in_degree, k_avg_cum_out_degree;
    std::vector<double> k_avg_cum_in_degree_err, k_avg_cum_out_degree_err;
    

    // Get global 'evenness_values' and 'mass_values' (only need to do this once)
    std::vector<double> mass_values;
    std::vector<double> evenness_values; // global population of splitters
    initialize_flows(G);
    G.compute_stationary_flows(Rand, (p.pbc_flag == 1)); // Compute flows once

    for(const auto& node : G.nodes) 
    {
        mass_values.push_back(node.mass);
        if(node.outnbrs.size() == 2) 
        {
            double evenness = 2.0 * std::min(node.outwgs[0], node.outwgs[1]);
            evenness_values.push_back(evenness);
        }
    }

    // We don't need the global means, but good to have
    double mass_mean, mass_var, evenness_mean, evenness_var;
    myfun::med_var(mass_values, mass_mean, mass_var);
    myfun::med_var(evenness_values, evenness_mean, evenness_var);

    
    std::vector<std::vector<std::size_t>> out_dist_hist;
    std::vector<std::vector<std::size_t>> in_dist_hist;
    G.compute_all_distance_histograms(max_k_for_analysis, out_dist_hist, in_dist_hist);
    


    // Start the loop
    for(std::size_t k : k_values)
    {
        std::cout << "\n--- Analyzing k = " << k << " ---" << std::endl;
        
        // Compute the k-th power matrix
        G.compute_weighted_power(k);

        // Define vectors for this k
        std::vector<double> k_in_degrees(G.nodes.size(), 0.0);
        std::vector<double> k_out_degrees(G.nodes.size(), 0.0);
        std::vector<double> t_fractions;
        std::vector<double> mass_pair_products, evenness_pair_products;
        std::vector<double> mass_values_i, mass_values_j;
        std::vector<double> evenness_values_i, evenness_values_j;

        const double epsilon = 1e-15; // Tolerance for non-zero check

        for(int i = 0; i < G.weights_matrix_powered.outerSize(); ++i) 
        {
            for(Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(G.weights_matrix_powered, i); it; ++it)
            {
                // Only count entries that are physically significant
                if (std::abs(it.value()) > epsilon) 
                {
                    std::size_t j = it.col();
                    k_out_degrees[i] += 1.0;
                    k_in_degrees[j] += 1.0;
                    t_fractions.push_back(it.value());

                    // Collect pairs for correlation
                    mass_values_i.push_back(G.nodes[i].mass);
                    mass_values_j.push_back(G.nodes[j].mass);
                    mass_pair_products.push_back(G.nodes[i].mass * G.nodes[j].mass);
                    
                    if((G.nodes[i].outnbrs.size() == 2) && (G.nodes[j].outnbrs.size() == 2))
                    {
                        double evenness_i = 2.0 * std::min(G.nodes[i].outwgs[0], G.nodes[i].outwgs[1]);
                        double evenness_j = 2.0 * std::min(G.nodes[j].outwgs[0], G.nodes[j].outwgs[1]);
                        
                        evenness_values_i.push_back(evenness_i);
                        evenness_values_j.push_back(evenness_j);
                        evenness_pair_products.push_back(evenness_i * evenness_j);
                    }
                }
            }
        } // end of matrix iteration

        // Save distributions for this k 
        std::string k_suffix = "_k=" + std::to_string(k) + ".dat";
        double min_val, max_val;
        
        auto k_degree_hist = myfun::buildHistogram(k_in_degrees, p.nbins, min_val, max_val, delta);
        myfun::guardaHistograma(k_degree_hist, min_val, delta, output_dir + "order_k_in_degree" + k_suffix);
        
        k_degree_hist = myfun::buildHistogram(k_out_degrees, p.nbins, min_val, max_val, delta);
        myfun::guardaHistograma(k_degree_hist, min_val, delta, output_dir + "order_k_out_degree" + k_suffix);
        
        auto k_t_frac_hist = myfun::buildHistogram(t_fractions, p.nbins, min_val, max_val, delta);
        myfun::guardaHistograma(k_t_frac_hist, min_val, delta, output_dir + "t_fractions_dist" + k_suffix);

        // Calculate and store summary statistics for this k 
        double mean, var, sem;

        // In-degree
        myfun::med_var(k_in_degrees, mean, var);
        sem = std::sqrt(var/k_in_degrees.size());
        k_avg_in_degree.push_back(mean);
        k_avg_in_degree_err.push_back(sem);
        std::cout << "Avg k-in-degree = " << mean << " +/- " << sem << std::endl;

        // Out-degree
        myfun::med_var(k_out_degrees, mean, var);
        sem = std::sqrt(var/k_out_degrees.size());
        k_avg_out_degree.push_back(mean);
        k_avg_out_degree_err.push_back(sem);
        std::cout << "Avg k-out-degree = " << mean << " +/- " << sem << std::endl;

        // Mass correlation
        double corr_mass = Pearson_correlation(mass_pair_products, mass_values_i, mass_values_j);
        k_corr_mass.push_back(corr_mass);
        std::cout << "Mass Correlation = " << corr_mass << std::endl;

        // Evenness Correlation
        double corr_evenness = Pearson_correlation(evenness_pair_products, evenness_values_i, evenness_values_j);
        k_corr_evenness.push_back(corr_evenness);
        std::cout << "Evenness Correlation = " << corr_evenness << std::endl;

        // cumulative distances
        std::cout << "Computing cumulative degrees for k <= " << k << "..." << std::endl;

        std::vector<double> cum_in_deg_vec, cum_out_deg_vec;
        cum_in_deg_vec.reserve(G.nodes.size());
        cum_out_deg_vec.reserve(G.nodes.size());

        for(std::size_t i = 0; i < G.nodes.size(); ++i)
        {
            std::size_t cum_out_count = 0;
            std::size_t cum_in_count = 0;
            
            for(std::size_t d = 1; d <= k; ++d) 
            {
                cum_out_count += out_dist_hist[i][d];
                cum_in_count += in_dist_hist[i][d];
            }
            cum_out_deg_vec.push_back(static_cast<double>(cum_out_count));
            cum_in_deg_vec.push_back(static_cast<double>(cum_in_count));
        }

        // Cumulative In-degree stats
        myfun::med_var(cum_in_deg_vec, mean, var);
        sem = std::sqrt(var/cum_in_deg_vec.size());
        k_avg_cum_in_degree.push_back(mean);
        k_avg_cum_in_degree_err.push_back(sem);

        // Cumulative Out-degree stats
        myfun::med_var(cum_out_deg_vec, mean, var);
        sem = std::sqrt(var/cum_out_deg_vec.size());
        k_avg_cum_out_degree.push_back(mean);
        k_avg_cum_out_degree_err.push_back(sem);

    
    } 


    // Save all "vs. k" data 
    std::cout << "\n--- K-Order Analysis Complete. Saving summary files. ---" << std::endl;
    
    
    // Save summary files using the helper function
    save_k_data(output_dir + "k_avg_in_degree.dat", k_values, k_avg_in_degree, k_avg_in_degree_err);
    save_k_data(output_dir + "k_avg_out_degree.dat", k_values, k_avg_out_degree, k_avg_out_degree_err);
    save_k_data(output_dir + "k_corr_mass.dat", k_values, k_corr_mass, std::vector<double>(k_values.size(), 0.0)); // No error for corr
    save_k_data(output_dir + "k_corr_evenness.dat", k_values, k_corr_evenness, std::vector<double>(k_values.size(), 0.0)); // No error for corr
    save_k_data(output_dir + "k_avg_CUM_in_degree.dat", k_values, k_avg_cum_in_degree, k_avg_cum_in_degree_err);
    save_k_data(output_dir + "k_avg_CUM_out_degree.dat", k_values, k_avg_cum_out_degree, k_avg_cum_out_degree_err);

    
    
    
    // We compute the stationary flows for p.nreals realizations of the disorder (weights):
    std::vector<double> mass_series, mass_1_series, mass_2_series;
    std::vector<double> flow_series, big_in_flow_series, big_out_flow_series, small_in_flow_series, small_out_flow_series; 
    double minflow, maxflow;
    bool use_pbc = (p.pbc_flag == 1);
    
    for(std::size_t e = 0; e < p.nreals; ++e)
    {
        
        initialize_flows(G); 
        if(p.randomize_flag == 1) fix_disorder(G);
        // G.check_mass_conservation();
        G.compute_stationary_flows(Rand, use_pbc);
        
        //We collect flow-rate data at pores (junctions):
        for(auto& pore : G.nodes) if(!pore.is_boundary)
        {
            mass_series.push_back(pore.mass);
            if(pore.state == 1) mass_1_series.push_back(pore.mass);
            if(pore.state == 2) mass_2_series.push_back(pore.mass);
        }
        
        //We collect flow-rate data at tubes:
        for(auto& tube : G.links) 
        {
            if((!G.nodes[tube.in].is_boundary)&&(!G.nodes[tube.out].is_boundary))
            {
                double flow = G.nodes[tube.in].mass*tube.weight;
                flow_series.push_back(flow);
                if((G.nodes[tube.in].innbrs.size() == 2)) big_in_flow_series.push_back(flow);
                else small_in_flow_series.push_back(flow);
                if((G.nodes[tube.out].outnbrs.size() == 2)) big_out_flow_series.push_back(flow);
                else small_out_flow_series.push_back(flow);  
            }
                
        }         
    }

    // Calculate the mean and variance of the raw data
    double flow_mean, flow_variance;
    // Esta primera llamada era solo para imprimir, pero se sobrescribía
    myfun::med_var(flow_series, flow_mean, flow_variance); 
    std::cout << "Calculated from data: mean = " << flow_mean << ", variance = " << flow_variance << std::endl;

    // Call Python script to generate the analytical curve
    std::cout << "Executing Python script to generate analytical fit..." << std::endl; 
    
    std::string  analytical_base_path = output_dir + "analytical";
    
    // --- ¡ESTA ES LA PARTE IMPORTANTE! ---
    // Aquí calculabas la media y varianza de 'mass_series'
    myfun::med_var(flow_series, flow_mean, flow_variance, 99.5);
    
    // Y aquí se las pasabas a Python
    std::string command = "python calculate_fit.py " + std::to_string(flow_mean) + " " + std::to_string(flow_variance) + " " + analytical_base_path; 
    
    // Execute the command
    int result = system(command.c_str());
    if (result != 0) 
        std::cerr << "Warning: Python script may have failed." << std::endl;
   

    

    // Compute correspoinding histograms
    std::string tube_hist_filename = "tube_flow_dist.dat";
    auto histogram = myfun::buildHistogram(mass_series, p.nbins, minflow, maxflow, delta);
    myfun::guardaHistograma(histogram, minflow, delta, output_dir + "pore_flow_dist.dat");
    histogram = myfun::buildHistogram(mass_1_series, p.nbins, minflow, maxflow, delta);
    myfun::guardaHistograma(histogram, minflow, delta, output_dir + "pore_1_flow_dist.dat");
    histogram = myfun::buildHistogram(mass_2_series, p.nbins, minflow, maxflow, delta);
    myfun::guardaHistograma(histogram, minflow, delta, output_dir + "pore_2_flow_dist.dat");
    histogram = myfun::buildHistogram(flow_series, p.nbins, minflow, maxflow, delta);
    myfun::guardaHistograma(histogram, minflow, delta, output_dir + tube_hist_filename);
    histogram = myfun::buildHistogram(big_in_flow_series, p.nbins, minflow, maxflow, delta);
    myfun::guardaHistograma(histogram, minflow, delta, output_dir + "big_in_tube_flow_dist.dat");
    histogram = myfun::buildHistogram(big_out_flow_series, p.nbins, minflow, maxflow, delta);
    myfun::guardaHistograma(histogram, minflow, delta, output_dir + "big_out_tube_flow_dist.dat");
    histogram = myfun::buildHistogram(small_in_flow_series, p.nbins, minflow, maxflow, delta);
    myfun::guardaHistograma(histogram, minflow, delta, output_dir + "small_in_tube_flow_dist.dat");
    histogram = myfun::buildHistogram(small_out_flow_series, p.nbins, minflow, maxflow, delta);
    myfun::guardaHistograma(histogram, minflow, delta, output_dir + "small_out_tube_flow_dist.dat");

    // Checking mass conservation: we display total inlet mass vs. total outlet mass
    double inlet_mass, outlet_mass;
    inlet_mass = outlet_mass = 0.0;
    for(auto& v : G.inlet) inlet_mass += G.nodes[v].mass;
    for(auto& v : G.outlet) outlet_mass += G.nodes[v].mass;
    std::cout << "Total inlet flow = " << inlet_mass << ", Total outlet flow = " << outlet_mass << std::endl;
    std::cout << "Graph has " << G.inlet.size() << " inlets." << std::endl;

    // Create the Gnuplot script for fitting
    myfun::med_var(flow_series, flow_mean, flow_variance, 99.5);
    create_all_plots(output_dir, tube_hist_filename, "pore_flow_dist.dat", "analytical", flow_mean, flow_variance);


    // We compute the fractions distribution
    std::vector<double> fractions;
    double fracmin, fracmax;
    for(auto& link : G.links ) fractions.push_back(link.weight);
    histogram = myfun::buildHistogram(fractions, p.nbins, fracmin, fracmax, delta);
    myfun::guardaHistograma(histogram, fracmin, delta, output_dir + "fractions_dist.dat");



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
        if(G.nodes[i].outnbrs.size() == 1)
        {
            G.nodes[i].outwgs[0] = 1.0;
            G.links[G.nodes[i].outlnks[0]].weight =  G.nodes[i].outwgs[0];
        }
        
        if(G.nodes[i].outnbrs.size()>2) 
        {
            double wgssum = 0.0;
            for(auto& weight : G.nodes[i].outwgs)
            {
                weight = Rand();
                wgssum += weight;
            }
            for(std::size_t j = 0; j < G.nodes[i].outnbrs.size(); ++j)
            {
                G.nodes[i].outwgs[j] /= wgssum;
                G.links[G.nodes[i].outlnks[j]].weight =  G.nodes[i].outwgs[j];   
            }   

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
    std::cerr << "Usage: " << argv[0] << " <filename> [n_realizations] [n_bins] [randomize_flag] [coords_file] [pbc_flag]" << std::endl;
    std::cerr << "  <filename>       : Path to the links file." << std::endl;
    std::cerr << "  [n_realizations] : (Default 1) Number of realizations." << std::endl;
    std::cerr << "  [n_bins]         : (Default 50) Number of histogram bins." << std::endl;
    std::cerr << "  [randomize_flag] : (Default 0) 0 for quenched, 1 for annealed." << std::endl;
    std::cerr << "  [coords_file]    : (Default \"\") Path to coordinates file. Use \"\" to skip." << std::endl;
    std::cerr << "  [pbc_flag]       : (Default 0) 0 for non-PBC, 1 for PBC solver." << std::endl;
    std::exit(1);
   }

   p.filename = argv[1];
   if(argc > 2) p.nreals = std::atoi(argv[2]);
   if(argc > 3) p.nbins = std::atoi(argv[3]);
   if(argc > 4) p.randomize_flag = std::atoi(argv[4]);
   if(argc > 5) p.coords_file = argv[5];
   if(argc > 6) p.pbc_flag = std::atoi(argv[6]); 

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

