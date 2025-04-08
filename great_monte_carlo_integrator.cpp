#define NMCS 100000000
#define NBINS 50

#include "myfunctions.h"
#include <vector>
#include <fstream>
#include <cmath>
#include <string>
#include <random>

// RNG plano entre 0 y 1
double Rand(void)
{
    // Usamos variables estaticas para inicializar el generador y la distribucion solo una vez
    static std::random_device rd;
    static std::mt19937 gen(rd()); // Mersenne Twister, es god
    static std::uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(gen);
}

int main(){

    std::vector<double> weights(NMCS), y_coords(NMCS), vels(NMCS), shears(NMCS), half_widths(NMCS), flux(NMCS);
    std::vector<double> histo_a(NBINS), histo_y(NBINS), histo_vels(NBINS), histo_shear(NBINS);
    double y, v, s, a, q, alpha;
    double max, min;
    int i;
    std::string filename;


    //Inicializamos todos los pesos a 1. Realmente no deberiamos usar pesos, pero mi funcion "histograma sin pesos" devuelve conteos enteros
    // en vez de doubles ya normalizados y me da palo cambiarla

    for(i=0; i<NMCS; i++) weights[i] = 1.0;

    
    alpha = 1.0;
    for(i=0; i<NMCS; i++)
    {
        a = 0.4*Rand()+ 0.8; //the half-width of the Poiseuille tube 
        y = (2.0*Rand()-1.0)*a; //the transverse coordinate, sampled uniformly
        v = alpha*(a*a-y*y); // eulerian velocity
        s = 2.0*alpha*std::fabs(y); //shear 

        //guardamos los valores en vectores:
        half_widths[i] = a;
        y_coords[i] = y;
        vels[i] = v;
        shears[i] = s;
    }

    histo_a = myfun::weightedHistogram(half_widths, weights, NBINS, min, max);
    myfun::guardaWeightedHistograma(histo_a, min, max, NBINS, "histograma_a.txt");
    histo_y = myfun::weightedHistogram(y_coords, weights, NBINS, min, max);
    myfun::guardaWeightedHistograma(histo_y, min, max, NBINS, "histograma_y.txt");
    histo_vels = myfun::weightedHistogram(vels, weights, NBINS, min, max);
    myfun::guardaWeightedHistograma(histo_vels, min, max, NBINS, "histograma_v.txt");
    histo_shear = myfun::weightedHistogram(shears, weights, NBINS, min, max);
    myfun::guardaWeightedHistograma(histo_shear, min, max, NBINS, "histograma_s.txt");


    
     //Creamos el histograma 2D de modulo de la velocidad y shear
     double min1, min2, max1, max2;
     std::vector<double> shear_given_v(NBINS,0.0), sh_error_given_v(NBINS,0.0), true_means(NBINS*NBINS, 0.0), errors(NBINS*NBINS, 0.0);
     min1 = 1.0e-3;
     max1 = *std::max_element(vels.begin(), vels.end());
     min2 = 1.0e-4;
     max2 = *std::max_element(shears.begin(), shears.end());
     auto histoVelShear = myfun::twoDLogHistogram(vels, shears, weights, true_means, errors, NBINS, NBINS, min1, max1, min2, max2);
     if(min1*min2==0.0) std::cout<<"liadon historico";
     std::cout<<"histograma creado\n";
 
     //Calculamos el shear promedio por velocidad y lo guardamos en un vector
    
     double v_i, s_i, logfv, logfs, fv, fs;
 
     v_i = min1;
     logfv = std::log(max1/min1)/static_cast<double>(NBINS); 
     logfs = std::log(max2/min2)/static_cast<double>(NBINS);
     fv = std::exp(logfv);
     fs = std::exp(logfs);
     
     for(i = 0; i<NBINS; i++)
     {
         s_i = min2;
         double row_sum, row_mean, row_mean2;
         row_sum = row_mean = row_mean2 = 0.0;
         for(std::size_t j = 0; j<NBINS; ++j)
         {
             row_sum += histoVelShear[i*NBINS+j];
             row_mean += true_means[i*NBINS+j];
             row_mean2 += errors[i*NBINS+j];
         }
         if(row_sum!=0.0)
         {
             row_mean/= row_sum;
             row_mean2/= row_sum;
             shear_given_v[i] = row_mean;  
             sh_error_given_v[i] = std::sqrt((row_mean2-row_mean*row_mean)/row_sum);
         }
         else shear_given_v[i] = sh_error_given_v[i] = 0.0;
         
         for(std::size_t j = 0; j<NBINS; ++j) 
         {
             if(row_sum!=0.0) histoVelShear[i*NBINS+j]/= (row_sum*v_i*(fv-1.0)*s_i*(fs-1.0));
             s_i *= fs;
         }
         v_i*=fv;
         
     }
 
 
     //Guardamos
     filename = "shear_given_v.dat";
     std::ofstream outFile(filename);
     if(!outFile.is_open()) 
     {
         std::cerr << "Error: No se pudo abrir el archivo " << filename << " para escribir.\n";
         return 0;
     }
     double bincenter, x_i = min1;
     double logf = std::log(max1/min1)/static_cast<double>(NBINS);
     double sqf = std::exp(0.5*logf);
     for(i = 0; i < NBINS; ++i)
     {
         bincenter = sqf*x_i;
         outFile << bincenter << "\t" << shear_given_v[i] << "\t" << sh_error_given_v[i] << "\n";
         x_i = bincenter*sqf;
     }
        
     outFile.close();
     std::cout << "Shear given vel ponderado y logaritmico guardado en " << filename << "\n";
     



    return 0;

}