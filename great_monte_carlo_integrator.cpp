#define NMCS 100000
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


    //Inicializamos todos los pesos a 1. Realmente no deberiamos usar pesos, pero mi funcion "histograma sin pesos" devuelve conteos enteros
    // en vez de doubles ya normalizados y me da palo cambiarla

    for(i=0; i<NMCS; i++) weights[i] = 1.0;

    
    alpha = 1.0;
    for(i=0; i<NMCS; i++)
    {
        a = Rand()+ 0.5; //the half-width of the Poiseuille tube 
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


    




    return 0;

}