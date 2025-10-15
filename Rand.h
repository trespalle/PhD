#ifndef RAND_H
#define RAND_H

#include <random>

class RNG
{
    std::mt19937 gen; // Mersenne Twister, god RNG
    std::uniform_real_distribution<double> uniform01; // U(0,1) RNG based on Mersenne Twister

    double x, y;
    char flag = 0;

public:
    //constructor:
    RNG(unsigned int seed = 123456789): gen(seed), uniform01(0.0, 1.0) {}

    //public methods:
    double operator()() //U(0,1) RNG
    {
        return uniform01(gen);
    }

    double exp(double theta = 1.0)
    {
        return -theta*std::log((*this)());
    }

    double norm(double mu=0, double sigma=1.0)
    {
        double u, v, t, prefactor;
        do
        {
            u = 2.0*(*this)() -1.0;
            v = 2.0*(*this)() -1.0;
            t = u*u + v*v;
        } while (t >= 1.0 || t == 0.0);
        

        if(flag==0)
        {
            prefactor = std::sqrt(-2.0*std::log(t)/t);
            x = prefactor*u;
            y = prefactor*v;
            ++flag;
            return  sigma*x + mu;   
        }
        else 
        {
            flag = 0;
            return sigma*y + mu;
        }

    }

};

#endif




