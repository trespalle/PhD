#include "fvCFD.H"
#include "Time.H"
#include "argList.H"
#include "fvSchemes.H"
#include "fvSolution.H"
#include "fvcGrad.H"      // Incluye la definición de fvc::grad

#include "myfunctions.h"

#include <vector>
#include <fstream>
#include <cmath>
#include <string>

int main(int argc, char *argv[])
{
    //Procesamos argumentos con la forma estandar de OpenFOAM
    Foam::argList::validArgs.append("case");
    Foam::argList args(argc, argv);

    //Creamos el objeto runTime para manejar la malla y los tiempos
    Foam::Time runTime
    (
        Foam::Time::controlDictName,  // normalmente "controlDict"
        args,                         // objeto argList
        /* setRootCase = */ false,
        /* setInstance = */ false,
        /* readOption = */ Foam::IOobject::MUST_READ
    );


    //Creamos la malla desde "constant/"
    Foam::fvMesh mesh 
    (
        Foam::IOobject
        (
            Foam::fvMesh::defaultRegion,
            runTime.constant(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    );


    //Forzamos el tiempo a 1047
    runTime.setTime(1047, 0);

    //Leemos la velocidad U en la carpeta "1047"
    Foam::volVectorField U
    (
        Foam::IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            Foam::IOobject::MUST_READ,
            Foam::IOobject::NO_WRITE
        ),
        mesh
    );


    //Calculamos el gradiente de U
    Foam::tmp<Foam::GeometricField<Foam::tensor, Foam::fvPatchField, Foam::volMesh> > gradU = Foam::fvc::grad(U);





    //Preparamos vector para guardar la velocidad X de cada celda
    std::vector<double> vels, shearXY, volumes, jacww, jacwz, jaczw, jaczz;
    vels.reserve(mesh.nCells());
    shearXY.reserve(mesh.nCells());
    volumes.reserve(mesh.nCells());
    jacww.reserve(mesh.nCells());
    jacwz.reserve(mesh.nCells());
    jaczw.reserve(mesh.nCells());
    jaczz.reserve(mesh.nCells());
    

    for(Foam::label cellI = 0; cellI<mesh.nCells(); ++cellI)
    {
        double velx2 = U[cellI].x()*U[cellI].x();
        double vely2 = U[cellI].y()*U[cellI].y();
        double velmag2 = velx2 + vely2; //modulo al cuadrado del campo de velocidades en el centroide de la celda
        double velmag = std::sqrt(velmag2);

        const Foam::tensor& gradUi = gradU()[cellI];
        double delxx = gradUi.xx();
        double delxy = gradUi.xy();
        double delyx = gradUi.yx();
        double delyy = gradUi.yy();



        double cos2 = velx2/velmag2;
        double sin2 = 1.0 - cos2;
        double sincos = U[cellI].x()*U[cellI].y()/velmag2;
        double alpha = delxy + delyx;
        double beta = delyy - delxx;
        double delww = delxx*cos2 + delyy*sin2 + alpha*sincos;
        double delwz = delxy*cos2 - delyx*sin2 + beta*sincos;
        double delzw = delyx*cos2 - delxy*sin2 + beta*sincos;
        double delzz = delyy*cos2 + delxx*sin2 - alpha*sincos;


        double volume = mesh.V()[cellI]*1E10;

        vels.push_back(velmag);
        volumes.push_back(volume);
        shearXY.push_back(delxx);
        jacww.push_back(delww);
        jacwz.push_back(delwz);
        jaczw.push_back(delzw);
        jaczz.push_back(delzz);

        //gradU[cellI](0,1) es la parcial de vx respecto de y (etc)
    }

    
    

    //Creamos un histograma con N bins
    double minVal, maxVal;
    std::size_t N = 50;
    auto histoU = myfun::weightedHistogram(vels, volumes, N, minVal, maxVal);

    
    //Guardamos resultados en un fichero
    std::string filename = "histograma_U.dat";
    myfun::guardaWeightedHistograma(histoU, minVal, maxVal, N, filename);


    //Hacemos lo mismo con los volumenes y shears:
    auto histoVol = myfun::weightedHistogram(volumes, volumes, N, minVal, maxVal);
    filename = "histograma_Vol.dat";
    myfun::guardaWeightedHistograma(histoVol, minVal, maxVal, N, filename);


    auto histoShearXY = myfun::weightedHistogram(shearXY, volumes, N, minVal, maxVal);
    filename = "histograma_shear_XY.dat";
    myfun::guardaWeightedHistograma(histoShearXY, minVal, maxVal, N, filename);

    auto histoShearWW = myfun::weightedHistogram(jacww, volumes, N, minVal, maxVal);
    filename = "histograma_shear_WW.dat";
    myfun::guardaWeightedHistograma(histoShearWW, minVal, maxVal, N, filename);

    auto histoShearWZ = myfun::weightedHistogram(jacwz, volumes, N, minVal, maxVal);
    filename = "histograma_shear_WZ.dat";
    myfun::guardaWeightedHistograma(histoShearWZ, minVal, maxVal, N, filename);

    auto histoShearZW = myfun::weightedHistogram(jaczw, volumes, N, minVal, maxVal);
    filename = "histograma_shear_ZW.dat";
    myfun::guardaWeightedHistograma(histoShearZW, minVal, maxVal, N, filename);

    auto histoShearZZ = myfun::weightedHistogram(jaczz, volumes, N, minVal, maxVal);
    filename = "histograma_shear_ZZ.dat";
    myfun::guardaWeightedHistograma(histoShearZZ, minVal, maxVal, N, filename);

    //Creamos el histograma 2D de modulo de la velocidad y shear WZ
    double min1, min2, max1, max2;
    auto histoVelShearWZ = myfun::twoDHistogram(vels, jacwz, volumes, N, N, min1, max1, min2, max2);
    filename = "histograma2D_vel_shear_WZ.dat";
    myfun::guarda2DHistograma(histoVelShearWZ, min1, max1, min2, max2, N, N, filename);

    //Calculamos el shear promedio por velocidad y lo guardamos en un vector
    //Ademas, renormalizamos el histograma 2D para que cada una de sus filas sea una dist. de prob. condicional (normalizada, por tanto)
   
    std::vector<double> shear_given_v(N,0.0), sh_error_given_v(N,0.0), cond_to_v_histo(N*N, 0.0);
    double sigma;

    for(std::size_t i = 0; i<N; ++i)
    {
        double row_sum, row_mean, row_mean2;
        row_sum = row_mean = row_mean2 = 0.0;
        for(std::size_t j = 0; j<N; ++j)
        {
            sigma = fabs(min2 + (max2-min2)*(static_cast<double>(j)+0.5)/N);
            row_sum += histoVelShearWZ[i*N+j];
            row_mean += histoVelShearWZ[i*N+j]*sigma;
            row_mean2 += histoVelShearWZ[i*N+j]*sigma*sigma;
            cond_to_v_histo[i*N+j] = histoVelShearWZ[i*N+j];
        }
        row_mean/= row_sum;
        row_mean2/= row_sum;
        shear_given_v[i] = row_mean;  
        sh_error_given_v[i] = std::sqrt((row_mean2-row_mean*row_mean)/N);
        for(std::size_t j = 0; j<N; ++j) cond_to_v_histo[i*N+j]/= row_sum;
    }


    //Guardamos

    filename = "shear_given_v.dat";
    std::ofstream outFile(filename);
    if(!outFile.is_open()) 
    {
        std::cerr << "Error: No se pudo abrir el archivo " << filename << " para escribir.\n";
        return 0;
    }
    for(std::size_t i = 0; i < N; ++i)
        outFile << min1 + (max1-min1)*(static_cast<double>(i)+0.5)/N << "\t" << shear_given_v[i] << "\t" << sh_error_given_v[i] << "\n";

    outFile.close();
    std::cout << "Histograma ponderado guardado en " << filename << "\n";

    filename = "conditioned_2D_histo.dat";
    myfun::guarda2DHistograma(cond_to_v_histo, min1, max1, min2, max2, N, N, filename);

    
    return 0;
    
}

