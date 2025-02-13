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
        shearXY.push_back(std::fabs(delxy));
        jacww.push_back(std::fabs(delww));
        jacwz.push_back(std::fabs(delwz));
        jaczw.push_back(std::fabs(delzw));
        jaczz.push_back(std::fabs(delzz));

        //gradU[cellI](0,1) es la parcial de vx respecto de y (etc)
    }

    
    

    //Creamos un histograma con N bins
    double minVal, maxVal;
    std::size_t N = 40;
    minVal = 1.0e-7;
    maxVal = *std::max_element(vels.begin(), vels.end());
    std::cout << "max y min son " << minVal <<" "<< maxVal << std::endl;
    auto histoU = myfun::weightedLogHistogram(vels, volumes, N, minVal, maxVal);

    
    //Guardamos resultados en un fichero
    std::string filename = "histograma_U.dat";
    myfun::guardaWeightedLogHistograma(histoU, minVal, maxVal, N, filename);


    //Hacemos lo mismo con los shears:
    
    minVal = 1.0e-4;
    maxVal = *std::max_element(shearXY.begin(), shearXY.end());
    auto histoShearXY = myfun::weightedLogHistogram(shearXY, volumes, N, minVal, maxVal);
    filename = "histograma_shear_XY.dat";
    myfun::guardaWeightedLogHistograma(histoShearXY, minVal, maxVal, N, filename);

    minVal = 1.0e-4;
    maxVal = *std::max_element(jacww.begin(), jacww.end());
    auto histoShearWW = myfun::weightedLogHistogram(jacww, volumes, N, minVal, maxVal);
    filename = "histograma_shear_WW.dat";
    myfun::guardaWeightedLogHistograma(histoShearWW, minVal, maxVal, N, filename);

    minVal = 1.0e-4;
    maxVal = *std::max_element(jacwz.begin(), jacwz.end());
    auto histoShearWZ = myfun::weightedLogHistogram(jacwz, volumes, N, minVal, maxVal);
    filename = "histograma_shear_WZ.dat";
    myfun::guardaWeightedLogHistograma(histoShearWZ, minVal, maxVal, N, filename);

    minVal = 1.0e-4;
    maxVal = *std::max_element(jaczw.begin(), jaczw.end());
    auto histoShearZW = myfun::weightedLogHistogram(jaczw, volumes, N, minVal, maxVal);
    filename = "histograma_shear_ZW.dat";
    myfun::guardaWeightedLogHistograma(histoShearZW, minVal, maxVal, N, filename);

    minVal = 1.0e-4;
    maxVal = *std::max_element(jaczz.begin(), jaczz.end());
    auto histoShearZZ = myfun::weightedLogHistogram(jaczz, volumes, N, minVal, maxVal);
    filename = "histograma_shear_ZZ.dat";
    myfun::guardaWeightedLogHistograma(histoShearZZ, minVal, maxVal, N, filename);

    //Creamos el histograma 2D de modulo de la velocidad y shear WZ
    double min1, min2, max1, max2;
    std::vector<double> shear_given_v(N,0.0), sh_error_given_v(N,0.0), true_means(N*N, 0.0), errors(N*N, 0.0);
    min1 = 1.0e-7;
    max1 = *std::max_element(vels.begin(), vels.end());
    min2 = 1.0e-4;
    max2 = *std::max_element(jacwz.begin(), jacwz.end());
    auto histoVelShearWZ = myfun::twoDLogHistogram(vels, jacwz, volumes, true_means, errors, N, N, min1, max1, min2, max2);
    if(min1*min2==0.0) std::cout<<"liadon historico";
    std::cout<<"histograma creado\n";

    //Calculamos el shear promedio por velocidad y lo guardamos en un vector
   
    double v_i, s_i, logfv, logfs, fv, fs;

    v_i = min1;
    logfv = std::log(max1/min1)/static_cast<double>(N); 
    logfs = std::log(max2/min2)/static_cast<double>(N);
    fv = std::exp(logfv);
    fs = std::exp(logfs);
    
    for(std::size_t i = 0; i<N; ++i)
    {
        s_i = min2;
        double row_sum, row_mean, row_mean2;
        row_sum = row_mean = row_mean2 = 0.0;
        for(std::size_t j = 0; j<N; ++j)
        {
            row_sum += histoVelShearWZ[i*N+j];
            row_mean += true_means[i*N+j];
            row_mean2 += errors[i*N+j];
        }
        if(row_sum!=0.0)
        {
            row_mean/= row_sum;
            row_mean2/= row_sum;
            shear_given_v[i] = row_mean;  
            sh_error_given_v[i] = std::sqrt((row_mean2-row_mean*row_mean)/row_sum);
        }
        else shear_given_v[i] = sh_error_given_v[i] = 0.0;
        
        for(std::size_t j = 0; j<N; ++j) 
        {
            if(row_sum!=0.0) histoVelShearWZ[i*N+j]/= (row_sum*v_i*(fv-1.0)*s_i*(fs-1.0));
            s_i *= fs;
        }
        v_i*=fv;
        
    }


    //Guardamos
    filename = "conditioned_2D_histo.dat";
    myfun::guarda2DLogHistograma(histoVelShearWZ, min1, max1, min2, max2, N, N, filename);
    filename = "shear_given_v.dat";
    std::ofstream outFile(filename);
    if(!outFile.is_open()) 
    {
        std::cerr << "Error: No se pudo abrir el archivo " << filename << " para escribir.\n";
        return 0;
    }
    double bincenter, x_i = min1;
    double logf = std::log(max1/min1)/static_cast<double>(N);
    double sqf = std::exp(0.5*logf);
    for(std::size_t i = 0; i < N; ++i)
    {
        bincenter = sqf*x_i;
        outFile << bincenter << "\t" << shear_given_v[i] << "\t" << sh_error_given_v[i] << "\n";
        x_i = bincenter*sqf;
    }
       
    outFile.close();
    std::cout << "Shear given vel ponderado y logaritmico guardado en " << filename << "\n";
    
    
    
    return 0;
    
}

