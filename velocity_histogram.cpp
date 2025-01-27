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


        double cos2 = velx2/velmag2;
        double sin2 = 1.0 - cos2;
        double sincos = U[cellI].x()*U[cellI].y()/velmag2;
        double alpha = gradU()[cellI](0,1) + gradU()[cellI](1,0);
        double beta = gradU()[cellI](1,1) - gradU()[cellI](0,0);
        double delww = gradU()[cellI](0,0)*cos2 + gradU()[cellI](1,1)*sin2 + alpha*sincos;
        double delwz = gradU()[cellI](0,1)*cos2 - gradU()[cellI](1,0)*sin2 + beta*sincos;
        double delzw = gradU()[cellI](1,0)*cos2 - gradU()[cellI](0,1)*sin2 + beta*sincos;
        double delzz = gradU()[cellI](1,1)*cos2 + gradU()[cellI](0,0)*sin2 - alpha*sincos;


        double volume = mesh.V()[cellI]*1E10;
        double sigma = gradU()[cellI](0,1);
        vels.push_back(velmag);
        volumes.push_back(volume);
        shearXY.push_back(sigma);
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
    
    return 0;
    
}

