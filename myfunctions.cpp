// myfunctions.cpp
#include "myfunctions.h"

namespace myfun {

    std::vector<std::size_t> buildHistogram(
        const std::vector<double>& data,
        std::size_t numBins,
        double& minValOut, 
        double& maxValOut
    )
    {
        if(data.empty()) throw std::runtime_error("buildHistogram: El vector de datos está vacío!");
        if(numBins == 0) throw std::runtime_error("buildHistogram: numBins no puede ser 0!");

        // Encontrar min y max
        minValOut = *std::min_element(data.begin(), data.end());
        maxValOut = *std::max_element(data.begin(), data.end());

        std::vector<std::size_t> histogram(numBins, 0);

        // Evitar división por cero si todos los valores son iguales
        if(minValOut == maxValOut)
        {
            histogram[0] = data.size();
            return histogram;
        }

        // Llenar el histograma
        for(auto val : data)
        {
            double ratio = (val - minValOut)/(maxValOut - minValOut);
            std::size_t binIndex = static_cast<std::size_t>(ratio * numBins);
            if(binIndex == numBins) binIndex--;
            histogram[binIndex]++;
        }

        return histogram;
    }

    std::vector<double> weightedHistogram(
        const std::vector<double>& data,
        const std::vector<double>& weights,
        std::size_t numBins,
        double& minValOut, 
        double& maxValOut
    )
    {
        if(data.empty()) throw std::runtime_error("weightedHistogram: El vector de datos está vacío!");
        if(weights.empty()) throw std::runtime_error("weightedHistogram: El vector de pesos está vacío!");
        if(data.size() != weights.size()) throw std::runtime_error("weightedHistogram: Los vectores de datos y pesos deben tener el mismo tamaño!");

        if(numBins == 0) throw std::runtime_error("weightedHistogram: numBins no puede ser 0!");

        // Encontrar min y max
        minValOut = *std::min_element(data.begin(), data.end());
        maxValOut = *std::max_element(data.begin(), data.end());

        std::vector<double> histogram(numBins, 0.0);

        // Evitar división por cero si todos los valores son iguales
        if(minValOut == maxValOut)
        {
            histogram[0] = data.size();
            return histogram;
        }

        // Llenar el histograma ponderado
        for(std::size_t i = 0; i < data.size(); ++i)
        {
            double ratio = (data[i] - minValOut)/(maxValOut - minValOut);
            std::size_t binIndex = static_cast<std::size_t>(ratio * numBins);
            if(binIndex == numBins) binIndex--;
            histogram[binIndex] += weights[i];
        }

        // Normalización
        double delta = (maxValOut - minValOut)/numBins;
        double suma_weights = 0.0;
        for(auto w : weights) suma_weights += w;
        for(auto &bin : histogram) bin /= (suma_weights * delta);

        return histogram;
    }

    std::vector<double> twoDHistogram(
        const std::vector<double>& data1,
        const std::vector<double>& data2,
        const std::vector<double>& weights,
        std::size_t numBins1,
        std::size_t numBins2,
        double& min1, 
        double& max1,
        double& min2,
        double& max2 
    )
    {
        if(data1.empty()) throw std::runtime_error("twoDHistogram: El vector data1 está vacío!");
        if(data2.empty()) throw std::runtime_error("twoDHistogram: El vector data2 está vacío!");
        if(weights.empty()) throw std::runtime_error("twoDHistogram: El vector de pesos está vacío!");
        if(data1.size() != weights.size()) throw std::runtime_error("twoDHistogram: Los vectores data1 y weights no son compatibles!");

        if (numBins1 == 0 || numBins2 == 0) throw std::runtime_error("twoDHistogram: numBins no puede ser 0!");

        // Encontrar min y max para ambos conjuntos de datos
        min1 = *std::min_element(data1.begin(), data1.end());
        max1 = *std::max_element(data1.begin(), data1.end());
        min2 = *std::min_element(data2.begin(), data2.end());
        max2 = *std::max_element(data2.begin(), data2.end());

        std::vector<double> histogram(numBins1 * numBins2, 0.0);

        // Evitar división por cero si todos los valores son iguales
        if((min1 == max1) || (min2 == max2))
        {
            histogram[0] = data1.size();
            return histogram;
        }

        // Llenar el histograma 2D
        for(std::size_t i = 0; i < data1.size(); ++i) 
        {
            double ratio1 = (data1[i] - min1)/(max1 - min1);
            std::size_t binIndex1 = static_cast<std::size_t>(ratio1 * numBins1);
            if(binIndex1 == numBins1) binIndex1--;

            double ratio2 = (data2[i] - min2)/(max2 - min2);
            std::size_t binIndex2 = static_cast<std::size_t>(ratio2 * numBins2);
            if(binIndex2 == numBins2) binIndex2--;

            histogram[binIndex1 * numBins2 + binIndex2] += weights[i];
        }

        // Normalización
        double delta1 = (max1 - min1)/numBins1;
        double delta2 = (max2 - min2)/numBins2;
        double suma_weights = 0.0;
        for(auto w : weights) suma_weights += w;
        for(auto &bin : histogram) bin /= (suma_weights * delta1 * delta2);

        return histogram;
    }

    void guardaHistograma(
        const std::vector<std::size_t>& histogram,
        double minVal,
        double maxVal,
        std::size_t numBins,
        const std::string& filename
    )
    {
        std::ofstream outFile(filename);
        if(!outFile.is_open()) 
        {
            std::cerr << "Error: No se pudo abrir el archivo " << filename << " para escribir.\n";
            return;
        }

        outFile << "# binIndex\tbinCenter\tcount\n";

        for(std::size_t i = 0; i < histogram.size(); ++i)
        {
            double binCenter = minVal + (maxVal - minVal) * ((static_cast<double>(i) + 0.5)/numBins);
            outFile << i << "\t" << binCenter << "\t" << histogram[i] << "\n";
        }

        outFile.close();
        std::cout << "Histograma guardado en " << filename << "\n";
    }

    void guardaWeightedHistograma(
        const std::vector<double>& histogram,
        double minVal,
        double maxVal,
        std::size_t numBins,
        const std::string& filename
    )
    {
        std::ofstream outFile(filename);
        if(!outFile.is_open()) 
        {
            std::cerr << "Error: No se pudo abrir el archivo " << filename << " para escribir.\n";
            return;
        }

        outFile << "# binIndex\tbinCenter\tcount\n";

        for(std::size_t i = 0; i < histogram.size(); ++i)
        {
            double binCenter = minVal + (maxVal - minVal) * ((static_cast<double>(i) + 0.5)/numBins);
            outFile << i << "\t" << binCenter << "\t" << histogram[i] << "\n";
        }

        outFile.close();
        std::cout << "Histograma ponderado guardado en " << filename << "\n";
    }

    void med_var(
        const std::vector<double>& data,
        double& media,
        double& varianza
    )
    {
        double media_acumulada = 0.0;
        double M2 = 0.0;
        std::size_t N = data.size();

        if (N == 0) throw std::invalid_argument("med_var: El vector de datos está vacío.");
        if (N == 1) throw std::invalid_argument("med_var: Se requiere al menos dos elementos para calcular la varianza.");

        for(std::size_t i = 0; i < N; ++i)
        {
            // Algoritmo de Welford para estabilidad numérica
            double delta = data[i] - media_acumulada;
            media_acumulada += delta / static_cast<double>(i + 1);
            M2 += delta * (data[i] - media_acumulada);
        }

        media = media_acumulada;
        varianza = M2 / static_cast<double>(N - 1);
    }

    void guarda2DHistograma(
        const std::vector<double>& histogram,
        double min1,
        double max1,
        double min2,
        double max2,
        std::size_t numBins1,
        std::size_t numBins2,
        const std::string& filename
    )
    {
        //Chequiar el tamano
        if(histogram.size() != numBins1*numBins2)
        {
            std::cerr << "Error: histogram size does not match de number of bins";
            return;
        }
        
        std::ofstream outFile(filename);
        if(!outFile.is_open()) 
        {
            std::cerr << "Error: No se pudo abrir el archivo " << filename << " para escribir.\n";
            return;
        }

        outFile << "# binIndex1\tbinCenter1\tbinIndex2\tbinCenter2\tcount\n";

        for(std::size_t i = 0; i < numBins1; ++i)
        {
            double binCenter1 = min1 + (max1 - min1) * ((static_cast<double>(i) + 0.5)/numBins1);
            for(std::size_t j = 0; j<numBins2; ++j)
            {
                double binCenter2 = min2 + (max2 - min2) * ((static_cast<double>(j) + 0.5)/numBins2);
                outFile << i << "\t" << binCenter1 << "\t" << j << "\t" << binCenter2 << "\t" << histogram[i*numBins2 + j] << "\n";
            }
            
        }

        outFile.close();
        std::cout << "Histograma 2D ponderado guardado en " << filename << "\n";
    }



} // namespace myfun
