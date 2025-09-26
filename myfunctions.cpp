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
        maxValOut = 4.5;

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
            if(binIndex < numBins) histogram[binIndex]++;
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

    std::vector<double> weightedLogHistogram(
        const std::vector<double>& data,
        const std::vector<double>& weights,
        std::size_t numBins,
        double minValOut, 
        double maxValOut
    )
    {
        if(data.empty()) throw std::runtime_error("weightedHistogram: El vector de datos está vacío!");
        if(weights.empty()) throw std::runtime_error("weightedHistogram: El vector de pesos está vacío!");
        if(data.size() != weights.size()) throw std::runtime_error("weightedHistogram: Los vectores de datos y pesos deben tener el mismo tamaño!");

        if(numBins == 0) throw std::runtime_error("weightedHistogram: numBins no puede ser 0!");
        

        std::vector<double> histogram(numBins, 0.0);

        // Evitar división por cero si todos los valores son iguales
        if(minValOut == maxValOut)
        {
            histogram[0] = data.size();
            return histogram;
        }

        if(minValOut<=0.0)
        {
            histogram[0] = data.size();
            return histogram;
        }

        double logf, f;

        logf = std::log(maxValOut/minValOut)/static_cast<double>(numBins);
        f = std::exp(logf);

        // Llenar el histograma logaritmico ponderado
        double suma_weights = 0.0;
        for(std::size_t i = 0; i < data.size(); ++i)
        {
            double ratio = std::log(data[i]/minValOut)/logf;
            if(ratio>=0.0)
            {
                std::size_t binIndex = static_cast<std::size_t>(ratio);
                if(binIndex >= numBins) binIndex = numBins-1;
                histogram[binIndex] += weights[i];
                suma_weights += weights[i];
            }  
        }

        // Normalización
        double delta_i, x_i;
        x_i = minValOut;
        for(std::size_t i=0; i<numBins; i++) 
        {
            delta_i = x_i*(f-1.0);
            histogram[i]/=(suma_weights*delta_i);
            x_i *= f;
        }
        return histogram;
    }

    void guardaWeightedLogHistograma(
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

        outFile << "# binCenter\tcount\n";

        double x_i = minVal;
        double logf, sqf;
        
        logf = std::log(maxVal/minVal)/static_cast<double>(numBins);
        sqf = std::exp(0.5*logf);
    


        for(std::size_t i = 0; i < histogram.size(); ++i)
        {
            double binCenter = x_i*sqf;
            outFile << binCenter << "\t" << histogram[i] << "\n";
            x_i = binCenter*sqf;
        }

        outFile.close();
        std::cout << "Histograma logaritmico ponderado guardado en " << filename << "\n";
    }

    std::vector<double> twoDLogHistogram(
        const std::vector<double>& data1,
        const std::vector<double>& data2,
        const std::vector<double>& weights,
        std::vector<double>& true_means,
        std::vector<double>& errors,
        std::size_t numBins1,
        std::size_t numBins2,
        double min1, 
        double max1,
        double min2,
        double max2 
    )
    {
        if(data1.empty()) throw std::runtime_error("twoDHistogram: El vector data1 está vacío!");
        if(data2.empty()) throw std::runtime_error("twoDHistogram: El vector data2 está vacío!");
        if(weights.empty()) throw std::runtime_error("twoDHistogram: El vector de pesos está vacío!");
        if(data1.size() != weights.size()) throw std::runtime_error("twoDHistogram: Los vectores data1 y weights no son compatibles!");
        if (data2.size() != data1.size())
            throw std::runtime_error("twoDHistogram: data2 and data1 must have the same size!");
        if (data2.size() != weights.size())
            throw std::runtime_error("twoDHistogram: data2 and weights must have the same size!");
        if (numBins1 == 0 || numBins2 == 0) throw std::runtime_error("twoDHistogram: numBins no puede ser 0!");
        std::cout << "numbins1 = " << numBins1 << ", " << "numbins2 = " << numBins2 << "\n";
        

        std::cout << "Las velocidades maxima y minima son " << max1 << " y " << min1 << "\n"; 
        std::cout << "Los shears maximo y minimo son " << max2 << " y " << min2 << "\n";
        std::vector<double> histogram(numBins1 * numBins2, 0.0);

        // Evitar división por cero si todos los valores son iguales
        if((min1 == max1) || (min2 == max2))
        {
            histogram[0] = data1.size();
            return histogram;
        }

        double logf1, logf2;
        double eps = 1.0E-30;
        logf1 = std::log(max1/min1)/static_cast<double>(numBins1);
        logf2 = std::log(max2/min2)/static_cast<double>(numBins2);
        std::cout<<"logf1 = "<<logf1<<"logf2 ="<<logf2<<"\n";
       
        // Llenar el histograma 2D
        for(std::size_t i = 0; i < data1.size(); ++i) 
        {
            double ratio1 = std::log(std::fabs(data1[i]/(min1+eps)))/logf1;
            if(ratio1>=0.0) 
            {
                std::size_t binIndex1 = static_cast<std::size_t>(ratio1);
                if(binIndex1 >= numBins1) binIndex1  = numBins1-1;

                double ratio2 = std::log(std::fabs(data2[i]/(min2+eps)))/logf2;
                if(ratio2>=0.0) 
                {
                    std::size_t binIndex2 = static_cast<std::size_t>(ratio2);
                    if(binIndex2 >= numBins2) binIndex2 = numBins2-1;

                    histogram[binIndex1 * numBins2 + binIndex2] += weights[i];
                    true_means[binIndex1 * numBins2 + binIndex2] += weights[i]*data2[i];
                    errors[binIndex1 * numBins2 + binIndex2] += weights[i]*data2[i]*data2[i];
                }
                 
            }
            
        }

        
        return histogram;
    }

    void guarda2DLogHistograma(
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

        double v_i, s_i;
        double logf1, sqf1, logf2, sqf2;
        logf1 = std::log(max1/min1)/static_cast<double>(numBins1);
        sqf1 = std::exp(0.5*logf1);
        logf2 = std::log(max2/min2)/static_cast<double>(numBins2);
        sqf2 = std::exp(0.5*logf2);

        v_i = min1;
        for(std::size_t i = 0; i < numBins1; ++i)
        {
            s_i = min2;
            double binCenter1 = v_i*sqf1;
            for(std::size_t j = 0; j<numBins2; ++j)
            {
                double binCenter2 = s_i*sqf2;
                outFile << i << "\t" << binCenter1 << "\t" << j << "\t" << binCenter2 << "\t" << histogram[i*numBins2 + j] << "\n";
                s_i = binCenter2*sqf2;
            }
            v_i = binCenter1*sqf1;  
        }

        outFile.close();
        std::cout << "Histograma 2D logaritmico ponderado guardado en " << filename << "\n";
    }





} // namespace myfun
