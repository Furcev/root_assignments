#include <fstream>
#include <TMath.h>
#include <TCanvas.h>
#include <iostream>
#include <TH1D.h>
#include <TMinuit.h>
#include <TF1.h>
#include <TFile.h>

double chiSquareValue;   
double signalCounts[100];  
double noiseCounts[100];
double noiseErrors[100];   // Ошибки для шума
double binCenters[100];    // Центры бинов  
double signalErrors[100];  

double Gaussian(double x, double *params) {
    return params[0] * TMath::Gaus(x, params[1], params[2], true);
}


double CombinedFunction(double *xValues, double *params) {
    return (params[0] / (sqrt(2 * M_PI) * params[2])) * TMath::Gaus(xValues[0], params[1], params[2]) + params[3];
}

// Функция для вычисления хи-квадрат
void CalculateChiSquare(int &nParams, double *gradients, double &chiSquare, double *params, int iflag) {
    const int numBins = 100;
    double chiSquared = 0.0;

    for (int i = 0; i < numBins; i++) {
        double deltaSignal = (signalCounts[i] - Gaussian(binCenters[i], params) - params[3]) / signalErrors[i];
        double deltaNoise = (noiseCounts[i] - params[3]) / noiseErrors[i];
        chiSquared += deltaSignal * deltaSignal + deltaNoise * deltaNoise;
    }

    chiSquareValue = chiSquared / numBins;
    chiSquare = chiSquared;
}

// Чтение данных из файла и заполнение гистограммы
void ReadHistogramFromFile(const std::string &filePath, TH1D *histogram) {
    std::ifstream dataFile(filePath);
    double value;

    while (dataFile >> value) {
        histogram->Fill(value);
    }
}


void ex10(){

    auto *signalHistogram = new TH1D("Signal Histogram", "Signal", 100, 500, 600);
    auto *noiseHistogram = new TH1D("Noise Histogram", "Noise", 100, 500, 600);
    

    ReadHistogramFromFile("./data_1.dat", signalHistogram);
    ReadHistogramFromFile("./data_2.dat", noiseHistogram);
    

    for (int i = 0; i < 100; i++) {
        signalCounts[i] = signalHistogram->GetBinContent(i + 1);
        noiseCounts[i] = noiseHistogram->GetBinContent(i + 1);
        binCenters[i] = i + 500; // Центры бинов
        signalErrors[i] = (signalCounts[i] == 0) ? sqrt(3.09) : sqrt(signalCounts[i]);
        noiseErrors[i] = (noiseCounts[i] == 0) ? sqrt(3.09) : sqrt(noiseCounts[i]);
    }
    

    TMinuit *minuit = new TMinuit(4);
    minuit->SetFCN(CalculateChiSquare);

    // Определение ошибок
    double arglist[10];
    arglist[0] = 1;
    int errorFlag = 0;
    minuit->mnexcm("SET ERR", arglist, 1, errorFlag);

    // Установка начальных значений и шагов
    static double initialValues[4] = {1.0, 530.0, 10.0, 1.0};
    static double steps[4] = {0.1, 0.1, 0.1, 0.1};

    for (int i = 0; i < 4; i++) {
        minuit->mnparm(i, (i == 3 ? "constant background" : (i == 2 ? "sigma" : (i == 1 ? "mean" : "amplitude"))),
                        initialValues[i], steps[i], 0, 0, errorFlag);
    }

    // Минимизация
    arglist[0] = 500; // максимальное количество итераций
    arglist[1] = 0.1; // точность
    minuit->mnexcm("MIGRAD", arglist, 2, errorFlag);


    double minValue, edm, errDef;
    int numVars, numParams, status;
    minuit->mnstat(minValue, edm, errDef, numVars, numParams, status);
    
    std::cout << "Minimum value: " << minValue << std::endl;


    double amplitude, mean, sigma, background;
    double amplitudeError, meanError, sigmaError, backgroundError;
    
    for (int i = 0; i < 4; i++) {
        minuit->GetParameter(i, (i == 0 ? amplitude : (i == 1 ? mean : (i == 2 ? sigma : background))), 
                             (i == 0 ? amplitudeError : (i == 1 ? meanError : (i == 2 ? sigmaError : backgroundError))));
    }


    auto *canvas = new TCanvas("Canvas", "Fitting Results", 1200, 600);
    auto *fitFunctionSignal = new TF1("Fit Function Signal", CombinedFunction, 500, 600, 4);
    fitFunctionSignal->SetParameters(amplitude, mean, sigma, background);

    auto *fitFunctionNoise = new TF1("Fit Function Noise", "[0]", 500, 600);
    fitFunctionNoise->SetParameter(0, background);


    canvas->Divide(2, 1); 
    canvas->cd(1); 
    signalHistogram->SetLineWidth(1);
    signalHistogram->Draw("E");
    fitFunctionSignal->Draw("SAME");

    canvas->cd(2); 
    noiseHistogram->SetLineWidth(1); 
    noiseHistogram->Draw("E");
    fitFunctionNoise->Draw("SAME");

    std::cout << "Chi-squared: " << chiSquareValue << std::endl;
    std::cout << "Estimated number of events: " << static_cast<int>(amplitude)
              << " ± " << static_cast<int>(amplitudeError) << std::endl;


    auto *outputFile = new TFile("fitting_results.root", "RECREATE");
    signalHistogram->Write();
    noiseHistogram->Write();
    fitFunctionSignal->Write();
    fitFunctionNoise->Write();
    outputFile->Close();
}