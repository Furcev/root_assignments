#include <fstream>
#include <TMath.h>
#include <TCanvas.h>
#include <iostream>
#include <TH1D.h>
#include <TMinuit.h>
#include <TF1.h>
#include <TFile.h>

double chiSquare;
double signalCounts[100];  
double noiseCounts[100];
double noiseErrors[100];
double signalErrors[100]; 
double xValues[100];


void ReadHistogramFromFile(const std::string &filePath, TH1D *histogram) {
    std::ifstream dataFile(filePath);
    double value;

    while (dataFile >> value) {
        histogram->Fill(value);
    }
}

double gaussianFunction(double inputX, double *parameters) {
    return parameters[0] / (sqrt(2 * M_PI) * parameters[2]) * TMath::Gaus(inputX, parameters[1], parameters[2]);
}

double combinedFunction(double* inputX, double *parameters) {
    return parameters[0] / (sqrt(2 * M_PI) * parameters[2]) * TMath::Gaus(inputX[0], parameters[1], parameters[2]) + parameters[3];
}

void CalculateChiSquare(int &numParams, double *gin, double &result, double *parameters, int flag) {
    const int numberOfBins = 100;

    double totalChiSquare = 0;
    chiSquare = 0;
    for (int i = 0; i < numberOfBins; i++) {
        double delta1 = TMath::Poisson(signalCounts[i], gaussianFunction(xValues[i], parameters) + parameters[3]);
        double delta2 = TMath::Poisson(noiseCounts[i], parameters[3]);
        totalChiSquare += -2.0 * log(delta1) - 2.0 * log(delta2);
    }
    chiSquare = totalChiSquare;
    result = totalChiSquare;
}

void loadHistogramData(std::string fileName, TH1D* histogram) {
    std::ifstream inputFile(fileName);
    double inputValue;
    if (inputFile.is_open()) {
        while (inputFile >> inputValue) {
            histogram->Fill(inputValue);
        }
    inputFile.close();
}
}
void ex11(){
    auto *signalHistogram = new TH1D("Signal Histogram", "Signal", 100, 500, 600);
     auto *noiseHistogram = new TH1D("Noise Histogram", "Noise", 100, 500, 600);

    ReadHistogramFromFile("./data_1.dat", signalHistogram);
    ReadHistogramFromFile("./data_2.dat", noiseHistogram);
    
    auto *canvas = new TCanvas("canvas", "", 1200, 600);
    canvas->Divide(2, 1);
    canvas->cd(1);
    signalHistogram->Draw("e");

    for (int i = 0; i < 100; i++) {
        signalCounts[i] = signalHistogram->GetBinContent(i + 1);
        noiseCounts[i] = noiseHistogram->GetBinContent(i + 1);
    }
    
    for (int i = 0; i < 100; i++) xValues[i] = i + 500;

    TMinuit *minuit = new TMinuit(4);
    minuit->SetFCN(CalculateChiSquare);

    double argList[10];
    argList[0] = 1;
    int errorFlag = 0;
    minuit->mnexcm("SET ERR", argList, 1, errorFlag);

    static double initialValues[4] = {0.1, 550., 10., 10.};
    static double steps[4] = {0.1, 0.1, 0.1, 0.1};
    
    minuit->mnparm(0, "amplitude", initialValues[0], steps[0], 0, 0, errorFlag);
    minuit->mnparm(1, "mean", initialValues[1], steps[1], 0, 0, errorFlag);
    minuit->mnparm(2, "sigma", initialValues[2], steps[2], 0, 0, errorFlag);
    minuit->mnparm(3, "constant", initialValues[3], steps[3], 0, 0, errorFlag);

    argList[0] = 500;
    argList[1] = 0.1;
    minuit->mnexcm("MIGRAD", argList, 2, errorFlag);

    double minValue, edm, errorDef;
    int numberOfParams, numberOfFixes, status;
    minuit->mnstat(minValue, edm, errorDef, numberOfParams, numberOfFixes, status);
    minuit->mnprin(3, minValue);
    
    double amplitude, mean, sigma, background;
    double ampError, meanError, sigmaError, constError;
    minuit->GetParameter(0, amplitude, ampError);
    minuit->GetParameter(1, mean, meanError);
    minuit->GetParameter(2, sigma, sigmaError);
    minuit->GetParameter(3, background, constError);

    auto fitFunction1 = new TF1("fitFunc1", "combinedFunction", 500, 600, 4);
    fitFunction1->SetParameter(0, amplitude);
    fitFunction1->SetParameter(1, mean);
    fitFunction1->SetParameter(2, sigma);
    fitFunction1->SetParameter(3, background);
    fitFunction1->Draw("Same");

    auto fitFunction2 = new TF1("fitFunc2", "[3]", 500, 600);
    fitFunction2->SetParameter(3, background);
    canvas->cd(2);
    noiseHistogram->Draw("E");
    fitFunction2->Draw("Same");

    std::cout << "Chi-squared: " << chiSquare << std::endl;
    std::cout << "Estimated number of events: " << static_cast<int>(amplitude)
              << " ± " << static_cast<int>(ampError);

    TFile* outputFile = new TFile("fitting_results.root", "RECREATE");
    signalHistogram->Write();
    noiseHistogram->Write();
    fitFunction1->Write();
    fitFunction2->Write();
    outputFile->Close();
}