#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TFile.h>
#include <iostream>
#include <fstream>

const int N_PARAMETERS = 4;

void loadHistogram(TH1D *histogram, const TString& filename) {
    std::ifstream inputFile(filename.Data(), std::ios::in);
    
    if (inputFile.is_open()) 
    {
        double value; 
        while (inputFile >> value)
            histogram->Fill(value);
        
        inputFile.close();
    } 
    else 
    {
        std::cerr << "Error: File " << filename << " cannot be opened!\n";
        exit(-1);
    }
}

double modelFunction(double *x, double *params) {
    return params[0] * TMath::Poisson(x[0], 0.2) 
         + params[1] * TMath::BreitWigner(x[0], params[2], params[3]);
}

void analyze_data() {
  
    gStyle->SetOptFit(1111);
    
    TH1D *dataHistogram = new TH1D("dataHist", "Data Distribution from task10Nov.dat;Value;Counts", 100, 0, 10);

    loadHistogram(dataHistogram, "task10Nov.dat");


    TF1 *fitFunction = new TF1("fitFunc", modelFunction, 0, 10, N_PARAMETERS);
    

    fitFunction->SetParameter(0, 1); 
    fitFunction->SetParLimits(0, 0, 1000);
    
    fitFunction->SetParameter(1, 1); 
    fitFunction->SetParLimits(1, 0, 1000);
    
    fitFunction->SetParameter(2, 5); 
    fitFunction->SetParLimits(2, 0, 10);
    
    fitFunction->SetParameter(3, 0.5); 
    fitFunction->SetParLimits(3, 0.01, 5);

    TCanvas *canvas = new TCanvas("canvas", "Fit Results", 1200, 600);
    

    dataHistogram->Draw("HIST");


    TFitResultPtr fitResult = dataHistogram->Fit(fitFunction, "SMPE", "", 0, 10);


    double chiSquare = fitFunction->GetChisquare();
    
    if (chiSquare <= 150) 
    {
        std::cout << "Fit successful with chi2 = " << chiSquare << std::endl;
        fitFunction->Draw("same");
    } 
    else 
    {
        std::cout << "Fit failed to meet chi2 criterion (chi2 = " << chiSquare << ")" << std::endl;
    }


    canvas->SaveAs("fit_results.png");
}