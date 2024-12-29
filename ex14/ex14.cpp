#include <iostream>
#include <fstream>
#include <TCanvas.h>
#include <TH1F.h>
#include <TVirtualFFT.h>

constexpr int dataSize = 500;

double computeTimeAtThreshold(TH1F* histogram, double thresholdRatio) {
    double peakValue = histogram->GetMaximum();
    double thresholdValue = thresholdRatio * peakValue;

    for (int binIndex = 1; binIndex <= histogram->GetNbinsX(); ++binIndex) {
        if (histogram->GetBinContent(binIndex) >= thresholdValue) {
            return histogram->GetBinCenter(binIndex);
        }
    }
    return -1;
}

void loadHistogramData(const std::string& filePath, TH1F* histogram) {
    std::ifstream inputFile(filePath);
    double value;
    int index = 0;

    if (inputFile.is_open()) {
        while (inputFile >> value && index < histogram->GetNbinsX()) {
            histogram->SetBinContent(index + 1, value);
            ++index;
        }
    inputFile.close();
}
}
void ex14() {    
    auto signalCanvas = new TCanvas("signalAnalysis", "Signal Analysis", 1200, 600);
    signalCanvas->Divide(2, 1);

    signalCanvas->cd(1);
    TH1F* rawSignalHistogram = new TH1F("rawSignal", "Raw Signal", dataSize, 0., 500.);
    loadHistogramData("dataFFT.dat", rawSignalHistogram);
    rawSignalHistogram->Draw();

    signalCanvas->cd(2);
    TH1* fftHistogram = rawSignalHistogram->FFT(nullptr, "MAG"); 
    fftHistogram->SetTitle("FFT Magnitude");
    fftHistogram->Draw();

    signalCanvas->cd(3);
    Double_t* realPart = new Double_t[dataSize];
    Double_t* imaginaryPart = new Double_t[dataSize];
    TVirtualFFT* currentFFT = TVirtualFFT::GetCurrentTransform();
    currentFFT->GetPointsComplex(realPart, imaginaryPart);

    for (int idx = 0; idx < dataSize; ++idx) {
        if (idx > 10 && idx < 490) { 
            realPart[idx] = 0;
            imaginaryPart[idx] = 0;
        }
    }

    int dims = dataSize; 
    TVirtualFFT* inverseFFT = TVirtualFFT::FFT(1, &dims, "C2R");
    inverseFFT->SetPointsComplex(realPart, imaginaryPart);
    inverseFFT->Transform();
    TH1* filteredHistogram = TH1::TransformHisto(inverseFFT, nullptr, "RE"); 

    auto filteredSignalHistogram = new TH1F("filteredSignal", "Filtered Signal", dataSize, 0, 500);
    for (int idx = 0; idx < dataSize; idx++) {
        filteredSignalHistogram->SetBinContent(idx + 1, filteredHistogram->GetBinContent(idx + 1) / dataSize); 
    }
    filteredSignalHistogram->Draw();

    double timeThreshold = computeTimeAtThreshold(rawSignalHistogram, 0.2);
    std::cout << "Time at which the signal reached 20% of maximum amplitude: " << timeThreshold << " ns" << std::endl;

    signalCanvas->Update();
}