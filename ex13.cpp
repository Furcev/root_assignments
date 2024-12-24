#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"

void analyze_data() {

    TH1D *hist = new TH1D("hist", "Histogram of m3piJPSI", 41, 3.0, 3.2); // 41 bins from 3 GeV to 3.2 GeV


    std::ifstream file("m3piJPSI_cut.dat");
    double value;
    while (file >> value) {
        hist->Fill(value);
    }
    file.close();


    TCanvas *canvas = new TCanvas("canvas", "Canvas", 800, 600);
    hist->Draw();

    // Определяем функцию Брейт-Вигнера + линейная функция
    TF1 *bw_linear = new TF1("bw_linear", "[0]*TMath::BreitWigner(x, [1], [2]) + [3]*x + [4]", 3.0, 3.2);
    bw_linear->SetParameters(1, 3.0969, 0.092, 0.1, 0); 


    //hist->Fit(bw_linear, "R");


    double mass = bw_linear->GetParameter(1);
    double width = bw_linear->GetParameter(2);
    double slope = bw_linear->GetParameter(3);
    double intercept = bw_linear->GetParameter(4);


    std::cout << "Mass: " << mass << " GeV" << std::endl;
    std::cout << "Width: " << width << " GeV" << std::endl;
    std::cout << "Slope: " << slope << std::endl;
    std::cout << "Intercept: " << intercept << std::endl;


    TF1 *double_gauss = new TF1("double_gauss", "[0]*TMath::Gaus(x, [1], [2]) + [3]*TMath::Gaus(x, [4], [5])", 3.0, 3.2);
    double_gauss->SetParameters(1, 3.0969, 0.01, 0.5, 3.1, 0.02); 


    hist->Fit(double_gauss, "R");


    double mean1 = double_gauss->GetParameter(1);
    double sigma1 = double_gauss->GetParameter(2);
    double mean2 = double_gauss->GetParameter(4);
    double sigma2 = double_gauss->GetParameter(5);


    std::cout << "Mean1: " << mean1 << " GeV" << std::endl;
    std::cout << "Sigma1: " << sigma1 << " GeV" << std::endl;
    std::cout << "Mean2: " << mean2 << " GeV" << std::endl;
    std::cout << "Sigma2: " << sigma2 << " GeV" << std::endl;


    canvas->cd();
    hist->Draw();
    //bw_linear->Draw("same");
    //double_gauss->Draw("same");

    canvas->SaveAs("fit_results.png"); 
}