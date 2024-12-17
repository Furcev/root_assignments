#include <iostream>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TMath.h>
#include <TGraphErrors.h>

using namespace TMath;

void analyze_photon_candidates(const char* inputFileName, const char* outputDir) {
    // Открываем входной файл
    TFile *inputFile = TFile::Open(inputFileName);
    TTree *tree = (TTree*)inputFile->Get("MyTree");

    // Переменные для чтения данных из дерева
    Int_t numPhotons;
    Float_t energyPhoton[300], thetaPhoton[300], phiPhoton[300];

    tree->SetBranchAddress("nph", &numPhotons);
    tree->SetBranchAddress("eph", energyPhoton);
    tree->SetBranchAddress("thetaph", thetaPhoton);
    tree->SetBranchAddress("phiph", phiPhoton);

    // Создаем выходной файл и гистограммы
    TFile *outputFile = new TFile(Form("%s/output.root", outputDir), "RECREATE");
    
    // Гистограммы для полярного и азимутального углов
    TH1F *histPolarAngle = new TH1F("histPolarAngle", "Polar Angle Distribution; Polar Angle (rad); Number of Candidates", 36, 0, TMath::Pi());
    TH1F *histAzimuthalAngle = new TH1F("histAzimuthalAngle", "Azimuthal Angle Distribution; Azimuthal Angle (rad); Number of Candidates", 36, 0, 2 * TMath::Pi());

    Int_t totalCandidates = 0;

    // Цикл по событиям
    Long64_t totalEntries = tree->GetEntries();
    for (Long64_t i = 0; i < totalEntries; i++) {
        tree->GetEntry(i);
        
        TLorentzVector photon1(0, 0, 0, 0), photon2(0, 0, 0, 0);
        
        // Цикл по парам фотонов
        for (Int_t j = 0; j < (numPhotons - 1); j++) {
            for (Int_t k = j + 1; k < numPhotons; k++) {
                photon1.SetPx(energyPhoton[j] * Sin(thetaPhoton[j]) * Cos(phiPhoton[j]));
                photon1.SetPy(energyPhoton[j] * Sin(thetaPhoton[j]) * Sin(phiPhoton[j]));
                photon1.SetPz(energyPhoton[j] * Cos(thetaPhoton[j]));
                photon1.SetE(energyPhoton[j]);

                photon2.SetPx(energyPhoton[k] * Sin(thetaPhoton[k]) * Cos(phiPhoton[k]));
                photon2.SetPy(energyPhoton[k] * Sin(thetaPhoton[k]) * Sin(phiPhoton[k]));
                photon2.SetPz(energyPhoton[k] * Cos(thetaPhoton[k]));
                photon2.SetE(energyPhoton[k]);

                Float_t invariantMass = (photon1 + photon2).M(); // Вычисление инвариантной массы
                
                // Проверка инвариантной массы
                if (invariantMass >= 0.1 && invariantMass <= 0.2) {
                    // Заполнение гистограмм
                    histPolarAngle->Fill(thetaPhoton[j]);
                    histPolarAngle->Fill(thetaPhoton[k]);
                    histAzimuthalAngle->Fill(phiPhoton[j]);
                    histAzimuthalAngle->Fill(phiPhoton[k]);
                    totalCandidates++;
                }
            }
        }
    }

    std::cout << "Total candidates found: " << totalCandidates << std::endl;

    // Создание графиков с ошибками
    TGraphErrors *graphPolarAngle = new TGraphErrors();
    TGraphErrors *graphAzimuthalAngle = new TGraphErrors();

    for (int bin = 1; bin <= histPolarAngle->GetNbinsX(); bin++) {
        double x = histPolarAngle->GetBinCenter(bin);
        double y = histPolarAngle->GetBinContent(bin);
        double error = TMath::Sqrt(y); // Ошибка по квадратному корню из количества
        graphPolarAngle->SetPoint(bin - 1, x, y);
        graphPolarAngle->SetPointError(bin - 1, 0, error);
    }

    for (int bin = 1; bin <= histAzimuthalAngle->GetNbinsX(); bin++) {
        double x = histAzimuthalAngle->GetBinCenter(bin);
        double y = histAzimuthalAngle->GetBinContent(bin);
                double error = TMath::Sqrt(y); // Ошибка по квадратному корню из количества
        graphAzimuthalAngle->SetPoint(bin - 1, x, y);
        graphAzimuthalAngle->SetPointError(bin - 1, 0, error);
    }

    // Создание канваса и рисование графиков
    TCanvas *canvas = new TCanvas("canvas", "Photon Candidates Analysis", 1200, 600);
    canvas->Divide(2, 1);

    canvas->cd(1);
    graphPolarAngle->SetMarkerStyle(21); // Квадратные маркеры
    graphPolarAngle->SetMarkerSize(0.5);
    graphAzimuthalAngle->SetLineColor(kRed);
    graphPolarAngle->SetTitle("Polar Angle Distribution; Polar Angle (rad); Number of Candidates");
    graphPolarAngle->Draw("AP");
    
    canvas->cd(2);
    graphAzimuthalAngle->SetMarkerStyle(21); // Квадратные маркеры
    graphAzimuthalAngle->SetMarkerSize(0.5);
    graphAzimuthalAngle->SetLineColor(kRed);
    graphAzimuthalAngle->SetTitle("Azimuthal Angle Distribution; Azimuthal Angle (rad); Number of Candidates");
    graphAzimuthalAngle->Draw("AP");

    // Добавление легенды
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(graphPolarAngle, "Candidates in Polar Angle", "p");
    legend->AddEntry(graphAzimuthalAngle, "Candidates in Azimuthal Angle", "p");
    legend->Draw();

    // Сохранение графиков в файл
    canvas->SaveAs(Form("%s/pion_candidates_analysis.png", outputDir));
}