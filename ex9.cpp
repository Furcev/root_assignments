#include <iostream>
#include <TTree.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <vector>

void analyze_photons(){
    TFile *inputFile = TFile::Open("newroot.root", "READ");

    TFile *outputFile = new TFile("new.root", "RECREATE");
    
    TTree *tree = (TTree *)inputFile->Get("MyTree");

    int numPhotons;
    float energyPhoton[40], thetaPhoton[40], phiPhoton[40];

    std::vector<float> thetaPhotonGraph(18, 0);
    std::vector<float> phiPhotonGraph(36, 0);

    tree->SetBranchAddress("nph", &numPhotons);
    tree->SetBranchAddress("eph", energyPhoton);
    tree->SetBranchAddress("thetaph", thetaPhoton);
    tree->SetBranchAddress("phiph", phiPhoton);

    // График для theta
    TCanvas *canvasTheta = new TCanvas("canvasTheta", "Candidates vs #theta", 600, 500);
    TGraphErrors *graphTheta = new TGraphErrors();
    graphTheta->SetTitle("Candidates vs #theta");

    // График для phi
    TCanvas *canvasPhi = new TCanvas("canvasPhi", "Candidates vs #phi", 600, 500);
    TGraphErrors *graphPhi = new TGraphErrors();
    graphPhi->SetTitle("Candidates vs #phi");

    // Главный цикл обработки событий
    for (int i = 0; i < tree->GetEntries(); i++) {
        tree->GetEntry(i);
        for (int j = 0; j < numPhotons - 1; j++) {
            for (int k = j + 1; k < numPhotons; k++) {
                TLorentzVector photon1(energyPhoton[j], 0, 0, energyPhoton[j]);
                photon1.SetTheta(thetaPhoton[j]);
                photon1.SetPhi(phiPhoton[j]);

                TLorentzVector photon2(energyPhoton[k], 0, 0, energyPhoton[k]);
                photon2.SetTheta(thetaPhoton[k]);
                photon2.SetPhi(phiPhoton[k]);

                TLorentzVector photonSum = photon1 + photon2;

                TVector3 vector1(sin(thetaPhoton[j]) * cos(phiPhoton[j]), sin(thetaPhoton[j]) * sin(phiPhoton[j]), cos(thetaPhoton[j]));
                TVector3 vector2(sin(thetaPhoton[k]) * cos(phiPhoton[k]), sin(thetaPhoton[k]) * sin(phiPhoton[k]), cos(thetaPhoton[k]));
                float angle = vector1.Angle(vector2);

                double invariantMass = sqrt(2. * energyPhoton[j] * energyPhoton[k] * (1. - cos(angle)));
                if (0.1 <= invariantMass && invariantMass <= 0.2) {
                    int thetaIndex = photonSum.Theta() * 18 / TMath::Pi();
                    int phiIndex = (photonSum.Phi() + TMath::Pi()) * 18 / TMath::Pi();
                    thetaPhotonGraph[thetaIndex]++;
                    phiPhotonGraph[phiIndex]++;
                }
            }
        }
    }

    for (int i = 0; i < 18; i++) {
        graphTheta->SetPoint(i, i * 10, thetaPhotonGraph[i]);
        graphTheta->SetPointError(i, 0, sqrt(thetaPhotonGraph[i]));
    }
    graphTheta->SetMarkerColor(kGreen); 
    graphTheta->SetMarkerStyle(21);
    graphTheta->SetMarkerSize(0.5); 

    for (int i = 0; i < 36; i++) {
        graphPhi->SetPoint(i, i * 10, phiPhotonGraph[i]);
        graphPhi->SetPointError(i, 0, sqrt(phiPhotonGraph[i]));
    }
    graphPhi->SetMarkerColor(kMagenta); 
    graphPhi->SetMarkerStyle(21);
    graphPhi->SetMarkerSize(0.5); 

    canvasTheta->cd();
    graphTheta->Draw("APC"); 
    canvasTheta->BuildLegend(0.48, 0.8, 0.9, 0.9);

    canvasPhi->cd();
    graphPhi->Draw("APC"); 
}