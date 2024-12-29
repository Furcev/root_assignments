#include <iostream>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TMath.h>

using namespace TMath;

void analyze_events() {

    TFile *inputFile = TFile::Open("newtree.root");
    TTree *tree = (TTree*)inputFile->Get("MyTree");


    int numPhotons;
    float energyPhoton[300], thetaPhoton[300], phiPhoton[300]; 
    
    tree->SetBranchAddress("nph", &numPhotons);
    tree->SetBranchAddress("eph", energyPhoton);
    tree->SetBranchAddress("thetaph", thetaPhoton);
    tree->SetBranchAddress("phiph", phiPhoton);


    TFile *outputFile = new TFile("new.root", "RECREATE");
    TH1F *histInvariantMass = new TH1F("histInvariantMass", " Mass of p0; Mass (MeV/c^{2}); Number of Events", 100, 100, 200);
    TH1F *histAngleBetweenPhotons = new TH1F("histAngleBetweenPhotons", "Angle Between #pi^{0} Candidates; Angle (rad); Number of Events", 100, 0, TMath::Pi());

    int totalCandidates = 0;


    Long64_t totalEntries = tree->GetEntries();
    for (Long64_t i = 0; i < totalEntries; i++) {
        tree->GetEntry(i);
        
        Float_t invariantMassArray[300], angleArray[300]; 
        int foundCandidates = 0;

        TLorentzVector photon1(0, 0, 0, 0), photon2(0, 0, 0, 0);
        
        for (int j = 0; j < (numPhotons - 1); j++) {
            for (int k = j + 1; k < numPhotons; k++) {
                photon1.SetPx(energyPhoton[j] * Sin(thetaPhoton[j]) * Cos(phiPhoton[j]));
                photon1.SetPy(energyPhoton[j] * Sin(thetaPhoton[j]) * Sin(phiPhoton[j]));
                photon1.SetPz(energyPhoton[j] * Cos(thetaPhoton[j]));
                photon1.SetE(energyPhoton[j]);

                photon2.SetPx(energyPhoton[k] * Sin(thetaPhoton[k]) * Cos(phiPhoton[k]));
                photon2.SetPy(energyPhoton[k] * Sin(thetaPhoton[k]) * Sin(phiPhoton[k]));
                photon2.SetPz(energyPhoton[k] * Cos(thetaPhoton[k]));
                photon2.SetE(energyPhoton[k]);

                float invariantMass = (photon1 + photon2).M(); 
                float angle = photon1.Vect().Angle(photon2.Vect()); 


                if (invariantMass >= 0.1 && invariantMass <= 0.2) {
                    invariantMassArray[foundCandidates] = invariantMass;
                    angleArray[foundCandidates] = angle;
                    foundCandidates++;
                }
            }
        }
        

        if (foundCandidates == 2) {
            histInvariantMass->Fill(invariantMassArray[0] * 1000);
            histInvariantMass->Fill(invariantMassArray[1] * 1000);
            histAngleBetweenPhotons->Fill(angleArray[0]);
            histAngleBetweenPhotons->Fill(angleArray[1]);
            totalCandidates++;
        }
    }

    std::cout << "candidate: " << totalCandidates << std::endl;

    TCanvas *canvas = new TCanvas("canvas", "#pi^{0} Candidates Analysis", 1200, 600);
    canvas->Divide(2, 1);

    canvas->cd(1);
    histAngleBetweenPhotons->Draw();
    histAngleBetweenPhotons->Write();
    
    canvas->cd(2);
    histInvariantMass->Draw();
    histInvariantMass->Write();


}