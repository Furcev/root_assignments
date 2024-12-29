#include <TH1D.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <iostream>
#include <TF1.h>
#include <TLorentzRotation.h>
#include <TCanvas.h>


double det_radius = 30; 
double det_length = 50;  
double avg_dist = 6; 
double tot_energy = 1020; 
double pion_m = 139.57039; 
double kaon_m = 497.611; 


void fillHistograms(const unsigned int evtCount, TH1D* hPionLen, TH1D* hPionTheta, TH1D* hPionPhi) {
    TRandom3 rng;
    rng.SetSeed(0);
    for (unsigned int i = 0; i < evtCount; ++i) {
        TLorentzRotation rot;
        double kaonE = (tot_energy - kaon_m) * rng.Rndm() - kaon_m;
        TLorentzVector kVec(0., 0., 0., kaon_m);
        double kMom = TMath::Sqrt(kaonE * kaonE - kaon_m * kaon_m);
        double azAngle = 2 * TMath::Pi() * rng.Rndm();
        double cosTheta = 2 * rng.Rndm() - 1;
        double sinTheta = TMath::Sqrt(1 - cosTheta * cosTheta);
        double pathLen = -TMath::Log(rng.Rndm()) * avg_dist;
        hPionLen->Fill(pathLen);
        //P+
        double pPionE = (kaon_m - pion_m) * rng.Rndm() - pion_m;
        double pPionMom = TMath::Sqrt(pPionE * pPionE - pion_m * pion_m);
        double pPionPhi = 2 * TMath::Pi() * rng.Rndm();
        double pPionCosTheta = 2 * rng.Rndm() - 1;
        double pPionSinTheta = TMath::Sqrt(1 - pPionCosTheta * pPionCosTheta);
        //P-
        double mPionE = pion_m - pPionE;
        double mPionMom = TMath::Sqrt(mPionE * mPionE - pion_m * pion_m);
        double mPionPhi = 2 * TMath::Pi() - pPionPhi;
        double mPionCosTheta = TMath::Cos(TMath::Pi() / 2 + TMath::ACos(pPionCosTheta));
        double mPionSinTheta = TMath::Sqrt(1 - mPionCosTheta * mPionCosTheta);

        TLorentzVector pPionVec(
            pPionMom * pPionSinTheta * TMath::Cos(pPionPhi),
            pPionMom * pPionSinTheta * TMath::Sin(pPionPhi),
            pPionMom * pPionCosTheta,
            pPionE
        );
        TLorentzVector mPionVec(mPionMom * mPionSinTheta * TMath::Cos(mPionPhi),mPionMom * mPionSinTheta * TMath::Sin(mPionPhi),mPionMom * mPionCosTheta, mPionE);
        pPionVec = rot.VectorMultiplication(pPionVec);
        mPionVec = rot.VectorMultiplication(mPionVec);

        hPionPhi->Fill(azAngle);
        hPionLen->Fill(pathLen);
        hPionTheta->Fill(TMath::ACos(cosTheta));
    }
}

void generateHistograms(unsigned int evtCount) {
    TH1D* hPionLen = new TH1D("hPionLen", "Pion Length Distribution", 50, 0, TMath::Sqrt(det_radius * det_radius + det_length * det_length / 4));
    TH1D* hPionTheta = new TH1D("hPionTheta", "Pion Theta Distribution", 50, 0, TMath::Pi() / 2);
    TH1D* hPionPhi = new TH1D("hPionPhi", "Pion Phi Distribution", 50, 0, TMath::Pi());
    fillHistograms(evtCount, hPionLen, hPionTheta, hPionPhi);
    TCanvas* canvas = new TCanvas("canvas", "Histogram Overview", 1200, 700);
    canvas->Divide(2, 2); 
    canvas->cd(1);
     hPionPhi->Draw();
    canvas->cd(3);
    hPionTheta->Draw();
    canvas->cd(4);
    hPionLen->Draw();
}