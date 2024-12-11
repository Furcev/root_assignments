#include "TBranch.h"
#include <TTree.h>
#include <TFile.h>
#include <ROOT/TSeq.hxx>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <iostream>

void analyzePhotonData()
{
    TString inputDataFile = "m3pimc.root"; 
    TString outputDataFile = "filtered_photon_tree.root"; 

    TFile *inputData = TFile::Open(inputDataFile); 
    TTree *dataTree = (TTree*)inputData->Get("h10"); 

    TFile *outputData = new TFile(outputDataFile, "recreate"); 
    TTree *selectedTree = dataTree->CopyTree("Isrfilter==1&&chi2_3p<30"); 


    selectedTree->SetBranchStatus("*", 0);
    

    for (auto branch : {"nph", "eph", "phiph", "thetaph"}) {
        selectedTree->SetBranchStatus(branch, 1);
    }

    TTree *finalDataTree = selectedTree->CloneTree();
    finalDataTree->SetName("FilteredPhotons");
    
    finalDataTree->Print();
    finalDataTree->Write();

    TH1F *energyHistogram = new TH1F("energyHistogram", "Photon-energy; MeV; count", 100, 0, 9);
    TF1 *fitFunction = new TF1("fitFunction", "[0]/x^[1]", 0, 9);

    fitFunction->SetParameter(1, 3);

    TCanvas *canvas = new TCanvas();
    finalDataTree->Draw("eph>>energyHistogram", "", "", finalDataTree->GetEntries(), 0);
    
    canvas->SetLogy();

    energyHistogram->Fit(fitFunction);
    energyHistogram->Write();

    Float_t maxPhotonEnergy = finalDataTree->GetMaximum("eph"); 
    Float_t minPhotonEnergy = finalDataTree->GetMinimum("eph"); 

    std::cout << std::endl
              << "Maximum : " << maxPhotonEnergy << " MeV" << std::endl
              << "Minimum : " << minPhotonEnergy << " MeV" << std::endl
              << std::endl;
}