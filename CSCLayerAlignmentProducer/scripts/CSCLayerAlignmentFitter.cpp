#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#include "TMath.h"
#include "TMinuit.h"
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TF1.h"

// Global variables for storing fit results and errors
std::vector<double> fitResults = {0., 0., 0., 0.};
std::vector<double> fitErrors = {0., 0., 0., 0.};

// Global variables for TTree data
TTree *dataTree;
Float_t residual, localPoint[3], globalPoint[3];
Long64_t numEvents;

// Function to calculate the log of a pure Gaussian
double logPureGaussian(double residual, double center, double sigma) {
    sigma = fabs(sigma);
    static const double cgaus = 0.5 * log(2. * M_PI);
    return (-pow(residual - center, 2) * 0.5 / sigma / sigma) - cgaus - log(sigma);
}

// Function to calculate the residual
double calculateResidual(double deltaX, double deltaY, double deltaPhiZ, double trackX, double trackY, double radius) {
    return deltaX - (trackX / radius - 3. * pow(trackX / radius, 3)) * deltaY - trackY * deltaPhiZ;
}

// Minuit function for fitting
void fitFunction(int &numParams, double *gradients, double &functionValue, double *params, int flag) {
    const double deltaX = params[0];
    const double deltaY = params[1];
    const double deltaPhiZ = params[2];
    const double sigma = params[3];

    functionValue = 0.;
    for (Long64_t i = 0; i < numEvents; i++) {
        dataTree->GetEntry(i);
        double trackX = localPoint[0];
        double trackY = localPoint[1];
        double radius = sqrt(globalPoint[0] * globalPoint[0] + globalPoint[1] * globalPoint[1]);
        double residualPeak = calculateResidual(deltaX, deltaY, deltaPhiZ, trackX, trackY, radius);
        functionValue += -1. * logPureGaussian(residual, residualPeak, sigma);
    }
}

// Function to perform the fit using TMinuit
void performFit(bool fitDeltaX, bool fitDeltaY, bool fitDeltaPhiZ) {
    TMinuit minuit(4);
    minuit.SetFCN(fitFunction);
    double initialParams[4] = {0., 0., 0., 0.5};
    minuit.DefineParameter(0, "deltaX", initialParams[0], 0.1, 0, 0);
    minuit.DefineParameter(1, "deltaY", initialParams[1], 0.1, 0, 0);
    minuit.DefineParameter(2, "deltaPhiZ", initialParams[2], 0.001, 0, 0);
    minuit.DefineParameter(3, "sigma", initialParams[3], 0.01, 0, 0);
    minuit.FixParameter(3); // Fix sigma

    if (!fitDeltaX) minuit.FixParameter(0);
    if (!fitDeltaY) minuit.FixParameter(1);
    if (!fitDeltaPhiZ) minuit.FixParameter(2);

    double argList[10] = {0};
    int errorFlag = 0;

    argList[0] = 0.5;
    minuit.mnexcm("SET ERR", argList, 1, errorFlag);

    argList[0] = 2;
    minuit.mnexcm("SET STR", argList, 1, errorFlag);

    argList[0] = 50000;
    minuit.mnexcm("MIGRAD", argList, 1, errorFlag);

    if (errorFlag != 0) {
        std::cout << "Fit failed, trying again with MIGRAD" << std::endl;
        minuit.mnexcm("MIGRAD", argList, 1, errorFlag);
    }

    Double_t minVal, errDist, errDef;
    Int_t numFreeParams, numParamsTotal, fitStatus;
    minuit.mnstat(minVal, errDist, errDef, numFreeParams, numParamsTotal, fitStatus);

    if (fitStatus != 3) {
        minuit.mnexcm("HESSE", argList, 0, errorFlag);
    }

    for (int i = 0; i < 3; i++) {
        double value, error;
        minuit.GetParameter(i, value, error);
        fitResults[i] = value;
        fitErrors[i] = error;
    }
}

int main() {
    // Configuration parameters
    const char* inputFileName = "/eos/user/m/mkizilov/crab3_out/ME11_ana/2023-12-13/Muon/2023-12-13/231213_200534/0000/merged_file.root";
    const char* treeName = "ME11ana/Inner_Prop";
    const char* residualName = "RdPhi";
    const char* outputPrefix = "CSC_layer_al_2023-12_13_1DOF_2022D";

    // Cuts for filtering the TTree
    const char* baseCuts = "muon_pt > 30 && abs(RdPhi) < 999 && has_fidcut";
    const char* positiveChargeCut = "muon_charge > 0";
    const char* negativeChargeCut = "muon_charge < 0";

    // Fit parameters
    bool fitDeltaX = true;
    bool fitDeltaY = false;
    bool fitDeltaPhiZ = false;

    // Configuration for layer or chamber level fit
    bool fitByLayer = true;
    int numCuts = 2;

    int maxLayer = fitByLayer ? 7 : 2;

    // Load input TTree
    TFile inputFile(inputFileName);
    TTree *inputTree = static_cast<TTree*>(inputFile.Get(treeName));
    TFile tmpFile("tmp1.root", "RECREATE");
    TTree *filteredTree = inputTree->CopyTree(baseCuts);
    inputFile.Close();

    for (int cutLevel = 1; cutLevel < numCuts; cutLevel++) {
        if (cutLevel > 1) {
            int totalEntries = filteredTree->GetEntries();
            int halfEntries = totalEntries / 2;
            filteredTree = filteredTree->CloneTree(halfEntries);
        }

        std::ofstream outputFile(Form("%s.csv", outputPrefix));

        for (int region = -1; region < 2; region += 2) {
            for (int chamber = 0; chamber < 36; chamber++) {
                for (int layer = 1; layer < maxLayer; layer++) {
                    int detectorId = region * (chamber + 101);

                    TFile tmpFile2("tmp2.root", "RECREATE");
                    TTree *detectorTree;
                    if (fitByLayer) {
                        detectorTree = filteredTree->CopyTree(Form("rechit_detId==%d && prop_location[4] == %d", detectorId, layer));
                    } else {
                        detectorTree = filteredTree->CopyTree(Form("rechit_detId==%d", detectorId));
                    }

                    if (detectorTree->GetEntries() < 28) continue;

                    TH1F hist("hist", "Histogram", 100, -20, 20);
                    detectorTree->Project("hist", residualName, "");
                    TF1 fitFunc("fitFunc", "gaus", -2, 2);
                    fitFunc.SetParLimits(1, -2, 2);
                    fitFunc.SetParLimits(2, 0, 2);
                    hist.Fit("fitFunc", "R");
                    float fitMean = fitFunc.GetParameter(1);
                    float fitStdDev = fitFunc.GetParameter(2);

                    dataTree = detectorTree->CopyTree(Form("RdPhi <= (%f + (1.6*%f)) && RdPhi >= (%f - (1.6*%f))", fitMean, fitStdDev, fitMean, fitStdDev));

                    if (dataTree->GetEntries() == 0) {
                        outputFile << detectorId << ", " << 0 << ", " << 0 << ", " << 0 << ", " << 0 << ", " << 0 << ", " << 0 << ", " << 0 << "\n";
                        continue;
                    }

                    dataTree->SetBranchAddress("RdPhi", &residual);
                    dataTree->SetBranchAddress("prop_LP", localPoint);
                    dataTree->SetBranchAddress("prop_GP", globalPoint);

                    numEvents = dataTree->GetEntries();
                    performFit(fitDeltaX, fitDeltaY, fitDeltaPhiZ);

                    double deltaX = fitResults[0];
                    double deltaY = fitResults[1];
                    double deltaPhiZ = fitResults[2];

                    if (fitByLayer) {
                        outputFile << detectorId << layer << ", " << deltaX << ", " << deltaY << ", " << 0.0 << ", " << 0.0 << ", " << 0.0 << ", " << deltaPhiZ << ", " << numEvents << "\n";
                    } else {
                        outputFile << detectorId << ", " << deltaX << ", " << deltaY << ", " << 0.0 << ", " << 0.0 << ", " << 0.0 << ", " << deltaPhiZ << ", " << numEvents << "\n";
                    }
                }
            }
        }
        outputFile.close();
    }
}
