#include <iostream>
#include <fstream>
#include <string>

#include "TMath.h"
#include "TMinuit.h"
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TF1.h"

std::vector<double> mResult = {0., 0., 0., 0.};
std::vector<double> mError = {0., 0., 0., 0.};

TTree *tt;
//Float_t mResidual, mTrackX, mTrackY, mR;
Float_t mResidual, mLocalPoint[3], mGlobalPoint[3];
Long64_t mEvents;

double MuonResidualsFitter_logPureGaussian(double residual, double center, double sigma) {
  sigma = fabs(sigma);
  static const double cgaus = 0.5 * log( 2.*M_PI );
  return (-pow(residual - center, 2) *0.5 / sigma / sigma) - cgaus - log(sigma);
}

double getResidual(double delta_x, double delta_y, double delta_phiz, double track_x, double track_y, double R) {
  return delta_x - (track_x/R - 3.*pow(track_x/R, 3)) * delta_y - track_y * delta_phiz;
}





void MuonResiduals3DOFFitter_FCN(int &npar, double *gin, double &fval, double *par, int iflag) {
  const double dx = par[0];
  const double dy = par[1];
  const double dphiz = par[2];
  const double sig = par[3];

  fval = 0.;
  for (Long64_t i=0;i<mEvents;i++) {
    tt->GetEntry(i);
    double residual = mResidual; double trackX = mLocalPoint[0]; double trackY = mLocalPoint[1]; double R = pow(pow(mGlobalPoint[0],2) + pow(mGlobalPoint[1],2),0.5);
    double residpeak = getResidual(dx, dy, dphiz, trackX, trackY, R);
    fval += -1.*MuonResidualsFitter_logPureGaussian(residual, residpeak, sig);
  }
}

void doFit(bool doDx, bool doDy, bool doDphiz) {
  TMinuit mfit(4);
  mfit.SetFCN(MuonResiduals3DOFFitter_FCN);
  double par[4] = {0., 0., 0., 0.5};
  mfit.DefineParameter(0, "dx", par[0], 0.1, 0, 0);
  mfit.DefineParameter(1, "dy", par[1], 0.1, 0, 0);
  mfit.DefineParameter(2, "dphiz", par[2], 0.001, 0, 0);
  mfit.DefineParameter(3, "sig", par[3], 0.01, 0, 0);
  mfit.FixParameter(3);
  if (!doDx) mfit.FixParameter(0);
  if (!doDy) mfit.FixParameter(1);
  if (!doDphiz) mfit.FixParameter(2);

  double arglist[10];
  int ierflg;
  int smierflg;

  for (int i = 0;  i < 10;  i++) arglist[i] = 0.;
  arglist[0] = 0.5;
  ierflg = 0;
  smierflg = 0;
  mfit.mnexcm("SET ERR", arglist, 1, ierflg);
  for (int i = 0;  i < 10;  i++) arglist[i] = 0.;
  arglist[0] = 2;
  ierflg = 0;
  mfit.mnexcm("SET STR", arglist, 1, ierflg);

  bool try_again = false;
  for (int i = 0;  i < 10;  i++) arglist[i] = 0.;
  arglist[0] = 50000;
  ierflg = 0;
  mfit.mnexcm("MIGRAD", arglist, 1, ierflg);
  if (ierflg != 0) try_again = true;
  if (try_again){
    std::cout << "try again" << std::endl;
    for (int i = 0;  i < 10;  i++) arglist[i] = 0.;
    arglist[0] = 50000;
    mfit.mnexcm("MIGRAD", arglist, 1, smierflg);
  }

  Double_t fmin, fedm, errdef;
  Int_t npari, nparx, istat;
  mfit.mnstat(fmin, fedm, errdef, npari, nparx, istat);
  if (istat != 3) {
    for (int i = 0;  i < 10;  i++) arglist[i] = 0.;
    ierflg = 0;
    mfit.mnexcm("HESSE", arglist, 0, ierflg);
  }
  for (int i = 0;  i < 3;  i++){
    double v,e;
    mfit.GetParameter(i,v,e);
    mResult[i] = v;
    mError[i] = e;
  }
}

int main() {
  //////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////
  //Input root file name
  const char* input_name = "/eos/user/m/mkizilov/crab3_out/ME11_ana/2023-12-13/Muon/2023-12-13/231213_200534/0000/merged_file.root";
  //Tree name ***Make sure to use correct one***
  const char* tree_name = "ME11ana/Inner_Prop";   //"analyzer/ME11Seg_Prop" or "ME11ana/Inner_Prop" for example
  const char* Rdphi_name = "RdPhi";
  //Will only change the name of the output csv file
  const char* outname_prefix = "CSC_layer_al_2023-12_13_1DOF_2022D";

  //Cuts on full tree in first cloning step
  const char* cuts = "muon_pt > 30 && abs(RdPhi) < 999 && has_fidcut"; //n_ME11_segment == 1
  const char* positive_cut = "muon_charge > 0";
  const char* negative_cut = "muon_charge < 0";
  //Option to turn on or off 3 dof alignments and layer level vs chamber level
  bool doDx = true;
  bool doDy = false;
  bool doDphiz = false;
  //Make sure this one matches the tree!!!!!!! CSCs do not do by layer!
  bool byLayer = true;
  //Number of cuts to fit on, each cut reduces the total sample remaining by 2 (n/2, 3n/4, 7n/8, 15n/16, ...)
  //2 is the base value
  int nCuts = 2;
  //////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////

  int max_layer = 2;
  if(byLayer) max_layer = 7;

  TFile *tf = new TFile(input_name);
  TTree *tmpTr = (TTree*)tf->Get(tree_name);
  EqualizeChargesInPtBins(cutEn, cuts, "output_file");
  //Create tmp files to stop memory errors, basic cuts on full Tree
  TFile* tmpTF = new TFile("tmp1.root","recreate");
  std::cout << "Copying Tree" << std::endl;
  TTree *cutEn = tmpTr->CopyTree(Form(cuts));
  std::cout << "Copied" << std::endl;
  std::cout << "Closing input" << std::endl;
  tf->Close();

  //If testing nEntries performance, will slim here
  for (int nCut = 1; nCut < nCuts; nCut++) {
    std::cout << "Starting cut number " << nCut << std::endl;
    std::cout << "Current number of entries = " << cutEn->GetEntries() << std::endl;
    if(nCut == 1){
      std::cout << "nCut 1, using full file" << std::endl;
    }
    else{
      std::cout << "Slimming tmp file" << std::endl;
      int total_entries = cutEn->GetEntries();
      int nCut_entries = int(total_entries/2.0);
      std::cout << "nCut = " << nCut << std::endl;
      std::cout << "total_entries = " << total_entries << std::endl;
      std::cout << "nCut_entries = " << nCut_entries << std::endl;
      std::cout << "Taking first half" << std::endl;
      cutEn = cutEn->CloneTree(nCut_entries);
    }

    //Finished base cloning, create output CSV "GE11_out_CutNumber*.csv"
    std::cout << "New number of entries = " << cutEn->GetEntries() << std::endl;
    std::ofstream myfile;
    std::cout << "Creating CSV file " << Form("%s.csv", outname_prefix) << std::endl;
    myfile.open (Form("%s.csv", outname_prefix));
    //std::cout << "Creating CSV file " << Form("%s_out_CutNumber%d.csv", outname_prefix, nCut) << std::endl;
    //myfile.open (Form("%s_out_CutNumber%d.csv", outname_prefix, nCut));
    double dx, dy, dz, dphix, dphiy, dphiz;
    int detNum;
    dz = 0.0; dphix = 0.0; dphiy = 0.0;



    //Loop over every region/chamber/layer*
    std::cout << "Starting Chamber loop" << std::endl;
    for (int j = -1; j < 2; j = j + 2){             // Region loop
      for (int i = 0; i<36;i++){                    // Chamber loop
        for (int k = 1; k < max_layer; k++){        // Layer loop
          detNum = j*(i+101);
          std::cout << "at chamber " << detNum << " and layer " << k << std::endl;

          //Copy Tree only using hits on the current Detector
          TFile* tmpTF = new TFile("tmp2.root","recreate");
          std::cout <<"About to copy tree" << std::endl;
          TTree* tt_tmp;
          
          if(byLayer){
            tt_tmp = cutEn->CopyTree(Form("rechit_detId==%d && prop_location[4] == %d",detNum, k));		//Only fits 1 chamber at a time (det_id)
          }
          else{
            tt_tmp = cutEn->CopyTree(Form("rechit_detId==%d", detNum));
          }
          std::cout << "Entries are on chamber are " << tt_tmp->GetEntries() << std::endl;

          if (tt_tmp->GetEntries()<28){continue;}
          muon_charge = tt_tmp->GetVector("muon_charge");
          //New hist of RdPhi to get STD and MEAN
          TH1F *h1 = new TH1F("h1", "h1 title", 100, -20, 20);
          tt_tmp->Project("h1", "RdPhi", "");

          //Fit RdPhi to get STD and MEAN
          TF1 f1 = TF1("f1", "gaus", -2, 2);
          f1.SetParLimits(1, -2, 2);
          f1.SetParLimits(2, 0, 2);
          h1->Fit("f1", "R");
          float fitMean = f1.GetParameter(1);
          float fitStd = f1.GetParameter(2);

          //Copy only RdPhi within 1.6sigma of mean
          std::cout << "Starting small copy" << std::endl;
          tt = tt_tmp->CopyTree(Form("RdPhi <= (%f + (1.6*%f)) && RdPhi >= (%f - (1.6*%f))", fitMean, fitStd, fitMean, fitStd));

	  
          //If there are no events on the chamber it is skipped
          if (tt->GetEntries() == 0){
            myfile << detNum << ", " << 0 << ", " << 0 << ", " << 0 << ", " << 0 << ", " << 0 << ", " << 0 << ", " << 0 << "\n";
            delete tmpTF;
            delete tt_tmp;
            continue;
          }

          //Get variables for running Fit
          tt->SetBranchAddress("RdPhi", &mResidual);
          tt->SetBranchAddress("prop_LP", &mLocalPoint);
          tt->SetBranchAddress("prop_GP", &mGlobalPoint);


          mEvents = tt->GetEntries();
          doFit(doDx, doDy, doDphiz);
          dx = mResult[0];
          dy = mResult[1];
          dphiz = mResult[2];

          //Save alignment solutions to csv
          if(byLayer){
            myfile << detNum << k << ", " << dx << ", " << dy << ", " << dz << ", " << dphix << ", " << dphiy << ", " << dphiz << ", " << mEvents << "\n";
          }
          else{
            myfile << detNum << ", " << dx << ", " << dy << ", " << dz << ", " << dphix << ", " << dphiy << ", " << dphiz << ", " << mEvents << "\n";
          }
          delete tt_tmp;
          delete tmpTF;
        }
      }
    }
    myfile.close();
  }
  //tf->Close();
}
