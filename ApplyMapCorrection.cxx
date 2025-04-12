// Authors: Stefano Piacentini, Giorgio Dho
// Date: 2025-04-12
// Description: This program applies a pixel-by-pixel map correction to tracks of a tree in a ROOT file.
// Usage: ./ApplyMapCorrection.exe <input file> <path to map> <output file>
// Compile with:
//     g++ cygno-analyzer/Analyzer.cxx ApplyMapCorrection.cxx -o run.exe `root-config --libs --cflags` -lSpectrum

#include <iostream>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <string>
#include <numeric>
#include <stdexcept>
#include <algorithm>
#include <filesystem>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "cygno-analyzer/Analyzer.h" // Not needed, but included for completeness

using namespace std;

void ScIndicesElem(int nSc, UInt_t npix, float* sc_redpixID, int &nSc_red, vector<int>& B, vector<int>& E);

int main(int argc, char** argv) {

    bool verbose = false;

    // Check if the correct number of arguments is provided
    if(argc < 4) {
        cerr<<"No file path!\nSuggested use: ./progname.exe <input file> <path to map> <output file>"<<endl; 
        exit(EXIT_FAILURE);
    } else if(argc > 4) {
        cerr<<"Too many arguments!\nSuggested use: ./progname.exe <input file> <path to map> <output file>"<<endl; 
        exit(EXIT_FAILURE);
    }

    // Get the file name from the command line argument
    string namein  = argv[1];
    string namemap = argv[2];
    string nameout = argv[3];

    // Check if the input and output file names are valid
    if (nameout == namein) {
        cerr<<"Output file name is the same as input file name: this is not safe! Exiting..."<<endl;
        exit(EXIT_FAILURE);
    }

    // Check if the input file exists
    if (!std::filesystem::exists(namein)) {
        cerr<<"Input file does not exist: "<<namein<<endl;
        exit(EXIT_FAILURE);
    }

    // Check if the map file exists
    if (!std::filesystem::exists(namemap)) {
        cerr<<"Map file does not exist: "<<namemap<<endl;
        exit(EXIT_FAILURE);
    }

    filesystem::copy(namein, nameout, filesystem::copy_options::overwrite_existing);
    // Check if the output file was created successfully
    if (!std::filesystem::exists(nameout)) {
        cerr<<"Output file was not created: "<<nameout<<endl;
        exit(EXIT_FAILURE);
    }

    // Open the copied output file fo the addition of the new branch
    TFile* f = TFile::Open(Form("%s",nameout.c_str()), "update");
    TTree* tree = (TTree*)f->Get("Events");

    // Open the map file
    TFile* f_map = TFile::Open(Form("%s",namemap.c_str()), "read");
    TH2F* hmap = (TH2F*)f_map->Get("hman");

    int Nbinsx = hmap->GetNbinsX();
    int Nbinsy = hmap->GetNbinsY();

    // DEBUG
    //cout<<"Map size: "<<Nbinsx<<" "<<Nbinsy<<endl;
    //cout<<"Map min: "<<hmap->GetXaxis()->GetXmin()<<" "<<hmap->GetYaxis()->GetXmin()<<endl;
    //cout<<"Map max: "<<hmap->GetXaxis()->GetXmax()<<" "<<hmap->GetYaxis()->GetXmax()<<endl;
    //cout<<"Map bin size: "<<hmap->GetXaxis()->GetBinWidth(1)<<" "<<hmap->GetYaxis()->GetBinWidth(1)<<endl;

    int binwidth  = hmap->GetXaxis()->GetBinWidth(1);
    int binheight = hmap->GetYaxis()->GetBinWidth(1);
    


    //define variables required
    unsigned int nSc     = 0;
    int          nSc_red = 0;
    UInt_t       Nredpix = 0;
    int run;

    int nbinsx = 2304;
    int nbinsy = 2304;
    int nmax   = nbinsx * nbinsy;
    int nscmax = 400000; //100 tracks per image as maximum 

    vector<float> sc_redpixID(nmax, -2.);
    vector<int>   XPix(nmax, -2);
    vector<int>   YPix(nmax, -2);
    vector<float> ZPix(nmax, -2.);
    vector<float> integral(nscmax, -2.);

    //Link the variables to the tree branches
    tree->SetBranchAddress("run",&run);   
    tree->SetBranchAddress("nSc",&nSc);
    tree->SetBranchAddress("sc_redpixIdx",sc_redpixID.data());
    tree->SetBranchAddress("nRedpix",&Nredpix);
    tree->SetBranchAddress("redpix_ix",XPix.data());
    tree->SetBranchAddress("redpix_iy",YPix.data());
    tree->SetBranchAddress("redpix_iz",ZPix.data());
    tree->SetBranchAddress("sc_integral",integral.data());

    //////////Analysis Variables//////////////
    vector<int> BeginScPix;
    vector<int> EndScPix;

    vector<float>* integral_mapcorr;
    TBranch *sc_integral_mapcorr = tree->Branch("sc_integral_mapcorr",&integral_mapcorr);

    
    int counter = 0;
    //Cycle on events
    cout<<"This run has "<<tree->GetEntries()<<" entries"<<endl;
    for(int k=0;k<tree->GetEntries();k++) {

        tree->GetEntry(k);
        if(integral_mapcorr->size()!=0) integral_mapcorr->clear();

        //for reduced pixels:
        ScIndicesElem(nSc, Nredpix ,sc_redpixID.data(), nSc_red, BeginScPix,EndScPix);
        
        for(int clu=0;clu<nSc;clu++) {
            
            if(sc_redpixID[clu]!=-1) {

                vector<float> ZPixCorr;

                for(int i=BeginScPix[clu];i<EndScPix[clu];i++) {
                    float corr = 1.0;
                    if(run<59253) { // fixed orientation from Run5 on
                        corr = hmap->GetBinContent((int)((nbinsy-YPix[i])/binheight), (int)((nbinsx-XPix[i])/binwidth));
                    } else {
                        corr = hmap->GetBinContent((int)((nbinsx-XPix[i])/binwidth), (int)((nbinsy-YPix[i])/binheight));
                    }
                    ZPixCorr.push_back(ZPix[i]/corr);
                }

                float corrint = accumulate(ZPixCorr.begin(), ZPixCorr.end(), 0.0);

                //DEBUG
                //float newint  = accumulate(ZPix.begin()+BeginScPix[clu], ZPix.begin()+EndScPix[clu], 0.0);
                //float corr = corrint/newint;

                integral_mapcorr->push_back(corrint);

                counter++;
            } else {
                integral_mapcorr->push_back(integral[clu]);
            }

        }

        sc_integral_mapcorr->Fill();

    }

    // Save the modified tree to the output file
    f->cd();
    tree->Write("", TObject::kOverwrite);
    f->Close();
    f_map->Close();

    return 0;

}

// Functions
void ScIndicesElem(int nSc, UInt_t npix, float* sc_redpixID, int &nSc_red, vector<int>& B, vector<int>& E) {
  B.clear();
  E.clear();

  vector<float> sc_redpix_start;

  int parcount=0;

  for(int i=0; i<nSc; i++){
    if(sc_redpixID[i]>=0)  sc_redpix_start.push_back(sc_redpixID[i]);
  }

  nSc_red = sc_redpix_start.size();

  sc_redpix_start.push_back(npix);

  for(int i=0;i<sc_redpix_start.size()-1;i++){
    B.push_back(sc_redpix_start[i]);
    E.push_back(sc_redpix_start[i+1]);
    //std::cout<<B[i]<<" "<<E[i]<<endl; // DEBUG
  }

  sc_redpix_start.clear();

  return;

}