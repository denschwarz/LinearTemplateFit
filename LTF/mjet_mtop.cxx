/*
   Linear template fit of top quark mass to data from boosted top decays
*/
// -------------------------------------------------------------------- //

#include <iostream>

#include "LTF/LTF.h"
#include "LTF/LTF_ROOTTools.h"

#include <TStyle.h>
#include <TH2D.h>
#include <TFile.h>
#include <TSystem.h>

void AddCovMatrix(LTF &ltf, TH2D* ch, string name);

int mjet_mtop() {

   using namespace std;

   gStyle->SetOptStat(0);
   gSystem->Load("libLTF.so");
   TH1D::AddDirectory(false);
   TH1::SetDefaultSumw2(false);

   // open the file with the templates and the data
   TFile* f = new TFile("data/UnfoldingCombination_pseudo173_5.root", "READ");
   // content:   
   /*
   data;
   cov_total;
   mc_mtop1695;
   mc_mtop1715;
   mc_mtop1725;
   mc_mtop1735;
   mc_mtop1755;
   cov_Stat;
   cov_Theo;
   cov_Model;
   cov_Exp;
   cov_stat_mtop1695;
   cov_stat_mtop1715;
   cov_stat_mtop1725;
   cov_stat_mtop1735;
   cov_stat_mtop1755; 
   */

   // get the data hist
   TH1D* data = (TH1D*)f->Get("data");

   // the values of mt used for the templates
   const vector<double> mtvals{169.5, 171.5, 172.5, 173.5, 175.5}; // template reference points

   // fill the binning
   int nbins = data->GetNbinsX();
   vector<double> bins;
   for (int i=0; i<nbins+1; ++i){
      bins.push_back(data->GetXaxis()->GetXbins()->At(i));
      //cout << "bin [" << i << "] = " << bins[i] << endl;
   }

   TH1D* mc_mtop1695 = (TH1D*)f->Get("mc_mtop1695");
   TH1D* mc_mtop1715 = (TH1D*)f->Get("mc_mtop1715"); // this is the MC used to obtain the pseudo-data 
   TH1D* mc_mtop1725 = (TH1D*)f->Get("mc_mtop1725");
   TH1D* mc_mtop1735 = (TH1D*)f->Get("mc_mtop1735"); 
   TH1D* mc_mtop1755 = (TH1D*)f->Get("mc_mtop1755");

   TH2D* cov_tot  = (TH2D*)f->Get("cov_total");
   TH2D* cov_stat = (TH2D*)f->Get("cov_Stat");
   TH2D* cov_exp  = (TH2D*)f->Get("cov_Exp");
   TH2D* cov_mod  = (TH2D*)f->Get("cov_Model");
   TH2D* cov_theo = (TH2D*)f->Get("cov_Theo");



   // ------------------------------------------------ //
   // ---  Do linear template fit
   // ------------------------------------------------ //
   // --- instantiate LTF object
   LTF ltf;
   ltf.SetGamma(vector<double>{1});
   ltf.UseNuisanceParameters(false);
   ltf.UseLogNormalUncertainties(false);

   // --- add templates 
   ltf.AddTemplate(mtvals[0], nbins, mc_mtop1695->GetArray()+1); 
   ltf.AddTemplate(mtvals[1], nbins, mc_mtop1715->GetArray()+1); // this is the MC used to obtain the pseudo-data 
   ltf.AddTemplate(mtvals[2], nbins, mc_mtop1725->GetArray()+1); 
   ltf.AddTemplate(mtvals[3], nbins, mc_mtop1735->GetArray()+1); 
   ltf.AddTemplate(mtvals[4], nbins, mc_mtop1755->GetArray()+1); 

   // stat uncertainty in templates: get diagonal elements from cov matrix
   // note: stored on the template is actually the diaganol element of the theory uncertainty
   const vector<TString> names{"cov_stat_mtop1695", "cov_stat_mtop1715", "cov_stat_mtop1725", "cov_stat_mtop1735", "cov_stat_mtop1755"};
   // sanity check  
   if (mtvals.size() != names.size()){
      cerr << "Not the same size of templates and template errors (mtvals and names)." << endl;
      exit(-1);
   }
   // fill the diagonal elements
   for (int i=0; i<mtvals.size(); ++i){
     vector<double> stat_unc_temp_sq;
     TH2D* cov_mcstat = (TH2D*) f->Get(names[i]); 
     for (int j=1; j<cov_mcstat->GetNbinsX()+1; ++j){
       stat_unc_temp_sq.push_back(cov_mcstat->GetBinContent(j,j)); 
       //cout << "stat unc sq of bin " << j-1 << " = " << stat_unc_temp_sq[j-1] << endl;
     }
     ltf.AddTemplateErrorSquared("statY", mtvals[i], nbins, stat_unc_temp_sq.data(), 0.); // set template error dY   
   }


   // --- initialize data
   ltf.SetData( nbins, data->GetArray()+1);
   //ltf.AddUncorrelatedErrorSquared("stat.", nbins, data->GetSumw2()->GetArray()+1);

   // hand over the covariance matrix, 
   //AddCovMatrix(ltf, cov_tot, "total");
   AddCovMatrix(ltf, cov_stat, "stat");
   AddCovMatrix(ltf, cov_exp,  "exp");
   AddCovMatrix(ltf, cov_mod,  "mod");
   AddCovMatrix(ltf, cov_theo, "theo");

   // remove one bin from the fit, because it is a normalised measurement
   // todo: adjust binning
   ltf.SetFitRange(1, 4);

   // for ( const auto& s : shiftsnuisance ) ltf.AddError("",N,s->GetArray()+1,1.);
   // for ( const auto& s : shifts ) ltf.AddError("",N,s->GetArray()+1,1.);
   // for ( const auto& s : shifts ) ltf.AddError("",N,s->GetArray()+1,0.5,LTF::Uncertainty::External);

   LTF::LiTeFit fit = ltf.DoLiTeFit();
   fit.PrintFull();
   //fit.DoIterativeFitNewton(6,0.6,2,1);
   //fit.DoIterativeFitTaylor();
   //fit.PrintFull();

   // for plotting: need to adjust binning (remove first bin, same as done in SetFitRange)
   bins.erase(bins.begin(), bins.begin()+1);

   // plot the results
   LTF_ROOTTools::plotLiTeFit(fit,bins, "norm. cross section", "m_{t} [GeV]", "Jet mass [GeV]");
   
   return 0;

}

void AddCovMatrix(LTF &ltf, TH2D* ch, string name)
{
   using namespace std;
   int nbinsx = ch->GetNbinsX();
   int nbinsy = ch->GetNbinsY();

   vector< vector<double> > cov;
   vector<double> cols;

   for (int i=1; i<nbinsx+1; ++i){
      for (int j=1; j<nbinsy+1; ++j){
         cols.push_back(ch->GetBinContent(i, j));
         //cout << "cov[" << i << ", " << j << "] = " << ch->GetBinContent(i, j) << endl;
      }
      cov.push_back(cols);
      cols.clear();
   }

   ltf.AddError(name, cov);

}


//! ------------------------------------------------------------------------ //
//! main function
int main() {

   return mjet_mtop();

}
