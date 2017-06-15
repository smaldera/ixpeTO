#include <iostream>
#include <TGraphErrors.h>                 // ci serve per istanziare grafici
#include <TMultiGraph.h>
#include <TAxis.h>                        // ci serve per manipolare gli assi dei grafici
#include <TLegend.h>
#include <TCanvas.h>                      // ci serve per disegnare i grafici
#include <TF1.h>                          // ci serve per scrivere le funzioni con cui fittare i grafici
#include "TVirtualFitter.h"               // ci serve come supporto per i fit
#include <TMath.h>
#include <string>
#include <TColor.h>
#include <TFormula.h>

using namespace std;

int mod_factor_vs_energy(){
  
  //Set di dati SIM
  const int ndata = 15;
  double energy[ndata]    = {2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0}; //[keV]
  double modulationFactor[ndata]   = {0.0850774926623, 0.167350250183, 0.244692798033, 0.306767319749, 0.348545691649, 0.383639843787, 0.407141488272, 0.435440247693, 0.461647508867, 0.490078282126, 0.518083903771, 0.54971168885, 0.558821519856, 0.57587527119, 0.590974684076};
  double thr1[ndata]  = {13, 8, 8, 8, 8, 5, 7, 5, 5, 5, 5, 5, 6, 5, 6};
  double thr2[ndata]  = { 8, 6, 5, 6, 6, 8, 5, 6, 5, 5, 8, 6, 5, 7, 9};
  double sThr1[ndata] = {0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.0001};
  double sThr2[ndata] = {0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001};
  //double sThr1[ndata] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  //double sThr2[ndata] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  
  double sEnergy[ndata]            = {0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001};
  
  double sModulationFactor[ndata]  = {0.00170015150827, 0.00183604733826, 0.00213269155356, 0.0025051795073, 0.00294589142676, 0.00341292162013, 0.00392178818043, 0.00446024618596, 0.00504053474788, 0.00559475017463, 0.00616594364876, 0.00680620474263, 0.00736353817818, 0.00466800320014, 0.00871803072181};



  //double sModulationFactor[ndata]  = {0.00156784835237, 0.00158063463164, 0.00174058556032, 0.0019778551829, 0.00228659626583, 0.00262167092106, 0.00299710083462, 0.00339441986975, 0.00382744839992, 0.0042501403814, 0.00469398539612, 0.00520288994111, 0.00564490618984, 0.00673535684944};

 
  //Set di dati VERI
  //xpol_ 2735, 2747, 2769
  const int ndataV = 3;
  double energyV[ndataV]            = {4.5, 6.4, 3.7}; //[keV]
  double modulationFactorV[ndataV]  = {0.384469553564, 0.46221075426, 0.331093023025};
  double thr1V[ndataV]  = {10, 6, 9};
  double thr2V[ndataV]  = { 7, 9, 20};
  double sThr1V[ndataV] = {0.00001, 0.00001, 0.00001};
  double sThr2V[ndataV] = {0.00001, 0.00001, 0.00001};
  double sEnergyV[ndataV]           = {0.00001, 0.00001, 0.00001};
  double sModulationFactorV[ndataV] = {0.00473813359058, 0.00423196298627, 0.00246242356854};
  //double sModulationFactorV[ndataV] = {0.00363470146419, 0.00320874764941, 0.00192225993575}; 

  
  
  //Tela
  TCanvas *c = new TCanvas("c","Modulation factor vs Energy", 200, 10, 600, 400);
  c->cd();
  
  //Grafico
  TGraphErrors *graph = new TGraphErrors(ndata, energy, modulationFactor, sEnergy, sModulationFactor);
  graph -> SetMarkerStyle(20);
  graph -> SetMarkerColor(1);
  graph -> SetMarkerSize(0.7);

  TGraphErrors *graphV = new TGraphErrors(ndataV, energyV, modulationFactorV, sEnergyV, sModulationFactorV);
  graphV -> SetMarkerStyle(22);
  graphV -> SetMarkerColor(kCyan+1);//(kOrange+1);
  graphV -> SetMarkerSize(1.3);  

  //-----------------------------------------------------------------------//  

  
  //Fit
  double min = 2.0;
  double max = 8.0;

  TF1 *fit4 = new TF1("fit4","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", min, max); 
  fit4 -> SetLineColor(kRed);
  graph->Fit(fit4,"MR");
  cout << "Chi^2: " << fit4->GetChisquare() << ", number of DoF: " << fit4->GetNDF() << " (Probability: " << fit4->GetProb() << ")." << endl;

  TF1 *fit5 = new TF1("fit","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*x*x*x*x*x", min, max); 
  fit5 -> SetLineColor(kSpring-6);
  graph->Fit(fit5,"MR+");
  cout << "Chi^2: " << fit5->GetChisquare() << ", number of DoF: " << fit5->GetNDF() << " (Probability: " << fit5->GetProb() << ")." << endl;

  double s0 = fit5 -> GetParError(0);
  double s1 = fit5 -> GetParError(1);
  double s2 = fit5 -> GetParError(2);
  double s3 = fit5 -> GetParError(3);
  double s4 = fit5 -> GetParError(4);
  double s5 = fit5 -> GetParError(5);

  //TFormula *smu = new TFormula("smu","sqrt(s0*s0+s1*s1*x*x+s2*s2*x*x*x*x+s3*s3*x*x*x*x*x*x+s4*s4*x*x*x*x*x*x*x*x+s5*s5*x*x*x*x*x*x*x*x*x*x)");
  
  TFormula *smu = new TFormula("smu","sqrt([0]*[0]+[1]*[1]*x*x+[2]*[2]*x*x*x*x+[3]*[3]*x*x*x*x*x*x+[4]*[4]*x*x*x*x*x*x*x*x+[5]*[5]*x*x*x*x*x*x*x*x*x*x)");
  smu -> SetParameter(0,s0);
  smu -> SetParameter(1,s1);
  smu -> SetParameter(2,s2);
  smu -> SetParameter(3,s3);
  smu -> SetParameter(4,s4);
  smu -> SetParameter(5,s5);
    
  cout << "2769 -> sim: mu = fit5(3.7 keV) = " << fit5->Eval(3.7) << " +- " << smu->Eval(3.7) << endl;
  cout << "2735 -> sim: mu = fit5(4.5 keV) = " << fit5->Eval(4.5) << " +- " << smu->Eval(4.5) << endl;
  cout << "2747 -> sim: mu = fit5(6.4 keV) = " << fit5->Eval(6.4) << " +- " << smu->Eval(6.4) << endl;
  
  TLegend *legend = new TLegend(.1,.75,.25,.9,"Legend");
  legend -> AddEntry(graph, "sim data", "p");
  legend -> AddEntry(graphV, "real data", "p");
  legend -> AddEntry(fit4, "poly4", "l");
  legend -> AddEntry(fit5, "poly5", "l");
  
  TMultiGraph *mg = new TMultiGraph();
  mg -> Add(graph);
  mg -> Add(graphV);
  mg -> SetTitle("Modulation factor vs Energy (sim + real data)");
  mg -> Draw("AP");
  mg -> GetXaxis() -> SetTitle("Energy [keV]");
  mg -> GetYaxis() -> SetTitle("Modulation factor");
  gPad -> Modified();
  
  legend -> Draw("same");


  //1st threshold vs energy
  TCanvas *c1 = new TCanvas("c1","First threshold vs Energy", 200, 10, 600, 400);
  c1->cd();
  
  TGraphErrors *graph1 = new TGraphErrors(ndata, energy, thr1, sEnergy, sThr1);
  graph1-> SetMarkerStyle(34);
  graph1-> SetMarkerColor(1);
  graph1-> SetMarkerSize(1.2);

  TGraphErrors *graph1V = new TGraphErrors(ndataV, energyV, thr1V, sEnergyV, sThr1V);
  graph1V -> SetMarkerStyle(22);
  graph1V -> SetMarkerColor(kCyan+1);
  graph1V -> SetMarkerSize(1.2);

  TLegend *legend1 = new TLegend(.1,.75,.25,.9,"Legend");
  legend1 -> AddEntry(graph1, "sim data", "p");
  legend1 -> AddEntry(graph1V, "real data", "p");
  
  TMultiGraph *mg1 = new TMultiGraph();
  mg1 -> Add(graph1);
  mg1 -> Add(graph1V);
  mg1 -> SetTitle("1st pass moments analysis vs Energy");
  mg1 -> Draw("AP");
  mg1 -> SetMinimum(4.);
  mg1 -> SetMaximum(18.);
  mg1 -> GetXaxis() -> SetTitle("Energy [keV]");
  mg1 -> GetYaxis() -> SetTitle("1st threshold");
  gPad -> Modified();

  legend1 -> Draw("same");



  //2nd threshold vs energy
  TCanvas *c2 = new TCanvas("c2","Second threshold vs Energy", 200, 10, 600, 400);
  c2->cd();
  
  TGraphErrors *graph2 = new TGraphErrors(ndata, energy, thr2, sEnergy, sThr2);
  graph2-> SetMarkerStyle(34);
  graph2-> SetMarkerColor(1);
  graph2-> SetMarkerSize(1.2);

  TGraphErrors *graph2V = new TGraphErrors(ndataV, energyV, thr2V, sEnergyV, sThr2V);
  graph2V -> SetMarkerStyle(22);
  graph2V -> SetMarkerColor(kCyan+1);
  graph2V -> SetMarkerSize(1.2);

  TLegend *legend2 = new TLegend(.1,.75,.25,.9,"Legend");
  legend2 -> AddEntry(graph2, "sim data", "p");
  legend2 -> AddEntry(graph2V, "real data", "p");
  
  TMultiGraph *mg2 = new TMultiGraph();
  mg2 -> Add(graph2);
  mg2 -> Add(graph2V);
  mg2 -> SetTitle("2nd pass moments analysis vs Energy");
  mg2 -> Draw("AP");
  mg2 -> SetMinimum(4.);
  mg2 -> SetMaximum(21.);
  mg2 -> GetXaxis() -> SetTitle("Energy [keV]");
  mg2 -> GetYaxis() -> SetTitle("2nd threshold");
  gPad -> Modified();
  
  legend2 -> Draw("same");


  
  /*
  if(fit->GetChisquare() < chiquadro95[fit->GetNDF()])
    cout << "Chi quadro accettabile" << endl;
  else
    cout << "Chi quadro non accettabile" << endl;
  */
  
    /*
  // intersezione tra rette "fit" e "fit1"
  double xInt  = (fit1->GetParameter(1) - fit->GetParameter(1)) / (fit->GetParameter(0) - fit1->GetParameter(0));
  double sxInt = x * ((fit1->GetParError(1) + fit->GetParError(1)) / (fit1->GetParameter(1) - fit->GetParameter(1)) + (fit1->GetParError(0) + fit->GetParError(0))/(fit->GetParameter(0) - fit1->GetParameter(0)));

  cout << "Intersezione tra le due rette: " << "\n\txInt = (" << xInt << "+-" << sxInt << ") u.d.m." << endl;
  
  */


  return 0;
  
}
