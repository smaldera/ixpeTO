
{
  TFile *f1=new TFile("modFact_zero5_5_5_noCuts_phi2.root","open");
  TGraph *gDiff_phi2= f1->Get("mu_diff");
  gDiff_phi2->SetTitle("Rel. difference #mu_{max} - #mu_{std}, no cuts, (zero supp=5)");
  gDiff_phi2->GetYaxis()->SetTitle("(#mu_{max}-#mu_{std})/mu_{std}");
  
 TFile *f2=new TFile("modFact_zero5_5_5_noCuts_phi1.root","open");
 TGraph *gDiff_phi1= f2->Get("mu_diff");

 gDiff_phi1->SetMarkerStyle(4);
 

 TCanvas *c1 = new TCanvas("c1","",0);
 gDiff_phi2->Draw("ap");
 gDiff_phi1->Draw("p");

 TLegend *leg = new TLegend(0.5,0.7,0.89,0.89);
 leg->AddEntry( gDiff_phi2 ,"phi2  ","lp");
  leg->AddEntry( gDiff_phi1 ,"phi1  ","lp");
  leg->Draw();

  ////////////////////////////////////////
  // plot con qCuts


  
 TFile *f1qc=new TFile("modFact_zero5_5_5_qcuts_phi2.root","open");
 TGraph *gDiff_phi2_qc= f1qc->Get("mu_diff");
  gDiff_phi2_qc->SetTitle("Rel. difference #mu_{max} - #mu_{std}, quality cuts, (zero supp=5)");
 gDiff_phi2_qc->GetYaxis()->SetTitle("(#mu_{max}-#mu_{std})/mu_{std}");
  
 TFile *f2qc=new TFile("modFact_zero5_5_5_qcuts_phi1.root","open");
 TGraph *gDiff_phi1_qc= f2qc->Get("mu_diff");

 gDiff_phi1_qc->SetMarkerStyle(4);
 

 TCanvas *c2 = new TCanvas("c2","",0);
 gDiff_phi2_qc->Draw("ap");
 gDiff_phi1_qc->Draw("p");

 TLegend *leg_qc = new TLegend(0.5,0.7,0.89,0.89);
 leg_qc->AddEntry( gDiff_phi2_qc ,"phi2  ","lp");
  leg_qc->AddEntry( gDiff_phi1_qc ,"phi1  ","lp");
  leg_qc->Draw();

  


  


}
