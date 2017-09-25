// compile with:
//c++ -Wl,--no-as-needed   readSimTree_Simo.cc -o readSimTree_Simo.exe `root-config --libs --cflags`

//reads the root ttree resulting from the ixpe simulation.
// for each events draws 3 views of the energy deposit and a xy view of the different tracks 
// a canvas is saved for each event.

// usage:
// ./readSimTree_Simo.exe inputfile.root  outfile.root  n_max_events (optional)
// example:
//./readSimTree_Simo.exe  sim_3KeV_-175deg.root prova.root 30



#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <math.h>


#include "TApplication.h"
#include <TCanvas.h>
#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TLegend.h>
#include <TH2F.h>
#include <TMarker.h>
#include <TF1.h>
#include <TRandom2.h>


using namespace std;


int trovaColore( TString process ){

  int red=2;
  int blu =4;
  int fucsia=6;
  int green=3;

  if (process=="phot") return blu;
  else if (process=="eIoni") return red;
  else if (process=="phot_auger") return green;
  else if (process=="Undefined") return fucsia;
 
  return 0;
}








int main(int argc, char **argv){

   if (argc<3 ){
      cout<<" usage:  ./readSimTree_Simo.exe inputfile.root  outfile.root  n_max_events (optional, default =100, <=0 --> all events)  "<<endl;
      return 0;
  }

   cout<<"max (2,4)="<<max (2,4)<<endl;
   cout<<"min (2,4)="<<min(2,4)<<endl;

   TRandom2 *RndGaussiana = new TRandom2;
   RndGaussiana->SetSeed( int (time (NULL)*4)  );
 //gauss= RndGaussiana->Gaus(s,sigma);
   TFile *fout = new  TFile(argv[2],"recreate");
   TFile *f = new TFile(argv[1],"open");

   int maxEvents=100;
   if (argc==4 )   maxEvents=atoi(argv[3]);
   cout<<endl<<"======================================================="<<endl<<endl;
   cout<<" running on : "<<argv[1]<<endl;
   cout<<" outfile = "<<argv[2]<<endl;
   cout<<" max events = "<<maxEvents<<endl<<endl;
   cout<<"======================================================="<<endl<<endl;
  
  
  TTree *tree= (TTree *)f->Get("ixpe");
  
  TTree *tree2=(TTree *)f->Get("ixpe2");

  double sigmaBlur=0.04;
  
  Int_t    NumEvents, TrackID,StepNum,FirstStep,LastStep;
  Double_t PosX, PosY, PosZ,EDep, PosXnorm, PosYnorm, PosZnorm ; 
  //Char_t StepProcName[100];
  Char_t TrackProcName[100];
  int nEntries=tree->GetEntries();
  int nEntries2=tree2->GetEntries();
  int n_Evts=0;
  cout<<nEntries2<<" events found in TTree ixpe2"<<endl;
  cout<<nEntries<<" entries in TTree ixpe (gas cell hits)"<<endl; 

  
  tree->SetBranchAddress("NumEvents",       &NumEvents     );     
  tree->SetBranchAddress("PosX",       &PosX     );     
  tree->SetBranchAddress("PosY",       &PosY     );     
  tree->SetBranchAddress("PosZ",       &PosZ     );     
  tree->SetBranchAddress("EDep",       &EDep     );     
  tree->SetBranchAddress("TrackID",    &TrackID     );     
  //tree->SetBranchAddress("StepProcName", &StepProcName);     
  tree->SetBranchAddress("TrackProcName", &TrackProcName);     
  tree->SetBranchAddress("StepNum", &StepNum);     
  tree->SetBranchAddress("LastStep", &LastStep);     
  tree->SetBranchAddress("FirstStep", &FirstStep);     
  

  
  
  Double_t x_conv, y_conv, z_conv, PePhi; 
  tree2->SetBranchAddress("AbsPosX",       &x_conv    );     
  tree2->SetBranchAddress("AbsPosY",       &y_conv    );     
  tree2->SetBranchAddress("AbsPosZ",       &z_conv    );     
  tree2->SetBranchAddress("PePhi",        &PePhi    );     

    
  TH2F *h_xy=new TH2F("h_xy","x-y",1000,-1,1.5,1000,-1,1);
  TH2F *h_xz=new TH2F("h_xz","x-z",1000,-1,1,1000,-1,1 );
  TH2F *h_yz=new TH2F("h_yz","y-z",1000,-1,1,1000,-1,1 );

  TH1F *h_long=new TH1F("h_long","long",1000,-1,1);
  TH1F *h_longAuger=new TH1F("h_longAuger","longAuger",1000,-1,1);
  TH1F *h_longAll=new TH1F("h_longAll","long",1000,-1,1);

  TH1F *h_long_blur=new TH1F("h_long_blur","long",1000,-1,1);
  TH1F *h_longAuger_blur=new TH1F("h_longAuger_blur","longAuger",1000,-1,1);
  TH1F *h_longAll_blur=new TH1F("h_longAll_blur","long",1000,-1,1);

 

  
  TH1F *h_longLin=new TH1F("h_longLin","longLin",1000,-1,1);
 
  h_xy->GetXaxis()->SetTitle("X-x_{conv}"); h_xy->GetYaxis()->SetTitle("Y-y_{conv}");
  h_xz->GetXaxis()->SetTitle("X-x_{conv}"); h_xz->GetYaxis()->SetTitle("Z-z_{conv}");
  h_yz->GetXaxis()->SetTitle("Y-y_{conv}"); h_yz->GetYaxis()->SetTitle("Z-z_{conv}");
 


  TCanvas *c0=new TCanvas("c0","",0);
  c0->Divide(2,2);
    
  Int_t prevEventId=-10;

  double maxX=-1000;
  double minX=1000;
  double maxY=-1000;
  double minY=1000;
  double maxZ=-1000;
  double minZ=1000;

  vector<Double_t> x[15];
  vector<Double_t> y[15];
  vector<Double_t> z[15];
  
  //vector<Double_t> eDep;

  TString processName[15];
  double prev_x=0;
  double prev_y=0;
  
  tree2->GetEntry(0);

  
  int maxTrackID=0;
  double dist=0;
  double distAuger=0;
  int n_pheDep=0;
  
  for (int i=0; i< nEntries; i++){
      tree->GetEntry(i);
      
      //cout<<"i = "<<i<<" NumEvents = "<<NumEvents<<" x = "<<PosX<<" stepProcName = "<<TrackProcName<<endl;
      if (prevEventId==NumEvents || prevEventId==-10 ||  n_pheDep==0 ){ // non cabiare evento!!! 

	PosXnorm=PosX-x_conv;
	PosYnorm=PosY-y_conv;
	PosZnorm=PosZ-z_conv;

	if (prev_x!=0){
 	     if  ( processName[TrackID]=="phot_auger"){
	          distAuger= distAuger -  pow(   (pow( (PosX-prev_x),2)+ pow( (PosY-prev_y),2) ),0.5);
		  h_longAuger->Fill(distAuger,EDep);
		  h_longAll->Fill(distAuger,EDep);

		  double distAugerBlur=RndGaussiana->Gaus(distAuger,sigmaBlur);
		  h_longAuger_blur->Fill(distAugerBlur,EDep);
		  h_longAll_blur->Fill(distAugerBlur,EDep);


		  
	     }
	     if ( processName[TrackID]=="phot"){
	          dist= dist+  pow(   (pow( (PosX-prev_x),2)+ pow( (PosY-prev_y),2) ),0.5);
		  h_long->Fill(dist,EDep);
		  h_longAll->Fill(dist,EDep);
		  //cout<<" dist= "<<dist<<" Edep = "<<EDep<<endl;
		  
		  double distBlur=RndGaussiana->Gaus(dist,sigmaBlur);
		  h_long_blur->Fill(distBlur,EDep);
		  h_longAll_blur->Fill(distBlur,EDep);
		  n_pheDep++;
	     }

	} else{
	  dist=0;
	  distAuger=0;
	  
	}
	
	

	double distL=pow(   (pow( PosX,2)+ pow( PosY,2) ),0.5);
	h_longLin->Fill(distL,EDep);
	
	//cout<<" trackID= "<<TrackID<<"  processName[TrackID]="<<  processName[TrackID]<<" numStep = "<<StepNum<<" firstStep= "<<FirstStep<<" Last = "<<LastStep<<endl;
	


	
	h_xy->Fill((PosXnorm),(PosYnorm),EDep);
	h_xz->Fill(PosXnorm,PosZnorm,EDep);
	h_yz->Fill(PosYnorm,PosZnorm,EDep);

	if (TrackID>maxTrackID) maxTrackID=TrackID;
	
	if (TrackID<15){
	  x[TrackID].push_back(PosXnorm);
	  y[TrackID].push_back(PosYnorm);
	  z[TrackID].push_back(PosZnorm); 
	  processName[TrackID]=TrackProcName;

	}
	
	if (PosXnorm>maxX) {maxX=PosXnorm;}
	if (PosXnorm<minX) {minX=PosXnorm;}
	if (PosYnorm>maxY) {maxY=PosYnorm;}
	if (PosYnorm<minY) {minY=PosYnorm;}
	if (PosZnorm>maxZ) {maxZ=PosZnorm;}
	if (PosZnorm<minZ) {minZ=PosZnorm;}

	prev_x=PosX;
	prev_y=PosY;

	
      }  
      else{
	//TCanvas *c0=new TCanvas("c0","",0);
	if (n_Evts>maxEvents  && maxEvents>0 ) break;
	if (n_Evts %1000==0) cout<<"read "<<n_Evts<<" events... "<<endl;
	//cout<<" ===> trovato nuovo evento! salvo... precedente   i = "<<i<<" NumEvents = "<<NumEvents<<" prevEvent = "<<prevEventId<<" n_Evts= "<<n_Evts<<endl;
	

	//TMarker *conv_p=new TMarker(x_conv,y_conv,20);
	TMarker *conv_p=new TMarker(0.,0.,20);
	
	//conv_p->SetMarkerSyle(20);
	double  m=(tan(PePhi));
	TF1 *fdir=new TF1("fdir","[0]*x",-0.2,0.2);
	fdir->SetParameter(0,m);
	fdir->SetLineColor(1);
	fdir->SetLineWidth(1);
	fdir->SetLineStyle(2);

	c0->cd(1);
	
	
	h_xy->GetXaxis()->SetRangeUser(min(minX,minY)-0.02, max(maxX,maxY)+0.02);
	h_xy->GetYaxis()->SetRangeUser(min(minX,minY)-0.02, max(maxX,maxY)+0.02 );
	h_xy->Draw("box");
	conv_p->Draw("samep");
	fdir->Draw("samel");
	
	h_xz->GetXaxis()->SetRangeUser(min(minX,minZ)-0.02, max(maxX,maxZ)+0.02);
	h_xz->GetYaxis()->SetRangeUser(min(minX,minZ)-0.02, max(maxX,maxZ)+0.02);
	
	c0->cd(3);
	h_xz->Draw("box");
	conv_p->Draw("samep");
	h_yz->GetXaxis()->SetRangeUser(min(minY,minZ)-0.02, max(maxY,maxZ)+0.02);
	h_yz->GetYaxis()->SetRangeUser(min(minY,minZ)-0.02, max(maxY,maxZ)+0.02);
	
	c0->cd(2);
	h_yz->Draw("box");
	conv_p->Draw("samep");
	
      //==============================================
	//cout<<" maxTrackID = "<<maxTrackID<<endl;
	c0->cd(4);

	double x_scala[2]={minX-0.05,maxX+0.05};
	double y_scala[2]={minY-0.05,maxY+0.05};
	//	double z_scala[2]={-0.2.,0.2};
	
	// creo tracce  per diverse particelle....
	TGraph *g_scala=new TGraph(2,x_scala,y_scala);
	//TGraph2D *g_scala=new TGraph2D(2,x_scala,y_scala,z_scala);

	g_scala->GetXaxis()->SetTitle("X");
	g_scala->GetYaxis()->SetTitle("Y");

	g_scala->Draw("ap");
	TGraph *tracce_xy[15];
	//TGraph2D *tracce_xy[15];

	TLegend *leg = new TLegend(0.4,0.6,0.89,0.89);

	char nome[10];
	int gia_eIoni=0;
	for (int j=1; j<maxTrackID; j++){
	  sprintf(nome,"track_%d",j);
	  tracce_xy[j]=new TGraph(x[j].size(),&x[j][0],&y[j][0]);
	  //tracce_xy[j]=new TGraph2D(x[j].size(),&x[j][0],&y[j][0],&z[j][0]);
	  tracce_xy[j]->SetName(nome);
	  
	  tracce_xy[j]->SetMarkerStyle(7);
	  int colore=trovaColore(processName[j]);
	  tracce_xy[j]->SetMarkerColor(colore);  
	  tracce_xy[j]->SetLineColor(colore);
	  tracce_xy[j]->SetLineWidth(3);
	  tracce_xy[j]->Draw("lp");

	  if (colore==2){ 
	    if (gia_eIoni==0) leg->AddEntry( tracce_xy[j] ,processName[j],"l");
	    gia_eIoni=1;
	  } else if (colore!=6) leg->AddEntry( tracce_xy[j] ,processName[j],"l");
	  
	} //end for trackid 
	leg->Draw();
	conv_p->Draw();
	fdir->Draw("samel");
		
	sprintf(nome,"ev_%d",	n_Evts+1);
	//	cout<<"nome=  "<<nome<<endl;
	c0->SetName(nome);
	c0->Update();
	c0->Draw();

	fout->cd();
	//c0->Write(); // non salvo il canvas!!!!
	
	n_Evts++;
	//cout<<"       xconv old === "<<x_conv<<endl;
	if (n_Evts>nEntries2) break;
	tree2->GetEntry(n_Evts); // prende le grandezze dell'evento successivo!!!!!! brutto!!!!
	
		
	h_xy->Reset();
	h_xz->Reset();
	h_yz->Reset();
	//h_long->Reset();
	dist=0;
	distAuger=0;
	
	prev_x=0;
	prev_y=0;
	n_pheDep=0;
	
	for (int j=0; j<15; j++){
	
	  x[j].clear();
	  y[j].clear();
	  z[j].clear();
	  if (j<maxTrackID ){
	
	    // delete tracce_xy[j];
	
	  }
	}
	
	

	delete g_scala;
	
	
	maxX=-1000;
	minX=1000;
	maxY=-1000;
	minY=1000;
	maxZ=-1000;
	minZ=1000;
        //cout<<"       xconv new === "<<x_conv<<endl;
	//cout<<"       xconv new ttree1 === "<<PosX<<endl;

	//x_conv=PosX;
	//y_conv=PosY;
	//z_conv=PosZ;
	
	
	
	PosXnorm=PosX-x_conv;
	PosYnorm=PosY-y_conv;
	PosZnorm=PosZ-z_conv;
	//delete c0;
	h_xy->Fill(PosXnorm,PosYnorm,EDep);
	h_xz->Fill(PosXnorm,PosZnorm,EDep);
	h_yz->Fill(PosYnorm,PosZnorm,EDep);
	maxTrackID=0;



	if (TrackID>maxTrackID) maxTrackID=TrackID;
	
	if (TrackID<15){
	  x[TrackID].push_back(PosX);
	  y[TrackID].push_back(PosY);
	  processName[TrackID]=TrackProcName;
	}
	//cout<<" -------------- "<<endl<<endl;
      }	//end else
      
      
      prevEventId=NumEvents;
      
  }

  fout->cd();
  h_long->Write();
  h_longAuger->Write();
  h_longAll->Write();

  h_long_blur->Write();
  h_longAuger_blur->Write();
  h_longAll_blur->Write();
 
  h_longLin->Write();
  cout<<endl<<"closing input root file... ";
  f->Close();
  cout<<"  done!"<<endl;

  cout<<"closing out root file... ";
  
  
  fout->Close();
  cout<<"  done!"<<endl;
  cout<<endl<<" simulated ph.electon tracks saved in "<<argv[2]<<endl<<endl;

  
}
