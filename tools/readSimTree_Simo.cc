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



   TFile *fout = new  TFile(argv[2],"recreate");
   TFile *f = new TFile(argv[1]);

   int maxEvents=100;
   if (argc==4 )   maxEvents=atoi(argv[3]);
   cout<<endl<<"======================================================="<<endl<<endl;
   cout<<" running on : "<<argv[1]<<endl;
   cout<<" outfile = "<<argv[2]<<endl;
   cout<<" max events = "<<maxEvents<<endl<<endl;
   cout<<"======================================================="<<endl<<endl;
  
  
  TTree *tree= (TTree *)f->Get("ixpe");
  TTree *tree2=(TTree *)f->Get("ixpe2");

  Int_t    NumEvents, TrackID;
  Double_t PosX, PosY, PosZ,EDep ;  
  //Char_t StepProcName[100];
  Char_t TrackProcName[100];
  
  
  
  tree->SetBranchAddress("NumEvents",       &NumEvents     );     
  tree->SetBranchAddress("PosX",       &PosX     );     
  tree->SetBranchAddress("PosY",       &PosY     );     
  tree->SetBranchAddress("PosZ",       &PosZ     );     
  tree->SetBranchAddress("EDep",       &EDep     );     
  tree->SetBranchAddress("TrackID",    &TrackID     );     
  //tree->SetBranchAddress("StepProcName", &StepProcName);     
  tree->SetBranchAddress("TrackProcName", &TrackProcName);     
  
  
  
  Double_t x_conv, y_conv, z_conv; 
  tree2->SetBranchAddress("AbsPosX",       &x_conv    );     
  tree2->SetBranchAddress("AbsPosY",       &y_conv    );     
  tree2->SetBranchAddress("AbsPosZ",       &z_conv    );     
  

    
  TH2F *h_xy=new TH2F("h_xy","x-y",300,-0.3,0.3,300,-0.3,0.3);
  TH2F *h_xz=new TH2F("h_xz","x-z",300,-0.3,0.3,300,-0.3,0.3 );
  TH2F *h_yz=new TH2F("h_yz","y-z",300,-0.3,0.3,300,-0.3,0.3 );
  
  h_xy->GetXaxis()->SetTitle("X-x_{conv}"); h_xy->GetYaxis()->SetTitle("Y-y_{conv}");
  h_xz->GetXaxis()->SetTitle("X-x_{conv}"); h_xz->GetYaxis()->SetTitle("Z-z_{conv}");
  h_yz->GetXaxis()->SetTitle("Y-y_{conv}"); h_yz->GetYaxis()->SetTitle("Z-z_{conv}");

  
  int nEntries=tree->GetEntries();
  int nEntries2=tree2->GetEntries();
  int n_Evts=0;

  cout<<nEntries2<<" events found in TTree ixpe2"<<endl;
  cout<<nEntries<<" entries in TTree ixpe (gas cell hits)"<<endl; 
  
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
  
  tree2->GetEntry(0);

  int maxTrackID=0;
  
  for (int i=0; i< nEntries; i++){
      tree->GetEntry(i);
      
      // cout<<"i = "<<i<<" NumEvents = "<<NumEvents<<" x = "<<PosX<<" stepProcName = "<<TrackProcName<<endl;
      if (prevEventId==NumEvents || prevEventId==-10 ){

	PosX=PosX-x_conv;
	PosY=PosY-y_conv;
	PosZ=PosZ-z_conv;
	
	h_xy->Fill((PosX),(PosY),EDep);
	h_xz->Fill(PosX,PosZ,EDep);
	h_yz->Fill(PosY,PosZ,EDep);

	if (TrackID>maxTrackID) maxTrackID=TrackID;
	
	if (TrackID<15){
	  x[TrackID].push_back(PosX);
	  y[TrackID].push_back(PosY);
	  z[TrackID].push_back(PosZ); 
	  processName[TrackID]=TrackProcName;

	}
	
	if (PosX>maxX) {maxX=PosX;}
	if (PosX<minX) {minX=PosX;}
	if (PosY>maxY) {maxY=PosY;}
	if (PosY<minY) {minY=PosY;}
	if (PosZ>maxZ) {maxZ=PosZ;}
	if (PosZ<minZ) {minZ=PosZ;}

	
      }  
      else{
	//TCanvas *c0=new TCanvas("c0","",0);
	if (n_Evts>maxEvents  && maxEvents>0 ) break;
	if (n_Evts %1000==0) cout<<"read "<<n_Evts<<" events... "<<endl;
	
	tree2->GetEntry(n_Evts);

	//TMarker *conv_p=new TMarker(x_conv,y_conv,20);
	TMarker *conv_p=new TMarker(0.,0.,20);
	
	//conv_p->SetMarkerSyle(20);

	c0->cd(1);
	
	
	h_xy->GetXaxis()->SetRangeUser(minX-0.02,maxX+0.02);
	h_xy->GetYaxis()->SetRangeUser(minY-0.02,maxY+0.02);
	h_xy->Draw("box");
	conv_p->Draw("samep");

	h_xz->GetXaxis()->SetRangeUser(minX-0.02,maxX+0.02);
	h_xz->GetYaxis()->SetRangeUser(minZ-0.02,maxZ+0.02);
	
	c0->cd(3);
	h_xz->Draw("box");
	conv_p->Draw("samep");
	h_yz->GetXaxis()->SetRangeUser(minY-0.02,maxY+0.02);
	h_yz->GetYaxis()->SetRangeUser(minZ-0.02,maxZ+0.02);
	
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

	
	sprintf(nome,"ev_%d",prevEventId);
	//	cout<<"nome=  "<<nome<<endl;
	c0->SetName(nome);
	c0->Update();
	c0->Draw();

	fout->cd();
	c0->Write();
	  
	n_Evts++;

	if (n_Evts>nEntries2) break;
	tree2->GetEntry(n_Evts); // prende le grandezze dell'evento successivo!!!!!! brutto!!!!
	
		
	h_xy->Reset();
	h_xz->Reset();
	h_yz->Reset();

	
	
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
	 
	//delete c0;
	h_xy->Fill(PosX,PosY,EDep);
	h_xz->Fill(PosX,PosZ,EDep);
	h_yz->Fill(PosY,PosZ,EDep);
	maxTrackID=0;



	if (TrackID>maxTrackID) maxTrackID=TrackID;
	
	if (TrackID<15){
	  x[TrackID].push_back(PosX);
	  y[TrackID].push_back(PosY);
	  processName[TrackID]=TrackProcName;
	}
	
	
	
      }	//end eels
      
      
      prevEventId=NumEvents;
      
  }

  cout<<endl<<"closing input root file... ";
  f->Close();
  cout<<"  done!"<<endl;

  cout<<"closing out root file... ";
  fout->Close();
  cout<<"  done!"<<endl;
  cout<<endl<<" simulated ph.electon tracks saved in "<<argv[2]<<endl<<endl;

  
}
