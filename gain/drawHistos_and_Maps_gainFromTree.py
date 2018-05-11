
import ROOT



nomefile="out_dataFF_0.08_new_readHisto_dev3.root"
nomefileOut="plots_mappe_dataFF_0.08_dev3.root"

infile =ROOT.TFile(nomefile)

tree=infile.Get("tree")


print "entries = ",tree.GetEntries()

h_mean=      ROOT.TH1F ("h_mean","",500,0,250)
h_meanG=     ROOT.TH1F ("h_meanG","",500,0,250)
h_meanCic=   ROOT.TH1F ("h_meanCic","",500,0,250)
h_meanCic3=  ROOT.TH1F ("h_meanCic3","",500,0,250)
h_gaussAve=  ROOT.TH1F ("h_gaussAve","",500,0,250)
h_landauAve= ROOT.TH1F ("h_landauAve","",500,0,250)


h2mean=ROOT.TH2F("h2mean","mean",300,0,300,352,0,352)
h2meanG=ROOT.TH2F("h2meanG","meanG",300,0,300,352,0,352)
h2meanCic=ROOT.TH2F("h2meanCic","meanCic",300,0,300,352,0,352)
h2meanCic3=ROOT.TH2F("h2meanCic3","meanCic3",300,0,300,352,0,352)
h2gaussAve=ROOT.TH2F("h2gaussAve","gaussAve",300,0,300,352,0,352)
h2landauAve=ROOT.TH2F("h2landauAve","landauAve",300,0,300,352,0,352)


h_mean.SetTitle(nomefile[:-5])


h_mean.SetLineColor(2)
h_meanG.SetLineColor(4)
h_meanCic.SetLineColor(1)
h_meanCic3.SetLineColor(6)
h_gaussAve.SetLineColor(7)
h_mean.SetLineWidth(2)
h_meanG.SetLineWidth(2)
h_meanCic.SetLineWidth(2)
h_meanCic3.SetLineWidth(2)
h_gaussAve.SetLineWidth(2)
h_landauAve.SetLineWidth(2)



padding=5

CUT='pix_x>'+str(padding)+' && pix_x<300-'+str(padding)+' && pix_y>'+str(padding)+' && pix_y<352-'+str(padding)

print "CUT=",CUT

tree.Draw("mean>>h_mean", CUT )
tree.Draw("meanG>>h_meanG", CUT )
tree.Draw("meanCic>>h_meanCic", CUT )
tree.Draw("meanCic3>>h_meanCic3",  CUT)
tree.Draw("gaussAve>>h_gaussAve", CUT )
tree.Draw("landauAve>>h_landauAve", CUT)



tree.Draw("pix_y:pix_x>>h2mean", "mean")
tree.Draw("pix_y:pix_x>>h2meanG", "meanG")
tree.Draw("pix_y:pix_x>>h2meanCic", "meanCic")
tree.Draw("pix_y:pix_x>>h2meanCic3", "meanCic3")
tree.Draw("pix_y:pix_x>>h2gaussAve", "gaussAve")
tree.Draw("pix_y:pix_x>>h2landauAve", "landauAve")

risMean=0
risMeanG=0 
risMeanCic=0
risMeanCic3=0
risGaussAve=0
risLAndauAve=0

if h_mean.GetMean()!=0:
   risMean=h_mean.GetRMS()/h_mean.GetMean()
if h_meanG.GetMean()!=0:   
   risMeanG=h_meanG.GetRMS()/h_meanG.GetMean()
if h_meanCic.GetMean()!=0:
   risMeanCic=h_meanCic.GetRMS()/h_meanCic.GetMean()
if h_meanCic3.GetMean()!=0:
   risMeanCic3=h_meanCic3.GetRMS()/h_meanCic3.GetMean()
if h_gaussAve.GetMean()!=0:
   risGaussAve=h_gaussAve.GetRMS()/h_gaussAve.GetMean()
if h_landauAve.GetMean()!=0:
   risLandauAve=h_landauAve.GetRMS()/h_landauAve.GetMean()

c1=ROOT.TCanvas("c1","",0)


h_mean.Draw()
h_meanG.Draw("sames")
h_meanCic.Draw("sames")
h_meanCic3.Draw("sames")
h_gaussAve.Draw("sames")
h_landauAve.Draw("sames")


leg=ROOT.TLegend(0.7,0.94,0.99,0.5)
print "{:10.4f}".format(risMean)
label='mean; ris= '+str("{:2.2f}".format(risMean*100))+' %'
leg.AddEntry(h_mean,label,"l")
label='meanG; ris= '+str("{:2.2f}".format(risMeanG*100))+' %'
leg.AddEntry(h_meanG,label,"l")

label='meanCic; ris= '+str("{:2.2f}".format(risMeanCic*100))+' %'
leg.AddEntry(h_meanCic,label,"l")

label='meanCic3; ris= '+str("{:2.2f}".format(risMeanCic3*100))+' %'
leg.AddEntry(h_meanCic3,label,"l")

label='gaussAve; ris= '+str("{:2.2f}".format(risGaussAve*100))+' %'
leg.AddEntry(h_gaussAve,label,"l")

label='landauAve; ris= '+str("{:2.2f}".format(risLandauAve*100))+' %'
leg.AddEntry(h_landauAve,label,"l")
leg.Draw()



c2=ROOT.TCanvas("c2","",0)
c2.cd()
h2pixelN=infile.Get("h2pixelN")
h2pixelN.SetTitle(nomefile[:-5])
h2pixelN.Draw("colz")



#nomefile_histos='plotHistos_'+nomefile[:-5]+'.png'
#c1.SaveAs(nomefile_histos)

#nomefile_entries='plotNentries_'+nomefile[:-5]+'.png'
#c2.SaveAs(nomefile_entries)


#########################################3
# mappe 2D

c3=ROOT.TCanvas("c3","",0)
c3.Divide(3,2)
c3.cd(1)
h2mean.Draw("colz")
c3.cd(2)
h2meanG.Draw("colz")
c3.cd(3)
h2meanCic.Draw("colz")
c3.cd(4)
h2meanCic3.Draw("colz")
c3.cd(5)
h2gaussAve.Draw("colz")
c3.cd(6)
h2landauAve.Draw("colz")


outfile=ROOT.TFile(nomefileOut,"recreate")
c3.Write()
c2.Write()
c1.Write()

h2mean.Write()
h2meanG.Write()
h2meanCic.Write()
h2meanCic3.Write()
h2gaussAve.Write()
h2landauAve.Write()
h_mean.Write()
h_meanG.Write()
h_meanCic.Write()
h_meanCic3.Write()
h_gaussAve.Write()
h_landauAve.Write()

outfile.Close()
