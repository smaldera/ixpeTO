# small script to fit the modulation factor vs energy.
# imput files, energies and name of output files are hardcoded!!!




import ROOT
import array as arr


filenames=['outSimo2new_sim3Kev.root','outSimo2new_sim5Kev.root','outSimo2new_sim7Kev.root','outSimo2new_sim9Kev.root']
energy=arr.array("d",[3,5,7,9])

modFactors=arr.array("d",[0.]*4)
err_mu=arr.array("d",[0.]*4)
err_e=arr.array("d",[0.]*4)

outFile = ROOT.TFile("outFile.root","recreate")



for i in range (0, len(energy)):



    print "\n\n i= ",i," processing file ",filenames[i], " E = ",energy[i]
    


    inRootFile=ROOT.TFile(filenames[i])
    tree=inRootFile.Get("reconT")
    name="phi_"+str(energy[i])+"KeV"
    c1=ROOT.TCanvas(name,name,0)   
    h=ROOT.TH1F("h","",400,-4,4)

    tree.Draw("trk_phi*TMath::DegToRad() >> h","", "")
    #tree.Scan("trk_phi")

    h.Draw()
    #valore = raw_input('continue?')  

    fmod=ROOT.TF1("fmod","[0]+[1]*(pow(cos(x+[2]),2))",-3.1415,3.1415)
    fmod.SetLineColor(4)
    fmod.SetParLimits(0,0,1e9)
    fmod.SetParLimits(1,0,1e9)
    
    h.Fit("fmod","MER")
    outFile.cd()
    c1.Write()
    
    A=fmod.GetParameter(0)
    B=fmod.GetParameter(1)
    phi0=fmod.GetParameter(2)

    Aerr=fmod.GetParError(0)
    Berr=fmod.GetParError(1)    
    modFactors[i]=B/(2.*A+B)

    dmu_dA=-2.*B/((2.*A+B)**2)
    dmu_dB=2.*A/((2.*A+B)**2)

    err_mu[i]=(  ((dmu_dA*Aerr)**2)+((dmu_dB*Berr)**2))**0.5
    

    print "modFact= ", modFactors[i]
    #valore = raw_input('continue?')    
    h.Delete()
    #c1.Delete()
    fmod.Delete()


c=ROOT.TCanvas("modFact","modFact",0)
g= ROOT.TGraphErrors(4,energy,modFactors,err_e,err_mu)
g.SetTitle("mod. factor vs energy")
g.SetMarkerStyle(20)
g.SetMarkerColor(2)
g.SetLineColor(2)
g.SetLineWidth(2)
g.GetXaxis().SetTitle("[KeV]")
g.GetYaxis().SetTitle("mod. factor")

g.Draw("ap")

c.Write()


outFile.Close()



