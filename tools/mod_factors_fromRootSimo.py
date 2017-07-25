# small script to fit the modulation factor vs energy.
# imput files, energies and name of output files are hardcoded!!!




import ROOT
import array as arr

path='/home/maldera/FERMI/Xipe/rec/data/sims/'
filenames=['sim_2keV_zero_5_th1_5_th2_5.root', 'outSimo2new_sim3Kev.root', 'sim_4keV_zero_5_th1_5_th2_5.root','outSimo2new_sim5Kev.root','sim_6keV_zero_5_th1_5_th2_5.root', 'outSimo2new_sim7Kev.root', 'sim_8keV_zero_5_th1_5_th2_5.root'   ,'outSimo2new_sim9Kev.root']


energy=arr.array("d",[2,3,4,5,6,7,8,9])
PA_cut=arr.array('d',[1200,2500,3400,4600])



modFactors=arr.array("d",[0.]*8)
err_mu=arr.array("d",[0.]*8)
err_e=arr.array("d",[0.]*8)

outFile = ROOT.TFile("outFile.root","recreate")



for i in range (0, len(energy)):

     
    filename=path+filenames[i]
    print "\n\n i= ",i," processing file ",filename, " E = ",energy[i]
    


    inRootFile=ROOT.TFile(filename)
    tree=inRootFile.Get("reconT")
    name="phi_"+str(energy[i])+"KeV"
    c1=ROOT.TCanvas(name,name,0)   
    h=ROOT.TH1F("h","",400,-4,4)

    #cut='trk_size>26 && trk_pulse_height > '+str(PA_cut[i])
    cut='trk_size>26 && trk_pulse_height > 0'
    
    print 'cut = ',cut 

    hShape=ROOT.TH1F("hShape","",1000,0,1)
    tree.Draw("trk_mom2trans/trk_mom2long >> hShape",cut,"")
    #ntot=hShape.GetEntries()
    ntot=tree.GetEntries()
   
    shapeCut=0
    for jj in range (0,1000):
        
        integral=hShape.Integral(1000-jj,1000)
        print "x0 = ",1000-jj," integral = ",integral/ntot
        if (integral/ntot >0.2):
            binCut=1000-jj
            shapeCut=hShape.GetBinCenter(binCut)
            break
            
        
    print "shapecut = ", shapeCut
   
    cut=cut+' &&  trk_mom2trans/trk_mom2long <  '+str(shapeCut)
    
    print "cut con shape = ",cut
    #valore = raw_input('continue?')
    tree.Draw("trk_phi*TMath::DegToRad() >> h",cut, "")
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
g= ROOT.TGraphErrors(8,energy,modFactors,err_e,err_mu)
g.SetTitle("mod. factor vs energy")
g.SetMarkerStyle(20)
g.SetMarkerColor(2)
g.SetLineColor(2)
g.SetLineWidth(2)
g.GetXaxis().SetTitle("[KeV]")
g.GetYaxis().SetTitle("mod. factor")

g.Draw("ap")

c.Write()
g.SetName("mu_std_qcut")
g.Write()
outFile.Close()



