import ROOT


fileDqmGainMap=ROOT.TFile('102_0000057_data_reconGainMapDqm.root','open')
#fileSimoGainMap=ROOT.TFile('102_0000057_data_reconGainMapSimo.root','open')
fileSimoGainMap=ROOT.TFile('102_0000057_data_reconGainMapSimo_mean0.20.root','open')


reconTdqm=fileDqmGainMap.Get("reconT")
reconTsimo=fileSimoGainMap.Get("reconT")


hPHA=ROOT.TH1F("hPHA","",1000,0,10000)
hPIdqm=ROOT.TH1F("hPIdqm","",1000,0,10000)
hPIsimo=ROOT.TH1F("hPIsimo","",1000,0,10000)

hPIdqm.SetLineColor(4)
hPIsimo.SetLineColor(2)

#CUT=""
CUT="TRK_ABSX>(-7.48+0.7) && TRK_ABSX<(7.48-0.7)&& TRK_ABSY<(7.46-0.7)&&TRK_ABSY>(-7.46+0.7) && NUM_CLU<30&& ROI_SIZE<1000"
#CUT="TRK_ABSX>(-2) && TRK_ABSX<(2)&& TRK_ABSY<(2)&&TRK_ABSY>(-2)"

reconTdqm.Draw("TRK_PHA>>hPHA",CUT)
reconTdqm.Draw("TRK_PI>>hPIdqm",CUT)
reconTsimo.Draw("TRK_PI>>hPIsimo",CUT)


hPHA.Draw()
hPIdqm.Draw("sames")
hPIsimo.Draw("sames")

leg=ROOT.TLegend(0.7,0.94,0.99,0.5)
leg.AddEntry(hPHA,"PHA","l")
leg.AddEntry(hPIdqm,"PI dqm","l")
leg.AddEntry(hPIsimo,"PI simo","l")
leg.Draw()
