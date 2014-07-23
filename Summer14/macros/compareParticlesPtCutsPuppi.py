#! /usr/bin/env python                                                                                                                                                                                         
import os
import glob
import math
import array
import sys
import time

import ROOT

ROOT.gROOT.ProcessLine(".L ~/tdrstyle.C");
ROOT.setTDRStyle();
#ROOT.gStyle.SetOptStat(1111)


from optparse import OptionParser
parser = OptionParser()
parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
parser.add_option('-r','--region',action="store",type="string",dest="region",default="all")
(options, args) = parser.parse_args()


    
def makeEff(f,var):
    hnum = f.Get('puppi/h'+var+'_good_puppi').Clone('hnum'+var)
    hden = f.Get('gen/h'+var+'_good_gen')
    if (hnum.GetNbinsX()!= hden.GetNbinsX()):
        rebin = hnum.GetNbinsX()/hden.GetNbinsX()
        hnum.Rebin(rebin)
    else:
        rebin = 2
        if 'pt' in var:
            rebin = 5
        hnum.Rebin(rebin)
        hden.Rebin(rebin)
            
    hnum.Divide(hnum,hden,1.,1.)
    return hnum


def makeRealJetsFraction(f,var):
    var = var.replace('gen','')
    hnum = f.Get('puppi/h'+var+'_good_puppi').Clone('hnumf'+var)
    hden = f.Get('puppi/h'+var+'_puppi')
    if (hnum.GetNbinsX()!= hden.GetNbinsX()):
        rebin = hnum.GetNbinsX()/hden.GetNbinsX()
        hnum.Rebin(rebin)
    else:
        rebin = 2
        if 'pt' in var:
            rebin = 10
        hnum.Rebin(rebin)
        hden.Rebin(rebin)
    hnum.Divide(hnum,hden,1.,1.)
    return hnum


def makeEffVsRealJetsFraction(fnames):

    g= ROOT.TGraph()
    for i,fname in enumerate(fnames):
        f = ROOT.TFile.Open(fname);
        h     = f.Get("puppi/hpt_puppi").Clone("h")
        hgood = f.Get("puppi/hpt_good_puppi").Clone("hgood")
        hgen  = f.Get("gen/hpt_gen").Clone("hgen")
        eff  = hgood.GetSumOfWeights()/hgen.GetSumOfWeights()
        frac = hgood.GetSumOfWeights()/h.GetSumOfWeights()
        g.SetPoint(i, frac, eff)
    g.SetMinimum(0)
    g.SetMaximum(1)    
    g.GetHistogram().GetXaxis().SetTitle("fraction of real jets")
    g.GetHistogram().GetYaxis().SetTitle("jet efficiency")
    return g



region = options.region

fnames = []

if (region == 'central'):
    fnames = ['../bin/histograms/histograms_ttbar_r08_minpartpt0.0GeV_minpt25_central.root',
              '../bin/histograms/histograms_ttbar_r08_minpartpt0.5GeV_minpt25_central.root',
              '../bin/histograms/histograms_ttbar_r08_minpartpt1.0GeV_minpt25_central.root',
              '../bin/histograms/histograms_ttbar_r08_minpartpt1.5GeV_minpt25_central.root'
              ]
elif (region == 'forward'):
    fnames = ['../bin/histograms/histograms_ttbar_r08_minpartpt0.0GeV_minpt25_fwd.root',
              '../bin/histograms/histograms_ttbar_r08_minpartpt0.5GeV_minpt25_fwd.root',
              '../bin/histograms/histograms_ttbar_r08_minpartpt1.0GeV_minpt25_fwd.root',
              '../bin/histograms/histograms_ttbar_r08_minpartpt1.5GeV_minpt25_fwd.root'
              ]
else:
    fnames = ['../bin/histograms/histograms_ttbar_r08_minpartpt0.0GeV_minpt25.root',
              '../bin/histograms/histograms_ttbar_r08_minpartpt0.5GeV_minpt25.root',
              '../bin/histograms/histograms_ttbar_r08_minpartpt1.0GeV_minpt25.root',
              '../bin/histograms/histograms_ttbar_r08_minpartpt1.5GeV_minpt25.root'
              ]



ptcuts = [0.0, 0.5, 1.0, 1.5]
colors = [1,2,4,8]

histograms = { 'heta':['#eta','events',2], 
               'hptgen':['p_{T} (GeV)','events',5],
               'hpt_response':['p_{T} - gen p_{T} (GeV)','events',4],
               'hm_response':['m - gen m (GeV)','events',4]
              }


leg = ROOT.TLegend(0.7,0.7,0.93,0.9);
leg.SetBorderSize(0);
leg.SetFillStyle(0);

gptmean= ROOT.TGraphErrors()
gmassmean= ROOT.TGraphErrors()
gptrms= ROOT.TGraphErrors()
gmassrms= ROOT.TGraphErrors()
gptfittedmean= ROOT.TGraphErrors()
gmassfittedmean= ROOT.TGraphErrors()
gptsigma = ROOT.TGraphErrors()
gmasssigma = ROOT.TGraphErrors() 

gptfittedmean.SetMarkerStyle(24)
gmassfittedmean.SetMarkerStyle(24)
gptsigma.SetMarkerStyle(24)
gmasssigma.SetMarkerStyle(24)


func = ROOT.TF1("func","gaus",3)

n = 0 
for hname in histograms:
    c = ROOT.TCanvas(hname.replace('h',''),hname.replace('h',''))

    for i,fname in enumerate(fnames):
        f = ROOT.TFile.Open(fname)
        h = f.Get('puppi/'+hname+'_puppi')
        h.GetXaxis().SetTitle(histograms[hname][0])
        h.GetYaxis().SetTitle(histograms[hname][1])
        h.Rebin(histograms[hname][2])
        h.SetLineColor(colors[i])

        if (n==0):
            leg.AddEntry(h,'p_{T,part} > %.1f GeV'%ptcuts[i])
            
        if ('response' in hname):
            func.SetRange(h.GetMean()-2*h.GetRMS(),h.GetMean()+2*h.GetRMS())
            h.Fit("func","RN")
            if ('pt' in hname):
                gptmean.SetPoint(i,ptcuts[i],h.GetMean())
                gptmean.SetPointError(i,0,h.GetMeanError())
                gptrms.SetPoint(i,ptcuts[i], h.GetRMS())
                gptrms.SetPointError(i,0, h.GetRMSError())
                gptfittedmean.SetPoint(i,ptcuts[i], func.GetParameter(1))
                gptfittedmean.SetPointError(i,0, func.GetParError(1))
                gptsigma.SetPoint(i,ptcuts[i], func.GetParameter(2))
                gptsigma.SetPointError(i,0, func.GetParError(2))
            else:
                gmassmean.SetPoint(i,ptcuts[i],h.GetMean())
                gmassmean.SetPointError(i,0,h.GetMeanError())
                gmassrms.SetPoint(i,ptcuts[i], h.GetRMS())
                gmassrms.SetPointError(i,0, h.GetRMSError())
                gmassfittedmean.SetPoint(i,ptcuts[i], func.GetParameter(1))
                gmassfittedmean.SetPointError(i,0, func.GetParError(1))
                gmasssigma.SetPoint(i,ptcuts[i], func.GetParameter(2))
                gmasssigma.SetPointError(i,0, func.GetParError(2))

        if (i==0):
            c.cd()
            print c.GetName()
            h.Draw()
        else:
            c.cd()
            h.Draw("sames")
    leg.Draw()
    n= n+1
    
    for typ in '.png','.pdf':
        c.SaveAs(c.GetName()+typ)


# eff vs real jets fraction
c = ROOT.TCanvas()
g = makeEffVsRealJetsFraction(fnames)
g.Draw("ap")
for typ in '.png','.pdf':
    c.SaveAs("efficiency_vs_realfraction"+typ)


# eff and real jets fractions
#for var in ['ptgen','eta']:
for var in ['ptgen']:
    c  = ROOT.TCanvas('c_eff_'+var,'c_eff_'+var)
    cc = ROOT.TCanvas('c_frac_'+var,'c_frac_'+var)
    heff  = []
    hfrac = []
    for i,fname in enumerate(fnames):
        f = ROOT.TFile.Open(fname)
        heff.append(makeEff(f,var))
        heff[i].SetLineColor(colors[i])
        hfrac.append(makeRealJetsFraction(f,var))
        hfrac[i].SetLineColor(colors[i])
        if 'pt' in var:
            heff[i].GetXaxis().SetRangeUser(25,300)
            heff[i].GetXaxis().SetRangeUser(25,300)
            hfrac[i].GetXaxis().SetRangeUser(25,300)
            hfrac[i].GetXaxis().SetRangeUser(25,300)
        heff[i].GetYaxis().SetRangeUser(0,1.4)
        heff[i].GetYaxis().SetRangeUser(0,1.4)
        hfrac[i].GetYaxis().SetRangeUser(0,1.4)
        hfrac[i].GetYaxis().SetRangeUser(0,1.4)

        if (i==0):
            c.cd()
            heff[i].Draw()
            leg.Draw()
            cc.cd()
            hfrac[i].Draw()
            leg.Draw()
        else:
            c.cd()
            heff[i].Draw("same")
            cc.cd()
            hfrac[i].Draw("same")


c.SaveAs("efficiency.png")



# mean and rms vs cut value
cptmean = ROOT.TCanvas("pt_mean_vs_cut","pt_mean_vs_cut")
gptmean.GetHistogram().GetXaxis().SetTitle("pT cut (GeV)")
gptmean.GetHistogram().GetYaxis().SetTitle("<#DeltapT> (GeV)")
gptmean.SetMinimum(-15)
gptmean.SetMaximum(15)
gptmean.Draw("apl")
gptfittedmean.Draw("plsame")

cptrms = ROOT.TCanvas("pt_rms_vs_cut","pt_rms_vs_cut")
gptrms.GetHistogram().GetXaxis().SetTitle("pT cut (GeV)")
gptrms.GetHistogram().GetYaxis().SetTitle("RMS(#DeltapT) (GeV)")
gptrms.SetMinimum(5)
gptrms.SetMaximum(20)
gptrms.Draw("apl")
gptsigma.Draw("plsame")

cmassmean = ROOT.TCanvas("mass_mean_vs_cut","mass_mean_vs_cut")
gmassmean.GetHistogram().GetXaxis().SetTitle("pT cut (GeV)")
gmassmean.GetHistogram().GetYaxis().SetTitle("<#Deltamass> (GeV)")
gmassmean.SetMinimum(-15)
gmassmean.SetMaximum(15)
gmassmean.Draw("apl")
gmassfittedmean.Draw("plsame")

cmassrms = ROOT.TCanvas("mass_rms_vs_cut","mass_rms_vs_cut")
gmassrms.GetHistogram().GetXaxis().SetTitle("pT cut (GeV)")
gmassrms.GetHistogram().GetYaxis().SetTitle("RMS(#Deltamass) (GeV)")
gmassrms.SetMinimum(5)
gmassrms.SetMaximum(20)
gmassrms.Draw("apl")
gmasssigma.Draw("plsame")

for canv in cptmean, cptrms, cmassmean, cmassrms:
    for typ in '.png','.pdf':
        canv.SaveAs(canv.GetName()+typ)



# response vs pT cut in different nPU bins
ptmean = {} # {nPU:[cut,mean, meanerr]}
ptrms  = {} # {nPU:[cut,mean, meanerr]}

h = []
rebin = 10
for i,fname in enumerate(fnames):
    f = ROOT.TFile.Open(fname)
    h.append( f.Get('puppi/hpt_response_vs_npu_puppi') )
    print h[i].GetName()
    j = 0 # pu bin
    px = (h[i].ProjectionX('px')).Clone('px') # projection on x 
    for bin in range(1,h[i].GetNbinsX()+1,rebin): # rebin 
         firstbin = bin
         lastbin  = bin+(rebin-1)
         py       = (h[i].ProjectionY('py',firstbin,lastbin)).Clone('py')
         mean     = py.GetMean()
         meanerr  = py.GetMeanError()
         rms      = py.GetRMS()
         rmserr   = py.GetRMSError()
         if (py.GetEntries()>0):
             npu  = (px.GetXaxis().GetBinUpEdge(lastbin) + px.GetXaxis().GetBinLowEdge(firstbin))/2
             npuerr = (px.GetXaxis().GetBinUpEdge(lastbin) - px.GetXaxis().GetBinLowEdge(firstbin))/2
             if npu not in ptmean:
                 ptmean[npu] = []
                 ptrms[npu] = []
             ptmean[npu].append([ptcuts[i],mean,meanerr])
             ptrms[npu].append([ptcuts[i],rms,rmserr])
             
cmean = ROOT.TCanvas("cmean","cmean")
crms = ROOT.TCanvas("crms","crms")
gmean = []
grms = []
ngraphs = 0
leg2 = ROOT.TLegend(0.2,0.2,0.5,0.5)
leg2.SetBorderSize(0)
leg2.SetFillStyle(0)
for npu in ptmean:
    if (npu < 60 ):
        gmean.append(ROOT.TGraphErrors())
        grms.append(ROOT.TGraphErrors())
        for i,el in enumerate(ptmean[npu]):
            gmean[ngraphs] .SetPoint(i,ptmean[npu][i][0], ptmean[npu][i][1])
            gmean[ngraphs] .SetPointError(i,0, ptmean[npu][i][2])
            grms[ngraphs] .SetPoint(i,ptrms[npu][i][0], ptrms[npu][i][1])
            grms[ngraphs] .SetPointError(i,0, ptrms[npu][i][2])
        
        gmean[ngraphs].SetLineColor(ngraphs+1)    
        gmean[ngraphs].SetMarkerColor(ngraphs+1)    
        grms[ngraphs].SetLineColor(ngraphs+1)    
        grms[ngraphs].SetMarkerColor(ngraphs+1)    
        leg2.AddEntry(gmean[ngraphs],"nPU=[%d,%d]"%(npu-5,npu+5),"L")
        if (ngraphs == 0):
            cmean.cd()
            gmean[ngraphs].SetMinimum(-0.5)
            gmean[ngraphs].SetMaximum(0.5)
            gmean[ngraphs].GetHistogram().GetXaxis().SetTitle("p_{T,part} cut (GeV)")
            gmean[ngraphs].GetHistogram().GetYaxis().SetTitle("< (pT-gen pT) /gen pT>")
            gmean[ngraphs].Draw("apl")
            crms.cd()
            grms[ngraphs].SetMinimum(-0.5)
            grms[ngraphs].SetMaximum(0.5)
            grms[ngraphs].GetHistogram().GetXaxis().SetTitle("p_{T,part} cut (GeV)")
            grms[ngraphs].GetHistogram().GetYaxis().SetTitle("RMS (pT-gen pT) /gen pT")
            grms[ngraphs].Draw("apl")
        else: 
            cmean.cd()
            gmean[ngraphs] .Draw("plsame")
            crms.cd()
            grms[ngraphs] .Draw("plsame")
        ngraphs = ngraphs+1

cmean.cd()
leg2.Draw()

crms.cd()
leg2.Draw()


raw_input('ok?')
