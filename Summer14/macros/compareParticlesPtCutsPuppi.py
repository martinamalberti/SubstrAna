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


def closestToZero(a):
    tmp = 9999
    ind = -1
    for k in range(0, len(a)):
        #print a[k][1], tmp
        if (abs(a[k][1]) < tmp ):
            tmp = abs(a[k][1])
            ind = k
    cut = a[ind][0]
    #print 'best cut = ', cut, '  diff = ', abs(a[ind][1])  
    return cut


def getZeros(a):
    x = []
    for i in range(0, len(a)-1):
        if ( (a[i][1]*a[i+1][1]<0) or (a[i][1]>0 and a[i+1][1]>0 and i ==  len(a)-2) ):
            gtemp = ROOT.TGraph()
            gtemp.SetPoint(0,a[i][0],a[i][1])
            gtemp.SetPoint(1,a[i+1][0],a[i+1][1])
            ftemp = ROOT.TF1("ftemp","[0]+[1]*x")
            gtemp.Fit("ftemp","Q")

            if (a[i][0]>a[i+1][0]):
                x1 = a[i+1][0]
                x2 = a[i][0]
            else:
                x2 = a[i+1][0]
                x1 = a[i][0]

            xtemp = - ftemp.GetParameter(0)/ftemp.GetParameter(1)
            
          
            if ((a[i][1]*a[i+1][1]<0) and  xtemp > x1 and xtemp < x2):
                x.append(xtemp)
            elif ( a[i][1]>0 and a[i+1][1]>0 and i ==  len(a)-2):
                x.append(xtemp)

    if (len(x)<1):
        return closestToZero(a)
    else:
        return x[-1]


region = options.region

fnames = []

if (region == 'central'):
    fnames = ['../scripts/hQCD300to470_MinNeutralPt00_central/histograms.root',
              '../scripts/hQCD300to470_MinNeutralPt05_central/histograms.root',
              '../scripts/hQCD300to470_MinNeutralPt10_central/histograms.root',
              '../scripts/hQCD300to470_MinNeutralPt15_central/histograms.root',
              '../scripts/hQCD300to470_MinNeutralPt20_central/histograms.root',
              '../scripts/hQCD300to470_MinNeutralPt25_central/histograms.root',
              ]
elif (region == 'forward'):
    fnames = ['../scripts/hQCD300to470_MinNeutralPt00_forward/histograms.root',
              '../scripts/hQCD300to470_MinNeutralPt05_forward/histograms.root',
              '../scripts/hQCD300to470_MinNeutralPt10_forward/histograms.root',
              '../scripts/hQCD300to470_MinNeutralPt15_forward/histograms.root',
              '../scripts/hQCD300to470_MinNeutralPt20_forward/histograms.root',
              '../scripts/hQCD300to470_MinNeutralPt25_forward/histograms.root',
              ]

if (region == 'central'):
    fnames = ['../scripts/hQCD300to470_MinNeutralPt00_JetPt25to100_central/histograms.root',
              '../scripts/hQCD300to470_MinNeutralPt05_JetPt25to100_central/histograms.root',
              '../scripts/hQCD300to470_MinNeutralPt10_JetPt25to100_central/histograms.root',
              '../scripts/hQCD300to470_MinNeutralPt15_JetPt25to100_central/histograms.root',
              '../scripts/hQCD300to470_MinNeutralPt20_JetPt25to100_central/histograms.root',
              '../scripts/hQCD300to470_MinNeutralPt25_JetPt25to100_central/histograms.root',
              ]
elif (region == 'forward'):
    fnames = ['../scripts/hQCD300to470_MinNeutralPt00_JetPt25to100_forward/histograms.root',
              '../scripts/hQCD300to470_MinNeutralPt05_JetPt25to100_forward/histograms.root',
              '../scripts/hQCD300to470_MinNeutralPt10_JetPt25to100_forward/histograms.root',
              '../scripts/hQCD300to470_MinNeutralPt15_JetPt25to100_forward/histograms.root',
              '../scripts/hQCD300to470_MinNeutralPt20_JetPt25to100_forward/histograms.root',
              '../scripts/hQCD300to470_MinNeutralPt25_JetPt25to100_forward/histograms.root',
              ]



ptcuts = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5]
colors = [1,2,4,5,6,8]

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
            leg.AddEntry(h,'p_{T,neutr} > %.1f GeV'%ptcuts[i])
            
        if ('response' in hname):
            func.SetRange(h.GetMean()-2*h.GetRMS(),h.GetMean()+2*h.GetRMS())
            h.Fit("func","QRN")
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
#for typ in '.png','.pdf':
#    c.SaveAs("efficiency_vs_realfraction"+typ)


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


#c.SaveAs("efficiency.png")



# mean and rms vs cut value
leg2 = ROOT.TLegend(0.7,0.7,0.93,0.9);
leg2.SetBorderSize(0);
leg2.SetFillStyle(0);
leg2.AddEntry(gptrms,"rms","LP")
leg2.AddEntry(gptsigma,"fitted #sigma","LP")

cptmean = ROOT.TCanvas("pt_mean_vs_cut","pt_mean_vs_cut")
cptmean.SetGridx()
cptmean.SetGridy()
gptmean.GetHistogram().GetXaxis().SetTitle("pT cut (GeV)")
gptmean.GetHistogram().GetYaxis().SetTitle("<#DeltapT> (GeV)")
gptmean.SetMinimum(-50)
gptmean.SetMaximum(50)
gptmean.Draw("apl")
gptfittedmean.Draw("plsame")

cptrms = ROOT.TCanvas("pt_rms_vs_cut","pt_rms_vs_cut")
gptrms.GetHistogram().GetXaxis().SetTitle("pT cut (GeV)")
gptrms.GetHistogram().GetYaxis().SetTitle("RMS(#DeltapT) (GeV)")
gptrms.SetMinimum(5)
gptrms.SetMaximum(50)
gptrms.Draw("apl")
gptsigma.Draw("plsame")
leg2.Draw()

cmassmean = ROOT.TCanvas("mass_mean_vs_cut","mass_mean_vs_cut")
cmassmean.SetGridx()
cmassmean.SetGridy()
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
leg2.Draw()

for canv in cptmean, cptrms, cmassmean, cmassrms:
    for typ in '.png','.pdf':
        canv.SaveAs(canv.GetName()+typ)



# response vs pT cut in different nPU bins
mbestcut = {}
ptbestcut = {}
gmbestcut = ROOT.TGraphErrors()
gptbestcut = ROOT.TGraphErrors()

rebin = 10
for var in 'pt','m':
    h = []
    response   = {} # {nPU:[cut, mean, meanerr]}
    resolution = {} # {nPU:[cut, rms, rmserr]}

    for i,fname in enumerate(fnames):
        f = ROOT.TFile.Open(fname)
        h.append( f.Get('puppi/h'+var+'_response_vs_npv_puppi') )
        j = 0 # j -> pu bin
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
                if npu not in response:
                    response[npu] = []
                    resolution[npu] = []
                response[npu].append([ptcuts[i],mean,meanerr])
                resolution[npu].append([ptcuts[i],rms,rmserr])

    cresponse = ROOT.TCanvas("cresponse"+var,"cresponse"+var)
    cresolution = ROOT.TCanvas("cresolution"+var,"cresolution"+var)
    cresponse.SetGridx()
    cresponse.SetGridy()

    gresponse = []
    gresolution = []
    nbinspu = 0
    leg3 = ROOT.TLegend(0.2,0.2,0.5,0.5)
    leg3.SetBorderSize(0)
    leg3.SetFillStyle(0)
    
    for npu in response:
        if (npu > 10 and npu < 50 ):
            gresponse.append(ROOT.TGraphErrors())
            gresolution.append(ROOT.TGraphErrors())
            for k,el in enumerate(response[npu]):
                gresponse[nbinspu] .SetPoint(k,response[npu][k][0], response[npu][k][1])
                gresponse[nbinspu] .SetPointError(k,0, response[npu][k][2])
                gresolution[nbinspu] .SetPoint(k,resolution[npu][k][0], resolution[npu][k][1])
                gresolution[nbinspu] .SetPointError(k,0, resolution[npu][k][2])
        
            gresponse[nbinspu].SetLineColor(nbinspu+1)    
            gresponse[nbinspu].SetMarkerColor(nbinspu+1)    
            gresolution[nbinspu].SetLineColor(nbinspu+1)    
            gresolution[nbinspu].SetMarkerColor(nbinspu+1)    
            leg3.AddEntry(gresponse[nbinspu],"nPV=[%d,%d]"%(npu-rebin/2,npu+rebin/2),"L")

            if (nbinspu == 0):
                if (var == 'pt'):
                    gresponse[nbinspu].SetMinimum(-0.3)
                    gresponse[nbinspu].SetMaximum(0.3)
                    gresponse[nbinspu].GetHistogram().GetXaxis().SetTitle("p_{T,neutr} cut (GeV)")
                    gresponse[nbinspu].GetHistogram().GetYaxis().SetTitle("< (pT-gen pT) /gen pT>")
                    gresolution[nbinspu].SetMinimum(0.0)
                    gresolution[nbinspu].SetMaximum(0.3)
                    gresolution[nbinspu].GetHistogram().GetXaxis().SetTitle("p_{T,neutr} cut (GeV)")
                    gresolution[nbinspu].GetHistogram().GetYaxis().SetTitle("RMS (pT-gen pT) /gen pT")
                else:
                    gresponse[nbinspu].SetMinimum(-30)
                    gresponse[nbinspu].SetMaximum(30)
                    gresponse[nbinspu].GetHistogram().GetXaxis().SetTitle("p_{T,neutr} cut (GeV)")
                    gresponse[nbinspu].GetHistogram().GetYaxis().SetTitle("< (m-gen m) > (GeV)")
                    gresolution[nbinspu].SetMinimum(0)
                    gresolution[nbinspu].SetMaximum(30)
                    gresolution[nbinspu].GetHistogram().GetXaxis().SetTitle("p_{T,neutr} cut (GeV)")
                    gresolution[nbinspu].GetHistogram().GetYaxis().SetTitle("RMS(m - gen m) (GeV)")
                cresponse.cd()
                gresponse[nbinspu].Draw("apl")
                cresolution.cd()
                gresolution[nbinspu].Draw("apl")
            else: 
                cresponse.cd()
                gresponse[nbinspu] .Draw("plsame")
                cresolution.cd()
                gresolution[nbinspu] .Draw("plsame")


            if (var=='pt'):
                #ptbestcut[npu] = closestToZero(response[npu])
                ptbestcut[npu] = getZeros(response[npu])
                print var, npu, ptbestcut[npu]
                gptbestcut.SetPoint(nbinspu, npu , ptbestcut[npu] )
                gptbestcut.SetPointError(nbinspu, 0. , 0.25 )
            else:
                print "*** MASS "
                #mbestcut[npu] = closestToZero(response[npu])
                mbestcut[npu] = getZeros(response[npu])
                uffa =  getZeros(response[npu])
                print "closest to zero:      ", npu, mbestcut[npu]
                print "intersection with zero: ", npu, uffa
                gmbestcut.SetPoint(nbinspu, npu , mbestcut[npu],)
                gmbestcut.SetPointError(nbinspu, 0. , 0.25 )

            nbinspu = nbinspu+1

    cresponse.cd()
    leg3.Draw()

    cresolution.cd()
    leg3.Draw()

    
    for typ in '.png','.pdf':
        cresponse.SaveAs(var+"_response_vs_npu"+typ)
        cresolution.SaveAs(var+"_resolution_vs_npu"+typ)

cbest=ROOT.TCanvas("cbest","cbest",700,700)
gmbestcut.SetMarkerColor(2)
gmbestcut.SetMarkerStyle(25)
gmbestcut.SetLineColor(2)

gptbestcut.SetMarkerColor(4)
gptbestcut.SetLineColor(4)

gmbestcut.GetHistogram().GetXaxis().SetTitle("n_{PV}")
gmbestcut.GetHistogram().GetYaxis().SetTitle("optimal p_{T, neutr} cut (GeV)")
gmbestcut.GetHistogram().GetYaxis().SetRangeUser(-1.,9.)
gmbestcut.Draw("ap")
gptbestcut.Draw("psame")
ffit = ROOT.TF1("ffit","pol1")
ffit.SetLineColor(2)
ffit.SetLineStyle(2)
gmbestcut.Fit("ffit")


leg4 = ROOT.TLegend(0.2,0.7,0.43,0.9);
leg4.SetBorderSize(0);
leg4.SetFillStyle(0);
leg4.AddEntry(gmbestcut,"mass","P")
leg4.AddEntry(gptbestcut,"pt","P")
leg4.Draw()

for typ in '.png','.pdf':
    cbest.SaveAs("optimal_ptNeutralCut"+typ)

raw_input('ok?')
