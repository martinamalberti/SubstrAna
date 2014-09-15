#! /usr/bin/env python
import os
import glob
import math
import array
import sys
import time

import ROOT

import CMS_lumi, tdrstyle

#set the tdr style
tdrstyle.setTDRStyle()

#change the CMS_lumi variables (see CMS_lumi.py)
CMS_lumi.lumi_13TeV = ""
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Simulation Preliminary"

#ROOT.gROOT.ProcessLine(".L ~/tdrstyle.C");
#ROOT.setTDRStyle();
#ROOT.gStyle.SetPadRightMargin(0.03);
##ROOT.gStyle.SetPadRightMargin(0.16);
ROOT.gStyle.SetLegendFont(42)
ROOT.gStyle.SetTextFont(42)

############################################
#            Job steering                  #
############################################
from optparse import OptionParser

parser = OptionParser()

parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
parser.add_option('-i','--input',action="store",type="string",dest="input",default="histograms.root")
parser.add_option('-o','--outdir',action="store",type="string",dest="outdir",default="plots")
parser.add_option('--nPU',action="store",type="int",dest="nPU",default=40)
parser.add_option('-r',action="store",type="float",dest="radius",default=0.8)
parser.add_option('--minPt',action="store",type="float",dest="minPt",default=25.)
parser.add_option('--maxPt',action="store",type="float",dest="maxPt",default=300.)
parser.add_option('--minEta',action="store",type="float",dest="minEta",default=0.)
parser.add_option('--maxEta',action="store",type="float",dest="maxEta",default=2.5)
parser.add_option('--sample',action="store",type="string",dest="sample",default="QCD")

(options, args) = parser.parse_args()



# cms preliminary 
cmsprel = ROOT.TLatex(0.20,0.96,("CMS Simulation Preliminary, #sqrt{s} = 13 TeV"))
cmsprel.SetNDC()
cmsprel.SetTextSize(0.03)

mytextsize = 0.04
latex0 = ROOT.TLatex(0.20,0.89,("%s"%(options.sample)))
latex0.SetNDC()
latex0.SetTextSize(mytextsize)
latex1 = ROOT.TLatex(0.20,0.84,("Anti-kT (R=%.1f)"%(options.radius)))
latex1.SetNDC()
latex1.SetTextSize(mytextsize)
latex2 = ROOT.TLatex(0.20,0.79,("<n_{PU}> = "+str(options.nPU)))
latex2.SetNDC()
latex2.SetTextSize(mytextsize)
latex3 = ROOT.TLatex(0.20,0.74,("%.0f GeV < p_{T} < %.0f GeV "%(options.minPt,options.maxPt)))
latex3.SetNDC()
latex3.SetTextSize(mytextsize)
latex4 = ROOT.TLatex(0.20,0.69,("%.1f  < |#eta| < %.1f "%(options.minEta,options.maxEta)))
if options.minEta == 0:
    latex4 = ROOT.TLatex(0.20,0.69,("|#eta| < %.1f "%(options.maxEta)))
latex4.SetNDC()
latex4.SetTextSize(mytextsize)



def makeMassVsPu(h2, rebin):
    h2.RebinX(rebin)
    pfx = h2.ProfileX().Clone('pfx')
    return pfx


def makeResponseVsPu(h2, hname, rebin):

    hmean = ROOT.TH1F(hname,hname,100/rebin,0,100)
    hrms  = ROOT.TH1F(hname+'2',hname+'2',100/rebin,0,100)
    hmeanfit   = ROOT.TH1F(hname+'3',hname+'3',100/rebin,0,100)
    hsigmafit  = ROOT.TH1F(hname+'4',hname+'4',100/rebin,0,100)

    fitfun = ROOT.TF1("fitfun","gaus")

    j = 0
    px = (h2.ProjectionX('px')).Clone('px')
    for bin in range(1,h2.GetNbinsX()+1,rebin): # rebin                                                                                                                                                                      
        firstbin = bin
        lastbin  = bin+(rebin-1)
        py       = (h2.ProjectionY('py',firstbin,lastbin)).Clone('py')
        mean     = py.GetMean()
        meanerr  = py.GetMeanError()
        rms      = py.GetRMS()
        rmserr   = py.GetRMSError()
        
        fitfun.SetRange(py.GetMean()-2*py.GetRMS(),py.GetMean()+2*py.GetRMS())
        py.Fit("fitfun","QRN")
        meanfit    = fitfun.GetParameter(1)
        meanfiterr = fitfun.GetParError(1)
        sigmafit     = fitfun.GetParameter(2)
        sigmafiterr  = fitfun.GetParError(2)


        if (py.GetEntries()>0):
            j = j + 1
            x  = (px.GetXaxis().GetBinUpEdge(lastbin) + px.GetXaxis().GetBinLowEdge(firstbin))/2
            ex = (px.GetXaxis().GetBinUpEdge(lastbin) - px.GetXaxis().GetBinLowEdge(firstbin))/2
            
            thisbin=hmean.FindBin(x)
            hmean.Fill(x,mean)
            hrms.Fill(x,rms)
            hmean.SetBinError(thisbin,meanerr)
            hrms.SetBinError(thisbin,rmserr)

            hmeanfit.Fill(x,meanfit)
            hsigmafit.Fill(x,sigmafit)
            hmeanfit.SetBinError(thisbin,meanfiterr)
            hsigmafit.SetBinError(thisbin,sigmafiterr)

    return hmean, hrms, hmeanfit, hsigmafit


    
if __name__ == '__main__':

    filename = options.input
    outdir   = options.outdir
    try:
        os.mkdir(outdir)
    except:
        print 'Cannot create output directory: directory already exists'
        #sys.exit()

    yaxisTitle = 'arbitrary units'


    algos = ['GEN', 'PF+PUPPI' , 'PF', 'PF+CHS', 'PF(Cleansing)', 'PF+CHS(Const.Sub.)']
    #algos = ['GEN', 'PF+PUPPI' , 'PF', 'PF+CHS', 'PF+CHS(Const.Sub.)']

    histos  = {'GEN' : ['gen/hm_leadjet_gen'],
               'PF+PUPPI' : ['puppi/hm_leadjet_puppi'],
               'PF'  : ['pf/hm_leadjet_pf'],
               'PF+CHS' : ['pfchs/hm_leadjet_pfchs'],
               'PF(Cleansing)' : ['pf/hmclean_leadjet_pf'],
               'PF+CHS(Const.Sub.)' : ['pfchs/hmconst_leadjet_pfchs'],
               }
     
    var = 'm'

    styles = {} # color, linestyle, line width
    styles['GEN'] = [ROOT.kBlack, 1, 2]
    styles['PF+PUPPI'] = [ROOT.kGreen+1, 1, 2]
    styles['PF'] = [ROOT.kBlue, 1, 2]
    styles['PF+CHS'] = [ROOT.kMagenta, 1, 2]
    styles['PF(Cleansing)'] = [ROOT.kOrange, 1, 2]
    styles['PF+CHS(Const.Sub.)'] = [ROOT.kCyan, 1, 2]

    
    # -- open file
    f = ROOT.TFile.Open(filename);

    # -- canvas
    #c = ROOT.TCanvas('%s_leadjet'%var,'%s_leadjet'%var,700,700)
    #cresponse = ROOT.TCanvas('%s_response_leadjet'%var,'%s_response_leadjet'%var,700,700)
    c = ROOT.TCanvas('%s_leadjet'%var,'%s_leadjet'%var,1000,800)
    cresponse = ROOT.TCanvas('%s_response_leadjet'%var,'%s_response_leadjet'%var,1000,800)

    # -- legend                                               
    leg1 = ROOT.TLegend(0.64,0.59,0.97,0.93);
    leg1.SetBorderSize(0);
    leg1.SetFillStyle(0);

    leg2 = ROOT.TLegend(0.68,0.30,0.97,0.91);
    leg2.SetBorderSize(0);
    leg2.SetFillStyle(0);
    leg2.SetTextSize(0.032)

    leg4 = ROOT.TLegend(0.64,0.61,0.97,0.93);
    leg4.SetBorderSize(0);
    leg4.SetFillStyle(0);

    # -- now plot
    
    # -- mass 1D distribution and mass resolution
    nre = 2
    nrer = 1
    i = 0
    j = 0

    func = ROOT.TF1('func','gaus',3)
    hmean   = ROOT.TH1F('hmean','hmean',5,0,5)
    hrms    = ROOT.TH1F('hrms','hrms',5,0,5)
    hsigma  = ROOT.TH1F('hsigma','hsigma',5,0,5)

    for algo in algos:

        h = f.Get(histos[algo][0])
        h.Rebin(nre)
        h.GetXaxis().SetTitle('m (GeV)')
        h.GetYaxis().SetTitle(yaxisTitle)
        h.GetYaxis().SetTitleOffset(1.6)
        h.SetLineColor(styles[algo][0])
        h.SetLineStyle(styles[algo][1])
        h.SetLineWidth(styles[algo][2])
        ymax = h.GetMaximum()

        leg1.AddEntry(h,algo,'L')

        if (i == 0):
            c.cd()
            if ('W' in options.sample):
                h.GetYaxis().SetRangeUser(0,ymax*1.4)
            else:
                h.GetYaxis().SetRangeUser(0,ymax*1.5)
            h.Draw()
        else:
            c.cd()
            h.Draw("same")

        i = i+1


        # response
        hr = f.Get(histos[algo][0].replace('_leadjet','_response_leadjet'))
        hr.Rebin(nrer)
        hr.GetXaxis().SetTitle('m_{reco} - m_{gen}(GeV)')
        hr.GetYaxis().SetTitle(yaxisTitle)
        hr.GetYaxis().SetTitleOffset(1.6)
        hr.SetLineColor(styles[algo][0])
        hr.SetLineStyle(styles[algo][1])
        hr.SetLineWidth(styles[algo][2])
        yrmax = hr.GetMaximum()
        
        if algo != 'GEN':
            
            func.SetRange(hr.GetMean()-2*hr.GetRMS(),hr.GetMean()+2*hr.GetRMS())
            hr.Fit("func","RNQ")
            mean    = hr.GetMean()
            meanerr = hr.GetMeanError()
            rms     = hr.GetRMS()
            rmserr  = hr.GetRMSError()
            sigma   = func.GetParameter(2)
            sigmaerr = func.GetParError(2)
            
            if (j==0): 
                print 'ALGO        <m-m_{gen}>(GeV)    RMS(GeV)    SIGMA(GeV)'
                print ('%s'%algo).ljust(18), ('    %.1f    %.1f    %.1f'%(mean, rms, sigma))

            hmean.Fill(algo,mean)
            hmean.SetBinError(j+1,meanerr)
            
            hrms.Fill(algo,rms)
            hrms.SetBinError(j+1,rmserr)
            
            hsigma.Fill(algo,sigma)
            hsigma.SetBinError(j+1,sigmaerr)
            
            #legentry = '#splitline{%s}{#splitline{<#Deltam>=%.1f GeV}{RMS=%.1f GeV}}'%(algo,mean,rms)
            #leg2.AddEntry(hr,legentry,'L')
            leg2.AddEntry(hr,algo,'L')
            leg2.AddEntry(0,'<#Deltam>=%.1f GeV'%mean,'')
            leg2.AddEntry(0,'RMS=%.1f GeV'%rms,'')
            leg2.AddEntry(0,'','')

            leg4.AddEntry(hr,algo,'L')

            if (j == 0):
                cresponse.cd()
                hr.GetYaxis().SetRangeUser(0, yrmax*1.4)
                hr.Draw("")
            else:
                cresponse.cd()
                hr.Draw("same")
                
            j = j+1

            
    # mass, mass response, mass resolution vs npu
    rebin = 5
    hmvspu = {}
    hmeanvspu = {}
    hrmsvspu = {}
    for algo in algos:
        h2 = f.Get(histos[algo][0].replace('_leadjet','_vs_npu_leadjet'))
        hmvspu[algo]=makeMassVsPu(h2,rebin)

        h2r =  f.Get(histos[algo][0].replace('_leadjet','_response_vs_npu_leadjet'))
        h,hh, hfit, hhfit = makeResponseVsPu(h2r,algo,rebin)
        hmeanvspu[algo]=h 
        hrmsvspu[algo]=hh

    cmvspu = ROOT.TCanvas('%s_vs_npu'%var,'%s_vs_npu'%var,1000,800)    
    cmeanvspu = ROOT.TCanvas('%s_response_vs_npu'%var,'%s_bias_vs_npu'%var,1000,800)    
    crmsvspu = ROOT.TCanvas('%s_resolution_vs_npu'%var,'%s_resolution_vs_npu'%var,1000,800)    

    for h in hmvspu, hmeanvspu, hrmsvspu:
        for algo in h:
            #h[algo].GetXaxis().SetTitle('n_{PU}')
            h[algo].GetXaxis().SetTitle('n_{PV}')
            #
            h[algo].GetYaxis().SetTitleOffset(1.3)
            h[algo].SetLineColor(styles[algo][0])
            h[algo].SetLineStyle(styles[algo][1])
            h[algo].SetLineWidth(styles[algo][2])
            h[algo].SetMarkerColor(styles[algo][0])
            h[algo].SetMarkerStyle(20)
            h[algo].SetMarkerSize(1.5)

    i = 0
    for algo in algos:
        if (i==0):
            cmvspu.cd()
            hmvspu[algo].GetYaxis().SetTitle('<m> (GeV)')
            hmvspu[algo].GetXaxis().SetRangeUser(15,50)
            if ("W" in options.sample):
                hmvspu[algo].GetYaxis().SetRangeUser(30,150)
            else:
                hmvspu[algo].GetYaxis().SetRangeUser(0,60)
            hmvspu[algo].Draw("l")

            cmeanvspu.cd()
            hmeanvspu[algo].GetYaxis().SetTitle('<m_{reco}-m_{gen}> (GeV)')
            hmeanvspu[algo].GetXaxis().SetRangeUser(15,50)
            hmeanvspu[algo].GetYaxis().SetRangeUser(-10,50)
            hmeanvspu[algo].Draw()

            crmsvspu.cd()
            hrmsvspu[algo].GetYaxis().SetTitle('RMS(m_{reco}-m_{gen}) (GeV)')
            hrmsvspu[algo].GetXaxis().SetRangeUser(15,50)
            hrmsvspu[algo].GetYaxis().SetRangeUser(0,35)
            hrmsvspu[algo].Draw()
        else:
            cmvspu.cd()
            hmvspu[algo].Draw('same')
            cmeanvspu.cd()
            hmeanvspu[algo].Draw('same')
            crmsvspu.cd()
            hrmsvspu[algo].Draw('same')
        i = i + 1


    
    # add text
    for canvas in c, cmvspu:
        canvas.cd()
        #cmsprel.Draw()
        latex0.Draw()
        latex1.Draw()
        latex2.Draw()
        latex3.Draw()
        latex4.Draw()
        leg1.Draw()
        CMS_lumi.CMS_lumi(canvas, 4, 0)

    for canvas in cresponse, cmeanvspu, crmsvspu:
        canvas.cd()
        #cmsprel.Draw()
        latex0.Draw()
        latex1.Draw()
        latex2.Draw()
        latex3.Draw()
        latex4.Draw()
        CMS_lumi.CMS_lumi(canvas, 4, 0)
    
    for canvas in cmeanvspu, crmsvspu: 
        canvas.cd()
        leg4.Draw()
        CMS_lumi.CMS_lumi(canvas, 4, 0)

    cresponse.cd()
    #leg2.SetTextAlign(12)
    leg2.Draw()
    CMS_lumi.CMS_lumi(cresponse, 4, 0)
        

    # mass vs algo
    ROOT.gStyle.SetPadRightMargin(0.15);
    ROOT.gStyle.SetErrorX(0.5);
    cmean=ROOT.TCanvas('%s_response_summary_leadjet'%var,'%s_response_summary_leadjet'%var, 1200,800)
    hmean.GetYaxis().SetRangeUser(-10,40)
    hmean.GetYaxis().SetTitle('<m_{reco} - m_{gen}> (GeV)')
    hmean.SetLineColor(ROOT.kRed)
    hmean.SetMarkerColor(ROOT.kRed)
    hmean.SetMarkerStyle(20)
    hmean.Draw('e')
    #cmsprel.Draw()
    latex0.Draw()
    latex1.Draw()
    latex2.Draw()
    latex3.Draw()
    latex4.Draw()
    CMS_lumi.CMS_lumi(cmean, 4, 0)

    # resolution vs algo
    cresol=ROOT.TCanvas('%s_resolution_summary_leadjet'%var,'%s_resolution_summary_leadjet'%var, 1200,800)
    #ROOT.gStyle.SetPadRightMargin(0.08);
    #ROOT.gStyle.SetErrorX(0.5);
    hrms.GetYaxis().SetRangeUser(0,35)
    hrms.GetYaxis().SetTitle('mass resolution (GeV)')
    hrms.SetLineColor(ROOT.kRed)
    hrms.SetMarkerColor(ROOT.kRed)
    hrms.SetMarkerStyle(20)
    hrms.Draw('e')
    hsigma.SetLineColor(ROOT.kRed)
    hsigma.SetMarkerColor(ROOT.kRed)
    hsigma.SetMarkerStyle(24)
    hsigma.Draw('esame')
    leg3 = ROOT.TLegend(0.67,0.75,0.90,0.9);
    leg3.SetBorderSize(0);
    leg3.SetFillStyle(0);
    leg3.AddEntry(hrms,'RMS','PL')
    leg3.AddEntry(hsigma,'fitted #sigma','PL')
    leg3.Draw()
    #cmsprel.Draw()
    latex0.Draw()
    latex1.Draw()
    latex2.Draw()
    latex3.Draw()
    latex4.Draw()
    CMS_lumi.CMS_lumi(cresol, 4, 0)

    raw_input('ok?')




    for typ in '.png','.pdf','.root':
        c.SaveAs(outdir+"/"+c.GetName()+typ)
        cresponse.SaveAs(outdir+"/"+cresponse.GetName()+typ)
        cmean.SaveAs(outdir+"/"+cmean.GetName() +typ)
        cresol.SaveAs(outdir+"/"+cresol.GetName() +typ)
        cmvspu.SaveAs(outdir+"/"+cmvspu.GetName()+typ)
        cmeanvspu.SaveAs(outdir+"/"+cmeanvspu.GetName()+typ)
        crmsvspu.SaveAs(outdir+"/"+crmsvspu.GetName()+typ)

