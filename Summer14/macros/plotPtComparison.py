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
ROOT.gStyle.SetPadRightMargin(0.03);
#ROOT.gStyle.SetPadRightMargin(0.16);
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

# text
latex1 = ROOT.TLatex(0.20,0.89,("%s jets, Anti-kT (R=%.1f)"%(options.sample,options.radius)))
latex1.SetNDC()
latex1.SetTextSize(0.03)
latex2 = ROOT.TLatex(0.20,0.84,("<n_{PU}> = "+str(options.nPU)))
latex2.SetNDC()
latex2.SetTextSize(0.03)
latex3 = ROOT.TLatex(0.20,0.79,("%.0f GeV < p_{T} < %.0f GeV "%(options.minPt,options.maxPt)))
latex3.SetNDC()
latex3.SetTextSize(0.03)
latex4 = ROOT.TLatex(0.20,0.74,("%.1f  < |#eta| < %.1f "%(options.minEta,options.maxEta)))
if options.minEta == 0:
    latex4 = ROOT.TLatex(0.20,0.74,("|#eta| < %.1f "%(options.maxEta)))
latex4.SetNDC()
latex4.SetTextSize(0.03)
    
if __name__ == '__main__':

    filename = options.input
    outdir   = options.outdir
    try:
        os.mkdir(outdir)
    except:
        print 'Cannot create output directory: directory already exists'
        #sys.exit()


    algos = ['GEN', 'PF+PUPPI' , 'PF', 'PF+CHS']

    histos  = {'GEN' : ['gen/hpt_leadjet_gen'],
               'PF+PUPPI' : ['puppi/hptcorr_leadjet_puppi'],
               'PF'  : ['pf/hptcorr_leadjet_pf'],
               'PF+CHS' : ['pfchs/hptcorr_leadjet_pfchs'],
               }
         
    var = 'ptcorr'

    styles = {} # color, linestyle, line width, marker style
    styles['GEN'] = [ROOT.kBlack, 1, 2, 20]
    styles['PF+PUPPI'] = [ROOT.kGreen+1, 1, 2, 21]
    styles['PF'] = [ROOT.kBlue, 1, 2, 22]
    styles['PF+CHS'] = [ROOT.kMagenta, 1, 2, 23]
    
    # -- open file
    f = ROOT.TFile.Open(filename);

    # -- canvas
    #c = ROOT.TCanvas('%s_leadjet'%var,'%s_leadjet'%var,700,700)
    #cresponse = ROOT.TCanvas('%s_response_leadjet'%var,'%s_response_leadjet'%var,700,700)
    c = ROOT.TCanvas('%s_leadjet'%var,'%s_leadjet'%var,1000,800)
    cresponse = ROOT.TCanvas('%s_response_leadjet'%var,'%s_response_leadjet'%var,1000,800)

    # -- legend                                               
    leg1 = ROOT.TLegend(0.7,0.7,0.97,0.92);
    leg1.SetBorderSize(0);
    leg1.SetFillStyle(0);

    leg2 = ROOT.TLegend(0.60,0.68,0.95,0.95);
    leg2.SetBorderSize(0);
    leg2.SetFillStyle(0);

    # -- now plot
    
    # -- pt 1D distribution and pt resolution
    nre = 5
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
        h.GetXaxis().SetRangeUser(options.minPt, options.maxPt)
        h.GetXaxis().SetTitle('p^{T} (GeV)')
        h.GetYaxis().SetTitle('events')
        h.GetYaxis().SetTitleOffset(1.6)
        h.SetLineColor(styles[algo][0])
        h.SetLineStyle(styles[algo][1])
        h.SetLineWidth(styles[algo][2])
        #h.SetMarkerStyle(styles[algo][3])
        #h.SetMarkerColor(styles[algo][0])
        ymax = h.GetMaximum()

        leg1.AddEntry(h,algo,'L')

        if (i == 0):
            c.cd()
            h.GetYaxis().SetRangeUser(0,ymax*1.5)
            h.Draw("hp")
        else:
            c.cd()
            h.Draw("hpsame")

        i = i+1


        hr = f.Get(histos[algo][0].replace('_leadjet','_response_leadjet'))
        hr.Rebin(nrer)
        hr.GetXaxis().SetTitle('p^{T} - p^{T}_{gen}(GeV)')
        hr.GetYaxis().SetTitle('events')
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
                print 'ALGO        <#DeltaPT>(GeV)    RMS(GeV)    SIGMA(GeV)'
            print ('%s'%algo).ljust(18), ('    %.1f    %.1f    %.1f'%(mean, rms, sigma))

            hmean.Fill(algo,mean)
            hmean.SetBinError(j+1,meanerr)

            hrms.Fill(algo,rms)
            hrms.SetBinError(j+1,rmserr)

            hsigma.Fill(algo,sigma)
            hsigma.SetBinError(j+1,sigmaerr)

            legentry = '#splitline{%s}{<#Deltap^{T}>=%.1f GeV, RMS=%.1f GeV}'%(algo,mean,rms)
            leg2.AddEntry(hr,legentry,'L')
            
            if (j == 0):
                cresponse.cd()
                hr.GetYaxis().SetRangeUser(0, yrmax*1.5)
                hr.Draw()
            else:
                cresponse.cd()
                hr.Draw("same")

            j = j+1

    c.cd()
    cmsprel.Draw()
    latex1.Draw()
    latex2.Draw()
    latex3.Draw()
    latex4.Draw()
    leg1.Draw()
    
    cresponse.cd()
    cmsprel.Draw()
    latex1.Draw()
    latex2.Draw()
    latex3.Draw()
    latex4.Draw()
    leg2.Draw()
    
    raw_input('ok?')




    for typ in '.png','.pdf','.root':
        c.SaveAs(outdir+"/"+c.GetName()+typ)
        cresponse.SaveAs(outdir+"/"+cresponse.GetName()+typ)
    
