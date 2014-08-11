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
parser.add_option('--groomer',action="store",type="string",dest="groomer",default="trim")

(options, args) = parser.parse_args()



# cms preliminary 
cmsprel = ROOT.TLatex(0.20,0.96,("CMS Simulation Preliminary, #sqrt{s} = 13 TeV"))
cmsprel.SetNDC()
cmsprel.SetTextSize(0.03)

# text
mytextsize = 0.037
latex1 = ROOT.TLatex(0.20,0.89,("%s jets, Anti-kT (R=%.1f)"%(options.sample, options.radius)))
latex1.SetNDC()
latex1.SetTextSize(mytextsize)
latex2 = ROOT.TLatex(0.20,0.84,("<n_{PU}> = "+str(options.nPU)))
latex2.SetNDC()
latex2.SetTextSize(mytextsize)
latex3 = ROOT.TLatex(0.20,0.79,("%.0f GeV < p_{T} < %.0f GeV "%(options.minPt,options.maxPt)))
latex3.SetNDC()
latex3.SetTextSize(mytextsize)
latex4 = ROOT.TLatex(0.20,0.74,("%.1f  < |#eta| < %.1f "%(options.minEta,options.maxEta)))
if options.minEta == 0:
    latex4 = ROOT.TLatex(0.20,0.74,("|#eta| < %.1f "%(options.maxEta)))
latex4.SetNDC()
latex4.SetTextSize(mytextsize)



def makeMassVsPu(h2, rebin):
    h2.RebinX(rebin)
    pfx = h2.ProfileX().Clone('pfx')
    return pfx


def makeResponseVsPu(h2, hname, rebin):

    hmean = ROOT.TH1F(hname,hname,100/rebin,0,100)
    hrms  = ROOT.TH1F(hname+'2',hname+'2',100/rebin,0,100)

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

        if (py.GetEntries()>0):
            j = j + 1
            x  = (px.GetXaxis().GetBinUpEdge(lastbin) + px.GetXaxis().GetBinLowEdge(firstbin))/2
            ex = (px.GetXaxis().GetBinUpEdge(lastbin) - px.GetXaxis().GetBinLowEdge(firstbin))/2
            
            thisbin=hmean.FindBin(x)
            hmean.Fill(x,mean)
            hrms.Fill(x,rms)
            hmean.SetBinError(thisbin,meanerr)
            hrms.SetBinError(thisbin,rmserr)

    return hmean, hrms


    
if __name__ == '__main__':

    filename = options.input
    outdir   = options.outdir
    try:
        os.mkdir(outdir)
    except:
        print 'Cannot create output directory: directory already exists'
        #sys.exit()

    
    var      = 'm'+options.groomer
    varsafe  = var+'safe'
    varlabel = 'm_{%s}'%options.groomer

    algos = ['GEN', 'PF+PUPPI' , 'PF', 'PF+CHS']

    histos  = {'GEN' : ['gen/h%s_leadjet_gen'%var],
               'PF+PUPPI' : ['puppi/h%s_leadjet_puppi'%varsafe],
               'PF'  : ['pf/h%s_leadjet_pf'%varsafe],
               'PF+CHS' : ['pfchs/h%s_leadjet_pfchs'%varsafe],
               }
     
    styles = {} # color, linestyle, line width
    styles['GEN'] = [ROOT.kBlack, 1, 2]
    styles['PF+PUPPI'] = [ROOT.kGreen+1, 1, 2]
    styles['PF'] = [ROOT.kBlue, 1, 2]
    styles['PF+CHS'] = [ROOT.kMagenta, 1, 2]
    
    # -- open file
    f = ROOT.TFile.Open(filename);

    # -- canvas
    #c = ROOT.TCanvas('%s_leadjet'%var,'%s_leadjet'%var,700,700)
    #cresponse = ROOT.TCanvas('%s_response_leadjet'%var,'%s_response_leadjet'%var,700,700)
    c = ROOT.TCanvas('%s_leadjet'%varsafe,'%s_leadjet'%varsafe,1000,800)
    cresponse = ROOT.TCanvas('%s_response_leadjet'%varsafe,'%s_response_leadjetsafe'%var,1000,800)

    # -- legend                                               
    leg1 = ROOT.TLegend(0.7,0.7,0.97,0.92);
    leg1.SetBorderSize(0);
    leg1.SetFillStyle(0);

    leg2 = ROOT.TLegend(0.60,0.68,0.95,0.95);
    leg2.SetBorderSize(0);
    leg2.SetFillStyle(0);

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

        print histos[algo][0]
        h = f.Get(histos[algo][0])
        h.Rebin(nre)
        h.GetXaxis().SetTitle('m (GeV)')
        h.GetYaxis().SetTitle('events')
        h.GetYaxis().SetTitleOffset(1.6)
        h.SetLineColor(styles[algo][0])
        h.SetLineStyle(styles[algo][1])
        h.SetLineWidth(styles[algo][2])
        ymax = h.GetMaximum()

        leg1.AddEntry(h,algo,'L')

        if (i == 0):
            c.cd()
            h.GetYaxis().SetRangeUser(0,ymax*1.7)
            h.Draw()
        else:
            c.cd()
            h.Draw("same")

        i = i+1


        hr = f.Get(histos[algo][0].replace('_leadjet','_response_leadjet'))
        hr.Rebin(nrer)
        hr.GetXaxis().SetTitle('m - m_{gen}(GeV)')
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
                print 'ALGO        <m-m_{gen}>(GeV)    RMS(GeV)    SIGMA(GeV)'
            print ('%s'%algo).ljust(18), ('    %.1f    %.1f    %.1f'%(mean, rms, sigma))

            hmean.Fill(algo,mean)
            hmean.SetBinError(j+1,meanerr)

            hrms.Fill(algo,rms)
            hrms.SetBinError(j+1,rmserr)

            hsigma.Fill(algo,sigma)
            hsigma.SetBinError(j+1,sigmaerr)

            legentry = '#splitline{%s}{<#Deltam>=%.1f GeV, RMS=%.1f GeV}'%(algo,mean,rms)
            leg2.AddEntry(hr,legentry,'L')
            
            if (j == 0):
                cresponse.cd()
                hr.GetYaxis().SetRangeUser(0, yrmax*1.5)
                hr.Draw()
            else:
                cresponse.cd()
                hr.Draw("same")

            j = j+1

                
    # add text
    for canvas in c, cresponse:
        canvas.cd()
        CMS_lumi.CMS_lumi(canvas, 4, 0)
        #        cmsprel.Draw()
        latex1.Draw()
        latex2.Draw()
        latex3.Draw()
        latex4.Draw()

    c.cd()
    leg1.Draw()

    cresponse.cd()
    leg2.Draw()

    raw_input('ok?')




    for typ in '.png','.pdf','.root':
        c.SaveAs(outdir+"/"+c.GetName()+typ)
        cresponse.SaveAs(outdir+"/"+cresponse.GetName()+typ)

