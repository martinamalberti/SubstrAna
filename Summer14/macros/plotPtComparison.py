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
parser.add_option('--all',action="store",type="int",dest="all",default=1)

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

    
if __name__ == '__main__':

    filename = options.input
    outdir   = options.outdir
    try:
        os.mkdir(outdir)
    except:
        print 'Cannot create output directory: directory already exists'
        #sys.exit()

    yaxisTitle = 'arbitrary units'

    algos = ['GEN', 'PF+PUPPI' , 'PF', 'PF+CHS', 'PF+SK(a=0.3)']
    
    #histos  = {'GEN' : ['gen/hpt_leadjet_gen'],
    #           'PF+PUPPI' : ['puppi/hptcorr_leadjet_puppi'],
    #           'PF'  : ['pf/hptcorr_leadjet_pf'],
    #           'PF+CHS' : ['pfchs/hptcorr_leadjet_pfchs'],
    #           }
         
    #var = 'ptcorr'

    #histos  = {'GEN' : ['gen/hpt_gen'],
    #           'PF+PUPPI' : ['puppi/hptcorrphil_puppi'],
    #           'PF'  : ['pf/hptcorrphil_pf'],
    #           'PF+CHS' : ['pfchs/hptcorrphil_pfchs'],
    #           }

    #histos  = {'GEN' : ['gen/hpt_leadjet_gen'],
    #           'PF+PUPPI' : ['puppi/hptcorrphil_leadjet_puppi'],
    #           'PF'  : ['pf/hptcorrphil_leadjet_pf'],
    #           'PF+CHS' : ['pfchs/hptcorrphil_leadjet_pfchs'],
    #           }


    var = 'pt'

    #histos  = {'GEN' : ['gen/h%s_gen'%var],
    #           'PF+PUPPI' : ['puppi/h%s_puppi'%var],
    #           'PF'  : ['pf/h%s_pf'%var],
    #           'PF+CHS' : ['pfchs/h%s_pfchs'%var],
    #           'PF+SOFTKILLER' : ['softkiller/h%s_softkiller'%var]}

    histos  = {'GEN' : ['gen/hpt_gen'],
               'PF+PUPPI' : ['puppi/hptraw_puppi'],
               'PF'  : ['pf/hpt_pf'],
               'PF+CHS' : ['pfchs/hpt_pfchs'],
               'PF+SK(a=0.3)' : ['softkiller/hptraw_softkiller']}
         


    styles = {} # color, linestyle, line width, marker style
    styles['GEN'] = [ROOT.kBlack, 1, 2, 20]
    styles['PF+PUPPI'] = [ROOT.kGreen+1, 1, 2, 21]
    styles['PF'] = [ROOT.kBlue, 1, 2, 22]
    styles['PF+CHS'] = [ROOT.kMagenta, 1, 2, 24]
    styles['PF+SK(a=0.3)'] = [ROOT.kRed, 1, 2, 25]
    
    # -- open file
    f = ROOT.TFile.Open(filename);

    # -- canvas
    c = ROOT.TCanvas('%s'%var,'%s'%var,1000,800)
    cresponse = ROOT.TCanvas('%s_response'%var,'%s_response'%var,1000,800)

    # -- legend                                               
    leg1 = ROOT.TLegend(0.7,0.68,0.97,0.92);
    leg1.SetBorderSize(0);
    leg1.SetFillStyle(0);

    leg2 = ROOT.TLegend(0.70,0.50,0.93,0.92);
    leg2.SetBorderSize(0);
    leg2.SetFillStyle(0);
    leg2.SetTextSize(0.030);

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
        h.GetYaxis().SetTitle(yaxisTitle)
        h.GetYaxis().SetTitleOffset(1.6)
        h.SetLineColor(styles[algo][0])
        h.SetLineStyle(styles[algo][1])
        h.SetLineWidth(styles[algo][2])
        h.SetMarkerStyle(styles[algo][3])
        h.SetMarkerColor(styles[algo][0])
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


        #hr = f.Get(histos[algo][0].replace('_leadjet','_response_leadjet'))
        #if  (options.all):
        hr = f.Get(histos[algo][0].replace('_','_response_'))
        hr.Rebin(2*nrer)
        hr.GetXaxis().SetTitle('p^{T} - p^{T}_{gen}(GeV)')
        hr.GetYaxis().SetTitle(yaxisTitle)
        hr.GetXaxis().SetTitleOffset(1.2)
        hr.GetYaxis().SetTitleOffset(1.6)
        hr.SetLineColor(styles[algo][0])
        hr.SetLineStyle(styles[algo][1])
        hr.SetLineWidth(styles[algo][2])
        hr.SetMarkerStyle(styles[algo][3])
        hr.SetMarkerColor(styles[algo][0])
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

            #legentry = '#splitline{%s}{#splitline{<#Deltap^{T}>=%.1f GeV}{RMS=%.1f GeV}}'%(algo,mean,rms)
            #leg2.AddEntry(hr,legentry,'L')
            leg2.AddEntry(hr,algo,'PL')
            leg2.AddEntry(0,'<#Deltap^{T}>=%.1f GeV'%mean,'')
            leg2.AddEntry(0,'RMS=%.1f GeV'%rms,'')
            leg2.AddEntry(0,'','')

            if (j == 0):
                cresponse.cd()
                #hr.GetYaxis().SetRangeUser(0, yrmax*1.5)
                hr.DrawNormalized("hp")
                yrmax = hr.GetMaximum()
                hr.GetYaxis().SetRangeUser(0., yrmax*1.5)
                hr.DrawNormalized("hp")
            else:
                cresponse.cd()
                hr.DrawNormalized("hpsame")

            j = j+1

    c.cd()
    CMS_lumi.CMS_lumi(c, 4, 0)
    #cmsprel.Draw()
    latex0.Draw()
    latex1.Draw()
    latex2.Draw()
    latex3.Draw()
    latex4.Draw()
    leg1.Draw()
    
    cresponse.cd()
    CMS_lumi.CMS_lumi(cresponse, 4, 0)
    #cmsprel.Draw()
    latex0.Draw()
    latex1.Draw()
    latex2.Draw()
    latex3.Draw()
    latex4.Draw()
    leg2.Draw()
    
    raw_input('ok?')




    for typ in '.png','.pdf','.root':
        if (options.all):
            c.SaveAs((outdir+"/"+c.GetName()+typ).replace('_leadjet',''))
            cresponse.SaveAs((outdir+"/"+cresponse.GetName()+typ).replace('_leadjet',''))
        else:
            c.SaveAs(outdir+"/"+c.GetName()+typ)
            cresponse.SaveAs(outdir+"/"+cresponse.GetName()+typ)
    
