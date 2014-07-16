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

############################################
#            Job steering                  #
############################################
from optparse import OptionParser

parser = OptionParser()

parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
parser.add_option('-i','--input',action="store",type="string",dest="input",default="outtre.root")
parser.add_option('-o','--outdir',action="store",type="string",dest="outdir",default="plots")
parser.add_option('--nPU',action="store",type="int",dest="nPU",default=40)
parser.add_option('-r',action="store",type="float",dest="radius",default=0.8)
parser.add_option('--minPt',action="store",type="float",dest="minPt",default=25.)
parser.add_option('--maxPt',action="store",type="float",dest="maxPt",default=300.)
parser.add_option('--minEta',action="store",type="float",dest="minEta",default=0.)
parser.add_option('--maxEta',action="store",type="float",dest="maxEta",default=2.5)

(options, args) = parser.parse_args()



# cms preliminary 
cmsprel = ROOT.TLatex(0.20,0.96,("CMS Simulation Preliminary, #sqrt{s} = 13 TeV"))
cmsprel.SetNDC()
cmsprel.SetTextSize(0.03)


# text
latex1 = ROOT.TLatex(0.20,0.90,("Anti-kT (R=%.1f)"%(options.radius)))
latex1.SetNDC()
latex1.SetTextSize(0.03)
latex2 = ROOT.TLatex(0.20,0.85,("<n_{PU}> = "+str(options.nPU)))
latex2.SetNDC()
latex2.SetTextSize(0.03)
latex3 = ROOT.TLatex(0.20,0.80,("%.0f GeV < p_{T} < %.0f GeV "%(options.minPt,options.maxPt)))
latex3.SetNDC()
latex3.SetTextSize(0.03)
latex4 = ROOT.TLatex(0.20,0.75,("%.1f  < |#eta| < %.1f "%(options.minEta,options.maxEta)))
if options.minEta == 0:
    latex4 = ROOT.TLatex(0.20,0.75,("|#eta| < %.1f "%(options.maxEta)))
latex4.SetNDC()
latex4.SetTextSize(0.03)

    
if __name__ == '__main__':

    filename = options.input
    outdir   = options.outdir
    try:
        os.mkdir(outdir)
    except:
        print 'Cannot create output directory: directory already exists'
        sys.exit()

    docmssw = False
    #docmssw = True
    



    algos = ['GEN', 'PUPPI' , 'PF', 'PFCHS', 'PF(Cleansing)', 'PFCHS(Const.Sub.)']
        
    histos  = {'GEN' : ['gen/hm_leadjet_gen'],
               'PUPPI' : ['puppi/hm_leadjet_puppi'],
               'PF'  : ['pf/hm_leadjet_pf'],
               'PFCHS' : ['pfchs/hm_leadjet_pfchs'],
               'PF(Cleansing)' : ['pf/hmclean_leadjet_pf'],
               'PFCHS(Const.Sub.)' : ['pfchs/hmconst_leadjet_pfchs']
               }

     
    styles = {} # color, linestyle, line width
    styles['GEN'] = [ROOT.kBlack, 1, 2]
    styles['PUPPI'] = [ROOT.kGreen+1, 1, 2]
    styles['PF'] = [ROOT.kBlue, 1, 2]
    styles['PFCHS'] = [ROOT.kMagenta, 1, 2]
    styles['PF(Cleansing)'] = [ROOT.kOrange, 1, 2]
    styles['PFCHS(Const.Sub.)'] = [ROOT.kCyan, 1, 2]
    
    # -- open file
    f = ROOT.TFile.Open(filename);

    # -- canvas
    c = ROOT.TCanvas('mass_leadjet','mass_leadjet',500,500)
    cresponse = ROOT.TCanvas('mass_response_leadjet','mass_response_leadjet',500,500)

    # -- legend                                               
    leg1 = ROOT.TLegend(0.65,0.6,0.98,0.9);
    leg1.SetBorderSize(0);
    leg1.SetFillStyle(0);

    leg2 = ROOT.TLegend(0.65,0.6,0.98,0.9);
    leg2.SetBorderSize(0);
    leg2.SetFillStyle(0);

    # -- now plot
    
    # mass distribution
    nre = 2
    i = 0
    j = 0
    for algo in algos:

        h = f.Get(histos[algo][0])
        h.Rebin(nre)
        h.GetXaxis().SetTitle('mass (GeV)')
        h.GetYaxis().SetTitle('events')
        h.GetYaxis().SetTitleOffset(1.6)
        h.SetLineColor(styles[algo][0])
        h.SetLineStyle(styles[algo][1])
        h.SetLineWidth(styles[algo][2])
        ymax = h.GetMaximum()

        leg1.AddEntry(h,algo,'L')

        if (i == 0):
            c.cd()
            h.GetYaxis().SetRangeUser(0,ymax*1.4)
            h.Draw()
        else:
            c.cd()
            h.Draw("same")

        i = i+1


        hresponse = f.Get(histos[algo][0].replace('_leadjet','_response_leadjet'))
        hresponse.Rebin(nre)
        hresponse.GetXaxis().SetTitle('mass - gen mass(GeV)')
        hresponse.GetYaxis().SetTitle('events')
        hresponse.GetYaxis().SetTitleOffset(1.6)
        hresponse.SetLineColor(styles[algo][0])
        hresponse.SetLineStyle(styles[algo][1])
        hresponse.SetLineWidth(styles[algo][2])
        yrmax = hresponse.GetMaximum()

        if algo != 'GEN':

            leg2.AddEntry(hresponse,algo,'L')

            if (j == 0):
                cresponse.cd()
                hresponse.GetYaxis().SetRangeUser(0, yrmax*1.3)
                hresponse.Draw()
            else:
                cresponse.cd()
                hresponse.Draw("same")

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

