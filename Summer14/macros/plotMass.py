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
ROOT.gStyle.SetPadLeftMargin(0.16);

############################################
#            Job steering                  #
############################################
from optparse import OptionParser

parser = OptionParser()

parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
parser.add_option('-i','--input',action="store",type="string",dest="input",default="outtre.root")
parser.add_option('-o','--outdir',action="store",type="string",dest="outdir",default="plots")
parser.add_option('--nPU',action="store",type="int",dest="nPU",default=22)
parser.add_option('-r',action="store",type="float",dest="radius",default=0.5)
parser.add_option('--minPt',action="store",type="float",dest="minPt",default=25.)

(options, args) = parser.parse_args()

############################################################
def makeKinComparisonPlots(f, hname, types, plotAttributes, styles, outdir):

    # canvas
    c = ROOT.TCanvas(plotAttributes[0],plotAttributes[0],700,700);

    # legend
    leg = ROOT.TLegend(0.7,0.7,0.93,0.9);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);

    h = {}
    ymax = -1
    for typ, suff in types.iteritems():
        h[typ] = f.Get(suff+'/'+hname+'_'+suff)
        tmp = h[typ].GetMaximum()
        if tmp > ymax:
            ymax = tmp

    n = 0        
    for typ in types:
        if (h[typ]):
            leg.AddEntry(h[typ], typ, "l")
            h[typ].SetXTitle(plotAttributes[1])
            h[typ].SetYTitle(plotAttributes[2])
            h[typ].SetLineColor(styles[typ][0])
            h[typ].SetLineStyle(styles[typ][1])
            h[typ].SetLineWidth(styles[typ][2])
            h[typ].GetYaxis().SetTitleOffset(1.3)
            c.cd()
            if (n==0):
                h[typ].GetYaxis().SetRangeUser(0,ymax*1.3)
                h[typ].Draw()
            else:
                h[typ].Draw('same')
            n = n + 1    

    # text
    latex1 = ROOT.TLatex(0.20,0.89,("Anti-kT (R=%.1f)"%(options.radius)))
    latex1.SetNDC()
    latex1.SetTextSize(0.03)
    latex2 = ROOT.TLatex(0.20,0.84,("n_{PU} = "+str(options.nPU)))
    latex2.SetNDC()
    latex2.SetTextSize(0.03)
    latex3 = ROOT.TLatex(0.20,0.79,("p_{T} > %.0f GeV "%options.minPt))
    latex3.SetNDC()
    latex3.SetTextSize(0.03)
    
    latex1.Draw()
    latex2.Draw()
    latex3.Draw()
    leg.Draw()
    
    c.SaveAs(outdir+"/"+c.GetName()+".png");
    c.SaveAs(outdir+"/"+c.GetName()+".pdf");


def makeResponseComparisonPlots(f, hname, types, plotAttributes, styles, outdir):

    # canvas
    c = ROOT.TCanvas(plotAttributes[0],plotAttributes[0],700,700);

    # legend
    leg = ROOT.TLegend(0.7,0.7,0.93,0.9);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);

    h = {}
    ymax = -1
    for typ, suff in types.iteritems():
        h[typ] = f.Get(suff+'/'+hname+'_'+suff)
        tmp = h[typ].GetMaximum()
        if (tmp > ymax and typ !='GEN'):
            ymax = tmp
            
    n = 0        
    for typ in types:
        if (h[typ] and typ != 'GEN' ):
            leg.AddEntry(h[typ], typ, "l")
            h[typ].SetXTitle(plotAttributes[1])
            h[typ].SetYTitle(plotAttributes[2])
            h[typ].SetLineColor(styles[typ][0])
            h[typ].SetLineStyle(styles[typ][1])
            h[typ].SetLineWidth(styles[typ][2])
            h[typ].GetYaxis().SetTitleOffset(1.3)
            c.cd()
            if (n==0):
                h[typ].GetYaxis().SetRangeUser(0,ymax*1.3)
                h[typ].Draw()
            else:
                h[typ].Draw('same')
        n = n + 1    


    
    # text
    latex1 = ROOT.TLatex(0.20,0.89,("Anti-kT (R=%.1f)"%(options.radius)))
    latex1.SetNDC()
    latex1.SetTextSize(0.03)
    latex2 = ROOT.TLatex(0.20,0.84,("n_{PU} = "+str(options.nPU)))
    latex2.SetNDC()
    latex2.SetTextSize(0.03)
    latex3 = ROOT.TLatex(0.20,0.79,("p_{T} > %.0f GeV "%options.minPt))
    latex3.SetNDC()
    latex3.SetTextSize(0.03)
    
    latex1.Draw()
    latex2.Draw()
    latex3.Draw()
    leg.Draw()
    
    c.SaveAs(outdir+"/"+c.GetName()+".png");
    c.SaveAs(outdir+"/"+c.GetName()+".pdf");



def makeResponseVsNpu(f, types, styles, outdir):

    for var in 'm', 'mraw':
        c1 = ROOT.TCanvas(var+'_responsemean_vs_npu',var+'_responsemean_vs_npu',700,700);
        c2 = ROOT.TCanvas(var+'_responserms_vs_npu',var+'_responserms_vs_npu',700,700);
        h = []
        graphmean = []
        graphrms = []

        n = 0
        for typ, suff in types.iteritems():

            h.append( f.Get(suff+'/h'+var+'_response_vs_npu_'+suff))
            
            graphmean.append(ROOT.TGraphErrors())
            graphrms.append(ROOT.TGraphErrors())
            graphmean[n].SetLineColor(styles[typ][0])
            graphrms[n].SetLineColor(styles[typ][0])
            graphmean[n].SetLineWidth(2)
            graphrms[n].SetLineWidth(2)

            j = 0
            for bin in range(1,h[n].GetNbinsX(),5):
                px = (h[n].ProjectionY('px',bin,bin+4)).Clone('px')
                mean = px.GetMean()
                meanerr = px.GetMeanError()
                rms = px.GetRMS()
                rmserr = px.GetRMSError()
                if (px.GetEntries()>0):
                    j = j + 1 
                    graphmean[n].SetPoint(j,bin+2.5,mean)
                    graphmean[n].SetPointError(j,2.5,meanerr)
                    graphrms[n].SetPoint(j,bin+2.5, rms)
                    graphrms[n].SetPointError(j,2.5,rmserr)
            print graphmean[n].GetN()
            if (n == 0):
                c1.cd()
                graphmean[n].Draw("ap")
                c2.cd()
                graphrms[n].Draw("ap")

            else:
                c1.cd()
                graphmean[n].Draw("p*same")
                c2.cd()
                graphrms[n].Draw("p*same")
            n = n + 1
            
        # save plots                                                                                                                                                                                                                  
        for p in '.pdf', '.png':
            c1.SaveAs(outdir+'/'+c1.GetName()+p)
            c2.SaveAs(outdir+'/'+c2.GetName()+p)


if __name__ == '__main__':

    filename = options.input
    outdir   = options.outdir
    try:
        os.mkdir(outdir)
    except:
        print 'Cannot create output directory: directory already exists'
        sys.exit()

    docmssw = False
    
    types = {'GEN':'gen','PUPPI':'puppi','PFlow':'pf','PFlowCHS':'pfchs'}
    if (docmssw):
        types = {'GEN':'gen','PUPPI':'puppi','PFlow':'pf','PFlowCHS':'pfchs','PF-CMSSW':'pfcmssw'}

    histograms = {'hmraw'          : ['mraw','raw mass (GeV)','events',1],
                  'hm'             : ['m','mass (GeV)','events',1],
                  'hmtrim'         : ['mtrim','trimmed mass (GeV)','events',1],
                  'hmtrimsafe'     : ['mtrimsafe','trimmed mass (GeV)','events',1],
                  'hmclean'        : ['mclean','clean mass (GeV)','events',1],
                  'hmconst'        : ['mconst','const subtracted mass (GeV)','events',1],

                  'hmraw_response'  : ['mraw_response','raw mass - gen mass(GeV)','events',1],
                  'hm_response'     : ['m_response','mass - gen mass(GeV)','events',1],
                  'hmtrim_response' : ['mtrim_response','trimmed mass - gen mass(GeV)','events',1],
                  'hmtrimsafe_response' : ['mtrimsafe_response','trimmed mass - gen mass(GeV)','events',1],
                  'hmclean_response' : ['mclean_response','cleansed mass - gen mass(GeV)','events',1],
                  'hmconst_response' : ['mconst_response','const subtracted mass - gen mass(GeV)','events',1],
                  }


    styles = {} # color, linestyle, line width
    styles['GEN'] = [ROOT.kBlack, 1, 2]
    styles['PUPPI'] = [ROOT.kGreen+1, 1, 2]
    styles['PFlow'] = [ROOT.kBlue, 1, 2]
    styles['PFlowCHS'] = [ROOT.kMagenta, 1, 2]
    styles['PF-CMSSW'] = [ROOT.kOrange, 1, 2]

    
    # -- open file
    f = ROOT.TFile.Open(filename);

    # -- rebin histograms
    for typ, suff in types.iteritems():
        for hname,plotAttributes in histograms.iteritems():
            h = f.Get(suff+'/'+hname+'_'+suff)
            h.Rebin(plotAttributes[3])

    # -- make plots 

    # kin. distributions        
    #for histogram, plotAttributes in histograms.iteritems():
    #    if ('pu' not in histogram and 'good' not in histogram and 'response' not in histogram):
    #        makeKinComparisonPlots(f, histogram, types, plotAttributes, styles, options.outdir)
    
    # response plots
    #for histogram, plotAttributes in histograms.iteritems():
    #    if ('response' in histogram):
    #        makeResponseComparisonPlots(f, histogram, types, plotAttributes, styles, options.outdir)

    

    makeResponseVsNpu(f, types, styles, options.outdir)
   


    raw_input('ok?')

        
        
