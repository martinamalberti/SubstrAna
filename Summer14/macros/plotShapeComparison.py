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
parser.add_option('--minPt',action="store",type="float",dest="minPt",default=200.)
parser.add_option('--maxPt',action="store",type="float",dest="maxPt",default=600.)
parser.add_option('--minEta',action="store",type="float",dest="minEta",default=0.)
parser.add_option('--maxEta',action="store",type="float",dest="maxEta",default=2.5)

(options, args) = parser.parse_args()



# cms preliminary 
cmsprel = ROOT.TLatex(0.20,0.96,("CMS Simulation Preliminary, #sqrt{s} = 13 TeV"))
cmsprel.SetNDC()
cmsprel.SetTextSize(0.03)

# text
latex1 = ROOT.TLatex(0.20,0.89,("Anti-kT (R=%.1f)"%(options.radius)))
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



def makeShapeVsPu(h2, rebin):
    h2.RebinX(rebin)
    pfx = h2.ProfileX().Clone('pfx')
    return pfx
    
if __name__ == '__main__':

    filename = options.input
    outdir   = options.outdir
    try:
        os.mkdir(outdir)
    except:
        print 'Cannot create output directory: directory already exists'
        #sys.exit()

    docmssw = False
    #docmssw = True
    

    #algos = ['GEN', 'PUPPI' , 'PF', 'PFCHS', 'PUPPI(soft-drop,#beta=2)', 'PFCHS(soft-drop,#beta=2)','PFCHS(soft-drop,#beta=2)']
    algos = ['GEN', 'PF+PUPPI' , 'PF', 'PF+CHS']
        
    histos  = {'GEN' : ['gen/htau21_leadjet_gen'],
               'PF+PUPPI' : ['puppi/htau21_leadjet_puppi'],
               'PF'  : ['pf/htau21_leadjet_pf'],
               'PF+CHS' : ['pfchs/htau21_leadjet_pfchs'],
               #'PUPPI(soft-drop,#beta=2)' : ['puppi/htau21_softdrop_leadjet_puppi'],
               #'PF(soft-drop,#beta=2)' : ['pf/htau21_softdrop_leadjet_pf'],
               #'PFCHS(soft-drop,#beta=2)' : ['pfchs/htau21_softdrop_leadjet_pfchs'],
               }

     
    styles = {} # color, linestyle, line width
    styles['GEN'] = [ROOT.kBlack, 1, 2]
    styles['PF+PUPPI'] = [ROOT.kGreen+1, 1, 2]
    styles['PF'] = [ROOT.kBlue, 1, 2]
    styles['PF+CHS'] = [ROOT.kMagenta, 1, 2]
    styles['PF+PUPPI(soft-drop,#beta=2)'] = [ROOT.kGreen+2, 1, 2]
    styles['PF(soft-drop,#beta=2)'] = [ROOT.kBlue+2, 1, 2]
    styles['PF+CHS(soft-drop,#beta=2)'] = [ROOT.kMagenta+2, 1, 2]


    # -- open file
    f = ROOT.TFile.Open(filename);

    # -- canvas
    c = ROOT.TCanvas('tau21_leadjet','tau21_leadjet',500,500)

    # -- legend                                               
    leg1 = ROOT.TLegend(0.7,0.7,0.97,0.92);
    leg1.SetBorderSize(0);
    leg1.SetFillStyle(0);

    # -- now plot
    
    # -- tau21 1D distribution and mass resolution

    nre = 10
    i = 0

    for algo in algos:

        h = f.Get(histos[algo][0])
        h.Rebin(nre)
        h.GetXaxis().SetTitle('#tau_{2}/#tau_{1}')
        h.GetYaxis().SetTitle('events')
        h.GetYaxis().SetTitleOffset(1.5)
        h.SetLineColor(styles[algo][0])
        h.SetLineStyle(styles[algo][1])
        h.SetLineWidth(styles[algo][2])
        ymax = h.GetMaximum()

        leg1.AddEntry(h,algo,'L')

        if (i == 0):
            c.cd()
            #h.GetYaxis().SetRangeUser(0,ymax*1.5)
            h.GetYaxis().SetRangeUser(0,ymax*3.5)
            h.Draw()
        else:
            c.cd()
            h.Draw("same")

        i = i+1

    # tau21 vs npu
    hvspu = {}
    for algo in algos:
        h2 = f.Get(histos[algo][0].replace('_leadjet','_vs_npu_leadjet'))
        hvspu[algo]=makeShapeVsPu(h2,5)

    cvspu = ROOT.TCanvas('tau21_vs_npu','tau21_vs_npu',1000,800)    
    

    i = 0
    for algo in algos:
        hvspu[algo].GetXaxis().SetTitle('N_{PU}')
        hvspu[algo].GetYaxis().SetTitleOffset(1.6)
        hvspu[algo].SetLineColor(styles[algo][0])
        hvspu[algo].SetLineStyle(styles[algo][1])
        hvspu[algo].SetLineWidth(styles[algo][2])
        hvspu[algo].SetMarkerColor(styles[algo][0])
        hvspu[algo].SetMarkerStyle(20)
        if (i==0):
            cvspu.cd()
            ROOT.gStyle.SetErrorX(0.5)
            hvspu[algo].GetYaxis().SetTitle('<#tau_{2}/#tau_{1}>')
            hvspu[algo].GetXaxis().SetRangeUser(25,55)
            hvspu[algo].GetYaxis().SetRangeUser(0,1.3)
            hvspu[algo].Draw("")
        else:
            cvspu.cd()
            hvspu[algo].Draw('same')
        i = i + 1


    
    # add text
    for canvas in c, cvspu:
        canvas.cd()
        cmsprel.Draw()
        latex1.Draw()
        latex2.Draw()
        latex3.Draw()
        latex4.Draw()
        leg1.Draw()


    raw_input('ok?')




    for typ in '.png','.pdf','.root':
        c.SaveAs(outdir+"/"+c.GetName()+typ)
        cvspu.SaveAs(outdir+"/tau21_vs_npu_leadjet"+typ)
        

