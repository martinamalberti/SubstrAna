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

ROOT.gStyle.SetLegendFont(42)
ROOT.gStyle.SetTextFont(42)
ROOT.gStyle.SetMarkerSize(1.5)

############################################
#            Job steering                  #
############################################
from optparse import OptionParser

parser = OptionParser()

parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
parser.add_option('-i','--input',action="store",type="string",dest="input",default="histograms.root")
parser.add_option('-o','--outdir',action="store",type="string",dest="outdir",default="plots")
parser.add_option('--nPU',action="store",type="int",dest="nPU",default=40)
parser.add_option('-r',action="store",type="float",dest="radius",default=0.4)
parser.add_option('--minPt',action="store",type="float",dest="minPt",default=50.)
parser.add_option('--maxPt',action="store",type="float",dest="maxPt",default=200.)
parser.add_option('--minEta',action="store",type="float",dest="minEta",default=0.)
parser.add_option('--maxEta',action="store",type="float",dest="maxEta",default=2.5)
parser.add_option('--sample',action="store",type="string",dest="sample",default="Pythia QCD")

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

 
def makeResponseVsPu(h2, hname, rebin):

    hmean = ROOT.TH1F(hname+'_mean',hname+'_mean',100/rebin,0,100)
    hrms  = ROOT.TH1F(hname+'_rms',hname+'_rms',100/rebin,0,100)
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
        #py.Fit("fitfun","QRN")
        #meanfit    = fitfun.GetParameter(1)
        #meanfiterr = fitfun.GetParError(1)
        #sigmafit     = fitfun.GetParameter(2)
        #sigmafiterr  = fitfun.GetParError(2)

        if (py.GetEntries()>0):
            j = j + 1
            x  = (px.GetXaxis().GetBinUpEdge(lastbin) + px.GetXaxis().GetBinLowEdge(firstbin))/2
            ex = (px.GetXaxis().GetBinUpEdge(lastbin) - px.GetXaxis().GetBinLowEdge(firstbin))/2

            thisbin=hmean.FindBin(x)
            if ('pt' in h2.GetName()):
                hmean.Fill(x,mean)
                hrms.Fill(x,rms/(1+mean))
                hmean.SetBinError(thisbin,meanerr)
                hrms.SetBinError(thisbin,rmserr/(1+mean))
            else:
                hmean.Fill(x,mean)
                hrms.Fill(x,rms)
                hmean.SetBinError(thisbin,meanerr)
                hrms.SetBinError(thisbin,rmserr)
    
            #hmeanfit.Fill(x,meanfit)
            #hsigmafit.Fill(x,sigmafit)
            #hmeanfit.SetBinError(thisbin,meanfiterr)
            #hsigmafit.SetBinError(thisbin,sigmafiterr)

    #return hmean, hrms, hmeanfit, hsigmafit
    return hmean, hrms



if __name__ == '__main__':

    region = 'central'
    if (options.minEta == 2.5 and options.maxEta == 5):
        region = 'forward'

    filenames = ['../scripts/hQCD_Pt-80to120_a02_%s/histograms.root'%region,
                 '../scripts/hQCD_Pt-80to120_a03_%s/histograms.root'%region,
                 '../scripts/hQCD_Pt-80to120_a04_%s/histograms.root'%region,
                 '../scripts/hQCD_Pt-80to120_a06_%s/histograms.root'%region]

    #filenames = ['../scripts/hQCD_Pt-300to470_a02_%s/histograms.root'%region,
    #             '../scripts/hQCD_Pt-300to470_a03_%s/histograms.root'%region,
    #             '../scripts/hQCD_Pt-300to470_a04_%s/histograms.root'%region,
    #             '../scripts/hQCD_Pt-300to470_a06_%s/histograms.root'%region]
 
    outdir   = options.outdir
    try:
        os.mkdir(outdir)
    except:
        print 'Cannot create output directory: directory already exists'
        #sys.exit()

    yaxisTitle = 'arbitrary units'

    algos = ['GEN','PF+SK(a=0.2)' , 'PF+SK(a=0.3)', 'PF+SK(a=0.4)', 'PF+SK(a=0.6)']
    #algos = ['GEN', 'PF+SK(a=0.2)' , 'PF+SK(a=0.4)', 'PF+SK(a=0.6)', 'PF+SK(a=0.2) - SafeSub.' , 'PF+SK(a=0.4) - SafeSub.', 'PF+SK(a=0.6) - SafeSub.']

    histos  = {'GEN' : ['gen/hptraw_gen','gen/hmraw_gen'],
               'PF+SK(a=0.2)' : ['softkiller/hptraw_softkiller','softkiller/hmraw_softkiller'],
               'PF+SK(a=0.3)' : ['softkiller/hptraw_softkiller','softkiller/hmraw_softkiller'],
               'PF+SK(a=0.4)' : ['softkiller/hptraw_softkiller','softkiller/hmraw_softkiller'],
               'PF+SK(a=0.6)' : ['softkiller/hptraw_softkiller','softkiller/hmraw_softkiller']
               #'PF+SK(a=0.2) - SafeSub.' : ['softkiller/hpt_softkiller'],
               #'PF+SK(a=0.4) - SafeSub.' : ['softkiller/hpt_softkiller'],
               #'PF+SK(a=0.6) - SafeSub.' : ['softkiller/hpt_softkiller'],
               }

    gridsize = [0.2,0.3,0.4,0.6]
         
    styles = {} # color, linestyle, line width, marker style
    styles['GEN'] = [ROOT.kBlack, 1, 2, 24]
    styles['PF+SK(a=0.2)'] = [ROOT.kOrange+1, 1, 2, 20]
    styles['PF+SK(a=0.3)'] = [ROOT.kAzure, 1, 2, 21]
    styles['PF+SK(a=0.4)'] = [ROOT.kMagenta, 1, 2, 22]
    styles['PF+SK(a=0.6)'] = [ROOT.kRed, 1, 2, 23]
    styles['PF+SK(a=0.2) - SafeSub.'] = [ROOT.kOrange+1, 2, 2, 24]
    styles['PF+SK(a=0.3) - SafeSub.'] = [ROOT.kAzure, 2, 2, 24]
    styles['PF+SK(a=0.4) - SafeSub.'] = [ROOT.kMagenta, 2, 2, 24]
    styles['PF+SK(a=0.6) - SafeSub.'] = [ROOT.kRed, 2, 2, 24]

    

    # -- 1D distributions and resolution
    nre  = 1
    nrer = 1

    func = ROOT.TF1('func','gaus',3)

    f = []
    for fname in filenames:
        f.append(ROOT.TFile.Open(fname))
         
    for var in 'ptraw','mraw':

        i = 0
        j = 0

        hmeanvsnpu = []
        hrmsvsnpu = []

        # -- canvas
        c = ROOT.TCanvas('%s'%var,'%s'%var,1000,800)
        cresponse = ROOT.TCanvas('%s_response'%var,'%s_response'%var,1000,800)
        
        cmeanvsarea = ROOT.TCanvas('%s_mean_vs_area'%var,'%s_mean_vs_area'%var,1000,800)
        crmsvsarea = ROOT.TCanvas('%s_rms_vs_area'%var,'%s_rms_vs_area'%var,1000,800)
        
        cmeanvsnpu = ROOT.TCanvas('%s_mean_vs_npu'%var,'%s_mean_vs_npu'%var,1000,800)
        crmsvsnpu = ROOT.TCanvas('%s_rms_vs_npu'%var,'%s_rms_vs_npu'%var,1000,800)

        # -- legend                                               
        leg1 = ROOT.TLegend(0.7,0.68,0.97,0.92);
        leg1.SetBorderSize(0);
        leg1.SetFillStyle(0);

        leg1a = ROOT.TLegend(0.7,0.68,0.97,0.92);
        leg1a.SetBorderSize(0);
        leg1a.SetFillStyle(0);

        leg2 = ROOT.TLegend(0.70,0.50,0.93,0.92);
        leg2.SetBorderSize(0);
        leg2.SetFillStyle(0);
        leg2.SetTextSize(0.030);

        # -- graphs
        gmean  = ROOT.TGraphErrors()
        gmean.SetName('gmean%s'%var)
        grms   = ROOT.TGraphErrors()
        grms.SetName('grms%s'%var)
        gsigma = ROOT.TGraphErrors()
        gsigma.SetName('gsigma%s'%var)

        for algo in algos:
            # -- get histograms from file
            ifile = i-1
            if algo=='GEN':
                ifile = 0
            if 'pt' in var :
                h = f[ifile].Get(histos[algo][0])
                h.Rebin(2)
                h.GetXaxis().SetRangeUser(options.minPt, options.maxPt)
                h.GetXaxis().SetTitle('p^{T} (GeV)')

                hr = f[ifile].Get(histos[algo][0].replace('_','_response_'))
                hr.Rebin(2)
                hr.GetXaxis().SetTitle('p^{T} - p^{T}_{gen}(GeV)')

                h2r = f[ifile].Get(histos[algo][0].replace('_','_response_vs_npu_'))

            else:
                h = f[ifile].Get(histos[algo][1])
                h.Rebin(nre)
                h.GetXaxis().SetRangeUser(0,60)
                h.GetXaxis().SetTitle('m (GeV)')

                hr = f[ifile].Get(histos[algo][1].replace('_','_response_'))
                hr.Rebin(nrer)
                hr.GetXaxis().SetRangeUser(-50,50)
                hr.GetXaxis().SetTitle('m - m_{gen}(GeV)')

                h2r = f[ifile].Get(histos[algo][1].replace('_','_response_vs_npu_'))


            # 1-D distributions
            h.GetYaxis().SetTitle(yaxisTitle)
            h.GetYaxis().SetTitleOffset(1.6)
            h.SetLineColor(styles[algo][0])
            h.SetLineStyle(styles[algo][1])
            h.SetLineWidth(styles[algo][2])
            h.SetMarkerStyle(styles[algo][3])
            h.SetMarkerColor(styles[algo][0])
            ymax = h.GetMaximum()

            leg1.AddEntry(h,algo,'PL')
            c.cd()

            if (i == 0):
                h.GetYaxis().SetRangeUser(0,ymax*1.7)
                h.Draw("hp")
            else:
                h.Draw("hpsame")

            i = i+1


            # response 
            if algo != 'GEN':

                # 2-D ( response vs npv)
                htmp1, htmp2 = makeResponseVsPu(h2r, '%s_algo_%d'%(var,j), 5)
                hmeanvsnpu.append(htmp1)
                hrmsvsnpu.append(htmp2)
                
                for hh in hr, hmeanvsnpu[j], hrmsvsnpu[j]:
                    hh.GetXaxis().SetTitleOffset(1.0)
                    hh.GetYaxis().SetTitleOffset(1.3)
                    hh.GetYaxis().SetTitle(yaxisTitle)
                    hh.SetLineColor(styles[algo][0])
                    hh.SetLineStyle(styles[algo][1])
                    hh.SetLineWidth(styles[algo][2])
                    hh.SetMarkerStyle(styles[algo][3])
                    hh.SetMarkerColor(styles[algo][0])

                leg1a.AddEntry(hr,algo,'PL')
            
                func.SetRange(hr.GetMean()-2*hr.GetRMS(),hr.GetMean()+2*hr.GetRMS())
                hr.Fit("func","RNQ")
                mean    = hr.GetMean()
                meanerr = hr.GetMeanError()
                rms     = hr.GetRMS()
                rmserr  = hr.GetRMSError()
                sigma   = func.GetParameter(2)
                sigmaerr = func.GetParError(2)
            
                if (j == 0): 
                    print 'ALGO        <#Delta>(GeV)    RMS(GeV)    SIGMA(GeV)'
                print ('%s'%algo).ljust(18), ('    %.1f    %.1f    %.1f'%(mean, rms, sigma))

                gmean.SetPoint(j,gridsize[j],mean)
                grms.SetPoint(j,gridsize[j],rms)
                gsigma.SetPoint(j, gridsize[j],sigma)
                
                leg2.AddEntry(hr,algo,'PL')
                if 'pt' in var: 
                    leg2.AddEntry(0,'<#Deltap^{T}>=%.1f GeV'%mean,'')
                else:
                    leg2.AddEntry(0,'<#Deltam>=%.1f GeV'%mean,'')
                leg2.AddEntry(0,'RMS=%.1f GeV'%rms,'')
                leg2.AddEntry(0,'','')

                if (j == 0):

                    cresponse.cd()
                    hr.DrawNormalized("hp")
                    yrmax = hr.GetMaximum()
                    hr.GetYaxis().SetRangeUser(0, yrmax*1.5)
                    hr.DrawNormalized("hp")

                    cmeanvsnpu.cd()
                    cmeanvsnpu.SetGridx()
                    cmeanvsnpu.SetGridy()                    
                    hmeanvsnpu[j].GetXaxis().SetRangeUser(0,55)
                    hmeanvsnpu[j].GetXaxis().SetTitle('n_{PV}')
                    if 'pt' in var:
                        hmeanvsnpu[j].GetYaxis().SetRangeUser(-0.2, 0.2)
                        hmeanvsnpu[j].GetYaxis().SetTitle('(p^{T}-p^{T}_{gen})/p^{T}_{gen}')
                    else:
                        hmeanvsnpu[j].GetYaxis().SetRangeUser(-20, 20)
                        hmeanvsnpu[j].GetYaxis().SetTitle('m - m_{gen}(GeV)')
                    hmeanvsnpu[j].Draw('ep')

                    crmsvsnpu.cd()
                    crmsvsnpu.SetGridx()
                    crmsvsnpu.SetGridy()                    
                    hrmsvsnpu[j].GetXaxis().SetRangeUser(0,55)
                    hrmsvsnpu[j].GetXaxis().SetTitle('n_{PV}')
                    if 'pt' in var:
                        hrmsvsnpu[j].GetYaxis().SetRangeUser(0, 0.3)
                        hrmsvsnpu[j].GetYaxis().SetTitle('RMS(p^{T}/p^{T}_{gen})')
                    else:
                        hrmsvsnpu[j].GetYaxis().SetRangeUser(0, 20)
                        hrmsvsnpu[j].GetYaxis().SetTitle('RMS(m - m_{gen}) (GeV)')
                    hrmsvsnpu[j].Draw('ep')
                else:
                    cresponse.cd()
                    hr.DrawNormalized("hpsame")

                    cmeanvsnpu.cd()
                    hmeanvsnpu[j].Draw('epsame')

                    crmsvsnpu.cd()
                    hrmsvsnpu[j].Draw('epsame')
                    
                j = j+1
            

        cmeanvsarea.cd()
        cmeanvsarea.SetGridx()
        cmeanvsarea.SetGridy()
        gmean.GetHistogram().SetXTitle('SoftKiller grid spacing a')
        if 'pt' in var:
            gmean.GetHistogram().SetYTitle('<#Delta p_{T}> (GeV)')
        else:
            gmean.GetHistogram().SetYTitle('<#Deltam> (GeV)')
        gmean.GetHistogram().GetYaxis().SetRangeUser(-30,30)
        gmean.SetMarkerStyle(33)
        gmean.Draw("apl")
        
        crmsvsarea.cd()
        crmsvsarea.SetGridx()
        crmsvsarea.SetGridy()
        grms.GetHistogram().SetXTitle('SoftKiller grid spacing a')
        if 'pt' in var:
            grms.GetHistogram().SetYTitle('RMS(p_{T}-p_{T,gen}) (GeV)')
            grms.GetHistogram().GetYaxis().SetRangeUser(0,20)
        else:
            grms.GetHistogram().SetYTitle('RMS(m-m_{gen}) (GeV)')
            grms.GetHistogram().GetYaxis().SetRangeUser(0,20)
        grms.SetMarkerStyle(33)
        grms.Draw("apl")
        
    
        for can in c, cresponse, cmeanvsarea, crmsvsarea,  cmeanvsnpu, crmsvsnpu:
            CMS_lumi.CMS_lumi(can, 4, 0)
            latex0.Draw()
            latex1.Draw()
            latex2.Draw()
            latex3.Draw()
            latex4.Draw()
            
        c.cd()
        leg1.Draw()
        
        cresponse.cd()
        leg2.Draw()
    
        cmeanvsnpu.cd()
        leg1a.Draw()

        crmsvsnpu.cd()
        leg1a.Draw()


        for typ in '.png','.pdf','.C', '.root':
            for can in c, cresponse, cmeanvsarea, crmsvsarea, cmeanvsnpu, crmsvsnpu:
                can.SaveAs(outdir+'/'+can.GetName()+typ)
