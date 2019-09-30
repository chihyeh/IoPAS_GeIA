#========Data-MC comparison for Electron=====#
import ROOT
import os
import math
from ROOT import TFile, TH1F, gDirectory, TCanvas, TPad, TProfile,TGraph, TGraphAsymmErrors,TGraphErrors, TF1
from ROOT import TH1D, TH1, TH1I, TH2Poly, TLine
from ROOT import gStyle
from ROOT import gROOT
from ROOT import TStyle
from ROOT import TLegend
from ROOT import TMath
from ROOT import TPaveText
from ROOT import TLatex
from array import array
from scipy.integrate import quad
import numpy as np
import scipy.integrate
#=======================
#Concentration(#/cm^3)
def Number_of_Atom(Avalanche_region,Concentration):
    return (Avalanche_region) * (Concentration)

#============================================================Studies at there!!================================================================
#################################
atlasStyle= TStyle("ATLAS","Atlas style")

icol=0
atlasStyle.SetFrameBorderMode(icol)
atlasStyle.SetCanvasBorderMode(icol)
atlasStyle.SetPadBorderMode(icol)
atlasStyle.SetPadColor(icol)
atlasStyle.SetCanvasColor(icol)
atlasStyle.SetStatColor(icol)
#atlasStyle.SetFillColor(icol)

# set the paper & margin sizes
atlasStyle.SetPaperSize(20,26)
atlasStyle.SetPadTopMargin(0.05)
#atlasStyle.SetPadRightMargin(0.05)
atlasStyle.SetPadRightMargin(0.12)
atlasStyle.SetPadBottomMargin(0.16)
atlasStyle.SetPadLeftMargin(0.16)
# use large fonts
#Int_t font=72.
font=42
tsize=0.04
atlasStyle.SetTextFont(font)


atlasStyle.SetTextSize(tsize)
atlasStyle.SetLabelFont(font,"x")
atlasStyle.SetTitleFont(font,"x")
atlasStyle.SetLabelFont(font,"y")
atlasStyle.SetTitleFont(font,"y")
atlasStyle.SetLabelFont(font,"z")
atlasStyle.SetTitleFont(font,"z")

atlasStyle.SetLabelSize(tsize,"x")
atlasStyle.SetTitleSize(tsize,"x")
atlasStyle.SetLabelSize(tsize,"y")
atlasStyle.SetTitleSize(tsize,"y")
atlasStyle.SetLabelSize(tsize,"z")
atlasStyle.SetTitleSize(tsize,"z")


#use bold lines and markers
atlasStyle.SetMarkerStyle(20)
atlasStyle.SetMarkerSize(1.2)
atlasStyle.SetHistLineWidth(2)
atlasStyle.SetLineStyleString(2,"[12 12]") # postscript dashes
#get rid of X error bars and y error bar caps
#atlasStyle.SetErrorX(0.001)

#do not display any of the standard histogram decorations
atlasStyle.SetOptTitle(0)
#atlasStyle.SetOptStat(1111)
#atlasStyle.SetOptStat()
atlasStyle.SetOptFit()
#atlasStyle.SetOptFit(0)

# put tick marks on top and RHS of plots
atlasStyle.SetPadTickX(1)
atlasStyle.SetPadTickY(1)

#gStyle.SetPadTickX(1)
#gStyle.SetPadTickY(1)

atlasStyle.SetOptFit(0)
#################################
gROOT.Reset()
gROOT.SetStyle("ATLAS")
gROOT.ForceStyle()
#################################

xarray=array("f",[])
yarray=array("f",[])
xarrayerror=array("f",[])
yarrayerror=array("f",[])

yarray=array("f",[])
yarray1=array("f",[])
yarray2=array("f",[])
yarray3=array("f",[])

Avalanche_region = 2.8 * math.pow(10,-5) #(cm)
Avalanche_region1 = 3 * math.pow(10,-5) #(cm)

for i in range(0,17):
    xarrayerror.append(0)
    yarrayerror.append(0)
    xarray.append(math.pow(10,10+i*0.5)) #cm^3
    yarray.append(Number_of_Atom(Avalanche_region,math.pow(math.pow(10,10+i*0.5),0.33)))
    yarray1.append(Number_of_Atom(Avalanche_region1,math.pow(math.pow(10,10+i*0.5),0.33)))
    print str(Number_of_Atom(Avalanche_region1,math.pow(math.pow(10,10+i*0.5),0.33)))
    print str(Number_of_Atom(Avalanche_region,math.pow(math.pow(10,10+i*0.5),0.33)))

c = TCanvas("c1", "c1",0,0,500,500)


gStyle.SetOptFit()
gr = TGraphErrors(17,xarray,yarray,xarrayerror,yarrayerror)
gr1 = TGraphErrors(17,xarray,yarray1,xarrayerror,yarrayerror)

gr.SetLineStyle(1)
gr.SetMarkerColor(2)
gr.SetMarkerStyle(8)
gr.SetLineWidth(2)

gr1.SetLineStyle(1)
gr1.SetMarkerColor(3)
gr1.SetMarkerStyle(8)
gr1.SetLineWidth(2)

gr1.SetTitle(";Concentration(1/cm^{3}) ;# of impurities ")
gr1.GetHistogram().SetMaximum(30)
gr1.GetHistogram().SetMinimum(0)
gr1.GetXaxis().SetLimits(math.pow(10,10+0*0.5),math.pow(10,10+16*0.5))

leg1=TLegend(0.2,0.7,0.5,0.9)
leg1.SetFillColor(0)
leg1.SetFillStyle(0)
leg1.SetTextSize(0.05)
leg1.SetBorderSize(0)
leg1.SetTextFont(22)
leg1.AddEntry("","GeIA group","")
leg1.AddEntry("","Avalanche region:","")
leg1.AddEntry(gr,str(Avalanche_region)+"cm(Square)","lp")
leg1.AddEntry(gr1,str(Avalanche_region1)+"cm(Circle)","lp")

line = TLine(math.pow(10,15),0,math.pow(10,15),Number_of_Atom(Avalanche_region,math.pow(math.pow(10,15),0.33)))
line1 = TLine(math.pow(10,16),0,math.pow(10,16),Number_of_Atom(Avalanche_region,math.pow(math.pow(10,16),0.33)))
line3 = TLine(math.pow(10,15),0,math.pow(10,16),0)
line3.SetLineWidth(10)

gr1.Draw("ALP")
gr.Draw("LPsame")
line.Draw("Lsame")
line1.Draw("Lsame")
c.SetLogx()
#c.SetLogy()
c.Draw()
leg1.Draw()
c.Print("Results/Concentration_problem.pdf")


















































