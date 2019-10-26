#========Data-MC comparison for Electron=====#
import ROOT
import os
import math
import Gain4_Singal_threshold
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
from numpy import log as ln
#=======================
#Gain factor exploration
#Principle: The electrons will go under the "impact ionizations" and give out more electrons.
#2 is the "a factor of two for give out the twice of the electrons
#Avlanche region is the region that can produce the "impact ionizations"
#The man free path is the length that how long "a electron" will keep running.
#(1) Option of electron/Hole: 0 is the electron and 1 is the hole - They have the different behaviors
#Option=0 #0 is the electron and 1 is the hole
#print 'Option: '+str(Option)
#(2) Temperature is the one of the important parameters in this studies
#Temperature=77 #4-77
#print 'temperature: '+str(Temperature)
#(3) Electric field is another one of the important parameters in this studies also
#Electric_field=100 #1000-100000(Lower-bound:15)
#print 'Electric_field: '+str(Electric_field*100)
#(4) Charge Constant
CC = 1.6 * math.pow(10,-19)
#(5) Boltzman COnstantKB = 8.617 * math.pow(10,-5) eV/K
KB = 8.617 * math.pow(10,-5)
#(6) Electron mass with the unit of eV/(cm2/s2)
Electron_mass = (0.51*math.pow(10,6)) / (9*math.pow(10,20))
#(7) P-type Impurities concentration
#Impurities_Ionization_Concentration = math.pow(10,10)
Impurities_Neutral_Concentration = 2*math.pow(10,15)
#(8) Crystal_Volume (cm^3)
#======================
#The fraction of the ionizations and neutral impurities ( K>120 )
# Three regions model #https://slideplayer.com/slide/4784387/
def Desity_of_impurities(Temperature,I_C):
    array1 = array("f",[])
    Charge_carrier_density = I_C * ( math.pow(2.718,0.0106*((1/(2*KB*120)- (1/(2*KB*Temperature))) )))
    if(Charge_carrier_density<=1):
        a = 1
    else:
        a = Charge_carrier_density
            #   print 'I_C: '+str(I_C)
            # print 'Charge_carrier_density: '+str(a)
#   print 'I_C/Charge_carrier_density: '+str(I_C/a)
    array1.append(a)
    return array1
#################################
def Raw_Singal_threshold(G,Nb):
# Formula Ns > ( 3 ((Nb*G)^(1/2)) /G )
    Ns = 3 * (math.pow( Nb/G , 0.5))
    if(G==1):
        print 'G==1: '+str(Ns*3)
        

    return Ns*3

def BKG_estimation(Volume,T,Charged_impurites):
    return (0.1*Volume*Desity_of_impurities(T,math.pow(10,Charged_impurites))[0])/(100) # #number/(cm^3*(mu*s))


#print str(Raw_Singal_threshold(100,BKG_estimation(7.5*7.5*3,4,10))) + "eV"
#print str(Raw_Singal_threshold(100,BKG_estimation(7.5*7.5*3,77,10))) + "eV"


'''
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
'''
xarray=array("f",[])
xarrayerror=array("f",[])
yarrayerror=array("f",[])
Esat = array("f",[])

yarray=array("f",[])
yarray1=array("f",[])
yarray2=array("f",[])
yarray3=array("f",[])

for i in range(4,120):
    xarray.append(1000/i)
    yarray.append(Desity_of_impurities(i,math.pow(10,10))[0])
    yarray1.append(Desity_of_impurities(i,math.pow(10,10))[0])
    xarrayerror.append(0)
    yarrayerror.append(0)
c = TCanvas("c1", "c1",0,0,500,500)

gStyle.SetOptFit()
gr = TGraphErrors(116,xarray,yarray,xarrayerror,yarrayerror)
gr1 = TGraphErrors(116,xarray,yarray1,xarrayerror,yarrayerror)

gr.SetLineStyle(1)
gr.SetMarkerColor(2)
gr.SetMarkerStyle(8)
gr.SetLineWidth(2)

gr1.SetLineStyle(1)
gr1.SetMarkerColor(3)
gr1.SetMarkerStyle(8)
gr1.SetLineWidth(2)

gr.SetTitle(";#frac{1000}{T}(#frac{1}{K}) ;Charge carrier density(cm^{-3}) ")
gr.GetXaxis().SetRangeUser(0,300)
gr.GetYaxis().SetRangeUser(0,math.pow(10,10))
#gr.GetHistogram().SetMaximum(6)
#gr.GetHistogram().SetMinimum(-2)
#gr.GetXaxis().SetLimits(0,1000)
gr.GetXaxis().CenterTitle()
gr.GetYaxis().CenterTitle()
#gr.GetXaxis().SetTitleOffset(1)
#gr.GetYaxis().SetTitleOffset(1)

#================================
gr.GetYaxis().SetTitleSize(0.04)
gr.GetXaxis().SetTitleSize(0.04)
gr.GetXaxis().SetLabelSize(0.04)
gr.GetYaxis().SetLabelSize(0.04)
gr.GetXaxis().SetLabelFont(22)
gr.GetYaxis().SetLabelFont(22)
gr.GetXaxis().SetTitleColor(1)
gr.GetYaxis().SetTitleColor(1)

leg1=TLegend(0.5,0.7,0.7,0.9)
leg1.SetFillColor(0)
leg1.SetFillStyle(0)
leg1.SetTextSize(0.05)
leg1.SetBorderSize(0)
leg1.SetTextFont(22)
leg1.AddEntry("","GeIA group","")
leg1.AddEntry("","#propto e^{-#frac{E}{2k_{B}T}}","")


line_120K =TLine(1000/120-2,0,1000/120+2,Desity_of_impurities(120,math.pow(10,10))[0])
line_77K =TLine(1000/77-2,0,1000/77+2,Desity_of_impurities(77,math.pow(10,10))[0])
line_4K =TLine(1000/4-2,0,1000/4+2,Desity_of_impurities(4,math.pow(10,10))[0])
line_120K.SetLineWidth(2)
line_77K.SetLineWidth(2)
line_4K.SetLineWidth(2)



t =  TLatex(1000/120,Desity_of_impurities(120,math.pow(10,10))[0],"120K:10^{10}(cm^{-3})");
t1 =  TLatex(1000/77,Desity_of_impurities(77,math.pow(10,10))[0]/10,"77K:"+str(round(Desity_of_impurities(77,math.pow(10,10))[0]/math.pow(10,9),1))+"#times 10^{9}(cm^{-3})")
t2 =  TLatex(1000/4-75,Desity_of_impurities(4,math.pow(10,10))[0],"4K:"+str(round(Desity_of_impurities(4,math.pow(10,10))[0]/math.pow(10,3),1))+"#times 10^{3}(cm^{-3})")

c.Draw()
gr.Draw("ALP")
leg1.Draw()
line_120K.Draw("same")
line_77K.Draw("same")
line_4K.Draw("same")
t.Draw("same")
t1.Draw("same")
t2.Draw("same")

c.SetLogy()
c.Print("Impurity_concentration.pdf")














































