#========Data-MC comparison for Electron=====#
import ROOT
import os
import math
import Gain_Function
from ROOT import TFile, TH1F, gDirectory, TCanvas, TPad, TProfile,TGraph, TGraphAsymmErrors,TGraphErrors, TF1
from ROOT import TH1D, TH1, TH1I, TH2Poly
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

#======================
#The fraction of the ionizations and neutral impurities ( K>120 )
# Three regions model #https://slideplayer.com/slide/4784387/
def Desity_of_impurities(Temperature,I_C):
    array1 = array("f",[])
    Charge_carrier_density = I_C * ( math.pow(2.718,0.01*((1/(2* KB*120)- (1/(2* KB*Temperature))) )))
    if(Charge_carrier_density<=1):
        a = 1
    else:
        a = Charge_carrier_density
    print 'I_C: '+str(I_C)
    print 'Charge_carrier_density: '+str(a)
    print 'I_C/Charge_carrier_density: '+str(I_C/a)
    #Neutral_impurities_density = I_C - Charge_carrier_density
    #print 'Charge_carrier_density: '+str(Charge_carrier_density)
    #print 'Neutral_impurities_density: '+str(Neutral_impurities_density)
    array1.append(a)
    return array1
#Desity_of_impurities(30)

def Ionization_impurities_scattering(I_C,Temperature):
    Hole1 = ( (1.59*math.pow(10,20))/(I_C) )* math.pow(( math.log( (2.8*math.pow(10,20))/(I_C) )  ),-1)
    Hole2 = ( (1.02*math.pow(10,21))/(I_C) )* math.pow(( math.log( (3.44*math.pow(10,18))/(I_C) )  ),-1)
    #print 'Ionization_impurities_mobility: '+str(Hole1+Hole2)
    return (Hole1+Hole2)*(math.log(math.pow(77,2))/math.pow(77,1.5))*(math.pow(Temperature,1.5)/math.log(math.pow(Temperature,2)))

def Phonon_impurities_scattering(T):
    #print 'Phonon_impurities_mobility: '+str(7.77*math.pow(10,7)/math.pow(T,1.5))
    return 7.77*math.pow(10,7)/math.pow(T,1.5)

def Neutral_impurities_scattering(I_C,I_D,T):
    print 'Neutral_impurities_concentration: '+str(I_C)
    print 'Forzen: '+str(I_D)
    Fixed_Neutral_impurities = ( ( 2.31*math.pow(10,18) + 2.36*math.pow(10,20) )/(I_C+I_D) )
    Real_Neutral_impurities = 0.82 * Fixed_Neutral_impurities * (0.228*math.pow(T,0.5) + 0.976*math.pow(T,-0.5))
    #print '0.82 * Fixed_Neutral_impurities: '+str(0.82 * Fixed_Neutral_impurities)
    #print 'Neutral_impurities_mobility: '+str(Real_Neutral_impurities)
    return Real_Neutral_impurities

def Total_Mobility(array_Mobility):
    Mobility = 0
    for i in range(len(array_Mobility)):
        Mobility = Mobility + (1/array_Mobility[i])
    #print 'Tem_Mobility: '+str(1/Mobility)
    return (1/Mobility)

def Real_Mobility(Tem_Mobility,Electric_field): #The Saturation condition(Two-region model)
    array1 = array("f",[])
    E_Sat = math.pow(10,7)/Tem_Mobility  #Find out the Saturation
    Mobility = (math.pow(10,7))*(1+(E_Sat/Electric_field)) / (Electric_field)
    #print 'Great!'
    if(Electric_field<E_Sat):
        #print 'Velocity: '+str(Tem_Mobility*Electric_field)
        array1.append(Tem_Mobility*Electric_field)
        array1.append(Tem_Mobility)
        array1.append(E_Sat)
        print 'Real_Mobility: '+str(Tem_Mobility)

    else:
        array1.append(math.pow(10,7))
        array1.append(Mobility)
        array1.append(E_Sat)
        print 'Real_Mobility: '+str(Tem_Mobility)
    return array1


#array1 = array("f",[])
#array1.append(Ionization_impurities_scattering(Desity_of_impurities(Temperature,Impurities_Ionization_Concentration)[0]))
#array1.append(Phonon_impurities_scattering(Temperature))
#array1.append(Neutral_impurities_scattering(Impurities_Neutral_Concentration,Temperature))
#Using_Mobility = Real_Mobility(Total_Mobility(array1),Electric_field)
#MFP_Ele = Gain_Function.Mean_free_path(Using_Mobility[0],float(Gain_Function.Relaxation_time(Using_Mobility[1],0.21,CC)))
#MFP_Hol = Gain_Function.Mean_free_path(Using_Mobility[0],float(Gain_Function.Relaxation_time(Using_Mobility[1],0.12,CC)))

#Gain_Function.Ionization_rate(0,MFP_Ele,MFP_Hol,0.01,Electric_field,Temperature)

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
Esat = array("f",[])

yarray=array("f",[])
yarray1=array("f",[])
yarray2=array("f",[])
yarray3=array("f",[])

for i in range(14,18):
    for Electric_field in range(5,110,5):
        Temperature = 4
        array_Mob = array("f",[])
        array_Mob.append(Ionization_impurities_scattering(Desity_of_impurities(Temperature,math.pow(10,i))[0],Temperature))
        array_Mob.append(Phonon_impurities_scattering(Temperature))
        array_Mob.append(Neutral_impurities_scattering(Impurities_Neutral_Concentration,math.pow(10,i)/(Desity_of_impurities(Temperature,math.pow(10,i))[0]),Temperature))
        print 'Ionization_impurities_scattering: '+str(array_Mob[0])
        print 'Phonon_impurities_scattering: '+str(array_Mob[1])
        print 'Neutral_impurities_scattering: '+str(array_Mob[2])
        Using_Mobility = Real_Mobility(Total_Mobility(array_Mob),Electric_field)
        MFP_Ele = Gain_Function.Mean_free_path(Using_Mobility[0],float(Gain_Function.Relaxation_time(Using_Mobility[1],0.21,CC)))
        MFP_Hol = Gain_Function.Mean_free_path(Using_Mobility[0],float(Gain_Function.Relaxation_time(Using_Mobility[1],0.12,CC)))
        Pre_selection = Gain_Function.Ionization_rate(1,MFP_Ele,MFP_Hol,0.01,Electric_field,Temperature,math.pow(10,i))[4]
        if(Electric_field==10):
            Esat.append(Using_Mobility[2])
        if(i==14):
            if(Pre_selection<=0):
                yarray.append(0)
            else:
                yarray.append(Pre_selection)
        if(i==15):
            if(Pre_selection<=0):
                yarray1.append(0)
            else:
                yarray1.append(Pre_selection)
        if(i==16):
            if(Pre_selection<=0):
                yarray2.append(0)
            else:
                yarray2.append(Pre_selection)
        if(i==17):
            if(Pre_selection<=0):
                yarray3.append(0)
            else:
                yarray3.append(Pre_selection)

        xarray.append(Electric_field)
        xarrayerror.append(0)
        yarrayerror.append(0)


c = TCanvas("c1", "c1",0,0,500,500)

gStyle.SetOptFit()
gr = TGraphErrors(21,xarray,yarray,xarrayerror,yarrayerror)
gr1 = TGraphErrors(21,xarray,yarray1,xarrayerror,yarrayerror)
gr2 = TGraphErrors(21,xarray,yarray2,xarrayerror,yarrayerror)
gr3 = TGraphErrors(21,xarray,yarray3,xarrayerror,yarrayerror)

gr.SetLineStyle(1)
gr.SetMarkerColor(2)
gr.SetMarkerStyle(8)
gr.SetLineWidth(2)
gr1.SetLineStyle(1)
gr1.SetMarkerColor(3)
gr1.SetMarkerStyle(8)
gr1.SetLineWidth(2)
gr2.SetLineStyle(1)
gr2.SetMarkerColor(4)
gr2.SetMarkerStyle(8)
gr2.SetLineWidth(2)
gr3.SetLineStyle(1)
gr3.SetMarkerColor(5)
gr3.SetMarkerStyle(8)
gr3.SetLineWidth(2)

gr1.SetTitle(";E(V/cm) ;Ionization Rate(/cm) ")
gr.GetHistogram().SetMaximum(3000)
gr.GetHistogram().SetMinimum(0)
gr.GetXaxis().SetLimits(0,10)
gr.GetXaxis().CenterTitle()
gr.GetYaxis().CenterTitle()
#================================
gr.GetYaxis().SetTitleSize(0.02)
gr.GetXaxis().SetTitleSize(0.02)
gr.GetXaxis().SetLabelSize(0.02)
gr.GetYaxis().SetLabelSize(0.02)
gr.GetXaxis().SetLabelFont(22)
gr.GetYaxis().SetLabelFont(22)
gr.GetXaxis().SetTitleColor(1)
gr.GetYaxis().SetTitleColor(1)

leg1=TLegend(0.2,0.6,0.5,0.8)
leg1.SetFillColor(0)
leg1.SetFillStyle(0)
leg1.SetTextSize(0.05)
leg1.SetBorderSize(0)
leg1.SetTextFont(22)
leg1.AddEntry("","GeIA group","")
leg1.AddEntry(gr,"10^{14}:"+str(Esat[0]),"lp")
leg1.AddEntry(gr1,"10^{15}"+str(Esat[1]),"lp")
leg1.AddEntry(gr2,"10^{16}"+str(Esat[2]),"lp")
leg1.AddEntry(gr3,"10^{17}"+str(Esat[3]),"lp")

gr1.Draw("ALP")
gr.Draw("LPsame")
gr2.Draw("LPsame")
gr3.Draw("LPsame")
#c.SetLogx()
#c.SetLogy()
c.Draw()
leg1.Draw()
c.Print("Impurity_IR_4K_16_to_19.pdf")














































