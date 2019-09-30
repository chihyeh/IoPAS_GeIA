#========Data-MC comparison for Electron=====#
import ROOT
import os
import math
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
#=======================
#Gain factor exploration
#Principle: The electrons will go under the "impact ionizations" and give out more electrons.
#2 is the "a factor of two for give out the twice of the electrons
#Avlanche region is the region that can produce the "impact ionizations"
#The man free path is the length that how long "a electron" will keep running.
#(1) Option of electron/Hole: 0 is the electron and 1 is the hole - They have the different behaviors
Option=0 #0 is the electron and 1 is the hole
#print 'Option: '+str(Option)
#(2) Temperature is the one of the important parameters in this studies
Temperature=77 #4-77
#print 'temperature: '+str(Temperature)
#(3) Electric field is another one of the important parameters in this studies also
Electric_field=100000 #1000-100000(Lower-bound:15)
#print 'Electric_field: '+str(Electric_field*100)
#(4) Charge Constant
CC = 1.6 * math.pow(10,-19)
#(5) Boltzman COnstantKB = 8.617 * math.pow(10,-5)
KB = 8.617 * math.pow(10,-5)
#(6) Electron mass with the unit of eV/(cm2/s2)
Electron_mass = (0.51*math.pow(10,6)) / (9*math.pow(10,20))
#==============================================================Avalanche================================================================
def integrand(x):
    return math.pow(2.718,x)
def f(A,B,E,x):
    return A * np.exp(-B/(Ex))
#==============================================================Avalanche================================================================
def Ionization_rate(Option,(MFP_El),(MFP_Ho),ionization_energy,E_x,Temperature):
    array1 = array("f",[])
    MFP_E = MFP_El
    MFP_H = MFP_Ho
    #print 'MFP_E: '+str(MFP_E)
    #print 'MFP_H: '+str(MFP_H)
    #s={n,p}
    if(Option==0):
        A_s = (1/MFP_E)
    #   print 'A_s: '+str(A_s)
    if(Option==1):
        A_s = (1/MFP_H)
#print 'A_s: '+str(A_s)
#===========
    B_s = ionization_energy * A_s
#print 'B_s: '+str(B_s)
    B_n = (ionization_energy)*(1/MFP_E)
    B_p = (ionization_energy)*(1/MFP_H)
#   print 'B_n: '+str(B_n)
#   print 'B_p: '+str(B_p)
    Z_E = 1.0 + (B_n/E_x)*math.pow(2.718,-B_n/E_x)+ (B_p/E_x)*math.pow(2.718,-B_p/E_x)
#   print 'Z_E: '+str(Z_E)
#return (A_s/Z_E)*math.pow(2.718,(-B_s/E_x))
    array1.append(A_s)
    array1.append(Z_E)
    array1.append(B_s)
    array1.append(E_x)
    array1.append((A_s/Z_E)*math.pow(2.718,(-B_s/E_x)))
    return array1
    #return 1.55* math.pow(10,7) * math.pow(2.718,-1.56*math.pow(10,6)/E_x)
#def Multuplication():

#==============================================================Avalanche================================================================
#Approximately, since there is no data for the electric_field below 10000 V/cm
#def Avalanche_region(Temperature,Electric_field):

#==============================================================Mean_free_path=================================================================
def Effective_mass_factor(Option,Temperature):#The factors of the electron and the hole
    if(Option==0):
        Effective_Mass_Of_The_Electron = 0.25 + 5 * math.pow( 10, -5) * Temperature
        EM = Effective_Mass_Of_The_Electron
    if(Option==1):
        Effective_Mass_Of_The_Hole = 0.34 + 5 * math.pow( 10, -4) * Temperature
        EM = Effective_Mass_Of_The_Hole
    #print 'Effective Mass(EM) factor: '+str(EM)
    return(EM)
#Temperature Dependence of Silicon Carrier E ective Masses with Application to FemtosecondRe ectivity Measurements(Paper)
#==========================================
def Mobility(Electric_field_1,Temperature):#Get the Mobility under the different electric fields and temperatures
    #E_Saturate = 50 + (450/63)*(Temperature-4)
    if(Temperature==77):
        E_Saturate = 2000
    if(Temperature==4):
        E_Saturate = 1000
    if(Electric_field_1>E_Saturate):
        Drift_Velocity = math.pow(10,7)
    if(Electric_field_1<=E_Saturate):
        #Constant = math.pow(10,7) / math.pow(E_Saturate,0.4798198+(0.5204107-0.4798198)/73*(Temperature-4))
#Drift_Velocity = Constant * math.pow(Electric_field_1,0.4798198+(0.5204107-0.4798198)/73*(Temperature-4))
        Drift_Velocity = 8996611.9*math.pow(Electric_field_1,0.5229874)
    #print 'Drift_Velocity: '+str(Drift_Velocity)
    #print 'E_Saturate: '+str(E_Saturate)
    Mobility = (Drift_Velocity)*(1+(Electric_field_1/E_Saturate)) / (Electric_field_1)
    #print 'Mobility: '+str(Mobility*0.0001)
# print 'Drift_Velocity: '+ str(Drift_Velocity)
    print 'Mobility: '+str(Mobility)
#return math.pow(10,7)/Electric_field_1
    return Mobility
#==========================================1
def Relaxation_time(Mobility,Effective_mass_factor,Charge_constant):#The time that the electron will bump into other electrons for the first time
    Effective_mass = Effective_mass_factor* Electron_mass
    #print 'Electron_mass: '+str(Electron_mass)
    #print 'Relaxation_time: '+str((Mobility * 9.1 * Effective_mass_factor * math.pow(10,-31)/(1.6*math.pow(10,-19)*10000)))
    return (Mobility * 9.1 * Effective_mass_factor * math.pow(10,-31)/(1.6*math.pow(10,-19)*10000))
#==========================================
def Thermal_Velocity(KB,Temperature,Effective_mass_factor):#Velocity that is caused by the thermal.
    Effective_mass = (Effective_mass_factor) * Electron_mass
    #print 'Thermal_Velocity: '+str(math.pow( (3*KB*Temperature)/Effective_mass ,0.5))
    return (math.pow( (3*KB*Temperature)/Effective_mass ,0.5))
#==========================================#==========================================#========================
def Mean_free_path(Thermal_Velocity,Relaxation_time):#The length that the particle will bump into another electrons for the first time.
    #print 'Mean_free_path: '+str((Thermal_Velocity)*(Relaxation_time))
    return (math.pow(10,7))*(Relaxation_time)
#==============================================================Gain=================================================================
def Gain(Avalanche_region,Mean_free_path):
    #print 'Avalanche_region/Mean_free_path: '+str(Avalanche_region/Mean_free_path)
    return math.pow(2,Avalanche_region/Mean_free_path)
#==========================================#==========================================#========================
#MFP_Ele = Mean_free_path(Thermal_Velocity(KB,Temperature,Effective_mass_factor(0,Temperature)),Relaxation_time(Mobility(Electric_field,Temperature), Effective_mass_factor(0,Temperature), CC))
#MFP_Hol = Mean_free_path(Thermal_Velocity(KB,Temperature,Effective_mass_factor(1,Temperature)),Relaxation_time(Mobility(Electric_field,Temperature), Effective_mass_factor(1,Temperature), CC))

f = lambda x : math.pow(2.718,-x**2)
scipy.integrate.quad(f, 0, 1)
(0.7468241328124271, 8.291413475940725e-15)




'''
def Gain_Final(Temperature,Electric_field):
    Eff_mass_factor = Effective_mass_factor(0,Temperature)
    Mob = Mobility(Electric_field,Temperature)
    MFP=Mean_free_path(Thermal_Velocity(KB,Temperature,Eff_mass_factor),Relaxation_time(Mob, Eff_mass_factor, CC))
    #print 'Gain: '+str(Gain(Avalanche_region(Temperature),MFP))
    #return Gain(Avalanche_region(Temperature,Electric_field),MFP)
    return Avalanche_region(Temperature,Electric_field)*100 * Ionization_rate_for_electron(Electric_field)
'''
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


c = TCanvas("c1", "c1",0,0,500,500)

for i in range(1,14):
    xarrayerror.append(0)
    yarrayerror.append(0)
    xarray.append(i*10000)
    MFP_Ele = Mean_free_path(Thermal_Velocity(KB,4,0.12),Relaxation_time(Mobility(i*10000,4), 0.12, CC))
    MFP_Hol = Mean_free_path(Thermal_Velocity(KB,4,0.21),Relaxation_time(Mobility(i*10000,4), 0.21, CC))
    
    Array = (Ionization_rate(0,MFP_Ele,MFP_Hol,3,i*10000,4))
    yarray.append((Array[0]/Array[1])*math.pow(2.718,(-Array[2]/Array[3])))
    #print str((Array[0]/Array[1])*math.pow(2.718,(-Array[2]/Array[3])))
    Array1 = (Ionization_rate(1,MFP_Ele,MFP_Hol,3,i*10000,4))
    yarray1.append((Array1[0]/Array1[1])*math.pow(2.718,(-Array1[2]/Array1[3])))


for i in range(1,13):
    MFP_Ele = Mean_free_path(Thermal_Velocity(KB,77,0.12),Relaxation_time(Mobility(i*10000,77), 0.12, CC))
    MFP_Hol = Mean_free_path(Thermal_Velocity(KB,77,0.21),Relaxation_time(Mobility(i*10000,77), 0.21, CC))
    Array = (Ionization_rate(0,MFP_Ele,MFP_Hol,3,i*10000,77))
    yarray2.append((Array[0]/Array[1])*math.pow(2.718,(-Array[2]/Array[3])))
    #print str((Array[0]/Array[1])*math.pow(2.718,(-Array[2]/Array[3])))
    Array1 = (Ionization_rate(1,MFP_Ele,MFP_Hol,3,i*10000,77))
    yarray3.append((Array1[0]/Array1[1])*math.pow(2.718,(-Array1[2]/Array1[3])))


'''
for i in range(1,1000000):
    Temperature = 77
    xarrayerror.append(0)
    yarrayerror.append(0)
    xarray.append(1*i)
    yarray.append(Gain_Final(Temperature,1*i))
    print 'Gain: '+str(Gain_Final(Temperature,1*i))
    Temperature1 = 4
    yarray1.append(Gain_Final(Temperature1,1*i))
    print 'Gain: '+str(Gain_Final(Temperature1,1*i))



for i in range(20000,25000):
    Temperature = 4 + (0.001 * 73 * i)
    Electric_field = 10
    xarrayerror.append(0)
    yarrayerror.append(0)
    xarray.append(Temperature)
    yarray.append(Gain_Final(Temperature,Electric_field))
    yarray1.append(Gain_Final(Temperature,Electric_field*10))
    print str(Gain_Final(Temperature,Electric_field))
    print str(Gain_Final(Temperature,Electric_field*10))
'''

gStyle.SetOptFit()
gr = TGraphErrors(12,xarray,yarray,xarrayerror,yarrayerror)
gr1 = TGraphErrors(12,xarray,yarray1,xarrayerror,yarrayerror)
gr2 = TGraphErrors(12,xarray,yarray2,xarrayerror,yarrayerror)
gr3 = TGraphErrors(12,xarray,yarray3,xarrayerror,yarrayerror)

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
gr.GetHistogram().SetMaximum(10000)
gr.GetHistogram().SetMinimum(0)
gr.GetXaxis().SetLimits(0,130000)
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

leg1=TLegend(0.2,0.7,0.5,0.9)
leg1.SetFillColor(0)
leg1.SetFillStyle(0)
leg1.SetTextSize(0.05)
leg1.SetBorderSize(0)
leg1.SetTextFont(22)
leg1.AddEntry("","GeIA group","")
leg1.AddEntry(gr,"Ge-Electron(4K)","lp")
leg1.AddEntry(gr1,"Ge-Hole(4K)","lp")
leg1.AddEntry(gr2,"Ge-Electron(77K)","lp")
leg1.AddEntry(gr3,"Ge-Hole(77K)","lp")

gr1.Draw("ALP")
gr.Draw("LPsame")
gr2.Draw("LPsame")
gr3.Draw("LPsame")
c.SetLogx()
#c.SetLogy()
c.Draw()
leg1.Draw()
c.Print("IR_77K_4K.pdf")


















































