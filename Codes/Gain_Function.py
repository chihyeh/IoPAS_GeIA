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
def Ionization_rate(Option,(MFP_El),(MFP_Ho),ionization_energy,E_x,Temperature,Impurity_Concentration):
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
    print 'Ionization_rate: '+str((A_s/Z_E)*math.pow(2.718,(-B_s/E_x)))
    return array1
    #return 1.55* math.pow(10,7) * math.pow(2.718,-1.56*math.pow(10,6)/E_x)
#==========================================1
def Relaxation_time(Mobility,Effective_mass_factor,Charge_constant):#The time that the electron will bump into other electrons for the first time
    #Effective_mass = Effective_mass_factor* Electron_mass
    #print 'Electron_mass: '+str(Electron_mass)
    #print 'Relaxation_time: '+str((Mobility * 9.1 * Effective_mass_factor * math.pow(10,-31)/(1.6*math.pow(10,-19)*10000)))
    return (Mobility * 9.1 * Effective_mass_factor * math.pow(10,-31)/(1.6*math.pow(10,-19)*10000))
#==========================================
def Thermal_Velocity(KB,Temperature,Effective_mass_factor):#Velocity that is caused by the thermal.
    Effective_mass = (Effective_mass_factor) * Electron_mass
    #print 'Thermal_Velocity: '+str(math.pow( (3*KB*Temperature)/Effective_mass ,0.5))
    return (math.pow( (3*KB*Temperature)/Effective_mass ,0.5))
#==========================================#==========================================#========================
def Mean_free_path(Velocity,Relaxation_time):#The length that the particle will bump into another electrons for the first time.
    #print 'Mean_free_path: '+str((Thermal_Velocity)*(Relaxation_time))
    return (Velocity)*(Relaxation_time)
#==============================================================Gain=================================================================
def Gain(Avalanche_region,Mean_free_path):
    #print 'Avalanche_region/Mean_free_path: '+str(Avalanche_region/Mean_free_path)
    return math.pow(2,Avalanche_region/Mean_free_path)
#==========================================#==========================================#========================


















































