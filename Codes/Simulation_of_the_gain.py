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
import decimal
#=======================
#Gain factor exploration
#Principle: The electrons will go under the "impact ionizations" and give out more electrons.
#2 is the "a factor of two for give out the twice of the electrons
#Avlanche region is the region that can produce the "impact ionizations"
#The man free path is the length that how long "a electron" will keep running.
#==========================================
#(1) Option of electron/Hole: 0 is the electron and 1 is the hole - They have the different behaviors
Option=0 #0 is the electron and 1 is the hole
print 'Option: '+str(Option)
#(2) Temperature is the one of the important parameters in this studies
Temperature=77 #4-77
print 'temperature: '+str(Temperature)
#(3) Electric field is another one of the important parameters in this studies also
Electric_field=1000 #5-10000(Lower-bound:15)
print 'Electric_field: '+str(Electric_field)
#(4) Charge Constant
CC = 1.6 * math.pow(10,-19)
#(5) Boltzman COnstant
KB = 8.617 * math.pow(10,-5)
#(6) Electron mass with the unit of eV/(m2/s2)
Electron_mass = (0.51*math.pow(10,6)) / (9*math.pow(10,16))
#==============================================================Avalanche======================================================================
#Approximately, since there is no data for the electric_field below 10000 V/cm
def Avalanche_region(Temperature):
    Slope = (25) / (73)
    Region_length = 5 + Temperature * Slope
    print 'Avalanche_region_length: '+str(Region_length * math.pow(10,-6))
    return Region_length * math.pow(10,-6)
#==============================================================Mean_free_path=================================================================
def Effective_mass_factor(Option,Temperature):#The factors of the electron and the hole
    if(Option==0):
        Effective_Mass_Of_The_Electron = 0.25 + 5 * math.pow( 10, -5) * Temperature
        EM = Effective_Mass_Of_The_Electron
    if(Option==1):
        Effective_Mass_Of_The_Hole = 0.34 + 5 * math.pow( 10, -4) * Temperature
        EM = Effective_Mass_Of_The_Hole
    print 'Effective Mass(EM) factor: '+str(EM)
    return(EM)
#Temperature Dependence of Silicon Carrier E ective Masses with Application to FemtosecondRe ectivity Measurements(Paper)
#==========================================
def Drift_Velocity_at_10V(Option,Temperature):#The one point of the drift velocity for fitting
    if(Option==0):
        Tem_term = 4*math.pow(Temperature,-1)
        pow_term = 1.66
    print 'Tem_term: '+str(Tem_term)
    if(Option==1):
        Tem_term = 77*math.pow(Temperature,-1)
        pow_term = 2.33
            #print 'Tem_term: '+str(Tem_term)
            #print 'pow_term: '+str(pow_term)
            #print 'Drift_Velocity_at_10V: '+str(math.pow(10,6) * math.pow(Tem_term,pow_term))
    return math.pow(10,6) * math.pow(Tem_term,pow_term)
#=======================
def A_Constant(Vd):#Get the one parameter of the function
    print 'A_Constant: '+ str((Vd - math.pow(10,7)) / math.pow(-9990,2))
    return (Vd - math.pow(10,7)) / math.pow(-9990,2)
#==========================================
def Mobility(A_Constant,Electric_field,Temperature):#Get the Mobility under the different electric fields and temperatures
    if(Electric_field>100 and Temperature ==8):
        Drift_Velocity = 10000000
    if(Electric_field>100 and Temperature ==77):
        Drift_Velocity = 10000000
    else:
        Drift_Velocity = A_Constant * math.pow(Electric_field-10000,2) + math.pow(10,7)
    Mobility = Drift_Velocity / Electric_field
    print 'Drift_Velocity: '+ str(Drift_Velocity)
    print 'Mobility: '+ str(Mobility)
    return Mobility
#==========================================
def Relaxation_time(Mobility,Effective_mass_factor,Charge_constant):#The time that the electron will bump into other electrons for the first time
    Effective_mass = Effective_mass_factor* Electron_mass * Charge_constant
    print 'Relaxation_time: '+str((Mobility/10000 * Effective_mass) / (Charge_constant))
    return (Mobility/10000 * Effective_mass) / (Charge_constant)
#==========================================
def Thermal_Velocity(KB,Temperature,Effective_mass_factor):#Velocity that is caused by the thermal.
    Effective_mass = (Effective_mass_factor) * Electron_mass
    print 'Thermal_Velocity: '+str(math.pow( (3*KB*Temperature)/Effective_mass ,0.5))
    return math.pow( (3*KB*Temperature)/Effective_mass ,0.5)
#==========================================#==========================================#========================
def Mean_free_path(Thermal_Velocity,Relaxation_time):#The length that the particle will bump into another electrons for the first time.
    print 'Mean_free_path: '+str((Thermal_Velocity)*(Relaxation_time))
    return (Thermal_Velocity)*(Relaxation_time)
#==============================================================Gain=================================================================
def Gain(Avalanche_region,Mean_free_path):
    print 'Avalanche_region/Mean_free_path: '+str(Avalanche_region/Mean_free_path)
    return math.pow(2,Avalanche_region/Mean_free_path)
#==========================================#==========================================#========================
#Try
Eff_mass_factor = Effective_mass_factor(0,Temperature)
Mob = Mobility(A_Constant(Drift_Velocity_at_10V(0,Temperature)),Electric_field,Temperature)
MFP=Mean_free_path(Thermal_Velocity(KB,Temperature,Eff_mass_factor),Relaxation_time(Mob, Eff_mass_factor, CC))
print 'Gain: '+str(Gain(Avalanche_region(Temperature),MFP))
#============================================================Studies at there!!================================================================
#(1) The same E => What's G?
xarray=array("f",[])
yarray=array("f",[])

#for Temperature1 in range(4,78):
#  Eff_mass_factor = Effective_mass_factor(0,Temperature1)
#  Mob = Mobility(A_Constant(Drift_Velocity_at_10V(0,Temperature1)),Electric_field)
#  MFP=Mean_free_path(Thermal_Velocity(KB,Temperature1,Eff_mass_factor),Relaxation_time(Mob, Eff_mass_factor, CC))
#  print 'Gain_fix_Electric_field: '+str(Gain(Avalanche_region(Temperature1),MFP))

#for Electric_field1 in range(1,10):
#   Eff_mass_factor = Effective_mass_factor(0,Temperature)
#   Mob = Mobility(A_Constant(Drift_Velocity_at_10V(0,Temperature)),Electric_field1*10)
#    MFP=Mean_free_path(Thermal_Velocity(KB,Temperature,Eff_mass_factor),Relaxation_time(Mob, Eff_mass_factor, CC))
#    print 'Gain_fix_temperature: '+str(Gain(Avalanche_region(Temperature),MFP))










































