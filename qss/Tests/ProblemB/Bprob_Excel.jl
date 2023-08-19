
using XLSX
using formalqss
#using formalqss
#using formalqssAB
using BenchmarkTools
#= using StaticArrays
using BSON =#
using Plots
include("/home/mongi/projects/formalqss/Tests/ProblemB/BProblem.jl")
#include("d://BProblem.jl") 
  function test11(ft::Float64,solver::Tuple{Val{solType}, Val{V}},absTol,relTol)where {V,solType} 
    

  solmliqss2=QSS_Solve(probB11(),solver,dQmin=absTol,dQrel=relTol,finalTime=ft)
  solmliqss2Interp=solInterpolated(solmliqss2,0.01)
    u1, u2 = -0.06805253887502617, -0.5248046039821167
    λ1, λ2 = -3.523632227874817, -0.32636777212518275
   
    c1, c2 = -4.487356192869495, 2.4873561928694947
    xp1, xp2 = 0.0, 0.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 


     ("B11",(er1+er2)/2,solmliqss2.totalSteps,solmliqss2.simulStepCount)
   
end
 
function test12(ft::Float64,solver::Tuple{Val{solType}, Val{V}},absTol,relTol)where {V,solType} 
  solmliqss2=QSS_Solve(probB12(),solver,dQmin=absTol,dQrel=relTol,finalTime=ft)
  solmliqss2Interp=solInterpolated(solmliqss2,0.01)
    u1, u2 = -0.06805253887502617, -0.5248046039821167
    λ1, λ2 = -3.523632227874817, -0.32636777212518275
   
    c1, c2 = -5.636348488851935, 2.6363484888519344
    xp1, xp2 = 1.1102230246251565e-16, 0.9999999999999999
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
   
    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
   

     ("B12",(er1+er2)/2,solmliqss2.totalSteps,solmliqss2.simulStepCount)
   
end

function test13(ft::Float64,solver::Tuple{Val{solType}, Val{V}},absTol,relTol)where {V,solType} 
    #BSON.@load "sysB/relaxedAE_sysB13_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    solmliqss2=QSS_Solve(probB13(),solver,dQmin=absTol,dQrel=relTol,finalTime=ft)
    solmliqss2Interp=solInterpolated(solmliqss2,0.01)
    u1, u2 = -0.06805253887502617, -0.5248046039821167
    λ1, λ2 = -3.523632227874817, -0.32636777212518275
   
    c1, c2 = -27.467202112518294, 5.467202112518295
    xp1, xp2 = 3.552713678800501e-15, 20.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
 
 
    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
   
     ("B13",(er1+er2)/2,solmliqss2.totalSteps,solmliqss2.simulStepCount)
   
end

function test14(ft::Float64,solver::Tuple{Val{solType}, Val{V}},absTol,relTol)where {V,solType} 
    #BSON.@load "sysB/relaxedAE_sysB14_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    solmliqss2=QSS_Solve(probB14(),solver,dQmin=absTol,dQrel=relTol,finalTime=ft)
    solmliqss2Interp=solInterpolated(solmliqss2,0.01)
    u1, u2 = -0.06805253887502617, -0.5248046039821167
    λ1, λ2 = -3.523632227874817, -0.32636777212518275
   
    c1, c2 = -4.498846115829319, 2.488846115829319
    xp1, xp2 = 0.0, 0.009999999999999998
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
 
   
    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
  

     ("B14",(er1+er2)/2,solmliqss2.totalSteps,solmliqss2.simulStepCount)
end

function test15(ft::Float64,solver::Tuple{Val{solType}, Val{V}},absTol,relTol)where {V,solType} 
    #BSON.@load "sysB/relaxedAE_sysB15_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    solmliqss2=QSS_Solve(probB15(),solver,dQmin=absTol,dQrel=relTol,finalTime=ft)
    solmliqss2Interp=solInterpolated(solmliqss2,0.01)
    u1, u2 = -0.06805253887502617, -0.5248046039821167
    λ1, λ2 = -3.523632227874817, -0.32636777212518275
   
    c1, c2 = -1.2576052870427024, -1.7423947129572968
    xp1, xp2 = -2.0, 0.9999999999999993
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
 
    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 


     ("B15",(er1+er2)/2,solmliqss2.totalSteps,solmliqss2.simulStepCount)
   
end


function test21(ft::Float64,solver::Tuple{Val{solType}, Val{V}},absTol,relTol)where {V,solType} 
   # BSON.@load "sysB/relaxedAE_sysB21_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    solmliqss2=QSS_Solve(probB21(),solver,dQmin=absTol,dQrel=relTol,finalTime=ft)
    solmliqss2Interp=solInterpolated(solmliqss2,0.01)
    u1, u2 = -0.03593893348744506, -1.739061066512555
    λ1, λ2 = -6.856244266050219, -0.04375573394978005
   
    c1, c2 = -2.6293605409677845, 0.6293605409677845
    xp1, xp2 = 0.0, 0.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2

    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 


     ("B21",(er1+er2)/2,solmliqss2.totalSteps,solmliqss2.simulStepCount)
   
end

function test22(ft::Float64,solver::Tuple{Val{solType}, Val{V}},absTol,relTol)where {V,solType} 
    #BSON.@load "sysB/relaxedAE_sysB22_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    solmliqss2=QSS_Solve(probB22(),solver,dQmin=absTol,dQrel=relTol,finalTime=ft)
   
    solmliqss2Interp=solInterpolated(solmliqss2,0.01)
    u1, u2 = -0.03593893348744506, -1.739061066512555
    λ1, λ2 = -6.856244266050219, -0.04375573394978005
   
    c1, c2 = -3.6504623356016257, 0.6504623356016255
    xp1, xp2 = 0.0, 1.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
 

    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 


     ("B22",(er1+er2)/2,solmliqss2.totalSteps,solmliqss2.simulStepCount)
   
end

function test23(ft::Float64,solver::Tuple{Val{solType}, Val{V}},absTol,relTol)where {V,solType} 
   # BSON.@load "sysB/relaxedAE_sysB23_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    solmliqss2=QSS_Solve(probB23(),solver,dQmin=absTol,dQrel=relTol,finalTime=ft)
    solmliqss2Interp=solInterpolated(solmliqss2,0.01)
    u1, u2 = -0.03593893348744506, -1.739061066512555
    λ1, λ2 = -6.856244266050219, -0.04375573394978005
   
    c1, c2 = -23.0513964336446, 1.0513964336445958
    xp1, xp2 = -1.4210854715202004e-14, 20.000000000000007
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
 
    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 

     ("B23",(er1+er2)/2,solmliqss2.totalSteps,solmliqss2.simulStepCount)
   
end

function test24(ft::Float64,solver::Tuple{Val{solType}, Val{V}},absTol,relTol)where {V,solType} 
    #BSON.@load "sysB/relaxedAE_sysB24_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    solmliqss2=QSS_Solve(probB24(),solver,dQmin=absTol,dQrel=relTol,finalTime=ft)
    solmliqss2Interp=solInterpolated(solmliqss2,0.01)
    u1, u2 = -0.03593893348744506, -1.739061066512555
    λ1, λ2 = -6.856244266050219, -0.04375573394978005
   
    c1, c2 = -2.6395715589141227, 0.6295715589141229
    xp1, xp2 = 0.0, 0.009999999999999998
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
 

    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
 

     ("B24",(er1+er2)/2,solmliqss2.totalSteps,solmliqss2.simulStepCount)
   
end

function test25(ft::Float64,solver::Tuple{Val{solType}, Val{V}},absTol,relTol)where {V,solType} 
    #BSON.@load "sysB/relaxedAE_sysB25_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    solmliqss2=QSS_Solve(probB25(),solver,dQmin=absTol,dQrel=relTol,finalTime=ft)
    solmliqss2Interp=solInterpolated(solmliqss2,0.01)
    u1, u2 = -0.03593893348744506, -1.739061066512555
    λ1, λ2 = -6.856244266050219, -0.04375573394978005
   
    c1, c2 = -2.4761484322014207, -0.5238515677985796
    xp1, xp2 = -2.0, 1.0000000000000002
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2

    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
 

     ("B25",(er1+er2)/2,solmliqss2.totalSteps,solmliqss2.simulStepCount)
   
end

function test31(ft::Float64,solver::Tuple{Val{solType}, Val{V}},absTol,relTol)where {V,solType} 
   # BSON.@load "sysB/relaxedAE_sysB31_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    solmliqss2=QSS_Solve(probB31(),solver,dQmin=absTol,dQrel=relTol,finalTime=ft)
    solmliqss2Interp=solInterpolated(solmliqss2,0.01)
    u1, u2 = -0.02252283136759291, -1.6649771686324073
    λ1, λ2 = -6.909908674529628, -0.3400913254703717
   
    c1, c2 = -2.6362707559197682, 0.6362707559197684
    xp1, xp2 = 0.0, 0.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
  
    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
 
     ("B31",(er1+er2)/2,solmliqss2.totalSteps,solmliqss2.simulStepCount)
end

function test32(ft::Float64,solver::Tuple{Val{solType}, Val{V}},absTol,relTol)where {V,solType} 
    #BSON.@load "sysB/relaxedAE_sysB32_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    solmliqss2=QSS_Solve(probB32(),solver,dQmin=absTol,dQrel=relTol,finalTime=ft)
    solmliqss2Interp=solInterpolated(solmliqss2,0.01)
    u1, u2 = -0.02252283136759291, -1.6649771686324073
    λ1, λ2 = -6.909908674529628, -0.3400913254703717
   
    c1, c2 = -3.6499836676620214, 0.6499836676620214
    xp1, xp2 = 1.1102230246251565e-16, 1.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
   
    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
  

     ("B32",(er1+er2)/2,solmliqss2.totalSteps,solmliqss2.simulStepCount)
   
end

function test33(ft::Float64,solver::Tuple{Val{solType}, Val{V}},absTol,relTol)where {V,solType} 
   # BSON.@load "sysB/relaxedAE_sysB33_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    solmliqss2=QSS_Solve(probB33(),solver,dQmin=absTol,dQrel=relTol,finalTime=ft)
    solmliqss2Interp=solInterpolated(solmliqss2,0.01)
    u1, u2 = -0.02252283136759291, -1.6649771686324073
    λ1, λ2 = -6.909908674529628, -0.3400913254703717
   
    c1, c2 = -22.910528990764828, 0.910528990764828
    xp1, xp2 = 1.7763568394002505e-15, 20.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
  
    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
 

     ("B33",(er1+er2)/2,solmliqss2.totalSteps,solmliqss2.simulStepCount)
   
end

function test34(ft::Float64,solver::Tuple{Val{solType}, Val{V}},absTol,relTol)where {V,solType} 
   # BSON.@load "sysB/relaxedAE_sysB34_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    solmliqss2=QSS_Solve(probB34(),solver,dQmin=absTol,dQrel=relTol,finalTime=ft)
    solmliqss2Interp=solInterpolated(solmliqss2,0.01)
    u1, u2 = -0.02252283136759291, -1.6649771686324073
    λ1, λ2 = -6.909908674529628, -0.3400913254703717
   
    c1, c2 = -2.646407885037191, 0.636407885037191
    xp1, xp2 = 8.673617379884035e-19, 0.01
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2

    
    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
  

     ("B34",(er1+er2)/2,solmliqss2.totalSteps,solmliqss2.simulStepCount)
   
end

function test35(ft::Float64,solver::Tuple{Val{solType}, Val{V}},absTol,relTol)where {V,solType} 
    #BSON.@load "sysB/relaxedAE_sysB35_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    solmliqss2=QSS_Solve(probB35(),solver,dQmin=absTol,dQrel=relTol,finalTime=ft)
    solmliqss2Interp=solInterpolated(solmliqss2,0.01)
    u1, u2 = -0.02252283136759291, -1.6649771686324073
    λ1, λ2 = -6.909908674529628, -0.3400913254703717
   
    c1, c2 = -2.432293802791496, -0.5677061972085039
    xp1, xp2 = -2.0, 1.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2

    #= saveAbsoluteError(solmliqss2Interp,1,x1)  
    saveAbsoluteError(solmliqss2Interp,2,x2)  =#
    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
   #=  er1a=getIntervalError(solmliqss2Interp,1,x1,0.0,20.0)  
    er2a=getIntervalError(solmliqss2Interp,2,x2,0.0,20.0)  
    er1b=getIntervalError(solmliqss2Interp,1,x1,20.0,1000.0)  
    er2b=getIntervalError(solmliqss2Interp,2,x2,20.0,1000.0)  =# 

     ("B35",(er1+er2)/2,solmliqss2.totalSteps,solmliqss2.simulStepCount)
   
end


function test41(ft::Float64,solver::Tuple{Val{solType}, Val{V}},absTol,relTol)where {V,solType} 
    #BSON.@load "sysB/relaxedAE_sysB41_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    solmliqss2=QSS_Solve(probB41(),solver,dQmin=absTol,dQrel=relTol,finalTime=ft)
    solmliqss2Interp=solInterpolated(solmliqss2,0.01)
    u1, u2 = -0.1144790629155213, -0.13539593708447864
    λ1, λ2 = -10.841674966758294, -9.168325033241706
   
    c1, c2 = -60.754387290570186, 58.754387290570186
    xp1, xp2 = 0.0, 0.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2

   #=  saveAbsoluteError(solmliqss2Interp,1,x1)  
    saveAbsoluteError(solmliqss2Interp,2,x2)  =#

    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
   #=  er1a=getIntervalError(solmliqss2Interp,1,x1,0.0,20.0)  
    er2a=getIntervalError(solmliqss2Interp,2,x2,0.0,20.0)  
    er1b=getIntervalError(solmliqss2Interp,1,x1,20.0,1000.0)  
    er2b=getIntervalError(solmliqss2Interp,2,x2,20.0,1000.0)  
 =#
     ("B41",(er1+er2)/2,solmliqss2.totalSteps,solmliqss2.simulStepCount)
   
end

function test42(ft::Float64,solver::Tuple{Val{solType}, Val{V}},absTol,relTol)where {V,solType} 
   # BSON.@load "sysB/relaxedAE_sysB42_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    solmliqss2=QSS_Solve(probB42(),solver,dQmin=absTol,dQrel=relTol,finalTime=ft)
    solmliqss2Interp=solInterpolated(solmliqss2,0.01)
    u1, u2 = -0.1144790629155213, -0.13539593708447864
    λ1, λ2 = -10.841674966758294, -9.168325033241706
   
    c1, c2 = -67.22743560509412, 64.22743560509412
    xp1, xp2 = 0.0, 1.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2


    #= saveAbsoluteError(solmliqss2Interp,1,x1)  
    saveAbsoluteError(solmliqss2Interp,2,x2)  =#
    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
   #=  er1a=getIntervalError(solmliqss2Interp,1,x1,0.0,20.0)  
    er2a=getIntervalError(solmliqss2Interp,2,x2,0.0,20.0)  
    er1b=getIntervalError(solmliqss2Interp,1,x1,20.0,1000.0)  
    er2b=getIntervalError(solmliqss2Interp,2,x2,20.0,1000.0)   =#

     ("B42",(er1+er2)/2,solmliqss2.totalSteps,solmliqss2.simulStepCount)
   
end


function test43(ft::Float64,solver::Tuple{Val{solType}, Val{V}},absTol,relTol)where {V,solType} 
    #BSON.@load "sysB/relaxedAE_sysB43_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    solmliqss2=QSS_Solve(probB43(),solver,dQmin=absTol,dQrel=relTol,finalTime=ft)

    solmliqss2Interp=solInterpolated(solmliqss2,0.01)
    u1, u2 = -0.1144790629155213, -0.13539593708447864
    λ1, λ2 = -10.841674966758294, -9.168325033241706
   
    c1, c2 = -190.21535358104902, 168.21535358104902
    xp1, xp2 = 0.0, 20.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
 
  

  #=   saveAbsoluteError(solmliqss2Interp,1,x1)  
    saveAbsoluteError(solmliqss2Interp,2,x2)  =#
    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
   #=  er1a=getIntervalError(solmliqss2Interp,1,x1,0.0,20.0)  
    er2a=getIntervalError(solmliqss2Interp,2,x2,0.0,20.0)  
    er1b=getIntervalError(solmliqss2Interp,1,x1,20.0,1000.0)  
    er2b=getIntervalError(solmliqss2Interp,2,x2,20.0,1000.0)   =#

    ("B43",(er1+er2)/2,solmliqss2.totalSteps,solmliqss2.simulStepCount)
   
end

 function test44(ft::Float64,solver::Tuple{Val{solType}, Val{V}},absTol,relTol)where {V,solType} 
    #BSON.@load "sysB/relaxedAE_sysB44_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    solmliqss2=QSS_Solve(probB44(),solver,dQmin=absTol,dQrel=relTol,finalTime=ft)
    solmliqss2Interp=solInterpolated(solmliqss2,0.01)
    u1, u2 = -0.1144790629155213, -0.13539593708447864
    λ1, λ2 = -10.841674966758294, -9.168325033241706
   
    c1, c2 = -60.81911777371541, 58.809117773715414
    xp1, xp2 = 0.0, 0.01
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2


  #=   saveAbsoluteError(solmliqss2Interp,1,x1)  
    saveAbsoluteError(solmliqss2Interp,2,x2)  =#
    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
  #=   er1a=getIntervalError(solmliqss2Interp,1,x1,0.0,20.0)  
    er2a=getIntervalError(solmliqss2Interp,2,x2,0.0,20.0)  
    er1b=getIntervalError(solmliqss2Interp,1,x1,20.0,1000.0)  
    er2b=getIntervalError(solmliqss2Interp,2,x2,20.0,1000.0)  
   =#
     ("B44",(er1+er2)/2,solmliqss2.totalSteps,solmliqss2.simulStepCount)
    
   
end
 
 function test45(ft::Float64,solver::Tuple{Val{solType}, Val{V}},absTol,relTol)where {V,solType} 
    #BSON.@load "sysB/relaxedAE_sysB45_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    solmliqss2=QSS_Solve(probB45(),solver,dQmin=absTol,dQrel=relTol,finalTime=ft)
    solmliqss2Interp=solInterpolated(solmliqss2,0.01)
    u1, u2 = -0.1144790629155213, -0.13539593708447864
    λ1, λ2 = -10.841674966758294, -9.168325033241706
   
    c1, c2 = 28.21275861611204, -31.18056545715832
    xp1, xp2 = -2.0, 1.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2

    


  #=   saveAbsoluteError(solmliqss2Interp,1,x1)  
    saveAbsoluteError(solmliqss2Interp,2,x2)  =#
    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
   #=  er1a=getIntervalError(solmliqss2Interp,1,x1,0.0,20.0)  
    er2a=getIntervalError(solmliqss2Interp,2,x2,0.0,20.0)  
    er1b=getIntervalError(solmliqss2Interp,1,x1,20.0,1000.0)  
    er2b=getIntervalError(solmliqss2Interp,2,x2,20.0,1000.0)   =#

     ("B45",(er1+er2)/2,solmliqss2.totalSteps,solmliqss2.simulStepCount)
end
function test51(ft::Float64,solver::Tuple{Val{solType}, Val{V}},absTol,relTol)where {V,solType} 
   # BSON.@load "sysB/relaxedAE_sysB51_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    solmliqss2=QSS_Solve(probB51(),solver,dQmin=absTol,dQrel=relTol,finalTime=ft)
    solmliqss2Interp=solInterpolated(solmliqss2,0.01)
    u1, u2 = -8.73522174738572, -7.385745994549763
    λ1, λ2 = -10.841674966758294, -9.168325033241706
   
    c1, c2 = 11.68712513430148, -13.68712513430148
    xp1, xp2 = 0.0, 0.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2

    #= saveAbsoluteError(solmliqss2Interp,1,x1)  
    saveAbsoluteError(solmliqss2Interp,2,x2)  =#

    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
   #=  er1a=getIntervalError(solmliqss2Interp,1,x1,0.0,20.0)  
    er2a=getIntervalError(solmliqss2Interp,2,x2,0.0,20.0)  
    er1b=getIntervalError(solmliqss2Interp,1,x1,20.0,1000.0)  
    er2b=getIntervalError(solmliqss2Interp,2,x2,20.0,1000.0)  
 =#
     ("B51",(er1+er2)/2,solmliqss2.totalSteps,solmliqss2.simulStepCount)
   
   
end
 function test52(ft::Float64,solver::Tuple{Val{solType}, Val{V}},absTol,relTol)where {V,solType} 
   # BSON.@load "sysB/relaxedAE_sysB52_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    solmliqss2=QSS_Solve(probB52(),solver,dQmin=absTol,dQrel=relTol,finalTime=ft)
    solmliqss2Interp=solInterpolated(solmliqss2,0.01)
    u1, u2 = -8.73522174738572, -7.385745994549763
    λ1, λ2 = -10.841674966758294, -9.168325033241706
   
    c1, c2 = 17.16017344882542, -20.16017344882542
    xp1, xp2 = 0.0, 1.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2


  #=   saveAbsoluteError(solmliqss2Interp,1,x1)  
    saveAbsoluteError(solmliqss2Interp,2,x2)  =#
    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
   #=  er1a=getIntervalError(solmliqss2Interp,1,x1,0.0,20.0)  
    er2a=getIntervalError(solmliqss2Interp,2,x2,0.0,20.0)  
    er1b=getIntervalError(solmliqss2Interp,1,x1,20.0,1000.0)  
    er2b=getIntervalError(solmliqss2Interp,2,x2,20.0,1000.0)   =#

     ("B52",(er1+er2)/2,solmliqss2.totalSteps,solmliqss2.simulStepCount)
   
end
function test53(ft::Float64,solver::Tuple{Val{solType}, Val{V}},absTol,relTol)where {V,solType} 
   # BSON.@load "sysB/relaxedAE_sysB53_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    solmliqss2=QSS_Solve(probB53(),solver,dQmin=absTol,dQrel=relTol,finalTime=ft)
    solmliqss2Interp=solInterpolated(solmliqss2,0.01)
    u1, u2 = -8.73522174738572, -7.385745994549763
    λ1, λ2 = -10.841674966758294, -9.168325033241706
   
    c1, c2 = 121.14809142478035, -143.14809142478035
    xp1, xp2 = 0.0, 20.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
 
    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
     ("B53",(er1+er2)/2,solmliqss2.totalSteps,solmliqss2.simulStepCount)
   
end
function test54(ft::Float64,solver::Tuple{Val{solType}, Val{V}},absTol,relTol)where {V,solType} 
    #BSON.@load "sysB/relaxedAE_sysB54_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    solmliqss2=QSS_Solve(probB54(),solver,dQmin=absTol,dQrel=relTol,finalTime=ft)
    solmliqss2Interp=solInterpolated(solmliqss2,0.01)
    u1, u2 = -8.73522174738572, -7.385745994549763
    λ1, λ2 = -10.841674966758294, -9.168325033241706
    c1, c2 = 11.74185561744672, -13.75185561744672
    xp1, xp2 = 0.0, 0.01
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
    ("B54",(er1+er2)/2,solmliqss2.totalSteps,solmliqss2.simulStepCount)
   
end
function test55(ft::Float64,solver::Tuple{Val{solType}, Val{V}},absTol,relTol)where {V,solType} 
    #BSON.@load "sysB/relaxedAE_sysB55_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
    solmliqss2=QSS_Solve(probB55(),solver,dQmin=absTol,dQrel=relTol,finalTime=ft)
    solmliqss2Interp=solInterpolated(solmliqss2,0.01)
    u1, u2 = -8.73522174738572, -7.385745994549763
    λ1, λ2 = -10.841674966758294, -9.168325033241706
   
    c1, c2 = 15.67811643831823, -18.67811643831823
    xp1, xp2 = -2.0, 1.0
    x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
    x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
    er1=getError(solmliqss2Interp,1,x1)  
    er2=getError(solmliqss2Interp,2,x2) 
      ("B55",(er1+er2)/2,solmliqss2.totalSteps,solmliqss2.simulStepCount)
   
end 



function nmliqsstest()
     absTol=1e-6
     relTol=1e-3
     ft=10.0
     solver=nmliqss2()
          
          #= 
          result11=test11(ft,solver,absTol,relTol)
          result12=test12(ft,solver,absTol,relTol)
          result13=test13(ft,solver,absTol,relTol)
          result14=test14(ft,solver,absTol,relTol)
          result15=test15(ft,solver,absTol,relTol)
          result21=test21(ft,solver,absTol,relTol)
          result22=test22(ft,solver,absTol,relTol)
          result23=test23(ft,solver,absTol,relTol)
          result24=test24(ft,solver,absTol,relTol)
          result25=test25(ft,solver,absTol,relTol)
          result31=test31(ft,solver,absTol,relTol)
          result32=test32(ft,solver,absTol,relTol)
          result33=test33(ft,solver,absTol,relTol)
          result34=test34(ft,solver,absTol,relTol)
          result35=test35(ft,solver,absTol,relTol)
          result41=test41(ft,solver,absTol,relTol)
          result42=test42(ft,solver,absTol,relTol) =#
          result43=test43(ft,solver,absTol,relTol)
          @show result43
         #=  result44=test44(ft,solver,absTol,relTol) 
          result45=test45(ft,solver,absTol,relTol) 
          result51=test51(ft,solver,absTol,relTol)
          result52=test52(ft,solver,absTol,relTol)
          result53=test53(ft,solver,absTol,relTol)
          result54=test54(ft,solver,absTol,relTol) 
          result55=test55(ft,solver,absTol,relTol) 


          XLSX.openxlsx("nmLiqss2______.xlsx", mode="w") do xf
          sheet = xf[1]
          sheet["A1"] = "nmLiqss2_"
          sheet["A4"] = collect(("Sys","error","totalSteps","simul_steps"))
          sheet["A5"] = collect(result11)
          sheet["A6"] = collect(result12)
          sheet["A7"] = collect(result13)
          sheet["A8"] = collect(result14)
          sheet["A9"] = collect(result15)
          sheet["A10"] = collect(result21)
          sheet["A11"] = collect(result22)
          sheet["A12"] = collect(result23)
          sheet["A13"] = collect(result24)
          sheet["A14"] = collect(result25)

          sheet["A15"] = collect(result31)
          sheet["A16"] = collect(result32)
          sheet["A17"] = collect(result33)
          sheet["A18"] = collect(result34)
          sheet["A19"] = collect(result35)

          sheet["A20"] = collect(result41)
          sheet["A21"] = collect(result42)
          sheet["A22"] = collect(result43)
          sheet["A23"] = collect(result44) 
          sheet["A24"] = collect(result45)
          sheet["A25"] = collect(result51)
          sheet["A26"] = collect(result52)
          sheet["A27"] = collect(result53)
          sheet["A28"] = collect(result54)
          sheet["A29"] = collect(result55)
          end   =#
end
nmliqsstest() 
############################################"

#= 
function test53only(ft::Float64,solver::Tuple{Val{solType}, Val{V}},absTol,relTol)where {V,solType} 
     # BSON.@load "sysB/relaxedAE_sysB53_mliqss2_errorE-3_step0.1_ft1000_Interp.bson" solmliqss2Interp
     solmliqss2=QSS_Solve(probB53(),solver,dQmin=absTol,dQrel=relTol,finalTime=ft)
     solmliqss2Interp=solInterpolated(solmliqss2,0.01)
     save_Sol(solmliqss2,"x1",1;xlims=(30.0,40.0),ylims=(-0.05,0.05))
     u1, u2 = -8.73522174738572, -7.385745994549763
     λ1, λ2 = -10.841674966758294, -9.168325033241706
     
     c1, c2 = 121.14809142478035, -143.14809142478035
     xp1, xp2 = 0.0, 20.0
     x1(t)=c1*u1*exp(λ1*t)+c2*u2*exp(λ2*t)+xp1
     x2(t)=c1*exp(λ1*t)+c2*exp(λ2*t)+xp2
     
     er1=getError(solmliqss2Interp,1,x1)  
     er2=getError(solmliqss2Interp,2,x2) 
          ("$(solmliqss2.algName)",(er1+er2)/2,solmliqss2.totalSteps,solmliqss2.simulStepCount)
  end

  
function globaltest()
          function liqsstestOneSys()
               ft=1000.0
               solver=liqss2()
               test53only(ft,solver,absTol,relTol)
          end
          #@time
          resultliqss2=liqsstestOneSys()
          function mliqsstestOneSys()
               ft=1000.0
               solver=mliqss2()
               test53only(ft,solver,absTol,relTol)
          end
          resultmliqss2=mliqsstestOneSys()
          function nliqsstestOneSys()
               ft=1000.0
               solver=nliqss2()
               test53only(ft,solver,absTol,relTol) 
          end
          #@time
          resultnliqss2=nliqsstestOneSys()
          function nmliqsstestOneSys()
               ft=1000.0
               solver=nmliqss2()
               test53only(ft,solver,absTol,relTol)
          end
          resultnmliqss2=nmliqsstestOneSys()


          XLSX.openxlsx("sys53 all solvers______.xlsx", mode="w") do xf
               sheet = xf[1]
               sheet["A1"] = "sys53 all solvers_"
               sheet["A4"] = collect(("solver","error","totalSteps","simul_steps"))
               sheet["A5"] = collect(resultliqss2)
               sheet["A6"] = collect(resultmliqss2)
               sheet["A7"] = collect(resultnliqss2)
               sheet["A8"] = collect(resultnmliqss2)
               
      
               end  

end

globaltest()
 =#


#=   
 function liqsstest()
     ft=1000.0
     solver=liqss2()
     test53(ft,solver,absTol,relTol)
     #= result22=test53(ft,solver,absTol,relTol)
     @show result22 =#
     result11=test11(ft,solver,absTol,relTol)
     result12=test12(ft,solver,absTol,relTol)
     result13=test13(ft,solver,absTol,relTol)
     result14=test14(ft,solver,absTol,relTol)
     result15=test15(ft,solver,absTol,relTol)
     result21=test21(ft,solver,absTol,relTol)
     result22=test22(ft,solver,absTol,relTol)
     result23=test23(ft,solver,absTol,relTol)
     result24=test24(ft,solver,absTol,relTol)
     result25=test25(ft,solver,absTol,relTol)
     result31=test31(ft,solver,absTol,relTol)
     result32=test32(ft,solver,absTol,relTol)
     result33=test33(ft,solver,absTol,relTol)
     result34=test34(ft,solver,absTol,relTol)
     result35=test35(ft,solver,absTol,relTol)
     result41=test41(ft,solver,absTol,relTol)
     result42=test42(ft,solver,absTol,relTol)
     result43=test43(ft,solver,absTol,relTol)
     result44=test44(ft,solver,absTol,relTol) 
     result45=test45(ft,solver,absTol,relTol) 
     result51=test51(ft,solver,absTol,relTol)
     result52=test52(ft,solver,absTol,relTol)
     result53=test53(ft,solver,absTol,relTol)
     result54=test54(ft,solver,absTol,relTol) 
     result55=test55(ft,solver,absTol,relTol) 


     XLSX.openxlsx("Liqss2______.xlsx", mode="w") do xf
     sheet = xf[1]
     sheet["A1"] = "Liqss2_"
     sheet["A4"] = collect(("Sys","error","totalSteps","simul_steps"))
     sheet["A5"] = collect(result11)
     sheet["A6"] = collect(result12)
     sheet["A7"] = collect(result13)
     sheet["A8"] = collect(result14)
     sheet["A9"] = collect(result15)
     sheet["A10"] = collect(result21)
     sheet["A11"] = collect(result22)
     sheet["A12"] = collect(result23)
     sheet["A13"] = collect(result24)
     sheet["A14"] = collect(result25)

     sheet["A15"] = collect(result31)
     sheet["A16"] = collect(result32)
     sheet["A17"] = collect(result33)
     sheet["A18"] = collect(result34)
     sheet["A19"] = collect(result35)

     sheet["A20"] = collect(result41)
     sheet["A21"] = collect(result42)
     sheet["A22"] = collect(result43)
     sheet["A23"] = collect(result44) 
     sheet["A24"] = collect(result45)
     sheet["A25"] = collect(result51)
     sheet["A26"] = collect(result52)
     sheet["A27"] = collect(result53)
     sheet["A28"] = collect(result54)
     sheet["A29"] = collect(result55)
     end  
end
 =#