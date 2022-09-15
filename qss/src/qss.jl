module qss
using RuntimeGeneratedFunctions
using StaticArrays
using SymEngine
using ExprTools
using MacroTools: postwalk, prewalk, @capture
#using TaylorSeries
import Base.:-
import Base.:+
import Base.:*
#use_show_default(true)
#= using LinearAlgebra: norm
import LinearAlgebra: norm =#

RuntimeGeneratedFunctions.init(@__MODULE__)


#this section belongs to taylorseries subcomponent
import Base: ==, +, -, *, /, ^

import Base: iterate, size, eachindex, firstindex, lastindex,
    eltype, length, getindex, setindex!, axes, copyto!

import Base:  sqrt, exp, log, sin, cos, sincos, tan,
    asin, acos, atan, sinh, cosh, tanh, atanh, asinh, acosh,
    zero, one, zeros, ones, isinf, isnan, iszero,
    convert, promote_rule, promote, show,
    real, imag, conj, adjoint,
    rem, mod, mod2pi, abs, abs2,
   
    power_by_squaring,
    rtoldefault, isfinite, isapprox, rad2deg, deg2rad
    #end of taylorseries section of imports


    # list of public (API) to the user, not between files as those are linked as if in one file
    export SimSettings,QSS_Problem,QSS_Solve ,  qss1,qss2,qss3,liqss,liqss2,saveat
   #=  ,Taylor0,AbstractSeries,resize_coeffs1!,constant_term  ,linear_polynomial, nonlinear_polynomial,
    normalize_taylor,taylor_expand,update! =#
    export  @NLodeProblem,save_prob_to_model,QSS_Solve_from_model
   # export mimicqss,@stub1,stub2,f1
    export Taylor0,mulT,createT,addsub,negateT,subsub,subadd,subT,addT,muladdT,mulsub,divT
#export section of ts subcomponent # to be deleted later after testing
#= export getcoeff, derivative, integrate, differentiate,#integrate!,#differentiate!,evaluateX,
    evaluate, evaluate!, inverse, set_taylor1_varname,
      displayBigO, use_show_default,
    get_order =#
     #end ts exports

       #include section of ts subcomponent
       include("ownTaylor/parameters.jl")  
       include("ownTaylor/constructors.jl") 
       include("ownTaylor/conversion.jl")        
       include("ownTaylor/arithmetic.jl")
       include("ownTaylor/arithmeticT.jl")
       include("ownTaylor/functions.jl")
       include("ownTaylor/functionsT.jl")
       include("ownTaylor/power.jl")
       include("ownTaylor/powerT.jl")
       include("ownTaylor/auxiliary.jl") 
       include("ownTaylor/calculus.jl")    
       include("ownTaylor/other_functions.jl")
       include("ownTaylor/evaluate.jl")
       
       
       #end ts includes
    #Utils
    include("Utils/SimUtils.jl") 
    #QSSFamily/common/
   # include("Common/SimSettings.jl")
    include("Common/QSSNLProblemHelper.jl")
    include("Common/QSSNLProblem.jl")
   # include("Common/QSSNLProblem2.jl")
    include("Common/QSS_data.jl")
    include("Common/Scheduler.jl")
    include("Common/QSS_modelworkshop.jl")
#= 
    #normal: only f parsed
   # include("NL/normal/NL_QSS2_Integrator.jl")
    include("NL/normal/NL_QSS2_Integratortimeit.jl")
    include("NL/normal/NL_QSS2.jl")
 =#
  #inclusive: f derf derderf parsed
  include("NL/model.jl")
       include("NL/normal/NL_QSS2_Integrator.jl")
       include("NL/normal/NL_LiQSS_Integrator.jl")
     #  include("NL/normal/NL_QSS_integrate_from_model.jl")
      # include("NL/normal/NL_LiQSS_integrate_from_model.jl")
       
       # include("NL/normal/NL_QSS2_Integrator0-timeit.jl")
        
      # include("NL/normal/NL_QSS2_Integrator-timeit.jl")
    #  include("NL/normal/NL_QSS2_Integrator-fbang.jl")
       include("NL/normal/NL_QSS2.jl")
#= include("../../diffequFunc.jl")
include("mimicqss.jl") =#

#include("ownTaylor/testOwnTaylor/testOperators.jl")

end # module
