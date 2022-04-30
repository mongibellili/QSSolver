module qss


    using StaticArrays 
    #using TimerOutputs
    # list of public (API) to the user, not between files as those are linked as if in one file
    export ModelSettings,QSS_simulate,QSS_simGenerate, QSS_integrate , plotX, qss1,qss2,qss3,saveat,liqss,liqss2
    
    #export @make_simulator, QSS_simulator, ModelSettings
    #export computeStates
#=     include("QSSFamily/make-sim-set.jl")
    include("QSSFamily/create-sim-set.jl") =#
    
    include("common/ModelSettings.jl")
    include("QSSFamily/QSS_modelworkshop.jl")
    include("QSSFamily/QSS_simulator.jl")
    include("QSSFamily/QSS_integrator.jl")
    include("QSSFamily/QSS1.jl")
    include("QSSFamily/Scheduler.jl")
    include("common/SimUtils.jl")
    


#include("common/ModelSettings.jl")
      # include("common/QSS_DataBase.jl")
   # include("QSSFamily/QSS_simulator.jl")
       #include("QSSFamily/Quantizer.jl")
    


end # module
