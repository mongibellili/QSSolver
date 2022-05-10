module qss


    using StaticArrays 
    #using TimerOutputs
    # list of public (API) to the user, not between files as those are linked as if in one file
    export ProblemSettings,QSS_Problem,QSS_Solve ,  qss1,qss2,qss3,liqss,liqss2,saveat
   
    #Utils
    include("Utils/SimUtils.jl")
    include("Utils/ProblemSettings.jl")
    
    #QSSFamily/common/
    include("QSSFamily/Common/QSS_Problem.jl")
    include("QSSFamily/Common/QSS_simulator.jl")
    include("QSSFamily/Common/Scheduler.jl")


    #QSSFamily/EasyQSS/
    include("QSSFamily/EasyQSS/Easy_QSS1_integrator.jl")
    include("QSSFamily/EasyQSS/Easy_QSS1.jl")
    include("QSSFamily/EasyQSS/Easy_QSS2_integrator.jl")
    include("QSSFamily/EasyQSS/Easy_QSS2.jl")
    include("QSSFamily/EasyQSS/Easy_QSS_modelworkshop.jl")


    #QSSFamily/FullQSS/
    include("QSSFamily/FullQSS/QSS_modelworkshop.jl")
    include("QSSFamily/FullQSS/QSS1_integrator.jl")
    include("QSSFamily/FullQSS/QSS1.jl")

    
 
    


#include("common/ModelSettings.jl")
      # include("common/QSS_DataBase.jl")
   # include("QSSFamily/QSS_simulator.jl")
       #include("QSSFamily/Quantizer.jl")
    


end # module
