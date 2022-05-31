module qss

    using SymEngine
    using StaticArrays 
    #using TimerOutputs
    # list of public (API) to the user, not between files as those are linked as if in one file
    export ProblemSettings,QSS_Problem,QSS_Solve ,  qss1,qss2,qss3,liqss,liqss2,saveat
    export @odeProblem
    #Utils
    include("Utils/SimUtils.jl")
    
    
    #QSSFamily/common/
    include("QSSFamily/Common/ProblemSettings.jl")
    include("QSSFamily/Common/QSSProblem.jl")
    include("QSSFamily/Common/QSS_Problem.jl")
    include("QSSFamily/Common/QSS_data.jl")
    include("QSSFamily/Common/Scheduler.jl")
    include("QSSFamily/Common/QSS_modelworkshop.jl")

    #QSSFamily/EasyQSS/
    include("QSSFamily/EasyQSS/Easy_QSS1_integrator.jl")
    include("QSSFamily/EasyQSS/Easy_QSS1.jl")
    include("QSSFamily/EasyQSS/Easy_QSS2_integrator.jl")
    include("QSSFamily/EasyQSS/Easy_QSS2.jl")
    
   

#= 
    #QSSFamily/FullQSS/
  
    include("QSSFamily/FullQSS/QSS1_integrator.jl")
    include("QSSFamily/FullQSS/QSS1.jl")

    
 
      #QSSFamily/MediumQSS/
      
      include("QSSFamily/MediumQSS/Medium_QSS1_integrator.jl")
      include("QSSFamily/MediumQSS/Medium_QSS1.jl")
  
      #QSSFamily/MediumQSS/ No jacobian
      
      include("QSSFamily/MediumQSS-NoJac/NoJac-Medium_QSS1_integrator.jl")
      include("QSSFamily/MediumQSS-NoJac/NoJac-Medium_QSS1.jl") =#
    
      #QSSFamily/NormalQSS/ jacobian not provided

      include("QSSFamily/NormalQSS-JacNotProvided/Normal_QSS1_integrator.jl")
      include("QSSFamily/NormalQSS-JacNotProvided/Normal_QSS1.jl")


end # module
