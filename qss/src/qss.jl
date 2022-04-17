module qss


    using Plots;gr()
    # list of public (API) to the user, not between files as those are linked as if in one file
    export ModelSettings, QSS_simulator ,QSS_simulate, plotX


    include("common/ModelSettings.jl")

    include("common/QSS_DataBase.jl")
    include("QSSFamily/QSS_simulator.jl")
    include("QSSFamily/QSS_integrator.jl")
    include("QSSFamily/QSS_modelworkshop.jl")
    include("QSSFamily/Quantizer.jl")
    include("QSSFamily/QSS1.jl")
    include("QSSFamily/Scheduler.jl")
    include("common/SimUtils.jl")


end # module
