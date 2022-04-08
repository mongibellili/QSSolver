
mutable struct SimulationSettings    
        ft :: Float64   
        dQmin ::Float64    
        dQrel ::Float64
        order :: Int        
        #solver :: SD_Solver  

    function SimulationSettings(ft:: Float64, dQmin ::Float64, dQrel ::Float64,order :: Int )
        settings = new()
        settings.ft = ft
        settings.dQmin = dQmin
        settings.dQrel = dQrel
        settings.order=order
        #settings.solver = solver
        settings
    end
end
  #@enum SD_Solver "qss1" "qss2" "qss3" "liqss"


  