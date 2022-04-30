macro make_simulator(settingstates)
    esc(quote struct ModelSettings
            initConditions ::SVector{$settingstates,Float64} 
            jacobian::SMatrix{$settingstates,$settingstates,Float64}
            finalTime:: Float64   
            initialTime :: Float64    
            dQmin ::Float64    
            dQrel ::Float64  
            solver ::Val{T} where {T}
            savetimeincrement::Float64
            end
        struct QSS_simulator
            
           
            savedVars::SVector{$settingstates,Array{Float64}}
            savedTimes::Array{Float64}
            quantum :: MVector{$settingstates,Float64} 
            x :: MVector{$settingstates*(1+1),Float64}
            q :: MVector{$settingstates*(1+1),Float64} 
            tx ::  MVector{$settingstates,Float64} 
            tq :: MVector{$settingstates,Float64}
            nextStateTime :: MVector{$settingstates,Float64} 
            minTimeValue::MVector{1,Float64} 
            minIndex :: MVector{1,Int}
            states::Int
            order::Int
            end
    end)
end
