
function integrateState(j::Int,order::Int,elapsed::Float64,x::  Vector{Float64} )
   0
    if order ==1
        x=(x[2*j-1])+ elapsed * (x[2*j])
       
       # push!(x[2*j-1],(x[2*j-1])+ elapsed * (x[2*j]))
    else
        println("keep coding next orders")
    end
    
end

function minPosRoot(coeff::Vector{Float64}, order ::Int)
    mpr=-1
    if order==1
        if coeff[2] == 0 
            mpr = Inf
        else 
            mpr = -coeff[1] / coeff[2];
            #println(mpr)
        end
        if mpr < 0
            mpr = Inf
        end
    end
    return mpr
end

function plotX(simulator :: QSS_simulator)
    x=simulator.qssData.x
    t=simulator.qssTime.tx
    states=simulator.settings.states
    order=simulator.settings.order

   for i = 1:states
      # println("steps to solve eq$i= ",length(x[(order+1)*i-1]))
        #println(length(t[i]))
        display(plot!(t[i],x[(order+1)*i-order]))#,xlims = (0.12,0.15),ylims = (1.23,1.26)))
        #, xlims = (0,0.00005),ylims = (1,1.0001)))
        #display(plot!(t[i],x[(order+1)*i])) #can't plot derivative cuz length does not equal time vector length
    end
    
    readline()
end