
function integrateState(j::Int,order::Int,elapsed::Float64,x:: Vector{Array{Float64}} )
    if order ==1
        newX=last(x[2*j-1])+ elapsed * last(x[2*j])
        push!(x[2*j-1],newX)
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
       println(length(x[(order+1)*i-1]))
        println(length(t[i]))
        display(plot!(t[i],x[(order+1)*i-1]))
    end
    readline()
end