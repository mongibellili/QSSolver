function f(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        createT(q[2], cache[1])
        #= none:1 =#
        for k = 2:1
            #= none:1 =#
            (cache[k]).coeffs .= 0.0
        end
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        subT(negateT(q[1], cache[2]), q[2], cache[1])
        #= none:1 =#
        for k = 2:2
            #= none:1 =#
            (cache[k]).coeffs .= 0.0
        end
        #= none:1 =#
        return nothing
    end
end
function zcf(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    return nothing
end
function eventf(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    return nothing
end
function getDependencies()
return qssspecialized.QSS_data{2, 0, 0}([1.0, 0.0], [0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], Float64[], 10, StaticArraysCore.SVector{2, Int64}[[0, 2], [1, 2]], [Int64[], Int64[]], StaticArraysCore.SVector{2, Int64}[], StaticArraysCore.SVector{0, Int64}[], StaticArraysCore.SVector{2, SymEngine.Basic}[])
end
