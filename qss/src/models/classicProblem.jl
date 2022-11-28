function f(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        createT(q[2], cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        subT(negateT(q[1], cache[2]), q[2], cache[1])
        #= none:1 =#
        return nothing
    end
end
function twoVarSys1(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        createT(q[2], cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        subT(negateT(q[1], cache[2]), q[2], cache[1])
        #= none:1 =#
        return nothing
    end
end

function twoVarSys12(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        mulT(0.01, q[2], cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        subadd(mulT(-100.0, q[1], cache[2]), mulT(100.0, q[2], cache[3]), 2020.0, cache[1])
        #= none:1 =#
        return nothing
    end
end
function twoVarSys13(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        subadd(negateT(q[1], cache[2]), q[2], 0.2, cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        subadd(q[1], q[2], 1.2, cache[1])
        #= none:1 =#
        return nothing
    end
end

function twoVarSys14(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        subadd(negateT(q[1], cache[2]), mulT(10.0, q[2], cache[3]), 0.2, cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        subadd(mulT(10.0, q[1], cache[2]), q[2], 1.2, cache[1])
        #= none:1 =#
        return nothing
    end
end
function twoVarSys16(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        subadd(mulT(-1.01, q[1], cache[2]), mulT(100.0, q[2], cache[3]), 1000.2, cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        subadd(mulT(10.1, q[1], cache[2]), mulT(1.001, q[2], cache[3]), 0.2, cache[1])
        #= none:1 =#
        return nothing
    end
end
function twoVarSys15(j::Int, q::Vector{Taylor0{Float64}}, d::Vector{Float64}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if j == 1
        #= none:1 =#
        addT(mulT(-4.0, q[1], cache[2]), mulT(0.5, q[2], cache[3]), 1.0, cache[1])
        #= none:1 =#
        return nothing
    else
        #= none:1 =#
        subadd(q[1], mulT(0.25, q[2], cache[2]), 0.25, cache[1])
        #= none:1 =#
        return nothing
    end
end
