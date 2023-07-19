function ICadvection1000V2(j::Int)
    if j == 0
        return nothing
    elseif 1 <= j <= 333
        return 1.0
    elseif 334 <= j <= 1000
        return 0.0
    end
end
function advection1000V2(i::Int, q::Vector{Taylor0{Float64}}, t::Taylor0{Float64}, cache::Vector{Taylor0{Float64}})
    if i == 0
        return nothing
    elseif i == 1000
        subadd(mulT(120.0, q[999], cache[2]), mulT(120.0, q[1000], cache[3]), mulTT(1000.0, q[1000], q[1000], subT(1.0, q[1000], cache[6]), cache[4], cache[5]), cache[1])
        return nothing
    elseif i == 1
        addT(mulT(-120.0, q[1], cache[2]), mulT(10.0, q[2], cache[3]), mulTT(1000.0, q[1], q[1], subT(1.0, q[1], cache[6]), cache[4], cache[5]), cache[1])
        return nothing
    elseif 2 <= i <= 999
        subadd(muladdT(10.0, q[i + 1], mulT(110.0, q[i - 1], cache[3]), cache[2]), mulT(120.0, q[i], cache[4]), mulTT(1000.0, q[i], q[i], subT(1.0, q[i], cache[7]), cache[5], cache[6]), cache[1])
        return nothing
    end
end
@resumable function jacadvection1000V2(i::Int)
    if i == 0
        return nothing
    elseif i == 1000
        #= none:3 =# @yield 1000
        #= none:4 =# @yield 999
    elseif i == 1
        #= none:6 =# @yield 2
        #= none:7 =# @yield 1
    elseif 2 <= i <= 999
        #= none:9 =# @yield i - 1
        #= none:10 =# @yield i + 1
        #= none:11 =# @yield i
    end
end
@resumable function SDadvection1000V2(i::Int)
    if 2 <= i - 1 <= 999
        #= none:1 =# @yield i - 1
    end
    if 2 <= i + 1 <= 999
        #= none:1 =# @yield i + 1
    end
    if 2 <= i <= 999
        #= none:1 =# @yield i
    end
    if i == 0
        return nothing
    elseif i == 999
        #= none:3 =# @yield 1000
    elseif i == 2
        #= none:5 =# @yield 1
    elseif i == 1000
        #= none:7 =# @yield 1000
    elseif i == 1
        #= none:9 =# @yield 1
    end
    return nothing
end
function mapadvection1000V2(cache::MVector{1, Int}, i::Int, j::Int)
    if i == 0
        return nothing
    elseif i == 1000
        if j == 0
            return nothing
        elseif j == 999
            cache[1] = 1
            return nothing
        elseif j == 1000
            cache[1] = 2
            return nothing
        end
    elseif i == 1
        if j == 0
            return nothing
        elseif j == 1
            cache[1] = 1
            return nothing
        elseif j == 2
            cache[1] = 2
            return nothing
        end
    elseif 2 <= i <= 999
        if j == 0
            return nothing
        elseif j == i - 1
            cache[1] = 1
            return nothing
        elseif j == i
            cache[1] = 2
            return nothing
        elseif j == i + 1
            cache[1] = 3
            return nothing
        end
    end
end
function jacDimensionadvection1000V2(i::Int)
    if i == 0
        return nothing
    elseif i == 1000
        return 2
    elseif i == 1
        return 2
    elseif 2 <= i <= 999
        return 3
    end
end
function prob_advection1000V2()
    prob = savedNLODEContProblem(Val(2), Val(1000), Val(0), Val(0), 7, ICadvection1000V2, advection1000V2, jacadvection1000V2, SDadvection1000V2, mapadvection1000V2, jacDimensionadvection1000V2)
end
