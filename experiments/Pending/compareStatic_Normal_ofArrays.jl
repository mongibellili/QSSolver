using StaticArrays
using BenchmarkTools

############contain normal array
struct normalvectorofArray
    v1 = Vector{Array{Float64}}(undef, 2)
    v1
end
struct staticMvectorofArray
    mv1 = MVector{2,Array{Float64}}([], [])
    # println(typeof(mv1))
    mv1
end
struct staticvectorofArray
    sv1 = SVector{2,Array{Float64}}([], [])
    # println(typeof(mv1))
    sv1
end
display(normalvectorofArray());
println();
display(staticMvectorofArray());
println();
display(staticvectorofArray());
println();


function createMutatevectorofArray()
    nv = normalvectorofArray()
    nv[1] = Array{Float64}[]
    # println(typeof(nv))
    push!(nv[1], 1.2)
    nv
end
function createMutateMvectorofArray()
    sv = staticMvectorofArray()
    #sv[1]=Array{Float64}[]
    # println(typeof(sv))
    push!(sv[1], 1.2)
    sv
end
function createMutateStaticvectorofArray()
    sv = staticvectorofArray()
    # sv[1]=Array{Float64}[]
    # println(typeof(sv))
    push!(sv[1], 1.2)
    sv
end


#createMutatevectorofArray()
#createMutateMvectorofArray()
#createMutateStaticvectorofArray()
#@btime createMutatevectorofArray()
#@btime createMutateMvectorofArray()
#@btime createMutateStaticvectorofArray()
