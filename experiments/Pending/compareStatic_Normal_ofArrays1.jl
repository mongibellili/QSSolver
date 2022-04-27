using StaticArrays
using BenchmarkTools
##########contain float64
function normalvector()
 v1 =  Vector{Float64}(undef,2)
 v1
end
function staticMvector()
    mv1 =  MVector{2,Float64}(1.0,2.0)
   # println(typeof(mv1))
    mv1
end
function staticvector()
    sv1 =  SVector{2,Float64}(1.0,2.0)
   # println(typeof(mv1))
    sv1
end
#=
@btime normalvector()
@btime staticMvector()
@btime staticvector()
function createMutatevector()
    nv=normalvector()
   # println(typeof(nv))
    nv[1]=1.2
    nv
end
function createMutateMvector()
    sv=staticMvector()
   # println(typeof(sv))
    sv[1]=1.2
    sv
end
function createMutateStaticvector()
    sv=staticvector()
   # println(typeof(sv))
    sv[1]=1.2
    sv
end
#println(createMutatevector())
#println(createMutateMvector())
#println(createMutateStaticvector()) #illegal
=#










############contain normal array
function normalvectorofArray()
    v1 =  Vector{Array{Float64}}(undef,2)
    v1
   end
   function staticMvectorofArray()
       mv1 =  MVector{2,Array{Float64}}([],[])
      # println(typeof(mv1))
       mv1
   end
   function staticvectorofArray()
       sv1 =  SVector{2,Array{Float64}}([],[])
      # println(typeof(mv1))
       sv1
   end
   
 #  display(normalvectorofArray());println()
  # display(staticMvectorofArray());println()
  # display(staticvectorofArray());println()


   function createMutatevectorofArray()
    nv=normalvectorofArray()
    nv[1]=Array{Float64}[]
   # println(typeof(nv))
    push!(nv[1],1.2)
    nv
end
function createMutateMvectorofArray()
    sv=staticMvectorofArray()
    #sv[1]=Array{Float64}[]
   # println(typeof(sv))
    push!(sv[1],1.2)
    sv
end
function createMutateStaticvectorofArray()
    sv=staticvectorofArray()
   # sv[1]=Array{Float64}[]
   # println(typeof(sv))
   push!(sv[1],1.2)
    sv
end

@btime createMutatevectorofArray()
@btime createMutateMvectorofArray()
@btime createMutateStaticvectorofArray()



























