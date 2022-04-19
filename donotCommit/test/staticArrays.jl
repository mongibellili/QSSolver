using StaticArrays

v1 =  SVector{Array{Float64}}
display(v1);println()
v1[1]=Array{Float64}[]
#=SA[1, 2, 3]     isa SVector{3,Int}
SA_F64[1, 2, 3] isa SVector{3,Float64}
SA_F32[1, 2, 3] isa SVector{3,Float32}
SA[1 2; 3 4]     isa SMatrix{2,2,Int}
SA_F64[1 2; 3 4] isa SMatrix{2,2,Float64}

# Create an SVector using various forms, using constructors, functions or macros
v1 = SVector(1, 2, 3)
display(v1);println()
v1.data === (1, 2, 3) # SVector uses a tuple for internal storage
v2 = SVector{3,Float64}(1, 2, 3) # length 3, eltype Float64
v3 = @SVector [1, 2, 3]
v4 = @SVector [i^2 for i = 1:10] # arbitrary comprehensions (range is evaluated at global scope)
v5 = zeros(SVector{3}) # defaults to Float64
v6 = @SVector zeros(3)
v7 = SVector{3}([1, 2, 3]) # Array conversions must specify size

# Can get size() from instance or type
size(v1) == (3,)
size(typeof(v1)) == (3,)
=#