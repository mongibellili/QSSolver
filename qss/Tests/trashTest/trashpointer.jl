

#= resppi = pointer(Vector{Float64}(undef, 2))
unsafe_store!(resppi, -1.0, 1);unsafe_store!(resppi, -1.0, 2)
resTupi=(unsafe_load(resppi,1),unsafe_load(resppi,2))
@show resTupi =#

f=(1.0,12.3,6.5)
g=max(f...,9.3)
@show g