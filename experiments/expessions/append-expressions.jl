code = Expr(:block)
v=quote
    function  $(Symbol(:f_, 1))(q::Vector{Taylor0{Float64}},d::Vector{Float64}, t::Taylor0{Float64},cache::Vector{Taylor0{Float64}})
 end
 push!(code.args,v)
 v=quote
    if j==1
 end
# push!(code.args, Meta.parse("if j==1"))
push!(code.args,v)
push!(code.args, equs[1])
v=quote
end
end
push!(code.args,v)
@show code