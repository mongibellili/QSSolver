using StaticArrays
macro make_struct(struct_name, schema...)
    fields=[:($(entry.args[1])::$(entry.args[2])) for entry in schema]
    esc(quote struct $struct_name
        $(fields...)
        end
    end)
end
println(@macroexpand @make_struct sim (x, Int) (y, Float64) (wRegister, MVector{2,Float64}))


#= v1=@MVector[1.1,2.6]
p=sim(1,2.0,v1)
display(p) =#
function problem()
    @make_struct sim (x, Int) (y, Float64) (wRegister, MVector{2,Float64})
end
function solve()
    v1=@MVector[1.1,2.6]
    p=sim(1,2.0,v1)
    display(p)
end


problem()

#solve()