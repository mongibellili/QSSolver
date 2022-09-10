using InteractiveUtils
import Base.:-
#import Base.:+

Base.:-(x::Float64, y::Float64,z::Float64) = sub_float(sub_float(x, 0.0),z)
# the above won't be called because the ast is represented differently and -(x,y,z) would never exist
# when called (after ast expression change) : error sub_float not defined
a=2.2
b=3.3
c=4.4
d=5.5
expr=:(a-b-c)

dump(expr)
#println(@code_lowered a+b+c)#%1 = a + b %2 = %1 + c  %3 = Core.tuple(Base.:+, %2) %4 = Core._apply_iterate(Base.iterate, Base.afoldl, %3, xs)
#println(@code_lowered a-b-c-d)#Base.sub_float(x, y)
#println(@code_typed a+b+c+d)
#=                 %1 = Base.add_float(a, b)::Float64
            │   %2 = Base.add_float(%1, c)::Float64
            │   %3 = Core.getfield(xs, 1)::Float64
            │   %4 = Base.add_float(%2, %3)::Float64 =#
# println(@code_typed a-b-c-d) #%1 = Base.sub_float(x, y)::Float64
#@show a-b-c
#println(@code_llvm a-b-c) #
#=             define double @julia_-_143(double %0, double %1) #0 {
            top:
            %2 = fsub double %0, %1 =#
#println(@code_llvm a+b+c+d)
#=         define double @"julia_+_162"(double %0, double %1, double %2, double %3) #0 {
        top:
        ;  @ operators.jl:655 within `+` @ float.jl:399
        %4 = fadd double %0, %1
        %5 = fadd double %4, %2
        ;  @ operators.jl:655 within `+`
        ; ┌ @ operators.jl:612 within `afoldl`
        ; │┌ @ float.jl:399 within `+`
            %6 = fadd double %5, %3 =#
        #    println(@code_native a-b-c-d) #
#=             vsubsd  %xmm1, %xmm0, %xmm0
            retq
            nopw    %cs:(%rax,%rax) =#























x = 2.1
y=1.1
z = 3.1
xx = 1.1
function substract(x::Float64,y::Float64,z::Float64,xx::Float64)
    sum=x-y-z-xx
    return sum
end
#@btime add(x,y,z,xx) #100.314 ns (2 allocations: 160 bytes)
#47.179 ns (1 allocation: 80 bytes)  if add-add activated
#@show substract(x,y,z,xx)