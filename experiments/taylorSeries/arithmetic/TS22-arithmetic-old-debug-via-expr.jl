import Base.:-
import Base.:+
import Base.:*
using MacroTools: prewalk, postwalk, @capture
#a=2.2;b=3.3;c=4.4;d=5.5;e=2.2;f=3.3;i=4.4;j=5.5


(-)(a, b, c)=(-)((-)(a,b),c)
(+)(a, b, c)=(+)((+)(a,b),c)
(addsub)(a, b, c)=(-)((+)(a,b),c)
(subadd)(a, b, c)=(+)((-)(a,b),c)

macro changeSubstractionAST(ex)
  Base.remove_linenums!(ex)
   dump(ex; maxdepth=8)
  esc(changeSubstractionASTFunc(ex))
# changeSubstractionASTFunc(ex)
end

function changeSubstractionASTFunc(ex)
    prewalk(ex) do x
     # @show x
     if x isa Expr && x.head == :call && x.args[1]==:- && length(x.args)==3 
       # @show x.args[1]
        if x.args[2] isa Expr && x.args[2].head == :call && x.args[2].args[1]==:- && length(x.args[2].args)==3 
          push!(x.args,x.args[3])#  args[4]=c
          x.args[3]=x.args[2].args[3]# args[3]=b
          x.args[2]=x.args[2].args[2]# args[2]=a
         
        elseif x.args[2] isa Expr && x.args[2].head == :call && x.args[2].args[1]==:+ && length(x.args[2].args)==3 
          push!(x.args,x.args[3])#  args[4]=c
          x.args[3]=x.args[2].args[3]# args[3]=b
          x.args[2]=x.args[2].args[2]# args[2]=a
          x.args[1]=:addsub # £ µ § ~....
         end
     elseif x isa Expr && x.head == :call && x.args[1]==:+ && length(x.args)==3 
         if x.args[2] isa Expr && x.args[2].head == :call && x.args[2].args[1]==:- && length(x.args[2].args)==3 
            push!(x.args,x.args[3])#  args[4]=c
            x.args[3]=x.args[2].args[3]# args[3]=b
            x.args[2]=x.args[2].args[2]# args[2]=a
            x.args[1]=:subadd#:µ  # £  § ....
           end
     end
     return x
      
    end
  end
  expr=:(a-b+c-d)
 # dump(expr)
 res= @changeSubstractionAST :(3*b) 
 @show res
#= expr=:(a+b+c+d+e+f+h*i+j)
dump(expr) =#
#=  head: Symbol call
args: Array{Any}((9,))
  1: Symbol +
  2: Symbol a
  3: Symbol b
  4: Symbol c
  5: Symbol d
  6: Symbol e
  7: Symbol f
  8: Expr
    head: Symbol call
    args: Array{Any}((3,))
      1: Symbol *
      2: Symbol h
      3: Symbol i
  9: Symbol j =#


 #=  expr=:(a+b+c+d+e+f+h*i-j)
  #dump(expr)
prewalk( expr) do x  
     @show x 
end =#



#= expr=:(a-b*c+d-e-f+h*i+j)
dump(expr) =#
#= head: Symbol call
  args: Array{Any}((4,))
    1: Symbol +
    2: Expr
      head: Symbol call
      args: Array{Any}((3,))
        1: Symbol -
        2: Expr
          head: Symbol call
          args: Array{Any}((3,))
            1: Symbol -
            2: Expr
              head: Symbol call
              args: Array{Any}((3,))
                1: Symbol +
                2: Expr
                3: Symbol d
            3: Symbol e
        3: Symbol f
    3: Expr
      head: Symbol call
      args: Array{Any}((3,))
        1: Symbol *
        2: Symbol h
        3: Symbol i
    4: Symbol j =#





















































#= 
function substract(x::Float64,y::Float64,z::Float64,xx::Float64)
    sum=x-y-z-xx
    return sum
end =#
#@btime add(x,y,z,xx) #100.314 ns (2 allocations: 160 bytes)
#47.179 ns (1 allocation: 80 bytes)  if add-add activated
#@show substract(x,y,z,xx)