using StaticArrays
#v1 = @SVector randn(Float64, 40)
#display(v1)
#A = rand(3);
#display(A)
#= using Base.Cartesian
display(@macroexpand @nloops 2 i A begin
    r += @nref 2 A i
end) =#
#= a=Float64(0)

if a==0
    println("a is zero")
end =#

#= struct DataSvector
    x::SVector{N,Float64} where {N} 
    y::Int
  end

n=5
v1 = @SVector zeros(n) 
p = DataSvector(v1, 5)  
display(p) =#
#= arr=[[],[],[]]
#= arr=[]
push!(arr,[]) =#
t=tuple(arr...)
display(t) =#
#display(sign(5))
#= s = 0
for i = 1:10
     s += i
end
println(s) =#
#= function test()
    s = 0
    for i = 1:10
        global s += i
    end
    s
end
s=2
println(test())
println(s) =#

                negEv_disArr=@SVector zeros(D) #  discrete
              negEv_conArr=@SVector zeros(T)  #  continous
              for j=1:length(negEvExp.args)  # j coressponds the number of statements under one negEvent
                 negEvLHS=negEvExp.args[j].args[1]
                 negEvRHS=negEvExp.args[j].args[2]
                 #dump(negEvLHS)# symbol at the LHS
                 basicLHS = convert(Basic, negEvLHS)
                 
                 discVarpositionArray = indexin(basicLHS, d)
                 
                 if !(discVarpositionArray[1] === nothing)
                  negEv_disArr=setindex(negEv_disArr,negEvRHS,discVarnegitionArray[1])
                 else # lhs is not a disc var 
                      conVarpositionArray = indexin(basicLHS, u)
                      if !( conVarpositionArray[1]=== nothing)
                          negEv_conArr=setindex(negEv_conArr,negEvRHS,conVarnegitionArray[1])
                      else
                              println("LHS is neither a cont nor a discr var!!")
                      end
                  end
              end
