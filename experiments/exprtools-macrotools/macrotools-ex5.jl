

using MacroTools: prewalk, postwalk, @capture
  macro removeConstantTerm(ex)
    Base.remove_linenums!(ex)
   # dump(ex.args[1].args[1]; maxdepth=8)
    removeConstantTermFunc(ex)
    removeConstantTermFunc2(ex)
  end
function removeConstantTermFunc(ex)
 # counter=0
  prewalk(ex) do x
    if x isa Expr && x.head == :call && x.args[1] in (:*,:/)  # head (call) arg1 (*,/) arg2 (expr,number)  arg3 (expr,number)  arg4 (expr,number) ...
      #mul or divide by q should be ignored
      constantTerm=true
      for i=2:length(x.args)
        
        if (x.args[i] isa Expr && x.args[i].head == :ref && x.args[i].args[1] == :q)  #q*... 
          constantTerm=false
          #@show x.args[i]
           for j=2:length(x.args)
              if j!=i
               if x.args[j] isa Expr   # prevent any future modification to this expression
                  x.args[j].head=:closed
                end
              end
            end#end j loop
 
        elseif x.args[i] isa Expr && x.args[i].head == :call  # another expr(...+....)
          constantTerm=false
        end#end if
      end#end for i loop
      if constantTerm
        x.args[2]=0   #  nullify one of them is sufficient 
      end
      return x
      
    elseif x isa Expr && x.head == :call && x.args[1] in (:+,:-)  #head (call) arg1 (*,/) arg2 (expr,number)  arg3 (expr,number)   possible arg4 and 5 ...
        for i=2:length(x.args)
            if (x.args[i] isa Expr && x.args[i].head == :ref && x.args[i].args[1] != :q)
              x.args[i]=0 
            elseif x.args[i] isa Expr
              removeConstantTermFunc(x.args[i])
            else
              x.args[i]=0 
            end
        end
        return x
    else
       return x
    end
  end
end


function removeConstantTermFunc2(ex)
  
  postwalk(ex) do x
    if x isa Expr && x.head == :closed 
      x.head=:call

    end
      return x
    
  end
end

using BenchmarkTools
function test()
  res=@removeConstantTerm quote 
   # q[1]+4*2                          # q[1] + 0 * 0
   # (7+2)*q[1]                        #(7 + 2) * q[1]
    #7+2*q[1] +5.2+q[2]*3-5            # (0 + 2 * q[1] + 0 + q[2] * 3) - 0
   # q[1]*q[2]+5+3*q[1]+2.2*1+3.2      #q[1] * q(2) + 0 + 3 * q[1] + 0 * 0 + 0
   #q[1]*q[2]*q[1]+3+4*2*3*4           # q[1] * q(2) * q[1] + 0 + 0 * 0 * 3 * 4
  # 4*2*q[1]                           #4 * 2 * q[1] 
  #test parentheses 
  #5+6*(q[1]+3)                         #0 + 6 * (q[1] + 0)
 # 5+6*(q[1]+3+5*2*q[1]+4.2*d[1])       #0 + 6 * (q[1] + 0 + 5 * 2 * q[1] + 0 * d[1])
 #test d[]
 #2+d[2]                                  #0 + 0
 # q[1]+2*d[1]+d[2]                       #q[1] + 0 * d[1] + 0
 5+6*(q[1]+3+5*2*q[1]+4.2*d[1]+q[1]*d[1])  #0 + 6 * (q[1] + 0 + 5 * 2 * q[1] + 0 * d[1] + q[1] * d(1))
 #still not tested for other types of expression and heads such as function call 
  end
  @show res
end

#@btime
 test()
 #=  x=[0.0]
  q=[1,2]
  @show (eval(res)) =#