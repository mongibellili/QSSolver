#using ExprTools
using MacroTools: prewalk, postwalk, @capture
  macro removeConstantTerm(ex)
    Base.remove_linenums!(ex)
    dump(ex.args[1]; maxdepth=8)
    postwalk(ex) do x
      if x isa Expr && x.head == :call && x.args[1] in (:+,:-) 
          if (@capture(x.args[3], n_)&& n isa Number)
            x.args[3]=0  
            return x
          elseif @capture(x.args[2], n_)&& n isa Number 
            x.args[2]=0  
            return x
          else
            return x
          end
      else
         return x
      end
    end
  end

  res=@removeConstantTerm quote 
    7+2*q[1] +5.2+q[2]*3-5
  end
  @show res

  x=[0.0]
  q=[1,2]
  @show (eval(res))