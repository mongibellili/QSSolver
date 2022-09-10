using TaylorSeries
using InteractiveUtils
using BenchmarkTools
using MacroTools: prewalk, postwalk, @capture
import Base.:-
import Base.:+
use_show_default(true)
#add functionality to arithmetic.jl inside TSfor (f, fc) in ((:+, :(add!)), (:-, :(subst!)))
## Addition and substraction ##
#this section until the macro should be for non-taylor types (floats)
#all taylor ops produced by the macro should be caught in the taylor section below
#(-)(a, b, c)=(-)((-)(a,b),c)
#(+)(a, b, c)=(+)((+)(a,b),c)
#= (addsub)(a, b, c)=(-)((+)(a,b),c)
(subsub)(a, b, c)=(-)((-)(a,b),c)
#(addsub)(cache,a, b, c)=(-)((+)(cache,a,b),c)# cache is at begining unlike the expression 
                                          #because taylor func flipped them: oldRes becomes at begining of exprss which makes sense.
(subadd)(a, b, c)=(+)((-)(a,b),c)
(subT)(a, b,c)=(-)(a,b)
(addT)(a, b,c)=(+)(a,b) =#


macro changeAST(ex)
  Base.remove_linenums!(ex)
  #dump(ex; maxdepth=16)
  twoInOne(ex)
 
  #esc(ex.args[1])# return result only , remove this for debugging
  @show ex.args[1]

  return nothing
end

function twoInOne(ex)
  cachexpr_lengthtracker = Expr(:mongi)
 
   i = 1 #index of cache
   
  prewalk(ex) do x
    # @show x
    #############################minus sign#############################################
    if x isa Expr && x.head == :call && x.args[1] == :- && length(x.args) == 3
      # @show x.args[1]
      if x.args[2] isa Expr && x.args[2].head == :call && x.args[2].args[1] == :- && length(x.args[2].args) == 3
        if i <= ex.args[2]
          push!(x.args, x.args[3])#  args[4]=c
          x.args[3] = x.args[2].args[3]# args[3]=b
          x.args[2] = x.args[2].args[2]# args[2]=a
          x.args[1] = :subsub
          # push!(x.args,ex.args[2])  # add cache for real
          push!(cachexpr_lengthtracker.args,:b)
          cachexpr = Expr(:ref, :cache)
          push!(cachexpr.args,length(cachexpr_lengthtracker.args))
          #cachexpr.args[2] = index[i]
          push!(x.args, cachexpr)
          i = i + 1
        end
        # @show x.args
      elseif x.args[2] isa Expr && x.args[2].head == :call && x.args[2].args[1] == :+ && length(x.args[2].args) == 3
        if i <= ex.args[2]
          push!(x.args, x.args[3])#  args[4]=c
          x.args[3] = x.args[2].args[3]# args[3]=b
          x.args[2] = x.args[2].args[2]# args[2]=a
          x.args[1] = :addsub # £ µ § ~....        
          push!(cachexpr_lengthtracker.args,:b)
          cachexpr = Expr(:ref, :cache)
          push!(cachexpr.args,length(cachexpr_lengthtracker.args))
          #cachexpr.args[2] = index[i]
          push!(x.args, cachexpr)
          i = i + 1
        end
      elseif x.args[2] isa Expr && x.args[2].head == :call && x.args[2].args[1] == :* && length(x.args[2].args) == 3
        #mulsub not implemented yet
        if i <= ex.args[2]
          push!(x.args, x.args[3])#  args[4]=c
          x.args[3] = x.args[2].args[3]# args[3]=b
          x.args[2] = x.args[2].args[2]# args[2]=a
          x.args[1] = :mulsub # £ µ § ~....
          push!(cachexpr_lengthtracker.args,:b)
          cachexpr1 = Expr(:ref, :cache)
          push!(cachexpr1.args,length(cachexpr_lengthtracker.args))
          #cachexpr.args[2] = index[i]
          push!(x.args, cachexpr1)
          i = i + 1
          push!(cachexpr_lengthtracker.args,:b)
          cachexpr2 = Expr(:ref, :cache)   #multiply needs two caches
          push!(cachexpr2.args,length(cachexpr_lengthtracker.args))
          #cachexpr.args[2] = index[i]
          push!(x.args, cachexpr2)
          i=i+1
        end

        # elseif  (x.args[2] isa Expr && x.args[2].head != :call) || !(x.args[2] isa Expr)
      else
        if i <= ex.args[2]
          x.args[1] = :subT  # symbol changed cuz avoid type taylor piracy    
          push!(cachexpr_lengthtracker.args,:b)
          cachexpr = Expr(:ref, :cache)
          push!(cachexpr.args,length(cachexpr_lengthtracker.args))
          #cachexpr.args[2] = index[i]
          push!(x.args, cachexpr)
          i = i + 1
        end


      end
    elseif x isa Expr && x.head == :call && x.args[1] == :- && length(x.args) == 2
      if i <= ex.args[2]
        x.args[1] = :negateT  # symbol changed cuz avoid type taylor piracy
        push!(cachexpr_lengthtracker.args,:b)
        cachexpr = Expr(:ref, :cache)
        push!(cachexpr.args,length(cachexpr_lengthtracker.args))
        #cachexpr.args[2] = index[i]
        push!(x.args, cachexpr)
        i = i + 1
      end


      ############################### plus sign#######################################
    elseif x isa Expr && x.head == :call && x.args[1] == :+ && length(x.args) == 3
      if x.args[2] isa Expr && x.args[2].head == :call && x.args[2].args[1] == :- && length(x.args[2].args) == 3
        if i <= ex.args[2]
          push!(cachexpr_lengthtracker.args,:b)
          cachexpr = Expr(:ref, :cache)
          push!(cachexpr.args,length(cachexpr_lengthtracker.args))
          #cachexpr.args[2] = index[i]
          push!(x.args, cachexpr)
          i = i + 1
          push!(x.args, x.args[3])#  args[4]=c
          x.args[3] = x.args[2].args[3]# args[3]=b
          x.args[2] = x.args[2].args[2]# args[2]=a
          x.args[1] = :subadd#:µ  # £  § ....
        end
      elseif x.args[2] isa Expr && x.args[2].head == :call && x.args[2].args[1] == :* && length(x.args[2].args) == 3
        if i <= ex.args[2]
          #muladd not implemented yet
          push!(x.args, x.args[3])#  args[4]=c
          x.args[3] = x.args[2].args[3]# args[3]=b
          x.args[2] = x.args[2].args[2]# args[2]=a
          x.args[1] = :muladd#:µ  # £  § ....
          push!(cachexpr_lengthtracker.args,:b)
          cachexpr1 = Expr(:ref, :cache)
          push!(cachexpr1.args,length(cachexpr_lengthtracker.args))
          #cachexpr.args[2] = index[i]
          push!(x.args, cachexpr1)
          i = i + 1
          push!(cachexpr_lengthtracker.args,:b)
          cachexpr2 = Expr(:ref, :cache)   #multiply needs two caches
          push!(cachexpr2.args,length(cachexpr_lengthtracker.args))
          #cachexpr.args[2] = index[i]
          push!(x.args, cachexpr2)
          i=i+1
        end


      else
        if i <= ex.args[2]
          x.args[1] = :addT
          push!(cachexpr_lengthtracker.args,:b)
          cachexpr = Expr(:ref, :cache)
          push!(cachexpr.args,length(cachexpr_lengthtracker.args))
          #cachexpr.args[2] = index[i]
          push!(x.args, cachexpr)
          i = i + 1
        end


      end
    elseif x isa Expr && x.head == :call && x.args[1] == :+ && (4 <= length(x.args) <= 9)
      if i <= ex.args[2]
        x.args[1] = :addT
        push!(cachexpr_lengthtracker.args,:b)
        cachexpr = Expr(:ref, :cache)
        push!(cachexpr.args,length(cachexpr_lengthtracker.args))
        #cachexpr.args[2] = index[i]
        push!(x.args, cachexpr)
        i = i + 1
      end


      ############################### multiply sign#######################################
      #never happens :#  elseif x isa Expr && x.head == :call && x.args[1]==:* && length(x.args)==3 to get addmul or submul

    elseif x isa Expr && x.head == :call && (x.args[1] == :*) && (3 <= length(x.args) <= 7)
      if i <= ex.args[2]
        x.args[1] = :mulT
        push!(cachexpr_lengthtracker.args,:b)
          cachexpr1 = Expr(:ref, :cache)
          push!(cachexpr1.args,length(cachexpr_lengthtracker.args))
          #cachexpr.args[2] = index[i]
          push!(x.args, cachexpr1)
          i = i + 1
          push!(cachexpr_lengthtracker.args,:b)
          cachexpr2 = Expr(:ref, :cache)   #multiply needs two caches
          push!(cachexpr2.args,length(cachexpr_lengthtracker.args))
          #cachexpr.args[2] = index[i]
          push!(x.args, cachexpr2)
          i=i+1
      end


    end
    return x

  end
end



res = @changeAST 2.0*3.0+a[1]*b[1], 10#-e-f+i) #10 or any number has to be hardcoded...you cannot pass any symbol cuz this is parse time
                        #you can have a whole different integrator with different cache size (specialized on via Val())
@show res
