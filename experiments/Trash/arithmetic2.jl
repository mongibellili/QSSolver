function twoInOne2(ex)
    cachexpr_lengthtracker = Expr(:mongi)
   
     i = 1 #index of cache
    #=  if ex.args[1] isa Number # case where lhs of eq is a number needed to change to taylor to be used for qss
      ex.args[1]=Expr(:call, :addT,ex.args[1],0.0)
      cachexpr = Expr(:ref, :cache)   #multiply needs two caches
       push!(cachexpr.args,1)
       push!( ex.args[1].args, cachexpr)
     
      return ex
     else =#
       prewalk(ex) do x
      # @show x
   #############################minus sign#############################################
      if x isa Expr && x.head == :call && x.args[1] == :- && length(x.args) == 3
        # @show x.args[1]
        if x.args[2] isa Expr && x.args[2].head == :call && x.args[2].args[1] == :- && length(x.args[2].args) == 3
          if i <= ex.args[2]  #do not want to change if mul is deprived of last cache (mul requires 2) and mul has priority add and substr
            push!(x.args, x.args[3])#  args[4]=c
            x.args[3] = x.args[2].args[3]# args[3]=b
            x.args[2] = x.args[2].args[2]# args[2]=a
            x.args[1] = :subsub
            # push!(x.args,ex.args[2])  # add cache for real
            push!(cachexpr_lengthtracker.args,:b) #b is anything ...not needed, just to fill vector
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
          #mulsub 
          if i < ex.args[2]  #less than cuz two caches needed
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
            push!(x.args, x.args[3])#  args[4]=c
            x.args[3] = x.args[2].args[3]# args[3]=b
            x.args[2] = x.args[2].args[2]# args[2]=a
            x.args[1] = :subadd#:µ  # £  § ....
            push!(cachexpr_lengthtracker.args,:b)
            cachexpr = Expr(:ref, :cache)
            push!(cachexpr.args,length(cachexpr_lengthtracker.args))
            #cachexpr.args[2] = index[i]
            push!(x.args, cachexpr)
            i = i + 1
          end
        elseif x.args[2] isa Expr && x.args[2].head == :call && x.args[2].args[1] == :* && length(x.args[2].args) == 3
        #muladd 
           if i < ex.args[2]
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
        if i < ex.args[2]
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
  
    end#end prewalk
    #= if ex.args[1] isa Number # case where lhs of eq is a number needed to change to taylor to be used for qss
      ex.args[1]=Expr(:call, :addT,ex.args[1],0.0)
      cachexpr = Expr(:ref, :cache)   #multiply needs two caches
       push!(cachexpr.args,1)
       push!( ex.args[1].args, cachexpr)
     
      return ex
    end =#
  #end #end if number else prewalk
  end