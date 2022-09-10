
using TaylorSeries
using StaticArrays
using InteractiveUtils
using BenchmarkTools
using MacroTools: prewalk, postwalk, @capture
import Base.:-
import Base.:+
use_show_default(true)

(addsub)(a, b, c)=(-)((+)(a,b),c)
(subsub)(a, b, c)=(-)((-)(a,b),c)
(subadd)(a, b, c)=(+)((-)(a,b),c)# this should not exist cuz addsub has 3 args + cache

  function addTfoldl(op, res, bs...)
      l = length(bs)
     
      #i =  0;  l == i && return a  # already checked
      l == 1 && return res 
      i =  2; l == i && return op(res, bs[1],bs[end])
      i =  3;  l == i && return op(res, bs[1], bs[2],bs[end])
      i =  4;  l == i && return op(op(res, bs[1], bs[2],bs[end]), bs[3],bs[end])
      i =  5;  l == i && return op(op(res, bs[1], bs[2],bs[end]), bs[3],bs[4],bs[end])
      i =  6;  l == i && return op(op(op(res, bs[1], bs[2],bs[end]), bs[3],bs[4],bs[end]),bs[5],bs[end])
      i =  7;  l == i && return op(op(op(res, bs[1], bs[2],bs[end]), bs[3],bs[4],bs[end]),bs[5],bs[6],bs[end])
      i =  8;  l == i && return op(op(op(op(res, bs[1], bs[2],bs[end]), bs[3],bs[4],bs[end]),bs[5],bs[6],bs[end]),bs[7],bs[end])
      i =  9;y=op(op(op(op(res, bs[1], bs[2],bs[end]), bs[3],bs[4],bs[end]),bs[5],bs[6],bs[end]),bs[7],bs[8],bs[end]);  l == i && return y
      #last i creates y and passes it to the next loop: note that from this point, op will allocate
      # either keep increasing i manually or warn the users . (inserting -0 will avoid allocation lol!)
       for i in i:l-1 # -1 cuz last one is cache
          y = op(y, bs[i],bs[end])
      end
      return y 
  end
    function addT(a, b, c,d, xs...)
      if length(xs)!=0# if ==0 case below where d==cache
     
        addTfoldl(addT, addT( addT(a,b,c,xs[end]) , d ,xs[end] ), xs...)
      end

    end

    function mulTfoldl(op, res, bs...)
      l = length(bs)
    #  @show l
      #println("entered multfoldl")
      #i =  0;  l == i && return a  # already checked
      i =  2;  l == i && return res 
      i =  3;  l == i && return              op(res, bs[1],bs[end-1],bs[end])
      i =  4; l == i && return           op(op(res, bs[1],bs[end-1],bs[end]),bs[2],bs[end-1],bs[end])
      i =  5;  l == i && return        op(op(op(res, bs[1],bs[end-1],bs[end]),bs[2],bs[end-1],bs[end]),bs[3],bs[end-1],bs[end])

      #the macro never creates long mulT in the first place. cuz it starts allocates (tested)
      #= if l==5
       
      # @show bs
        return     op(op(op(op(res, bs[1],bs[end]),bs[2],bs[end]),bs[3],bs[end]),bs[4],bs[end])
      end =#

      #since bs[end][1](cache1 has previous result then i can use it directly op(cache1,b,cache)) shortlines :)
     #=  i =  5;  l == i &&  return     op(op(op(op(res, bs[1],bs[end]),bs[2],bs[end]),bs[3],bs[end]),bs[4],bs[end])
      i =  6; @show i; l == i &&  return op(op(op(op(op(res, bs[1],bs[end]), bs[2],bs[end]), bs[3],bs[end]),bs[4],bs[end]),bs[5],bs[end])
      i =  7; @show i; l == i &&  return op(op(op(op(op(op(res, bs[1],bs[end]), bs[2],bs[end]), bs[3],bs[end]),bs[4],bs[end]),bs[5],bs[end]),bs[6],bs[end])
      i =  8; @show i; l == i && return op(op(op(op(op(op(op(res, bs[1],bs[end]), bs[2],bs[end]), bs[3],bs[end]),bs[4],bs[end]),bs[5],bs[end]),bs[6],bs[end]),bs[7],bs[end])
      i =  9;@show i;y=op(op(op(op(op(op(op(op(res, bs[1],bs[end]), bs[2],bs[end]), bs[3],bs[end]),bs[4],bs[end]),bs[5],bs[end]),bs[6],bs[end]),bs[7],bs[end]),bs[8],bs[end]);  l == i && return y
      #last i creates y and passes it to the next loop: note that from this point, op will allocate
      # either keep increasing i manually or warn the users . (inserting -0 will avoid allocation lol!)
       for i in i:l-1 # -1 cuz last one is cache
        @show i;
          y = op(y, bs[i],bs[end])
      end
     # println("at end")
      return y  =#
    end
    function mulT(a, b, c, xs...)
     # println("entered mult(a b c)")
      if length(xs)>1# length(xs)==0 will never happen, length(xs)==1 is the case far-below
     #   println("fuck multi xs= ",xs[end])
        mulTfoldl(mulT, mulT( mulT(a,b,xs[end-1],xs[end]) , c,xs[end-1] ,xs[end] ), xs...)
      end

    end
 
  macro changeAST(ex)
    Base.remove_linenums!(ex)
    # dump(ex; maxdepth=8)
    twoInOne(ex)
   
  # dump( ex.args[1])
   #= @show (ex.args[1])# return  
  return nothing =#
  esc(ex.args[1])# return 
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
          #mulsub 
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
        #muladd 
           if i <= ex.args[2]
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
  
  ##########################Taylor section##############################""
 ##########################Taylor section##############################""

#all these functions can be written in terms of each other without cost (I think). i will test later.

  # all methods should have new names . no type piracy!!!
  function addsub(a::Taylor1{T}, b::Taylor1{T},c::Taylor1{T},cache::Taylor1{T}) where {T<:Number}
    #v = similar(a.coeffs)
    
    #@show cache
   # @__dot__ cache.coeffs = addsub(cache.coeffs,a.coeffs, b.coeffs,c.coeffs)
   @__dot__ cache.coeffs = addsub(a.coeffs, b.coeffs,c.coeffs)
    return cache
  end
  function addsub(a::Taylor1{T}, b::Taylor1{T},c::T,cache::Taylor1{T}) where {T<:Number}
    #v = similar(a.coeffs)
    
    #@show cache
   # @__dot__ cache.coeffs = addsub(cache.coeffs,a.coeffs, b.coeffs,c.coeffs)
            cache.coeffs.=b.coeffs  
            cache[0]=cache[0]-c
   @__dot__ cache.coeffs = (+)(cache.coeffs, a.coeffs)
    #return 0#$T(v, a.order)
    return cache
  end
  function addsub(a::T, b::Taylor1{T},c::Taylor1{T},cache::Taylor1{T}) where {T<:Number}
            cache.coeffs.=b.coeffs  
            cache[0]=a+ cache[0]
           # @show a
   @__dot__ cache.coeffs = (-)(cache.coeffs, c.coeffs)
    return cache
  end
  addsub(a::Taylor1{T}, b::T,c::Taylor1{T},cache::Taylor1{T}) where {T<:Number}=addsub(b, a ,c,cache)#addsub(b::T, a::Taylor1{T} ,c::Taylor1{T},cache::Taylor1{T}) where {T<:Number}
 
  function addsub(a::T, b::T,c::Taylor1{T},cache::Taylor1{T}) where {T<:Number}
    #cache.coeffs.=-c.coeffs #not tested
    #cache.coeffs.=0.0
    @__dot__ cache.coeffs = (-)(c.coeffs)
   # @__dot__ cache.coeffs=(-)(cache.coeffs, c.coeffs)#   ==-c
    cache[0]=a+b+ cache[0] 
    return cache
end



function addsub(c::Taylor1{T},a::T, b::T,cache::Taylor1{T}) where {T<:Number}
  cache.coeffs.=c.coeffs 
  cache[0]=a+ cache[0]-b
  return cache
end
 addsub(a::T, c::Taylor1{T},b::T,cache::Taylor1{T}) where {T<:Number}= addsub(c,a, b,cache) 
 
 function addsub(c::T,a::T, b::T,cache::Taylor1{T}) where {T<:Number}
  ##[[[[[[[[[!!!!!!!!!!!!!!!!!]]]]]]]]]cache[1]...should be emptied
  #cache.coeffs.=0.0   # since C is number this cant be an op in the middle (all in middle op return taylor)
                      # still the cache can be full from another completely diff equation
                      # emptying the cache does not affect c (in this case c is not cache)
                      #i think emptying the cashe once before the equation is cheaper than emptying it 
                      #for every op
  cache[0]=a+c-b
  return cache
end

###########"negate###########
function negateT(a::Taylor1{T},cache::Taylor1{T}) where {T<:Number}
  @__dot__ cache.coeffs = (-)(a.coeffs)
  return cache
end
#= function negateT(b::T,cache::Taylor1{T}) where {T<:Number} # no need since ast of -Number is direct
  cache.coeffs.=0.0
  cache[0]=-b
  return cache
end =#
#################################################subsub########################################################""


  function subsub(a::Taylor1{T}, b::Taylor1{T},c::Taylor1{T},cache::Taylor1{T}) where {T<:Number}
    
    @__dot__ cache.coeffs = subsub(a.coeffs, b.coeffs,c.coeffs)
    
    return cache
  end

  function subsub(a::Taylor1{T}, b::Taylor1{T},c::T,cache::Taylor1{T}) where {T<:Number}
    #@show a
   cache.coeffs.=b.coeffs  
   cache[0]=cache[0]+c     #(a-b-c=a-(b+c))
  @__dot__ cache.coeffs = (-)(a.coeffs, cache.coeffs)
    return cache
  end
  subsub(a::Taylor1{T}, b::T,c::Taylor1{T},cache::Taylor1{T}) where {T<:Number}=subsub(a, c,b,cache) 
   #------!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!----------------
 #


 #      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 #--------!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-------
  function subsub(a::T,b::Taylor1{T}, c::Taylor1{T},cache::Taylor1{T}) where {T<:Number}
   # @show c
# bug: this is the first time cache is used upfront..!!!cache not empty
#cache.coeffs.=0.0
@__dot__ cache.coeffs = (-)(b.coeffs) # like before since a is number modifying the cache is fine
  # @__dot__ cache.coeffs=(-)(cache.coeffs, b.coeffs)
   cache[0]=a+ cache[0]
  # @show a
@__dot__ cache.coeffs = (-)(cache.coeffs, c.coeffs)
return cache



   #=  cache.coeffs.=-b.coeffs #not tested
    cache[0]=a+cache[0]     #(a-b-c=a-(b+c))
   @__dot__ cache.coeffs = (-)(cache.coeffs, c.coeffs)
     return cache =#
   end


   function subsub(a::Taylor1{T}, b::T,c::T,cache::Taylor1{T}) where {T<:Number}
    cache.coeffs.=a.coeffs  # if first op then needed...else cache is already  a
                                #maybe later have firstOp for every op
    cache[0]=cache[0]-b-c     
     return cache
   end

   function subsub( a::T,b::Taylor1{T},c::T,cache::Taylor1{T}) where {T<:Number}
    #cache.coeffs.=-a.coeffs  #not tested
   # cache.coeffs.=0.0 #]]]]]]]]]]]]]]]]]]]]]]$$$$$$$$$$$$$$$$********************
    @__dot__ cache.coeffs = (-)(b.coeffs)
   # @__dot__ cache.coeffs=(-)(cache.coeffs, a.coeffs)
    cache[0]=cache[0]+a-c     
     return cache
   end
   subsub( b::T,c::T,a::Taylor1{T},cache::Taylor1{T}) where {T<:Number}=subsub( b,a,c,cache) 


   function subsub( a::T,b::T,c::T,cache::Taylor1{T}) where {T<:Number}
    #cache.coeffs.=0.0 # a is number, ok to empty cache
     #i think emptying the cashe once before the equation is cheaper than emptying it 
                      #for every op
    cache[0]=a-b-c    
     return cache
   end


   #################################subadd####################################""
  #= function subadd(a::Taylor1{T}, b::Taylor1{T},c::Taylor1{T},cache::Taylor1{T}) where {T<:Number}

   @__dot__ cache.coeffs = subadd(a.coeffs, b.coeffs,c.coeffs)
    #return 0#$T(v, a.order)
    return cache
  end =#

  #subadd(a,b,c)::Taylor1{T} where {T<:Number}=addsub(a,c,b)   #not tested
  function subadd(a::P,b::Q,c::R,cache::Taylor1{T}) where {P,Q,R <:Union{Taylor1,Number},T<:Number}
    addsub(a,c,b,cache)  
end

  ##################################""subT################################
  function subT(a::Taylor1{T}, b::Taylor1{T},cache::Taylor1{T}) where {T<:Number}
   @__dot__ cache.coeffs = (-)(a.coeffs, b.coeffs)
    return cache
  end
  function subT(a::Taylor1{T}, b::T,cache::Taylor1{T}) where {T<:Number}
    #= for i=1:length(cache.coeffs)
        cache.coeffs[i]=a.coeffs[i] 
    end =#
    @__dot__ cache.coeffs = a.coeffs  # not needed for in middle ops
   cache[0]=cache[0]-b    
    #@__dot__ cache.coeffs = (-)(a.coeffs, b.coeffs)
     return cache
   end
   function subT(a::T, b::Taylor1{T},cache::Taylor1{T}) where {T<:Number}
    #cache.coeffs.=-b.coeffs  #not tested
    #cache.coeffs.=0.0
    @__dot__ cache.coeffs = (-)(b.coeffs) # a is number ok to modify cache
   # @__dot__ cache.coeffs=(-)(cache.coeffs, b.coeffs)
   cache[0]=cache[0]+a    
    #@__dot__ cache.coeffs = (-)(a.coeffs, b.coeffs)
     return cache
   end
   function subT( a::T,b::T,cache::Taylor1{T}) where {T<:Number}
    #cache.coeffs.=0.0 # a is number ok to empty
     #i think emptying the cashe once before the equation is cheaper than emptying it 
                        #for every op
    cache[0]=a-b   
    #@show cache
     return cache
   end
function addT(a::Taylor1{T}, b::Taylor1{T},cache::Taylor1{T}) where {T<:Number}
  
        @__dot__ cache.coeffs = (+)(a.coeffs, b.coeffs)
                return cache
end
function addT(a::Taylor1{T}, b::T,cache::Taylor1{T}) where {T<:Number}
  #@show a==cache
    cache.coeffs.=a.coeffs  #not needed for inmiddle ops
    cache[0]=cache[0]+b 
          return cache
end
addT(b::T,a::Taylor1{T},cache::Taylor1{T}) where {T<:Number}=addT(a, b,cache) 

function addT( a::T,b::T,cache::Taylor1{T}) where {T<:Number}
  #cache.coeffs.=0.0 # a is number ok to empty
   #i think emptying the cashe once before the equation is cheaper than emptying it 
                      #for every op
  cache[0]=a+b   
  #@show cache
   return cache
 end

 #add Three vars a b c
function addT(a::Taylor1{T}, b::Taylor1{T},c::Taylor1{T},cache::Taylor1{T}) where {T<:Number}
 # println("add abc a=",a)
 # println("add abc cache=",cache)
  
 # println("equal cache=",cache==a)

@__dot__ cache.coeffs = (+)(a.coeffs, b.coeffs,c.coeffs)
  return cache
end   
function addT(a::Taylor1{T}, b::Taylor1{T},c::T,cache::Taylor1{T}) where {T<:Number}
  
  cache.coeffs.=a.coeffs  #not needed for in middle ops
    cache[0]=cache[0]+c 
  @__dot__ cache.coeffs = (+)(cache.coeffs, b.coeffs)
    return cache
  end   
  addT(a::Taylor1{T},c::T, b::Taylor1{T},cache::Taylor1{T}) where {T<:Number}=addT(a, b,c,cache) 
  addT(c::T, a::Taylor1{T}, b::Taylor1{T},cache::Taylor1{T}) where {T<:Number}=addT(a, b,c,cache)

  function addT(a::Taylor1{T}, b::T,c::T,cache::Taylor1{T}) where {T<:Number}
  
    cache.coeffs.=a.coeffs  #not needed for inmiddle ops
      cache[0]=cache[0]+c+b
      return cache
    end 
    addT( b::T,a::Taylor1{T},c::T,cache::Taylor1{T}) where {T<:Number}=addT(a, b,c,cache) 
    addT( b::T,c::T,a::Taylor1{T},cache::Taylor1{T}) where {T<:Number}=addT(a, b,c,cache) 

    function addT( a::T,b::T,c::T,cache::Taylor1{T}) where {T<:Number}
      #cache.coeffs.=0.0 #a is number ok to empty
       #i think emptying the cashe once before the equation is cheaper than emptying it 
                      #for every op
      cache[0]=a+b+c    
       return cache
     end


#########################mul###########################

function mulT(a::Taylor1{T}, b::T,cache1::Taylor1{T},cache2::Taylor1{T}) where {T<:Number}
  fill!(cache2.coeffs, b)
  @__dot__ cache1.coeffs = a.coeffs * cache2.coeffs  ##fixed broadcast dimension mismatch
#maybe clear cache2 and see if it allocs
  return cache1
end
mulT(a::T,b::Taylor1{T}, cache1::Taylor1{T},cache2::Taylor1{T}) where {T<:Number} = mulT(b , a,cache1,cache2)

function mulT(a::Taylor1{T}, b::Taylor1{T},cache1::Taylor1{T},cache2::Taylor1{T}) where {T<:Number}
 for k in eachindex(a)
    #println("inloop a= ",a)
   # 
      @inbounds cache2[k] = a[0] * b[k]
      @inbounds for i = 1:k
        cache2[k] += a[i] * b[k-i]
      end
     # println("inloop cache= ",cache)
  end

 @__dot__ cache1.coeffs = cache2.coeffs
#maybe clear cache2 and see if it allocs
  return cache1
end
function mulT(a::T, b::T,cache1::Taylor1{T},cache2::Taylor1{T}) where {T<:Number}
  #cache1 needs to be clean: a clear here will alloc...this is the only "mul" that wants a clean cache
  #sometimes the cache can be dirty at only first position: it will not throw an assert error!!!
 # cache[1].coeffs.=0.0 #it causes two allocs for mul(cst,cst,a,b)
  cache1[0]=a*b  
       return cache1
end
function muladd(a::P,b::Q,c::R,cache1::Taylor1{T},cache2::Taylor1{T}) where {P,Q,R <:Union{Taylor1,Number},T<:Number}
      addT(mulT(a, b,cache1,cache2),c,cache1)
 end

 function mulsub(a::P,b::Q,c::R,cache1::Taylor1{T},cache2::Taylor1{T}) where {P,Q,R <:Union{Taylor1,Number},T<:Number}
  subT(mulT(a, b,cache1,cache2),c,cache1)
end



a =   [ Taylor1([1.0,1.1,1.2],2),Taylor1([1.3,1.4,1.5],2)]
b =   [ Taylor1([2.0,2.1,2.2],2),Taylor1([2.3,2.4,2.5],2)]
c =   [ Taylor1([3.0,3.1,3.2],2),Taylor1([3.3,3.4,3.5],2)]
d =    [Taylor1([4.0,4.1,4.2],2),Taylor1([4.3,4.4,4.5],2)]
#cache of size 12: we can later have the integrator specialize based on the prefered cache size via Val
#between equations cache has to be cleared...we can add an index to tell which cache is dirty
cache=[Taylor1([0.0,0.0,0.0],2),Taylor1([0.0,0.0,0.0],2),Taylor1([0.0,0.0,0.0],2),Taylor1([0.0,0.0,0.0],2),Taylor1([0.0,0.0,0.0],2),Taylor1([0.0,0.0,0.0],2),Taylor1([0.0,0.0,0.0],2),Taylor1([0.0,0.0,0.0],2),Taylor1([0.0,0.0,0.0],2),Taylor1([0.0,0.0,0.0],2),Taylor1([0.0,0.0,0.0],2),Taylor1([0.0,0.0,0.0],2)]

function fastArith(a::Vector{Taylor1{Float64}},b::Vector{Taylor1{Float64}},c::Vector{Taylor1{Float64}},d::Vector{Taylor1{Float64}},cache::Vector{Taylor1{Float64}})
  sum=@changeAST 2.0*2.0*a[2]*a[1],6
 #sum=@changeAST 2.0*3.0+a[1]*b[1],6
end
function normalArith(a::Vector{Taylor1{Float64}},b::Vector{Taylor1{Float64}},c::Vector{Taylor1{Float64}},d::Vector{Taylor1{Float64}},cache::Vector{Taylor1{Float64}})
  sum= 2.0*2.0*a[2]*a[1]
end
#=  @show normalArith(a,b,c,d,cache)
@show fastArith(a,b,c,d,cache) =#

#@btime normalArith(a,b,c,d,cache) 
#@btime fastArith(a,b,c,d,cache) 


#########stress testing exhaustive testing#########


function boundaryTestSub(a::Vector{Taylor1{Float64}},b::Vector{Taylor1{Float64}},c::Vector{Taylor1{Float64}},d::Vector{Taylor1{Float64}},cache::Vector{Taylor1{Float64}})
 # @assert 0.0-a[1]==@changeAST 0.0-a[1],12
 #@assert 0.0-a[1]-1.0==@changeAST 0.0-a[1]-1.0,12
 # @assert 1.2-2.1-a[1]-2.04-b[2]-3.2-c[1]-d[1]-1.12-3.2-d[2]==@changeAST 1.2-2.1-a[1]-2.04-b[2]-3.2-c[1]-d[1]-1.12-3.2-d[2],12
end
#boundaryTestSub(a,b,c,d,cache)
function benchmarkboundarySub(a::Vector{Taylor1{Float64}},b::Vector{Taylor1{Float64}},c::Vector{Taylor1{Float64}},d::Vector{Taylor1{Float64}},cache::Vector{Taylor1{Float64}})
#@changeAST 0.0-a[1],12 #15.648 ns (0 allocations: 0 bytes)
#0.0-a[1] #44.508 ns (1 allocation: 80 bytes)
#@changeAST 0.0-a[1]-1.0,12 #16.099 ns (0 allocations: 0 bytes)
#0.0-a[1]-1.0 #80.048 ns (2 allocations: 160 bytes)
 #@changeAST 1.2-2.1-a[1]-2.04-b[2]-3.2-c[1]-d[1]-1.12-3.2-d[2],12 #89.111 ns (0 allocations: 0 bytes)
# 1.2-2.1-a[1]-2.04-b[2]-3.2-c[1]-d[1]-1.12-3.2-d[2] #394.097 ns (9 allocations: 720 bytes)
end
#@btime benchmarkboundarySub(a,b,c,d,cache)
function closeEnough(sum1::Taylor1{Float64},sum2::Taylor1{Float64}) 
  isclose=true
  for i=0:length(sum1.coeffs)-1
    if abs(sum1[i]-sum2[i])>1e-12
      isclose=false
      break
    end
  end
  return isclose
end
function closeEnough(sum1::Float64,sum2::Taylor1{Float64}) 
  return abs(sum1-sum2[0])<1e-12
end
function stressTestSub(a::Vector{Taylor1{Float64}},b::Vector{Taylor1{Float64}},c::Vector{Taylor1{Float64}},d::Vector{Taylor1{Float64}},cache::Vector{Taylor1{Float64}})
  #test 2 vars
  for term1 in(1.0,a[1])
    for term2 in(1.0,b[1])
     if term1==1.0 && term2==1.0  
      @assert term1-term2==(@changeAST term1-term2,12)[0]
     else
        @assert term1-term2==@changeAST term1-term2,12
     end
    end
  end
  #test 3 vars
  for term1 in(d[1],1.2,a[1])
      for term2 in(1.33,b[1])
        for term3 in(2.07,c[1]) 
          if term1==1.2 && term2==1.33 && term3==2.07
            clearCache(cache)
            @assert term1-term2-term3==(@changeAST term1-term2-term3,12)[0]
           # @show cache
          else
            clearCache(cache)

            @assert closeEnough(term1-term2-term3,@changeAST term1-term2-term3,12)
             
          end
        end
      end
  end
    #test 4 vars
    for term1 in(d[1],1.2,a[1])
      for term2 in(1.33,b[1])
        for term3 in(2.07,c[1]) 
          for term4 in(1.07,d[1]) 
              if term1==1.2 && term2==1.33 && term3==2.07 && term4==1.07
                clearCache(cache)
                @assert term1-term2-term3-term4==(@changeAST term1-term2-term3-term4,12)[0]
              # @show cache
              else
                clearCache(cache)
                @assert closeEnough(term1-term2-term3-term4,@changeAST term1-term2-term3-term4,12)
              end
           end
        end
      end
  end
      #test 5 vars
      for term1 in(d[1],1.2,a[1])
        for term2 in(1.33,b[1])
          for term3 in(2.07,c[1]) 
            for term4 in(1.07,d[1]) 
              for term5 in(0.7,a[2])
                  if term1==1.2 && term2==1.33 && term3==2.07 && term4==1.07 && term5==0.7
                    clearCache(cache)
                    @assert term1-term2-term3-term4-term5==(@changeAST term1-term2-term3-term4-term5,12)[0]
                  # @show cache
                  else
                    clearCache(cache)
                    @assert closeEnough(term1-term2-term3-term4-term5,@changeAST term1-term2-term3-term4-term5,12)
                  end
                end
             end
          end
        end
     end
  
end
#a deep investigation may avoid the need to clean cache
function clearCache(cache::Vector{Taylor1{T}}) where {T<:Number}
    for i=1:length(cache)
      cache[i].coeffs.=0.0
    end
end
stressTestSub(a,b,c,d,cache)





