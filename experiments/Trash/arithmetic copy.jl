
#= using TaylorSeries
#using InteractiveUtils
#using BenchmarkTools
using MacroTools: prewalk, postwalk, @capture
import Base.:-
import Base.:+
use_show_default(true) =#
#add functionality to arithmetic.jl inside TSfor (f, fc) in ((:+, :(add!)), (:-, :(subst!)))
## Addition and substraction ##
#this section until the macro should be for non-taylor types (floats)
#all taylor ops produced by the macro should be caught in the taylor section below
#(-)(a, b, c)=(-)((-)(a,b),c)
#(+)(a, b, c)=(+)((+)(a,b),c)
(addsub)(a, b, c)=(-)((+)(a,b),c)
(subsub)(a, b, c)=(-)((-)(a,b),c)
#(addsub)(cache,a, b, c)=(-)((+)(cache,a,b),c)# cache is at begining unlike the expression 
                                          #because taylor func flipped them: oldRes becomes at begining of exprss which makes sense.
(subadd)(a, b, c)=(+)((-)(a,b),c)# this should not exist cuz addsub has 3 args + cache
#(subadd)(cache,a, b, c)=(+)((-)(cache,a,b),c) 
#= afoldl(op, a) = a # afoldl  required or import it ..later can cause redefine warning 
function afoldl(op, a, bs...)
    l = length(bs)
  #  @show l
    i =  0; y = a;            l == i && return y
    #@nexprs 31 i -> (y = op(y, bs[i]); l == i && return y)
    i =  1; y = op(y, bs[i]); l == i && return y
    i =  2; y = op(y, bs[i]); l == i && return y
    i =  3; y = op(y, bs[i]); l == i && return y
    i =  4; y = op(y, bs[i]); l == i && return y
end
  @eval begin
      (-)(a, b, c, xs...) = afoldl(-, (-)((-)(a,b),c), xs...) # should never be here since julia ast does not produce
      (+)(a, b, c, xs...) = afoldl(+, (+)((+)(a,b),c), xs...)
  end =#

 #=  macro changeAST(ex,cacheT)
    Base.remove_linenums!(ex)
    # dump(cacheT; maxdepth=8)
    twoInOne(ex)
   esc(twoInOne2(ex,cacheT))# adds cacheT to each call op
  end =#
  macro changeAST(ex)
    Base.remove_linenums!(ex)
    # dump(ex; maxdepth=8)
    twoInOne(ex)
    twoInOne2(ex)# adds cacheT to each call op
    esc(ex.args[1])# return  only
  end
  function twoInOne2(ex)
    prewalk(ex) do x
     # @show x
     if x isa Expr && x.head == :call && (x.args[1]==:- ) && length(x.args)>=3 
      push!(x.args,ex.args[2])  
      x.args[1]=:subT  # symbol changed cuz avoid type taylor piracy
      #push!(x.args,ex.args[1].args[2])# use this when exprs used for testing
    elseif x isa Expr && x.head == :call && (x.args[1]==:-  ) && length(x.args)==2  #negation
      push!(x.args,ex.args[2])  # add cache for real
     x.args[1]=:negateT  # symbol changed cuz avoid type taylor piracy
    elseif x isa Expr && x.head == :call && (x.args[1]==:+  ) && length(x.args)>=3 
      push!(x.args,ex.args[2])  
      x.args[1]=:addT
      #push!(x.args,ex.args[1].args[2])# use this when exprs used for testing cahce changes pos when used expression from outside
     
    elseif x isa Expr && x.head == :call && (x.args[1]==:* ) && length(x.args)>=3 
      push!(x.args,ex.args[2])  
      x.args[1]=:mulT
      #push!(x.args,ex.args[1].args[2])# use this when exprs used for testing
     end
    
     return x
    end
  end
  function twoInOne(ex)
    prewalk(ex) do x
     # @show x
     if x isa Expr && x.head == :call && x.args[1]==:- && length(x.args)==3 
       # @show x.args[1]
        if x.args[2] isa Expr && x.args[2].head == :call && x.args[2].args[1]==:- && length(x.args[2].args)==3 
          push!(x.args,x.args[3])#  args[4]=c
          x.args[3]=x.args[2].args[3]# args[3]=b
          x.args[2]=x.args[2].args[2]# args[2]=a
          x.args[1]=:subsub
          push!(x.args,ex.args[2])  # add cache
        elseif x.args[2] isa Expr && x.args[2].head == :call && x.args[2].args[1]==:+ && length(x.args[2].args)==3 
          push!(x.args,x.args[3])#  args[4]=c
          x.args[3]=x.args[2].args[3]# args[3]=b
          x.args[2]=x.args[2].args[2]# args[2]=a
          x.args[1]=:addsub # £ µ § ~....
          push!(x.args,ex.args[2])
         end
     elseif x isa Expr && x.head == :call && x.args[1]==:+ && length(x.args)==3 
         if x.args[2] isa Expr && x.args[2].head == :call && x.args[2].args[1]==:- && length(x.args[2].args)==3 
            push!(x.args,x.args[3])#  args[4]=c
            x.args[3]=x.args[2].args[3]# args[3]=b
            x.args[2]=x.args[2].args[2]# args[2]=a
            x.args[1]=:subadd#:µ  # £  § ....
            push!(x.args,ex.args[2])
          end
     end
     return x
      
    end
  end



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
  cache.coeffs.=0.0
  cache[0]=a+c-b
  return cache
end

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
@__dot__ cache.coeffs = (-)(b.coeffs)
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
    cache.coeffs.=a.coeffs  
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
    cache.coeffs.=0.0
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
###########"negate###########
function negateT(a::Taylor1{T},cache::Taylor1{T}) where {T<:Number}
  @__dot__ cache.coeffs = (-)(a.coeffs)
  return cache
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
    @__dot__ cache.coeffs = a.coeffs
   cache[0]=cache[0]-b    
    #@__dot__ cache.coeffs = (-)(a.coeffs, b.coeffs)
     return cache
   end
   function subT(a::T, b::Taylor1{T},cache::Taylor1{T}) where {T<:Number}
    #cache.coeffs.=-b.coeffs  #not tested
    #cache.coeffs.=0.0
    @__dot__ cache.coeffs = (-)(b.coeffs)
   # @__dot__ cache.coeffs=(-)(cache.coeffs, b.coeffs)
   cache[0]=cache[0]+a    
    #@__dot__ cache.coeffs = (-)(a.coeffs, b.coeffs)
     return cache
   end

function addT(a::Taylor1{T}, b::Taylor1{T},cache::Taylor1{T}) where {T<:Number}
        @__dot__ cache.coeffs = (+)(a.coeffs, b.coeffs)
                return cache
end
function addT(a::Taylor1{T}, b::T,cache::Taylor1{T}) where {T<:Number}
    cache.coeffs.=a.coeffs  
    cache[0]=cache[0]+b 
          return cache
end
addT(b::T,a::Taylor1{T},cache::Taylor1{T}) where {T<:Number}=addT(a, b,cache) 


function addT(a::Taylor1{T}, b::Taylor1{T},c::Taylor1{T},cache::Taylor1{T}) where {T<:Number}
@__dot__ cache.coeffs = (+)(a.coeffs, b.coeffs,c.coeffs)
  return cache
end   
function addT(a::Taylor1{T}, b::Taylor1{T},c::T,cache::Taylor1{T}) where {T<:Number}
  cache.coeffs.=a.coeffs  #not tested
    cache[0]=cache[0]+c 
  @__dot__ cache.coeffs = (+)(cache.coeffs, b.coeffs)
    return cache
  end   
  addT(a::Taylor1{T},c::T, b::Taylor1{T},cache::Taylor1{T}) where {T<:Number}=addT(a, b,c,cache) 
  addT(c::T, a::Taylor1{T}, b::Taylor1{T},cache::Taylor1{T}) where {T<:Number}=addT(a, b,c,cache)

  function addT(a::Taylor1{T}, b::T,c::T,cache::Taylor1{T}) where {T<:Number}
    cache.coeffs.=a.coeffs  
      cache[0]=cache[0]+c+b
      return cache
    end 
    addT( b::T,a::Taylor1{T},c::T,cache::Taylor1{T}) where {T<:Number}=addT(a, b,c,cache) 
    addT( b::T,c::T,a::Taylor1{T},cache::Taylor1{T}) where {T<:Number}=addT(a, b,c,cache) 

    function addT( a::T,b::T,c::T,cache::Taylor1{T}) where {T<:Number}
      cache.coeffs.=0.0
      cache[0]=a+b+c    
       return cache
     end


#four variables still not included constants
    function addT(a::Taylor1{T}, b::Taylor1{T},c::Taylor1{T},d::Taylor1{T},cache::Taylor1{T}) where {T<:Number}
@__dot__ cache.coeffs = (+)(a.coeffs, b.coeffs,c.coeffs,d.coeffs)
  return cache
end        
function addT(a::Taylor1{T}, b::Taylor1{T},c::Taylor1{T},d::Taylor1{T},e::Taylor1{T},cache::Taylor1{T}) where {T<:Number}
  #v = similar(a.coeffs)
#=    println("ab")
  @show $f
  @show a.coeffs
  @show b.coeffs
  @show cache.coeffs =#
@__dot__ cache.coeffs = (+)(a.coeffs, b.coeffs,c.coeffs,d.coeffs,e.coeffs)
 # @show cache.coeffs
#return 0#$T(v, a.order)
  return cache
end      
function addT(a::Taylor1{T}, b::Taylor1{T},c::Taylor1{T},d::Taylor1{T},e::Taylor1{T},f::Taylor1{T},cache::Taylor1{T}) where {T<:Number}
  #v = similar(a.coeffs)
#=    println("ab")
  @show $f
  @show a.coeffs
  @show b.coeffs
  @show cache.coeffs =#
@__dot__ cache.coeffs = (+)(a.coeffs, b.coeffs,c.coeffs,d.coeffs,e.coeffs,f.coeffs)
 # @show cache.coeffs
#return 0#$T(v, a.order)
  return cache
end   


####################################################
#= a =   [ Taylor1([1.0,1.1,1.2],2),Taylor1([1.3,1.4,1.5],2)]
b =   [ Taylor1([2.0,2.1,2.2],2),Taylor1([2.3,2.4,2.5],2)]
c =   [ Taylor1([3.0,3.1,3.2],2),Taylor1([3.3,3.4,3.5],2)]
d =    [Taylor1([4.0,4.1,4.2],2),Taylor1([4.3,4.4,4.5],2)]
cacheT=Taylor1([0.0,0.0,0.0],2)
function fastArith(a::Vector{Taylor1{Float64}},b::Vector{Taylor1{Float64}},c::Vector{Taylor1{Float64}},d::Vector{Taylor1{Float64}},cacheT::Taylor1{Float64})
    #sum=@changeAST a[1]-b[1]-c[1]-d[1]+a[2]-b[2],cacheT
    sum=@changeAST a[1]+b[1]+c[1]-4.0,cacheT
    return sum
end
function normalArith(a::Vector{Taylor1{Float64}},b::Vector{Taylor1{Float64}},c::Vector{Taylor1{Float64}},d::Vector{Taylor1{Float64}},cacheT::Taylor1{Float64})
  #sum=a[1]-b[1]-c[1]-d[1]+a[2]-b[2]
  sum= a[1]+b[1]+c[1]-4.0
  return sum
end
@show normalArith(a,b,c,d,cacheT)
@btime normalArith(a,b,c,d,cacheT)  =#
#@show fastArith(a,b,c,d,cacheT)
#@btime fastArith(a,b,c,d,cacheT) 

###################results######################
#expression                            normalArith                                             fast

#a+b+c+d                        143.642 ns (3 allocations: 240 bytes)              30.696 ns (0 allocations: 0 bytes)

#a+b+c-d                        136.106 ns (3 allocations: 240 bytes)              29.332 ns (0 allocations: 0 bytes)

#a-b+c-d                         134.348 ns (3 allocations: 240 bytes)              28.232 ns (0 allocations: 0 bytes)

#a-b-c-d                         136.262 ns (3 allocations: 240 bytes)              27.714 ns (0 allocations: 0 bytes)

#a[1]-b[1]-c[1]-d[1]             142.936 ns (3 allocations: 240 bytes)             32.863 ns (0 allocations: 0 bytes)

#a[1]-b[1]-c[1]-d[1]+a[2]-b[2]  142.936 ns (3 allocations: 240 bytes)               46.704 ns (0 allocations: 0 bytes)


#a[1]+b[1]+c[1]-4.0            136.052 ns (3 allocations: 240 bytes)           32.335 ns (0 allocations: 0 bytes)


























