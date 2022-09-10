
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
     #   println("fuck add xs= ",xs[end])
        addTfoldl(addT, addT( addT(a,b,c,xs[end]) , d ,xs[end] ), xs...)
      end

    end

    function mulTfoldl(op, res, bs...)
      l = length(bs)
      
      #i =  0;  l == i && return a  # already checked
      i =  1;  l == i && return res 
      i =  2;  l == i && return              op(res, bs[1],bs[end])
      i =  3;  l == i && return           op(op(res, bs[1],bs[end]),bs[2],bs[end])
      i =  4;  l == i && return        op(op(op(res, bs[1],bs[end]),bs[2],bs[end]),bs[3],bs[end])
      if l==5
       
        return     op(op(op(op(res, bs[1],bs[end]),bs[2],bs[end]),bs[3],bs[end]),bs[4],bs[end])
      end
     # i =  5;  l == i &&  return     op(op(op(op(res, bs[1],bs[end]),bs[2],bs[end]),bs[3],bs[end]),bs[4],bs[end])
      i =  6;  l == i &&  return op(op(op(op(op(res, bs[1],bs[end]), bs[2],bs[end]), bs[3],bs[end]),bs[4],bs[end]),bs[5],bs[end])
      i =  7;  l == i &&  return op(op(op(op(op(op(res, bs[1],bs[end]), bs[2],bs[end]), bs[3],bs[end]),bs[4],bs[end]),bs[5],bs[end]),bs[6],bs[end])
      i =  8; @show i; l == i && return op(op(op(op(op(op(op(res, bs[1],bs[end]), bs[2],bs[end]), bs[3],bs[end]),bs[4],bs[end]),bs[5],bs[end]),bs[6],bs[end]),bs[7],bs[end])
      i =  9;@show i;y=op(op(op(op(op(op(op(op(res, bs[1],bs[end]), bs[2],bs[end]), bs[3],bs[end]),bs[4],bs[end]),bs[5],bs[end]),bs[6],bs[end]),bs[7],bs[end]),bs[8],bs[end]);  l == i && return y
      #last i creates y and passes it to the next loop: note that from this point, op will allocate
      # either keep increasing i manually or warn the users . (inserting -0 will avoid allocation lol!)
       for i in i:l-1 # -1 cuz last one is cache
        @show i;
          y = op(y, bs[i],bs[end])
      end
     # println("at end")
      return y 
    end
    function mulT(a, b, c, xs...)
      
      if length(xs)!=0# if ==0 case below where c==cache
     #   println("fuck multi xs= ",xs[end])
        mulTfoldl(mulT, mulT( mulT(a,b,xs[end]) , c ,xs[end] ), xs...)
      end

    end
 
  macro changeAST(ex)
    Base.remove_linenums!(ex)
    # dump(ex; maxdepth=8)
    twoInOne(ex)
   twoInOne2(ex)# adds cache to each call op
   esc(ex.args[1])# return cache only
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
      #push!(x.args,ex.args[1].args[2])# use this when exprs used for testing
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
         # @show x.args
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
  cache.coeffs.=0.0   # since C is number this cant be an op in the middle (all in middle op return taylor)
                      # still the cache can be full from another completely diff equation
                      # emptying the cache does not affect c (in this case c is not cache)
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
    cache.coeffs.=0.0 # a is number, ok to empty cache
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
    @__dot__ cache.coeffs = a.coeffs
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

function addT(a::Taylor1{T}, b::Taylor1{T},cache::Taylor1{T}) where {T<:Number}
  
        @__dot__ cache.coeffs = (+)(a.coeffs, b.coeffs)
                return cache
end
function addT(a::Taylor1{T}, b::T,cache::Taylor1{T}) where {T<:Number}
  @show a==cache
    cache.coeffs.=a.coeffs  
    cache[0]=cache[0]+b 
          return cache
end
addT(b::T,a::Taylor1{T},cache::Taylor1{T}) where {T<:Number}=addT(a, b,cache) 

function addT( a::T,b::T,cache::Taylor1{T}) where {T<:Number}
  cache.coeffs.=0.0 # a is number ok to empty
  cache[0]=a+b   
  @show cache
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
      cache.coeffs.=0.0 #a is number ok to empty
      cache[0]=a+b+c    
       return cache
     end


#= #four variables still not included constants
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
 =#

#########################mul###########################
#########################mul###########################

#= function mulT(a::T, b::Taylor1{T},cache::Vector{Taylor1{T}}) where {T<:Number}
  
     @__dot__ cache[1].coeffs = a * b.coeffs #broadcast dimension mismatch
  return cache[1]
end =#
function mulT(a::Taylor1{T}, b::T,cache::Vector{Taylor1{T}}) where {T<:Number}
  
  #@__dot__ cache[1].coeffs = a.coeffs * b  ##broadcast dimension mismatch
  for i=0:length(a)-1
    cache[1][i]=a[i]*b
  end
  return cache[1]
end
mulT(a::T,b::Taylor1{T}, cache::Vector{Taylor1{T}}) where {T<:Number} = mulT(b , a,cache)

function mulT(a::Taylor1{T}, b::Taylor1{T},cache::Vector{Taylor1{T}}) where {T<:Number}
  #println(" b= ",b)
 #=  if a.order != b.order
      a, b = fixorder(a, b)
  end =#
  n=length(a.coeffs)
v=SVector{3,Float64}(a.coeffs)  # float64 can be changed later to T
 # println("taylor *taylor")
  
 # cache[2].coeffs.=0.0  # !!!!!!!!!bug
 # println(" cache= ",cache)
 # @show a===cache[1]   # true this is bad since modifying cache[1] will modify a in below loop !!!!
#@show a===cache[2]
 # println("after empty a= ",a)
 #println("isbits v= ",isbits(v))
 for k in eachindex(a)
    #println("inloop a= ",a)
   # 
      @inbounds cache[1][k] = v[1] * b[k]
      @inbounds for i = 1:k
        cache[1][k] += v[i+1] * b[k-i]
      end
     # println("inloop cache= ",cache)
  end
 # println("2Taylors cache[1]= ",cache[1])
# @__dot__ cache[1].coeffs = cache[2].coeffs
 #println(" cache= ",cache)
 #println(" cache1===cache2 ",cache[1]===cache[2])

 #even i this function has something simple like

 # cache[1][1] = b[1] #* b[k-i]

  return cache[1]
end
function mulT(a::T, b::T,cache::Vector{Taylor1{T}}) where {T<:Number}
 # println("number *number to provide taylor")
  #cache[1].coeffs.=0.0
  fill!(cache[1].coeffs, 0.0)
  cache[1][0]=a*b   
       return cache[1]
end






a =   [ Taylor1([1.0,1.1,1.2],2),Taylor1([1.3,1.4,1.5],2)]
b =   [ Taylor1([2.0,2.1,2.2],2),Taylor1([2.3,2.4,2.5],2)]
c =   [ Taylor1([3.0,3.1,3.2],2),Taylor1([3.3,3.4,3.5],2)]
d =    [Taylor1([4.0,4.1,4.2],2),Taylor1([4.3,4.4,4.5],2)]
cache=[Taylor1([0.0,0.0,0.0],2),Taylor1([0.0,0.0,0.0],2)]

function fastArith(a::Vector{Taylor1{Float64}},b::Vector{Taylor1{Float64}},c::Vector{Taylor1{Float64}},d::Vector{Taylor1{Float64}},cache::Vector{Taylor1{Float64}})
  
  # @show 3.2*1.3
  #=  sum=@changeAST a[1]*1.3,cacheT
   @assert sum==a[1]*1.3
   sum=@changeAST 1.3*a[1],cacheT
   @assert sum==1.3*a[1]
   sum=@changeAST a[1]*b[1],cacheT
   @assert sum==a[1]*b[1]
   sum=@changeAST 3.2*1.3,cacheT#-d[1]+a[2]-b[2],cacheT
   #@show sum
   sum=@changeAST a[1]*1.3*2.0,cacheT
   @assert sum==a[1]*1.3*2.0
   sum=@changeAST 1.3*a[1]*2.0,cacheT
   @assert sum==1.3*a[1]*2.0
   sum=@changeAST a[1]*b[1]*2.0,cacheT
   @assert sum==a[1]*b[1]*2.0 =#
   
  # sum=@changeAST a[1]*b[2]*b[1]*c[1],cacheT
 #sum=@changeAST a[1]+1.0+a[2]+1.0+1.0+b[1]+b[2]+c[2],cacheT
# sum=@changeAST a[1]*3.2*0.9*a[2],cacheT #i=2                  38.835 ns (0 allocations: 0 bytes)
 #sum=@changeAST a[1]*3.2*c[1]*0.9*a[2],cacheT #i=3             56.858 ns (0 allocations: 0 bytes)
# sum=@changeAST a[1]*3.2*1.3*c[1]*0.9*a[2],cacheT #i=4        70.905 ns (0 allocations: 0 bytes)
 sum=@changeAST 2.2*a[1]*c[2]*d[1]*a[2]*2.0*c[1]*2.0,cache #i=5    601.403 ns (16 allocations: 528 bytes)
 
 # sum=@changeAST 2.2-1.1-3.0,cacheT
  return sum
end
function normalArith(a::Vector{Taylor1{Float64}},b::Vector{Taylor1{Float64}},c::Vector{Taylor1{Float64}},d::Vector{Taylor1{Float64}},cacheT::Taylor1{Float64})
  #sum=a[1]-b[1]-c[1]-d[1]+a[2]-b[2]
 # sum= 1.1- a[1]-b[1]
  #sum=a[1]+b[1]+c[1]+d[1]+a[2]-b[2]+c[2]
#  sum=a[1]*b[2]*b[1]*c[1]
#=   sum=3.2*1.3
  sum=a[1]*1.3
  sum=1.3*a[1]
  sum=a[1]*b[1] =#
  sum=2.2*a[1]*c[2]*d[1]*a[2]*2.0*c[1]*2.0
  return sum
end
 #@show normalArith(a,b,c,d,cache[1])
#@btime normalArith(a,b,c,d,cache)  
#@show fastArith(a,b,c,d,cache)
@btime fastArith(a,b,c,d,cache) 







function showFastTaylor(x::Taylor1{Float64},T::Taylor1{Float64},cacheT::Taylor1{Float64})
  term2=1.0
  term1=1.0
  #for term1 in(1.0,T)
    #for term2 in(1.0,T)
     # for term3 in(1.0,T)   
        #@eval begin
          @show term1 
          @show term2 
        
          sum1=@changeAST term1 + term2 ,cacheT
          @show sum1
         sum2=term1 + term2 
        # @show sum2
           sum1===sum2 || @show sum2 
           
           term1=T
           term2=x
           sum1=@changeAST term1 + term2 ,cacheT
         sum2=term1 + term2 
        # @show sum2
           sum1==sum2 || @show sum2 

       # end
   #  end
   # end
 # end
  end
  r=Taylor1([1.1,2.2,3.3])
 # showFastTaylor(a[1],r,cache)