
using TaylorSeries
using StaticArrays
using InteractiveUtils
using TimerOutputs
using BenchmarkTools
using MacroTools: prewalk, postwalk, @capture
import Base.:-
import Base.:+
use_show_default(true)

(addsub)(a, b, c)=(-)((+)(a,b),c)
(subsub)(a, b, c)=(-)((-)(a,b),c)
(subadd)(a, b, c)=(+)((-)(a,b),c)# this should not exist cuz addsub has 3 args + cache[1]

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
       for i in i:l-1 # -1 cuz last one is cache[1]
          y = op(y, bs[i],bs[end])
      end
      return y 
  end
    function addT(a, b, c,d, xs...)
      if length(xs)!=0# if ==0 case below where d==cache[1]
     
        addTfoldl(addT, addT( addT(a,b,c,xs[end]) , d ,xs[end] ), xs...)
      end

    end

    function mulTfoldl(op, res, bs...)
      l = length(bs)
    #  @show l
      #println("entered multfoldl")
      #i =  0;  l == i && return a  # already checked
      i =  1;  l == i && return res 
      i =  2;  l == i && return              op(res, bs[1],bs[end])
      i =  3; l == i && return           op(op(res, bs[1],bs[end]),bs[2],bs[end])
      i =  4;  l == i && return        op(op(op(res, bs[1],bs[end]),bs[2],bs[end]),bs[3],bs[end])
      #= if l==5
       
      # @show bs
        return     op(op(op(op(res, bs[1],bs[end]),bs[2],bs[end]),bs[3],bs[end]),bs[4],bs[end])
      end =#

      #since bs[end][1](cache1 has previous result then i can use it directly op(cache1,b,cache)) shortlines :)
      i =  5;  l == i &&  return     op(op(op(op(res, bs[1],bs[end]),bs[2],bs[end]),bs[3],bs[end]),bs[4],bs[end])
      i =  6; @show i; l == i &&  return op(op(op(op(op(res, bs[1],bs[end]), bs[2],bs[end]), bs[3],bs[end]),bs[4],bs[end]),bs[5],bs[end])
      i =  7; @show i; l == i &&  return op(op(op(op(op(op(res, bs[1],bs[end]), bs[2],bs[end]), bs[3],bs[end]),bs[4],bs[end]),bs[5],bs[end]),bs[6],bs[end])
      i =  8; @show i; l == i && return op(op(op(op(op(op(op(res, bs[1],bs[end]), bs[2],bs[end]), bs[3],bs[end]),bs[4],bs[end]),bs[5],bs[end]),bs[6],bs[end]),bs[7],bs[end])
      i =  9;@show i;y=op(op(op(op(op(op(op(op(res, bs[1],bs[end]), bs[2],bs[end]), bs[3],bs[end]),bs[4],bs[end]),bs[5],bs[end]),bs[6],bs[end]),bs[7],bs[end]),bs[8],bs[end]);  l == i && return y
      #last i creates y and passes it to the next loop: note that from this point, op will allocate
      # either keep increasing i manually or warn the users . (inserting -0 will avoid allocation lol!)
       for i in i:l-1 # -1 cuz last one is cache[1]
        @show i;
          y = op(y, bs[i],bs[end])
      end
     # println("at end")
      return y 
    end
    function mulT(a, b, c, xs...)
     # println("entered mult(a b c)")
      if length(xs)!=0# if ==0 case below where c==cache[1]
     #   println("fuck multi xs= ",xs[end])
       @timeit "mulfold" mulTfoldl(mulT, mulT( mulT(a,b,xs[end]) , c ,xs[end] ), xs...)
     # @timeit "4 terms"  mulT(mulT( mulT(a,b,xs[end]) , c ,xs[end] ),xs[1],xs[end] )
      end

    end
  macro changeAST(ex)
    Base.remove_linenums!(ex)
    # dump(ex; maxdepth=8)
    twoInOne(ex)
   twoInOne2(ex)# adds cache[1] to each call op
  # @show ex.args[1]
   esc(ex.args[1])# return 
  end
  function twoInOne2(ex)
    prewalk(ex) do x
     # @show x
     if x isa Expr && x.head == :call && (x.args[1]==:- ) && length(x.args)>=3 
      push!(x.args,ex.args[2])  
      x.args[1]=:subT  # symbol changed cuz avoid type taylor piracy
      #push!(x.args,ex.args[1].args[2])# use this when exprs used for testing
    elseif x isa Expr && x.head == :call && (x.args[1]==:-  ) && length(x.args)==2  #negation
       push!(x.args,ex.args[2])  # add cache[1] for real
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
          push!(x.args,ex.args[2])  # add cache[1]
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
  function addsub(a::Taylor1{T}, b::Taylor1{T},c::Taylor1{T},cache::Vector{Taylor1{T}}) where {T<:Number}
    #v = similar(a.coeffs)
    
    #@show cache[1]
   # @__dot__ cache[1].coeffs = addsub(cache[1].coeffs,a.coeffs, b.coeffs,c.coeffs)
   @__dot__ cache[1].coeffs = addsub(a.coeffs, b.coeffs,c.coeffs)
    return cache[1]
  end
  function addsub(a::Taylor1{T}, b::Taylor1{T},c::T,cache::Vector{Taylor1{T}}) where {T<:Number}
    #v = similar(a.coeffs)
    
    #@show cache[1]
   # @__dot__ cache[1].coeffs = addsub(cache[1].coeffs,a.coeffs, b.coeffs,c.coeffs)
            cache[1].coeffs.=b.coeffs  
            cache[1][0]=cache[1][0]-c
   @__dot__ cache[1].coeffs = (+)(cache[1].coeffs, a.coeffs)
    #return 0#$T(v, a.order)
    return cache[1]
  end
  function addsub(a::T, b::Taylor1{T},c::Taylor1{T},cache::Vector{Taylor1{T}}) where {T<:Number}
            cache[1].coeffs.=b.coeffs  
            cache[1][0]=a+ cache[1][0]
           # @show a
   @__dot__ cache[1].coeffs = (-)(cache[1].coeffs, c.coeffs)
    return cache[1]
  end
  addsub(a::Taylor1{T}, b::T,c::Taylor1{T},cache::Vector{Taylor1{T}}) where {T<:Number}=addsub(b, a ,c,cache)#addsub(b::T, a::Taylor1{T} ,c::Taylor1{T},cache::Vector{Taylor1{T}}) where {T<:Number}
 
  function addsub(a::T, b::T,c::Taylor1{T},cache::Vector{Taylor1{T}}) where {T<:Number}
    #cache[1].coeffs.=-c.coeffs #not tested
    #cache[1].coeffs.=0.0
    @__dot__ cache[1].coeffs = (-)(c.coeffs)
   # @__dot__ cache[1].coeffs=(-)(cache[1].coeffs, c.coeffs)#   ==-c
    cache[1][0]=a+b+ cache[1][0] 
    return cache[1]
end



function addsub(c::Taylor1{T},a::T, b::T,cache::Vector{Taylor1{T}}) where {T<:Number}
  cache[1].coeffs.=c.coeffs 
  cache[1][0]=a+ cache[1][0]-b
  return cache[1]
end
 addsub(a::T, c::Taylor1{T},b::T,cache::Vector{Taylor1{T}}) where {T<:Number}= addsub(c,a, b,cache) 
 
 function addsub(c::T,a::T, b::T,cache::Vector{Taylor1{T}}) where {T<:Number}
  ##[[[[[[[[[!!!!!!!!!!!!!!!!!]]]]]]]]]cache[1][1]...should be emptied
  cache[1].coeffs.=0.0   # since C is number this cant be an op in the middle (all in middle op return taylor)
                      # still the cache[1] can be full from another completely diff equation
                      # emptying the cache[1] does not affect c (in this case c is not cache[1])
  cache[1][0]=a+c-b
  return cache[1]
end

###########"negate###########
function negateT(a::Taylor1{T},cache::Vector{Taylor1{T}}) where {T<:Number}
  @__dot__ cache[1].coeffs = (-)(a.coeffs)
  return cache[1]
end
#= function negateT(b::T,cache::Vector{Taylor1{T}}) where {T<:Number} # no need since ast of -Number is direct
  cache[1].coeffs.=0.0
  cache[1][0]=-b
  return cache[1]
end =#
#################################################subsub########################################################""


  function subsub(a::Taylor1{T}, b::Taylor1{T},c::Taylor1{T},cache::Vector{Taylor1{T}}) where {T<:Number}
    
    @__dot__ cache[1].coeffs = subsub(a.coeffs, b.coeffs,c.coeffs)
    
    return cache[1]
  end

  function subsub(a::Taylor1{T}, b::Taylor1{T},c::T,cache::Vector{Taylor1{T}}) where {T<:Number}
    #@show a
   cache[1].coeffs.=b.coeffs  
   cache[1][0]=cache[1][0]+c     #(a-b-c=a-(b+c))
  @__dot__ cache[1].coeffs = (-)(a.coeffs, cache[1].coeffs)
    return cache[1]
  end
  subsub(a::Taylor1{T}, b::T,c::Taylor1{T},cache::Vector{Taylor1{T}}) where {T<:Number}=subsub(a, c,b,cache) 
   #------!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!----------------
 #


 #      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 #--------!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-------
  function subsub(a::T,b::Taylor1{T}, c::Taylor1{T},cache::Vector{Taylor1{T}}) where {T<:Number}
   # @show c
# bug: this is the first time cache[1] is used upfront..!!!cache[1] not empty
#cache[1].coeffs.=0.0
@__dot__ cache[1].coeffs = (-)(b.coeffs) # like before since a is number modifying the cache[1] is fine
  # @__dot__ cache[1].coeffs=(-)(cache[1].coeffs, b.coeffs)
   cache[1][0]=a+ cache[1][0]
  # @show a
@__dot__ cache[1].coeffs = (-)(cache[1].coeffs, c.coeffs)
return cache[1]



   #=  cache[1].coeffs.=-b.coeffs #not tested
    cache[1][0]=a+cache[1][0]     #(a-b-c=a-(b+c))
   @__dot__ cache[1].coeffs = (-)(cache[1].coeffs, c.coeffs)
     return cache[1] =#
   end


   function subsub(a::Taylor1{T}, b::T,c::T,cache::Vector{Taylor1{T}}) where {T<:Number}
    cache[1].coeffs.=a.coeffs  
    cache[1][0]=cache[1][0]-b-c     
     return cache[1]
   end

   function subsub( a::T,b::Taylor1{T},c::T,cache::Vector{Taylor1{T}}) where {T<:Number}
    #cache[1].coeffs.=-a.coeffs  #not tested
   # cache[1].coeffs.=0.0 #]]]]]]]]]]]]]]]]]]]]]]$$$$$$$$$$$$$$$$********************
    @__dot__ cache[1].coeffs = (-)(b.coeffs)
   # @__dot__ cache[1].coeffs=(-)(cache[1].coeffs, a.coeffs)
    cache[1][0]=cache[1][0]+a-c     
     return cache[1]
   end
   subsub( b::T,c::T,a::Taylor1{T},cache::Vector{Taylor1{T}}) where {T<:Number}=subsub( b,a,c,cache) 


   function subsub( a::T,b::T,c::T,cache::Vector{Taylor1{T}}) where {T<:Number}
    cache[1].coeffs.=0.0 # a is number, ok to empty cache[1]
    cache[1][0]=a-b-c    
     return cache[1]
   end


   #################################subadd####################################""
  #= function subadd(a::Taylor1{T}, b::Taylor1{T},c::Taylor1{T},cache::Vector{Taylor1{T}}) where {T<:Number}

   @__dot__ cache[1].coeffs = subadd(a.coeffs, b.coeffs,c.coeffs)
    #return 0#$T(v, a.order)
    return cache[1]
  end =#

  #subadd(a,b,c)::Taylor1{T} where {T<:Number}=addsub(a,c,b)   #not tested
  function subadd(a::P,b::Q,c::R,cache::Vector{Taylor1{T}}) where {P,Q,R <:Union{Taylor1,Number},T<:Number}
    addsub(a,c,b,cache[1])  
end

  ##################################""subT################################
  function subT(a::Taylor1{T}, b::Taylor1{T},cache::Vector{Taylor1{T}}) where {T<:Number}
   @__dot__ cache[1].coeffs = (-)(a.coeffs, b.coeffs)
    return cache[1]
  end
  function subT(a::Taylor1{T}, b::T,cache::Vector{Taylor1{T}}) where {T<:Number}
    #= for i=1:length(cache[1].coeffs)
        cache[1].coeffs[i]=a.coeffs[i] 
    end =#
    @__dot__ cache[1].coeffs = a.coeffs
   cache[1][0]=cache[1][0]-b    
    #@__dot__ cache[1].coeffs = (-)(a.coeffs, b.coeffs)
     return cache[1]
   end
   function subT(a::T, b::Taylor1{T},cache::Vector{Taylor1{T}}) where {T<:Number}
    #cache[1].coeffs.=-b.coeffs  #not tested
    #cache[1].coeffs.=0.0
    @__dot__ cache[1].coeffs = (-)(b.coeffs) # a is number ok to modify cache[1]
   # @__dot__ cache[1].coeffs=(-)(cache[1].coeffs, b.coeffs)
   cache[1][0]=cache[1][0]+a    
    #@__dot__ cache[1].coeffs = (-)(a.coeffs, b.coeffs)
     return cache[1]
   end

function addT(a::Taylor1{T}, b::Taylor1{T},cache::Vector{Taylor1{T}}) where {T<:Number}
  
        @__dot__ cache[1].coeffs = (+)(a.coeffs, b.coeffs)
                return cache[1]
end
function addT(a::Taylor1{T}, b::T,cache::Vector{Taylor1{T}}) where {T<:Number}
  #@show a==cache[1]
    cache[1].coeffs.=a.coeffs  
    cache[1][0]=cache[1][0]+b 
          return cache[1]
end
addT(b::T,a::Taylor1{T},cache::Vector{Taylor1{T}}) where {T<:Number}=addT(a, b,cache) 

function addT( a::T,b::T,cache::Vector{Taylor1{T}}) where {T<:Number}
  cache[1].coeffs.=0.0 # a is number ok to empty
  cache[1][0]=a+b   
  @show cache[1]
   return cache[1]
 end

 #add Three vars a b c
function addT(a::Taylor1{T}, b::Taylor1{T},c::Taylor1{T},cache::Vector{Taylor1{T}}) where {T<:Number}
 # println("add abc a=",a)
 # println("add abc cache[1]=",cache[1])
  
 # println("equal cache[1]=",cache[1]==a)

@__dot__ cache[1].coeffs = (+)(a.coeffs, b.coeffs,c.coeffs)
  return cache[1]
end   
function addT(a::Taylor1{T}, b::Taylor1{T},c::T,cache::Vector{Taylor1{T}}) where {T<:Number}
  
  cache[1].coeffs.=a.coeffs  #not tested
    cache[1][0]=cache[1][0]+c 
  @__dot__ cache[1].coeffs = (+)(cache[1].coeffs, b.coeffs)
    return cache[1]
  end   
  addT(a::Taylor1{T},c::T, b::Taylor1{T},cache::Vector{Taylor1{T}}) where {T<:Number}=addT(a, b,c,cache) 
  addT(c::T, a::Taylor1{T}, b::Taylor1{T},cache::Vector{Taylor1{T}}) where {T<:Number}=addT(a, b,c,cache)

  function addT(a::Taylor1{T}, b::T,c::T,cache::Vector{Taylor1{T}}) where {T<:Number}
  
    cache[1].coeffs.=a.coeffs  
      cache[1][0]=cache[1][0]+c+b
      return cache[1]
    end 
    addT( b::T,a::Taylor1{T},c::T,cache::Vector{Taylor1{T}}) where {T<:Number}=addT(a, b,c,cache) 
    addT( b::T,c::T,a::Taylor1{T},cache::Vector{Taylor1{T}}) where {T<:Number}=addT(a, b,c,cache) 

    function addT( a::T,b::T,c::T,cache::Vector{Taylor1{T}}) where {T<:Number}
      cache[1].coeffs.=0.0 #a is number ok to empty
      cache[1][0]=a+b+c    
       return cache[1]
     end


#= #four variables still not included constants
    function addT(a::Taylor1{T}, b::Taylor1{T},c::Taylor1{T},d::Taylor1{T},cache::Vector{Taylor1{T}}) where {T<:Number}
@__dot__ cache[1].coeffs = (+)(a.coeffs, b.coeffs,c.coeffs,d.coeffs)
  return cache[1]
end        
function addT(a::Taylor1{T}, b::Taylor1{T},c::Taylor1{T},d::Taylor1{T},e::Taylor1{T},cache::Vector{Taylor1{T}}) where {T<:Number}
  #v = similar(a.coeffs)
#=    println("ab")
  @show $f
  @show a.coeffs
  @show b.coeffs
  @show cache[1].coeffs =#
@__dot__ cache[1].coeffs = (+)(a.coeffs, b.coeffs,c.coeffs,d.coeffs,e.coeffs)
 # @show cache[1].coeffs
#return 0#$T(v, a.order)
  return cache[1]
end      
function addT(a::Taylor1{T}, b::Taylor1{T},c::Taylor1{T},d::Taylor1{T},e::Taylor1{T},f::Taylor1{T},cache::Vector{Taylor1{T}}) where {T<:Number}
  #v = similar(a.coeffs)
#=    println("ab")
  @show $f
  @show a.coeffs
  @show b.coeffs
  @show cache[1].coeffs =#
@__dot__ cache[1].coeffs = (+)(a.coeffs, b.coeffs,c.coeffs,d.coeffs,e.coeffs,f.coeffs)
 # @show cache[1].coeffs
#return 0#$T(v, a.order)
  return cache[1]
end   
 =#

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
  #= n=length(a.coeffs)
v=SVector{n,Float64}(a.coeffs) =#  # float64 can be changed later to T
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
      @inbounds cache[2][k] = a[0] * b[k]
      @inbounds for i = 1:k
        cache[2][k] += a[i] * b[k-i]
      end
     # println("inloop cache= ",cache)
  end
 # println("2Taylors cache[1]= ",cache[1])
 @__dot__ cache[1].coeffs = cache[2].coeffs
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
  reset_timer!()
  # @show 3.2*1.3
  #=  sum=@changeAST a[1]*1.3,cache[1]
   @assert sum==a[1]*1.3
   sum=@changeAST 1.3*a[1],cache[1]
   @assert sum==1.3*a[1]
   sum=@changeAST a[1]*b[1],cache[1]
   @assert sum==a[1]*b[1]
   sum=@changeAST 3.2*1.3,cache[1]#-d[1]+a[2]-b[2],cache[1]
   #@show sum
   sum=@changeAST a[1]*1.3*2.0,cache[1]
   @assert sum==a[1]*1.3*2.0
   sum=@changeAST 1.3*a[1]*2.0,cache[1]
   @assert sum==1.3*a[1]*2.0
   sum=@changeAST a[1]*b[1]*2.0,cache[1]
   @assert sum==a[1]*b[1]*2.0 =#
   #sum=@changeAST b[1]*a[1],cache
  # sum=@changeAST a[1]*b[2]*b[1]*c[1],cache[1]
 #sum=@changeAST a[1]+1.0+a[2]+1.0+1.0+b[1]+b[2]+c[2],cache[1]
# sum=@changeAST a[1]*3.2*0.9*a[2],cache[1] #i=2                  38.835 ns (0 allocations: 0 bytes)
 #sum=@changeAST a[1]*3.2*c[1]*0.9*a[2],cache[1] #i=3             56.858 ns (0 allocations: 0 bytes)
 #sum=@changeAST a[1]*3.2*1.3*c[1]*0.9*a[2],cache #i=4        70.905 ns (0 allocations: 0 bytes)
 #sum=@changeAST a[1]*c[2]*d[1]*a[2]*c[1]*b[1]*b[2],cache #i=5    601.403 ns (16 allocations: 528 bytes)
 #sum=@changeAST 3.2*1.3*0.9*2.1*3.2*2.1*1.6,cache
sum=@changeAST 1.0*a[1]*2.5*1.0,cache
 #sum=@changeAST a[1]+c[2]+d[1]+2.0+c[1]+2.0+1.1+1.1+1.1+1.1+a[2],cache
 # sum=@changeAST 2.2-1.1-3.0,cache[1]
 #sum=@changeAST 3.2+1.3+1.2+a[1],cache
 print_timer()
  return sum
end
function normalArith(a::Vector{Taylor1{Float64}},b::Vector{Taylor1{Float64}},c::Vector{Taylor1{Float64}},d::Vector{Taylor1{Float64}},cache::Vector{Taylor1{Float64}})
  #sum=a[1]-b[1]-c[1]-d[1]+a[2]-b[2]
 # sum= 1.1- a[1]-b[1]
  #sum=a[1]+b[1]+c[1]+d[1]+a[2]-b[2]+c[2]
#  sum=a[1]*b[2]*b[1]*c[1]
#=   sum=3.2*1.3
  sum=a[1]*1.3
  sum=1.3*a[1]
  sum=a[1]*b[1] =#
  sum=1.0*a[1]*2.5*1.0
  #sum=b[1]*a[1]
  return sum
end
 #@show normalArith(a,b,c,d,cache)
#@btime normalArith(a,b,c,d,cache)  
#@show 
#@show 
fastArith(a,b,c,d,cache)
#@btime fastArith(a,b,c,d,cache) 




















#= 
function showFastTaylor(x::Taylor1{Float64},T::Taylor1{Float64},cache::Vector{Taylor1{Float64}})
  term2=1.0
  term1=1.0
  #for term1 in(1.0,T)
    #for term2 in(1.0,T)
     # for term3 in(1.0,T)   
        #@eval begin
          @show term1 
          @show term2 
        
          sum1=@changeAST term1 + term2 ,cache[1]
          @show sum1
         sum2=term1 + term2 
        # @show sum2
           sum1===sum2 || @show sum2 
           
           term1=T
           term2=x
           sum1=@changeAST term1 + term2 ,cache[1]
         sum2=term1 + term2 
        # @show sum2
           sum1==sum2 || @show sum2 

       # end
   #  end
   # end
 # end
  end
  r=Taylor1([1.1,2.2,3.3])
 # showFastTaylor(a[1],r,cache[1]) =#