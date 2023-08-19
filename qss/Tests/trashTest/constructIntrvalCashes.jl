using BenchmarkTools
using InteractiveUtils
using TimerOutputs
#= 
vec1=(6.5,2.3)
vec2=(2.23,6.3)
@btime vec=(vec1[1]::Float64,vec2[1]::Float64)::NTuple{2, Float64} =#
#@btime vec=(vec1[1]::Float64,vec1[2]::Float64,vec2[1]::Float64,vec2[2]::Float64)
#@btime vec=(vec1[1]::Float64,vec2[1]::Float64)
 #@show vec
#= function mergesmallTuples(vec1::NTuple{M, Float64},vec2::NTuple{N, Float64}) where {M,N}
    if length(vec1)==0
        vec=vec2
    elseif length(vec2)==0
        vec=vec1
    elseif length(vec1)==1
        if length(vec2)==1
            vec=(vec1[1]::Float64,vec2[1]::Float64)
        else #1 2
            vec=(vec1[1]::Float64,vec2[1]::Float64,vec2[2]::Float64)
        end
    else #2 x
        if length(vec2)==1
            vec=(vec1[1]::Float64,vec1[2]::Float64,vec2[1]::Float64)
        else #2 2
            vec=(vec1[1]::Float64,vec1[2]::Float64,vec2[1]::Float64,vec2[2]::Float64)
        end
    end
    vec
end
@btime mergesmallTuples(vec1,vec2) =#





#display(vec)
#= res=filter((x) -> x < 3.0 , vec)
display(typeof(res)) =#

constructIntrval0(tup::Tuple{},tup2::Tuple{})=((0,Inf),)

#merge has 1 element
function constructIntrval1(tup::Tuple{Float64},tup2::NTuple{0, Float64})
    res=((0,tup[1]),)
end
function constructIntrval1(tup::NTuple{0, Float64},tup2::NTuple{1, Float64})
    #constructIntrval(tup2,tup)
    res=((0,tup2[1]),)
end

#merge has 2 elements
function constructIntrval2(tup::Tuple{Float64},tup2::Tuple{Float64})
  h1=min(tup[1],tup2[1])
  h4=max(tup[1],tup2[1])
  res=((0.0,h1),(h4,Inf))
end

function constructIntrval2(tup::NTuple{2, Float64},tup2::NTuple{0, Float64})
    h1=min(tup[1],tup[2])
    h4=max(tup[1],tup[2])
    res=((0.0,h1),(h4,Inf))
end 

function constructIntrval2(tup2::NTuple{0, Float64},tup::NTuple{2, Float64})
    h1=min(tup[1],tup[2])
    h4=max(tup[1],tup[2])
    res=((0.0,h1),(h4,Inf))
end 
#constructIntrval(tup2::Tuple{},tup::NTuple{2, Float64})=constructIntrval(tup,tup2)

#merge has 3 elements
function constructIntrval3(tup::NTuple{2, Float64},tup2::NTuple{1,Float64})
  t1=tup[1];t2=tup[2];t3=tup2[1]
  #t12=min(t1,t2);t21=max(t1,t2)
  if t1<t2
    t12=t1;t21=t2
  else
    t12=t2;t21=t1
  end

  if t3<t12
   res=((0.0,t3),(t12,t21))
  else
    if t3<t21
      res=((0.0,t12),(t3,t21))
    else
      res=((0.0,t12),(t21,t3))
    end
  end
  res
end
function constructIntrval3(tup2::Tuple{Float64},tup::NTuple{2, Float64})
   # constructIntrval(tup,tup2)
   t1=tup[1];t2=tup[2];t3=tup2[1]
  #t12=min(t1,t2);t21=max(t1,t2)
  if t1<t2
    t12=t1;t21=t2
  else
    t12=t2;t21=t1
  end

  if t3<t12
   res=((0.0,t3),(t12,t21))
  else
    if t3<t21
      res=((0.0,t12),(t3,t21))
    else
      res=((0.0,t12),(t21,t3))
    end
  end
  res
end

#merge has 4 elements
function constructIntrval4(tup::NTuple{2, Float64},tup2::NTuple{2, Float64})
  t1=tup[1];t2=tup[2];t3=tup2[1];t4=tup2[2]
    h1=min(t1,t2,t3,t4);h4=max(t1,t2,t3,t4)
    if t1==h1
      if t2==h4
        res=((0,t1),(min(t3,t4),max(t3,t4)),(t2,Inf))
      elseif t3==h4
        res=((0,t1),(min(t2,t4),max(t2,t4)),(t3,Inf))
      else
        res=((0,t1),(min(t2,t3),max(t2,t3)),(t4,Inf))
      end

    elseif t2==h1
      if t1==h4
        res=((0,t2),(min(t3,t4),max(t3,t4)),(t1,Inf))
      elseif t3==h4
        res=((0,t2),(min(t1,t4),max(t1,t4)),(t3,Inf))
      else
        res=((0,t2),(min(t1,t3),max(t1,t3)),(t4,Inf))
      end
     
    elseif t3==h1
      if t1==h4
        res=((0,t3),(min(t2,t4),max(t2,t4)),(t1,Inf))
      elseif t2==h4
        res=((0,t3),(min(t1,t4),max(t1,t4)),(t2,Inf))
      else
        res=((0,t3),(min(t1,t2),max(t1,t2)),(t4,Inf))
      end
      
    else
      if t1==h4
        res=((0,t4),(min(t2,t3),max(t2,t3)),(t1,Inf))
      elseif t2==h4
        res=((0,t4),(min(t1,t3),max(t1,t3)),(t2,Inf))
      else
        res=((0,t4),(min(t1,t2),max(t1,t2)),(t3,Inf))
      end
    end
   
    res
 
end




vec1=(23.3,)
vec2=()

#display(@code_warntype constructIntrval(vec1,vec2))
#= @show constructIntrval0(vec1,vec2)
@btime constructIntrval0(vec1,vec2) =#
#@show ress

#= oi=NTuple{3, Float64}((2.0,3.2,2.2))
display(typeof(oi)) =#

#println(methods(constructIntrval))

function constructIntrval(vec1,vec2)
    #= vec1=(1.2,3.3)
    vec2=() =#
    lv1=length(vec1);lv2=length(vec2)
    lv=length(vec1)+length(vec2)
    ress=nothing
    if lv==0
        ress=constructIntrval0(vec1,vec2)
    elseif lv==1
      ress=constructIntrval1(vec1,vec2)
    elseif lv==2
      #= @timeit "2vers" =#  ress=constructIntrval2(vec1,vec2)
    elseif lv==3
        ress=constructIntrval3(vec1,vec2)
    else
        ress=constructIntrval4(vec1,vec2)
    end
    ress
end


#= @show test()
@btime test() =#

function testouter()
  reset_timer!()
vec1=(23.3,)
vec2=(3.3,)

#@show test(cache,vec1,vec2)
#@show cache
@timeit "unsafeload" constructIntrval(vec1,vec2)
#= @show test()
@btime test() =#
print_timer()
end


testouter()