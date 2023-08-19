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

#= constructIntrval0(tup::Tuple{},tup2::Tuple{})=((0.0,Inf),)

#merge has 1 element
function constructIntrval10(tup::Tuple{Float64},tup2::NTuple{0, Float64})
    res=((0.0,tup[1]),)
end =#
#= function constructIntrval01(tup::NTuple{0, Float64},tup2::NTuple{1, Float64})
    #constructIntrval(tup2,tup)
    res=((0.0,tup2[1]),)
end =#

#merge has 2 elements
#= function constructIntrval11(tup::Tuple{Float64},tup2::Tuple{Float64})
  h1=min(tup[1],tup2[1])
  h4=max(tup[1],tup2[1])
  res=((0.0,h1),(h4,Inf))
end =#

function constructIntrval2(cache::Vector{Vector{Float64}},t1::Float64,t2::Float64)
    h1=min(t1,t2)
    h4=max(t1,t2)
   # res=((0.0,h1),(h4,Inf))
    cache[1][1]=0.0;cache[1][2]=h1;
    cache[2][1]=h4;cache[2][2]=Inf;
    return nothing
end 

#= function constructIntrval02(tup2::NTuple{0, Float64},tup::NTuple{2, Float64})
    h1=min(tup[1],tup[2])
    h4=max(tup[1],tup[2])
    res=((0.0,h1),(h4,Inf))
end  =#
#constructIntrval(tup2::Tuple{},tup::NTuple{2, Float64})=constructIntrval(tup,tup2)

#merge has 3 elements
function constructIntrval3(cache::Vector{Vector{Float64}},t1::Float64,t2::Float64,t3::Float64)
  #t1=tup[1];t2=tup[2];t3=tup2[1]
  t12=min(t1,t2);t21=max(t1,t2)
 #=  if t1<t2
    t12=t1;t21=t2
  else
    t12=t2;t21=t1
  end =#

  if t3<t12
  # res=((0.0,t3),(t12,t21))
   cache[1][1]=0.0;cache[1][2]=t3;
   cache[2][1]=t12;cache[2][2]=t21;
  else
    if t3<t21
     # res=((0.0,t12),(t3,t21))
      cache[1][1]=0.0;cache[1][2]=t12;
      cache[2][1]=t3;cache[2][2]=t21;
    else
      #res=((0.0,t12),(t21,t3))
      cache[1][1]=0.0;cache[1][2]=t12;
      cache[2][1]=t21;cache[2][2]=t3;
    end
  end
  return nothing
end
#= function constructIntrval12(tup2::Tuple{Float64},tup::NTuple{2, Float64})
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
end =#

#merge has 4 elements
function constructIntrval4(cache::Vector{Vector{Float64}},t1::Float64,t2::Float64,t3::Float64,t4::Float64,)
 # t1=tup[1];t2=tup[2];t3=tup2[1];t4=tup2[2]
    h1=min(t1,t2,t3,t4);h4=max(t1,t2,t3,t4)
    if t1==h1
      if t2==h4
       # res=((0,t1),(min(t3,t4),max(t3,t4)),(t2,Inf))
        cache[1][1]=0.0;cache[1][2]=t1;
        cache[2][1]=min(t3,t4);cache[2][2]=max(t3,t4);
        cache[3][1]=t2;cache[3][2]=Inf;
      elseif t3==h4
        #res=((0,t1),(min(t2,t4),max(t2,t4)),(t3,Inf))
        cache[1][1]=0.0;cache[1][2]=t1;
        cache[2][1]=min(t2,t4);cache[2][2]=max(t2,t4);
        cache[3][1]=t3;cache[3][2]=Inf;
      else
       # res=((0,t1),(min(t2,t3),max(t2,t3)),(t4,Inf))
        cache[1][1]=0.0;cache[1][2]=t1;
        cache[2][1]=min(t3,t2);cache[2][2]=max(t3,t2);
        cache[3][1]=t4;cache[3][2]=Inf;
      end

    elseif t2==h1
      if t1==h4
        #res=((0,t2),(min(t3,t4),max(t3,t4)),(t1,Inf))
        cache[1][1]=0.0;cache[1][2]=t2;
        cache[2][1]=min(t3,t4);cache[2][2]=max(t3,t4);
        cache[3][1]=t1;cache[3][2]=Inf;
      elseif t3==h4
        #res=((0,t2),(min(t1,t4),max(t1,t4)),(t3,Inf))
        cache[1][1]=0.0;cache[1][2]=t2;
        cache[2][1]=min(t1,t4);cache[2][2]=max(t1,t4);
        cache[3][1]=t3;cache[3][2]=Inf;
      else
       # res=((0,t2),(min(t1,t3),max(t1,t3)),(t4,Inf))
        cache[1][1]=0.0;cache[1][2]=t2;
        cache[2][1]=min(t3,t1);cache[2][2]=max(t3,t1);
        cache[3][1]=t4;cache[3][2]=Inf;
      end
     
    elseif t3==h1
      if t1==h4
        #res=((0,t3),(min(t2,t4),max(t2,t4)),(t1,Inf))
        cache[1][1]=0.0;cache[1][2]=t3;
        cache[2][1]=min(t2,t4);cache[2][2]=max(t2,t4);
        cache[3][1]=t1;cache[3][2]=Inf;
      elseif t2==h4
        #res=((0,t3),(min(t1,t4),max(t1,t4)),(t2,Inf))
        cache[1][1]=0.0;cache[1][2]=t3;
        cache[2][1]=min(t1,t4);cache[2][2]=max(t1,t4);
        cache[3][1]=t2;cache[3][2]=Inf;
      else
        #res=((0,t3),(min(t1,t2),max(t1,t2)),(t4,Inf))
        cache[1][1]=0.0;cache[1][2]=t3;
        cache[2][1]=min(t1,t2);cache[2][2]=max(t1,t2);
        cache[3][1]=t4;cache[3][2]=Inf;
      end
      
    else
      if t1==h4
        #res=((0,t4),(min(t2,t3),max(t2,t3)),(t1,Inf))
        cache[1][1]=0.0;cache[1][2]=t4;
        cache[2][1]=min(t3,t2);cache[2][2]=max(t3,t2);
        cache[3][1]=t1;cache[3][2]=Inf;
      elseif t2==h4
       # res=((0,t4),(min(t1,t3),max(t1,t3)),(t2,Inf))
        cache[1][1]=0.0;cache[1][2]=t4;
        cache[2][1]=min(t3,t1);cache[2][2]=max(t3,t1);
        cache[3][1]=t2;cache[3][2]=Inf;
      else
        #res=((0,t4),(min(t1,t2),max(t1,t2)),(t3,Inf))
        cache[1][1]=0.0;cache[1][2]=t4;
        cache[2][1]=min(t1,t2);cache[2][2]=max(t1,t2);
        cache[3][1]=t3;cache[3][2]=Inf;
      end
    end
   
    return nothing
 
end




#= ress=constructIntrval0((),())

@btime ress=constructIntrval1(vec1,vec2)

    ress=constructIntrval2((3.3,),(2.3,))

    ress=constructIntrval3((3.2,26.3),(1.2,))

    ress=constructIntrval4((3.3,6.2),(0.2,0.6)) =#
#display(@code_warntype constructIntrval(vec1,vec2))
#= @show constructIntrval0(vec1,vec2)
@btime constructIntrval0(vec1,vec2) =#
#@show ress

#= oi=NTuple{3, Float64}((2.0,3.2,2.2))
display(typeof(oi)) =#

#println(methods(constructIntrval))

function constructIntrval(cache::Vector{Vector{Float64}},vec1::NTuple{M, Float64},vec2::NTuple{N, Float64})where{M,N}
    #= vec1=(1.2,3.3)
    vec2=() =#
    lv1=length(vec1);lv2=length(vec2)
    lv=length(vec1)+length(vec2)
   
    if lv==0
        #res=((0.0,Inf),)
        cache[1][1]=0.0;cache[1][2]=Inf;
    elseif lv==1
        if lv1==1
            #res=((0.0,vec1[1]),)
            cache[1][1]=0.0;cache[1][2]=vec1[1];
        else
           # res=((0.0,vec2[2]),)
            cache[1][1]=0.0;cache[1][2]=vec2[1];
        end
    elseif lv==2
        if lv1==1
             constructIntrval2(cache,vec1[1],vec2[1])
        elseif lv1==2
           constructIntrval2(cache,vec1[1],vec1[2])
        else
           constructIntrval2(cache,vec2[1],vec2[2])
        end
    elseif lv==3
        if lv1==1
           constructIntrval3(cache,vec1[1],vec2[1],vec2[2])
         else
             constructIntrval3(cache,vec1[1],vec1[2],vec2[1])
        
         end
    else
       constructIntrval4(cache,vec1[1],vec1[2],vec2[1],vec2[2])
    end
    return nothing
end


function testouter()
    reset_timer!()
vec1=(23.3,)
vec2=(3.3,)
cache=Vector{Vector{Float64}}(undef,4)
for i =1:4
  cache[i]=zeros(2)
end
#@show test(cache,vec1,vec2)
#@show cache
@timeit "unsafeload" constructIntrval(cache,vec1,vec2)
#= @show test()
@btime test() =#
print_timer()
end


testouter()


#= if res1<0 && res2<0 && res3<0 && res4<0
  #res=((0.0,Inf),)
  cache[1][1]=0.0;cache[1][2]=Inf;
elseif res1>0 && res2<0 && res3<0 && res4<0
      cache[1][1]=0.0;cache[1][2]=res1;
elseif res2>0 && res1<0 && res3<0 && res4<0
      cache[1][1]=0.0;cache[1][2]=res2;
elseif res3>0 && res2<0 && res1<0 && res4<0
      cache[1][1]=0.0;cache[1][2]=res3;
elseif res4>0 && res2<0 && res3<0 && res1<0
      cache[1][1]=0.0;cache[1][2]=res4;
elseif res1>0 && res2>0 && res3<0 && res4<0
    constructIntrval2(cache,res1,res2)
elseif res1>0 && res3>0 && res2<0 && res4<0
    constructIntrval2(cache,res1,res3)
elseif res1>0 && res3>0 && res2<0 && res4<0
    constructIntrval2(cache,res1,res3) =#