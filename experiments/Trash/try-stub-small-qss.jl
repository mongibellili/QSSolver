
using StaticArrays
using BenchmarkTools


#TP13

struct ModelSettings    
  initConditions ::SVector{2,Float64}
  jacobian::SMatrix{2,2,Float64}
  finalTime:: Float64   
  initialTime :: Float64  
  dQmin ::Float64    
  dQrel ::Float64
  order :: Int       
  solver ::Val{T} where {T}
end
function ModelSettings(initConditions ,jacobian ,finalTime; initialTime=0.0,dQmin=1e-6,dQrel=1e-3,order=1,solver=Val(1))
  return ModelSettings(initConditions, jacobian, finalTime,initialTime,dQmin,dQrel,order,solver)
end
qss1()=Val(1)
qss2()=Val(2)
qss3()=Val(3)
liqss()=Val(4)

mutable struct DataSvector
  x::MVector{2,Float64} 
  tx::MVector{2,Float64} 
  q::MVector{2,Float64} 
  tq::MVector{2,Float64} 
  quantum::MVector{2,Float64} 
end



function integrate(probSetting::ODEProblem)
  p1 = DataSvector(zeros(2),zeros(2),zeros(2),zeros(2),zeros(2))
  s=0.25* probSetting.jacobian[1,2]
  for k=1:1e+2
    for j = 1:2
      p1.quantum[j]+=0.01
      for i = 1:1e+1
        n=rand(1:2)
        p1.x[j] =p1.x[n]+i*0.5 +0.25* probSetting.jacobian[1,2]
        p1.q[j]=p1.x[j]+p1.quantum[j]
      end
      p1.tx[j]+=1.0
      p1.tq[j]=p1.tx[j]
    end
  end
  p1
end

#main : user should provide and code this
function qss_Svector()
  ft=5
  initConditions=@SVector[1.0,2.0]
  jacobian=@SMatrix[0.0 1.0;-1.0 -1.0 ]
  problem=ModelSettings(initConditions,jacobian,ft)
  #display(problem.jacobian[1,2])
  s=integrate(problem)
  #display(s)
end
#qss_Svector()
@btime qss_Svector()







#=
struct DataMvector
  x::MVector{5,Float64} 
  y::Int
end
struct DataNormal
  x::Vector{Float64}
  y::Int
end

function qss_Mvector()
  s=0.0
  #v1 =  @MVector[1.1,2.2,3.3,4.4,5.5]
  p1 = SimulatorDataMvector(v1, 5)
  for i = 1:1e+3
    for j = 1:5
      s=s+p1.x[j]    
    end
  end
  s
end
function qss_vector()
  s=0.0
  v2 = [1.1,2.2,3.3,4.4,5.5]
  p2 = SimulatorDataNormal(v2, 5)

  for i = 1:1e+3 for j = 1:5
        s=s+p2.x[j]
  end end
  s
end


=#




