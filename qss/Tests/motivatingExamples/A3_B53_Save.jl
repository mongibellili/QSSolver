

using formalqss
function test53()
       @saveNLodeProblem begin
        destination=(sysb53,"/home/unknown/formalqss/Tests/trash/testSave.jl")
        # u = [-1.0, -2.0]
u[1]=-1.0; u[2]=-2.0
        # discrete = [0.0]
         du[1] = -20.0*u[1]-80.0*u[2]+1600.0
         du[2] =1.24*u[1]-0.01*u[2]+0.2
     end      
  end
  test53()
  


 