using Plots;gr()
using StaticArrays
using qss
#=
in later versions, a parser will be added to let the user describe their problem in an easy way like this:
function initialvariables()
    x[1]=0 # state var1     use x[
    dx[1]=0 #derx1          use dx[     
    x[2]=0
    discreteContact=0       use discrete
    #for loop for many vars
end
function equations()
    m=1.0       use any naming except x[ or  dx[
    b=30.0
    g=9.8
    k=1e6  
    dx[1]=x[2]
    dx[2]=-g-(d/m)*(k*x[1]+contact*x[2])
end
function events()
    if x1+discreteContact >0
        discreteContact=expression...
    else
        x[1]=another expression...
    end
    if x2-x1+2 >0
        discreteContact=expression...
    else
        discreteContact=another expression...
    end
end
problem=Problem(initialVariables,equations,events)
settings=ProblemSettings(finalTime,initialTime,dQmin,dQrel,solver,savetimeincrement) # allow user to skip args and default values used
 sol=QSS_solve(problem,settings)
display(plot!(sol[1],sol[2][1])) # for the case of two vars
=#
# --------the code produced by the parser looks like:
function qssApproachInitialInside() 
    DEBUG=false
    #-------------------------------------state vars and SD matrix--------------------------
    #x[1](initTime)=10 ;x[2](initTime)=0 ; initTime if not explicitly given it is zero...parser produce this code
    initConditions=@SVector[10.0,0.0]
    m=1.0
    b=30.0
    g=9.8
    k=1e+5 
    #dx[1]=x[2]
    #dx[2]=-g-(d/m)*(k*x[1]+contact*x[2])
    # later a parser will find the coeff of the statevars x[] convert to the following code
    f1(q,d,tq,order)=q[order+1]   #later specify the arguments ...now they are any any any any
     #q1 is always =q[1], q2 in order 2 becomes q[3] formula is "order*index-order+1"
    f2(q,d,tq,order)= -g-(d[1]/m)*(k*q[1]+b*q[order+1])
   

#=     f1(q,d,t,order)=q[order+1]
    f2(q,d,t,order)=-q[1]-q[order+1] +1 =#
 
    #dx[2]=-g-(d/m)*(k*x[1]+contact*x[2])
    jacobian=SVector{2,Function}(f1,f2)  
    SD=@SMatrix[2 0;1 2] # ncolomns = number of der, number of rows = number of states
     #the SD matrix is build inside (at initialize integrator)---leter can be moved here just to be consistent with other dependence matrices

    #-------------------------------------input vars and multiplexer-------------------------
    # any term left from the equations (not a coeff) is an inputvar. if all inputvars do not depend on time the parser will create the following code
    inputVars=@SVector[0.0,-g] 
    # if one of the input var depends on time the parser will produce the following code
    #=u1(t)=t+1
    u2(t)=10-t*t+t +cos(t)
    inputVars=SVector{2,Function}(u1,u2) =# 
    #ID=@SVector[n m...] #add this feature later for now mapping goes j to j  #...[3 1] means inputvar1 influences derx3 and inputVar2 influences derx1...note that derx2 does not have an inpputvar
   
    #---------------------------------------discrete vars------------------------------
    # any variable starts with discrete in its name, the compiler creates a vector and adds its initial value, then constructs an svector out of the vector
    #discreteContact=0.0  -->
    discreteVars=@MVector[0.0]  
    
     #---------------------------------------zero crossing---------------------------
     #for any expression after the "if" in the events function the parser counts them and creates:
    function zc_1(q,d,t,order) 
       # println("q1= ",q[1])
        (sign(q[1]),q[1])# the same expression after "if" in eventhandlerFun
    end

   
    function zc_2(q,d,t,order) 
        #println(" in zc2")
        (sign(q[1]-2),q[1]-2)
    end 
     zcFunctions=SVector{2,Function}(zc_1,zc_2)  # this 1 (the number of zc) determines the number of columns of HZ and SZ
     SZ=@SMatrix[1 2;0 0]  #x1 influences the zc1 function...[2 0;1 2] means x1 influences zc2 and x2 influences both zc1 and zc2
    
 
    #---------------------------------------event handlers---------------------------

    function ev1(q,d,tq)
        if q[1] >0  #not tested  ...maybe expression=0
            #discreteVars=setindex(discreteVars,0.0,1)
            discreteVars[1]=0.0
        else
            discreteVars[1]=1.0
            #discreteVars=setindex(discreteVars,1.0,1)
        end
    end
    function ev2(q,d,tq)
#=         if q[1]-2 > 0  #not tested  ...maybe expression=0
            #discreteVars=setindex(discreteVars,0.0,1)
            discreteVars[1]=0.0
        else
            discreteVars[1]=0.0
            #discreteVars=setindex(discreteVars,1.0,1)
        end =#
    end
    eventHandlerFunctions=SVector{2,Function}(ev1,ev2)   # this 1 (the number of events) determines the number of rows of HZ and HD
    HD=@SMatrix[0 2] #ev1 influence only derx2...[0 2;1 2] means ev1 influences derx2 and ev2 influences both derx1 and derx2
    HZ=@SMatrix[0 0] #there is one event and it does not influence the only 1 zc...[1 0;1 2] means event1 influences only zc1 and event2 influences both zc1 and zc2    
    

   
   

    psettings = ProblemSettings(9.0,saveat(0.001),qss1())
    #problem=Problem(initialVariables,equations,events) #later with the parser, user enters this
    prob = QSS_Problem(initConditions,jacobian,inputVars,discreteVars,zcFunctions,eventHandlerFunctions,SD,SZ,HD,HZ) # this should be done under the hood after parsing
    sol=QSS_Solve(psettings,prob)
    display(plot!(sol[1],sol[2][1]))
    display(plot!(sol[1],sol[2][2])) 
    println("done") 
    readline()
      



    if DEBUG
        println("jacobain");display(jacobian);println()
        println("inputvars");display(inputVars);println()
        println("discrete vars");display(discreteVars);println()
        println("eventhandlerFunctions");display(eventHandlerFunctions);println()
        println("HZ");display(HZ);println()
        println("HD");display(HD);println()
        println("zcFunctions");display(zcFunctions);println()
        println("SZ");display(SZ);println()
    end   
end

qssApproachInitialInside()


# odeDiffEquPackage()





#ZH=@SVector[...] # does not exist cuz each zc func has an eventHandler and their indexing is the same ie  zc1 triggers ev1...zc(j) triggers ev(j)

    #jacobian=@SMatrix[0.0 1.0;a b]
#=     u1(t)=-g
    inputVars=SVector{2,Function}(u1) =#

   


       #inputJac=@SMatrix[0.0 0.0;0.0 1.0 ]
#=    inputVars=@SVector[0.0,1.0]
   settings = ModelSettings(initConditions,jacobian,5.7,saveat(0.1),qss1(),inputVars) =#

       #states=2
    #settings = ModelSettings(initConditions,jacobian,5.7,saveat(0.1),qss2())#do not call saveat to not save, i should fix when called with zero; it does not save at all
   # settings = ModelSettings(initConditions,jacobian,5.7)