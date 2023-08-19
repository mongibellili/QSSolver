function probB11()
    odeprob = @NLodeProblem begin
         #sys b11
        u = [-1.0, -2.0]
        
        du[1] = 0.15*u[1]+0.25*u[2]
        du[2] =-7.0*u[1]-4.0*u[2]    
    end
end


 function probB12()
    odeprob = @NLodeProblem begin
         #sys b11
        u = [-1.0, -2.0]
        
        du[1] = 0.15*u[1]+0.25*u[2]-0.25
        du[2] =-7.0*u[1]-4.0*u[2]+4.0    
    end
end


 function probB13()
    odeprob = @NLodeProblem begin
         #sys b13
        u = [-1.0, -2.0]
        
        du[1] = 0.15*u[1]+0.25*u[2]-5.0
        du[2] =-7.0*u[1]-4.0*u[2]+80.0
    end
  
end

 function probB14()
    odeprob = @NLodeProblem begin
         #sys b14
        u = [-1.0, -2.0]
        
        du[1] = 0.15*u[1]+0.25*u[2]-0.0025
        du[2] =-7.0*u[1]-4.0*u[2]+0.04
    end
   
end

function probB15()
    odeprob = @NLodeProblem begin         
        #sys b15
        u = [-1.0, -2.0]
        
        du[1] = 0.15*u[1]+0.25*u[2]+0.05
        du[2] =-7.0*u[1]-4.0*u[2]-10.0
       
    end
end

function probB21()
    odeprob = @NLodeProblem begin   
         #sys b21
        u = [-1.0, -2.0]
        
        du[1] = 0.1*u[1]+0.25*u[2]
        du[2] =-4.0*u[1]-7.0*u[2]
    end
end

function probB22()
    odeprob = @NLodeProblem begin    
         #sys b22
        u = [-1.0, -2.0]
        
        du[1] = 0.1*u[1]+0.25*u[2]-0.25
        du[2] =-4.0*u[1]-7.0*u[2]+7.0
    end   
end


function probB23()
    odeprob = @NLodeProblem begin    
         #sys b25
        u = [-1.0, -2.0]
        
        du[1] = 0.1*u[1]+0.25*u[2]-5.0
        du[2] =-4.0*u[1]-7.0*u[2]+140.0
    end
end

function probB24()
    odeprob = @NLodeProblem begin       
         #sys b24
        u = [-1.0, -2.0]
        
        du[1] = 0.1*u[1]+0.25*u[2]-0.0025
        du[2] =-4.0*u[1]-7.0*u[2]+0.07     
    end
end

function probB25()
    odeprob = @NLodeProblem begin        
         #sys b25
        u = [-1.0, -2.0]
        
        du[1] = 0.1*u[1]+0.25*u[2]-0.05
        du[2] =-4.0*u[1]-7.0*u[2]-1.0
    end
end

function probB31()
    odeprob = @NLodeProblem begin        
         #sys b31
        u = [-1.0, -2.0]
        
        du[1] = -0.25*u[1]+0.15*u[2]
        du[2] =-4.0*u[1]-7.0*u[2]
    end

end

function probB32()
    odeprob = @NLodeProblem begin
         
         #sys b31
        u = [-1.0, -2.0]
        
        du[1] = -0.25*u[1]+0.15*u[2]-0.15
        du[2] =-4.0*u[1]-7.0*u[2]+7.0
    end
end

function probB33()
    odeprob = @NLodeProblem begin        
         #sys b31
        u = [-1.0, -2.0]
        
        du[1] = -0.25*u[1]+0.15*u[2]-3.0
        du[2] =-4.0*u[1]-7.0*u[2]+140.0
    end
end


function probB34()
    odeprob = @NLodeProblem begin      
         #sys b31
        u = [-1.0, -2.0]
        
        du[1] = -0.25*u[1]+0.15*u[2]-0.0015
        du[2] =-4.0*u[1]-7.0*u[2]+0.07
    end
end

function probB35()
    odeprob = @NLodeProblem begin       
         #sys b35
        u = [-1.0, -2.0]
        
        du[1] = -0.25*u[1]+0.15*u[2]-0.65
        du[2] =-4.0*u[1]-7.0*u[2]-1.0
    end
end

function probB41()
    odeprob = @NLodeProblem begin
         #sys b41
        u = [-1.0, -2.0]
        
        du[1] = -0.01*u[1]+1.24*u[2]
        du[2] =-80.0*u[1]-20.0*u[2]    
    end  
end


 function probB42()
    odeprob = @NLodeProblem begin
         #sys b42
        u = [-1.0, -2.0]
        
        du[1] = -0.01*u[1]+1.24*u[2]-1.24
        du[2] =-80.0*u[1]-20.0*u[2]+20.0      
    end
   
end

function probB43()
    odeprob = @NLodeProblem begin
        u = [-1.0, -2.0]
        
        du[1] = -0.01*u[1]+1.24*u[2]-24.8
        du[2] =-80.0*u[1]-20.0*u[2]+400.0
    end 
  
end


function probB44()
    odeprob = @NLodeProblem begin
         #sys b44
        u = [-1.0, -2.0]
        
        du[1] = -0.01*u[1]+1.24*u[2]-0.0124
        du[2] =-80.0*u[1]-20.0*u[2]+0.2
    end
  
end

  function probB45()
    odeprob = @NLodeProblem begin         
        #sys b45
        u = [-1.0, -2.0]
        
        du[1] = -0.01*u[1]+1.24*u[2]-1.26
        du[2] =-80.0*u[1]-20.0*u[2]-140.0
    end
  
end

##########################################

 function probB51()
    odeprob = @NLodeProblem begin
         #sys b51
        u = [-1.0, -2.0]
        
        du[1] = -20.0*u[1]-80.0*u[2]
        du[2] =1.24*u[1]-0.01*u[2]      
    end  
end

function probB52()
    odeprob = @NLodeProblem begin
         #sys b52
        u = [-1.0, -2.0]
        
        du[1] = -20.0*u[1]-80.0*u[2]+80.0
        du[2] =1.24*u[1]-0.01*u[2]+0.01     
    end 
end

function probB53()
    odeprob = @NLodeProblem begin
         #sys b53
        u = [-1.0, -2.0]
        
        du[1] = -20.0*u[1]-80.0*u[2]+1600.0
        du[2] =1.24*u[1]-0.01*u[2]+0.2
    end  
end


function probB54()
    odeprob = @NLodeProblem begin
         #sys b54
        u = [-1.0, -2.0]
        
        du[1] = -20.0*u[1]-80.0*u[2]+0.8
        du[2] =1.24*u[1]-0.01*u[2]+0.0001
    end
end

function probB55()
    odeprob = @NLodeProblem begin         
        #sys b55
        u = [-1.0, -2.0]
        
        du[1] = -20.0*u[1]-80.0*u[2]+40.0
        du[2] =1.24*u[1]-0.01*u[2]+2.49
    end
end

