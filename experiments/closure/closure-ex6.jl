
 function parentfunction()
    
        x = 0
        function getx()
            x
        end
       function counter()
         x += 1
       end
    ()->(counter,getx)
 end


i2 = parentfunction() # i2 is the counter function
#println( i2().getx) #LoadError: type Tuple has no field getx
#println( i2.getx)#getx
println( i2.getx())#0  #this is the correct way
