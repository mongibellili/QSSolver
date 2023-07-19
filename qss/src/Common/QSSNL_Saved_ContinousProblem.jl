
struct savedNLODEContProblem{PRTYPE,T,Z,Y}<: NLODEProblem{PRTYPE,T,Z,Y} 
    prtype::Val{PRTYPE}
    a::Val{T}
    b::Val{Z}
    c::Val{Y}
    cacheSize::Int   
    initConditions::Function  
    eqs::Function
    jac::Function#Jacobian dependency..I have a der and I want to know which vars affect it...opposite of SD
    #jac::Vector{Vector{Int}}
    SD::Function#  I have a var and I want the der that are affected by it
    map::Function
    jacDim::Function
end
function getInitCond(prob::savedNLODEContProblem,i::Int)
    return prob.initConditions(i)
end
function saveNLodeProblemFunc(odeExprs::Expr,::Val{T},::Val{0},::Val{0},initCond::Dict{Union{Int,Expr},Float64},du::Symbol)where {T}
    if verbose println("SAVEnlodeprobfun only T= $T") end
    equs=Dict{Union{Int,Expr},Expr}()
    jac = Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}}()# set used because do not want to insert an existing varNum
    SD = Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}}()  # SD will be used in cretsdfun() to create the real SD in a function: this SD is not full SD: only extract partial SD from numbers
                                                              # symbol and expression will be stored like in Jac and later in cretsdfun() more work to get full SD 
    num_cache_equs=1#cachesize
    for argI in odeExprs.args
        #only diff eqs: du[]= number/one ref/call  
        if argI isa Expr &&  argI.head == :(=)  && argI.args[1] isa Expr && argI.args[1].head == :ref && argI.args[1].args[1]==du#&& ((argI.args[2] isa Expr && (argI.args[2].head ==:ref || argI.args[2].head ==:call ))||argI.args[2] isa Number)
            y=argI.args[1];rhs=argI.args[2]
            varNum=y.args[2] # order of variable
            if rhs isa Number # rhs of equ =number  
                equs[varNum]=:($((transformFSimplecase(:($(rhs))))))
            elseif rhs.head==:ref #rhs is only one var
                if rhs.args[1]==:q # check not needed
                extractJacDepNormal(varNum,rhs,jac ,SD ) 
                end
                equs[varNum ]=:($((transformFSimplecase(:($(rhs))))))
            else #rhs head==call...to be tested later for  math functions and other possible scenarios or user erros                 
                extractJacDepNormal(varNum,rhs,jac ,SD ) 
                temp=(transformF(:($(rhs),1))).args[2]  #number of caches distibuted   ...no need interpolation and wrap in expr....before was cuz quote....
                if num_cache_equs<temp 
                        num_cache_equs=temp
                end 
                equs[varNum]=rhs
            end 
        elseif @capture(argI, for counter_ in b_:niter_ loopbody__ end)
            specRHS=loopbody[1].args[2]
             extractJacDepLoop(b,niter,specRHS,jac ,SD )   
             temp=(transformF(:($(specRHS),1))).args[2]
                if num_cache_equs<temp 
                    num_cache_equs=temp
                end 
               equs[:(($b,$niter))]=specRHS            
        else#end of equations and user enter something weird...handle later
               # error("expression $x: top level contains only expressions 'A=B' or 'if a b' ")#wait until exclude events
       # end#end cases inside @capture
        end #end for x in 
      
    end #end for args #########################################################################################################

    fname= :f #default name
    path="./temp.jl" #default path
    if odeExprs.args[1] isa Expr && odeExprs.args[1].args[2] isa Expr && odeExprs.args[1].args[2].head == :tuple#if user chose path and name
        fname= odeExprs.args[1].args[2].args[1]
        path=odeExprs.args[1].args[2].args[2]
    end
    initCondfunction=createInitCondFun(initCond,fname)
    diffEqfunction=createContEqFun(equs,fname) # found in not saved cont prob file
    jacFunction=createJacFun(jac,fname)
    jacDimFunction=createJacDimensionFun(jac,fname) # found in not saved cont prob file
    SDFunction=createSDFun(SD,fname)
    mapFun=createMapFun(jac,fname) # found in not saved cont prob file
   

   s=" prob=savedNLODEContProblem(Val(2),Val($T),Val(0),Val(0),$num_cache_equs,$(Symbol(:IC,fname)),$fname,$(Symbol(:jac,fname)),$(Symbol(:SD,fname)),$(Symbol(:map,fname)),$(Symbol(:jacDimension,fname)))"
    myex1=Meta.parse(s)
    Base.remove_linenums!(myex1)
  
    def=Dict{Symbol,Any}()
    def[:head] = :function
    def[:name] = Symbol(:prob_,fname)  
    #def[:args] = []
    def[:body] = myex1
    probFun=combinedef(def)


    open(path, "a") do io    
        println(io,string(initCondfunction))  
        println(io,string(diffEqfunction)) 
        print(io,"@resumable ") 
        println(io,string(jacFunction)) 
        print(io,"@resumable ")     
        println(io,string(SDFunction))  
        println(io,string(mapFun)) 
        println(io,string(jacDimFunction)) 

        println(io,string(probFun)) 
      end

end


function createInitCondFun(initCond::Dict{Union{Int,Expr},Float64},funName::Symbol)
    s="if j==0 return 0.0\n"
    for i in initCond
        Base.remove_linenums!(i[1])
        Base.remove_linenums!(i[2])
        if i[1] isa Int
            s*="elseif j==$(i[1])   ;return $(i[2])\n"
        end
        if i[1] isa Expr
            s*="elseif $(i[1].args[1])<=j<=$(i[1].args[2])   ;return $(i[2])\n"
        end
    end
    s*=" end "
    myex1=Meta.parse(s)
    Base.remove_linenums!(myex1)
    def=Dict{Symbol,Any}()
    def[:head] = :function
    def[:name] = Symbol(:IC,funName)   
    def[:args] = [:(j::Int)]
    def[:body] = myex1  
    functioncode=combinedef(def)
   # @show functioncode;functioncode
end


function createJacFun(jac:: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}},funName::Symbol)
    ss="if i==0 return nothing\n"
    for dictElement in jac
    #=  Base.remove_linenums!(dictElement[1])
        Base.remove_linenums!(dictElement[2]) =#
       # counterJac=1
        if dictElement[1] isa Int
            ss*="elseif i==$(dictElement[1])  \n"
            for k in dictElement[2]
                ss*="@yield $k \n"
               # counterJac+=1
            end
           # ss*=" return nothing \n"
        elseif dictElement[1] isa Expr
            ss*="elseif $(dictElement[1].args[1])<=i<=$(dictElement[1].args[2])  \n"
            for k in dictElement[2]
                ss*="@yield $k \n"
               # counterJac+=1
            end
           # ss*=" return nothing \n"
        end     
    end
    ss*=" end \n"         
    myex1=Meta.parse(ss)
    Base.remove_linenums!(myex1)
    def1=Dict{Symbol,Any}() #any changeto Union{expr,Symbol}  ????
    def1[:head] = :function
    def1[:name] = Symbol(:jac,funName)  
    def1[:args] = [:(i::Int)]
    def1[:body] = myex1
    functioncode1=combinedef(def1)
end

function createSDFun(SD:: Dict{Union{Int,Expr},Set{Union{Int,Symbol,Expr}}},funName::Symbol)
    allEpxpr=Expr(:block)
    s="if i==0 return nothing\n"
    for dictElement in SD
        #it is needed when counting the influences to b:niter later..
        if dictElement[1] isa Int # key is an int
            s*="elseif i==$(dictElement[1])  \n"
            for k in dictElement[2]   #elments values 
                if k isa Int # a value is an int
                    s*="@yield $k \n"
                elseif k isa Expr # a value is an expression
                    s*="for j=$(k.args[1]):$(k.args[2])  @yield j end \n"           
                end
            end
        elseif dictElement[1] isa Expr # key is an expression
            for k in dictElement[2]
                sspecial="if $(dictElement[1].args[1])<=$k<=$(dictElement[1].args[2])  @yield $k end  \n"
                ex=Meta.parse(sspecial)
                Base.remove_linenums!(ex.args[2])
                push!(allEpxpr.args,ex)
            end
        end     
    end
    s*=" end \n"         
    myex=Meta.parse(s)
    Base.remove_linenums!(myex)
    push!(allEpxpr.args,myex)
    sEnd="return nothing"
    myex3=Meta.parse(sEnd)
    push!(allEpxpr.args,myex3)
    Base.remove_linenums!(allEpxpr)
    def=Dict{Symbol,Any}() # can be the same dict passed around
    def[:head] = :function
    def[:name] = Symbol(:SD,funName)  
    def[:args] = [:(i::Int)]
    def[:body] = allEpxpr 
    functioncode=combinedef(def)
end

