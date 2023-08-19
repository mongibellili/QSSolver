using qss
using XLSX
using BenchmarkTools
include("/home/mongi/projects/formalqss/Tests/ProblemC/typeCSelect.jl")
function solveProblem(prblem::Function,ft::Float64,solver::Tuple{Val{solType}, Val{V}},absTol,relTol)where {V,solType} 
    pr=prblem()
    odeprob=pr[1]
    x1=pr[2]
    x2=pr[3]
    timenmliqss=0.0
    #= absTol=1e-5
    relTol=1e-3 =#
    solnmliqss2=QSS_Solve(odeprob,solver,dQmin=absTol,dQrel=relTol,finalTime=ft,maxErr=100*relTol)
    #save_Sol(solnmliqss2)
    solnmliqss2Interp=solInterpolated(solnmliqss2,0.01)
    er1=getError(solnmliqss2Interp,1,x1)  
    er2=getError(solnmliqss2Interp,2,x2) 
    # timenmliqss=@belapsed QSS_Solve($odeprob,nmliqss2(),dQmin=$absTol,saveat=0.01,dQrel=$relTol,finalTime=100.0)
   #=  res= ("$(solnmliqss2.sysName)",relTol,absTol,(er1+er2)/2,solnmliqss2.totalSteps,solnmliqss2.simulStepCount,timenmliqss)
    @show res =#
    resnmliqss= ((er1+er2)/2,solnmliqss2.totalSteps,solnmliqss2.simulStepCount)
    #@show resnmliqss
    errr=(er1+er2)/2
    if errr>relTol*50
      println("name=$(solnmliqss2.sysName); relTol=$relTol absTol=$absTol error=$errr")
    end
    return resnmliqss
end


function mainTest(cfuns::Vector{Function},title::String)
   #=  absTols=[#= 1e-2, =#1e-3,1e-4,1e-5,1e-6 #=,1e-7,1e-8 =#]
    relTols=[#= 1e-1, =#1e-2,1e-3,1e-4,1e-5#= ,1e-6 =#] =#
    absTol=1e-3
    relTol=1e-2
  #=   tolColumn=[]
    for relTol in relTols
      for absTol in absTols
        if absTol<=relTol
        push!(tolColumn,(relTol,absTol))
        end
      end
    end =#
    ft=100.0
    solver=nmliqss1()
   
    numProblems=length(cfuns)
    numProblems%10!=0 && error("provide a set of multiple of 10")
    numDistinctJac=div(numProblems,10)
    cResults=Vector{Tuple{Float64,Float64,Float64}}(undef, numProblems)
    vecPerInput=Vector{Vector{Tuple{Float64,Float64,Float64}}}(undef, 10)
    for i=1:10
        vecPerInput[i]=[]
    end
    vecPerJac=Vector{Vector{Tuple{Float64,Float64,Float64}}}(undef, numDistinctJac)
    
    for i=1:numDistinctJac
        vecPerJac[i]=[]
    end
    for i=1:numProblems
      cResults[i]=solveProblem(cfuns[i],ft,solver,absTol,relTol) 
      if  i%10==0
        push!(vecPerInput[10],cResults[i])
      end
      for k=1:10
          if  i%10==k
              push!(vecPerInput[k],cResults[i])
          end
        
          if div(i,10)==k-1 && i%10!=0 || (div(i,10)==k && i%10==0) 
              push!(vecPerJac[k],cResults[i])
          end
          if div(i,10)==k+9 && i%10!=0 || (div(i,10)==10+k && i%10==0)
              push!(vecPerJac[10+k],cResults[i])
          end
          if div(i,10)==k+19 && i%10!=0 || (div(i,10)==20+k && i%10==0)
              push!(vecPerJac[20+k],cResults[i])
          end
      end
    end
  #@show cResults
   #Problemletter=["B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z","AA","AB","AC","AD","AE","AF","AG","AH","AI","AJ","AK","AL","AM","AN","AO","AP","AQ","AR","AS","AT","AU","AV","AW","AX","AY","AZ","BA","BB","BC","BD","BE","BF","BG","BH","BI","BJ","BK","BL","BM","BN","BO","BP","BQ","BR","BS","BT","BU","BV","BW","BX","BY","BZ"]#for excel not be confused with problem name

    #= XLSX.openxlsx("AllMetrics_LTI_$(title)_nmliqss2_$(relTol).xlsx", mode="w") do xf
    #  XLSX.openxlsx("LTI_testLette_r.xlsx", mode="w") do xf
      sheet = xf[1]
      sheet["A1"]="LTI_$(title)_nmliqss2_"
      sheet["A2"]="sub_system"
      sheet["B2"]="Error"
      sheet["C2"]="Steps"
      sheet["D2"]="simulSteps"
      for k=3:numProblems+2
          sheet["A$(k)"] ="$(cfuns[k-2])"
          sheet["B$(k)"] =cResults[k-2][1]
          sheet["C$(k)"] =cResults[k-2][2]
          sheet["D$(k)"] =cResults[k-2][3]
      end
    end =#
    Problemletter=["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z","AA","AB","AC","AD","AE","AF","AG","AH","AI","AJ","AK","AL","AM","AN","AO","AP","AQ","AR","AS","AT","AU","AV","AW","AX","AY","AZ","BA","BB","BC","BD","BE","BF","BG","BH","BI","BJ","BK","BL","BM","BN","BO","BP","BQ","BR","BS","BT","BU","BV","BW","BX","BY","BZ"]#for excel not be confused with problem name
    XLSX.openxlsx("AllMetrics__perInput_LTI_C_$(title)__$(relTol)_abs-3.xlsx", mode="w") do xf
      #  XLSX.openxlsx("LTI_testLette_r.xlsx", mode="w") do xf
        sheet = xf[1]
        sheet["A1"]="LTI_$(title)__e-2e-3_"
        for k=1:10
          sheet["$(Problemletter[4k-3])2"] =collect(("problem","error","totalSteps","simul_steps"))
      
            for i=1:numDistinctJac
            sheet["$(Problemletter[4k-3])$(i+2)"] =collect(("C$(k+(i-1)*10)",vecPerInput[k][i][1],vecPerInput[k][i][2],vecPerInput[k][i][3]))
            end
        end
      end
end
 cfuns=[C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15,C16,C17,C18,C19,C20,C21,C22,C23,C24,C25,C26,C27,C28,C29,C30,C31,C32,C33,C34,C35,C36,C37,C38,C39,C40,C41,C42,C43,C44,C45,C46,C47,C48,C49,C50,C51,C52,C53,C54,C55,C56,C57,C58,C59,C60,C61,C62,C63,C64,C65,C66,C67,C68,C69,C70,C71,C72,C73,C74,C75,C76,C77,C78,C79,C80,C81,C82,C83,C84,C85,C86,C87,C88,C89,C90,C91,C92,C93,C94,C95,C96,C97,C98,C99,C100,C101,C102,C103,C104,C105,C106,C107,C108,C109,C110,C111,C112,C113,C114,C115,C116,C117,C118,C119,C120,C121,C122,C123,C124,C125,C126,C127,C128,C129,C130,C131,C132,C133,C134,C135,C136,C137,C138,C139,C140,C141,C142,C143,C144,C145,C146,C147,C148,C149,C150,C151,C152,C153,C154,C155,C156,C157,C158,C159,C160,C161,C162,C163,C164,C165,C166,C167,C168,C169,C170,C171,C172,C173,C174,C175,C176,C177,C178,C179,C180,C181,C182,C183,C184,C185,C186,C187,C188,C189,C190,C191,C192,C193,C194,C195,C196,C197,C198,C199,C200,C201,C202,C203,C204,C205,C206,C207,C208,C209,C210,C211,C212,C213,C214,C215,C216,C217,C218,C219,C220]

 #cfuns=[C81,C82,C83,C84,C85,C86,C87,C88,C89,C90,C91,C92,C93,C94,C95,C96,C97,C98,C99,C100,C121,C122,C123,C124,C125,C126,C127,C128,C129,C130,C131,C132,C133,C134,C135,C136,C137,C138,C139,C140,C151,C152,C153,C154,C155,C156,C157,C158,C159,C160,C201,C202,C203,C204,C205,C206,C207,C208,C209,C210,C211,C212,C213,C214,C215,C216,C217,C218,C219,C220]
#cfuns=[C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15,C16,C17,C18,C19,C20,C21,C22,C23,C24,C25,C26,C27,C28,C29,C30]
mainTest( cfuns," order1_exacteA_analyticUpdateQ_only1Delta_qx_Inf_withoutSimulupdate") 



 #=cfuns1_20=[C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15,C16,C17,C18,C19,C20]
cfuns21_40=[C21,C22,C23,C24,C25,C26,C27,C28,C29,C30,C31,C32,C33,C34,C35,C36,C37,C38,C39,C40]
cfuns41_60=[C41,C42,C43,C44,C45,C46,C47,C48,C49,C50,C51,C52,C53,C54,C55,C56,C57,C58,C59,C60]
cfuns61_80=[C61,C62,C63,C64,C65,C66,C67,C68,C69,C70,C71,C72,C73,C74,C75,C76,C77,C78,C79,C80]
cfuns81_100=[C81,C82,C83,C84,C85,C86,C87,C88,C89,C90,C91,C92,C93,C94,C95,C96,C97,C98,C99,C100]
cfuns101_120=[C101,C102,C103,C104,C105,C106,C107,C108,C109,C110,C111,C112,C113,C114,C115,C116,C117,C118,C119,C120]
cfuns121_140=[C121,C122,C123,C124,C125,C126,C127,C128,C129,C130,C131,C132,C133,C134,C135,C136,C137,C138,C139,C140]
cfuns141_160=[C141,C142,C143,C144,C145,C146,C147,C148,C149,C150,C151,C152,C153,C154,C155,C156,C157,C158,C159,C160]
cfuns161_180=[C161,C162,C163,C164,C165,C166,C167,C168,C169,C170,C171,C172,C173,C174,C175,C176,C177,C178,C179,C180]
cfuns181_200=[C181,C182,C183,C184,C185,C186,C187,C188,C189,C190,C191,C192,C193,C194,C195,C196,C197,C198,C199,C200]
cfuns201_220=[C201,C202,C203,C204,C205,C206,C207,C208,C209,C210,C211,C212,C213,C214,C215,C216,C217,C218,C219,C220]
mainTest( cfuns1_20," cfuns1_20")
mainTest( cfuns21_40," cfuns21_40")
mainTest( cfuns41_60," cfuns41_60")
 mainTest( cfuns61_80," cfuns61_80")
mainTest( cfuns81_100," cfuns81_100")
mainTest(cfuns101_120," cfuns101_120")
mainTest( cfuns121_140," cfuns121_140")
mainTest(cfuns141_160," cfuns141_160")
mainTest( cfuns161_180," cfuns161_180")
mainTest( cfuns181_200," cfuns181_200")
mainTest( cfuns201_220," cfuns201_220") =#