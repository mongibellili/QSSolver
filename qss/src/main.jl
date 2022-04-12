
using qss
initialTime=0.1
finalTime=4.0
dQmin=0.0001
dQrel=0.005
order=1 # order2 means we will consider second derivatives
solver="qss1"
initConditions=[1.0,2.0]
jacobian=[0.0 1.0; -1.0 -1.0]

settings = ModelSettings(initialTime,finalTime,dQmin,dQrel,order,solver,initConditions,jacobian)
println(settings.solver)

simulator = QSS_simulator(settings)
QSS_simulate(simulator)
#display(simulator.qssData.x)
plotX(simulator)

    
    

