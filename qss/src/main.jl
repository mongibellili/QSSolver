
using qss
initialTime=0.1
finalTime=5.0
dQmin=0.01
dQrel=0.5
order=1 # order2 means we will consider second derivatives
solver="qss1"
initConditions=[1.0,2.0]
jacobian=[-2.5 20.0; 0.0 7.0]

settings = ModelSettings(initialTime,finalTime,dQmin,dQrel,order,solver,initConditions,jacobian)
println(settings.solver)

simulator = QSS_simulator(settings);
QSS_simulate(simulator);