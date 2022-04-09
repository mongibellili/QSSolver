
using qss
initialTime=0.0
finalTime=5.0
dQmin=0.01
dQrel=0.5
order=1
solver="qss1"
initConditions=[1.0,2.0]
jacobian=[10.0 20.0; 30.0 40.0]

settings = ModelSettings(initialTime,finalTime,dQmin,dQrel,order,solver,initConditions,jacobian)
println(settings.solver)

simulator = QSS_simulator(settings);
QSS_simulate(simulator);