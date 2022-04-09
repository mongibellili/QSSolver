mutable struct QSS_simulator
    settings :: ModelSettings
    qssData :: QSSdata
    qssTime :: QSStime
    qssModel :: QSSmodel
    function QSS_simulator(settings :: ModelSettings )
        simulator = new()
        simulator.settings=settings
        simulator.qssData= QSSdata(settings.states, settings.initConditions)
        simulator.qssTime= QSStime(settings.states, settings.initialTime)
        simulator.qssModel= QSSmodel(settings.jacobian)
        simulator
    end
end
function QSS_simulate(qssSimulator:: QSS_simulator)
    println("........started simulation.........")
    integrator = QSS_integrator(qssSimulator);
    INT_integrate(integrator, qssSimulator);
end
