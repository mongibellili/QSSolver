# struct that glues the other 3 structs settings, data, time, model
# information flow is done through it
mutable struct QSS_simulator
    settings :: ModelSettings
    qssData :: QSSdata
    qssTime :: QSStime
    qssModel :: QSSmodel
    function QSS_simulator(settings :: ModelSettings )
        simulator = new()
        simulator.settings=settings
        simulator.qssData= QSSdata(settings.states,settings.order)
        simulator.qssTime= QSStime(settings.states, settings.initialTime)
        simulator.qssModel= QSSmodel(settings.jacobian)
        simulator
    end
end
function QSS_simulate(qssSimulator:: QSS_simulator)
    #println("........started simulation.........")
    QSS_integrate(qssSimulator);
end
