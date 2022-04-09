mutable struct QSS_integrator
    settings :: ModelSettings
    qssData :: QSSdata
    qssTime :: QSStime
    qssModel ::QSSmodel
    function QSS_integrator(qssSimulator:: QSS_simulator )
        integrator = new()
        integrator.settings=qssSimulator.settings
        integrator.qssData= qssSimulator.qssData
        integrator.qssTime= qssSimulator.qssTime
        integrator.qssModel= qssSimulator.qssModel
        integrator
    end
end
function INT_integrate(qssI::QSS_integrator,qssSimulator:: QSS_simulator)
    println("........started integration.........")
end