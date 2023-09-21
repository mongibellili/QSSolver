set title "mliqss_test"
set ylabel "State Variables"
set xlabel "Time"
set xrange [0:1.0]
set grid
plot "/home/mongi/qss-solver/output/mliqss_test/x1.dat" with lines title "x1", "/home/mongi/qss-solver/output/mliqss_test/x2.dat" with lines title "x2" 
