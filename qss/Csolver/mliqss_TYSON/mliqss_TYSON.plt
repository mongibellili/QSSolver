set title "mliqss_TYSON"
set ylabel "State Variables"
set xlabel "Time"
set xrange [0:25]
set grid
plot "/home/mongi/qss-solver/output/mliqss_TYSON/C2.dat" with lines title "C2", "/home/mongi/qss-solver/output/mliqss_TYSON/CP.dat" with lines title "CP", "/home/mongi/qss-solver/output/mliqss_TYSON/M.dat" with lines title "M", "/home/mongi/qss-solver/output/mliqss_TYSON/pM.dat" with lines title "pM", "/home/mongi/qss-solver/output/mliqss_TYSON/Y.dat" with lines title "Y", "/home/mongi/qss-solver/output/mliqss_TYSON/yP.dat" with lines title "yP" 
