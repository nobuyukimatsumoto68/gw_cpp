set terminal pdfcairo enhanced color size 20cm,18cm font 'Helvetica, 34'
set output 'geff.pdf'

set logscale xy
set xrange [1000:0.1]
set yrange [50:130]

set xlabel "T [GeV]" # offset -4,0,0
set ylabel "g_{eff,SM}" # offset -1,0,0

plot "geff.dat" u 1:2 w l lw 3 notitle #  "g_{eff}(T)"

