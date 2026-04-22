# set terminal pdfcairo color size 20cm,20cm font 'Helvetica, 34'
set terminal pdfcairo color size 20cm,20cm font 'Helvetica, 34'

set output 'disc24.pdf'

# set xrange [-1.0:1.0]
# set yrange [-0.1:0.6]
set size ratio 1.0

# set ytics
set xtics -0.005, 0.0025, 0.005

# set lmargin at screen 0.15
# set rmargin at screen 0.95
set bmargin at screen 0.14
# set tmargin at screen 0.9

# set key opaque

set ylabel "~L‾ _0"
set xlabel "{/Symbol b}-{/Symbol b}_c"

set style fill transparent solid 0.6 # partial transparency
set style fill noborder # no separate top/bottom lines

set key at 0.0048, 0.2 spacing 1.2 font 'Helvetica, 29' # -0.0005, 0.38

# set title "~{m}{0.4\\136}_{ q}=0.4"

plot "../spline_cpp_eps2/L0_24c_mass0p1000_1.dat" ti "~{m}{0.4\\136}_{ q}=0.1" lc 7 pointtype 8, \
     "../spline_cpp_eps2/L0_24c_mass0p1000_1.dat" using 1:($2-$3):($2+$3) with filledcurves lc 7 ti "stat error", \
     "../spline_cpp_eps2/L0_24c_mass0p1000_2.dat" notitle lc 7 pointtype 8, \
     "../spline_cpp_eps2/L0_24c_mass0p1000_2.dat" using 1:($2-$3):($2+$3) with filledcurves lc 7 notitle, \
     "../spline_cpp_eps2/L0_24c_mass0p2000_1.dat" ti "~{m}{0.4\\136}_{ q}=0.2" lc 4 pointtype 4, \
     "../spline_cpp_eps2/L0_24c_mass0p2000_1.dat" using 1:($2-$3):($2+$3) with filledcurves lc 4 ti "stat error", \
     "../spline_cpp_eps2/L0_24c_mass0p2000_2.dat" notitle lc 4 pointtype 4, \
     "../spline_cpp_eps2/L0_24c_mass0p2000_2.dat" using 1:($2-$3):($2+$3) with filledcurves lc 4 notitle, \
     "../spline_cpp_eps2/L0_24c_mass0p3000_1.dat" ti "~{m}{0.4\\136}_{ q}=0.3" lc 3 pointtype 3, \
     "../spline_cpp_eps2/L0_24c_mass0p3000_1.dat" using 1:($2-$3):($2+$3) with filledcurves lc 3 ti "stat error", \
     "../spline_cpp_eps2/L0_24c_mass0p3000_2.dat" notitle lc 3 pointtype 3, \
     "../spline_cpp_eps2/L0_24c_mass0p3000_2.dat" using 1:($2-$3):($2+$3) with filledcurves lc 3 notitle, \
     "../spline_cpp_eps2/L0_24c_mass0p4000_1.dat" ti "~{m}{0.4\\136}_{ q}=0.4" lc 6 pointtype 6, \
     "../spline_cpp_eps2/L0_24c_mass0p4000_1.dat" using 1:($2-$3):($2+$3) with filledcurves lc 6 ti "stat error", \
     "../spline_cpp_eps2/L0_24c_mass0p4000_2.dat" notitle lc 6 pointtype 6, \
     "../spline_cpp_eps2/L0_24c_mass0p4000_2.dat" using 1:($2-$3):($2+$3) with filledcurves lc 6 notitle

# plot "Veff32_mq0.4_P.dat" ti "{/Symbol b}=11.04775>{/Symbol b}_c" lc 4 pointtype 4 w e, \
#      "Veff32_mq0.4.dat" ti "{/Symbol b}=11.04275={/Symbol b}_c" lc 3 pointtype 3 w e, \
#      "Veff32_mq0.4_M.dat" ti "{/Symbol b}=11.03775<{/Symbol b}_c" lc 6 pointtype 6 w e

