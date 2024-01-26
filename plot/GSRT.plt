set title "van Hove correlation of Li^+ ion"
set ylabel "4πr^2 G_s(r,t)"
set xlabel "r (Å)"

set xrange [:15]

set grid

plot "./dat/GSRT.dat" u 1:2 w l lw 3 lc "#BE3455" t "100 ps", \
     "./dat/GSRT.dat" u 1:3 w l lw 3 lc "#6868AB" t "500 ps", \
     "./dat/GSRT.dat" u 1:4 w l lw 3 lc "#939597" t "1000 ps", \
     "./dat/GSRT.dat" u 1:5 w l lw 3 lc "#F5DF4D" t "1500 ps", \
     "./dat/GSRT.dat" u 1:6 w l lw 3 lc "#34568B" t "2000 ps", \

