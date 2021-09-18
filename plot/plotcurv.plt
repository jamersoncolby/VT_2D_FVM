# plotcurv.plt
set term x11 font "-*-helvetica-medium-r-*-*-14-*-*-*-*-*-*-*"
set title "Curvilinear Mesh"
set nokey
set grid
set xlabel "x"
set ylabel "y"
m="plot/curvgrid.txt"
plot m using 1:2 with linespoints



