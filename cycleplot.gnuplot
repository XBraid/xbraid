#
# Type the following command in the run directory (e.g., in 'examples'):
#   gnuplot ../cycleplot.gnuplot
#
set yrange [:] reverse
set offsets 0, 0, 1, 1
plot 'braid.out.cycle' using ($1-$2) with linespoints pt 5
pause -1 "Hit any key to continue"
