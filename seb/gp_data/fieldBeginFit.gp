#!/usr/bin/env gnuplot5

fName = "'gradient2d_-44.dat'"

# functions depending on index
# linear g(x):

x0 = -0.19
fStr(N) = sprintf("f%d(x) = a%d * (x-x0)**2 + b%d", N, N, N)

# fits depending on index
fFitStr(N) = sprintf("fit [-0.19:-0.18] f%d(x) %s i %d u 3:1 via a%d, b%d", N, fName, N, N, N)

a = 8.
b = 0.1

# Set initial values
setVarStr(N) = sprintf("a%d = 8.; b%d = 0.1", N, N)

# Execute fits
do for [i=0:19] {
	eval( setVarStr(i) )
	eval( fStr(i) )
	eval( fFitStr(i) )
}

# Export fit parameters
set print "paramBegin.dat"
print "# a\tb"
expStr(N) = sprintf("print a%d, b%d", N, N, N)

do for [i=0:19] {
	eval(expStr(i))
}
print ""
unset print

# Plot
plotStr = "plot "
do for [i=0:19] {
    plotStr = plotStr . sprintf("f%d(x) t '', 'gradient2d_-44.dat' i %d u 3:1 w l t ''%s ", i, i, (i == 19) ? "" : ", ")
}

eval(plotStr)

# plot 'gradient2d_-44.dat' i 19 u 3:1 w l t 'Data',\
# 	g(x) t 'linear',\
#	f(x) t 'curve'

pause -1

