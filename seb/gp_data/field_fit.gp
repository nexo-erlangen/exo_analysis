#!/usr/bin/env gnuplot5

fName = "'gradient2d_-44.dat'"

# functions depending on index
# linear g(x):
gStr(N) = sprintf("g%d(x) = m%d * x + t%d", N, N, N)
# curve f(x):
fStr(N) = sprintf("f%d(x) = g%d(x) + A%d * sin(-(x - s%d)/b%d )", N, N, N, N, N)

# fits depending on index
# gFitStr(N) = sprintf("fit [-0.14:-0.06] g%d(x) %s i %d u 3:1 via m%d, t%d", N, fName, N, N, N)
fFitStr(N) = sprintf("fit [-0.14:-0.06] f%d(x) %s i %d u 3:1 via A%d, b%d, m%d, t%d", N, fName, N, N, N, N, N)

m = 0.1
t = 0.18
g(x) = m*x + t

A = 0.0008
b = 0.00260513
s = 1.
f(x) = g(x) + A*sin(-(x - s)/b )

# fit [-0.14:-0.06] g(x) 'gradient2d_-44.dat' i 19 u 3:1 via m, t
# fit [-0.14:-0.06] f(x) 'gradient2d_-44.dat' i 19 u 3:1 via A, b, s

# Set initial values
setVarStr(N) = sprintf("m%d = 0.1; t%d = 0.18; A%d = 0.0008; b%d = 0.0026; s%d = 1.", N, N, N, N, N)

# Execute fits
do for [i=0:19] {
	eval( setVarStr(i) )
	eval( gStr(i) )
	eval( fStr(i) )
	# eval( gFitStr(i) )
	eval( fFitStr(i) )
}

# Export fit parameters
set print "param.dat"
print "# m\tt\tA\tb\ts"
expStr(N) = sprintf("print m%d, t%d, A%d, b%d, s%d", N, N, N, N, N)

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

