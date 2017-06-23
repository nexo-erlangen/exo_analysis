#!/usr/bin/env gnuplot5

# m Fit
m(x) = m_c * x + m_d # exp(m_c*(x-m_d))
m_c = 1.
m_d = 0.175

m2(x) = m_c2 * x + m_d2
m_c2 = 1.
m_d2 = 0. #0.175

fit [0.:0.18] m2(x) 'param_out.dat' u 1:2 via m_c2
fit [0.18:0.19] m(x) 'param_out.dat' u 1:2 via m_c, m_d

# t Fit
t(x) = t_m * x # + t_t
t_m = 1.
#t_t = 0.1

fit t(x) 'param_out.dat' u 1:3 via t_m #, t_t

# A Fit
A(x) = -A_c*((x - A_d)**2) + A_d2 #-exp(A_c * (x - A_d))
# A_A = 0.001
A_c = 5
A_d = 0.175

A2_(x) = A_c2*(x - A_d) + A_d2
A2(x) = A_d2 / A_d * x
A_c2 = 0.
A_d2 = 0.

fit [0.:0.175] A2_(x) 'param_out.dat' u 1:4 via A_c2, A_d2
fit [0.175:0.19] A(x) 'param_out.dat' u 1:4 via A_c, A_d

set multiplot layout 2, 2 rowsfirst
plot 'param_out.dat' u 1:2 w l t 'm(r)',\
	m(x) t '',\
	m2(x) t ''
plot 'param_out.dat' u 1:3 w l t 't(r)',\
	t(x) t ''
plot 'param_out.dat' u 1:4 w l t 'A(r)',\
	A(x) t '',\
	A2(x) t 'A2(r)'
plot 'param_out.dat' u 1:5 w l t 'b(r)'
unset multiplot

print ""
print "c_m = ", m_c
print "d_m = ", m_d
print "c_m2 = ", m_c2
print "d_m2 = ", m_d2
print ""

print "m_t = ", t_m
print ""

# print "t_t = ", t_t
# print "A_A = ", A_A
print "c_A = ", A_c
print "d_A = ", A_d
print "c_A2 = ", A_c2
print "d_A2 = ", A_d2

pause -1

