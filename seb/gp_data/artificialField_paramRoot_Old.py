# Constants
FIELDRING_WIDTH = 9.652e-3
REFLECTORINNERRAD = 183.2356e-3

# Fit parameter
# Distance of field shaping rings
b = 16.8656e-3 # 0.00261

# Phase shift of S(z)
z_0 = -0.75 * b - 18.2753e-3 - FIELDRING_WIDTH/2. # -0.0028
# Center of h(r, z)
z0 = -0.19 
# Start of f(r, z)
z1 = -0.02 # b/2. - FIELDRING_WIDTH/2. - 18.2753e-3 
# Start of h(r, z)
z2 = -0.176
h_z0 = -0.19

# Parameter for f(r, z)
# =====================
# m(r)
c_m = 1.5
c_m2 = 0. # 0.0408885891201855 # -0.0177060260082422
d_m2 = 0. # 0.0099479324509953

s_m = REFLECTORINNERRAD - 0.005
s_m_ = (d_m2 - s_m * (c_m - c_m2)) / c_m + REFLECTORINNERRAD
d_m = c_m * (s_m_ - REFLECTORINNERRAD)

# t(r)
m_t = 1.00109952487923

# A(r)
t_h = 0.60000000
c_A = 5.50000000
d_A = REFLECTORINNERRAD - 0.01
# c_A2 = 0.000143566044902527 # -0.00292957115068489
# d_A2 = 0. # -3.19103731564881e-05

c_A2 = 0. # -0.00292957115068489
d_A2 = 0. # -3.19103731564881e-05

# Parameter for h(r, z)
# =====================
r0 = 0.163856
a_t = 50. # 136.273

# Parameter for g(r, z)
# =====================
A_gs = 0. # 0.01
s_g = 0.03 # 0.019

a_C = 3.e-2 
b_C = -0.15
c_C = -0.05 

