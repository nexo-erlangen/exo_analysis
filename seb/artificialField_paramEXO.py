# Constants
FIELDRING_WIDTH = 9.652e-3
REFLECTORINNERRAD = 183.2356e-3

# Fit parameter
# Distance of field shaping rings
b = 16.8656e-3 # 0.00261

# Phase shift of S(z)
z_0 = -0.75 * b - 18.2753e-3 - FIELDRING_WIDTH/2. + 1.3359507223967966e-3 # -0.0028
# Center of h(r, z)
z0 = -0.19 
# Start of f(r, z)
z1 = -0.02 # b/2. - FIELDRING_WIDTH/2. - 18.2753e-3 
# Start of h(r, z)
z2 = -0.187
h_z0 = -0.19

# Parameter for f(r, z)
# =====================
# m(r)
c_m = -1.593037
c_m2 = -0.001
d_m2 = 0. 

s_m = 0.172 # REFLECTORINNERRAD - 0.01
d_m = s_m *(c_m2 - c_m) + d_m2
# c_m2 = 0.

# c_m_ = (c_m * REFLECTORINNERRAD + d_m) / (REFLECTORINNERRAD - s_m)
# d_m_ = c_m * REFLECTORINNERRAD + d_m - c_m_ * REFLECTORINNERRAD

# c_m  = c_m_
# d_m  = d_m_

# t(r)
m_t = 1.00109952487923

# A(r)
t_h = 1.177416
c_A = 7.671374 
d_A = s_m 

c_A2 = -0.0003
d_A2 = c_A2 * d_A

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

