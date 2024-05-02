
import numpy as np
from numpy import array, sqrt, sum, repeat, zeros, newaxis, copy, pi, arange, linspace, clip, real, abs, ones, exp, argwhere, arctan
import matplotlib.pyplot as plt
plt.style.use('default')
from scipy import special

uM = 1e-6 #micro metre
nM = 1e-9 #nano metre
epsilon = 1*uM
LAMDA = 532*nM
W0 = 5*uM
Z_R = pi*1*(W0**2)/LAMDA
K0 = 2*pi/LAMDA


def radial_distance(plane_array_in):
    return sqrt(sum(plane_array_in**2, axis=2))

# print(radial_distance(plane))
N = 100
X = linspace(-50, 50, N)
Y = linspace(-50, 50, N)
Z = zeros(X.shape)
GLOBAL_EXTENT = [X[0], X[-1], Y[0], Y[-1]]

# Delta arrays for gradient
DX = []
DY = []
DZ = []
Delta_step = 1/N
# Z = zeros(n)
input_plane = []
for x in X:
    x_y = []
    DX_y = []
    DY_y = []
    DZ_y = []
    for y in Y:
        x_y_z = [x, y, 0]
        x_y.append(x_y_z)
        DX_y.append([Delta_step, 0, 0])
        DY_y.append([0, Delta_step, 0])
        DZ_y.append([0, 0, Delta_step])
    input_plane.append(x_y)
    DX.append(DX_y)
    DY.append(DY_y)
    DZ.append(DZ_y) 
input_plane = array(input_plane)
DX = array(DX)
DY = array(DY)
DZ = array(DZ)
# print(input_plane)


def get_radial_distance(plane_array_in):
    x = plane_array_in[:,:, 0].T
    y = plane_array_in[:,:, 1].T
    r = sqrt(x**2 + y**2)
    return r.T

def do_dot_product(plane_one, plane_two):
    x1, y1, z1 = plane_one[:,:,0].T, plane_one[:,:,1].T, plane_one[:,:,2].T
    x2, y2, z2 = plane_two[:,:,0].T, plane_two[:,:,1].T, plane_two[:,:,2].T
    dot_x = x1*x2
    dot_y = y1*y2
    dot_z = z1*z2
    return (dot_x + dot_y + dot_z).T    


def get_magnitude(plane_array_in):
    x = plane_array_in[:,:, 0].T
    y = plane_array_in[:,:, 1].T
    z = plane_array_in[:,:, 2].T
    r = sqrt(x**2 + y**2 + z**2)
    return r.T

def get_beam_width(plane_array_in, w0, zR):
    z = plane_array_in[:, :, 2]
    return (w0*sqrt(1 + (z/zR)**2))

def get_wavefront_curvature(plane_array_in, zR):
    z = plane_array_in[:, :, 2]
    return z*(1 + (zR/z)**2)

def get_gouy_phase(plane_array_in, zR):
    z = plane_array_in[:, :, 2]
    return arctan(z/zR)

def Hermite_n(x, n=0):
    ''' 
    x-> input variables
    n-> order of Hermite Polynomial (default = 0)
    '''
    h_n = special.hermite(n, monic=False)
    return h_n(x)

def get_refractive_index(plane_array_in, n0=1.3, A=sqrt(184615384.61538467)):
    # x = plane_array_in[:,:,0].T
    # y = plane_array_in[:,:,1].T
    # z = plane_array_in[:,:,2].T
    # r = sqrt(x**2 + y**2)
    # n = n0*(1 - 0.5*((A*r)**2))
    # return n.T

    r = get_radial_distance(plane_array_in)
    return n0*(1 - 0.5*(A*r)**2)

def plot_vector(vector_array_in, vector_name, lims = (-50, 50), cmap = 'viridis'):
    plt.figure()
    plt.title(f'{vector_name} in x')
    plt.imshow(vector_array_in[:, :, 0].T, extent=GLOBAL_EXTENT, origin='lower', cmap=cmap)
    plt.xlim(lims)
    plt.ylim(lims)
    plt.colorbar()

    plt.figure()
    plt.title(f'{vector_name} in y')
    plt.imshow(vector_array_in[:, :, 1].T, extent=GLOBAL_EXTENT, origin='lower', cmap=cmap)
    plt.xlim(lims)
    plt.ylim(lims)  
    plt.colorbar()

    plt.figure()
    plt.title(f'{vector_name} in z')
    plt.imshow(vector_array_in[:, :, 2].T, extent=GLOBAL_EXTENT, origin='lower', cmap=cmap)
    plt.xlim(lims)
    plt.ylim(lims)  
    plt.colorbar()

    plt.figure()
    plt.title(f"magnitude (|{vector_name}|)")
    plt.imshow(get_magnitude(vector_array_in).T, extent=GLOBAL_EXTENT, origin='lower',cmap=cmap )
    plt.xlim(lims)
    plt.ylim(lims)  
    plt.colorbar()
    plt.show()

def vectorized_gaussian_beam(plane_array_in, E_0_i=2, mode=(0, 0), refractive_index=1, l0 = LAMDA, w0 = W0):
    l, m = mode
    x = plane_array_in[:,:,0].T
    y = plane_array_in[:,:,1].T
    zR = (pi*refractive_index*w0**2)/l0
    r = get_radial_distance(plane_array_in).T
    w = get_beam_width(plane_array_in, w0=w0, zR=zR).T
    hx = Hermite_n(sqrt(2)*x/w, l)
    hy = Hermite_n(sqrt(2)*y/w, m)
    return (E_0_i*hx*hy*(w0/w)*exp(-(r/(1*w))**2)).T

def vectorized_psi(plane_array_in):
    x = plane_array_in[:, :, 0].T
    y = plane_array_in[:, :, 1].T
    z = plane_array_in[:, :, 2].T
    r_sq = x**2 + y**2
    R = get_wavefront_curvature(plane_array_in, zR=Z_R)
    return (K0*z + K0*(r_sq/(2*R)) - get_gouy_phase(plane_array_in, zR=Z_R)).T

def vectorized_k(plane_array_in):
    n = get_refractive_index(plane_array_in).T
    R = get_wavefront_curvature(plane_array_in, zR=Z_R).T
    x = plane_array_in[:,:,0].T
    y = plane_array_in[:,:,1].T
    z = plane_array_in[:,:,2].T
    kx = K0*x/R
    ky = K0*y/R
    kz = sqrt((K0*n)**2 - kx**2 - ky**2)
    return array([kx, ky, kz]).T

def vectorized_s_hat(plane_array_in):
    k_vec = vectorized_k(plane_array_in)
    
    k_vec_x = k_vec[:,:,0].T
    k_vec_y = k_vec[:,:,1].T
    k_vec_z = k_vec[:,:,2].T
    k_mag = get_magnitude(k_vec).T  

    s_hat_x = k_vec_x/k_mag
    s_hat_y = k_vec_y/k_mag
    s_hat_z = k_vec_z/k_mag
    return array([s_hat_x, s_hat_y, s_hat_z]).T

def vectorized_T(plane_array_in):
    
    n = get_refractive_index(plane_array_in).T
    s_hat = vectorized_s_hat(plane_array_in).T
    return (n*s_hat).T

def gradient_index(plane_array_in, n0=1.3, A=sqrt(184615384.61538467)):
    x = plane_array_in[:,:,0].T
    y = plane_array_in[:,:,1].T
    z = (plane_array_in[:,:,2]*0).T
    o = ones(z.shape)
    r = sqrt(x**2 + y**2) 
    a = argwhere(r >= 50*uM)
    o[a[:, 0], a[:, 1]] = 0
    return -(n0*(A**2)*array([x, y, z])*o).T
    
medium_in_plane = copy(input_plane*uM)
medium_in_plane[:,:, 2] = (160*uM*ones(medium_in_plane[:, :, 2].shape)).T

gauss_field = vectorized_gaussian_beam(medium_in_plane, mode=(0,0))

plt.figure()
im = plt.imshow(gauss_field.T, extent=GLOBAL_EXTENT, cmap='gist_rainbow', origin='lower')
plt.xlim(-20, 20)
plt.ylim(-20, 20)
ctr = plt.contour(input_plane[:,:, 0], input_plane[:, :, 1], gauss_field, 5, colors='black')
plt.clabel(ctr, inline='1')
cbr = plt.colorbar(im)
cbr.add_lines(ctr)


def step(plane_array_in, D_in, T_in, delta=epsilon):
    n = get_refractive_index(plane_array_in)
    n_e = repeat(n[:, :, newaxis], 3, axis=2)
    DdotT = repeat(do_dot_product(D_in, T_in)[:,:, newaxis], 3, axis=2)

    r_d = plane_array_in + (T_in/n_e)*(delta/4) + ((delta**2)/32)*((D_in/n_e) - T_in*(DdotT/n_e**3))
    r_p = plane_array_in + (T_in/n_e)*(delta/2) + ((delta**2)/8)*((D_in/n_e) - T_in*(DdotT/n_e**3))
    # D_p = gradient_index(r_p)
    # T_p = T_in + (delta/12)*(D_in + 4*gradient_index(r_d) + D_p)
    return r_d, r_p # D_p, T_p
    
    # rx, ry, rz = plane_array_in[:,:,0].T,plane_array_in[:,:,1].T,plane_array_in[:,:,2].T
    # T_x, T_y, T_z = T_in[:,:, 0].T, T_in[:,:,1].T, T_in[:,:,2].T
    # r_d_x = rx + (T_x/n)*(delta/2) + ((delta**2)/32)*(D)

def scalar_factor(plane_array_in):
    x = plane_array_in[:,:,0].T
    y = plane_array_in[:,:,1].T
    z = plane_array_in[:,:,2].T
    r = sqrt(x**2 + y**2)
    f = 1 - 0.01*(r-5*uM)/(5*uM)
    a = argwhere(r < 5*uM)
    f[a[:, 0], a[:,1]] = 1
    return f.T

def vectorized_propagate(plane_array_in, z_out = 260*uM, delta_s = epsilon):
    z_0 = plane_array_in[0,0,2]
    E_in = vectorized_gaussian_beam(plane_array_in)
    d_in = gradient_index(plane_array_in)
    t_in = vectorized_T(plane_array_in)
    psi_in = vectorized_psi(plane_array_in)
    T_in = copy(t_in)
    n_i = get_refractive_index(plane_array_in)
    u_hat_in = zeros(plane_array_in.shape)
    u_hat_in[:,:,1] = 1
    a_hat_in = zeros(plane_array_in.shape)
    a_hat_in[:,:,2] = 1
    a_hat_out = zeros(plane_array_in.shape)
    a_hat_out[:,:,2] = 1
    print(medium_in_plane[0,0, 2])
    c = 0
    while z_0 <= z_out:

        #r_i to r_prime
        t_in_copy = copy(t_in)
        r_double_prime, r_prime = step(plane_array_in, d_in, t_in, delta = delta_s)
        D_i, D_d_prime, D_prime = gradient_index(plane_array_in), gradient_index(r_double_prime), gradient_index(r_prime)
        T_prime = t_in + (delta_s/12)*(d_in + 4*D_d_prime + D_prime)
        # D_prime = gradient_index(r_prime)
        n_prime = get_refractive_index(r_prime)
        psi_in = psi_in + (K0*delta_s/4)*(n_i + 2*n_prime)
        f1 = -repeat((do_dot_product(u_hat_in, d_in)/n_i**2)[:,:, newaxis], 3, axis=2)*t_in
        f2 = -repeat((do_dot_product((u_hat_in + delta_s*f1/2), D_prime)/n_prime**2)[:,:, newaxis], 3,axis=2)*T_prime
        f3 = -repeat((do_dot_product((u_hat_in + delta_s*f2/2), D_prime)/n_prime**2)[:,:, newaxis], 3,axis=2)*T_prime
        
        # r_prime to r_i+1
        r_double_prime, plane_array_in = step(r_prime, D_prime, T_prime, delta = delta_s)
        d_in = gradient_index(plane_array_in)
        t_in = T_prime + (delta_s/12)*(gradient_index(r_prime) + 4*gradient_index(r_double_prime) + gradient_index(plane_array_in))
        n_i = get_refractive_index(plane_array_in)
        f4 = -repeat((do_dot_product((u_hat_in + delta_s*f3), d_in)/n_i**2)[:,:, newaxis], 3,axis=2)*t_in
        u_hat_in = u_hat_in + (delta_s/6)*(f1 + 2*f2 + 2*f3 + f4)
        psi_in = psi_in + (K0*delta_s/4)*get_refractive_index(plane_array_in)
        z_0 = plane_array_in[0,0,2]
        f_sc = scalar_factor(plane_array_in)
        E_in = f_sc*E_in*sqrt(5*do_dot_product(t_in_copy, a_hat_in)/do_dot_product(t_in, a_hat_out))
        # c += 1
        print(plane_array_in[0, 0, 2]*1e6)        
    # E_out = E_in*sqrt(5*do_dot_product(T_in, a_hat_in)/do_dot_product(t_in, a_hat_out))
    E_out = E_in
    E_final = repeat((E_out*exp(1j*psi_in))[:,:,newaxis], 3, axis=2)*u_hat_in
    # return z_0, u_hat_in, psi_in
    return E_final


E_OUT = vectorized_propagate(medium_in_plane, z_out=260*uM, delta_s=epsilon/10)

E_show = abs(real(E_OUT))
plot_vector(E_show, "E-field", lims=(-20, 20), cmap='gist_rainbow')