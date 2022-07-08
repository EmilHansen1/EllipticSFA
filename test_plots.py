# %% IMPORTS

from attr import NOTHING
from matplotlib.colors import LogNorm
import numpy as np
import matplotlib.pyplot as plt

# %% PARAMETERS AND FUNCTIONS

Up = 0.55
Ip = 0.579
omega = 0.057
N = 4
CEP = 0
ellip = np.pi/2
A0 = np.sqrt(2*Up)
phi = CEP
epsilon = ellip

def envelope(t):
    return np.sin(omega * t / (2 * N))**2

def e_field(t):
    return 2 * np.sqrt(Up) * envelope(t) * omega * np.array([np.sin(omega * t + CEP) * np.cos(ellip / 2), -np.cos(omega * t + CEP) * np.sin(ellip / 2), 0])

def a_field(t):
    return 2 * np.sqrt(Up) * envelope(t) * np.array([np.cos(omega * t + CEP) * np.cos(ellip / 2), np.sin(omega * t + CEP) * np.sin(ellip / 2), 0])



def p2(t, p_vec): 
    px, py, pz = p_vec
    return (px+A0*np.sin(omega*t/(2*N))**2 * np.cos(omega*t + phi)*np.cos(epsilon/2))**2 + (py + A0*np.sin(omega*t/(2*N))**2 * np.sin(omega*t+phi)*np.sin(epsilon/2))**2 + pz**2 

def S(t, p_vec): 
    px, py, pz = p_vec
    cp = np.cos(phi)
    sp = np.sin(phi)
    co = np.cos(omega * t)
    so = np.sin(omega * t)
    ce = np.cos(epsilon)
    se = np.sin(epsilon)
    con = np.cos(omega * t / N)
    son = np.sin(omega * t / N)

    front_fac = 1 / (64 * N**4 * omega - 80 * N**2*omega + 16*omega)
    term1 = 16 * (N-0.5) * A0**2 * N**2 * (cp * sp * co**2 + so * (cp**2 - 0.5)*co - sp*cp/2) * (N+0.5) * ce * con**2
    term2 = 4 * A0 * N * con * (((-4*cp**2 * ce + 2*ce) * co**2 + 4*cp*sp*ce*co*so + 2*cp**2*ce + N**2 - ce -1) *
            (N-0.5) * A0 * (N+0.5) * son - 8*(2*(N-0.5) * px * (so*cp + co*sp) * (N+0.5)*np.cos(epsilon/2) - 2*(N-0.5) * 
            (co*cp - so*sp)*py*(N+0.5)*np.sin(epsilon/2) + A0*(N-1) * (cp*sp*co**2 + so*(cp**2-0.5)*co - sp*cp/2)*(N+1)*ce)*N)
    
    term3 = -16*A0*N*son * (-4*(N-0.5) * (co*cp - so*sp) * px * (N+0.5)*np.cos(epsilon/2) - 4*(N-0.5)*py* 
            (so*cp + co*sp) * (N+0.5)*np.sin(epsilon/2) + A0*(N-1)*((-cp**2*ce + ce/2)*co**2 + cp*sp*ce*co*so + 
            cp**2 * ce/2 + N**2 - ce/4 - 1/4)*(N+1))

    term4 = 64 * (N-0.5) * (sp * (N**2-1)*co + cp*(N**2-1)*so + sp) * A0 * px * (N+0.5)*np.cos(epsilon/2)
    term5 = -64*(N-0.5)*A0*(cp*(N**2-1)*co + (1-N**2)*sp*so + cp) * py * (N+0.5) * np.sin(epsilon/2) 
    term6 = 16 * (N-0.5) * A0**2 *(N**2 - 3/2) *sp * cp * (N+0.5) * ce * co**2
    term7 = 16 * (N-0.5) * A0**2 * (N**2-3/2) * so * (cp**2 - 0.5) * (N+0.5) * ce * co
    term8 = -8*(A0**2 * ce * sp * (N**2-3/4) * cp - (3*(N-0.5)*(A0**2 + 16*px**2/3 + 16*py**2/3 + 16*pz**2/3)*omega*(N+0.5)*t)/2) * (N-1)*(N+1)

    return Ip * t + 0.5 * front_fac * (term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8)


p_vec = [0.3, 0.2, 0.]
t_list = np.linspace(0, 2*np.pi * N/omega, 100)
res = []

for ti in t_list: 
    int_list = np.linspace(0, ti, 1000)
    num_int = np.trapz([p2(tj, p_vec) for tj in int_list], int_list)
    res.append(num_int)

plt.plot(t_list, res)
plt.plot(t_list, S(t_list, p_vec))
plt.show()


def get_poly_coeffs(p):
    
    ec = np.exp(1j*CEP)
    ec2 = ec**2
    rtUp = np.sqrt(Up)
    px, py, pz = p
    A = np.sqrt(2*Up)
    C = np.cos(ellip/2)
    S = np.sin(ellip/2)
    S2 = S**2
    C2 = C**2

    #coeffs[0:3] = [(-px + 1j*py)*A0*ec/8, (px - 1j*py)*A0*ec/4, (-px + 1j*py)*A0*ec/8]
    #coeffs[N-1:N+4] += [Up/16, -Up/4, Ip + 3*Up/8 + (px**2 + py**2 + pz**2)/2, -Up/4, Up/16]
    #coeffs[2*N:2*N+3] += [-(px + 1j*py)*A0/(8*ec), (px + 1j*py)*A0/(4*ec), -(px + 1j*py)*A0/(8*ec)]
    if N == 2:
        coeffs = np.array([(C2 - S2)*Up*ec2/64, (S2 - C2)*Up*ec2/16, 3*(C2 - S2)*Up*ec2/32, # 0 - 2
            (S2 - C2)*Up*ec2/16 + (-C*px + 1j*S*py)*A*ec/8, (C2 + S2)*Up/32 + (C2 - S2)*Up*ec2/64 + (C*px - 1j*S*py)*A*ec/4, # 3 - 4
            (C2 + S2)*Up/8 + (-C*px + 1j*S*py)*A*ec/8, (px**2 + py**2 + pz**2)/2 + Ip + 3*(C2 + S2)*Up/16, # 5 - 6
            -(C2 + S2)*Up/8 - (C*px + 1j*S*py)*A/(8*ec), (C2 + S2)*Up/32 + (C2 - S2)*Up/(64*ec2) + (C*px + 1j*S*py)*A/(4*ec), # 7 - 8
            (-C2 + S2)*Up/(16*ec2) - (C*px + 1j*S*py)*A0/(8*ec), 3*(C2 - S2)*Up/(32*ec2), # 9 - 10
            (-C2 + S2)/(16*ec2), (C2 - S2)*Up/(64*ec2)], dtype=complex) # 10 - 12
    if N > 3:
        coeffs = np.zeros(4 * N + 4 + 1, dtype=complex)
        coeffs[0:5] = [(C2 - S2)*Up*ec2/64, (S2 - C2)*Up*ec2/16, 3*(C2 - S2)*Up*ec2/32, (S2 - C2)*Up*ec2/16, (C2 - S2)*Up*ec2/64]
        coeffs[N:N+5] = [0, (-C*px + 1j*S*py)*A*ec/8, (C*px - 1j*S*py)*A*ec/4, (-C*px + 1j*S*py)*A*ec/8, 0]
        coeffs[2*N:2*N+5] = [(C2 + S2)*Up/32, -(C2 + S2)*Up/8, (px**2 + py**2 + pz**2)/2 + Ip + 3*(C2 + S2)*Up/16, -(C2 + S2)*Up/8, (C2 + S2)*Up/32]
        coeffs[3*N:3*N+5] = [0, -(C*px + 1j*S*py)*A/(8*ec), (C*px + 1j*S*py)*A/(4*ec), -(C*px + 1j*S*py)*A/(8*ec), 0]
        coeffs[4*N:4*N+5] = [(C2 - S2)*Up/(64*ec2), (S2 - C2)*Up/(16*ec2), 3*(C2 - S2)*Up/(32*ec2), (S2 - C2)*Up/(16*ec2), (C2 - S2)*Up/(64*ec2)]
    return 0

def get_times(p):
    coeffs = get_poly_coeffs(p)
    z_roots = np.polynomial.polynomial.polyroots(coeffs)
    ts = 1j * N * np.log(z_roots) / omega
    ts = [t for t in ts if np.imag(t) > 0]
    ts = [t if np.real(t) > 0 else t + 2*np.pi*N/omega for t in ts]
    return ts

ts = get_times([0.2, 0.2, 0])


# %%

t_lst = np.linspace(0, 2 * np.pi * N / omega, 1000)
A_x, A_y, A_z = np.array([a_field(t) for t in t_lst]).T
E_x, E_y, E_z = np.array([e_field(t) for t in t_lst]).T
E_xs, E_ys, E_zs = np.array([e_field(np.real(t)) for t in ts]).T

# %%

plt.plot(omega*A_x, omega*A_y)
plt.plot(E_x, E_y)
plt.plot(E_xs, E_ys, 'ko')

# %%

plt.plot(t_lst, E_x)
plt.plot(np.real(ts), [e_field(np.real(t))[0] for t in ts], 'ko')

# %%

plt.plot(np.real(ts), np.imag(ts), 'ko')

# %%
t_re = np.linspace(0, 2*np.pi*N/omega, 300)
t_im = np.linspace(0, 70, 200)
S_mat = np.zeros((len(t_re), len(t_im)), dtype=complex)
for i, tr in enumerate(t_re):
    for j, it in enumerate(t_im):
        S_mat[i, j] += np.angle(np.exp(1j*S(tr + 1j*it, [0.2, 0.2, 0.0])))

 
# %%

plt.ylim(0, 100)
plt.imshow(np.real(np.flip(S_mat.T, 0)), extent=(0, 2*np.pi*N/omega, 0, 70), cmap='twilight', interpolation='bicubic', aspect='auto')
plt.plot(np.real(ts), np.imag(ts), 'ko')
# %%
