import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import orbit_info as oi
import tdp1
import tdp3
import standard_const as stdc
import ecc

earth_radius = 6378.0
earth_mu = 398600.0    #   km3/s2
thrust = 10e3
I_sp = 300   # s
g0 = stdc.g0   # m/s2
mass_0 = 2000.0
r1, r2, r3 = 436, 6083, 2529 # km
v1, v2, v3 = -7.340, -0.5125, 2.497   # km/s

def derivative(vec_sv, t):

    rx, ry, rz, vx, vy, vz, m = vec_sv
    r = np.array([rx, ry, rz])
    v = np.array([vx, vy, vz])
    norm_r = np.linalg.norm(r)
    norm_v = np.linalg.norm(v)
    ax, ay, az = -r*earth_mu/norm_r**3 + thrust*v/m/norm_v/1000.0
    dmdt = -thrust/I_sp/g0
    dydt= [vx, vy, vz, ax, ay, az, dmdt]
    return dydt

def orbit_prop(vec_sv, t):

    rx, ry, rz, vx, vy, vz, m = vec_sv
    r = np.array([rx, ry, rz])
    v = np.array([vx, vy, vz])
    norm_r = np.linalg.norm(r)
    norm_v = np.linalg.norm(v)
    ax, ay, az = -r*earth_mu/norm_r**3
    dmdt =0
    dydt = [vx, vy, vz, ax, ay, az, dmdt]
    return dydt


if __name__  == '__main__' :

  # Coasting period

    tspan1 = 89*60 # seconds
    dt1 = 0.1         # secs
    n_steps1 = int(np.ceil(tspan1 / dt1))
    r0 = [r1, r2, r3]
    v0 = [v1, v2, v3]
    m0 = [mass_0]

    print(oi.sv_to_oe(r0, v0))
    v_sv_op = r0 + v0 + m0
    path = np.zeros((n_steps1, 7))
    t1 = np.linspace(0, tspan1, n_steps1)
    path1 = odeint(orbit_prop, v_sv_op, t1)

    x_t1 = path1[:, 0]
    y_t1 = path1[:, 1]
    z_t1 = path1[:, 2]
    ecc1 = ecc.ecc(r0, v0)

    r_t1 = (x_t1 ** 2 + y_t1 ** 2 + z_t1 ** 2) ** 0.5
    rs1 = path1[-1]
    #   plt.plot(t1, x_t1)
    #   plt.plot(t1, y_t1)
    #   plt.plot(t1, z_t1)
    plt.plot(t1, r_t1, label='Radial distance' )
    plt.legend(loc='upper left')
    plt.xlabel('time(seconds)')
    plt.ylabel('Radial Distance (in km)')
    plt.show()
    tdp1.plot(path1)

    pre_burn_sv = path1[-1]  # state vector just before burn starts
    print('Pre-burn state vector: ', pre_burn_sv)
    print('Pre-burn Orbital Elements', oi.sv_to_oe(pre_burn_sv[0:3] , pre_burn_sv[3:6]))  #to check if coe remained same or not

#  Thrust on
    tspan = 120   # secs
    dt = 0.1        # secs
    n_steps = int(np.ceil(tspan/dt))
    t = np.linspace(0, tspan, n_steps)
    vec_sv0 = pre_burn_sv
    sol = odeint(derivative, vec_sv0, t)

    x_t = sol[:, 0]
    y_t = sol[:, 1]
    z_t = sol[:, 2]
    m_t = sol[:, 6]

    r_t = (x_t**2 + y_t**2 + z_t**2)**0.5
    rs = sol[-1]
    v1 = rs[3:6]
    ecc2 = ecc.ecc(rs[0:3], rs[3:6])
    plt.plot(t, x_t, label='R_x')
    plt.plot(t, y_t, label='R_y')
    plt.plot(t, z_t, label='R_z')
    plt.xlabel('time(seconds)')
    plt.ylabel('In km')
    plt.legend(loc='center')
    plt.show()

    plt.plot(t, r_t, label='Radial Distance')
    plt.xlabel('time(seconds)')
    plt.ylabel('In km')
    plt.legend(loc='center')
    plt.show()
    print('Post-burn State vector: ', rs[0:6])

    print('Post burn Orbital elements:', oi.sv_to_oe(rs[0:3], rs[3:6]))
    plt.plot(t, m_t, label='Mass')
    plt.xlabel('time(seconds)')
    plt.ylabel('In Kg')
    plt.legend(loc='upper right')
    plt.show()
    print('Post burn Mass:', rs[6])


    print('Fuel expended = ', 2000.0-rs[6],
          'Delta V non-impulsive (in m/s) =' ,
          (np.linalg.norm(v1) - np.linalg.norm(v0))*1000)

    # Coasting period 2
    sv2 = oi.sv_to_oe(rs[0:3], rs[3:6])
    tspan2 = 3 * 7000      # seconds
    dt2 = 1.0              # secs
    n_steps2 = 21000
    t2 = np.linspace(0, tspan2, n_steps2)
    path2 = odeint(orbit_prop, rs, t2)

    x_t2 = path2[:, 0]
    y_t2 = path2[:, 1]
    z_t2 = path2[:, 2]

    r_t2 = (x_t2 ** 2 + y_t2 ** 2 + z_t2 ** 2) ** 0.5
    rs2 = path2[-1]
    # plt.plot(t2, x_t2)
    # plt.plot(t2, y_t2)
    # plt.plot(t2, z_t2)
    plt.plot(t2, r_t2, label='Radial distance')
    plt.legend(loc='upper left')
    plt.xlabel('time(seconds)')
    plt.ylabel('Radial Distance (in km )')
    plt.show()

    tdp3.plot(path1, path, path2)

    apse_rot = np.arccos(np.dot(ecc1, ecc2)/np.linalg.norm(ecc1)/np.linalg.norm(ecc2))
    print('Apse line rotated by ', apse_rot*180/np.pi, ' degrees')
    def impulsive():
        ra1 = 6602.479985800101
        rp1 = 6601.7270536077
        ra2 = 9502.609203208358
        rp2 = 6603.060672045933

        h2= (2*earth_mu)**0.5 * (rp1*ra2)**0.5 / (rp1 + ra2)**0.5
        vp2 = h2/rp1
        delta_vp = vp2 - np.linalg.norm(v0)
        va2 = h2/ ra2
        h3 = (2*earth_mu)**0.5 * (rp2*ra2)**0.5 / (rp2 + ra2)**0.5
        va3= h3/ ra2
        delta_va = va3 - va2
        del_total = (delta_va + delta_vp) * 1000
        fuel_mass = mass_0 * (1 - np.exp(-del_total / g0 / I_sp))
        return fuel_mass, delta_va, delta_vp, del_total

    fuel_mass, delta_va, delta_vp, del_total = impulsive()
    print( 'Total Delta V for impulsive (in m/s)  = ', del_total ,
           ' Total fuel requirement (in KG) = ', fuel_mass )


