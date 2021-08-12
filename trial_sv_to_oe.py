import standard_const as stdc
import orbit_info as oi

r = [stdc.Re + 500, 0, 0]    # km
v = [0, 7.2, 0]    #km/s

oi.sv_to_oe(r,v)

