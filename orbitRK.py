'''
orbitRK.py 
implement rk4 to solve orbit problem
'''

from visual import *
from math import *

G = 6.67e-11
AU = 1.5e11
year = 365.25*24.*60.*60.
dt = 1.0e4

scene.autoscale = 1

sun = sphere(pos= (0,0,0), radius= 6.9e8, mass = 2e30, color=color.yellow)
earth = sphere(pos=(AU,0,0), radius = 400*6.4e6, mass= 6e24, color=color.blue)
earth.vel = vector(0,1.4*pi*AU/year,1.4*pi*AU/year)
earth.trail = curve(color = color.blue)
def acc(earthpos):
    r = mag(earthpos - sun.pos)
    return -G*sun.mass*(earthpos - sun.pos)/r**3
def rk4(earth):
    k1v = acc(earth.pos)*dt
    k1x = earth.vel*dt

    k2v = acc(earth.pos + k1x/2.0)*dt
    k2x = (earth.vel + k1v/2.0)*dt

    k3v = acc(earth.pos + k2x/2.0)*dt
    k3x = (earth.vel + k2v/2.0)*dt

    k4v = acc(earth.pos + k3x)*dt
    k4x = (earth.vel + k3v)*dt

    earth.vel += (k1v + 2.0*k2v + 2.0*k3v +k4v)/6.0
    earth.pos += (k1x + 2.0*k2x + 2.0*k3x +k4x)/6.0

ecc = 0
ecc_s = str("%.3f"%ecc)
v = str("%.3f"%mag(earth.vel))
vel_label = label(pos = earth.pos, text = "vel: "+v+"\necc: "+ecc_s)
varrow = arrow(pos=earth.pos, axis= earth.vel, leght=7e9, color=color.cyan)

while True:
    if scene.kb.keys:           #check if there is a keyboard event waiting
        s = scene.kb.getkey()   #obtain the keyboard entry
        if s == 'b':
            earth.vel *= 1.1
        if s == 'r':
            earth.vel *= 0.9
    rk4(earth)                  #call runge-kutta to calculate velocity and position

    earth.trail.append(pos=earth.pos)   #see path of earth
    varrow.pos = earth.pos
    varrow.axis = earth.vel*1.0e6
    
    r = mag(earth.pos - sun.pos)
    L= earth.mass*cross(earth.pos, earth.vel)       #calculate angular momentum
    E = 0.5*earth.mass*mag2(earth.vel) - G*earth.mass*sun.mass/r
    ecc = sqrt(1.0+2*mag2(L)*E/((G*sun.mass**2)*earth.mass**3))

    ecc_s = str("%.3f"%ecc)
    v = str("%.3f"%mag(earth.vel))
    vel_label.pos = earth.pos
    vel_label.text = "vel: "+v+"\necc: "+ecc_s
    rate(200)






    
    
                  
