
import pynbody
import numpy as np
import sys


"""


Extracts halo particles in a spherical region in a
Gadget snapshot. Takes in the snapshot path as 
first argument, 6D coordinates from halo finder, 
and writes the ASCII file in GalactICS units to 
the last argument.


"""

data = pynbody.load(sys.argv[1])
a = float(data.properties['time'])
X0 = float(sys.argv[2])
Y0 = float(sys.argv[3])
Z0 = float(sys.argv[4])
VX0 = float(sys.argv[5])
VY0 = float(sys.argv[6])
VZ0 = float(sys.argv[7])
out = sys.argv[8]


masses = []
posx   = []
posy   = []
posz   = []
velx   = []
vely   = []
velz   = []

i = 0
x_dm  = np.array(data.dm['pos'].T[0], dtype=np.float32)
y_dm  = np.array(data.dm['pos'].T[1], dtype=np.float32)
z_dm  = np.array(data.dm['pos'].T[2], dtype=np.float32)
vx_dm = np.array(data.dm['vel'].T[0], dtype=np.float32)
vy_dm = np.array(data.dm['vel'].T[1], dtype=np.float32)
vz_dm = np.array(data.dm['vel'].T[2], dtype=np.float32)
m_dm  = np.array(data.dm['mass'], dtype=np.float32)
print("Reading data...\n")
for x,y,z,vx,vy,vz,m in zip(x_dm,y_dm,z_dm, vx_dm, vy_dm, vz_dm, m_dm):
    r = np.array([x - X0, y - Y0, z - Z0])
    if (np.linalg.norm(r) < 250.):
        masses.append(m*4.302)
        posx.append(r[0]*a)
        posy.append(r[1]*a)
        posz.append(r[2]*a)
        velx.append(vx / (a * 100.))
        vely.append(vy / (a * 100.))
        velz.append(vz / (a * 100.))
    i = i + 1

print("Shifting arrays...\n")
masses = np.array(masses)
velx   = np.array(velx) - VX0/(100. * a)
vely   = np.array(vely) - VY0/(100. * a)
velz   = np.array(velz) - VZ0/(100. * a)


print("Writing Xhalo of length + " + str(len(posx))  + "...\n")
f = open(out, "w")
f.write(str(len(posx)) + " 0.000000\n")
for m,x,y,z,vx,vy,vz in zip(masses, posx, posy,\
                                posz, velx, vely, velz):
    f.write(str(float(m)) + " " + str(float(x)) + " " +\
            str(float(y)) + " " + str(float(z)) + \
            " " + str(float(vx)) + " " + str(float(vy)) +\
            " " + str(float(vz)) + "\n")

print("Done.")
f.close()
