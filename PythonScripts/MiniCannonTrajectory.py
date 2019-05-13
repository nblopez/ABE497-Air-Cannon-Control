#Python Conversion of Trajectory Calculation

#Required Modules
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

import time


#AeroDynFunction
def aerodyn(x_state, t_sim , K, g):
	#print(x_state)


	xp = [x_state[1], -K*x_state[1]*np.sqrt(x_state[1]**2+x_state[3]**2),
	    	x_state[3], -K*x_state[3]*np.sqrt(x_state[1]**2+x_state[3]**2)-g]

	return xp


start = time.time()

#Constant Values
g = float(input('g (m/s^2) [9.81] = ') or 9.81)
k = float(input('k [1.4] = ') or 1.4)
Mu_air = float(input('Viscosity of Air (mu) [1.82 * 10**-5] = ') or 1.82 * 10**-5)
Rho_air = float(input('Density of Air (rho, kg/m^3) [1.18] = ') or 1.18)


#Object Chars
Cd = float(input('Cd of Object [0.445] = ') or 0.445)
Obj_Mass = float(input('Mass of Object (g) [2.7] = ') or 2.7) / 1000
Obj_Diam = float(input('Diameter of Object (mm) [40] = ') or 40) / 1000

#Cannon Chars
Barrel_Length = float(input('Length of Cannon Barrel (m) [0.5] = ') or 0.5)
Barrel_Diam = float(input('Diameter of Cannon Barrel (mm) [42] = ') or 42)
Barrel_Thickness = float(input('Thickness of Cannon Barrel Wall (mm) [6] = ') or 6) / 1000

Res_Length = float(input('Length of Reservoir (m) [0.6] = ') or 0.6)
Res_Diam = float(input('Diameter of Reservoir (m) [0.1524] = ') or 0.1524)

Air_Vol_PostFire = Res_Diam**2 * np.pi/4 * Res_Length
Air_Vol_PreFire = Air_Vol_PostFire - (Barrel_Diam**2 * np.pi/4 * Barrel_Length)

#Cannon Setup
Barrel_Angle = float(input('Firing Angle (deg) [45] = ') or 45)
Barrel_Angle_rad = Barrel_Angle * np.pi/180
Barrel_Exit_Height = Barrel_Length * np.sin(Barrel_Angle_rad)


Res_Pressure_psi = float(input('Pressure in Reservoir (psi) [50] = ') or 50)
Res_Pressure_pa = Res_Pressure_psi * 6894.76

#Ambient Air Parameters
Temp_Outside_f = float(input('Ambient Air Temperature (f) [65] = ')  or 65)
Temp_Outside_K = (Temp_Outside_f - 32) * 5/9 + 273.15

#Initial Calculations
#Exit_Vel_Array = []
#Max_Dist_Check = 0.0


#Begin Traj Calcs

Pressure_PostFire = (Res_Pressure_pa * Air_Vol_PreFire**k)/Air_Vol_PostFire**k
Temp_PostFire = Temp_Outside_K * (Pressure_PostFire/Res_Pressure_pa)**((k-1)/k)

Work = ((Pressure_PostFire * Air_Vol_PostFire) - (Res_Pressure_pa * Air_Vol_PreFire))/(1-k)
print(Work)
Vel_Exit = np.sqrt((2/Obj_Mass)*Work)

Vy_Exit = Vel_Exit * np.sin(Barrel_Angle_rad)
Vx_Exit = Vel_Exit * np.cos(Barrel_Angle_rad)

#Cross-Sectional Area of Object
A_p = np.pi/4 * Obj_Diam**2

K = (Cd * Rho_air * A_p) / (2 * Obj_Mass)

Vel_Term = np.sqrt(g/K)


#Define Initial Conditions for odeint solver
t = np.linspace(0,10,1000)

init_cond = [0, Vx_Exit, Barrel_Exit_Height, Vy_Exit]

print(init_cond)

#Solve For Trajectory
sol = odeint(aerodyn, init_cond, t, args=(K,g))
#print(sol)
#np.savetxt('soln.csv', sol)

#Row in the solution where y goes below 0m
index = np.argmax(sol[:,2]<0)

x_dist = sol[:index,0]
x_vel  = sol[:index,1]
y_dist = sol[:index,2]
y_vel  = sol[:index,3]
time_index   = t[:index]


FinalTrajSlope = (x_dist[-1] - x_dist[-2])/(y_dist[-1] - y_dist[-2])

MaxDist = x_dist[-1] - y_dist[-2] * FinalTrajSlope
MaxHeight = max(y_dist)
AverageVel = np.sqrt(np.mean(abs(x_vel))**2 + np.mean(abs(y_vel))**2)
print("Projected Distance (m):       ", round(MaxDist, 2))
print("Projected Maximum Height (m): ", round(MaxHeight, 2))
print("Projected Flight Time (s):    ", round(np.amax(time_index), 2))
print("Projected Ave Velocity (m/s): ", round(AverageVel, 2))

'''
#Plot Trajectory
fig, (ax_dist, ax_vel) = plt.subplots(nrows=2)



ax_dist.plot(x_dist, y_dist, label='Trajectory')
ax_dist.set_title('Trajectory Plots')
ax_dist.set_xlabel('Horizontal Distance (m)')
ax_dist.set_ylabel('Vertical \n Distance (m)')
ax_dist.grid()

ax_vel.plot(time_index, x_vel, label='X Velocity')
ax_vel.plot(time_index, y_vel, label='Y Velocity')
ax_vel.set_xlabel('Flight Time (s)')
ax_vel.set_ylabel('Velocity (m/s)')
ax_vel.grid()
ax_vel.legend(loc = 'best')

plt.show()
'''
end = time.time()
print("Calculation Time (s):         ", round((end - start), 2))
