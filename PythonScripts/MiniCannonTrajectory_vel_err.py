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
g = 9.81
k = 1.4
Mu_air = 1.82 * 10**-5
Rho_air = 1.18


#Object Chars
Cd = 0.445
Obj_Mass = 2.7 / 1000
Obj_Diam = 40 / 1000

#Cannon Chars
Barrel_Length = 0.5
Barrel_Diam = 42 / 1000
Barrel_Thickness = 0.006

Res_Length = 0.6
Res_Diam = 0.1524

Air_Vol_PostFire = Res_Diam**2 * np.pi/4 * Res_Length
Air_Vol_PreFire = Air_Vol_PostFire - (Barrel_Diam**2 * np.pi/4 * Barrel_Length)

#Cannon Setup
Barrel_Angle = 32.0
Barrel_Angle_rad = Barrel_Angle * np.pi/180
Barrel_Exit_Height = Barrel_Length * np.sin(Barrel_Angle_rad)


Res_Pressure_psi = 302.0
Res_Pressure_pa = Res_Pressure_psi * 6894.76

#Ambient Air Parameters
Temp_Outside_f = 68.00
Temp_Outside_K = (Temp_Outside_f - 32) * 5/9 + 273.15

#Initial Calculations
#Exit_Vel_Array = []
#Max_Dist_Check = 0.0


#Begin Traj Calcs

Pressure_PostFire = (Res_Pressure_pa * Air_Vol_PreFire**k)/Air_Vol_PostFire**k
Temp_PostFire = Temp_Outside_K * (Pressure_PostFire/Res_Pressure_pa)**((k-1)/k)

Work = ((Pressure_PostFire * Air_Vol_PostFire) - (Res_Pressure_pa * Air_Vol_PreFire))/(1-k)
A_p = np.pi/4 * Obj_Diam**2
K = (Cd * Rho_air * A_p) / (2 * Obj_Mass)
Vel_Term = np.sqrt(g/K)

fig, (ax_dist, ax_vel) = plt.subplots(nrows=2)

#Define Initial Conditions for odeint solver
t = np.linspace(0,10,10000)
Vel_Exit_real = np.sqrt((2/Obj_Mass)*Work)


for fact_err in np.linspace(0.9,1.1,3):
	
	Vel_Exit = Vel_Exit_real * fact_err

	Vy = Vel_Exit * np.sin(Barrel_Angle_rad)
	Vx = Vel_Exit * np.cos(Barrel_Angle_rad)






	init_cond = [0, Vx, Barrel_Exit_Height, Vy]
	#print(init_cond)

	#Solve For Trajectory
	sol = odeint(aerodyn, init_cond, t, args=(K,g))
	#print(sol)
	#np.savetxt('soln.csv', sol)

	index = np.argmax(sol[:,2]<0)

	x_dist = sol[:index,0]
	x_vel  = sol[:index,1]
	y_dist = sol[:index,2]
	y_vel  = sol[:index,3]
	time_index   = t[:index]


	FinalTrajSlope = (x_dist[-1] - x_dist[-2])/(y_dist[-1] - y_dist[-2])

	MaxDist = x_dist[-1] - y_dist[-2] * FinalTrajSlope
	#print("Projected Distance @ %i m/s: "%(int(Vel_Exit)), MaxDist)
	#print("Projected Flight Time @ %i m/s: "%(int(Vel_Exit)), np.amax(time_index))


	plt_i = -1

	ax_dist.plot(x_dist, y_dist, label='Dist: %im'%(int(MaxDist)))
	#ax_vel.plot(time_index[:plt_i], x_vel[:plt_i], label='X Velocity: %i'%(int(Vel_Exit)))
	ax_vel.plot(time_index[:plt_i], y_vel[:plt_i], label='Dist %im'%(int(MaxDist)))
	
	
ax_dist.set_title('Trajectory Plots')
ax_dist.set_xlabel('Horizontal Distance (m)')
ax_dist.set_ylabel('Vertical \n Distance (m)')
ax_dist.legend(loc='best')
ax_dist.grid()


ax_vel.set_xlabel('Flight Time (s)')
ax_vel.set_ylabel('Velocity (m/s)')
ax_vel.grid()
ax_vel.legend(loc = 'best')

end = time.time()
print("Calculation Time: ", (end - start))

fig.tight_layout()
plt.show()




