#Python Conversion of Trajectory Calculation

#Required Modules
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time


#import csv



#AeroDynFunction
def aerodyn(x_state, t_sim , K, g):
	#print(x_state)

	
	xp = [x_state[1], -K*x_state[1]*np.sqrt(x_state[1]**2+x_state[3]**2),
	    	x_state[3], -K*x_state[3]*np.sqrt(x_state[1]**2+x_state[3]**2)-g]
 
	return xp
	
	


#Constant Values
g = 9.81
k = 1.4
Mu_air = 1.82 * 10**-5
Rho_air = 1.18


#Object Chars
Cd = 0.2
#Obj_Mass = 2.0
Obj_Diam = 0.25

#Cannon Chars
Barrel_Length = 4.0
Barrel_Diam = 0.254
Barrel_Thickness = 0.125

Res_Length = 3.0
Res_Diam = 1.0

Air_Vol_PostFire = Res_Diam**2 * np.pi/4 * Res_Length
Air_Vol_PreFire = Air_Vol_PostFire - (Barrel_Diam**2 * np.pi/4 * Barrel_Length)

#Cannon Setup
#Barrel_Exit_Height = 1.0
#Barrel_Angle = 32.0
#Barrel_Angle_rad = Barrel_Angle * np.pi/180

Res_Pressure_psi = 200.00
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






#Define Initial Conditions for odeint solver
t = np.linspace(0,40,500)

MaxDist_Check = 0
Angle_Check = 0
Obj_Mass_Check = 0.0

fig = plt.figure()
ax = fig.gca(projection='3d')

for Obj_Mass in np.linspace(0.05, 4, 20):

	Vel_Exit = np.sqrt((2/Obj_Mass)*Work)
	K = (Cd * Rho_air * A_p) / (2 * Obj_Mass)
	
	Vel_Term = np.sqrt(g/K)
	
	for Barrel_Angle in np.linspace(1,89,80):
		Barrel_Angle_rad = Barrel_Angle * np.pi/180
		Barrel_Exit_Height = Barrel_Angle_rad * Barrel_Length + 1
		
		Vy = Vel_Exit * np.sin(Barrel_Angle_rad)
		Vx = Vel_Exit * np.cos(Barrel_Angle_rad)

		

		init_cond = [0, Vx, Barrel_Exit_Height, Vy]
		#print(init_cond)

		#Solve For Trajectory
		sol = odeint(aerodyn, init_cond, t, args=(K,g))
		#print(sol)

		index = np.argmax(sol[:,2]<0)

		x_dist = sol[:index,0]
		x_vel  = sol[:index,1]
		y_dist = sol[:index,2]
		y_vel  = sol[:index,3]
		time   = sol[:index]

		FinalTrajSlope = (x_dist[-1] - x_dist[-2])/(y_dist[-1] - y_dist[-2])

		MaxDist = x_dist[-1] - y_dist[-2] * FinalTrajSlope
		#print(MaxDist)
		
		if MaxDist >= MaxDist_Check:
			MaxDist_Check = MaxDist
			Angle_Check = Barrel_Angle
			Obj_Mass_Check = Obj_Mass
			
		#Plot Trajectory
		ax.plot(x_dist, y_dist, Obj_Mass)
		
		
		


print("Max Distance (m): %f"%MaxDist_Check)
print("Optimal Angle (deg): %f"%Angle_Check)
print("Optimal Mass (kg): %f"%Obj_Mass_Check)

#print("Time Taken to Run Loop (s): %i"%(delt_time))

#plt.legend(loc='best')
#plt.grid()
plt.show()




