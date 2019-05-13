#Python Conversion of Trajectory Calculation

#Required Modules
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import csv
from time import sleep



#AeroDynFunction
def aerodyn(x_state, t_sim , K, g):
	#print(x_state)

	
	xp = [x_state[1], -K*x_state[1]*np.sqrt(x_state[1]**2+x_state[3]**2),
	    	x_state[3], -K*x_state[3]*np.sqrt(x_state[1]**2+x_state[3]**2)-g]
 
	return xp

def default_input():
	Cd = 0.4
	#Obj_Mass = 0.057
	Obj_Diam = 0.065
	return [Cd, Obj_Diam]
		
	
#Constant Values
g = 9.81
k = 1.4
Mu_air = 1.82 * 10**-5
Rho_air = 1.18

user_default = input('Press "y" and then "Enter" if you would like to enter object parameters? ')
if user_default == 'y':
	print("\nNOTE: All values are assumed to be greater than zero and numerical.\n")
	sleep(1)
	Cd = float(input('Please enter the Cd value for your object: '))
	#Obj_Mass = float(input('Please enter the mass for your object (kg): '))
	Obj_Diam = float(input('Please enter the diameter for your object (m): '))
	
	print("\nProceeding with calculations....")
	
	
else:
	print("\nUsing default values and proceeding with calculations....")
	[Cd, Obj_Diam] = default_input()


#Cannon Chars
Barrel_Length = 1.676
Barrel_Diam = 0.076
Barrel_Thickness = 0.006

Res_Length = 1.829
Res_Diam = 0.165

Air_Vol_PostFire = Res_Diam**2 * np.pi/4 * Res_Length
Air_Vol_PreFire = Air_Vol_PostFire - (Barrel_Diam**2 * np.pi/4 * Barrel_Length)

#Cannon Setup
Barrel_Angle = 45.0
Barrel_Angle_rad = Barrel_Angle * np.pi/180
Barrel_Exit_Height = Barrel_Angle_rad * Barrel_Length + 1

Res_Pressure_psi = 100.00
Res_Pressure_pa = Res_Pressure_psi * 6894.76

#Ambient Air Parameters
Temp_Outside_f = 68.00
Temp_Outside_K = (Temp_Outside_f - 32) * 5/9 + 273.15


#Begin Traj Calcs

Pressure_PostFire = (Res_Pressure_pa * Air_Vol_PreFire**k)/Air_Vol_PostFire**k
Temp_PostFire = Temp_Outside_K * (Pressure_PostFire/Res_Pressure_pa)**((k-1)/k)

Work = ((Pressure_PostFire * Air_Vol_PostFire) - (Res_Pressure_pa * Air_Vol_PreFire))/(1-k)

A_p = np.pi/4 * Obj_Diam**2

	



#Define Initial Conditions for odeint solver
t = np.linspace(0,40,10000)

MaxDist_Check = 0
Mass_Check = 0

for Obj_Mass in np.linspace(0.005,5,50):

	Vel_Exit = np.sqrt((2/Obj_Mass)*Work)
	K = (Cd * Rho_air * A_p) / (2 * Obj_Mass)
	Vel_Term = np.sqrt(g/K)

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
		Mass_Check = Obj_Mass
		
	#Plot Trajectory
	plt.plot(x_dist, y_dist)

print("Max Distance (m): %f"%MaxDist_Check)
print("Optimal Mass (kg): %f"%Mass_Check)

	

#plt.legend(loc='best')
plt.text(1,1, "Maximum Distance (m): ~%i \nOptimal Mass (kg): ~%.3f"%(MaxDist_Check,Mass_Check),bbox = dict(facecolor='gray', alpha = 0.75))
plt.title("Plot of Trajectories given variation in object mass")
plt.xlabel("Distance across ground from cannon (m)")
plt.ylabel("Height above ground (m)")
plt.grid()
plt.show()




