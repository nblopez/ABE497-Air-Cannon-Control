#Python Conversion of Trajectory Calculation

#Required Modules
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from os import system



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

Cd = 0.4

'''user_default = input('Press "y" and then "Enter" if you would like to enter object parameters? ')
if user_default == 'y':
	print("\nNOTE: All values are assumed to be greater than zero and numerical.\n")
	sleep(1)
	Cd = float(input('Please enter the Cd value for your object: '))
	Obj_Mass = float(input('Please enter the mass for your object (kg): '))
	Obj_Diam = float(input('Please enter the diameter for your object (m): '))
	
	print("\nProceeding with calculations....")
	
	
else:
	print("\nUsing default values and proceeding with calculations....")
	[Cd, Obj_Mass, Obj_Diam] = default_input()
'''
Barrel_Thickness = 0.006
Res_Pressure_psi = 100
Res_Pressure_pa = Res_Pressure_psi * 6894.76
#Define Initial Conditions for odeint solver
t = np.linspace(0,40,100)

MaxDist_Check = 0.0
Opt_Mass = 0.0
Opt_Angle = 0.0
Opt_Barrel_Length = 0.0
Opt_Barrel_Diam = 0.0
Opt_Res_Diam = 0.0
Opt_Res_Pressure = 0.0
Opt_Temp = 0.0
Opt_Obj_Diam = 0.0
Optimal_Vals = []

loop_num = 0
for Barrel_Length in np.linspace(10,15,2):

	Res_Length = Barrel_Length + 0.15	
	
	for Barrel_Diam in np.linspace(0.125, 0.25,5):
		Obj_Diam = Barrel_Diam
		A_p = np.pi/4 * Obj_Diam**2
		for Res_Diam in np.linspace(0.26,2,10):

			Air_Vol_PostFire = Res_Diam**2 * np.pi/4 * Res_Length
			Air_Vol_PreFire = Air_Vol_PostFire - (Barrel_Diam**2 * np.pi/4 * Barrel_Length)

		
		

			Pressure_PostFire = (Res_Pressure_pa * Air_Vol_PreFire**k)/Air_Vol_PostFire**k
			Work = ((Pressure_PostFire * Air_Vol_PostFire) - (Res_Pressure_pa * Air_Vol_PreFire))/(1-k)
			
			for Temp_Outside_f in np.linspace(0,100,5):
				
				Temp_Outside_K = (Temp_Outside_f - 32) * 5/9 + 273.15
				Temp_PostFire = Temp_Outside_K * (Pressure_PostFire/Res_Pressure_pa)**((k-1)/k)
				
				for Barrel_Angle in range(28,40):
				
					Barrel_Angle_rad = Barrel_Angle * np.pi/180
					Barrel_Exit_Height = Barrel_Length * np.sin(Barrel_Angle_rad)
			

							
					for Obj_Mass in np.linspace(0.02,10,10):
					
						K = (Cd * Rho_air * A_p) / (2 * Obj_Mass)
						Vel_Term = np.sqrt(g/K)
						Vel_Exit = np.sqrt((2/Obj_Mass)*Work)
					
						Vy = Vel_Exit * np.sin(Barrel_Angle_rad)
						Vx = Vel_Exit * np.cos(Barrel_Angle_rad)
					
					
						


						init_cond = [0, Vx, Barrel_Exit_Height, Vy]



						sol = odeint(aerodyn, init_cond, t, args=(K,g))


						index = np.argmax(sol[:,2]<0)

						x_dist = sol[:index,0]
						x_vel  = sol[:index,1]
						y_dist = sol[:index,2]
						y_vel  = sol[:index,3]
						time   = sol[:index]

						FinalTrajSlope = (x_dist[-1] - x_dist[-2])/(y_dist[-1] - y_dist[-2])

						MaxDist = x_dist[-1] - y_dist[-2] * FinalTrajSlope

						loop_num +=1
						if MaxDist >= MaxDist_Check:
							MaxDist_Check = MaxDist
							Opt_Mass = Obj_Mass
							Opt_Angle = Barrel_Angle
							Opt_Barrel_Length = Barrel_Length
							Opt_Barrel_Diam = Barrel_Diam
							Opt_Res_Diam = Res_Diam
							#Opt_Res_Pressure = Res_Pressure_psi
							Opt_Temp = Temp_Outside_f
							#Opt_Obj_Diam = Obj_Diam
							
							Optimal_Vals = [Opt_Mass, Opt_Angle, Opt_Barrel_Length, Opt_Barrel_Diam, Opt_Res_Diam, Opt_Temp]
						system('clear')
						print(loop_num, Optimal_Vals)

#Plot Trajectory
#plt.plot(x_dist, y_dist, label='%i Degrees'%Barrel_Angle)

print("Max Distance (m): %f"%MaxDist_Check)
print("Optimal Angle (deg): %f"%Angle_Check)
print("[Opt_Mass, Opt_Angle, Opt_Barrel_Length, Opt_Barrel_Diam, Opt_Res_Diam, Opt_Res_Pressure, Opt_Temp, Opt_Obj_Diam]")
print(Optimal_Vals)

	
'''
#plt.legend(loc='best')
plt.title("Plot of Trajectories given variation in launch angle")
plt.text(1,1, "Maximum Distance (m): ~%i \nOptimal Angle (deg): ~%.1f"%(MaxDist_Check,Angle_Check),bbox = dict(facecolor='gray', alpha = 0.75))
plt.xlabel("Distance across ground from cannon (m)")
plt.ylabel("Height above ground (m)")
plt.grid()
plt.show()
'''



