import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def heave_amp(t):
    return h0*np.sin(omega*t)

def pitch_amp(t):
    return theta0*np.cos(omega*t)

def heave_vel(t):
    return h0*omega*np.cos(omega*t)

def pitch_vel(t):
    return -1*theta0*omega*np.sin(omega*t)
   
def induced_vel(x, y, x0, y0, gamma):
    rr = (x-x0)**2 + (y-y0)**2 
    a = np.array([[0, 1], [-1, 0]])
    b = np.array([[x-x0], [y-y0]])
    ab = a@b
    return (gamma/(2*np.pi*rr))*ab

def rotate(theta):
    return np.array([
        [np.cos(theta), -np.sin(theta)],
        [np.sin(theta),  np.cos(theta)]
    ])

# Constants and Parameters
N = 5                # number of panels
c = 1.0              # chord length
dl = c/N            
dt = 0.02             # time step
T = 15               # total time
Nt = int(T / dt)      # number of time steps
U_inf = 1.0           # freestream velocity

# Initial geometry in local plate frame
x_endpoints = np.linspace(0, c, N+1)
z_endpoints = np.zeros_like(x_endpoints)
x_vortex = x_endpoints[:-1] + (c/N) * 0.25
x_colloc = x_endpoints[:-1] + (c/N) * 0.75
z_vortex = np.zeros(N)
z_colloc = np.zeros(N)

# Convert to (N, 2) coordinates
endpoint_coord_i = np.stack((x_endpoints, z_endpoints), axis=1)
vortex_coord_i = np.stack((x_vortex, z_vortex), axis=1)
colloc_coord_i = np.stack((x_colloc, z_colloc), axis=1)

# Define pitching axis (0.5c)
pitch_axis = np.array([0.5 * c, 0.0])

# Storage for positions over time
vortex_hist = []
colloc_hist = []
endpoint_history = []

# Kinematics
h0 = 0.5
theta0 = np.pi/30
reduced_freq = 3
non_dim_amp = h0/c
# omega = 2*np.pi*0.5
omega = (reduced_freq*U_inf)/c
# normal = np.array([0, 1]) # normal doesnt change with time in body frame
# tangent = np.array([1, 0]) # tangent also doesnt change with time in body frame
wake_pos_t = []
endpoint_history = []
gamma_wake = [] # contains the circulation strengths of each TEV
RHS = np.zeros(N+1)
gamma_body = np.zeros((Nt, N))
wake_hist = []
rho = 1.225
Lift = []
Drag = []
# Moment = []
time = []
print("kh = ", non_dim_amp*reduced_freq)
normal = np.array([0, 1]) # normal doesnt change with time
tangent = np.array([1, 0]) # tangent also doesnt change with time


for k in range(0,Nt):
    t = k*dt
    time.append(t)
    theta = -pitch_amp(t)
    h = heave_amp(t)
    x_vel = -U_inf
    z_vel = heave_vel(t)
    # x_shift = -U_inf * t
    
    # body frame velocities
    body_vel = rotate(theta) @ np.array([-x_vel, -z_vel])
    x_shift = 0
    R = rotate(theta)
    drag_at_j = 0

    # computing tangent and normal
    # normal = np.array([np.sin(theta), np.cos(theta)])
    # tangent = np.array([np.cos(theta), -np.sin(theta)])

    

    # Move to origin for rotation, then rotate, then move back
    endpoint_coord = (R @ (endpoint_coord_i - pitch_axis).T).T + pitch_axis
    vortex_coord = (R @ (vortex_coord_i - pitch_axis).T).T + pitch_axis
    colloc_coord = (R @ (colloc_coord_i - pitch_axis).T).T + pitch_axis

    # Add freestream x-shift and heave z-shift
    endpoint_coord[:, 0] += x_shift
    vortex_coord[:, 0] += x_shift
    colloc_coord[:, 0] += x_shift

    endpoint_coord[:, 1] += h
    vortex_coord[:, 1] += h
    colloc_coord[:, 1] += h

    endpoint_history.append(endpoint_coord.copy())
    vortex_hist.append(vortex_coord.copy())
    colloc_hist.append(colloc_coord.copy())

    # calculation of delta x and delta z (x and z distances between the trailing edge points)
    if k == 0:
        delta_x, delta_z = endpoint_coord[-1] - endpoint_coord_i[-1]
    else:
        delta_x, delta_z = endpoint_history[k][-1] - endpoint_history[k-1][-1]


    P = []
    # vortex release
    # number of wake vortices = k+1
    wake_x = endpoint_coord[5,0] - delta_x*0.3
    wake_z = endpoint_coord[5,1] - delta_z*0.3
    
    wake_pos_t.append([float(wake_x), float(wake_z)])
    
    # changing wake vortex positions based on wake rollup
    for abc in range(k):
        wake_pos_t[abc][0] += float((wake_vel_sum[abc,0] + U_inf)*dt)
        wake_pos_t[abc][1] += float(wake_vel_sum[abc,1]*dt)

    wake_hist.append([coord.copy() for coord in wake_pos_t])

    ind_vel = np.zeros((N, N+1, 2))
    infl_coeff = np.zeros((N+1, N+1)) 

    # induced velocity due to latest wake vortex of unit strength
    for p in range (N):
        vel_ind = induced_vel(colloc_coord[p,0], colloc_coord[p,1], wake_pos_t[-1][0], wake_pos_t[-1][1], 1)
        # print(u_ind, v_ind)
        ind_vel[p,N,0] = vel_ind[0, 0]
        ind_vel[p,N,1] = vel_ind[1, 0]
    
    # induced velocity due to other vortices of unit strength
    for v in range(0, N):
        for w in range(0,N):
            vel_ind = induced_vel(colloc_coord[v,0], colloc_coord[v,1], vortex_coord[w,0], vortex_coord[w,1], 1)
            ind_vel[v,w,0] = vel_ind[0, 0]
            ind_vel[v,w,1] = vel_ind[1, 0]

    # now calculate influence coefficients
    
    for v in range(0,N):
        for w in range(0,N+1):
            infl_coeff[v,w] = ind_vel[v,w,0]*normal[0] + ind_vel[v,w,1]*normal[1]

    for x in range(0, N+1):
        infl_coeff[N, x] = 1.0

    # influence coefficient matrix calculation is done
    # now calculate the RHS

    for c1 in range(0,N+1):
        if c1 != 5:
            RHS[c1] = -(body_vel[0]*normal[0] + body_vel[1]*normal[1]) 
            for prev_wake_vortex in range(0,k):
                wake_vel = induced_vel(colloc_coord[c1, 0], colloc_coord[c1, 1], wake_pos_t[prev_wake_vortex][0], wake_pos_t[prev_wake_vortex][1], gamma_wake[prev_wake_vortex]) # gamma k is calculated from the previous steps
                
                RHS[c1] += -1*(wake_vel[0,0]*normal[0] + wake_vel[1,0]*normal[1])
        else:
            if k == 0: # airfoil doesnt have circulation before 1st time step 
                RHS[c1] = 0
            else:
                RHS[c1] = sum(gamma_body[k-1])

    # RHS calculation done
    # now solve the equation

    ans = np.linalg.inv(infl_coeff) @ RHS
    
    for a in range(0, N+1):
        if a == N:
            gamma_wake.append(float(ans[a]))
        else:
            gamma_body[k, a] = ans[a]    

    # aerodynamic forces computation
    for j in range(N):
        gamma_sum = 0
        gamma_prev_sum = 0
        wwj = 0
        P.append(rho*(body_vel[0]*tangent[0] + body_vel[1]*tangent[1]))
        for something in range(0,k+1): # need to calculate the induced velocity due to wake vortices
            wake_vel = induced_vel(colloc_coord[j, 0], colloc_coord[j, 1], wake_pos_t[something][0], wake_pos_t[something][1], gamma_wake[something])
            P[j] += rho*(wake_vel[0,0]*tangent[0] + wake_vel[1,0]*tangent[1])
            wwj += wake_vel[1]
        
        P[j] *= (gamma_body[k, j]/dl)

        for ki in range(j+1):
            gamma_sum += gamma_body[k,ki]
            if k != 0:
                gamma_prev_sum += gamma_body[k-1,ki]
        
        time_derivative_gamma = (gamma_sum - gamma_prev_sum)/dt 
        term = rho*((wwj*gamma_body[k, j]) + (time_derivative_gamma*dl*np.sin(theta)))
        drag_at_j += term
        P[j] += rho*time_derivative_gamma
        
    Lift.append((sum(P)*dl*np.cos(theta))/(0.5*rho*U_inf*U_inf*c))
    Drag.append(drag_at_j/(0.5*rho*U_inf*U_inf*c))
    
    # wake rollup
    for wake_vort in range(k+1):
        if wake_vort == 0:
            wake_vel_sum = np.zeros((k+1, 2)) # sum of velocities induced by all vortices on a specific wake vortex 
        
        for n in range(N+1):
            if n == N:
                for other_wake_vort in range(k+1):
                    if other_wake_vort == wake_vort: continue
                    else:
                        vel = induced_vel(wake_pos_t[wake_vort][0], wake_pos_t[wake_vort][1], wake_pos_t[other_wake_vort][0], wake_pos_t[other_wake_vort][1], gamma_wake[other_wake_vort])
                        wake_vel_sum[wake_vort,0] += vel[0,0]
                        wake_vel_sum[wake_vort,1] += vel[1,0]
            else:
                # velocity induced on wake due to bound vortices
                vel = induced_vel(wake_pos_t[wake_vort][0], wake_pos_t[wake_vort][1], vortex_coord[n,0], vortex_coord[n,1], gamma_body[k, n])
                # print("vel=", vel)
                wake_vel_sum[wake_vort,0] += vel[0,0]
                wake_vel_sum[wake_vort,1] += vel[1,0]


    # print(ind_vel)
    # print(infl_coeff)
    # print (RHS)
    # print (gamma_body)
    # print(gamma_wake)
    # print(wake_vel_sum)
    # print(wake_pos_t)
    # print()
    # print("done")
    # print(Lift)

# print(wake_hist)   
# print(endpoint_history)
print("mean drag = ", sum(Drag[100:])/Nt)
print("mean lift = ", sum(Lift[100:])/Nt)
# print ("done")

fig, ax = plt.subplots()
ax.set_xlim(-0.5, 20)
ax.set_ylim(-2, 2)
ax.set_aspect(1.5)

wing_line, = ax.plot([], [], 'k-', lw=1.5)
wake_dots, = ax.plot([], [], 'r*', markersize=2)

def init():
    wing_line.set_data([], [])
    wake_dots.set_data([], [])
    return wing_line, wake_dots

def update(frame):
    # Update wing and wake vortices
    if frame < len(wake_hist):
        xdata = [pt[0] for pt in wake_hist[frame]]
        zdata = [pt[1] for pt in wake_hist[frame]]
        wake_dots.set_data(xdata, zdata)
        wing_x = [ep[0] for ep in endpoint_history[frame]]
        wing_z = [ep[1] for ep in endpoint_history[frame]]
        wing_line.set_data(wing_x, wing_z)
    return wing_line, wake_dots

ani = animation.FuncAnimation(
    fig, update, frames=Nt, init_func=init, blit=True, interval=50
)
# Optional: Save the animation as an MP4 or GIF
# ani.save("flapping_wake.mp4", fps=30, dpi=150)
# ani.save("flapping_wake.gif", writer='imagemagick', fps=30)

plt.title("Wake Rollup Visualization")
plt.xlabel("x")
plt.ylabel("z")
plt.grid(True)
plt.show()

# Cl plotting
plt.plot(time, Lift)
plt.title("Cl vs time")
plt.xlabel("time")
plt.ylabel("Cl")
plt.show()

# Cd plotting
plt.plot(time, Drag)
plt.title("Cd vs time")
plt.xlabel("time")
plt.ylabel("Cd")
plt.show()
