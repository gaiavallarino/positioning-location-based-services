from math import cos, sin
import numpy as np


def calcTrajectory(P0 , a_x, a_y, omegaz, epoch):
    # Initialize and compute velocities and delta positions in body frame

    pos = [P0]
    vx = [0]
    vy = [0]
    alpha = [0]
    for t in range(1,len(epoch)):
        deltat = (epoch[t] - (epoch[t-1]))
        vx.append(vx[t-1] + a_x[t] * deltat)  
        dx = vx[t] * deltat + (1/2) * a_x[t] * deltat**2
        ac = omegaz[t]*vx[t]
        resay = a_y[t] - ac
        vy.append(vy[t-1] + resay*deltat)
        dy = vy[t] * deltat + (1/2) * a_y[t] * deltat**2

        alpha.append(alpha[t-1] +omegaz[t] * deltat)
        R = np.matrix([ [ cos(alpha[t]), sin(alpha[t])],
                        [-sin(alpha[t]), cos(alpha[t])]])
        delta = np.array([[dx],[dy]])
        delta =  R * delta
        pos.append(np.array(pos[t-1] + delta).tolist())
    return pos

    # Implement a loop to compute X velocities and delta positions in body frame
    # y is constrained to zero because the cart is on a rail: only centrifugal, no skidding
    # clean apparent centrifugal from Y acceleration
    # compute, in case, skidding velocity and displacement
    
    # Implement a loop to compute for each epoch alpha, R(alpha), rotate Dx from body to
    # inertial and update intertial coordinates
    
    # Return the structure containing the
