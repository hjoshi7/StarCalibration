import os
import sys
import numpy
sys.path.append('/opt/homebrew/lib/python3.11/site-packages')
import astropy
from astropy.coordinates import SkyCoord, EarthLocation, ICRS, AltAz, solar_system_ephemeris, get_body
from astropy.time import Time
import numpy as np
import astropy.units as u
import pandas as pd
import math
from lmfit import Minimizer, Parameters, fit_report
from matplotlib import pyplot as plt

#System Coordi¡nates Conversions
def GpsToECEF(lat, long, alt):
    # alt in meters, out in km
    lat = lat * math.pi / 180
    long = long * math.pi / 180
    a = 6378.1
    b = 6356.8
    e = 1 - (b**2 / a**2)
    N = a / (math.sqrt(1.0 - (e * (math.sin(lat)**2))))
    cosLatRad = math.cos(lat)
    cosLongiRad = math.cos(long)
    sinLatRad = math.sin(lat)
    sinLongiRad = math.sin(long)
    x = (N + 0.001 * alt) * cosLatRad * cosLongiRad
    y = (N + 0.001 * alt) * cosLatRad * sinLongiRad
    z = ((b**2 / a**2) * N + 0.001 * alt) * sinLatRad
    return [x, y, z]

def GpstoENU(camera_position):
    # alt in meters
    # e, n, u returns eastward, northward and zenith distance from reference to i
    #lati, longi, alti, latr, longr, altr
    # Define the reference GPS coordinates Camera 1
    longr = camera_position.loc[0,'E']
    latr = camera_position.loc[0,'N']
    altr = camera_position.loc[0,'H'] #Radius Earth at sea level
    

    # Define the point GPS coordinates Camera 2
    longi= camera_position.loc[1,'E']
    lati = camera_position.loc[1,'N']
    alti = camera_position.loc[1,'H'] #Radius Earth at sea level    
    
    cosLatRad = math.cos(latr * math.pi / 180)
    cosLongRad = math.cos(longr * math.pi / 180)
    sinLatRad = math.sin(latr * math.pi / 180)
    sinLongRad = math.sin(longr * math.pi / 180)
    loci = GpsToECEF(lati, longi, alti)
    locr = GpsToECEF(latr, longr, altr)
    dx = loci[0] - locr[0]
    dy = loci[1] - locr[1]
    dz = loci[2] - locr[2]
    e = dx * (-sinLongRad) + dy * cosLongRad
    n = dx * (-sinLatRad) * cosLongRad + dy * (-sinLatRad) * sinLongRad + dz * cosLatRad
    u = dx * cosLatRad * cosLongRad + dy * cosLatRad * sinLongRad + dz * sinLatRad
    print(f"Eastward distance (e): {e:.3f} km, Northward distance (n): {n:.3f} km, Upward distance (u): {u:.3f} km")
    return e*1000, n*1000, u*1000

#camera parameters for the projection matrices
def configfile(camcnt, configfilename):
    # Read the Excel file into a pandas dataframe
    df = pd.read_excel(configfilename, sheet_name='camara_config')
    fx = [0,0]
    fy = [0,0]
    cx = [0,0]
    cy = [0,0]
    pitch = [0,0]
    roll = [0,0]
    yaw = [0,0]
    for i in range(camcnt):
        fx[i] = df.iloc[0, i+1]
        fy[i] = df.iloc[1, i+1]
        cx[i] = df.iloc[2, i+1]
        cy[i] = df.iloc[3, i+1]
        pitch[i] = df.iloc[4, i+1]
        roll[i] = df.iloc[5, i+1]
        yaw[i] = df.iloc[6, i+1]
    
    #origin_camera_system= pd.read_excel('ConfigFile.xlsx', sheet_name='origin_camera_system') #RO This may not be necessary
    camera_position=pd.read_excel(configfilename, sheet_name='camera_position')
    return fx, fy, cx, cy, pitch, roll, yaw,camera_position
    
#Matrix operations
def matrix_multiply(Cm,R):#this function multiplies two matrices in the given order
    result = [[sum(a*b for a,b in zip(Cm_row,R_col)) for R_col in zip(*R)] for Cm_row in Cm]
    return(result)

def invert_matrix(A):
    B = numpy.linalg.inv(A)
    return(B)

def projection_matrix(focusx,focusy,cx,cy,yaw,pitch,roll,tx,theight,ty):#this function gets the projection matrix P
    Cm = [[focusx, 0.0, cx],
         [0. ,focusy, cy],
         [0., 0.,1.]]

    R12 = [[math.cos(yaw),0 ,-math.sin(yaw)],
          [0, 1, 0],
          [math.sin(yaw),0,math.cos(yaw)]]

    R11 = [[1,0,0],
          [0,math.cos(pitch),math.sin(pitch)],
          [0,-math.sin(pitch),math.cos(pitch)]]

    R13 = [[math.cos(roll),math.sin(roll),0],
          [-math.sin(roll),math.cos(roll),0],
          [0, 0, 1]]
    T = [[tx],
        [theight],
        [ty]]
    RB = matrix_multiply(R11,R12)
    R = matrix_multiply(RB,R13)
    RT = -1*np.array(matrix_multiply(R,T))
    Rt = np.append(R,RT, axis=1)
    P = matrix_multiply(Cm,Rt)
    return(P) #Gets the projection matrix P  


def enu_to_cam_coordinate(dx, dy, dz, facing):
    # This routine is to translate x = E, y = N, z = zenith to camera coordinate system which is
    # x = baseline from pairing to reference, y = -zenith, z = normal to baseline
    # By default, assume N
    tx = -dx
    ty = -dz
    tz = dy
    if facing == 'S':
        tx = dx
        tz = -dy
    elif facing == 'W':
        tx = dy
        tz = -dx
    elif facing == 'E':
        tx = -dy
        tz = dx
    return tx, ty, tz #in meters


def get_orientation(dx,dy):
    facing = 'N' #by default
    phi = math.atan(abs(dy)/abs(dx))
    if dx < 0:
     if dy < 0:
          if phi <= math.pi/4:
               print("The system is facing North")
          else:
               print("The system is facing West")
               facing = 'W'
     else:
          if phi <= math.pi/4:
               print("The system is facing North")
          else:
               print("The system is facing East")
               facing = 'E'
    else:
        facing = 'S'
        if dy < 0:
            if phi <= math.pi/4:
                print("The system is facing South")
            else:
                print("The system is facing West")
                facing = 'W'
        else:
            if phi <= math.pi/4:
                print("The system is facing South")
            else:
                print("The system is facing East")
                facing = 'E'
    return(facing)

def arrange_translation_vector(dx, dy, dz, facing):
    # This routine is to translate x = E, y = N, z = zenith to camera coordinate system which is
    # x = baseline from pairing to reference, y = -zenith, z = normal to baseline
    # By default, assume N
    tx = dx
    ty = -dz
    tz = dy
    if facing == 'S':
        tx = -dx
        tz = -dy
    elif facing == 'W':
        tx = dy
        tz = dx
    elif facing == 'E':
        tx = -dy
        tz = -dx
    return [[tx], [ty], [tz],[1]]

def fundamental_matrix(P_R, P_P, T):
#Building the Fundamental matrix

    PP = matrix_multiply(numpy.transpose(P_P),invert_matrix(matrix_multiply(P_P,numpy.transpose(P_P))))
    f1 = numpy.array(matrix_multiply(P_R, T))
    f1mat = [[0,-f1.item(2),f1.item(1)],[f1.item(2),0,-f1.item(0)],[-f1.item(1),f1.item(0),0]]
    F = matrix_multiply(matrix_multiply(f1mat,P_R),PP)
    FF = numpy.divide(numpy.array(F), numpy.array(F).item(8)) #This is the final Fundamental matrix
    return (FF)


def project_to_cam(P, e, n, u, facing):
    Pos = []
    
    for i in range(len(e)):
        px,py,pz = enu_to_cam_coordinate(e[i],n[i],u[i],facing)
        Pos.append ([px, py, pz, 1]) #Nx4 matr of world positions
     
        
    proj = matrix_multiply(P, np.transpose(Pos))
    proj_x = np.divide(numpy.array(proj[0]),np.array(proj[2]))
    proj_y = np.divide(numpy.array(proj[1]),np.array(proj[2]))
    
    return(np.concatenate((proj_x,proj_y),axis=0)) # 2Nx1 vector of (xi,yi) 
      

#Reading matching pairs
def reading_matchingpairs(file_name):
        with open(file_name, 'r') as f:
        # initialize empty lists for each column
            dates = []
            times = []
            r_values = []
            p_values = []

            # iterate over each line in the file
            cnt = 0
            for line in f:
                # skip lines that do not start with a digit
                if not line[0].isdigit():
                    continue

                # split the line into individual values
                values = line.split()
         
                # convert the values to their appropriate types
                date = int(values[0])
                time = int(values[1])
                xr = (values[2])
                yr = (values[3])
                xp = (values[4])
                yp = (values[5])
                
                # append the values to their respective arrays
                dates.append(date)
                times.append(time)
                r_values.append(float(xr))
                r_values.append(float(yr))
                r_values.append(float(1))
                p_values.append(float(xp))
                p_values.append(float(yp))
                p_values.append(float(1))
                cnt = cnt+1


        ptsr = numpy.array(r_values)
        ptsp = numpy.array(p_values)
        return(ptsr.reshape((cnt,3)),ptsp.reshape((cnt,3)), cnt)
                
def reading_starpos(file_name):
    with open(file_name, 'r') as f:
        # initialize empty lists for each column
        stars = []
        dates = []
        times = []
        x_values = []
        y_values = []
        e_values = []
        n_values = []
        u_values = []

        line_counter = 0

        # iterate over each line in the file
        for line in f:
            # skip lines that do not start with a digit
            if line_counter == 0  or line_counter%2!=0:
                line_counter += 1
                continue
            # split the line into individual values
            values = line.split()
            
            # convert the values to their appropriate types
            star = values[1]
            date = int(values[2].split('T')[0].replace('-',''))
            time = int(values[2].split('T')[1].replace(':','').split('.')[0])

            x = float(values[3])
            y = float(values[4])
            e = float(values[5])
            n = float(values[6])
            u = float(values[7])

            # append the values to their respective arrays
            stars.append(star)
            dates.append(date)
            times.append(time)
            x_values.append(x)
            y_values.append(y)
            e_values.append(e)
            n_values.append(n)
            u_values.append(u)

            line_counter += 1
    return [x_values, y_values,e_values,n_values,u_values]



def enforce_epipolar_constraint(mpr,mpp,FF):
    epipolar_vector_cost = []        
        
    epiLine = matrix_multiply(mpr,FF)
    epiLine_proj = numpy.absolute(numpy.sum(numpy.multiply(mpp,epiLine),1))

    scale = numpy.square(epiLine)
    scale = numpy.sqrt(scale[:,0]+scale[:,1])
    res_err = numpy.array(numpy.divide(epiLine_proj,scale))
    
    return res_err
    
         
def enforce_projection_constraint(PR, PP, er, nr, ur, ep, np, up, facing ):
    
    xryr = project_to_cam(PR, er, nr, ur, facing)
    xpyp = project_to_cam(PP, ep, np, up, facing)
    
    return(xryr,xpyp)
    
#Initialize parameters for the optimization routine of lmfit    
def initialize_axes(yaw,pitch,roll,xr,yr,xp,yp,match_len):
    params = Parameters()
    params.add('yaw_r', value=yaw[0], min=yaw[0]-0.1, max=yaw[0]+0.1)
    params.add('yaw_p', value=yaw[1], min=yaw[1]-0.1, max=yaw[1]+0.1)
    params.add('pitch_r', value=pitch[0], min=pitch[0]-0.1, max=pitch[0]+0.1)
    params.add('pitch_p', value=pitch[1], min=pitch[1]-0.1, max=pitch[1]+0.1)
    params.add('roll_r', value=roll[0], min=roll[0]-0.1, max=roll[0]+0.1)
    params.add('roll_p', value=roll[1], min=roll[1]-0.1, max=roll[1]+0.1)
    
    zerovec = np.zeros(match_len)
    data = np.concatenate((xr,yr,xp,yp,zerovec), axis=0)
    
    return params, data

def optimization_func(pars,fx,fy,cx,cy,tx,ty,tz,facing,mptsr,mptsp,e_r,n_r,u_r,e_p,n_p,u_p, data=None):

    yawr, pitchr, rollr = pars['yaw_r'], pars['pitch_r'], pars['roll_r']
    yawp, pitchp, rollp = pars['yaw_p'], pars['pitch_p'], pars['roll_p']
    P_R=projection_matrix(fx[0],fy[0],cx[0],cy[0],yawr,pitchr,rollr,0,0,0)
    P_P=projection_matrix(fx[1],fy[1],cx[1],cy[1],yawp,pitchp,rollp,tx,ty,tz)
    T = [[tx],[ty], [tz], [1] ] 
    FF = fundamental_matrix(P_R, P_P, T)
    
    #Enforcing the epipolar constraint
    epip_err = enforce_epipolar_constraint(mptsr,mptsp,FF)
    #Enforcing the projection constraint
    #the output is the projected pixels x_r, y_r, x_p, y_p
    proj_r,proj_p = enforce_projection_constraint(P_R,P_P,e_r,n_r,u_r,e_p,n_p,u_p, facing)
    model = np.concatenate((proj_r,proj_p,epip_err),axis=0)
    if data is None:
        return model
    return (model-data)
    
def print_parameter_table(params):
    param_names = ['yaw_r', 'yaw_p', 'pitch_r', 'pitch_p', 'roll_r', 'roll_p']
    print('Parameter\tValue\tUncertainty\tBounds')
    for name in param_names:
        value = params[name].value
        uncertainty = params[name].stderr
        lower_bound = value - uncertainty
        upper_bound = value + uncertainty
        bound_str = f"[{lower_bound:.6f}, {upper_bound:.6f}]"
        print(f'{name}\t{value:.6f}\t{uncertainty:.6f}\t{bound_str}')

#Main part of the code
def main():
    
    ######Camera setup and config file
    camcnt = 2
    configfilename = 'ConfigFile.xlsx'
    # data files 
    matchingpairsfilename = 'matching_points.txt'
    starposfilename_r = 'starpos_r.txt'
    starposfilename_p = 'starpos_p.txt'
    
    #read initial camera parameterscamera positions for two cameras (reference and pair)
    fx, fy, cx, cy, pitch, roll, yaw, camera_position = configfile(camcnt,configfilename)

    #find the distance between cameras and the direction cameras are facing 
    x_dist,y_dist,z_dist= GpstoENU(camera_position)
    facing = get_orientation(x_dist,y_dist)  
    #find the translation vector
    tx,ty,tz = enu_to_cam_coordinate(x_dist,y_dist,z_dist,facing)

    ######Input vectors
    #Reading matching pairs --xr,yr,xp,yp-- 
    mpts_r, mpts_p, pts_cnt = reading_matchingpairs(matchingpairsfilename)
    #Reading xr,yr,xp,yp,E,N,U 
    x_r,y_r,e_r,n_r,u_r = reading_starpos(starposfilename_r)
    x_p,y_p,e_p,n_p,u_p = reading_starpos(starposfilename_p)

    #Optimization routine
    #Initialize parameters
    params,data = initialize_axes(yaw,pitch,roll,x_r,y_r,x_p,y_p,pts_cnt)
    #run optimization
    best_val = Minimizer(optimization_func, params, fcn_args=(fx,fy,cx,cy,tx,ty,tz,facing,mpts_r,mpts_p,e_r,n_r,u_r,e_p,n_p,u_p,), fcn_kws={'data': data})    
    best_out = best_val.leastsq()
    best_fit = optimization_func(best_out.params, fx,fy,cx,cy,tx,ty,tz,facing,mpts_r,mpts_p,e_r,n_r,u_r,e_p,n_p,u_p)
    
    print_parameter_table(best_out.params)

    #plt.plot(x, data, 'o', label='data'
    #plt.plot(x, fit, label='with analytical derivative')

    #Plotting the error
    residual_error = data - best_fit
    plt.plot(residual_error)
    plt.xlabel('Point index')
    plt.ylabel('Residual error')
    plt.title('Residual error plot')
    plt.show()
    
    #Printing the min, max, and average of residual error
    print(f"Min residual error: {np.min(residual_error)}")
    print(f"Max residual error: {np.max(residual_error)}")
    print(f"Average residual error: {np.mean(residual_error)}")

#Execution of the code
main()
