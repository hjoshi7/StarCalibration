""" """ import sys
import numpy
sys.path.append('/opt/homebrew/lib/python3.11/site-packages')
from astropy.coordinates import SkyCoord, EarthLocation, ICRS, AltAz, solar_system_ephemeris, get_body
from astropy.time import Time
import numpy as np
import astropy.units as u
import pandas as pd
import math
import matplotlib.pyplot as plt
import numpy as np

from lmfit import Minimizer, Parameters


    
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
    return e, n, u

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
    return tx*1000, ty*1000, tz*1000 #in meters


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
    
    for i in range(len(e)):
        px,py,pz = enu_to_cam_coordinate(e[i],n[i],u[i],facing)
        Pos = [[px], [py], [pz], [1]]
        Cam_pos = matrix_multiply(P, Pos)
        Camcorr[i][0] =  Cam_pos[0][0] / Cam_pos[2][0]
        Camcorr[i][1] =  Cam_pos[1][0] / Cam_pos[2][0]
    return Camcorr #Nx2 vector  (xi,yi) 
    
    

#Reading files
def reading_1(file_name):
        with open(file_name, 'r') as f:
        # initialize empty lists for each column
            dates_1 = []
            times_1 = []
            xr_values_1 = []
            yr_values_1 = []
            xp_values_1 = []
            yp_values_1 = []

            # iterate over each line in the file
            for line in f:
                # skip lines that do not start with a digit
                if not line[0].isdigit():
                    continue

                # split the line into individual values
                values = line.split()

                # convert the values to their appropriate types
                date = int(values[0])
                time = int(values[1])
                xr = int(values[2])
                yr = int(values[3])
                xp = int(values[4])
                yp = int(values[5])

                # append the values to their respective arrays
                dates_1.append(date)
                times_1.append(time)
                xr_values_1.append(xr)
                yr_values_1.append(yr)
                xp_values_1.append(xp)
                yp_values_1.append(yp)

        return [xr_values_1, yr_values_1, xp_values_1, yp_values_1]
                
def reading_2(file_name):
    with open(file_name, 'r') as f:
        # initialize empty lists for each column
        stars = []
        dates_2 = []
        times_2 = []
        xr_values_2 = []
        yr_values_2 = []
        xp_values_2 = []
        yp_values_2 = []
        e_values_2 = []
        n_values_2 = []
        u_values_2 = []

        line_counter = 0

        # iterate over each line in the file
        for line in f:
            # skip lines that do not start with a digit
            if line_counter == 0 or line_counter ==1:
                line_counter += 1
                continue
            # split the line into individual values
            values = line.split()

            # convert the values to their appropriate types
            star = values[0]
            date = int(values[1])
            time = int(values[2])
            xr = int(values[3])
            yr = int(values[4])
            xp = int(values[5])
            yp = int(values[6])
            e = float(values[7])
            n = float(values[8])
            u = float(values[9])

            # append the values to their respective arrays
            stars.append(star)
            dates_2.append(date)
            times_2.append(time)
            xr_values_2.append(xr)
            yr_values_2.append(yr)
            xp_values_2.append(xp)
            yp_values_2.append(yp)
            e_values_2.append(e)
            n_values_2.append(n)
            u_values_2.append(u)

            line_counter += 1
    return [xr_values_2, yr_values_2, xp_values_2, yp_values_2,e,n,u]



def enforce_epipolar_constraint(xryr,xpyp,FF):
    epipolar_vector_cost = []

    for i in range(len(xp_values)):
        M_1 = matrix_multiply(F, np.array(X_p[i]).reshape((3,1)))
        M_2 = matrix_multiply(np.array(X_r[i]).reshape(1,3),M_1)
        M = M_2[0][0]
        epipolar_vector_cost.append(M)
    return epipolar_vector_cost
    
         
def enforce_projection_constraint(PR, PP, er, nr, ur, ep, np, up ):
    num_rows = len(xr_pred)  # get the number of rows in the data

    xryr_pred = project_to_cam(PR, er, nr, ur, facing)
    xpyp_pred = project_to_cam(PP, ep, np, up, facing)
    # loop over the rows and compute the matrix for each row
    matrices = []
    for i in range(num_rows):
        args = np.array([xr[i] - xr_pred[i],
                         yr[i] - yr_pred[i],
                         xp[i] - xp_pred[i],
                         yp[i] - yp_pred[i]])
        matrix = np.reshape(args, (4, 1))
        matrices.append(matrix)

    # concatenate the matrices into a single 4n x 1 vector
    vector = np.concatenate(matrices, axis=0)

    return vector
    
def initialize_axes(yaw,pitch,roll):
    yaw_ref = yaw[0]
    yaw_proj = yaw[1]
    pitch_ref = pitch[0]
    pitch_proj = pitch[1]
    roll_ref = roll[0]
    roll_proj = roll[1]
    params = Parameters()
    params.add('yaw_ref', value=yaw[0], min=yaw[0]-0.1, max=yaw[0]+0.1)
    params.add('yaw_proj', value=yaw[1], min=yaw[1]-0.1, max=yaw[1]+0.1)
    params.add('pitch_ref', value=pitch[0], min=pitch[0]-0.1, max=pitch[0]+0.1)
    params.add('pitch_proj', value=pitch[1], min=pitch[1]-0.1, max=pitch[1]+0.1)
    params.add('roll_ref', value=roll[0], min=roll[0]-0.1, max=roll[0]+0.1)
    params.add('roll_proj', value=roll[1], min=roll[1]-0.1, max=roll[1]+0.1)
    return params

def optimization_func(xr_pred, yr_pred, xp_pred, yp_pred, xr, yr, xp, yp):
    #Enforcing the epipolar constraint
    costf = enforce_epipolar_constraint()
    #Enforcing the projection constraint
    costp = enforce_projection_constraint(xr_pred, yr_pred, xp_pred, yp_pred, xr, yr, xp, yp)

#Main part of the code
def main():
    
    ######Camera setup
    #read initial camera parameters and camera positions for two cameras (reference and pair)
    camcnt = 2
    configfilename = 'Patrick-Code/ConfigFile.xlsx'
    fx, fy, cx, cy, pitch, roll, yaw, camera_position = configfile(camcnt,configfilename)
    #find the distance between cameras and the direction cameras are facing 
    x_dist,y_dist,z_dist= GpstoENU(camera_position)
    facing = get_orientation(x_dist,y_dist)  
    #find the translation vector
    tx,ty,tz = enu_to_cam_coordinate(x_dist,y_dist,z_dist,facing)
    T = [tx,ty, tz, [1] ] 

    ######Input vectors
    #Reading xr,yr,xp,yp from FirstDummy.txt
    #reading_1('FirstDummy.txt')
    #Reading xr,yr,xp,yp,E,N,U from SecondDummy.txt
    #reading_2('SecondDummy1.txt')
    
    #xr_values = reading_1('FirstDummy.txt')[0]
    #yr_values = reading_1('FirstDummy.txt')[1]
    #xp_values = reading_1('FirstDummy.txt')[2]
    #yp_values = reading_1('FirstDummy.txt')[3]

    P_R=projection_matrix(fx[0],fy[0],cx[0],cy[0],yaw[0],pitch[0],roll[0],0,0,0)
    P_P=projection_matrix(fx[1],fy[1],cx[1],cy[1],yaw[1],pitch[1],roll[1],tx,ty,tz)
    FF = fundamental_matrix(P_R, P_P, T)
    print(FF)
    
    params = initialize_axes(yaw,pitch,roll)
    min_val = Minimizer(optimization_func, params, fcn_args=(x,), fcn_kws={'data': data})
    out = min_val.leastsq(Dfun=dfunc, col_deriv=1)
    fit2 = optimization_func(out.params, x)
    
    #Enforcing the epipolar constraint
    costf = enforce_fundamental_matrix()
    #Enforcing the projection constraint
    costp = enforce_projection_constraint(xr_pred, yr_pred, xp_pred, yp_pred, xr, yr, xp, yp)

#Execution of the code
main()