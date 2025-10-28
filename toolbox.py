# 2-layer QG Channel Model

from model_parameters import *
from scipy.sparse import csr_matrix
from scipy.sparse import diags
eps = 1e-4

def diff_x(Y_mat=Y-4,dx=d_lon*np.pi/180):
    # Central difference along the x direction (lon).
    # Cyclic boundary condition.
    A = np.zeros((Y_mat*X, Y_mat*X)) # 1-st order operating matrix.

    for i in range(Y_mat):
        for j in range(X):
            index = i*X + j
            jm1 = i*X + (j-1+X)%X
            jp1 = i*X + (j+1)%X
            # Fill the operating matrix.
            A[index,jm1] = -0.5
            A[index,jp1] = 0.5
    
    A = A/dx # Add resolution.
    return csr_matrix(A)

def diff_y(Y_mat=Y-4,dy=d_lat*np.pi/180):
    # Central difference along the y direction (lat).
    A = np.zeros((Y_mat*X, Y_mat*X)) # 1-st order operating matrix.

    for i in range(Y_mat):
        for j in range(X):
            index = i*X + j
            if (i==0) and (Y_mat==Y-4):    # Top boundary.
                ip1 = (i+1)*X+j
                A[index,ip1] = 0.5
            elif (i==Y_mat-1) and (Y_mat==Y-4):    # Bottom boundary.
                im1 = (i-1)*X+j
                A[index,im1] = -0.5
            elif 0<i<Y_mat-1:
                ip1 = (i+1)*X+j
                im1 = (i-1)*X+j
                A[index,ip1] = 0.5
                A[index,im1] = -0.5

    #A += eps*np.eye(Y*X) # Add perturbation.
    A = A/dy # Add resolution.
    return csr_matrix(A)

def diff_x2(Y_mat=Y-4,dx=d_lat*np.pi/180):
    # 2-nd order Central difference along the x direction (lon).
    A = np.zeros((Y_mat*X, Y_mat*X)) # 1-st order operating matrix.

    for i in range(Y_mat):
        for j in range(X):
            index = i*X + j
            jm1 = i*X + (j-1+X)%X
            jp1 = i*X + (j+1)%X
            # Fill the operating matrix.
            A[index,jm1] = 1.0
            A[index,jp1] = 1.0
            A[index,index] = -2.0
    
    #A += eps*np.eye(Y*X) # Add perturbation.
    A = A/(dx**2) # Add resolution.
    return csr_matrix(A)

def diff_y2(Y_mat=Y-4,dy=d_lat*np.pi/180):
    # 2nd order Central difference along the y direction (lat).
    A = np.zeros((Y_mat*X, Y_mat*X)) # 1-st order operating matrix.

    for i in range(Y_mat):
        for j in range(X):
            index = i*X + j
            if (i==0) and (Y_mat==Y-4):
                A[index,(i+1)*X+j] = 1.0
                A[index,index] = -2.0
            elif (i==Y_mat-1) and (Y_mat==Y-4):
                A[index,(i-1)*X+j] = 1.0
                A[index,index] = -2.0
            elif 0<i<Y_mat-1:
                A[index,(i-1)*X+j] = 1.0
                A[index,(i+1)*X+j] = 1.0
                A[index,index] = -2.0

    #A += eps*np.eye(Y*X) # Add perturbation.
    A = A/dy**2 # Add resolution.
    return csr_matrix(A)

def diff_xy(Y_mat=Y-4,dx=d_lat*np.pi/180,dy=d_lat*np.pi/180):
    # 1st central difference along the x and y direction.
    A = np.zeros((Y_mat*X, Y_mat*X)) # 1-st order operating matrix.

    for i in range(Y_mat):
        for j in range(X):
            index = i*X + j
            ip1_jp1 = (i+1)*X+(j+1)%X
            im1_jm1 = (i-1)*X+(j-1+X)%X
            ip1_jm1 = (i+1)*X+(j-1+X)%X
            im1_jp1 = (i-1)*X+(j+1)%X
            if (i==0) and (Y_mat==Y-4):
                A[index,ip1_jp1] = 1.0
                A[index,ip1_jm1] = -1.0
            elif (i==Y_mat-1) and (Y_mat==Y-4):
                A[index,im1_jm1] = 1.0
                A[index,im1_jp1] = -1.0
            elif 0<i<Y_mat-1:
                A[index,im1_jm1] = 1.0
                A[index,im1_jp1] = -1.0
                A[index,ip1_jp1] = 1.0
                A[index,ip1_jm1] = -1.0

    #A += eps*np.eye(Y*X) # Add perturbation.
    A = A/(4*dx*dy) # Add resolution.
    return csr_matrix(A)

def diff_x3(Y_mat=Y-4,dx=d_lat*np.pi/180):
    # 3-rd order Central difference along the x direction (lon).
    A = np.zeros((Y_mat*X, Y_mat*X)) # 1-st order operating matrix.

    for i in range(Y_mat):
        for j in range(X):
            index = i*X + j
            jm2 = i*X + (j-2+X)%X
            jm1 = i*X + (j-1+X)%X
            jp1 = i*X + (j+1)%X
            jp2 = i*X + (j+2)%X
            # Fill the operating matrix.
            A[index,jm2] = -1.0
            A[index,jm1] = 2.0
            A[index,jp1] = -2.0
            A[index,jp2] = 1.0
    
    #A += eps*np.eye(Y*X) # Add perturbation.
    A = A/(2*dx**3) # Add resolution.
    return csr_matrix(A)

def diff_y3(Y_mat=Y-4,dy=d_lat*np.pi/180):
    # 2nd order Central difference along the y direction (lat).
    A = np.zeros((Y_mat*X, Y_mat*X)) # 1-st order operating matrix.

    for i in range(Y_mat):
        for j in range(X):
            index = i*X + j
            im2 = (i-2)*X + j
            im1 = (i-1)*X + j
            ip1 = (i+1)*X + j
            ip2 = (i+2)*X + j

            if (i==0) and (Y_mat==Y-4):
                A[index,ip1] = -2.0
                A[index,ip2] = 1.0
            elif (i==1) and (Y_mat==Y-4):
                A[index,im1] = 2.0
                A[index,ip1] = -2.0
                A[index,ip2] = 1.0
            elif (i==Y_mat-2) and (Y_mat==Y-4):
                A[index,im2] = -1.0
                A[index,im1] = 2.0
                A[index,ip1] = -2.0
            elif (i==Y_mat-1) and (Y_mat==Y-4):
                A[index,im2] = -1.0
                A[index,im1] = 2.0
            elif 1<i<Y_mat-2:
                A[index,im2] = -1.0
                A[index,im1] = 2.0
                A[index,ip1] = -2.0
                A[index,ip2] = 1.0

    #A += eps*np.eye(Y*X) # Add perturbation.
    A = A/(2*dy**3) # Add resolution.
    return csr_matrix(A)

def diff_xy2(Y_mat=Y-4,dx=d_lat*np.pi/180,dy=d_lat*np.pi/180):
    # 1st derivative along x and 2nd derivative along y.
    A = np.zeros((Y_mat*X, Y_mat*X)) # 1-st order operating matrix.

    for i in range(Y_mat):
        for j in range(X):
            index = i*X + j
            ip1_jp1 = (i+1)*X+(j+1)%X
            im1_jm1 = (i-1)*X+(j-1+X)%X
            ip1_jm1 = (i+1)*X+(j-1+X)%X
            im1_jp1 = (i-1)*X+(j+1)%X
            i_jp1 = i*X + (j+1)%X
            i_jm1 = i*X + (j-1+X)%X
            if (i==0) and (Y_mat==Y-4):
                A[index,i_jm1] = 2.0
                A[index,i_jp1] = -2.0
                A[index,ip1_jp1] = 1.0
                A[index,ip1_jm1] = -1.0
            elif (i==Y_mat-1) and (Y_mat==Y-4):
                A[index,im1_jm1] = -1.0
                A[index,im1_jp1] = 1.0
                A[index,i_jm1] = 2.0
                A[index,i_jp1] = -2.0
            elif 0<i<Y_mat-1:
                A[index,im1_jm1] = -1.0
                A[index,im1_jp1] = 1.0
                A[index,i_jm1] = 2.0
                A[index,i_jp1] = -2.0
                A[index,ip1_jp1] = 1.0
                A[index,ip1_jm1] = -1.0

    #A += eps*np.eye(Y*X) # Add perturbation.
    A = A/(2*dx*dy**2) # Add resolution.
    return csr_matrix(A)

def diff_x2y(Y_mat=Y-4,dx=d_lat*np.pi/180,dy=d_lat*np.pi/180):
    # 1st derivative along x and 2nd derivative along y.
    A = np.zeros((Y_mat*X, Y_mat*X)) # 1-st order operating matrix.

    for i in range(Y_mat):
        for j in range(X):
            index = i*X + j
            ip1_jp1 = (i+1)*X+(j+1)%X
            im1_jm1 = (i-1)*X+(j-1+X)%X
            ip1_jm1 = (i+1)*X+(j-1+X)%X
            im1_jp1 = (i-1)*X+(j+1)%X
            ip1_j = (i+1)*X+j
            im1_j = (i-1)*X+j
            if (i==0) and (Y_mat==Y-4):
                A[index,ip1_j] = -2.0
                A[index,ip1_jp1] = 1.0
                A[index,ip1_jm1] = 1.0
            elif (i==Y_mat-1) and (Y_mat==Y-4):
                A[index,im1_jm1] = -1.0
                A[index,im1_jp1] = -1.0
                A[index,im1_j] = 2.0
            elif 0<i<Y_mat-1:
                A[index,im1_jm1] = -1.0
                A[index,im1_jp1] = -1.0
                A[index,im1_j] = 2.0
                A[index,ip1_j] = -2.0
                A[index,ip1_jp1] = 1.0
                A[index,ip1_jm1] = 1.0

    #A += eps*np.eye(Y*X) # Add perturbation.
    A = A/(2*dx**2*dy) # Add resolution.
    return csr_matrix(A)

def Vort():
    # Get the Laplacian operator in the spherical coordinate.
    # Output: an operating matrix of size (N,N).
    lat_all = get_lat_all()
    term_1 = diags((1/np.cos(lat_all)**2)) @ diff_x2()
    term_2 = -diags((np.tan(lat_all))) @ diff_y()
    #term_3 = -diags((1/np.tan(lat_all))) @ diff_y()
    term_4 = diff_y2()
    mat = diags((1/Corio()))@(term_1+term_2+term_4)/r**2
    return mat

def Laplacian():
    # Laplacian operator.
    lat_all = get_lat_all()
    term_1 = diags((1/np.cos(lat_all)**2)) @ diff_x2()
    term_2 = -diags((np.tan(lat_all))) @ diff_y()
    #term_3 = -diags((1/np.tan(lat_all))) @ diff_y()
    term_4 = diff_y2()
    mat = (term_1+term_2+term_4)/r**2
    return mat

def Corio(Y_mat=Y-4):
    # Get the Coriolis parameter of all the lat.
    # Output: vector of size N.
    lat_all = get_lat_all(Y_mat=Y_mat)
    f = 2*ome*np.sin(lat_all)
    return f

def get_lat_all(lat=lat,Y_mat=Y-4):
    # Get the lat matrix of size (Y,X) in the domain.
    # Flatten it into a vector of size N.
    # Unit: radian.
    if Y_mat==Y-4:
        lat_rad = lat[2:Y-2]*np.pi/180
    else:
        lat_rad = lat*np.pi/180
    lat_mat = np.tile(lat_rad,(X,1)).T
    lat_vect = lat_mat.flatten()
    return lat_vect

def reconstruct(vect,Y_mat=Y-4):
    # Convert the background flow to the correct size.
    if len(vect)!=Y*X:
        print("Wrong input in reconstruction!")
    mat = vect.reshape((Y,X))
    new_vect = mat[2:Y-2].reshape(-1)
    return new_vect

def Vort_y(Y_mat=Y-4):
    lat_all = get_lat_all(Y_mat=Y_mat)
    term_1 = diags(1/np.cos(lat_all)**2) @ diff_x2y(Y_mat=Y_mat)
    term_2 = - diags(np.tan(lat_all)) @ diff_y2(Y_mat=Y_mat)
    #term_3 = - diags(1/np.tan(lat_all)) @ diff_y2(Y_mat=Y_mat)
    term_4 = diff_y3(Y_mat=Y_mat)
    mat = diags(1/Corio(Y_mat=Y_mat)) @ ((term_1+term_2+term_4)/r**3)
    return mat

def Lap_y(Y_mat=Y-4):
    lat_all = get_lat_all(Y_mat=Y_mat)
    term_1 = diags(1/np.cos(lat_all)**2) @ diff_x2y(Y_mat=Y_mat)
    term_2 = - diags(np.tan(lat_all)) @ diff_y2(Y_mat=Y_mat)
    #term_3 = - diags(1/np.tan(lat_all)) @ diff_y2(Y_mat=Y_mat)
    term_4 = diff_y3(Y_mat=Y_mat)
    mat = (term_1+term_2+term_4)/r**3
    return mat

def Vort_x(Y_mat=Y-4):
    lat_all = get_lat_all(Y_mat=Y_mat)
    term_1 = diags((1/np.cos(lat_all)**2)) @ diff_x3(Y_mat=Y_mat)
    term_2 = -diags((np.tan(lat_all))) @ diff_xy(Y_mat=Y_mat)
    #term_3 = -diags((1/np.tan(lat_all))) @ diff_xy(Y_mat=Y_mat)
    term_4 = diff_xy2(Y_mat=Y_mat)
    mat = diags((1/(Corio(Y_mat=Y_mat)*np.cos(lat_all))))\
          @ ((term_1+term_2+term_4)/r**3)
    return mat

def Lap_x(Y_mat=Y-4):
    lat_all = get_lat_all(Y_mat=Y_mat)
    term_1 = diags((1/np.cos(lat_all)**2)) @ diff_x3(Y_mat=Y_mat)
    term_2 = -diags((np.tan(lat_all))) @ diff_xy(Y_mat=Y_mat)
    #term_3 = -diags((1/np.tan(lat_all))) @ diff_xy(Y_mat=Y_mat)
    term_4 = diff_xy2(Y_mat=Y_mat)
    mat = diags((1/np.cos(lat_all))) @ ((term_1+term_2+term_3+term_4)/r**3)
    return mat

def Vort_y_new(Y_mat=Y-4):
    lat_all = get_lat_all(Y_mat=Y_mat)
    term_1 = diags((1/np.cos(lat_all)**2)) @ diff_x2y(Y_mat=Y_mat)
    term_2 = -diags((np.tan(lat_all))) @ diff_y2(Y_mat=Y_mat)
    term_3 = -diags((1/np.tan(lat_all))) @ diff_y2(Y_mat=Y_mat)
    term_4 = diff_y3(Y_mat=Y_mat)
    old = diags((1/Corio(Y_mat=Y_mat))) @ ((term_1+term_2+term_3+term_4)/r**3)

    term_5 = diags((2/np.cos(lat_all)**3-\
            np.sin(lat_all)**(-2)*np.cos(lat_all)**(-1))) @ diff_x2(Y_mat=Y_mat)
    term_6 = -diags((np.tan(lat_all)/np.cos(lat_all))) @ diff_y(Y_mat=Y_mat)
    term_7 = diags((2*np.cos(lat_all)**2/np.sin(lat_all)**3+\
                     1/np.sin(lat_all))) @ diff_y(Y_mat=Y_mat)
    new = (term_5+term_6+term_7)/(2*ome*r**3)
    
    return old+new