""" In this module, we create functions that can be then used to create an animation
of the 2D shallow water equations.
"""
import numpy as np

def check_stability(dx, dy, dt, H, xdomain, ydomain, init_eta):
    """Checks whether the simulation will be stable with the given parameters and returns
    True if it will be stable and False otherwise. It also returns the maximum value for
    dt which is <checkval>
    
    Parameters:
        dx: step size x-direction
        dy: step size y-direction
        dt: time step
        H: Bathymetry function (should take in the y value first and then the x value)
        xdomain: A list of 2 elements, the first being the starting point of the
                x-interval and the second being the end point of the interval. E.g. the
                domain for the x values is [0,20] where 0m is the starting point and 20m is
                the end point.
        ydomain: similar to <xdomain>, but now for the y direction interval
        init_eta: array of length <len(xdomain[1]/dx)+2>x<len(ydomain[1]/dy)+2> containing
                  the values of eta at time t=0.
    """
    g = 9.8 # acceleration due to gravity in meters/second
    xvals = np.arange(xdomain[0], xdomain[1])
    yvals = np.arange(ydomain[0], ydomain[1])
    X,Y = np.meshgrid(xvals, yvals)
    hmax = np.max(H(Y,X) + init_eta)
    dxdymin = min(dx,dy)
    checkval = dxdymin/((2*g*hmax)**(1/2))
    if dt <= checkval:
        boolean = True
    else:
        boolean = False
    return boolean, checkval

def plus_and_minus(arr,n,j,k):
    """Calculates the terms needed in the function <uvh>. <val_p> is the u+ or v+ value
    from the lab manual and <val_m> is the u- or v- value from the lab manual. These
    are needed so that we can choose the correct indices to evaluate (eta+H) when we're
    calculating u*(eta+H) or v*(eta+H), correct meaning with respect to the upstream
    scheme. Also,
    
        if arr[n+1,j,k] > 0:
            val_p == arr[n+1,j,k] and val_m == 0
            (i.e. it "chooses" eta[n,j,k]+H(at x,y corresponding to j,k) in <uvh>)
        elif  arr[n+1,j,k] < 0:
            val_p == 0 and val_m == arr[n+1,j,k]
            (i.e. it "chooses" eta[n,j,k+1]+H(...) for u or eta[n,j+1,k]+H(...) for v)
        else:
            i.e. if arr[n+1,j,k] == 0
            val_p = val_m = 0
    
    Parameters:
        arr: This is either the u-array or the v-array.
        n: This is the time index we use for the eta+H terms and <n+1> is the index
            we use for the u and v terms.
        j: y-axis index
        k: x-axis index
    """
    val_p = (1/2)*(arr[n+1,j,k] + abs(arr[n+1,j,k])) # the u+ or v+ value
    val_m = (1/2)*(arr[n+1,j,k] - abs(arr[n+1,j,k])) # the u- or v- value
    return val_p, val_m



def uvh(n,j,k,index,u,v,eta,H,dx,dy):
    """Calculates the u*(eta+H) and v*(eta+H) values that show up in the calculation
    of the eta-star values. The value that is returned then gets used to calculate
    the derivative of u*(eta+H) and v*(eta+H) which in turn gets used in calculating
    eta-star values. 
    
    Note that when evaluating H at the (x,y) values corresponding to the j,k indices, 
    we use y = (j-1)*dy and x = (k-1)*dx instead of j and k since we chose our array 
    lengths in such a way that we can calcualte the eta values more conveniently but 
    also so that we get that the eta values in the array properly correspond to the 
    spatial domain.
    
    Parameters:
        n: time index for the eta star value. This ends up being the time index we use
            for the eta+H terms and <n+1> is the index we use for the u and v terms.
        j: y-axis index
        k: x-axis index
        index: Tells us whether we should calculate the u*(eta+H) or the v*(eta+H) value.
            <index>=1 tells us to calculate the u-array value and <index>=2 for the v-array
            value.
    """
    # if asking for u array:
    if index == 1:
        u_p, u_m = plus_and_minus(u,n,j,k)
        s = u_p*(eta[n,j,k] + H(dy*(j-1),dx*(k-1))) + u_m*(eta[n,j,k+1]+H(dy*(j-1),dx*k))
            
    # if asking for v array:        
    elif index == 2:
        v_p, v_m = plus_and_minus(v,n,j,k)
        s = v_p*(eta[n,j,k] + H(dy*(j-1),dx*(k-1))) + v_m*(eta[n,j+1,k]+H(dy*j,dx*(k-1)))
            
    return s



def eta_star(n,j,k,u,v,eta,H,dt,dx,dy):
    """Calculates eta-star. For the points in the interior, the function uses the
    back-difference derivative along with the upstream scheme to calculate the
    d(uh)/dx and d(vh)/dx terms. For the cases at the edges, the function uses
    the back-difference derivative as before, but instead of using the upstream 
    scheme to determine which indices to evaluate (eta+H) at, it simply uses the
    value of (eta+H) at [n,j,k].
    
    Parameters:
        n: time index for the eta star value. This ends up being the time index we use
            for the eta+H terms and <n+1> is the index we use for the u and v terms.
        j: y-axis index
        k: x-axis index
    """
    M, N = len(u[0]), len(u[0,0])
    # for d(uh)
    if (k != N-1 and k != 1): # If not at the left or right edge
        u_he = uvh(n,j,k,1,u,v,eta,H,dx,dy)
        u_hw = uvh(n,j,k-1,1,u,v,eta,H,dx,dy)
        duh = u_he - u_hw
    else:
        duh = (u[n+1,j,k] - u[n+1,j,k-1])*(eta[n,j,k] + H(j,k))
        
    # for d(vh)
    if (j != M-1 and j != 1): # If not at the top or bottom edge
        v_hn = uvh(n,j,k,2,u,v,eta,H,dx,dy)
        v_hs = uvh(n,j-1,k,2,u,v,eta,H,dx,dy)
        dvh = v_hn - v_hs
    else:
        dvh = (v[n+1,j,k] - v[n+1,j-1,k])*(eta[n,j,k] + H(j,k))
        
    eta_star_val = eta[n,j,k] - dt*((duh/dx) + (dvh/dy))
        
    return eta_star_val



def shallow_water_simulation(Lx, Ly, T, eps, dx, dy, dt, init_cond, H):
    """This function returns arrays with values that can be used to create a 2d 
    shallow water simulation. The simulation is only for a rectangular "tub" with
    dimensions <Lx>x<Ly>, where the edges of the "tub" are considered to be rigid
    walls (i.e. we will take the x-component of the velocities to be 0 at the left
    and right edges and the y-component of the velocities to be 0 at the top and
    bottom edges. Also, 
    
    This function returns 3 KxMxN arrays:
            u: array containing the x-components of the velocity
            v: array containing the y-components of the velocity
            eta: array containing the height of the water above the undisturbed 
                 water level
    K is the time dimension, M is the y dimension, and N is the x dimension. Also,
    note that if you want to plot any of the final arrays with the x values horizonatally
    and the y values vertically, you have to tranpose the spatial dimensions of the 
    array. E.g. for eta, you would need to plot eta[time].transpose().
    
    Note that eta[:,1:,1:] contains the values while eta[:,0,:] and eta[:,:,0]
    are set to None. eta is padded this way to make the calculations easier. Also,
    for <u> and <v>, since the eta array is created in such a way that we have eta
    values at the edge of the "tub", the values u[:,:,0] and v[:,0,:] don't
    correspond to any physical location inside or at the edges of the tub.
        
    This function also uses the outside functions <plus_and_minus>, <uvh>, and
    <eta_star>.
    
    Note: H (the bathymetry function) should be defined before calling this function.
    Also, it should take in the y value first and then the x value. I.e. H(y,x).
    
    Parameters:
        Lx: x-dimension of the tub
        Ly: y-dimension of the tub
        T: total duration of the simulation in seconds
        eps: smoothing parameter (should be small, i.e. much less than 1, and positive)
        dx: step size for x-direction
        dy: step size for y-direction
        dt: time step
        init_cond: array of 3 arrays, each with the same dimensions as the values that
                   are returned by thefunction, containing the values of u,v,η at time
                   t=0.
    """
    g = 9.8 # acceleration due to gravity in meters/second
    
    # Dimensions for the arrays:
    N = int(Lx/dx) + 2 # number of grid points, x-direction
    M = int(Ly/dy) + 2 # number of grid points, y-direction 
    K = int(T/dt) + 1 # number of grid points, in time
    
    # Intial arrays
    eta = np.zeros((K, M, N)) # holds the values of the height above sea level
    u = np.zeros((K, M, N)) # holds the x components of horizontal velocity
    v = np.zeros((K, M, N)) # holds the y components of horizontal velocity
    etas = np.zeros((M,N)) # holds the eta star values for the current iteration of the loop
    
    etas[0,:], etas[:,0] = None, None # sets the indices at which we don't have eta star as None
    eta[:,0,:], eta[:,:,0] = None, None # sets the indices at which we don't have eta as None
    u[:,0,:] = None # sets the indices at which we don't have u as None
    v[:,:,0] = None # sets the indices at which we don't have v as None
    
    # Plugging in the initial conditions:
    u[0,1:,:] = init_cond[0][1:,:]
    v[0,:,1:] = init_cond[1][:,1:]
    eta[0,1:,1:] = init_cond[2][1:,1:]
    
    # === For loop ===
    for n in range(0, K-1):
        # Calculating the u[n+1,j,k] values
        # Not including the left and right edges of the array since u is set to 0 there
        for j in range(1, M):
            for k in range(1, N-1):
                u[n+1,j,k] = u[n,j,k] - (dt/dx)*g*(eta[n,j,k+1] - eta[n,j,k])
        
        # Calculating the v[n+1,j,k] values
        # Not including the top and bottom edges of the array since v is set to 0 there
        for k in range(1, N):
            for j in range(1, M-1):
                v[n+1,j,k] = v[n,j,k] - (dt/dy)*g*(eta[n,j+1,k] - eta[n,j,k])
                
        # Calculating the eta star values:
        for j in range(1,M):
            for k in range(1, N):
                etas[j,k] = eta_star(n,j,k,u,v,eta,H,dt,dx,dy)
                
        # Calculating the eta[n+1,j,k] values
        # First, for the interior points (i.e. not on the edges):
        for j in range(2, M-1):
            for k in range(2, N-1):
                eta[n+1,j,k] = (1-eps)*etas[j,k] + \
                (eps/4)*(etas[j+1,k]+etas[j-1,k]+etas[j,k+1]+etas[j,k-1])
        
        # For the top and bottom edges, but not at the corners:
        for j in [1, M-1]:
            # Want to choose l=1 or l=-1 so that we have either etas[j+1,k] or etas[j-1,k]
            # but not both when calculating eta[n+1,j,k] since we dont have the j+1 value for
            # eta at j = M-1 and the j-1 value for eta at j=1:
            if j == 1:
                l = 1
            elif j == M-1:
                l = -1
            for k in range(2, N-1):
                eta[n+1,j,k] = (1-eps)*etas[j,k] + (eps/3)*(etas[j+l,k]+etas[j,k+1]+etas[j,k-1])
                
        # For the right and left edges, but not at the corners:
        for k in [1, N-1]:
            # Want to choose l=1 or l=-1 so that we have either etas[j,k+1] or etas[j,k-1]
            # but not both when calculating eta[n+1,j,k] since we dont have the k+1 value for
            # eta at k = N-1 and the k-1 value for eta at k=1:
            if k == 1:
                l = 1
            elif k == N-1:
                l = -1
            for j in range(2, M-1):
                eta[n+1,j,k] = (1-eps)*etas[j,k] + (eps/3)*(etas[j+1,k]+etas[j-1,k]+etas[j,k+l])
                
        # At the corners
        for j in [1, M-1]:
            # l is the same as it was before for j=1 and j=M-1
            if j == 1:
                l = 1
            elif j == M-1:
                l = -1
            for k in [1,N-1]:
                # m serves the same purpose that l did before for k=1 and k=N-1
                if k == 1:
                    m = 1
                elif k == N-1:
                    m = -1
                eta[n+1,j,k] = (1-eps)*etas[j,k] + (eps/2)*(etas[j+l,k]+etas[j,k+m])
                
    return u, v, eta