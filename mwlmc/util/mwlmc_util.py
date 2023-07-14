import numpy as np

__all__ = ['get_harmonic_flag', 'print_field_info', 'make_random_obs', 'radecpms_to_xyzv']

def get_harmonic_flag(l_array):
    ''' 
    Returns the harmonic flag given harmonic moments. The monopole is always on.
    
    Parameters
    ----------
    l_array : array-like integers 
       Selected spherical harmonics. 
    
    Returns
    -------
    flag : int
       Flag to put into mwlmc field functions.
    '''
    flag = int(sum([2**(l-1) for l in l_array]))
    return flag
  
def print_field_info(*args):
    ''' 
    Prints information on the fields.
    
    Parameters
    ----------
    fx, fy, fz : float
       Forces in x, y and z-direction.
    rho : float
       Density 
    Phi : float
       Potential
    '''    
    fx, fy, fz, dens, pot = args
    print(r'f_x = %.2f [km/s^2]' % fx)
    print(r'f_y = %.2f [km/s^2]' % fy)
    print(r'f_z = %.2f [km/s^2]' % fz)
    print(r'rho = %.2e [Msun/pc^3]' % dens)
    print(r'Phi = %8e [(km/s)^2]' % pot)
    return None
  
def make_random_obs(*obs, N=10, seed=265):
    ''' 
    Returns random draws from observations.
    
    Parameters
    ----------
    ra, dec : float 
       Ra and dec in deg.
    dhel : float
       Distance in kpc.
    pmra, pmdec : float
       Proper motions in mas/yr.
    vhel : float
       Heliocentric radial velocity in km/s.
    edhel : float
       Error of distance in kpc.
    epmra, epmdec : float
       Error of proper motions in mas/yr.
    evhel : float
       Error of radial velocity in km/s.
    N : int (default = 10)
       Number of random realisations of observables.
    seed : int (default = 265)
       Seed for the rng.
    
    Returns
    -------
    ra, dec, dhel, pmra, pmdec, vhel : array-like float
       
    '''        
    rng = np.random.default_rng(seed)
    ra, dec, dhel, pmra, pmdec, vhel, edhel, epmra, epmdec, evhel = obs
    ra = np.full(N, ra)
    dec = np.full(N, dec)
    dhel = rng.normal(dhel, scale=edhel, size=N)
    pmra = rng.normal(pmra, scale=epmra, size=N)
    pmdec = rng.normal(pmdec, scale=epmdec, size=N)
    vhel = rng.normal(vhel, scale=evhel, size=N)
    return ra, dec, dhel, pmra, pmdec, vhel
  
###### Coordinate transformations
k = 4.74047

AG = np.array([[-0.0548755604, 0.4941094279, -0.8676661490], 
               [-0.8734370902, -0.4448296300, -0.1980763734],
               [-0.4838350155, 0.7469822445, 0.4559837762]])

R_phirad = np.array([[-0.4776303088, -0.1738432154, 0.8611897727],
                     [0.510844589, -0.8524449229, 0.111245042],
                     [0.7147776536, 0.4930681392, 0.4959603976]])

def _M_UVW_pm(phi1, phi2):
    if not isinstance(phi1, np.ndarray):
        M = np.array([[np.cos(phi1)*np.cos(phi2), -np.sin(phi1),
                       -np.cos(phi1)*np.sin(phi2)],
                      [np.sin(phi1)*np.cos(phi2), np.cos(phi1), 
                       -np.sin(phi1)*np.sin(phi2)],
                      [np.sin(phi2), 0., np.cos(phi2)]])
        
    else:
        reshaping = (len(phi1), 3, 3)
        M = np.array([])
        for i in range(len(phi1)):
            M = np.append(M, np.array([[np.cos(phi1[i]) * np.cos(phi2[i]), 
                                        - np.sin(phi1[i]), - np.cos(phi1[i]) 
                                        * np.sin(phi2[i])],
                          [np.sin(phi1[i]) * np.cos(phi2[i]), 
                           np.cos(phi1[i]), - np.sin(phi1[i]) * 
                           np.sin(phi2[i])],
                          [np.sin(phi2[i]), 0., np.cos(phi2[i])]]))

        M = M.reshape(reshaping)
    return M
  
  
def _pmrapmdec_to_UVW(ra, dec, d, pmra, pmdec, vr):
    input_vec = np.array([vr, k*d*pmra, k*d*pmdec]).T
    ra, dec = ra*np.pi/180., dec*np.pi/180.
    M_mat = _M_UVW_pm(ra, dec)    
    
    if not isinstance(ra, np.ndarray):
        res = np.matmul(np.matmul(AG.transpose(), M_mat), input_vec)
    else:
        res = np.zeros((len(M_mat), 3))

        for ii in range(len(M_mat)):
            res[ii] = np.matmul(np.matmul(AG.transpose(), 
                                          M_mat[ii]), input_vec[ii])
            
        res = res.T    
    U, V, W = res
    return U, V, W
  
def _radec_to_lb(ra, dec, d):
    input_vec = np.array([d * np.cos(ra*np.pi/180.) * np.cos(dec*np.pi/180.), 
                          d * np.sin(ra*np.pi/180.) * np.cos(dec*np.pi/180.), 
                          d * np.sin(dec*np.pi/180.)])
    res = np.matmul(AG.transpose(), input_vec)
    
    res0 = res[0] # = d * cos(l) * cos(b)
    res1 = res[1] # = d * sin(l) * cos(b)
    res2 = res[2] # = d * sin(b)
    
    b = np.arcsin(res2 / d)
    cosl = res0 / d / np.cos(b)
    sinl = res1 / d / np.cos(b)
    l = np.arctan2(sinl, cosl)

    l = l * 180. / np.pi
    b = b * 180. / np.pi
    return l, b, d  

def _lbd_to_xyz(l, b, d):
    try: 
        X, Y, Z = np.array([d * np.cos(b*np.pi/180.) * np.cos(l*np.pi/180.), 
                            d * np.cos(b*np.pi/180.) * np.sin(l*np.pi/180.), 
                            d * np.sin(b*np.pi/180.)]).T
    except:
        X, Y, Z = np.array([d * np.cos(b*np.pi/180.) * np.cos(l*np.pi/180.), 
                            d * np.cos(b*np.pi/180.) * np.sin(l*np.pi/180.), 
                            d * np.sin(b*np.pi/180.)])
    return X, Y, Z
  
def _gal_to_galcenrect(X, Y, Z, U, V, W, 
                      x_sun=-8.249, z_sun=0., 
                      v_sun=[11.1, 245, 7.3]):
    d_GC = np.linalg.norm([x_sun, z_sun])
    costheta, sintheta = x_sun/d_GC, z_sun/d_GC
    xg_out = np.dot(np.array([[costheta, 0., -sintheta], [0.,1.,0.], 
                              [sintheta, 0., costheta]]), 
                    np.array([- X + d_GC, Y , np.sign(x_sun)*Z])).T    
    vg_out = np.dot(np.array([[costheta, 0., -sintheta], [0.,1.,0.], 
                              [sintheta, 0., costheta]]), 
                    np.array([- U, V, np.sign(x_sun) * W])).T + np.array(v_sun)
    return xg_out, vg_out
  
def radecpms_to_xyzv(ra, dec, d, pmra, pmdec, vr, 
                     x_sun=-8.249, z_sun=0., 
                     v_sun=[11.1, 245, 7.3]):
    ''' 
    Return cartesian coordinates from observational data
    
    Parameters
    ----------
    ra, dec : (array-like) float 
       Ra and dec in deg
    d : (array-like) float
       Distance in kpc
    pmra, pmdec : (array-like) float
       Proper motions in mas/yr
    vhel : (array-like) float
       Heliocentric radial velocity in km/s
    x_sun, z_sun: float (optional, defaults=-8.249, 0.)
       Position of Sun
    v_sun : array-like floats (optinal, default=[11.1, 245, 7.3])
    
    Returns
    -------
    x, y, z : (array-like) float
       Position in kpc
    vx, vy, vz : (array-like) float
       Velocity in km/s
    '''  
    U, V, W = _pmrapmdec_to_UVW(ra, dec, d, pmra, pmdec, vr)
    l, b, d = _radec_to_lb(ra, dec, d)
    X, Y, Z = _lbd_to_xyz(l, b, d)
    xg_out, vg_out = _gal_to_galcenrect(X, Y, Z, U, V, W, 
                                        x_sun=x_sun, z_sun=z_sun, 
                                        v_sun=v_sun)
    if not isinstance(ra, np.ndarray):
        x, y, z = xg_out
        vx, vy, vz = vg_out
    else:
        x, y, z = xg_out[:,0], xg_out[:,1], xg_out[:,2]
        vx, vy, vz = vg_out[:,0], vg_out[:,1], vg_out[:,2]    
    return x, y, z, vx, vy, vz   