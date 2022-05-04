def get_harmonic_flag(l_array):
    ''' 
    Return the harmonic flag given harmonic moments. The monopole is always on.
    
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
    fx, fy, fz, dens, pot = args
    print(r'f_x = %.2f [km/s^2]' % fx)
    print(r'f_y = %.2f [km/s^2]' % fy)
    print(r'f_z = %.2f [km/s^2]' % fz)
    print(r'rho = %.2e [Msun/pc^3]' % dens)
    print(r'Phi = %8e [(km/s)^2]' % pot)
    return None