import numpy as np

def errorFunction(z, a):
    """
    Compute exp(a) * erf(z) using various approximations.
    
    Parameters:
    z : array_like
        Input values (can be complex).
    a : float
        Exponent value.
        
    Returns:
    erf_val : ndarray
        The computed error function values.
    """
    z = np.asarray(z, dtype=complex)
    erf_val = np.sign(np.real(z)) * np.exp(a) - np.exp(a - z**2) * z / (np.sqrt(np.pi) * (z**2 + 0.5 / (1 + 1 / (1.5 + z**2))))
    erf_val = np.asarray(erf_val)
    # Patch Poles with Pade Approximation
    mask1 = np.abs(z + 2j) < 1.1
    if np.any(mask1):
        sub_z = z[mask1]
        erf_val[mask1] = np.exp(a) * (-18.564802j + 12.326185 * (2j + sub_z) + 13.9957261j * (2j + sub_z)**2 
                                     - 7.0912286 * (2j + sub_z)**3 - 2.1539997j * (2j + sub_z)**4 + 0.80057144 * (2j + sub_z)**5) \
                        / (1 - 2.6545518j * (2j + sub_z) - 2.9260194 * (2j + sub_z)**2 + 1.665267j * (2j + sub_z)**3 
                           + 0.48178231 * (2j + sub_z)**4 - 0.054052386j * (2j + sub_z)**5)

    mask2 = np.abs(z - 2j) < 1.1
    if np.any(mask2):
        sub_z = z[mask2]
        erf_val[mask2] = np.exp(a) * (18.564802j + 12.326185 * (-2j + sub_z) - 13.9957261j * (-2j + sub_z)**2 
                                     - 7.0912286 * (-2j + sub_z)**3 + 2.1539997j * (-2j + sub_z)**4 + 0.80057144 * (-2j + sub_z)**5) \
                        / (1 + 2.6545518j * (-2j + sub_z) - 2.9260194 * (-2j + sub_z)**2 - 1.665267j * (-2j + sub_z)**3 
                           + 0.48178231 * (-2j + sub_z)**4 + 0.054052386j * (-2j + sub_z)**5)

    mask3 = np.abs(z) < 1.4
    if np.any(mask3):
        sub_z = z[mask3]
        erf_val[mask3] = np.sign(np.real(sub_z)) * np.sqrt(np.exp(2 * a) - np.exp(2 * a - (1.2732395 + 0.14001229 * sub_z**2) 
                                                                 / (1 + 0.14001229 * sub_z**2) * sub_z**2))
    
    return erf_val
