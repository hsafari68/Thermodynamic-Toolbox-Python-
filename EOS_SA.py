def PR78_SA(Tc, Pc, w, MW, T, P, BIC, comp):
    
    # Import necessary libraries and functions
    import numpy as np
    from EOS import PR78
    from Mixing_Rule import Mixing_Rule
    
    # Defining the constants
    delta1 = -0.4142 
    delta2 = 2.4142 
    C = 1
    
    # Defining the reduced properties
    Tr = T / Tc
    Pr = P / Pc
    
    # Calculating EOS parameters
    out_EOS_SA = PR78(Tc, Pc, w, MW, Tr, Pr)
    b_SA = out_EOS_SA['b']
    a_SA = out_EOS_SA['a']
    
    # Approximating the properties using mixing rule
    Out_MR_SA = Mixing_Rule(a_SA, b_SA, comp, T, P, BIC)
    Bm_SA = Out_MR_SA['Bm']
    bm_SA = Out_MR_SA['bm']
    Am_SA = Out_MR_SA['Am']
    am_SA = Out_MR_SA['am']
    Si_SA = Out_MR_SA['Si']
    Si_SA = Out_MR_SA['Si']
    
    # Constructing the Z factor coefficents
    coef_SA = [1, C * Bm_SA - 1, \
               Am_SA - Bm_SA * (1 + C) - Bm_SA ** 2 * (1 + 2 * C), \
                   C * (Bm_SA ** 2 + Bm_SA ** 3) - Am_SA * Bm_SA]
        
    # Solving Z factor cubic equaion using numpy
    Z_SA = np.roots(coef_SA)
    
    # Checking the roots
    if sum(np.isreal(Z_SA)) == 1:
        Z_SA = float(Z_SA[np.isreal(Z_SA)].real)
    else:
        if sum(Z_SA > 0) <= 2:
            Z_SA = max(Z_SA)
        else:
            Z_SA = min(Z_SA)
    
    # Updating fugacities and n_v
    phi_SA = np.exp(b_SA * (Z_SA - 1) / bm_SA - np.log(Z_SA - Bm_SA) \
                     - 1 / (delta2 - delta1) * Am_SA / Bm_SA * \
                         (Si_SA / am_SA * 2 - b_SA / bm_SA) * \
                             np.log((Z_SA + delta2 * Bm_SA) / \
                                    (Z_SA + delta1 * Bm_SA)))   
    f_SA = P * comp * phi_SA
    return f_SA