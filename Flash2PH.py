def Flash2PH(CompName, z, T, P, BIC, c, Vp = 100 ):
    
    # Import necessary libraries and functions
    import numpy as np
    from EOS import PR78
    from K_value import K_Wilson
    from Stability_Analysis import Stability_Analysis
    from Mixing_Rule import Mixing_Rule
    from Volume_Derivative import Volume_Derivative
    from Select_Component import Select_Component
    
    # Inserting the component properties
    Props = Select_Component(CompName)
    
    # Converting the input data to array for saving time
    z = np.array(z)
    Tc = np.array(Props['Tc'])
    Pc = np.array(Props['Pc'])
    MW = np.array(Props['MW'])
    w =  np.array(Props['w'])
    BIC = np.array(BIC)
    c = np.array(c)
    PCH = np.array(Props['PCH'])
    
    # Defining the gas global constant
    R = 8.31447
    
    # Defining the constants
    delta1 = -0.4142 
    delta2 = 2.4142 
    C = 1
     
    # Normalizing the initial compositions    
    z = z / sum(z)

        
    # Defining the reduced properties
    Tr = T / Tc
    Pr = P / Pc
        
    # Calculating EOS parameters
    Out_EOS = PR78(Tc, Pc, w, MW, Tr, Pr) 
    b = Out_EOS['b']
    a = Out_EOS['a']
        
    # Initialization of n_v and K_val of Withson
    n_v = 0.5
    K_est = K_Wilson(w, Tr, Pr)

    # Initialization of liquid and vapor compositions
    K_l = K_v = K_est
    x = z / (1 + n_v * (K_l - 1))
    y = z / (1 + (1 - n_v) * (1 / K_v - 1))
    
    # Itrating the x and y
    Tol = 1
    Itr = 0
    while Tol > 1e-12:
        Itr += 1

        # Stability check for n_v>1 or n_v<0 
        if n_v > 1.0:
            Out_SA = Stability_Analysis(Tc, Pc, MW, w, BIC, z, T, P)
            K_est = Out_SA['Ki'][0]
        elif n_v < 0:
            Out_SA = Stability_Analysis(Tc, Pc, MW, w, BIC, z, T, P)
            K_est = Out_SA['Ki'][1]
        # Using mixing rule to calculate am, bm, Am and Bm of phases
        Out_MR_l = Mixing_Rule(a, b, x, T, P, BIC)
        Out_MR_v = Mixing_Rule(a, b, y, T, P, BIC)
        Bm_l = Out_MR_l['Bm']
        bm_l = Out_MR_l['bm']
        Am_l = Out_MR_l['Am']
        am_l = Out_MR_l['am']
        Si_l = np.array(Out_MR_l['Si'])
        Bm_v = Out_MR_v['Bm']
        bm_v = Out_MR_v['bm']      
        Am_v = Out_MR_v['Am']
        am_v = Out_MR_v['am']        
        Si_v = np.array(Out_MR_v['Si'])
        
        # Constructing the Z factor coefficents
        coef_l = [1, \
                  C * Bm_l - 1, \
                  Am_l - Bm_l * (1 + C) - Bm_l ** 2 * (1 + 2 * C), \
                      C * (Bm_l ** 2 + Bm_l ** 3) - Am_l * Bm_l]
        coef_v = [1, \
                  C * Bm_v - 1, \
                  Am_v - Bm_v * (1 + C) - Bm_v ** 2 * (1 + 2 * C), \
                      C * (Bm_v ** 2 + Bm_v ** 3) - Am_v * Bm_v]
    
        # Solving Z factor cubic equaion using numpy
        Z_l = np.roots(coef_l)
        Z_v = np.roots(coef_v)
        
        # Checking the roots
        if sum(np.isreal(Z_l)) == 1:
            Z_l = float(Z_l[np.isreal(Z_l)].real)
        else:
            if sum(Z_l > 0) <= 2:
                Z_l = max(Z_l)
            else:
                Z_l = min(Z_l)
        if sum(np.isreal(Z_v)) == 1:
            Z_v = float(Z_v[np.isreal(Z_v)].real)
        else:
            Z_v = max(Z_v)
        
        # Updating fugacities and n_v
        phi_l = np.exp(b * (Z_l - 1) / bm_l - np.log(Z_l - Bm_l) \
                         - 1 / (delta2 - delta1) * Am_l / Bm_l * \
                             (Si_l / am_l * 2 - b / bm_l) * \
                                 np.log((Z_l + delta2 * Bm_l) / \
                                        (Z_l + delta1 * Bm_l)))
        phi_v = np.exp(b * (Z_v - 1) / bm_v - np.log(Z_v - Bm_v) \
                         - 1 / (delta2 - delta1) * Am_v / Bm_v * \
                             (Si_v / am_v * 2 - b / bm_v) * \
                                 np.log((Z_v + delta2 * Bm_v) / \
                                          (Z_v + delta1 * Bm_v)))
        sum1 = sum(z * (K_est - 1) / (1 + n_v * (K_est - 1)))
        sum2 = -sum(z * (K_est - 1) ** 2 / (1 + n_v * (K_est - 1)) ** 2)
        n_v -= sum1 / sum2
        n_v_min = 1 / (1 - max(K_est))
        n_v_max = 1 / (1 - min(K_est))
        n_v = min(max(n_v, n_v_min), n_v_max)
        
        # Updating K_est
        f_l = P * x * phi_l
        f_v = P * y * phi_v
        if (Itr + 2) % 6 == 0:
            f_l_2 = f_l
            f_v_2 = f_v
        elif (Itr + 1) % 6 == 0:
            f_l_1 = f_l
            f_v_1 = f_v    
        if Itr % 6 == 0:
            b02 = sum(np.log(f_l / f_v) * np.log(f_l_2 / f_v_2))
            b01 = sum(np.log(f_l /f_v) * np.log(f_l_1 / f_v_1))
            b12 = sum(np.log(f_l_1 / f_v_1) * np.log(f_l_2 / f_v_2))
            b11 = sum(np.log(f_l_1 / f_v_1) * np.log(f_l_1 / f_v_1))
            b22 = sum(np.log(f_l_2 / f_v_2) * np.log(f_l_2 / f_v_2))
            mu_1 = (b02 * b12 - b01 * b22) / (b11 * b22 - b12 * b12)
            mu_2 = (b01 * b12 - b02 * b11) / (b11 * b22 - b12 * b12)
            K_est = np.exp(np.log(K_est) + \
                           (np.log(f_l / f_v) - mu_2 * \
                            np.log(f_l_1 / f_v_1)) / (1 + mu_1 + mu_2))      
        else:            
            K_est = K_est * f_l / f_v 
        
        # Updating liquid and vapor compositions
        x = z / (1 + n_v * (K_est - 1))
        y = z / (1 + (1 - n_v) * (1/K_est - 1))
        
        # Updating Stopping criteria tolerance
        Tol = sum((1 - f_l / f_v) ** 2)
    
    
    # Normalizing the compositions
    x = x / sum(x)
    y = y / sum(y)

    # Calculating liquid and vapor volumes
    V_l_EOS = Z_l * R * T / P
    V_v_EOS = Z_v * R * T / P
    
    # Applying volume shifts to volumes
    V_l_SH = V_l_EOS - sum(x * c)
    V_v_SH = V_v_EOS - sum(y * c)
    V = [V_l_SH, V_v_SH]
    
    # Calculating liquid and vapor molar densities
    DenM = [1 / V_l_SH, 1 / V_v_SH]
    
    # Calculating liquid and vapor densities
    Rho_l = DenM[0] * sum(x * MW)
    Rho_v = DenM[1] * sum(y * MW)
    
    # Calculating IFT using Parachor
    MW_l = sum(x * MW)
    MW_v = sum(y * MW)
    IFT_PCH = (sum(((x * Rho_l / MW_l * 1E-6) - (y * Rho_v / MW_v * 1E-6)) * PCH)) ** 4
    # Evaluating volume and pressure derivatives
    Out_VD = Volume_Derivative(V, a, am_l, Am_l, am_v, Am_v, b, bm_l, Bm_l, \
                               bm_v, Bm_v, Si_l, Si_v, Z_l, Z_v, P, T, BIC, \
                                   x, y, Vp)
    DVP = Out_VD['DVP']
    Vti = Out_VD['Vti']

    # Constructing and returning output dictionary   
    Flash2PH_Out = {'n_v':n_v, 'x':list(x), 'y':list(y), 'DenM':DenM, \
                    'Rho_l':Rho_l, 'Rho_v':Rho_v, 'IFT_PCH':IFT_PCH, \
                        'DVP':DVP, 'Vti':Vti}
    return Flash2PH_Out