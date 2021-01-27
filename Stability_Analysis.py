def Stability_Analysis(Tc, Pc, MW, w, BIC, z, T, P):
    
    # Inserting the necessary name spaces
    import numpy as np
    from EOS_SA import PR78_SA
    from K_value import K_Wilson
    
    # Defining the reduced properties
    Tr = T / Tc
    Pr = P / Pc
    
    # Calculating fugacity
    f_mix = PR78_SA(Tc, Pc, w, MW, T, P, BIC, z)
    
    # Approximating K_est
    K_est_SA = K_Wilson(w, Tr, Pr)
    
    # Calculating 2nd phase mole numbers
    Ki_v = K_est_SA
    Ki_l = K_est_SA
    itr = 0
    Err_v = 1
    Err_l = 1
    while Err_v > 1e-12 and Err_l > 1e-12:
        itr += 1
        
        # Initial guess
        Kin_v = Ki_v
        Kin_l = Ki_l
        Yi_v = z * Kin_v
        Yi_l = z / Kin_l
        
        # Summation of mole numbers
        Sv = sum(Yi_v)
        Sl = sum(Yi_l)
        
        # Normalizing the second phase mole numbers to calculate yi
        yi_v = Yi_v / Sv
        yi_l = Yi_l / Sl
        
        # Calculating fyi_v or fyi_l from EOS
        fyi_v = PR78_SA(Tc, Pc, w, MW, T, P, BIC, yi_v)
        fyi_l = PR78_SA(Tc, Pc, w, MW, T, P, BIC, yi_l)
        # Calculating  the fugasity ratio correction
        Ri_v = (f_mix / fyi_v) * (1 / Sv)
        Ri_l = (fyi_l / f_mix) * Sl
        if (itr + 1) % 5 == 0:
            Ri_v_nb = Ri_v
            Ri_l_nb = Ri_l
        
        # Checking convergenc
        Err_v = sum(Ri_v - 1) ** 2
        Err_l = sum(Ri_l - 1) ** 2
        
        # Updating K-values
        if itr % 5 == 0:
            b01_v = sum(np.log(Ri_v) * np.log(Ri_v_nb))
            b11_v = sum(np.log(Ri_v) * np.log(Ri_v))
            lambda_v = abs(b11_v / (b11_v - b01_v))
            b01_l = sum(np.log(Ri_l) * np.log(Ri_l_nb))
            b11_l = sum(np.log(Ri_l) * np.log(Ri_l))
            lambda_l = abs(b11_l / (b11_l - b01_l))
            Ki_v *= np.power(Ri_v, lambda_v)
            Ki_l *= np.power(Ri_l, lambda_l)
        else:
            Ki_v = Kin_v * Ri_v
            Ki_l = Kin_l * Ri_l
    
    # Creating results and checking stability
    TS_v = sum(np.power(np.log(Ki_v), 2))
    TS_l = sum(np.power(np.log(Ki_l), 2))
    
    # Stable phase
    if TS_l < 1e-4 and TS_v < 1e-4:
        Stability_Analysis_Output = {'I':1, 'Sl':Sl, 'Sv':Sv, 'TS_l':TS_l, 'TS_v':TS_v, \
                  'Ki':[Ki_l, Ki_v]}
    elif Sl - 1 <= 1e-4 and TS_v < 1e-4:
        Stability_Analysis_Output = {'I':1, 'Sl':Sl, 'Sv':Sv, 'TS_l':TS_l, 'TS_v':TS_v, \
                  'Ki':[Ki_l, Ki_v]} 
    elif TS_l < 1e-4 and Sv - 1 <= 1e-4:
        Stability_Analysis_Output = {'I':1, 'Sl':Sl, 'Sv':Sv, 'TS_l':TS_l, 'TS_v':TS_v, \
                  'Ki':[Ki_l, Ki_v]} 
    elif Sv - 1 <= 1e-4 and Sl - 1 <= 1e-4:
        Stability_Analysis_Output = {'I':1, 'Sl':Sl, 'Sv':Sv, 'TS_l':TS_l, 'TS_v':TS_v, \
                  'Ki':[Ki_l, Ki_v]}
    
    # Unstable phase
    if Sv - 1 > 1e-4 and TS_l > 1e-4:
        Stability_Analysis_Output = {'I':-1, 'Sl':Sl, 'Sv':Sv, 'TS_l':TS_l, 'TS_v':TS_v, \
                  'Ki':[Ki_l, Ki_v]}
    elif TS_l < 1e-4 and Sl - 1 > 1e-4:
        Stability_Analysis_Output = {'I':-1, 'Sl':Sl, 'Sv':Sv, 'TS_l':TS_l, 'TS_v':TS_v, \
                  'Ki':[Ki_l, Ki_v]}
    elif Sv - 1 > 1e-4 and Sl - 1 > 1e-4:
        Stability_Analysis_Output = {'I':-1, 'Sl':Sl, 'Sv':Sv, 'TS_l':TS_l, 'TS_v':TS_v, \
                  'Ki':[Ki_l, Ki_v]} 
    elif Sv - 1 > 1e-4 and Sl - 1 > 1e-4:
        Stability_Analysis_Output = {'I':-1, 'Sl':Sl, 'Sv':Sv, 'TS_l':TS_l, 'TS_v':TS_v, \
                  'Ki':[Ki_l, Ki_v]}
    elif Sv - 1 > 1e-4 and Sl - 1 <= 1e-4:
        Stability_Analysis_Output = {'I':-1, 'Sl':Sl, 'Sv':Sv, 'TS_l':TS_l, 'TS_v':TS_v, \
                  'Ki':[Ki_l, Ki_v]}
    elif Sv - 1 < 1e-4 and Sl - 1 > 1e-4:
        Stability_Analysis_Output = {'I':-1, 'Sl':Sl, 'Sv':Sv, 'TS_l':TS_l, 'TS_v':TS_v, \
                  'Ki':[Ki_l, Ki_v]}
    return Stability_Analysis_Output