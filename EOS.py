def PR78(Tc, Pc, w, MW, Tr, Pr):
    
    # Inserting necessary name spaces
    import numpy as np
     
    # Defining the gas global constant
    R = 8.31447
    
    # Inserting the EOS constants
    omega_a = 0.4572355289
    omega_b = 0.0777960739
    
    # Calculating EOS parameters
    ac = omega_a * R ** 2 * Tc ** 2 / Pc
    b = omega_b * R * Tc / Pc
    m = []
    for i in range(len(MW)):
        if MW[i] < 0.142285:
            m.append(0.37464 + 1.5422 * w[i] - 0.26992 * w[i] ** 2)
        else:
            m.append(0.379642 + 1.48503 * w[i] - 0.164423 * w[i] ** 2 + \
                     0.016666 * w[i] ** 3)
    m = np.array(m)
    alpha = (1 + m * (1 - Tr ** 0.5)) ** 2
    a = ac * alpha
        
    # Constructing and returning output dictionary
    EOS_Output = {'a':a, 'b':b}
    return EOS_Output