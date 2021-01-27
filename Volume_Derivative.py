def Volume_Derivative(V, a, am_l, Am_l, am_v, Am_v, b, bm_l, Bm_l, bm_v, \
                      Bm_v, Si_l, Si_v, Z_l, Z_v, P, T, BIC, x, y, Vp):

    # Import necessary libraries and functions
    import numpy as np
    import math
    
     # Defining the gas global constant
    R = 8.31447
        
    # Initializing the parameters
    V_l = np.array(V[0])
    V_v = np.array(V[1])
    
    # Fugacity Derivative equations with respect to Mole Number (Liquid phase)
    n_ll = Vp / V_l
    Db_nl = (1 / n_ll) * (-bm_l + b)
    Da_nl = (2 / n_ll) * (-am_l + Si_l)
    DB_nl = (1 / n_ll) * (-Bm_l + (P * b / (R * T)))
    DA_nl = (2 / n_ll) * (-Am_l + (P / (R * T) ** 2) * Si_l)
    DZ_nl = ((Bm_l - Z_l) * DA_nl + ((Am_l - 2 * Bm_l - 3 * Bm_l ** 2) + \
                                        2 * (3 * Bm_l + 1) * Z_l - Z_l ** 2) \
                 * DB_nl) / (3 * Z_l ** 2 - 2 * (1 - Bm_l) * Z_l + \
                             Am_l - 3 * Bm_l ** 2 - 2 * Bm_l)                          
    DG_nl = (2 * math.sqrt(2) / (Z_l + (1 - math.sqrt(2)) * \
                                 Bm_l) ** 2) * (Z_l * DB_nl - Bm_l * DZ_nl)
    DSig_nl = np.zeros(len(a))
    DE_nl = np.zeros((len(a), len(a)))
    DD_nl = np.zeros((len(a), len(a)))
    DC_nl = np.zeros((len(a), len(a)))
    KronDel = np.eye(len(a))
    DF_nl = np.zeros((len(a), len(a)))
    CC_l = (x * P) / (Z_l - Bm_l)
    EE_l = -Am_l / (2 * math.sqrt(2) * Bm_l) * (2 / am_l * Si_l - b / bm_l)
    GG_l = (Z_l + (1 + math.sqrt(2)) * Bm_l) / \
        (Z_l + (1 - math.sqrt(2)) * Bm_l)
    for k in range(len(a)):
        for i in range(len(a)):
            DSig_nl[i] = (1 / n_ll) * ((1 - BIC[i][k]) * \
                                       math.sqrt(a[i] * a[k])- Si_l[i])
            DE_nl[i][k] = (1 / (2 * math.sqrt(2) * Bm_l ** 2)) * \
                ((b[i] / bm_l) - (2 / am_l) * Si_l[i]) * \
                    (Bm_l * DA_nl[k] - Am_l * DB_nl[k]) + Am_l / \
                        (2 * math.sqrt(2) * Bm_l) * \
                            (((-b[i] / bm_l ** 2) * Db_nl[k]) - \
                             (2 / am_l ** 2) * ((am_l * DSig_nl[i]) - \
                                                Si_l[i] * Da_nl[k]))
            DD_nl[i][k] = (b[i] / bm_l) * (DZ_nl[k] - \
                                           ((Z_l - 1) / (bm_l * n_ll)) * \
                                               (b[k] - bm_l))
            DC_nl[i][k] = ( P / (Z_l - Bm_l) ** 2) * \
                (((Z_l - Bm_l) / n_ll) * (KronDel[i][k] - x[i]) - \
                 x[i] * (DZ_nl[k] - DB_nl[k]))                    
            DF_nl[i][k] = 1 / CC_l[i] * DC_nl[i][k] + DD_nl[i][k] + \
                math.log(GG_l, math.e) * DE_nl[i][k] + EE_l[i] / \
                    GG_l * DG_nl[k]
                
    # Fugacity Derivative equations with respect to Mole Number (Vapor phase)
    n_vv = Vp / V_v
    Db_nv = (1 / n_vv) * (-bm_v + b)
    Da_nv = (2 / n_vv) * (-am_v + Si_v)
    DB_nv = (1 / n_vv) * (-Bm_v + (P * b / (R * T)))
    DA_nv = (2 / n_vv) * (-Am_v + (P / (R * T) ** 2) * Si_v)
    DZ_nv = ((Bm_v - Z_v) * DA_nv + ((Am_v - 2 * Bm_v - 3 * Bm_v ** 2) + \
                                        2 * (3 * Bm_v + 1) * Z_v - Z_v ** 2) \
                 * DB_nv) / (3 * Z_v ** 2 - 2 * (1 - Bm_v) * Z_v + \
                             Am_v - 3 * Bm_v ** 2 - 2 * Bm_v)                                               
    DG_nv = (2 * math.sqrt(2) / (Z_v + (1 - math.sqrt(2)) * \
                                 Bm_v) ** 2) * (Z_v * DB_nv - Bm_v * DZ_nv)
    DSig_nv = np.zeros(len(a))
    DE_nv = np.zeros((len(a), len(a)))
    DD_nv = np.zeros((len(a), len(a)))
    DC_nv = np.zeros((len(a), len(a)))
    KronDel = np.eye(len(a))
    DF_nv = np.zeros((len(a), len(a)))
    CC_v = (y * P) / (Z_v - Bm_v)
    EE_v = -Am_v / (2 * math.sqrt(2) * Bm_v) * (2 / am_v * Si_v - b / bm_v)
    GG_v = (Z_v + (1 + math.sqrt(2)) * Bm_v) / \
        (Z_v + (1 - math.sqrt(2)) * Bm_v)
    for k in range(len(a)):
        for i in range(len(a)):
            DSig_nv[i] = (1 / n_vv) * ((1 - BIC[i][k]) * \
                                       math.sqrt(a[i] * a[k])- Si_v[i])
            DE_nv[i][k] = (1 / (2 * math.sqrt(2) * Bm_v ** 2)) * \
                ((b[i] / bm_v) - (2 / am_v) * Si_v[i]) * \
                    (Bm_v * DA_nv[k] - Am_v * DB_nv[k]) + Am_v / \
                        (2 * math.sqrt(2) * Bm_v) * \
                            (((-b[i] / bm_v ** 2) * Db_nv[k]) - \
                             (2 / am_v ** 2) * ((am_v * DSig_nv[i]) - \
                                                Si_v[i] * Da_nv[k]))                                  
            DD_nv[i][k] = (b[i] / bm_v) * (DZ_nv[k] - \
                                           ((Z_v - 1) / (bm_v * n_vv)) * \
                                               (b[k] - bm_v))
            DC_nv[i][k] = ( P / (Z_v - Bm_v) ** 2) * \
                (((Z_v - Bm_v) / n_vv) * (KronDel[i][k] - y[i]) - \
                 y[i] * (DZ_nv[k] - DB_nv[k]))
            DF_nv[i][k] = 1 / CC_v[i] * DC_nv[i][k] + DD_nv[i][k] + \
                math.log(GG_v, math.e) * DE_nv[i][k] + EE_v[i] / \
                    GG_v * DG_nv[k]
                    
    # Fugacity Derivative equations with respect to pressure (Liquid phase)
    DZ_pl = (Am_l * (Bm_l - Z_l) + Bm_l * \
             (Am_l - 2 * Bm_l - 3 * Bm_l ** 2 + 2 * (3 * Bm_l + 1) * Z_l - \
              Z_l ** 2))/(P * (3 * Z_l ** 2 - 2*(1 - Bm_l) * \
                               Z_l + Am_l - 3 * Bm_l ** 2 - 2 * Bm_l)) 
    DC_pl = x / (Z_l - Bm_l) ** 2 * (Z_l - Bm_l - P * \
                                     (DZ_pl - (bm_l / (R * T))))
    DD_pl = b / bm_l * DZ_pl
    DG_pl = 2 * math.sqrt(2) / (Z_l + (1 - math.sqrt(2)) * Bm_l) ** 2 * \
        (Z_l * bm_l / (R * T) - Bm_l * DZ_pl)
    DF_pl = 1/CC_l * DC_pl + DD_pl + EE_l / GG_l * DG_pl

    # Fugacity Derivative equations with respect to pressure (Vapor phase)
    DZ_pv = (Am_v * (Bm_v - Z_v) + Bm_v * \
             (Am_v - 2 * Bm_v - 3 * Bm_v ** 2 + 2 * (3 * Bm_v + 1) * Z_v - \
              Z_v ** 2))/(P * (3 * Z_v ** 2 - 2*(1 - Bm_v) * \
                               Z_v + Am_v - 3 * Bm_v ** 2 - 2 * Bm_v)) 
    DC_pv = y / (Z_v - Bm_v) ** 2 * (Z_v - Bm_v - P * \
                                     (DZ_pv - (bm_v / (R * T))))
    DD_pv = b / bm_v * DZ_pv
    DG_pv = 2 * math.sqrt(2) / (Z_v + (1 - math.sqrt(2)) * Bm_v) ** 2 * \
        (Z_v * bm_v / (R * T) - Bm_v * DZ_pv)
    DF_pv = 1/CC_v * DC_pv + DD_pv + EE_v / GG_v * DG_pv
    
    # Derivative of total fluid volume with respect to component mole
    COEFF = np.zeros((len(a) * len(a), len(a) * len(a)))
    CONST = np.zeros((len(a) * len(a)))
    for i in range(len(a)):
        for j in range(len(a)):
            COEFF[i][j] = DF_nl[i][j] + DF_nv[i][j]
            COEFF[i + len(a)][j + len(a)] =  DF_nl[i][j] + DF_nv[i][j]
            CONST[i * len(a) + j] = DF_nv[j][i]
    Dn_Nil = np.matmul(np.linalg.inv(COEFF), CONST)
    Dn_Nil = np.reshape(Dn_Nil, (len(a), len(a)), 'F')
    Dn_Niv = np.eye(len(a)) - Dn_Nil
    Vti = np.zeros(len(a))
    Sum1 = 0
    Sum2 = 0
    for k in range(len(a)):
        for i in range(len(a)):
            Sum1 += (V_l + (n_ll * R * T / P) * DZ_nl[i]) * Dn_Nil[i][k]
            Sum2 += (V_v + (n_vv * R * T / P) * DZ_nv[i]) * Dn_Niv[i][k]
        Vti[k] = (Sum1 + Sum2) * 1000
    Vti = list(Vti)
    
    # Derivative of total fluid volume with respect to pressure
    COEFF = np.zeros((len(a), len(a)))
    CONST = np.zeros(len(a))
    for i in range(len(a)):
        for j in range(len(a)):
            COEFF[i][j] = DF_nl[i][j] + DF_nv[i][j]
        CONST[i] = DF_pv[i] + DF_pl[i]
    Dn_Pil = np.matmul(np.linalg.inv(COEFF), CONST)
    Dn_Piv = -Dn_Pil
    Sum1 = 0
    Sum2 = 0
    for i in range(len(a)):
        Sum1 += (V_l + (n_ll * R * T / P) * DZ_nl[i]) * Dn_Pil[i] + \
            (n_ll * (R * T / P ** 2 * (P * DZ_pl - Z_l)))
        Sum2 += (V_v + (n_vv * R * T / P) * DZ_nv[i]) * Dn_Piv[i] + \
            (n_vv * (R * T / P ** 2 * (P * DZ_pv - Z_v)))
    DPl = Sum1
    DPv = Sum2
    DVP = DPl + DPv
    
    # Constructing and returning output dictionary
    Volume_Derivative_Output = {'DVP':DVP, 'Vti':Vti}
    return Volume_Derivative_Output