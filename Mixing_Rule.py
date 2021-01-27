def Mixing_Rule(a, b, comp, T, P, BIC):
    
    # Inserting necessary name spaces
    import math
    
    # Defining the gas global constant
    R = 8.31447
    
    # Mixing the properties
    bm = 0
    am = 0
    Si = []
    for i in range(len(a)):
        r = 0
        for j in range(len(a)):
            r += comp[j] * (1 - BIC[i][j]) * math.sqrt(a[j])
        Si.append(math.sqrt(a[i]) * r)
        am += comp[i] * Si[i]
        bm += b[i] * comp[i]
    Am = am * P / (R * T) ** 2
    Bm = bm * P / (R * T)
    
    # Constructing and returning output dictionary
    Mixing_Rule_Out = {'am':am, 'bm':bm, 'Am':Am, 'Bm':Bm, 'Si':Si}
    return Mixing_Rule_Out