def K_Wilson(w, Tr, Pr):
    
    # Inserting necessary libraries
    import numpy as np
    
    # Calculating K-value using Wilson correlation
    K_value_Output = (1 / Pr) * np.exp(5.37 * (1 + w) * (1 - 1 / Tr))
    
    # Returning output value
    return K_value_Output