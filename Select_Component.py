def Select_Component(CompName):
    # Inserting necessary name spaces
    import pandas as pd
    
    # Importing pure component data bank
    Pure_Data = pd.read_excel(r'PureDataProperties.xlsx', header = 0)
    
    # Selecting the properties of components
    MW = []
    Tc = []
    Pc = []
    Vc = []
    Zc = []
    w = []
    PCH = []
    for i in CompName:
        MW.append(Pure_Data[i][1])
        Tc.append(Pure_Data[i][2])
        Pc.append(Pure_Data[i][3])
        Vc.append(Pure_Data[i][4])
        Zc.append(Pure_Data[i][5])
        w.append(Pure_Data[i][6])
        PCH.append(Pure_Data[i][7])
    
    # Constructing and returning output dictionary
    Props = {'MW':MW, 'Tc':Tc, 'Pc':Pc, 'Vc':Vc, 'Zc':Zc, 'w':w, 'PCH':PCH}
    return Props