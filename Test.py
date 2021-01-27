from Flash2PH import Flash2PH

CName = ['CO2','nC10']
COMP = [0.8, 0.2]
TEMP = 313.15
PRESS = 6E+6
KIJ = [[0, 0.01312], [0.01312, 0]]
VSH = [0, 0]

Result = Flash2PH(CName, COMP, TEMP, PRESS, KIJ, VSH)

for key, value in Result.items():
    print(key, ': ', value)