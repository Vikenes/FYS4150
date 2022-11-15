import numpy as np 

def Z(beta):
    return 12 + 4*np.cosh(8*beta)

def avg_E(beta):
    return -32*np.sinh(8*beta) / Z(beta)

def avg_eps(beta):
    return avg_E(beta)/4 

def avg_E2(beta):
    return 256/Z(beta) * np.cosh(8*beta) 

def avg_M(beta):
    return 8/Z(beta) * (2 + np.exp(8*beta))

def avg_M2(beta):
    return 32/Z(beta) * (1 + np.exp(8*beta))

def CV(beta):
    return 1/4 * (avg_E2(beta) - avg_E(beta)**2)

def chi(beta):
    return 1/4 * (avg_M2(beta) - avg_M(beta)**2)

beta = 1 

# print(f'E :  {avg_E(beta):.4f}')
# print(f'E2: {avg_E2(beta):.4f}')
# print(f'M :  {avg_M(beta):.4f}')
# print(f'M2: {avg_M2(beta):.4f}')
# print(f'Cv :    {CV(beta):.6f}')
# print(f'Chi:  { chi(beta):.6f}')

# print(f'e: {avg_eps(1):.4f}')
