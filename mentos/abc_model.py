import pandas as pd
import numpy as np

mets = ['A','B','C','D','E','F','A_ext', 'D_ext', 'E_ext','F_ext']
internal_mets = [m for m in mets if 'ext' not in m]
rxns = ['R_{}'.format(i) for i in range(1,10)]

data = {'R_1': pd.Series({'A':1, 'A_ext': -1}),
       'R_2': pd.Series({'A':-1,'B':1}),
       'R_3': pd.Series({'A':-1,'C':1}),
       'R_4': pd.Series({'B':-1, 'D': 2, 'E': -1}),
       'R_5': pd.Series({'E':1, 'E_ext':-1}),
       'R_6': pd.Series({'B':-2, 'C':1, 'F': 1}),
       'R_7': pd.Series({'C':-1, 'D':1}),
       'R_8': pd.Series({'D': -1, 'D_ext': 1}),
       'R_9': pd.Series({'F':-1, 'F_ext':1})}
fullS = pd.DataFrame(data, columns=rxns, index=mets,dtype='int64').fillna(0)

biomass = fullS.columns.get_loc('R_8')  # Index of biomass reaction

A_uptake, E_uptake = 0, 4   # Index of uptake reactions

R = 8.3144598/1000.0 # ideal gas constant
T = 298.15           # standard temperature
n_A = 6.022e23       # Avogadro's number
V = 1e-15            # volume of cell in Liters
c_L = 1e-8           # lower bound of metabolite concentrations
c_U = 1e-3           # upper bound of metabolite concentrations
v_L = -100
v_U = 100
A_ext = 6            # Index of A_ext
E_ext = 8            # Index of E_ext
D_ext = 7            # Index of D_ext
F_ext = 9            # Index of F_ext
A,B,C,D,E,F = range(6) # Index of internal metabolites
lambda_x = 0.5
S = fullS.loc[internal_mets].as_matrix() 

external_mets = [met for met in fullS.index if 'ext' in met]

m,n = fullS.shape

metab = {}
reactions = {}
true_metab, true_reactions = {},{}

mu0 = pd.Series([0.0,-2,-2,-4.0,-2.,-6.0,0.0,-4.0,-2.0,-6.0], index=mets,dtype='float64')
deltaG0 = fullS.T.dot(mu0)


met_bounds = pd.Series({'A_ext':c_U, 'E_ext': c_U, 'F_ext': c_L, 'D_ext': c_L}, index=external_mets)

efflux = [fullS.index.get_loc(met) for met in met_bounds[met_bounds == c_L].index]
uptake = [fullS.index.get_loc(met) for met in met_bounds[met_bounds == c_U].index]
internal = [fullS.index.get_loc(met) for met in internal_mets]
mu_ext = mu0[external_mets] + R*T*met_bounds.apply(np.log)
external_free_energy =  (mu_ext['D_ext'] + mu_ext['F_ext']) - (mu_ext['A_ext'] + mu_ext['E_ext'])


