from src.grid.grid_model import Grid
import numpy as np



###########################################################################
# Configuração do 33 barras
###########################################################################

# Declaração do objeto grid
path = "data/33bus"
grid = Grid(path)

# Limite de fluxo ativo passante
grid.DLIN.loc[:,'PMAX'] = 97

# Ativação das GDs
grid.DBAR.loc[23, 'GD'] = 1
grid.DBAR.loc[24, 'GD'] = 1
grid.DBAR.loc[31, 'GD'] = 1

# # Ativação das baterias
# grid.DBAR.loc[23, 'BAT'] = 3
# grid.DBAR.loc[24, 'BAT'] = 3
# grid.DBAR.loc[31, 'BAT'] = 3

# # Ativação dos SOPs para o 33 barras
# grid.DLIN.loc[(7,20), 'ST'] = 2             # linha 7-20
# grid.DLIN.loc[(8,14), 'ST'] = 2             # linha 8-14
# grid.DLIN.loc[(11,21), 'ST'] = 2            # linha 11-21
# grid.DLIN.loc[(17,32), 'ST'] = 2            # linha 17-32
# grid.DLIN.loc[(24,28), 'ST'] = 2            # linha 24-28



###########################################################################
# Configuração do 94 barras
###########################################################################

# # Declaração do objeto grid
# path = "data/94bus"
# grid = Grid(path)
#
# # # Limite de fluxo ativo passante
# grid.DLIN.loc[:,'PMAX'] = 1
#
# # Ativação das GDs
# grid.DBAR.loc[71, 'GD'] = 1
# grid.DBAR.loc[79, 'GD'] = 1
# grid.DBAR.loc[28, 'GD'] = 1
#
# # # Ativação das baterias
# grid.DBAR.loc[71, 'BAT'] = 2
# grid.DBAR.loc[79, 'BAT'] = 2
# grid.DBAR.loc[28, 'BAT'] = 2
#
# # # Ativação dos SOPs para o 94 barras
# grid.DLIN.loc[(5,55), 'ST'] = 2    # caso 2
# grid.DLIN.loc[(7,60), 'ST'] = 2     # caso 3
# grid.DLIN.loc[(11,43), 'ST'] = 2    # caso 4
# grid.DLIN.loc[(12,72), 'ST'] = 2    # caso 5
# grid.DLIN.loc[(13,76), 'ST'] = 2    # caso 6
# grid.DLIN.loc[(14,18), 'ST'] = 2    # caso 7
# grid.DLIN.loc[(16,26), 'ST'] = 2    # caso 8
# grid.DLIN.loc[(20,83), 'ST'] = 2    # caso 9
# grid.DLIN.loc[(28,32), 'ST'] = 2    # caso 10
# grid.DLIN.loc[(29,39), 'ST'] = 2    # caso 11
# grid.DLIN.loc[(34,46), 'ST'] = 2    # caso 12
# grid.DLIN.loc[(40,42), 'ST'] = 2    # caso 13
# grid.DLIN.loc[(53,64), 'ST'] = 2    # caso 14



###########################################################################
# Execução do fluxo de potência ótimo
###########################################################################

''' opções opf
opf = 1 => socp_opf
'''

solver = 'ipopt'

# PERCENTUAL DE CORTES (MANUAL)
percc = 1

# PERCENTUAL DE PERDAS (AUTOMÁTICO)
percl = 1 - percc

options = {
    'percc': percc,
    'percl': percl,
    'solver':solver,
    'ipopt_maxiter':30000,
    'gurobi_time_limit':120
}

opf = 1
FOB, RBAR, RLIN = grid.run_opf(opf=opf, options=options)

FOB = FOB*grid.DBASE.SBASE[0]
if opf == 1:
    pos = RBAR.columns.get_loc('u')
    RBAR.insert(pos+1, 'v', np.sqrt(np.sqrt(2)*RBAR['u']))
RBAR.loc[:, 'pg':] *= grid.DBASE.SBASE[0]

FOB.to_csv(f"RFOB_{len(grid.DBAR)}.csv")
RBAR.to_csv(f"RBAR_{len(grid.DBAR)}.csv", index=True)
RLIN.to_csv(f"RLIN_{len(grid.DBAR)}.csv", index=True)





