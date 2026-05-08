import numpy as np


def create_busmap(dbar):

    busmap = {dbar.index.tolist()[k]: k for k in range(len(dbar))}

    return busmap


def calc_ybus(grid):

    ybus = np.zeros((len(grid.DBAR), len(grid.DBAR)), dtype=complex)

    for lin in grid.DLIN.index:
        ybus[grid.BSMP[lin[0]]][grid.BSMP[lin[1]]] += -(
                    grid.DLIN.ST[lin]*complex(grid.DLIN.G[lin], grid.DLIN.B[lin]) + grid.DLIN.BSH[lin])
        ybus[grid.BSMP[lin[1]]][grid.BSMP[lin[0]]] += -(
                    grid.DLIN.ST[lin]*complex(grid.DLIN.G[lin], grid.DLIN.B[lin]) + grid.DLIN.BSH[lin])
        ybus[grid.BSMP[lin[0]]][grid.BSMP[lin[0]]] += (
                grid.DLIN.ST[lin]*complex(grid.DLIN.G[lin], grid.DLIN.B[lin]) + grid.DLIN.BSH[lin])
        ybus[grid.BSMP[lin[1]]][grid.BSMP[lin[1]]] += (
                grid.DLIN.ST[lin]*complex(grid.DLIN.G[lin], grid.DLIN.B[lin]) + grid.DLIN.BSH[lin])

    return ybus


def create_connect(grid):

    cnct = {bus: [bus] for bus in grid.DBAR.index}

    for lin in grid.DLIN.index:
        if grid.DLIN.ST[lin] > 0:
            cnct[lin[0]].append(lin[1])
            cnct[lin[1]].append(lin[0])

    return cnct


def create_connect_conic(grid):

    cnct = {bus: [] for bus in grid.DBAR.index}

    for lin in grid.DLIN.index:
        if grid.DLIN.ST[lin] > 0:
            cnct[lin[0]].append(lin[1])
            cnct[lin[1]].append(lin[0])

    return cnct

    pass