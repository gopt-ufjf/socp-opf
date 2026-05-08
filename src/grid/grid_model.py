import numpy as np
import pandas as pd
import src.grid.grid_utils as gu
from src.opf import socp_opf


class Grid:

    def __init__(self, path):
        # Leitura dos arquivos csv
        self.DBAR = pd.read_csv(path + "/DBAR.CSV", sep=";", index_col=0, decimal=',')
        self.DLIN = pd.read_csv(path + "/DLIN.CSV", sep=";", index_col=['BUSK', 'BUSM'], decimal=',')
        self.DGER = pd.read_csv(path + "/DGER.CSV", sep=";", index_col=0, decimal=',')
        self.DTEN = pd.read_csv(path + "/DTEN.CSV", sep=";", index_col=0, decimal=',')
        self.DBASE = pd.read_csv(path + "/DBASE.CSV", sep=";", decimal=',')
        self.DSOP = pd.read_csv(path + "/DSOP.CSV", sep=";", index_col=0, decimal=',')
        self.DGD = pd.read_csv(path + "/DGD.CSV", sep=";", index_col=0, decimal=',')
        self.DGD = self.DGD.astype(float)
        self.DBAT = pd.read_csv(path + "/DBAT.CSV", sep=";", index_col=0, decimal=',')
        self.CURVAP = pd.read_csv(path + "/CURVAP.CSV", sep=";", index_col=0, decimal=',')
        self.CURVAQ = pd.read_csv(path + "/CURVAQ.CSV", sep=";", index_col=0, decimal=',')

        # Passagem dos parâmetros para P.U.
        self.DBAR.PL = self.DBAR.PL / self.DBASE.SBASE[0]
        self.DBAR.QL = self.DBAR.QL / self.DBASE.SBASE[0]
        self.DBAR.CL = self.DBAR.CL / self.DBASE.SBASE[0]
        self.DLIN[["R", "X"]] = self.DLIN[["R", "X"]] / self.DBASE.ZBASE[0]
        self.DLIN[["G", "B", "BSH"]] = self.DLIN[["G", "B", "BSH"]] * self.DBASE.ZBASE[0]
        self.DGER[["PMIN", "PMAX", "QMIN", "QMAX"]] = self.DGER[["PMIN", "PMAX", "QMIN", "QMAX"]] / self.DBASE.SBASE[0]
        self.DSOP.SMAX = self.DSOP.SMAX / self.DBASE.SBASE[0]
        self.DGD.PMAX = self.DGD.PMAX / self.DBASE.SBASE[0]
        self.DBAT.PMAX = self.DBAT.PMAX / self.DBASE.SBASE[0]
        self.DBAT.EMAX = self.DBAT.EMAX / self.DBASE.SBASE[0]  # 1 dia com 24 periodos (sem manipulação adicional)

        # Manipulações adicionais
        self.DBAR.QL = self.DBAR.QL - self.DBAR.CL

        # Linhas conectadas à subestação (SE única)
        self.DLINSE = self.DLIN[
            (self.DLIN.index.get_level_values(0) == self.DGER.index[0]) |
            (self.DLIN.index.get_level_values(1) == self.DGER.index[0])
        ]

        # Cálculo do número de períodos
        self.PER = np.arange(1, self.CURVAP.shape[1] + 1)

        # Normalização das curva de Geração
        max_gd = self.DGD.iloc[:, 1:].to_numpy().max()
        self.DGD.iloc[:, 1:] = self.DGD.iloc[:, 1:] / max_gd

        # Tratamento das baterias
        self.DBAT.SOCMIN = (self.DBAT.SOCMIN / 100) * self.DBAT.EMAX
        self.DBAT.SOCMAX = (self.DBAT.SOCMAX / 100) * self.DBAT.EMAX
        self.DBAT.SOCINIC = (self.DBAT.SOCINIC / 100) * (self.DBAT.SOCMAX - self.DBAT.SOCMIN) + self.DBAT.SOCMIN
        self.DBAT.EFFCH = self.DBAT.EFFCH / 100
        self.DBAT.EFFDCH = self.DBAT.EFFDCH / 100

        # Criação das estruturas de mapeamento
        self.BSMP = gu.create_busmap(self.DBAR)

        # Formação da matriz Ybus
        self.YBUS = gu.calc_ybus(self)


    def run_opf(self, opf=1, options=None):

        if opf == 1:
            FOB, RBAR, RLIN = socp_opf.run(self, options)

        return FOB, RBAR, RLIN
