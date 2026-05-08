import pyomo.environ as pyo
import numpy as np
from src.grid.grid_utils import create_connect_conic
from src.report.report import create_report
from time import time

def run(grid, options):
    t0 = time()

    ####################################################################################################################
    # Declaracão do modelo
    ####################################################################################################################

    model = pyo.ConcreteModel()



    ####################################################################################################################
    # Estruturas auxiliares necessárias
    ####################################################################################################################

    cnct = create_connect_conic(grid)
    ridx = grid.DLIN.index.tolist()
    tidx = ridx + [(m,k) for (k,m) in ridx]
    DLINSP = grid.DLIN[grid.DLIN.ST > 1]      # subconjunto de DLIN referente aos SOPs
    KSP = [idx[0] for idx in DLINSP.index]
    MSP = [idx[1] for idx in DLINSP.index]
    DBARGD = grid.DBAR[grid.DBAR.GD > 0]      # Subconjunto de DBAR referente às GDs
    DBARBAT = grid.DBAR[grid.DBAR.BAT > 0]    # Subconjunto de DBAR referente às baterias



    ####################################################################################################################
    # Declaracão das variáveis
    ####################################################################################################################

    # Tensões nodais (fase)
    def ph_bounds(model, bus, t):
        return (grid.DTEN.PHMIN[bus], grid.DTEN.PHMAX[bus])
    def ph_inic(model, bus, t):
        return grid.DBAR.PH[bus]
    model.ph = pyo.Var(grid.DBAR.index, grid.PER, bounds=ph_bounds, initialize=ph_inic)

    # sqrt(2)*V_k = u_k
    def u_bounds(model, bus, t):
        return ((grid.DTEN.VMIN[bus]**2/np.sqrt(2)), (grid.DTEN.VMAX[bus]**2/np.sqrt(2)))
    def u_inic(model, bus, t):
        # return 1/np.sqrt(2)
        return (grid.DTEN.VMAX[bus]**2/np.sqrt(2) + grid.DTEN.VMIN[bus]**2/np.sqrt(2))/2
    model.u = pyo.Var(grid.DBAR.index, grid.PER, bounds=u_bounds, initialize=u_inic)

    # R = V_k*V_m*cos(theta_km)
    def R_bounds(model):
        return (0, max(grid.DTEN.VMAX)**2)
    model.R = pyo.Var(ridx, grid.PER, domain=pyo.NonNegativeReals)

    # T = V_k*V_m*sen(theta_km)
    def T_bounds(model):
        return (-max(grid.DTEN.VMAX) ** 2, max(grid.DTEN.VMAX) ** 2)
    model.T = pyo.Var(tidx, grid.PER)

    # Geração na SE (ativa)
    def pg_bounds(model, bus, t):
        return (grid.DGER.PMIN[bus], grid.DGER.PMAX[bus])
    model.pg = pyo.Var(grid.DGER.index, grid.PER,  bounds=pg_bounds)

    # Geração na SE (reativa)
    def qg_bounds(model, bus, t):
        return (grid.DGER.QMIN[bus], grid.DGER.QMAX[bus])
    model.qg = pyo.Var(grid.DGER.index, grid.PER, bounds=qg_bounds)

    # Variáveis dos SOPs
    if len(DLINSP) != 0:
        model.psok = pyo.Var(KSP, grid.PER, initialize=1e-10)
        model.psom = pyo.Var(MSP, grid.PER, initialize=1e-10)
        model.qsok = pyo.Var(KSP, grid.PER, initialize=1e-10)
        model.qsom = pyo.Var(MSP, grid.PER, initialize=1e-10)
        model.psokl = pyo.Var(KSP, grid.PER, initialize=1e-10)
        model.psoml = pyo.Var(MSP, grid.PER, initialize=1e-10)

    # Variáveis das GDs (curtailment de GDs)
    if len(DBARGD) != 0:
        def cgd_bounds(model, bus, t):
            return (0, grid.DGD.PMAX[grid.DBAR.GD[bus]]*grid.DGD.loc[grid.DBAR.GD[bus], f'{t}'])
        model.cgd = pyo.Var(DBARGD.index, grid.PER, bounds=cgd_bounds)

    # Variáveis das baterias
    if len(DBARBAT)!= 0:
        def pbat_bounds(model, bus, t):
            return (0, grid.DBAT.PMAX[grid.DBAR.BAT[bus]])
        model.pch = pyo.Var(DBARBAT.index, grid.PER, bounds=pbat_bounds)
        model.pdch = pyo.Var(DBARBAT.index, grid.PER, bounds=pbat_bounds)
        def soc_bounds(model, bus, t):
            return (grid.DBAT.SOCMIN[grid.DBAR.BAT[bus]], grid.DBAT.SOCMAX[grid.DBAR.BAT[bus]])
        model.soc = pyo.Var(DBARBAT.index, grid.PER, bounds=soc_bounds)

    # Perdas e cortes
    model.perdas = pyo.Var()
    model.cortes = pyo.Var()



    ####################################################################################################################
    # Declaracão da função objetivo
    ####################################################################################################################

    def obj_rule(model):
        return options['percc']*model.cortes + options['percl']*model.perdas
    model.obj = pyo.Objective(rule=obj_rule, sense=pyo.minimize)



    ####################################################################################################################
    # Declaracão das restrições
    ####################################################################################################################

    # Recuperação das fases ############################################################################################

    def fases_rule(model,k,m,t):
        return -model.ph[k,t] + model.ph[m,t] + pyo.atan(model.T[k,m,t]/model.R[k,m,t]) == 0
    model.fases = pyo.Constraint(grid.DLIN.index, grid.PER, rule=fases_rule)

    # Cálculo das perdas e cortes ######################################################################################

    def calc_perdas_rule(model):
        perdas_aux = sum(
                        grid.DLIN.ST[lin]*(
                        grid.DLIN.G[lin]*(np.sqrt(2)*model.u[lin[0],t] + np.sqrt(2)*model.u[lin[1],t] - 2*model.R[lin,t])
                ) for lin in grid.DLIN.index for t in grid.PER
            )
        return model.perdas == perdas_aux
    model.calc_perdas = pyo.Constraint(rule=calc_perdas_rule)

    def calc_cortes_rule(model):
        return model.cortes == sum(model.cgd[bus, t] for bus in DBARGD.index for t in grid.PER)
    model.calc_cortes = pyo.Constraint(rule=calc_cortes_rule)


    # Balanço de potência (ativa) ######################################################################################

    def balp_rule(model, k, t):

        pgen = 0; pflux = 0
        psopk = 0; psopm = 0
        pgend = 0
        pbatch = 0; pbatdch = 0

        # Checagem da SE
        if k in grid.DGER.index:
            pgen += model.pg[k, t]

        # Checagem SOP
        if k in KSP:
            psopk += model.psok[k,t]
        if k in MSP:
            psopm += model.psom[k,t]

        # Checagem GD
        if k in DBARGD.index:
            pgend += grid.DGD.PMAX[grid.DBAR.GD[k]]*grid.DGD.loc[grid.DBAR.GD[k], f'{t}'] - model.cgd[k,t]

        # Checagem BAT
        if k in DBARBAT.index:
            pbatch += model.pch[k,t]
            pbatdch += model.pdch[k,t]

        # Fluxos passantes cônicos
        for m in cnct[k]:
            if (k,m) in ridx:
                pflux += (grid.YBUS[grid.BSMP[k]][grid.BSMP[m]].real*model.R[k,m,t] +
                          grid.YBUS[grid.BSMP[k]][grid.BSMP[m]].imag*model.T[k,m,t])
            else:
                pflux += (grid.YBUS[grid.BSMP[k]][grid.BSMP[m]].real*model.R[m,k,t] +
                          grid.YBUS[grid.BSMP[k]][grid.BSMP[m]].imag*model.T[k,m,t])

        return (
            pgen + psopk + psopm + pgend + pbatdch - pbatch - grid.CURVAP.loc[k, f'{t}']*grid.DBAR.PL[k] -
            np.sqrt(2)*grid.YBUS[grid.BSMP[k]][grid.BSMP[k]].real*model.u[k,t] - pflux == 0
        )

    model.balp = pyo.Constraint(grid.DBAR.index, grid.PER, rule=balp_rule)


    # Balanço de potência (reativa) ####################################################################################

    def balq_rule(model, k, t):

        qgen = 0; qflux = 0
        qsopk = 0; qsopm = 0

        # Checagem da SE
        if k in grid.DGER.index:
            qgen += model.qg[k,t]

        # Checagem SOP
        if k in KSP:
            qsopk += model.qsok[k,t]
        if k in MSP:
            qsopm += model.qsom[k,t]

        # Fluxos passantes cônicos
        for m in cnct[k]:
            if (k,m) in ridx:
                qflux += (grid.YBUS[grid.BSMP[k]][grid.BSMP[m]].imag*model.R[(k,m,t)] -
                          grid.YBUS[grid.BSMP[k]][grid.BSMP[m]].real*model.T[(k,m,t)])
            else:
                qflux += (grid.YBUS[grid.BSMP[k]][grid.BSMP[m]].imag*model.R[(m,k,t)] -
                          grid.YBUS[grid.BSMP[k]][grid.BSMP[m]].real*model.T[(k,m,t)])

        return (
            qgen + qsopk + qsopm - grid.CURVAQ.loc[k, f'{t}']*grid.DBAR.QL[k] +
            np.sqrt(2)*grid.YBUS[grid.BSMP[k]][grid.BSMP[k]].imag*model.u[k,t] + qflux == 0
        )
    model.balq = pyo.Constraint(grid.DBAR.index, grid.PER, rule=balq_rule)


    # Restrições cônicas - cone rotacionado ############################################################################

    def rotc_rule(model, k, m, t):
        return 2*model.u[k,t]*model.u[m,t] == model.R[k,m,t]**2 + model.T[k,m,t]**2
    model.rotc = pyo.Constraint(grid.DLIN.index, grid.PER, rule=rotc_rule)


    # Restrições cônicas - variável T ##################################################################################

    def restt_rule(model, k, m, t):
        return model.T[k,m,t] == - model.T[m,k,t]
    model.restt = pyo.Constraint(grid.DLIN.index, grid.PER, rule=restt_rule)


    # Restrições cônicas - fluxo passante ##############################################################################

    def fluxp1(model, k, m, t):
        return np.sqrt(2)*grid.DLIN.G[k,m]*model.u[k,t] - grid.DLIN.G[k,m]*model.R[k,m,t] - grid.DLIN.B[k,m]*model.T[k,m,t] >= -grid.DLIN.PMAX[k,m]
    model.fluxp1 = pyo.Constraint(grid.DLIN.index, grid.PER, rule=fluxp1)

    def fluxp2(model, k, m, t):
        return np.sqrt(2)*grid.DLIN.G[k,m]*model.u[k,t] - grid.DLIN.G[k,m]*model.R[k,m,t] - grid.DLIN.B[k,m]*model.T[k,m,t] <= grid.DLIN.PMAX[k,m]
    model.fluxp2 = pyo.Constraint(grid.DLIN.index, grid.PER, rule=fluxp2)

    def fluxp3(model, k, m, t):
        return np.sqrt(2)*grid.DLIN.G[k,m]*model.u[k,t] - grid.DLIN.G[k,m]*model.R[k,m,t] + grid.DLIN.B[k,m]*model.T[k,m,t] >= -grid.DLIN.PMAX[k,m]
    model.fluxp3 = pyo.Constraint(grid.DLIN.index, grid.PER, rule=fluxp3)

    def fluxp4(model, k, m, t):
        return np.sqrt(2)*grid.DLIN.G[k,m]*model.u[k,t] - grid.DLIN.G[k,m]*model.R[k,m,t] + grid.DLIN.B[k,m]*model.T[k,m,t] <= grid.DLIN.PMAX[k,m]
    model.fluxp4 = pyo.Constraint(grid.DLIN.index, grid.PER, rule=fluxp4)


    # Restrições dos SOPs ##############################################################################################

    def balpsp_rule(model, k, m, t):
        return model.psok[k,t] + model.psom[m,t] + model.psokl[k,t] + model.psoml[m,t] == 0
    model.balpsp = pyo.Constraint(DLINSP.index, grid.PER, rule=balpsp_rule)

    def spconk_rule(model, k, m, t):
        return (
                model.psok[k,t]**2 + model.qsok[k,t]**2 == 2*
                (model.psokl[k,t]/(np.sqrt(2)*grid.DSOP.PERDAS[DLINSP.ST[k,m]]))*
                (model.psokl[k,t]/(np.sqrt(2)*grid.DSOP.PERDAS[DLINSP.ST[k,m]]))
        )
    model.spconk = pyo.Constraint(DLINSP.index, grid.PER, rule=spconk_rule)

    def spconm_rule(model, k, m, t):
        return (
                model.psom[m,t]**2 + model.qsom[m,t]**2 == 2*
                (model.psoml[m,t]/(np.sqrt(2)*grid.DSOP.PERDAS[DLINSP.ST[k,m]]))*
                (model.psoml[m,t]/(np.sqrt(2)*grid.DSOP.PERDAS[DLINSP.ST[k,m]]))
        )
    model.spconm = pyo.Constraint(DLINSP.index, grid.PER, rule=spconm_rule)

    def splimk_rule(model, k, m, t):
        return (
                model.psok[k,t]**2 + model.qsok[k,t]**2 <= 2*
                (grid.DSOP.SMAX[DLINSP.ST[k, m]]/np.sqrt(2))*(grid.DSOP.SMAX[DLINSP.ST[k, m]]/np.sqrt(2))
        )
    model.splimk = pyo.Constraint(DLINSP.index, grid.PER, rule=splimk_rule)

    def splimm_rule(model, k, m, t):
        return (
                model.psom[m,t]**2 + model.qsom[m,t]**2 <= 2*
                (grid.DSOP.SMAX[DLINSP.ST[k, m]]/np.sqrt(2))*(grid.DSOP.SMAX[DLINSP.ST[k, m]]/np.sqrt(2))
        )
    model.splimm = pyo.Constraint(DLINSP.index, grid.PER, rule=splimm_rule)


    # Restrições das baterias ##########################################################################################

    # Complementariedade
    def comp_rule(modek, k, t):
        return model.pch[k,t]*model.pdch[k,t] == 0
    model.comp = pyo.Constraint(DBARBAT.index, grid.PER, rule=comp_rule)

    # Atualização do SOC
    def atsoc_rule(model, k, t):
        if t > 1:
            return (
                    model.soc[k,t] == model.soc[k,t-1] +
                    grid.DBAT.EFFCH[DBARBAT.BAT[k]]*model.pch[k,t-1] -
                    (1/grid.DBAT.EFFDCH[DBARBAT.BAT[k]])*model.pdch[k,t-1]
            )
        else:
            return model.soc[k, t] == grid.DBAT.SOCINIC[DBARBAT.BAT[k]]
    model.atsoc = pyo.Constraint(DBARBAT.index, grid.PER, rule=atsoc_rule)

    tf = time() - t0


    ####################################################################################################################
    # Resolução do modelo
    ####################################################################################################################

    solver = pyo.SolverFactory(options['solver'])
    if options['solver'] == 'ipopt':
        if options['ipopt_maxiter'] != None:
            solver.options['max_iter'] = options['ipopt_maxiter']

    try:
        results = solver.solve(model, tee=True)
    except ValueError:
        results = 'ERRO'


    ####################################################################################################################
    # Relatório de resolução
    ####################################################################################################################

    print()
    print()
    print('#' * 50)
    print('Relatório'.upper())
    print('#' * 50)
    print()
    print(f'Status: {results.solver.termination_condition}')
    print(f'Função Objetivo: {pyo.value(model.obj)*grid.DBASE.SBASE[0]} _W')
    print(f'Perdas: {pyo.value(model.perdas)*grid.DBASE.SBASE[0]} _W')
    print(f'Cortes: {pyo.value(model.cortes)*grid.DBASE.SBASE[0]} _W')
    print(f'Tempo de geração: {tf} s')
    print(f'Tempo de resolução: {results.solver.time} s')

    FOB, RBAR, RLIN = create_report(model)
    return FOB, RBAR, RLIN










