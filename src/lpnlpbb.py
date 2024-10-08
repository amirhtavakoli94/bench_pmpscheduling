#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 11:07:46 2020

@author: Sophie Demassey, Gratien Bonvin
"""

import gurobipy as gp
from gurobipy import GRB
# import primalheuristic as ph
import time
from hydraulics import HydraulicNetwork
import graphic
import csv
from pathlib import Path

# every unfeasible integer node X is discarded with a nogood cut: |x-X|>=1
# RECORDSOL = False:
# feasible integer nodes are discarded with nogoodcuts too and incumbent/cutoff are updated/checked outside of Gurobi
# => ObjBound is not a valid lower bound
# RECORDSOL = True:
# feasible integer nodes are updated with bound cut: obj >= (realcost-eps) * (1-X)
# which should invalidate the current MILP solution (when its cost is strictly lower than the real feasible solution)
# the real feasible solution is provided as a heuristic solution to Gurobi
# Problem 1: Gurobi 9.1 allows to set solutions at MIPNODE but not at MIPSOL
# so the solution found at MIPSOL must be recorded to be set at the next MIPNODE event
# Problem 2: even if Gurobi accepts the provided solution, it does not seem to directly update the cutoff value

OUTDIR = Path("../output/")
IISFILE = Path(OUTDIR, f'modeliis.ilp')


def _attach_callback_data(model, instance, modes):
    model._instance = instance
    model._nperiods = instance.nperiods()
    model._network = HydraulicNetwork(instance, model.Params.FeasibilityTol)

    model._incumbent = GRB.INFINITY
    model._callbacktime = 0
    model._solutions = []
    model._intnodes = {'unfeas': 0, 'feas': 0, 'adjust': 0}

    model._recordsol = (modes["solve"] == "RECORD")

    if modes["adjust"] != "NOADJUST":
        print("the primal heuristic based on time period adjustment is currently deactivated !")
    model._adjustmode = None
    model._adjusttime = time.time()
    model._adjust_solutions = []

    vs = [*model._svar.values(), *model._qvar.values()]
    model._lastsol = {'node': -1, 'cost': GRB.INFINITY, 'vars': vs}

    model._clonemodel = model.copy()

def mycallback(m, where):

    # STOP if UB-LB < tol
    if where == GRB.Callback.MIP:
        if m._incumbent - m.cbGet(GRB.Callback.MIP_OBJBND) < m.Params.MIPGap * m._incumbent:
            print('Stop early - ', m.Params.MIPGap * 100, '% gap achieved')
            m.terminate()

    elif m._recordsol and where == GRB.Callback.MIPNODE:
        if m._lastsol['cost'] < GRB.INFINITY:
            lastcost = m._lastsol['cost']
            m._lastsol['cost'] = GRB.INFINITY
            assert m._lastsol['node'] == m.cbGet(GRB.Callback.MIPNODE_NODCNT)
            assert len(m._lastsol['vars']) == len(m._lastsol['vals'])
            oldbest = m.cbGet(GRB.Callback.MIPNODE_OBJBST)
            m.cbSetSolution(m._lastsol['vars'], m._lastsol['vals'])
            objval = m.cbUseSolution()
            print(f"MIPNODE #{int(m.cbGet(GRB.Callback.MIPNODE_NODCNT))} "
                  f"oldbest = {oldbest} "
                  f"set solution #{m.cbGet(GRB.Callback.MIPNODE_SOLCNT)}: {lastcost} -> {objval}")
            if objval >= GRB.INFINITY and abs(oldbest-lastcost) > m.Params.MIPGapAbs:
                cloneandchecksolution(m, m._lastsol['vals'])
                print("if MILP feasible then the solution must violate a lazy cut")
        # if not m._rootlb:
        #    m._rootlb = m.cbGet(GRB.Callback.MIPNODE_OBJBND)

    # at an integer solution
    elif where == GRB.Callback.MIPSOL:

        # note that this integer solution may not be LP-optimal (e.g. if computed by a heuristic)
        costmipsol = m.cbGet(GRB.Callback.MIPSOL_OBJ)
        bestmipsol = m.cbGet(GRB.Callback.MIPSOL_OBJBST)
        currentlb = m.cbGet(GRB.Callback.MIPSOL_OBJBND)
        currentnode = int(m.cbGet(GRB.Callback.MIPSOL_NODCNT))
        fstring = f"MIPSOL #{currentnode}: {costmipsol} [{currentlb:.2f}, " \
                  f"{'inf' if m._incumbent == GRB.INFINITY else round(m._incumbent, 2)}]" #, best = {bestmipsol}"

        # prevent to evaluate again the same plan... (warning: no cut so the MILP solution is considered as feasible)
        # 1/ set by useSolution at the next MIPNODE (then MIPSOL is called again)
        # 2/ recomputed after adding the bound cut (sometimes (why not always ??) MIPSOL is called again) (=> the 10*)
        # 3/ due to multithreading delays (the incumbent update seems to be communicated but not the cutoff value)
        # if costmipsol > m._incumbent - 10*m.Params.MIPGapAbs:
        #    print(fstring + f" obj >= incumbent {m._incumbent}... SKIP (MILP solution accepted) ")
        #    if m._lastsol['cost'] < GRB.INFINITY:
        #        print("last solution has not yet been posted to the solver... let gurobi record this one instead")
        #        m._lastsol['cost'] = GRB.INFINITY
        #        assert costmipsol
        #        sol = [0 if m.cbGetSolution(svar) < 0.5 else 1 for svar in m._svar.values()]
        #        assert sol == m._lastsol['vals'][:len(sol)], "plan has changed between 2 consecutive MIPSOL ???"
        #    return
        if m._lastsol['node'] == currentnode:
            print(fstring + " same node: no way to cut this MILP solution, let gurobi save it and skip NLP computation")
            if m._lastsol['cost'] < GRB.INFINITY:
                assert m._lastsol['cost'] == m._incumbent
                print(f"last feasible solution {m._incumbent} has not yet been posted")
                # if abs(costmipsol - m._incumbent) < 10 * m.Params.MIPGapAbs:
                # print("obj = last ... SKIP (probably ?! same plan / almost same solution)")
                # m._lastsol = GRB.INFINITY
                # else:
                # print("obj > last ... SKIP (not same plan but should be fathomed later)")
                # sol = [0 if m.cbGetSolution(svar) < 0.5 else 1 for svar in m._svar.values()]
                # print(f"same plan ? {sol == m._lastsol['vals'][:len(sol)]}")
            return

        m._starttime = time.time()
        nogood_lastperiod = m._nperiods
        costrealsol = GRB.INFINITY

        inactive, activity = getplan(m)
        qreal, hreal, vreal, violperiod = m._network.extended_period_analysis(inactive, stopatviolation=True)

        if violperiod:
            m._intnodes['unfeas'] += 1
            nogood_lastperiod = violperiod
            print(fstring + f" t={violperiod}")
            addnogoodcut(m, _linearnorm(m._svar, nogood_lastperiod, activity), currentnode)

        else:
            m._intnodes['feas'] += 1
            costrealsol = solutioncost(m, activity, qreal)
            print(fstring + f" feasible: {costrealsol}")

            if costrealsol < costmipsol - m.Params.MIPGapAbs:
                print(f"mip solution cost {costmipsol} (heuristic, non-lp optimal?) > real cost {costrealsol}")
                # solvecvxmodelwithsolution(m, getrealsol(m, activity, qreal))

            linexpr = _linearnorm(m._svar, nogood_lastperiod, activity)
            if m._recordsol:
                print(f"bound cut {costrealsol}")
#                print('nogood:', str(linexpr))
                addboundcut(m, linexpr, costrealsol, currentnode)
            else:
                addnogoodcut(m, linexpr, currentnode)

        if costrealsol < m._incumbent:
            m._incumbent = costrealsol
            m._solutions.append({'plan': activity, 'cost': costrealsol, 'flows': qreal, 'volumes': vreal,
                                 'cpu': m.cbGet(GRB.Callback.RUNTIME), 'adjusted': (qreal is None)})
            gap = (m._incumbent - currentlb) / m._incumbent
            print(f'UPDATE INCUMBENT gap={gap * 100:.4f}%')

            if m._recordsol:
                m._lastsol['cost'] = costrealsol
                m._lastsol['vals'] = getrealsol(m, activity, qreal)
                m._lastsol['node'] = m.cbGet(GRB.Callback.MIPSOL_NODCNT)
            else:
                addcutoff(m, (1 - m.Params.MIPGap) * m._incumbent, currentnode)

        m._callbacktime += time.time() - m._starttime


def getrealsol(m, activity, qreal):
    solx = [activity[t][a] for (a, t) in m._svar]
    solq = [qreal[t][a] for (a, t) in m._qvar]
    # for (j,t) in m._hvar:
    #     sol.append(hreal[t][j])
    return [*solx, *solq]

def getplan(m):
    inactive = {t: set() for t in range(m._nperiods)}
    activity = {t: {} for t in range(m._nperiods)}
    for (a, t), svar in m._svar.items():
        if m.cbGetSolution(svar) < 0.5:
            inactive[t].add(a)
            activity[t][a] = 0
        else:
            activity[t][a] = 1
    return inactive, activity

def addnogoodcut(m, linnorm, n):
    m.cbLazy(linnorm >= 1)
    if m._clonemodel:
        c = clonelinexpr(m, linnorm)
        m._clonemodel.addConstr(c >= 1, name=f'nogood{n}')

def addboundcut(m, linnorm, costrealsol, n):
    cost = costrealsol #- m.Params.MIPGapAbs
    m.cbLazy(m._obj >= cost * (1 - linnorm))
    if m._clonemodel:
        c = clonelinexpr(m, linnorm)
        cobj = m._clonemodel.getObjective()
        m._clonemodel.addConstr(cobj >= cost * (1 - c), name=f'bound{n}')

def addcutoff(m, cutoff, n):
    m.cbLazy(m._obj <= cutoff)
    if m._clonemodel:
        cobj = m._clonemodel.getObjective()
        m._clonemodel.addConstr(cobj <= cutoff, name=f'cutoff{n}')


def clonelinexpr(m, linexpr):
    assert m._clonemodel
    sz = linexpr.size()
    coeffs = [linexpr.getCoeff(i) for i in range(sz)]
    vars = [m._clonemodel.getVarByName(linexpr.getVar(i).varname) for i in range(sz)]
    c = gp.LinExpr(linexpr.getConstant())
    c.addTerms(coeffs, vars)
    assert c.size() == sz and c.getConstant() == linexpr.getConstant()
    return c


def _parse_activity(horizon, svars):
    inactive = {t: set() for t in horizon}
    activity = {t: {} for t in horizon}
    for (a, t), svar in svars.items():
        if svar.x < 0.5:
            inactive[t].add(a)
            activity[t][a] = 0
        else:
            activity[t][a] = 1
    return inactive, activity


def _linearnorm(svars, last_period, activity):
    linexpr = gp.LinExpr()
    nbact = 0

    for t in range(last_period):
        for a, active in activity[t].items():
            # _filterwithsymmetry(activity[t], a):
            if active:
                linexpr.addTerms(-1.0, svars[a, t])
                nbact += 1
            else:
                linexpr.addTerms(1.0, svars[a, t])
    linexpr.addConstant(nbact)
    return linexpr


def solutioncost(cvxmodel, status, flows):
    return sum(status[t][a] * cvxmodel._svar[a, t].obj
               + flows[t][a] * cvxmodel._qvar[a, t].obj for (a, t) in cvxmodel._svar)


def lpnlpbb(cvxmodel, instance, modes, drawsolution=True):
    """Apply the LP-NLP Branch-and-bound using the convex relaxation model cvxmodel."""

    _attach_callback_data(cvxmodel, instance, modes)
    cvxmodel.params.LazyConstraints = 1

    cvxmodel.optimize(mycallback)

    if cvxmodel.status != GRB.OPTIMAL:
        print('Optimization was stopped with status %d' % cvxmodel.status)

    if cvxmodel._recordsol and cvxmodel.solcount == 0:
        return 0, {}, 0, cvxmodel.ObjBound, cvxmodel._solutions

    print("check gurobi best solution")
    solx =  [v.x for v in cvxmodel._svar.values()]
    print(solx)
    # cvxmodel._clonemodel.write("check.lp")
    solq =  [v.x for v in cvxmodel._qvar.values()]
    cloneandchecksolution(cvxmodel, [*solx, *solq])

    cost = 0
    plan = {}
    volume= 0
    lower_boun=0
    solu = 0
    if cvxmodel._solutions:
        bestsol = cvxmodel._solutions[-1]
        plan = bestsol['plan']
        flow = bestsol['flows']
        volume = bestsol['volumes']
        cost = cvxmodel._incumbent
        lower_boun= cvxmodel.ObjBound
        solu = cvxmodel._solutions
        if not flow:
            print('best solution found by the time-adjustment heuristic')
            inactive = {t: set(a for a, act in activity_t.items() if not act) for t, activity_t in plan.items()}
            flow, hreal, volume, nbviolations = cvxmodel._network.extended_period_analysis(inactive, stopatviolation=False)
            assert nbviolations, 'solution was time-adjusted and should be slightly unfeasible'
            cost = solutioncost(cvxmodel, plan, flow)
            lower_boun= cvxmodel.ObjBound
            print(f'real plan cost = {cost} / time adjustment cost = {cvxmodel._incumbent}')
        if drawsolution:
            graphic.pumps(instance, flow)
            graphic.tanks(instance, flow, volume)
    else:
        print(f'no solution found write IIS file {IISFILE}')
        cvxmodel.computeIIS()
        cvxmodel.write(IISFILE)

    return cost, plan, volume, lower_boun, solu


def solveconvex(cvxmodel, instance, drawsolution=True):
    """Solve the convex relaxation model cvxmodel."""
    # cvxmodel.params.SolutionLimit = 1

    cvxmodel.optimize()

    if cvxmodel.status != GRB.OPTIMAL:
        print('Optimization was stopped with status %d' % cvxmodel.status)

    costreal = 0
    plan = {}
    if cvxmodel.SolCount:
        inactive, activity = _parse_activity(instance.horizon(), cvxmodel._svar)
        net = HydraulicNetwork(instance, cvxmodel.Params.FeasibilityTol)
        qreal, hreal, vreal, nbviolations = net.extended_period_analysis(inactive, stopatviolation=False)
        print(f"real plan with {nbviolations} violations")
        plan = {t: {a: (0 if abs(q) < 1e-6 else 1) for a, q in qreal[t].items()} for t in qreal}
        costreal = solutioncost(cvxmodel, plan, qreal)

        if drawsolution:
            graphic.pumps(instance, qreal)
            graphic.tanks(instance, qreal, vreal)

        for a, s in cvxmodel._svar.items():
            print(f'{a}: {round(s.x)} {round(cvxmodel._qvar[a].x, 4)}')
    else:
        print(f"f write IIS in {IISFILE}")
        cvxmodel.computeIIS()
        cvxmodel.write(IISFILE)

    return costreal, plan


def recordandwritesolution(m, activity, qreal, filename):
    sol = getrealsol(m, activity, qreal)
    f = open(filename, 'a')
    write = csv.writer(f)
    write.writerow(sol)
    return sol


def cloneandchecksolution(m, vals):
    vars = m._lastsol['vars']
    assert len(vals) == len(vars)
    model = m._clonemodel if m._clonemodel else m.copy()
    for i, var in enumerate(vars):
        clonevar = model.getVarByName(var.varname)
        clonevar.lb = vals[i]
        clonevar.ub = vals[i]
    print("test real solution / cvx model (without internal cuts nor processing)")
    model.optimize()
    if model.status != GRB.OPTIMAL:
        print('Optimization was stopped with status %d' % model.status)

    model.terminate()
