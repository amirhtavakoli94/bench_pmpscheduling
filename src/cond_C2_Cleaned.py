# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 19:01:10 2023

@author: amirhossein.tavakoli
"""

import gurobipy as gp
from gurobipy import GRB
import outerapproximation as oa
from instance import Instance
import datetime as dt
import numpy as np
import math
from math import sqrt


def beta(q, coeff):
    return coeff[2] * q * abs(q) + coeff[1] * q + coeff[0]


def beta1(q, coeff):
    return -coeff[2] * q * abs(q) + (-coeff[1]) * q - coeff[0]


# !!! check round values
# noinspection PyArgumentList
def build_common_model(inst: Instance, Z, C, D, P0, P1, m_, t, mode_BT, oa_types, accuracy: float, envelop: float, oagap: float, config = [], dh = {}, Discrete = False):
    """Build the convex relaxation gurobi model."""

    milp = gp.Model('Pumping_Scheduling')
#    milp.params.NonConvex = 2
    milp.update()

    qvar = {}  # arc flow
    dhvar = {}  # arc hloss
    svar = {}  # arc status
    hvar = {}  # node head
    qexpr = {}  # node inflow

    nperiods = inst.nperiods()
    horizon = inst.horizon()


    for j in inst.junctions:
        hvar[j, t] = milp.addVar(name=f'hj({j},{t})')

    for j, res in inst.reservoirs.items():
        hvar[j, t] = milp.addVar(lb=res.head(t), ub=res.head(t), name=f'hr({j},{t})')

    for j, tank in inst.tanks.items():
        lbt = tank.head(tank.vinit) if t == 0 else D[j, t][0]
        ubt = tank.head(tank.vinit) if t == 0 else D[j, t][1]
        hvar[j, t] = milp.addVar(lb=lbt, ub=ubt, name=f'ht({j},{t})')
        hvar[j, t+1] = milp.addVar(lb=tank.head(tank.vmin), ub=tank.head(tank.vmax), name=f'ht({j},{t})')
        if t != nperiods-1:
            hvar[j, t+1] = milp.addVar(lb=D[j, t+1][0], ub=D[j, t+1][1], name=f'ht({j},{t})')
            
            
        if t>=1:
            
            milp.addConstr(hvar[j, t] <= tank.hmax[t])
            milp.addConstr(hvar[j, t] >= tank.hmin[t])
        

    milp.update()

    for (i, j), a in inst.arcs.items():

                    
        if a.control:

            if ((i, j), t) in P0:
                svar[(i, j), t] = milp.addVar(ub=0, vtype=GRB.BINARY, name=f'x({i},{j},{t})')
            elif ((i, j), t) in P1:
                svar[(i, j), t] = milp.addVar(lb=1, vtype=GRB.BINARY, name=f'x({i},{j},{t})')
            else:
                    svar[(i, j), t] = milp.addVar(vtype=GRB.BINARY, name=f'x({i},{j},{t})')

            
            qvar[(i, j), t] = milp.addVar(lb=-GRB.INFINITY, name=f'q({i},{j},{t})')
            dhvar[(i, j), t] = milp.addVar(lb=-GRB.INFINITY, name=f'H({i},{j},{t})')
            milp.addConstr(qvar[(i, j), t] <= Z[(i, j), t][1] * svar[(i, j), t], name=f'qxup({i},{j},{t})')
            milp.addConstr(qvar[(i, j), t] >= Z[(i, j), t][0] * svar[(i, j), t], name=f'qxlo({i},{j},{t})')
            dhmin = max(a.hlossval(a.qmin), hvar[i, t].lb - hvar[j, t].ub)
            dhmax = min(a.hlossval(a.qmax), hvar[i, t].ub - hvar[j, t].lb)
            
            milp.addConstr(dhvar[(i, j), t] <= dhmax * svar[(i, j), t], name=f'dhxup({i},{j},{t})')
            milp.addConstr(dhvar[(i, j), t] >= dhmin * svar[(i, j), t], name=f'dhxlo({i},{j},{t})')
            if (oa_types == 'oa_cuts' or 'partial_SOS'):
                milp.addConstr(dhvar[(i, j), t] <= hvar[i, t] - hvar[j, t] - C[(i, j), t][0] * (1-svar[(i, j), t]), name=f'dhhub({i},{j},{t})')
                milp.addConstr(dhvar[(i, j), t] >= hvar[i, t] - hvar[j, t] - C[(i, j), t][1] * (1-svar[(i, j), t]), name=f'dhhlo({i},{j},{t})')
            else:
                pass
        else:
            qvar[(i, j), t] = milp.addVar(lb=Z[(i, j), t][0], ub=Z[(i, j), t][1], name=f'q({i},{j},{t})')
            dhvar[(i, j), t] = milp.addVar(lb=a.hlossval(a.qmin), ub=a.hlossval(a.qmax), name=f'H({i},{j},{t})')
            svar[(i, j), t] = milp.addVar(vtype=GRB.BINARY, lb=1, name=f'x({i},{j},{t})')
            milp.addConstr(dhvar[(i, j), t] == hvar[i, t] - hvar[j, t], name=f'dhh({i},{j},{t})')            
            
             
    for j, tank in inst.tanks.items():
        hvar[j, nperiods] = milp.addVar(lb=tank.head(tank.vinit), ub=tank.head(tank.vmax), name=f'ht({j},T)')

    milp.update()

    # FLOW CONSERVATION

    for j in inst.nodes:
        qexpr[j, t] = gp.quicksum(qvar[a, t] for a in inst.inarcs(j)) \
                          - gp.quicksum(qvar[a, t] for a in inst.outarcs(j))
                          

    for j, junc in inst.junctions.items():
        milp.addConstr(gp.quicksum(qvar[a, t] for a in inst.inarcs(j))
                           - gp.quicksum(qvar[a, t] for a in inst.outarcs(j)) == junc.demand(t), name=f'fc({j},{t})')

    for j, tank in inst.tanks.items():
        milp.addConstr(hvar[j, t+1] - hvar[j, t] == inst.flowtoheight(tank) * qexpr[j, t], name=f'fc({j},{t})')
        
    for j, tank in inst.tanks.items():
        milp.addConstr(qexpr[j, t] <= tank.qtmax[t])
        milp.addConstr(qexpr[j, t] >= tank.qtmin[t])
        

    # CONVEXIFICATION OF HEAD-FLOW
    for (i, j), arc in inst.arcs.items():


        x = svar[(i, j), t] if arc.control else 1
        if oa_types == 'full_SOS':
            f_SOS(milp, arc, i, j, t, x, dhvar, hvar, qvar, Z, oagap, accuracy, envelop, True)
        elif oa_types == 'partial_SOS':
            f_SOS(milp, arc, i, j, t, x, dhvar, hvar, qvar, Z, oagap, accuracy, envelop, False)
        else:
            cutbelow, cutabove = oa.hlossoa(Z[(i, j), t][0], Z[(i, j), t][1]+0.000001, arc.hloss, (i, j), oagap, drawgraph=False)
            for n, c in enumerate(cutbelow):
                milp.addConstr(dhvar[(i, j), t] >= c[1] * qvar[(i, j), t] + c[0] * x, name=f'hpl{n}({i},{j},{t})')
            for n, c in enumerate(cutabove):
                milp.addConstr(dhvar[(i, j), t] <= c[1] * qvar[(i, j), t] + c[0] * x, name=f'hpu{n}({i},{j},{t})')
                
                
    # CONVEXIFICATION OF HEAD-FLOW
#    for (i, j), arc in inst.arcs.items():


#        x = svar[(i, j), t] if arc.control else 1
#        if oa_types == 'full_SOS':
#            f_SOS(milp, arc, i, j, t, x, dhvar, hvar, qvar, Z, oagap, accuracy, envelop, True)
#        elif oa_types == 'partial_SOS':
#            f_SOS(milp, arc, i, j, t, x, dhvar, hvar, qvar, Z, oagap, accuracy, envelop, False)
#        else:
#            lower = arc.qmin_if_on[t] if arc.control else arc.qmin_t[t]
#            upper = arc.qmax_if_on[t] if arc.control else arc.qmax_t[t]
#            cutbelow, cutabove = oa.hlossoa(lower, upper, arc.hloss, (i, j), oagap, drawgraph=False)
#            for n, c in enumerate(cutbelow):
#                milp.addConstr(dhvar[(i, j), t] >= c[1] * qvar[(i, j), t] + c[0] * x, name=f'hpl{n}({i},{j},{t})')
#            for n, c in enumerate(cutabove):
#                milp.addConstr(dhvar[(i, j), t] <= c[1] * qvar[(i, j), t] + c[0] * x, name=f'hpu{n}({i},{j},{t})')



    binarydependencies(inst, t, milp, svar, nperiods, horizon)
    
    

        
    
    if m_ in inst.arcs.keys():
        if mode_BT =='MIN':
            obj = qvar[m_, t]
            milp.setObjective(obj, GRB.MINIMIZE)
        elif mode_BT == 'MAX':
            obj = qvar[m_, t]
            milp.setObjective(obj, GRB.MAXIMIZE)

    elif m_ in inst.tanks.keys():
        if mode_BT == 'MIN':
            obj = qexpr[m_, t]
            milp.setObjective(obj, GRB.MINIMIZE)
        elif mode_BT == 'MAX':
            obj = qexpr[m_, t]
            milp.setObjective(obj, GRB.MAXIMIZE)    
    ### probably means the same thing
    for ctan, ctanks in inst.tanks_couples.items():
        if t>=1:
            hvar[(ctan, t)] = milp.addVar(lb= ctanks.dhmin[t], ub= ctanks.dhmax[t])
            milp.addConstr(hvar[(ctan, t)] == hvar[(ctan[0], t)] - hvar[(ctan[1], t)])
        
    # you should have a look if it is fully compatible with the above constraints
#    if t >= 1:
#        milp.addConstr(hvar['t5', t] - hvar['t6', t] <= D['T56', t][1] , name=f'fc({j},{t})')
#        milp.addConstr(hvar['t5', t] - hvar['t6', t] >= D['T56', t][0] , name=f'fc({j},{t})')
        

    
    if Discrete == False:
    #enforcing configurations pumps and all the arcs relevant for the probing(conditional bounds)
        setting_configuration(inst, milp, svar, config, t)
    elif Discrete == True:

        discretization_varc_bounds(inst, milp, hvar, dh, t)

    milp.update()

    milp._svar = svar
    milp._qvar = qvar
    milp._hvar = hvar
    milp._obj = obj


    return milp


def binarydependencies(inst, t, milp, svar, nperiods, horizon):
    # PUMP SWITCHING
    sympumps = inst.symmetries
    uniquepumps = inst.pumps_without_sym()
    print('symmetries:', uniquepumps)

    def getv(vdict, pump, t):
        return gp.quicksum(vdict[a, t] for a in sympumps) if pump == 'sym' else vdict[pump, t]


    # PUMP DEPENDENCIES
    if sympumps:

        for i, pump in enumerate(sympumps[:-1]):
#                milp.addConstr(ivar[pump, t] >= ivar[sympumps[i + 1], t], name=f'symi({t})')
            milp.addConstr(svar[pump, t] >= svar[sympumps[i + 1], t], name=f'symx({t})')

    if inst.dependencies:
        for s in inst.dependencies['p1 => p0']:
            milp.addConstr(svar[s[0], t] >= svar[s[1], t], name=f'dep1({t})')
        for s in inst.dependencies['p0 xor p1']:
#            milp.addConstr(svar[s[0], t] + svar[s[1], t] >= 1, name=f'dep2({t})')
            milp.addConstr(svar[s[0], t] + svar[s[1], t] <= 1, name=f'dep2({t})')
        for s in inst.dependencies['p0 = p1 xor p2']:
            milp.addConstr(svar[s[0], t] == svar[s[1], t] + svar[s[2], t], name=f'dep3({t})')
            # for s in inst.dependencies['p1 => not p0']:
            #    milp.addConstr(svar[s[0], t] + svar[s[1], t] <= 1, name=f'dep4({t})')



def f_SOS(milp, arc, i, j, t, x, dhvar, hvar, qvar, Z, oagap, accuracy, envelop, full_SOS):

    if arc.control:
        if full_SOS == True:              
            step= sqrt(accuracy /(arc.hloss[2]+abs(arc.hloss[1])))
                    
            numbe= math.floor((Z[(i, j), t][1]-arc.qmin)/step)+1
            if numbe <=5:
                numbe=numbe+3
            else:
                pass
                
            x_samples = np.linspace(arc.qmin, Z[(i, j), t][1], numbe)
            y_samples1 = beta1(x_samples, arc.hloss)+envelop
            y_samples2 = beta1(x_samples, arc.hloss)-envelop


# 1) Instantiate a new model
            dhmax=arc.dhmax
            dhmin=arc.dhmin

            x_ = milp.addVar(lb=arc.qmin, ub=Z[(i, j), t][1], vtype=GRB.CONTINUOUS)
            y = milp.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS)
            weights = milp.addVars(len(x_samples), lb=0, ub=1, vtype=GRB.CONTINUOUS)


            milp.addSOS(GRB.SOS_TYPE2, weights)

            milp.addConstr(gp.quicksum(weights) == 1)
            milp.addConstr(gp.quicksum([weights[i]*x_samples[i] for i in range(len(x_samples))]) == x_)
            milp.addConstr(gp.quicksum([weights[i]*y_samples1[i] for i in range(len(y_samples1))]) >= y)
            milp.addConstr(gp.quicksum([weights[i]*y_samples2[i] for i in range(len(y_samples2))]) <= y)
                
            milp.addConstr(x_ == qvar[(i, j), t])
            milp.addConstr( (hvar[j, t] - hvar[i, t])-y <= (-dhmin+arc.hloss[0])*(1-x))
            milp.addConstr( (hvar[j, t] - hvar[i, t])-y >= (-dhmax+arc.hloss[0])*(1-x))
        elif full_SOS == False:
            cutbelow, cutabove = oa.hlossoa(Z[(i, j), t][0], Z[(i, j), t][1]+0.00001, arc.hloss, (i, j), oagap, drawgraph=False)
#%#        print(f'{arc}: {len(cutbelow)} cutbelow, {len(cutabove)} cutabove')
            for n, c in enumerate(cutbelow):
                milp.addConstr(dhvar[(i, j), t] >= c[1] * qvar[(i, j), t] + c[0] * x, name=f'hpl{n}({i},{j},{t})')
            for n, c in enumerate(cutabove):
                milp.addConstr(dhvar[(i, j), t] <= c[1] * qvar[(i, j), t] + c[0] * x, name=f'hpu{n}({i},{j},{t})')
                            
    else:
            
#            step= 2* sqrt(oagap/arc.hloss[2])
        step= sqrt(accuracy /(arc.hloss[2]+0.01*abs(arc.hloss[1])))
        numbe= math.floor((Z[(i, j), t][1]-Z[(i, j), t][0])/step)+1
        if numbe <=3:
            numbe=numbe+2
        else:
            pass
        x_samples = np.linspace(Z[(i, j), t][0], Z[(i, j), t][1], numbe)
        y_samples1 = beta(x_samples, arc.hloss)+envelop
        y_samples2 = beta(x_samples, arc.hloss)-envelop



        x_ = milp.addVar(lb=Z[(i, j), t][0], ub=Z[(i, j), t][1], vtype=GRB.CONTINUOUS)
        y = milp.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS)
        weights = milp.addVars(len(x_samples), lb=0, ub=1, vtype=GRB.CONTINUOUS)


        milp.addSOS(GRB.SOS_TYPE2, weights)

        milp.addConstr(gp.quicksum(weights) == 1)
        milp.addConstr(gp.quicksum([weights[i]*x_samples[i] for i in range(len(x_samples))]) == x_)
        milp.addConstr(gp.quicksum([weights[i]*y_samples1[i] for i in range(len(y_samples1))]) >= y)
        milp.addConstr(gp.quicksum([weights[i]*y_samples2[i] for i in range(len(y_samples2))]) <= y)
                
        milp.addConstr(x_ == qvar[(i, j), t])
        milp.addConstr(y == (hvar[i, t] - hvar[j, t]))
        
        
def setting_configuration(instance, milp, svar, config, t):
    for a in instance.configurations_probing[config]['ON']:
        milp.addConstr(svar[a, t] == 1)
    for a in instance.configurations_probing[config]['OFF']:
        milp.addConstr(svar[a, t] == 0)
        
        
def discretization_varc_bounds(inst, milp, hvar, dh, t):
        
        
            if t >= 1:
                for i in inst.arcs_discretizing.keys():
                    if i in dh:
                    
                    
                    
                        milp.addConstr(hvar[i[0], t] - hvar[i[1], t] <= dh[i][1] )
                        milp.addConstr(hvar[i[0], t] - hvar[i[1], t] >= dh[i][0] )
        


