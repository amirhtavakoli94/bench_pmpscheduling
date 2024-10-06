# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 19:35:35 2022

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
#h_max_milp= BT_h.build_model_BT_h(instance, Z, C, D, P0, P1, K, t_, 'MILP', 'oa_cuts', accuracy=0.01, envelop=0.1, Minim= False, oagap=OA_GAP, arcvals=None)
def build_model_BT_h(inst: Instance, Z, C, D, P0, P1, k, t_, mode_h_BT , oa_types, accuracy, envelop, Minim, two_h ,oagap: float, arcvals=None):
    """Build the convex relaxation gurobi model."""

    milp = gp.Model('Pumping_Scheduling')
    milp.params.TimeLimit= 80
    
    milp.params.Presolve = 1


    qvar = {}  # arc flow
    dhvar = {}  # arc hloss
    svar = {}  # arc status
    ivar = {}  # pump ignition status
    hvar = {}  # node head
    qexpr = {}  # node inflow
    
    
    svar_new = {} #lifted variable
    hdiff_var = {}
    q_dem = {}
    sq_var = {}


    nperiods = inst.nperiods()
    horizon = inst.horizon()
    if t_>=4 and t_<= nperiods-3:
        TT= range(t_-4, t_+3)
    elif t_<= 3:
        TT= range(0, t_+4)
    elif nperiods-4 <= t_:
        TT= range(t_-7, nperiods)
    
        
        

    for t in TT:
        for (i, j), pump in inst.pumps.items():
            ivar[(i, j), t] = milp.addVar(vtype=GRB.BINARY, name=f'ik({i},{j},{t})')

        for j in inst.junctions:
            hvar[j, t] = milp.addVar(name=f'hj({j},{t})')

        for j, res in inst.reservoirs.items():
            hvar[j, t] = milp.addVar(lb=res.head(t), ub=res.head(t), name=f'hr({j},{t})')

        for j, tank in inst.tanks.items():
            lbt = tank.head(tank.vinit) if t == 0 else D[j, t][0]
            ubt = tank.head(tank.vinit) if t == 0 else D[j, t][1]
            hvar[j, t] = milp.addVar(lb=lbt, ub=ubt, name=f'ht({j},{t})')
            if t < nperiods-1:
                
                hvar[j, t+1]= milp.addVar(lb=D[j, t+1][0], ub=D[j, t+1][1], name=f'ht({j},{t})')
        milp.update()

        for (i, j), a in inst.arcs.items():

            if a.control:
                qvar[(i, j), t] = milp.addVar(lb=-GRB.INFINITY, name=f'q({i},{j},{t})')
                dhvar[(i, j), t] = milp.addVar(lb=-GRB.INFINITY, name=f'H({i},{j},{t})')
                if ((i, j), t) in P0:
                    svar[(i, j), t] = 0
                elif ((i, j), t) in P1:
                    svar[(i, j), t] = 1
                else:
                    if mode_h_BT == 'MILP':
                        svar[(i, j), t] = milp.addVar(vtype=GRB.BINARY, name=f'x({i},{j},{t})')
                    elif mode_h_BT == 'LP':
                        svar[(i, j), t] = milp.addVar(lb=0, ub=1, name=f'x({i},{j},{t})')

                # q_a=0 if x_a=0 otherwise in [qmin,qmax]
                milp.addConstr(qvar[(i, j), t] <= Z[(i, j), t][1] * svar[(i, j), t], name=f'qxup({i},{j},{t})')
                milp.addConstr(qvar[(i, j), t] >= Z[(i, j), t][0] * svar[(i, j), t], name=f'qxlo({i},{j},{t})')
                # dh_a = (h_i - h_j) * x_a
                dhmin = max(a.hlossval(a.qmin), hvar[i, t].lb - hvar[j, t].ub)
                dhmax = min(a.hlossval(a.qmax), hvar[i, t].ub - hvar[j, t].lb)
                milp.addConstr(dhvar[(i, j), t] <= dhmax * svar[(i, j), t], name=f'dhxup({i},{j},{t})')
                milp.addConstr(dhvar[(i, j), t] >= dhmin * svar[(i, j), t], name=f'dhxlo({i},{j},{t})')
#                milp.addConstr(dhvar[(i, j), t] <= hvar[i, t] - hvar[j, t] - a.dhmin * (1-svar[(i, j), t]), name=f'dhhub({i},{j},{t})')
                milp.addConstr(dhvar[(i, j), t] <= hvar[i, t] - hvar[j, t] - C[(i, j), t][0] * (1-svar[(i, j), t]), name=f'dhhub({i},{j},{t})')
#                milp.addConstr(dhvar[(i, j), t] >= hvar[i, t] - hvar[j, t] - a.dhmax * (1-svar[(i, j), t]), name=f'dhhlo({i},{j},{t})')
                milp.addConstr(dhvar[(i, j), t] >= hvar[i, t] - hvar[j, t] - C[(i, j), t][1] * (1-svar[(i, j), t]), name=f'dhhlo({i},{j},{t})')

            else:
                qvar[(i, j), t] = milp.addVar(lb=Z[(i, j), t][0], ub=Z[(i, j), t][1], name=f'q({i},{j},{t})')
                dhvar[(i, j), t] = milp.addVar(lb=a.hlossval(a.qmin), ub=a.hlossval(Z[(i, j), t][1]), name=f'H({i},{j},{t})')
                svar[(i, j), t] = milp.addVar(vtype=GRB.BINARY, lb=1, name=f'x({i},{j},{t})')
                milp.addConstr(dhvar[(i, j), t] == hvar[i, t] - hvar[j, t], name=f'dhh({i},{j},{t})')

    for j, tank in inst.tanks.items():
        hvar[j, nperiods] = milp.addVar(lb=tank.head(tank.vinit), ub=tank.head(tank.vmax), name=f'ht({j},T)')

    milp.update()

    # FLOW CONSERVATION
    for t in TT:
        for j in inst.nodes:
            qexpr[j, t] = gp.quicksum(qvar[a, t] for a in inst.inarcs(j)) \
                          - gp.quicksum(qvar[a, t] for a in inst.outarcs(j))

        for j, junc in inst.junctions.items():
            milp.addConstr(gp.quicksum(qvar[a, t] for a in inst.inarcs(j))
                           - gp.quicksum(qvar[a, t] for a in inst.outarcs(j)) == junc.demand(t), name=f'fc({j},{t})')

        for j, tank in inst.tanks.items():
            milp.addConstr(hvar[j, t+1] - hvar[j, t] == inst.flowtoheight(tank) * qexpr[j, t], name=f'fc({j},{t})')

    # MAX WITHDRAWAL AT RESERVOIRS
    for j, res in inst.reservoirs.items():
        if res.drawmax:
            milp.addConstr(res.drawmax >=
                           inst.flowtovolume() * gp.quicksum(qexpr[j, t] for t in TT), name=f'w({j})')

    # CONVEXIFICATION OF HEAD-FLOW
    for (i, j), arc in inst.arcs.items():
        for t in TT:
            x = svar[(i, j), t] if arc.control else 1
            if oa_types == 'full_SOS':
                f_SOS(milp, arc, i, j, t, x, dhvar, hvar, qvar, Z, oagap, accuracy, envelop, True)
            elif oa_types == 'partial_SOS':
                f_SOS(milp, arc, i, j, t, x, dhvar, hvar, qvar, Z, oagap, accuracy, envelop, False)
            else:
                cutbelow, cutabove = oa.hlossoa(Z[(i, j), t][0], Z[(i, j), t][1]+0.000001, arc.hloss, (i, j), oagap, drawgraph=False)
#%#        print(f'{arc}: {len(cutbelow)} cutbelow, {len(cutabove)} cutabove')
                for n, c in enumerate(cutbelow):
                    milp.addConstr(dhvar[(i, j), t] >= c[1] * qvar[(i, j), t] + c[0] * x, name=f'hpl{n}({i},{j},{t})')
                for n, c in enumerate(cutabove):
                    milp.addConstr(dhvar[(i, j), t] <= c[1] * qvar[(i, j), t] + c[0] * x, name=f'hpu{n}({i},{j},{t})')
                    


    for ctan, ctanks in inst.tanks_couples.items():
        for t in TT:
            if t >= 1:
                hvar[(ctan, t)] = milp.addVar(lb= ctanks.dhmin[t], ub= ctanks.dhmax[t])
                milp.addConstr(hvar[ctan[0], t] - hvar[ctan[1], t] <= ctanks.dhmax[t])
                milp.addConstr(hvar[ctan[0], t] - hvar[ctan[1], t] >= ctanks.dhmin[t])




    # to have svar related to the different cnfigurations
    lifting_svar(inst, milp, TT, svar)

###    strongdualityconstraints(inst, Z, C, D, milp, hvar, qvar, svar, dhvar, qexpr, horizon, True)

    binarydependencies(inst, milp, ivar, svar, nperiods, TT)

    ###### CONDITIONAL BOUNDS BASED ON ACTIVITIES OF CONTROL ARCS - FIRST PARTITION OF VAN ZYL NETWORK #####
    probing_confs(inst, milp, svar, qvar, hvar, dhvar, qexpr, TT, oagap)

#    if hasattr(inst, '_number_disc'):
#        dh = inst.disc_tank_couples(inst._number_disc)
#        discretization_vars_consts(inst, milp, svar, qvar, hvar, dhvar, dh, inst._number_disc, horizon, oagap)

    if getattr(inst, '_number_disc', None) is not None:
        
#@        dh = inst.disc_tank_couples(inst._number_disc)
        dh = inst.disc_tank_couples
        discretization_vars_consts(inst, milp, svar, qvar, hvar, dhvar, dh, inst._number_disc, TT, oagap)
#        print("The attribute exists and is not None.")
#    else:
    # This will execute if the attribute is either None or does not exist
#        print("The attribute is None or does not exist.")
        
        
    if arcvals:
        postbinarysolution(inst, arcvals, horizon, svar)
        
    

    if two_h == False:
        if type(k) is tuple:
            obj = hvar[k[0], t_] - hvar[k[1], t_]
        else:
            obj = hvar[k, t_]
    else:
        obj = hvar[k, t_]+ hvar[k, t_+1]

    if Minim == True:

        milp.setObjective(obj, GRB.MINIMIZE)
    else:
        milp.setObjective(obj, GRB.MAXIMIZE)
    milp.update()

    milp._svar = svar
    milp._ivar = ivar
    milp._qvar = qvar
    milp._hvar = hvar
    milp._Z= Z
    milp._obj = obj

    return milp


def binarydependencies(inst, milp, ivar, svar, nperiods, TT):
    # PUMP SWITCHING
    sympumps = inst.symmetries
    uniquepumps = inst.pumps_without_sym()
####    print('symmetries:', uniquepumps)

    def getv(vdict, pump, t):
        return gp.quicksum(vdict[a, t] for a in sympumps) if pump == 'sym' else vdict[pump, t]

    # !!! check the max ignition constraint for the symmetric group
    # !!! make ivar[a,0] = svar[a,0]


    # PUMP DEPENDENCIES


    if inst.dependencies:
        for t in TT:
            for s in inst.dependencies['p1 => p0']:
                milp.addConstr(svar[s[0], t] >= svar[s[1], t], name=f'dep1({t})')
            for s in inst.dependencies['p0 xor p1']:
#                milp.addConstr(svar[s[0], t] + svar[s[1], t] >= 1, name=f'dep2({t})')
                milp.addConstr(svar[s[0], t] + svar[s[1], t] <= 1, name=f'dep2({t})')
            for s in inst.dependencies['p0 = p1 xor p2']:
                milp.addConstr(svar[s[0], t] == svar[s[1], t] + svar[s[2], t], name=f'dep3({t})')
            # for s in inst.dependencies['p1 => not p0']:
            #    milp.addConstr(svar[s[0], t] + svar[s[1], t] <= 1, name=f'dep4({t})')


def postbinarysolution(inst, arcvals, TT, svar):
    assert arcvals
    for a in inst.varcs:
        for t in TT:
            v = arcvals.get((a, t))
            if v == 1:
                svar[a, t].lb = 1
            elif v == 0:
                svar[a, t].ub = 0


def postsolution(model, vals, precision=1e-6):
    i = 0
    for var in model._svar:
        var.lb = vals[i]
        var.ub = vals[i]
        i += 1
    for var in model._qvar:
        var.lb = vals[i] - precision
        var.ub = vals[i] + precision
        i += 1
#    for var in model._hvar:
#        var.lb = vals[i] - precision
#        var.ub = vals[i] + precision
#        i += 1



def f_SOS(milp, arc, i, j, t, x, dhvar, hvar, qvar, Z, oagap, accuracy, envelop, full_SOS):

#%#        print(f'{arc}: {len(cutbelow)} cutbelow, {len(cutabove)} cutabove')
#        for t in horizon:


    if arc.control:
        if full_SOS == True:              
            step= 0.5*sqrt(accuracy /(arc.hloss[2]+abs(arc.hloss[1])))
                    
            numbe= math.floor((arc.qmax-arc.qmin)/step)+1
            if numbe <=2:
                numbe=numbe+3
            else:
                pass
                
            x_samples = np.linspace(arc.qmin, arc.qmax, numbe)
            y_samples1 = beta1(x_samples, arc.hloss)+envelop
            y_samples2 = beta1(x_samples, arc.hloss)-envelop


# 1) Instantiate a new model
            dhmax=arc.dhmax
            dhmin=arc.dhmin

            x_ = milp.addVar(lb=arc.qmin, ub=arc.qmax, vtype=GRB.CONTINUOUS)
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
            cutbelow, cutabove = oa.hlossoa(Z[(i, j), t][0], Z[(i, j), t][1], arc.hloss, (i, j), oagap, drawgraph=False)
#%#        print(f'{arc}: {len(cutbelow)} cutbelow, {len(cutabove)} cutabove')
            for n, c in enumerate(cutbelow):
                milp.addConstr(dhvar[(i, j), t] >= c[1] * qvar[(i, j), t] + c[0] * x, name=f'hpl{n}({i},{j},{t})')
            for n, c in enumerate(cutabove):
                milp.addConstr(dhvar[(i, j), t] <= c[1] * qvar[(i, j), t] + c[0] * x, name=f'hpu{n}({i},{j},{t})')
                            





    else:
            
#            step= 2* sqrt(oagap/arc.hloss[2])
        step= 0.5*sqrt(accuracy /arc.hloss[2])
        numbe= math.floor((arc.qmax-arc.qmin)/step)+1
        if numbe <=2:
            numbe=numbe+2
        else:
            pass
        x_samples = np.linspace(arc.qmin, arc.qmax, numbe)
        y_samples1 = beta(x_samples, arc.hloss)+envelop
        y_samples2 = beta(x_samples, arc.hloss)-envelop



        x_ = milp.addVar(lb=arc.qmin, ub=arc.qmax, vtype=GRB.CONTINUOUS)
        y = milp.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS)
        weights = milp.addVars(len(x_samples), lb=0, ub=1, vtype=GRB.CONTINUOUS)


        milp.addSOS(GRB.SOS_TYPE2, weights)

        milp.addConstr(gp.quicksum(weights) == 1)
        milp.addConstr(gp.quicksum([weights[i]*x_samples[i] for i in range(len(x_samples))]) == x_)
        milp.addConstr(gp.quicksum([weights[i]*y_samples1[i] for i in range(len(y_samples1))]) >= y)
        milp.addConstr(gp.quicksum([weights[i]*y_samples2[i] for i in range(len(y_samples2))]) <= y)
                
        milp.addConstr(x_ == qvar[(i, j), t])
        milp.addConstr(y == (hvar[i, t] - hvar[j, t]))
        

def lifting_svar(instance, milp, horizon, svar):
    for conf, config_val in instance.configurations_probing.items():
        for t in horizon:
            svar[conf, t] = milp.addVar(lb=0, ub=1, vtype = GRB.BINARY)
            for K in config_val['ON']:
                milp.addConstr(svar[conf, t] <= svar[K, t])
            for K in config_val['OFF']:
                milp.addConstr(svar[conf, t] <= 1 - svar[K, t])
            milp.addConstr(svar[conf, t] >= gp.quicksum(svar[K, t] for K in config_val['ON'])+ gp.quicksum(1-svar[K, t] for K in config_val['OFF']) - len(config_val['ON'])- len(config_val['OFF']) - 1)
            
            
            
        
def probing_confs(inst, milp, svar, qvar, hvar, dhvar, qexpr, horizon, oagap):
    
    for conf_arc_keys, conf_arcs_values in inst.arcs_probing.items():
        for (i, j), arc in inst.arcs.items():
            if conf_arc_keys[1] == (i, j):
                for t in horizon:
#@                    dhvar[(i, j), t, conf_arc_keys[0]] = milp.addVar(lb=-abs(arc.hlossval(arc.qmin)), ub=abs(arc.hlossval(arc.qmax))) 
                    
                    dhvar[(i, j), t, conf_arc_keys[0]] = milp.addVar(lb=-10000, ub=10000) 
                    
    
    for conf_arc_keys, conf_arcs_values in inst.arcs_probing.items():
        for (i, j), arc in inst.arcs.items():
            if conf_arc_keys[1] == (i, j):
                for t in horizon:
                    qvar[(i, j), t, conf_arc_keys[0]] = milp.addVar(lb= arc.qmin, ub= conf_arcs_values.qmax_prob[t]) 
                    
                    milp.addConstr(qvar[(i, j), t, conf_arc_keys[0]] <= conf_arcs_values.qmax_prob[t] * svar[conf_arc_keys[0], t])
                    milp.addConstr(qvar[(i, j), t, conf_arc_keys[0]] >= conf_arcs_values.qmin_prob[t] * svar[conf_arc_keys[0], t])
                    
                
                    
    
    for conf_arc_keys, conf_arcs_values in inst.arcs_probing.items():
        for (i, j), arc in inst.arcs.items():
            if conf_arc_keys[1] == (i, j):
                for t in horizon:


  #                  cutbelow, cutabove = oa.hlossoa(max(conf_arcs_values.qmin_prob[t], temp_min), min(conf_arcs_values.qmax_prob[t]+10e-4, temp_max+10e-4), arc.hloss, (i, j), oagap, drawgraph=False)
                    if arc.control:
                        cutbelow, cutabove = oa.hlossoa(conf_arcs_values.qmin_prob[t], conf_arcs_values.qmax_prob[t]+10e-4, arc.hloss, (i, j), oagap, drawgraph=False)
                        
                        
                        for n, c in enumerate(cutbelow):
                            milp.addConstr(dhvar[(i, j), t, conf_arc_keys[0]]   >= c[1] * qvar[(i, j), t, conf_arc_keys[0]] + c[0] * svar[conf_arc_keys[0], t], name=f'hpl{n}({i},{j},{t})')
                        for n, c in enumerate(cutabove):
                            milp.addConstr(dhvar[(i, j), t, conf_arc_keys[0]]  <= c[1] * qvar[(i, j), t, conf_arc_keys[0]] + c[0] * svar[conf_arc_keys[0], t], name=f'hpu{n}({i},{j},{t})')
                    
                    
                    else:
                        cutbelow, cutabove = oa.hlossoa(conf_arcs_values.qmin_prob[t], conf_arcs_values.qmax_prob[t]+10e-4, arc.hloss, (i, j), oagap, drawgraph=False)
                    
#@                    print(f'{arc}: {len(cutbelow)} cutbelow, {len(cutabove)} cutabove')
                    for n, c in enumerate(cutbelow):
                        milp.addConstr(dhvar[(i, j), t, conf_arc_keys[0]]   >= c[1] * qvar[(i, j), t, conf_arc_keys[0]] + c[0] * svar[conf_arc_keys[0], t], name=f'hpl{n}({i},{j},{t})')
                    for n, c in enumerate(cutabove):
                        milp.addConstr(dhvar[(i, j), t, conf_arc_keys[0]]  <= c[1] * qvar[(i, j), t, conf_arc_keys[0]] + c[0] * svar[conf_arc_keys[0], t], name=f'hpu{n}({i},{j},{t})')



    unique_keys_ij = set(conf_arc_keys[1] for conf_arc_keys, _ in inst.arcs_probing.items())
    unique_keys_conf = set(conf_arc_keys[0] for conf_arc_keys, _ in inst.arcs_probing.items())
    
    
    new_dict_arcs_probing = {}

    for (conf, arc), times in inst.arcs_probing.items():
        for t in inst.horizon():
            new_key = (arc, t)  # Combine arc tuple and time into a new tuple key
            if new_key not in new_dict_arcs_probing:
                new_dict_arcs_probing[new_key] = {}
            new_dict_arcs_probing[new_key][conf] = [ inst.arcs_probing[conf, arc].qmin_prob[t],  inst.arcs_probing[conf, arc].qmax_prob[t]]
            
            
    new_dict_mein = {}

    for (conf, arc), times in inst.arcs_probing.items():
        for t in horizon:
            new_key = (arc, t)
            if new_key not in new_dict_mein:
                new_dict_mein[new_key] = set()
            new_dict_mein[new_key].add(conf)
            

    for (arc, t_), set_confs in new_dict_mein.items():
                (i, j) = arc
                
                
                milp.addConstr(qvar[(i, j), t_] == gp.quicksum(qvar[(i, j), t_, conf_arc_keys] for conf_arc_keys in set_confs ))
                
###                print((i, j, set_confs))

                if inst.arcs[arc].control:
                    milp.addConstr(dhvar[(i, j), t_]  <= gp.quicksum(dhvar[(i, j), t_, conf] for conf in set_confs))
                    milp.addConstr(dhvar[(i, j), t_]  >= gp.quicksum(dhvar[(i, j), t_, conf] for conf in set_confs))
#                    pass
                else:
                
                
                    milp.addConstr(dhvar[(i, j), t_] <= gp.quicksum(dhvar[(i, j), t_, conf] for conf in set_confs))
                    milp.addConstr(dhvar[(i, j), t_] >= gp.quicksum(dhvar[(i, j), t_, conf] for conf in set_confs))
            

    
    ###though these consts must be eclipsed by the below consts
    for j, tank in inst.tanks.items():
        milp.addConstr(qexpr[j, t] <= tank.qtmax[t])
        milp.addConstr(qexpr[j, t] >= tank.qtmin[t])    



    for j, tank in inst.tanks.items():
        for t in horizon:
            milp.addConstr(qexpr[j, t] <= gp.quicksum(svar[conf_tank_keys[0], t] * conf_tank_values.qexprmax_prob[t] for conf_tank_keys, conf_tank_values in inst.tanks_inflow_probing.items() if conf_tank_keys[1] == j) + tank.qtmax[t]*(1-gp.quicksum(svar[conf_tank_keys[0], t] for conf_tank_keys, conf_tank_values in inst.tanks_inflow_probing.items() if conf_tank_keys[1] == j   )))
            milp.addConstr(qexpr[j, t] >= gp.quicksum(svar[conf_tank_keys[0], t] * conf_tank_values.qexprmin_prob[t] for conf_tank_keys, conf_tank_values in inst.tanks_inflow_probing.items() if conf_tank_keys[1] == j) + tank.qtmin[t]*(1-gp.quicksum(svar[conf_tank_keys[0], t] for conf_tank_keys, conf_tank_values in inst.tanks_inflow_probing.items() if conf_tank_keys[1] == j   )))



def discretization_vars_consts(inst, milp, svar, qvar, hvar, dhvar, dh, number_disc, horizon, oagap):
    
    
    for t in horizon:
        if t >= 1:
            for i in inst.arcs_discretizing.keys():
                for n in range(0, number_disc[i]):
                    svar[i, f'discrete_{n}', t] = milp.addVar(vtype = GRB.BINARY) 
                
                
                for n in range(0, number_disc[i]):
                    milp.addConstr(hvar[i[0], t] - hvar[i[1], t] <= dh[i][f'discrete_{n}', t][1] * svar[i, f'discrete_{n}', t] + inst.tanks_couples[i].dhmax[t]*(1-svar[i, f'discrete_{n}', t]))
                    milp.addConstr(hvar[i[0], t] - hvar[i[1], t] >= dh[i][f'discrete_{n}', t][0] * svar[i, f'discrete_{n}', t] + inst.tanks_couples[i].dhmin[t]*(1-svar[i, f'discrete_{n}', t]))

    for t in horizon:
        if t >= 1:
            for i in inst.arcs_discretizing.keys():
                milp.addConstr(gp.quicksum(svar[i, f'discrete_{n}', t] for n in range(0, number_disc[i])) == 1)
                    

                    
    #all related methods should beforehand initialized 
    for t in horizon:
        if t >= 1:
            for i, arc in inst.arcs_discretizing.items():
                for ii in arc:
                    for jj in ii:
                
                        
                        for n in range(0, number_disc[i]):
                            qvar[jj, f'discrete_{n}', t] = milp.addVar(lb = inst.arcs[jj].qmin, ub = inst.arcs[jj].qmax)
                            milp.addConstr(qvar[jj, f'discrete_{n}', t] >= svar[i, f'discrete_{n}', t] * inst.arcs_discretization_initialization[(i, jj)].qtmin_disc[t, f'discrete_{n}'])
                            milp.addConstr(qvar[jj, f'discrete_{n}', t] <= svar[i, f'discrete_{n}', t] * inst.arcs_discretization_initialization[(i, jj)].qtmax_disc[t, f'discrete_{n}'])
                
                            dhvar[jj, f'discrete_{n}', t] = milp.addVar(lb=-abs(inst.arcs[jj].hlossval(inst.arcs[jj].qmin)), ub=abs(inst.arcs[jj].hlossval(inst.arcs[jj].qmax)))
                            
                        milp.addConstr(gp.quicksum(qvar[jj, f'discrete_{n}', t] for n in range(0, number_disc[i])) == qvar[jj, t])
                
    
                        

    #all related methods should beforehand initialized 
    for t in horizon:
        if t >= 1:
            for i, arcs in inst.arcs_discretizing.items():
                for ii in arcs:
                    for jj in ii:
                        
                        for n_ in range(0, number_disc[i]):
                    
                            cutbelow, cutabove = oa.hlossoa(inst.arcs_discretization_initialization[(i, jj)].qtmin_disc[t, f'discrete_{n_}'], inst.arcs_discretization_initialization[(i, jj)].qtmax_disc[t, f'discrete_{n_}'], inst.arcs[jj].hloss, jj, oagap, drawgraph=False)
 #@                           print(f'{arc}: {len(cutbelow)} cutbelow, {len(cutabove)} cutabove')
                            for n, c in enumerate(cutbelow):
                                milp.addConstr(dhvar[jj, f'discrete_{n_}', t]  >= c[1] * qvar[jj, f'discrete_{n_}', t] + c[0] * svar[i, f'discrete_{n_}', t])
                            for n, c in enumerate(cutabove):
                                milp.addConstr(dhvar[jj, f'discrete_{n_}', t]  <= c[1] * qvar[jj, f'discrete_{n_}', t] + c[0] * svar[i, f'discrete_{n_}', t])
                                
                        milp.addConstr(gp.quicksum(dhvar[jj, f'discrete_{n}', t] for n in range(0, number_disc[i]))   == dhvar[jj, t])
                                    
                
    for t in horizon:
        if t >= 1:
            for i in inst.arcs_discretizing.keys():
                milp.addConstr(gp.quicksum(svar[i, f'discrete_{n}', t] for n in range(0, number_disc[i])) == 1)
                
#2    for t in horizon:
#@        if t >= 1:
#@            for i, arcs in inst.arcs_discretizing.items():
#@                for ii in arcs:
#@                    for jj in ii:
#                        dhvar[jj, t, 'all_discrete'] = milp.addVar(lb = -2000, ub = 2000)
#2                        temp = 0
#2                        for n in range(0, number_disc[i]):
#2                            temp += dhvar[jj, f'discrete_{n}', t]     
                        
#2                        milp.addConstr(dhvar[jj, t] == temp)
                        
#@                        milp.addConstr(gp.quicksum(dhvar[jj, f'discrete_{n}', t] for n in range(0, number_disc[i])) - 0.01 <= dhvar[jj, t])
#                        milp.addConstr(gp.quicksum(dhvar[jj, f'discrete_{n}', t] for n in range(0, number_disc[i])) + 0.01 >= dhvar[jj, t])
    
#@    for t in horizon:
#@        if t >= 1:
#@            for i, arcs in inst.arcs_discretizing.items():
#@                for ii in arcs:
#@                    for jj in ii:
#@                        milp.addConstr(gp.quicksum(dhvar[jj, f'discrete_{n}', t] for n in range(0, number_disc[i])) - 0.01 <= dhvar[jj, t])
#@                        milp.addConstr(gp.quicksum(dhvar[jj, f'discrete_{n}', t] for n in range(0, number_disc[i])) + 0.01 >= dhvar[jj, t])