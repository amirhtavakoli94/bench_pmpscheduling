# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 18:53:11 2022

@author: amirhossein.tavakoli
"""

import gurobipy as gp
from gurobipy import GRB
import outerapproximation as oa
from instance import Instance
import datetime as dt
import math


# !!! check round values
# noinspection PyArgumentList
def build_model(inst: Instance, Z, C, D, P0, P1, sn_b_all, oagap: float, arcvals=None):
    """Build the convex relaxation gurobi model."""

    milp = gp.Model('Pumping_Scheduling')

    qvar = {}  # arc flow
    dhvar = {}  # arc hloss
    svar = {}  # arc status
    ivar = {}  # pump ignition status
    hvar = {}  # node head
    qexpr = {}  # node inflow
    
    lifted_svar = {}
    
    dh = inst.disc_tank_couples
    s3_var={}
    
#    svar_disc = {}
    
    
###    s3_var= {}


    nperiods = inst.nperiods()
    horizon = inst.horizon()

    for t in horizon:
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
            
        for j, tank in inst.tanks.items():
            lbt = tank.head(tank.vinit) if t == 0 else D[j, t][0]
            ubt = tank.head(tank.vinit) if t == 0 else D[j, t][1]


        milp.update()


        for (i, j), a in inst.arcs.items():

            if a.control:
                qvar[(i, j), t] = milp.addVar(lb=-GRB.INFINITY, name=f'q({i},{j},{t})')
                dhvar[(i, j), t] = milp.addVar(lb=-GRB.INFINITY, name=f'H({i},{j},{t})')
                if ((i, j), t) in P0:
                    svar[(i, j), t] = milp.addVar(lb=0, ub=0, name=f'x({i},{j},{t})')
                elif ((i, j), t) in P1:
                    svar[(i, j), t] = milp.addVar(lb=1, ub=1, name=f'x({i},{j},{t})')
                else:
                    svar[(i, j), t] = milp.addVar(vtype=GRB.BINARY, name=f'x({i},{j},{t})')
                # q_a=0 if x_a=0 otherwise in [qmin,qmax]
                milp.addConstr(qvar[(i, j), t] <= Z[(i, j), t][1] * svar[(i, j), t], name=f'qxup({i},{j},{t})')
                milp.addConstr(qvar[(i, j), t] >= Z[(i, j), t][0] * svar[(i, j), t], name=f'qxlo({i},{j},{t})')
                

                dhmin = max(a.hlossval(a.qmin), hvar[i, t].lb - hvar[j, t].ub)
                dhmax = min(a.hlossval(Z[(i, j), t][1]), hvar[i, t].ub - hvar[j, t].lb)
                
                milp.addConstr(dhvar[(i, j), t] <= dhmax * svar[(i, j), t], name=f'dhxup({i},{j},{t})')
                milp.addConstr(dhvar[(i, j), t] >= dhmin * svar[(i, j), t], name=f'dhxlo({i},{j},{t})')
                milp.addConstr(dhvar[(i, j), t] <= hvar[i, t] - hvar[j, t] - a.dhmin_t[t] * (1-svar[(i, j), t]), name=f'dhhub({i},{j},{t})')
                milp.addConstr(dhvar[(i, j), t] <= hvar[i, t] - hvar[j, t] - C[(i, j), t][0] * (1-svar[(i, j), t]), name=f'dhhub({i},{j},{t})')
                milp.addConstr(dhvar[(i, j), t] >= hvar[i, t] - hvar[j, t] - a.dhmax_t[t] * (1-svar[(i, j), t]), name=f'dhhlo({i},{j},{t})')
                milp.addConstr(dhvar[(i, j), t] >= hvar[i, t] - hvar[j, t] - C[(i, j), t][1] * (1-svar[(i, j), t]), name=f'dhhlo({i},{j},{t})')

            else:
                qvar[(i, j), t] = milp.addVar(lb=Z[(i, j), t][0], ub=Z[(i, j), t][1], name=f'q({i},{j},{t})')
                dhvar[(i, j), t] = milp.addVar(lb=a.hlossval(a.qmin), ub=a.hlossval(Z[(i, j), t][1]), name=f'H({i},{j},{t})')
                svar[(i, j), t] = milp.addVar(vtype=GRB.BINARY, lb=1, name=f'x({i},{j},{t})')
                milp.addConstr(dhvar[(i, j), t] == hvar[i, t] - hvar[j, t], name=f'dhh({i},{j},{t})')
                


    for j, tank in inst.tanks.items():
        hvar[j, nperiods] = milp.addVar(lb=tank.head(tank.vinit), ub=tank.head(tank.vmax), name=f'ht({j},T)')
        
        


    milp.update()
    
    
    # to have svar related to the different cnfigurations
    lifting_svar(inst, milp, horizon, svar, lifted_svar)
    
    milp.update()

    # FLOW CONSERVATION
    for t in horizon:
        for j in inst.nodes:
            qexpr[j, t] = gp.quicksum(qvar[a, t] for a in inst.inarcs(j)) \
                          - gp.quicksum(qvar[a, t] for a in inst.outarcs(j))

        for j, junc in inst.junctions.items():
            milp.addConstr(gp.quicksum(qvar[a, t] for a in inst.inarcs(j))
                           - gp.quicksum(qvar[a, t] for a in inst.outarcs(j)) == junc.demand(t), name=f'fc({j},{t})')

        for j, tank in inst.tanks.items():
            milp.addConstr(hvar[j, t+1] - hvar[j, t] == inst.flowtoheight(tank) * qexpr[j, t], name=f'fc({j},{t})')
            milp.addConstr(qexpr[j, t]<= Z[j, t][1]+0.01)
            milp.addConstr(qexpr[j, t]>= Z[j, t][0]-0.01)
            



    # MAX WITHDRAWAL AT RESERVOIRS
    for j, res in inst.reservoirs.items():
        if res.drawmax:
            milp.addConstr(res.drawmax >=
                           inst.flowtovolume() * gp.quicksum(qexpr[j, t] for t in horizon), name=f'w({j})')
            


    # CONVEXIFICATION OF HEAD-FLOW
    for (i, j), arc in inst.arcs.items():
        for t in horizon:
            #to do instead of Z, you can use arc.qmin_t or arc.qmax_t
            cutbelow, cutabove = oa.hlossoa(Z[(i, j), t][0], Z[(i, j), t][1]+0.00001, arc.hloss, (i, j), oagap, drawgraph=False)
###@            print(f'{arc}: {len(cutbelow)} cutbelow, {len(cutabove)} cutabove')
            x = svar[(i, j), t] if arc.control else 1
            for n, c in enumerate(cutbelow):
                milp.addConstr(dhvar[(i, j), t] >= c[1] * qvar[(i, j), t] + c[0] * x, name=f'hpl{n}({i},{j},{t})')
            for n, c in enumerate(cutabove):
                milp.addConstr(dhvar[(i, j), t] <= c[1] * qvar[(i, j), t] + c[0] * x, name=f'hpu{n}({i},{j},{t})')
                
                
    ###### CONDITIONAL BOUNDS BASED ON ACTIVITIES OF CONTROL ARCS - FIRST PARTITION OF VAN ZYL NETWORK #####
#@    probing_confs(inst, milp, svar, lifted_svar, qvar, hvar, dhvar, qexpr, horizon, oagap)
    ##### then discretization, adding constraints related to the difference of the two tanks levels
    
#@    discretization_vars_consts(inst, milp, svar, qvar, hvar, dhvar, dh, inst._number_disc, horizon, oagap)
                        
    







#    nested_dict = {}
#    for key, value in inst._tank_mir_cuts.items():
#    # Unpacking the original key
#        tank, set_arcs, ts, tf = key
#    # Creating the outer key based on ts and tf
#        outer_key = (ts, tf)
    
#    # Creating the middle layer if not exists
#        if outer_key not in nested_dict:
#            nested_dict[outer_key] = {}
    
#    # Checking if set_arcs already exists in the middle dictionary
#        if set_arcs not in nested_dict[outer_key]:
#            nested_dict[outer_key][set_arcs] = {}
    
#    # Assigning the value with tank as the most inner key
#        nested_dict[outer_key][set_arcs][tank] = value





#    hvar_ts = {}
#    hvar_tf = {}
#    for j, tank in inst.tanks.items():
#        for t in inst.horizon():
#            lbt = tank.head(tank.vinit) if t == 0 else D[j, t][0]
#            ubt = tank.head(tank.vinit) if t == 0 else D[j, t][1]
#            hvar_ts[j, t] = milp.addVar(lb=0, ub=ubt-lbt, name=f'hts({j},{t})')
#            hvar_tf[j, t] = milp.addVar(lb=0, ub=ubt-lbt, name=f'htf({j},{t})')


#    for ii, jj in nested_dict.items():
#        
#        for set_arcs, rest in jj.items():
#            
#           for ij, jk in rest.items():
#               bn = jk.x_star
#           b = bn - sum(mir_val.mu_s*(inst.tanks[tan].h_star[ii[0]] - inst.tanks[tan].hmin[ii[0]]) for tan, mir_val in rest.items()) + sum(mir_val.mu_f*(inst.tanks[tan].h_star[ii[1]] - inst.tanks[tan].hmin[ii[1]]) for tan, mir_val in rest.items()) 
#           
#           b_ = math.floor(b)
           
#           b_hat = b-b_
           
#           if b_hat > 0.0001:
               
#               milp.addConstr(s3_var[ii[0]]>= (1/b_hat)*((gp.quicksum(mir_val.mu_s*hvar_ts[j, ii[0]] for tan, mir_val in rest.items()))-gp.quicksum(mir_val.mu_f*hvar_ts[j, ii[1]] for tan, mir_val in rest.items())+b_+1))

                

#    for ii, jj in inst._tank_mir_cuts.items():
#        b = jj.x_star - sum(jj.mu_star*(inst.tanks[ii[0]].h_star[ii[2]] - inst.tanks[ii[0]].hmin[ii[2]]) )
                        






#1    for t_ in range(1, nperiods-1):
###        for j, tank in inst.tanks.items():
            
#####            b = instance
        
#1            b=S_n1[(t_,t_+1)][0]-sum(pi_1[j, t_][0]*(bt_h_p[f'ht({j},{t_})'][0]-D[j, t_][0]) for j, tank in inst.tanks.items())+sum(pi_1[j, t_][1]*(-bt_h_p[f'ht({j},{t_+1})'][0]+D[j, t_+1][1])for j, tank in inst.tanks.items())
#1            b_= math.floor(b)
#1            b_hat=b-b_
#1            if b_hat>0.0000001:
#1                    milp.addConstr(s1_var[t_]>= (1/b_hat)*((gp.quicksum(pi_1[j, t_][0]*hvar_ts[j, t_] for j, tank in inst.tanks.items()))-gp.quicksum(pi_1[j, t_][1]*hvar_tf[j, t_+1]for j, tank in inst.tanks.items())+b_+1))
    milp.update()
                    
                    


    for t in horizon:
        s3_var[t]=0
        
        
    milp.update()
    
#    aa=inst.arcs

#    for t in horizon:
#        for (i, j), arc in inst.arcs.items():
#            if arc.control:
#                if aa[(i, j)].type == 'FSD':
####                s3_var[t]= svar[(i, j), t] + s3_var[t]
#                    s3_var[t]= svar[(i, j), t] + s3_var[t]
                    
                    
    aa=inst.arcs
    
    
    # THE NUMBER PUMPS ACTIVE UNTIL CERTAIN TIME

    for t in horizon:
        for (i, j), arc in inst.arcs.items():
            if arc.control:
                if aa[(i, j)].type == 'FSD':
                    s3_var[t]= svar[(i, j), t] + s3_var[t]
                    

#    for t in horizon:
#        milp.addConstr(s_n_var[t] == s3_var[t])               

    # CARDINALITY INEQUALITY RELATED TO ACTIVITY OF THE PUMPS TO AND FROM CERTAIN TIME 
#    for t_ in range(0, nperiods):
#        milp.addConstr(gp.quicksum(s3_var[t1] for t1 in range(0, t_))>= math.ceil(sn_b_all[t_][0]))
#        milp.addConstr(gp.quicksum(s3_var[t1] for t1 in range(t_, nperiods))>= math.ceil(sn_b_all[t_][1]))
                    
                    
    if inst._tank_mir_cuts is None:
        pass
#    else:
#        if len(inst._tank_mir_cuts.items())>=1:
#            mir_cuts_implementation_n(milp, inst, D, s3_var, hvar)
                    
                    


####    strongdualityconstraints(inst, Z, C, D, milp, hvar, qvar, svar, dhvar, qexpr, horizon, True)

    binarydependencies(inst, milp, ivar, svar, nperiods, horizon)
    

    
    ###### CONDITIONAL BOUNDS BASED ON ACTIVITIES OF CONTROL ARCS - FIRST PARTITION OF VAN ZYL NETWORK #####
 #   probing_confs(inst, milp, svar, qvar, hvar, dhvar, qexpr, horizon, oagap)
    ##### then discretization, adding constraints related to the difference of the two tanks levels
    
 #   discretization_vars_consts(inst, milp, svar, qvar, hvar, dhvar, dh, inst._number_disc, horizon, oagap)
                        

    if arcvals:
        postbinarysolution(inst, arcvals, horizon, svar)

    obj = gp.quicksum(inst.eleccost(t) * (pump.power[0] * svar[k, t] + pump.power[1] * qvar[k, t])
                      for k, pump in inst.pumps.items() for t in horizon)

    milp.setObjective(obj, GRB.MINIMIZE)
    milp.update()

    milp._svar = svar
    milp._ivar = ivar
    milp._qvar = qvar
    milp._hvar = hvar
    milp._Z= Z
    milp._obj = obj

    return milp


def strongdualityconstraints(inst, Z, C, D, milp, hvar, qvar, svar, dhvar, qexpr, horizon, withoutz):
    print("#################  STRONG DUALITY: 5 gvar(pipe) + 10 gvar (pump)")
    # strong duality constraint: sum_a gvar[a,t] + sdexpr[t] <= 0
    gvar = {}    # arc component:    x_a * (\Phi_a(q_a) - \Phi_a(\phi^{-1}(h_a)) + h_a\phi^{-1}(h_a))
    sdexpr = {}  # node component:   sum_n (q_n * h_n)
    hqvar = {}   # tank component:   q_r * h_r
    for t in horizon:

        # McCormick's envelope of hq_rt = h_rt * q_rt = h_rt * (h_{r,t+1}-h_rt)/c
        for j, tank in inst.tanks.items():
            c = inst.flowtoheight(tank)
            (h0, h1) = (hvar[j, t], hvar[j, t + 1])
            if t==0:
            
                (l0, l1, u0, u1) = (h0.lb, D[j, t+1][0], h0.ub, D[j, t+1][1])
            elif t== inst.nperiods()-1:
                (l0, l1, u0, u1) = (D[j, t][0], h1.lb, D[j, t][1], h1.ub)
            else:
                (l0, l1, u0, u1) = (D[j, t][0], D[j, t+1][0], D[j, t][1], D[j, t+1][1])
            if l0 == u0:
                hqvar[j, t] = (h1 - l0) * l0 / c
            else:
                hqvar[j, t] = milp.addVar(lb=-GRB.INFINITY, name=f'hqt({j},{t})')
                inflow = {a: [inst.arcs[a].qmin, inst.arcs[a].qmax] for a in inst.inarcs(j)}
                outflow = {a: [inst.arcs[a].qmin, inst.arcs[a].qmax] for a in inst.outarcs(j)}
                print(f"inflow: {inflow}")
                print(f"outflow: {outflow}")
                lq = max(c * inst.inflowmin(j), c * Z[j, t][0], l1 - u0)
                uq = min(c * inst.inflowmax(j), c * Z[j, t][1], u1 - l0)
                # refining with a direction indicator variable
                if withoutz:
                    milp.addConstr(c * hqvar[j, t] >= l0 * (h1 - h0) + lq * (h0 - l0), name=f'hqlo({j},{t})')
                    milp.addConstr(c * hqvar[j, t] >= u0 * (h1 - h0) + uq * (h0 - u0), name=f'hqup({j},{t})')
                else:
                    zvar = milp.addVar(vtype=GRB.BINARY, name=f'z({j},{t})')
                    hzvar = milp.addVar(lb=0, ub=u0, name=f'hz({j},{t})')
                    milp.addConstr(hzvar <= u0 * zvar, name=f'hz1up({j},{t})')
                    milp.addConstr(hzvar >= l0 * zvar, name=f'hz1lo({j},{t})')
                    milp.addConstr(hzvar <= h0 - l0 * (1 - zvar), name=f'hz0up({j},{t})')
                    milp.addConstr(hzvar >= h0 - u0 * (1 - zvar), name=f'hz0lo({j},{t})')
                    milp.addConstr(c * hqvar[j, t] >= l0 * (h1 - h0) + lq * (hzvar - l0 * zvar), name=f'hqlo({j},{t})')
                    milp.addConstr(c * hqvar[j, t] >= u0 * (h1 - h0) + uq * (h0 - hzvar - u0 * (1 - zvar)), name=f'hqup({j},{t})')

        # sdexpr[t] = milp.addVar(lb=-GRB.INFINITY, name=f'sd({t})')
        sdexpr[t] = gp.quicksum(hqvar[j, t] for j in inst.tanks) \
            + gp.quicksum(junc.demand(t) * hvar[j, t] for j, junc in inst.junctions.items()) \
            + gp.quicksum(res.head(t) * qexpr[j, t] for j, res in inst.reservoirs.items())

        # OA for the convex function g_a >= Phi_a(q_a) - Phi_a(phi^{-1}(h_a)) + h_a * phi^{-1}(h_a)
        for (i,j), arc in inst.arcs.items():
            a = (i, j)
            gvar[a, t] = milp.addVar(lb=-GRB.INFINITY, name=f'g({i},{j},{t})')
            noacut = 10 if a in inst.pumps else 5
            for n in range(noacut):
                qstar = (arc.qmin + Z[(i,j), t][1]) * n / (noacut - 1)
                milp.addConstr(gvar[a, t] >= arc.hlossval(qstar) *
                               (qvar[a, t] - qstar * svar[a, t]) + qstar * dhvar[a, t], name=f'goa{n}({i},{j},{t})')

        milp.addConstr(gp.quicksum(gvar[a, t] for a in inst.arcs) + sdexpr[t] <= milp.Params.MIPGapAbs, name=f'sd({t})')


def binarydependencies(inst, milp, ivar, svar, nperiods, horizon):
    # PUMP SWITCHING
    sympumps = inst.symmetries
    uniquepumps = inst.pumps_without_sym()
    print('symmetries:', uniquepumps)

    def getv(vdict, pump, t):
        return gp.quicksum(vdict[a, t] for a in sympumps) if pump == 'sym' else vdict[pump, t]

    # !!! check the max ignition constraint for the symmetric group
    # !!! make ivar[a,0] = svar[a,0]
    for a in uniquepumps:
        rhs = 6 * len(sympumps) if a == 'sym' else 6 - svar[a, 0]
        milp.addConstr(gp.quicksum(getv(ivar, a, t) for t in range(1, nperiods)) <= rhs, name='ig')
        for t in range(1, nperiods):
            milp.addConstr(getv(ivar, a, t) >= getv(svar, a, t) - getv(svar, a, t - 1), name=f'ig({t})')
            if inst.tsduration == dt.timedelta(minutes=30) and t < inst.nperiods() - 1:
                # minimum 1 hour activity
                milp.addConstr(getv(svar, a, t + 1) + getv(svar, a, t - 1) >= getv(svar, a, t), name=f'ig1h({t})')

    # PUMP DEPENDENCIES
    if sympumps:
        for t in horizon:
            for i, pump in enumerate(sympumps[:-1]):
                milp.addConstr(ivar[pump, t] >= ivar[sympumps[i + 1], t], name=f'symi({t})')
                milp.addConstr(svar[pump, t] >= svar[sympumps[i + 1], t], name=f'symx({t})')

    if inst.dependencies:
        for t in horizon:
            for s in inst.dependencies['p1 => p0']:
                milp.addConstr(svar[s[0], t] >= svar[s[1], t], name=f'dep1({t})')
            for s in inst.dependencies['p0 xor p1']:
                milp.addConstr(svar[s[0], t] + svar[s[1], t] >= 1, name=f'dep2({t})')
            for s in inst.dependencies['p0 = p1 xor p2']:
                milp.addConstr(svar[s[0], t] == svar[s[1], t] + svar[s[2], t], name=f'dep3({t})')
            # for s in inst.dependencies['p1 => not p0']:
            #    milp.addConstr(svar[s[0], t] + svar[s[1], t] <= 1, name=f'dep4({t})')


def postbinarysolution(inst, arcvals, horizon, svar):
    assert arcvals
    for a in inst.varcs:
        for t in horizon:
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

def lifting_svar(instance, milp, horizon, svar, lifted_svar):
    

    for conf, config_val in instance.configurations_probing.items():
        for t in horizon:
            lifted_svar[conf, t] = milp.addVar(lb=0, ub=1, vtype = GRB.BINARY)
            for K in config_val['ON']:
                milp.addConstr(lifted_svar[conf, t] <= svar[K, t])
            for K in config_val['OFF']:
                milp.addConstr(lifted_svar[conf, t] <= 1 - svar[K, t])
            milp.addConstr(lifted_svar[conf, t] >= gp.quicksum(svar[K, t] for K in config_val['ON'])+ gp.quicksum(1-svar[K, t] for K in config_val['OFF']) - len(config_val['ON'])- len(config_val['OFF']) - 1)
    

def discretization_vars_consts(inst, milp, svar, qvar, hvar, dhvar, dh, number_disc, horizon, oagap):
    qvar_dict = {}
    svar_dict = {}
    for t in horizon:
        if t >= 1:
            for i in inst.arcs_discretizing.keys():
                for n in range(0, number_disc[i]):
                    svar_dict[i, f'discrete_{n}', t] = milp.addVar(vtype = GRB.BINARY) 
                
                
                for n in range(0, number_disc[i]):
                    milp.addConstr(hvar[i[0], t] - hvar[i[1], t] <= dh[i][f'discrete_{n}', t][1] * svar_dict[i, f'discrete_{n}', t] + inst.tanks_couples[i].dhmax[t]*(1-svar_dict[i, f'discrete_{n}', t]))
                    milp.addConstr(hvar[i[0], t] - hvar[i[1], t] >= dh[i][f'discrete_{n}', t][0] * svar_dict[i, f'discrete_{n}', t] + inst.tanks_couples[i].dhmin[t]*(1-svar_dict[i, f'discrete_{n}', t]))

    for t in horizon:
        if t >= 1:
            for i in inst.arcs_discretizing.keys():
                milp.addConstr(gp.quicksum(svar_dict[i, f'discrete_{n}', t] for n in range(0, number_disc[i])) == 1)
                    

                    
    #all related methods should beforehand initialized 
    for t in horizon:
        if t >= 1:
            for i, arc in inst.arcs_discretizing.items():
                for ii in arc:
                    for jj in ii:
                
                        
                        for n in range(0, number_disc[i]):
                            qvar_dict[jj, f'discrete_{n}', t] = milp.addVar(lb = inst.arcs[jj].qmin, ub = inst.arcs[jj].qmax)
                            milp.addConstr(qvar_dict[jj, f'discrete_{n}', t] >= svar_dict[i, f'discrete_{n}', t] * inst.arcs_discretization_initialization[(i, jj)].qtmin_disc[t, f'discrete_{n}'])
                            milp.addConstr(qvar_dict[jj, f'discrete_{n}', t] <= svar_dict[i, f'discrete_{n}', t] * inst.arcs_discretization_initialization[(i, jj)].qtmax_disc[t, f'discrete_{n}'])
                
                            dhvar[jj, f'discrete_{n}', t] = milp.addVar(lb=-abs(inst.arcs[jj].hlossval(inst.arcs[jj].qmin)), ub=abs(inst.arcs[jj].hlossval(inst.arcs[jj].qmax)))
                            
                        milp.addConstr(gp.quicksum(qvar_dict[jj, f'discrete_{n}', t] for n in range(0, number_disc[i])) == qvar[jj, t])
                
    
                        

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
                                milp.addConstr(dhvar[jj, f'discrete_{n_}', t]  >= c[1] * qvar_dict[jj, f'discrete_{n_}', t] + c[0] * svar_dict[i, f'discrete_{n_}', t])
                            for n, c in enumerate(cutabove):
                                milp.addConstr(dhvar[jj, f'discrete_{n_}', t]  <= c[1] * qvar_dict[jj, f'discrete_{n_}', t] + c[0] * svar_dict[i, f'discrete_{n_}', t])
                                
                        milp.addConstr(gp.quicksum(dhvar[jj, f'discrete_{n}', t] for n in range(0, number_disc[i]))   == dhvar[jj, t])


def probing_confs(inst, milp, svar, lifted_svar, qvar, hvar, dhvar, qexpr, horizon, oagap):
    
    qvar_probing = {}
    
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
                    qvar_probing[(i, j), t, conf_arc_keys[0]] = milp.addVar(lb= arc.qmin, ub= conf_arcs_values.qmax_prob[t]) 
                    
                    milp.addConstr(qvar_probing[(i, j), t, conf_arc_keys[0]] <= conf_arcs_values.qmax_prob[t] * lifted_svar[conf_arc_keys[0], t])
                    milp.addConstr(qvar_probing[(i, j), t, conf_arc_keys[0]] >= conf_arcs_values.qmin_prob[t] * lifted_svar[conf_arc_keys[0], t])
                    
                
                    
    
    for conf_arc_keys, conf_arcs_values in inst.arcs_probing.items():
        for (i, j), arc in inst.arcs.items():
            if conf_arc_keys[1] == (i, j):
                for t in horizon:


  #                  cutbelow, cutabove = oa.hlossoa(max(conf_arcs_values.qmin_prob[t], temp_min), min(conf_arcs_values.qmax_prob[t]+10e-4, temp_max+10e-4), arc.hloss, (i, j), oagap, drawgraph=False)
                    if arc.control:
                        cutbelow, cutabove = oa.hlossoa(conf_arcs_values.qmin_prob[t], conf_arcs_values.qmax_prob[t]+10e-4, arc.hloss, (i, j), oagap, drawgraph=False)
                        
                        
                        for n, c in enumerate(cutbelow):
                            milp.addConstr(dhvar[(i, j), t, conf_arc_keys[0]]   >= c[1] * qvar_probing[(i, j), t, conf_arc_keys[0]] + c[0] * lifted_svar[conf_arc_keys[0], t], name=f'hpl{n}({i},{j},{t})')
                        for n, c in enumerate(cutabove):
                            milp.addConstr(dhvar[(i, j), t, conf_arc_keys[0]]  <= c[1] * qvar_probing[(i, j), t, conf_arc_keys[0]] + c[0] * lifted_svar[conf_arc_keys[0], t], name=f'hpu{n}({i},{j},{t})')
                    
                    
                    else:
                        cutbelow, cutabove = oa.hlossoa(conf_arcs_values.qmin_prob[t], conf_arcs_values.qmax_prob[t]+10e-4, arc.hloss, (i, j), oagap, drawgraph=False)
                    
#@                    print(f'{arc}: {len(cutbelow)} cutbelow, {len(cutabove)} cutabove')
                    for n, c in enumerate(cutbelow):
                        milp.addConstr(dhvar[(i, j), t, conf_arc_keys[0]]   >= c[1] * qvar_probing[(i, j), t, conf_arc_keys[0]] + c[0] * lifted_svar[conf_arc_keys[0], t], name=f'hpl{n}({i},{j},{t})')
                    for n, c in enumerate(cutabove):
                        milp.addConstr(dhvar[(i, j), t, conf_arc_keys[0]]  <= c[1] * qvar_probing[(i, j), t, conf_arc_keys[0]] + c[0] * lifted_svar[conf_arc_keys[0], t], name=f'hpu{n}({i},{j},{t})')



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
                
                
                milp.addConstr(qvar[(i, j), t_] == gp.quicksum(qvar_probing[(i, j), t_, conf_arc_keys] for conf_arc_keys in set_confs ))
                
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



#    for j, tank in inst.tanks.items():
#        for t in horizon:
#            milp.addConstr(qexpr[j, t] <= gp.quicksum(svar[conf_tank_keys[0], t] * conf_tank_values.qexprmax_prob[t] for conf_tank_keys, conf_tank_values in inst.tanks_inflow_probing.items() if conf_tank_keys[1] == j))
#            milp.addConstr(qexpr[j, t] >= gp.quicksum(svar[conf_tank_keys[0], t] * conf_tank_values.qexprmin_prob[t] for conf_tank_keys, conf_tank_values in inst.tanks_inflow_probing.items() if conf_tank_keys[1] == j))    
            
            
    for j, tank in inst.tanks.items():
        for t in horizon:
            milp.addConstr(qexpr[j, t] <= gp.quicksum(lifted_svar[conf_tank_keys[0], t] * conf_tank_values.qexprmax_prob[t] for conf_tank_keys, conf_tank_values in inst.tanks_inflow_probing.items() if conf_tank_keys[1] == j) + tank.qtmax[t]*(1-gp.quicksum(lifted_svar[conf_tank_keys[0], t] for conf_tank_keys, conf_tank_values in inst.tanks_inflow_probing.items() if conf_tank_keys[1] == j   )))
            milp.addConstr(qexpr[j, t] >= gp.quicksum(lifted_svar[conf_tank_keys[0], t] * conf_tank_values.qexprmin_prob[t] for conf_tank_keys, conf_tank_values in inst.tanks_inflow_probing.items() if conf_tank_keys[1] == j) + tank.qtmin[t]*(1-gp.quicksum(lifted_svar[conf_tank_keys[0], t] for conf_tank_keys, conf_tank_values in inst.tanks_inflow_probing.items() if conf_tank_keys[1] == j   )))

def mir_cuts_implementation(milp, inst, D, s3_var, hvar):

    nested_dict = {}
    for key, value in inst._tank_mir_cuts.items():
    # Unpacking the original key
        tank, set_arcs, ts, tf = key
    # Creating the outer key based on ts and tf
        outer_key = (ts, tf)
    
    # Creating the middle layer if not exists
        if outer_key not in nested_dict:
            nested_dict[outer_key] = {}
    
    # Checking if set_arcs already exists in the middle dictionary
        if set_arcs not in nested_dict[outer_key]:
            nested_dict[outer_key][set_arcs] = {}
    
    # Assigning the value with tank as the most inner key
        nested_dict[outer_key][set_arcs][tank] = value





    hvar_ts = {}
    hvar_tf = {}
    
    hvar_ts_alternate = {}
    hvar_tf_alternate = {}
    for j, tank in inst.tanks.items():
        for t in inst.horizon():
            if t == 0:
                pass
            else:
                lbt = tank.head(tank.vinit) if t == 0 else D[j, t][0]
                ubt = tank.head(tank.vinit) if t == 0 else D[j, t][1]
                hvar_ts[j, t] = milp.addVar(lb=lbt-ubt, ub=ubt-lbt, name=f'hts({j},{t})')
                
       #         hvar_ts_alternate[j, t] = milp.addVar(lb=lbt-ubt, ub=0, name=f'hts({j},{t})')
                
       #         milp.addConstr(hvar_ts_alternate[j, t] == hvar[j, t] - D[j, t][1])
            
                milp.addConstr(hvar_ts[j, t] == hvar[j, t] - D[j, t][0])
            
                hvar_tf[j, t] = milp.addVar(lb=lbt-ubt, ub=ubt-lbt, name=f'htf({j},{t})')
                
       #         hvar_tf_alternate[j, t] =  milp.addVar(lb=lbt-ubt, ub=0, name=f'htf({j},{t})')
            
                milp.addConstr(hvar_tf[j, t] == hvar[j, t] - D[j, t][1])
                
       #         milp.addConstr(hvar_tf_alternate[j, t] == hvar[j, t] - D[j, t][1])
                
            


    for ii, jj in nested_dict.items():
        
        for set_arcs, rest in jj.items():
#            
           for ij, jk in rest.items():
               bn = jk.x_star
#           b = bn - sum(mir_val.mu_s*(inst.tanks[tan].h_star[ii[0]] - inst.tanks[tan].hmin[ii[0]]) for tan, mir_val in rest.items()) + sum(mir_val.mu_f*(inst.tanks[tan].h_star[ii[1]] - inst.tanks[tan].hmin[ii[1]]) for tan, mir_val in rest.items()) 
           
           b = bn + sum((mir_val.mu_s*(-inst.tanks[tan].h_star[ii[0]] + inst.tanks[tan].hmin[ii[0]])) for tan, mir_val in rest.items()) + sum((mir_val.mu_f*(-inst.tanks[tan].h_star[ii[1]] + inst.tanks[tan].hmax[ii[1]])) for tan, mir_val in rest.items()) 
           
           
#           b_alternate = bn + sum(abs(mir_val.mu_s*(inst.tanks[tan].h_star[ii[0]] - inst.tanks[tan].hmax[ii[0]])) for tan, mir_val in rest.items()) + sum(abs(mir_val.mu_f*(inst.tanks[tan].h_star[ii[1]] - inst.tanks[tan].hmax[ii[1]])) for tan, mir_val in rest.items()) 
#          


           
           b_ = math.floor(b)
           
#           b_alternate_ = math.floor(b_alternate)
           
           b_hat = b-b_
           
#           b_hat_alternate = b_alternate - b_alternate_
           
#@           if b_hat > 0.00001  :
               
#@                   milp.addConstr(s3_var[ii[0]]>= (-1/b_hat)*((gp.quicksum((mir_val.mu_s)*hvar_ts[j, ii[0]] for tan, mir_val in rest.items() if mir_val.mu_s<=0))+gp.quicksum((mir_val.mu_f)*hvar_tf[j, ii[1]] for tan, mir_val in rest.items() if mir_val.mu_f>=0))+b_+1)
                   
#@           if b_hat > 0.00001 and ii[0] !=3 :
               
#@                   milp.addConstr(s3_var[ii[0]]>= (-1/b_hat)*((gp.quicksum((mir_val.mu_s)*hvar_ts[j, ii[0]] for tan, mir_val in rest.items()))+gp.quicksum((mir_val.mu_f)*hvar_tf[j, ii[1]] for tan, mir_val in rest.items()))+b_+1)
                   
                   
 #                  milp.addConstr(s3_var[ii[0]]>= (1/b_hat_alternate)*((gp.quicksum(abs(mir_val.mu_s)*hvar_ts_alternate[j, ii[0]] for tan, mir_val in rest.items()))+gp.quicksum(abs(mir_val.mu_f)*hvar_tf_alternate[j, ii[1]] for tan, mir_val in rest.items()))+b_alternate_+1)
                   
# For items where mir_val.mu_s <= 0, sum (mir_val.mu_s * hvar_ts[j, ii[0]])
           sum_neg_mu_s = gp.quicksum(mir_val.mu_s * hvar_ts[j, ii[0]] for tan, mir_val in rest.items() if mir_val.mu_s <= 0)

# For items where mir_val.mu_f >= 0, sum (mir_val.mu_f * hvar_tf[j, ii[1]])
           sum_pos_mu_f = gp.quicksum(mir_val.mu_f * hvar_tf[j, ii[1]] for tan, mir_val in rest.items() if mir_val.mu_f >= 0)

# Full constraint
#@           constraint_expr = (-1 / b_hat) * (sum_neg_mu_s + sum_pos_mu_f) + b_ + 1
#@#@           milp.addConstr(s3_var[ii[0]] >= constraint_expr)
           if b_hat>= 0.01 and ii[0] <=23 and ii[0] >=12:
               milp.addConstr(s3_var[ii[0]] >= (-1 / b_hat) * (
    gp.quicksum((-mir_val.mu_s * hvar_ts[j, ii[0]]) if mir_val.mu_s <= 0 else 1000 for tan, mir_val in rest.items()) +
    gp.quicksum((-mir_val.mu_f * hvar_tf[j, ii[1]]) for tan, mir_val in rest.items() if mir_val.mu_f >= 0)
) + b_ + 1)

                   
                   

#    for ii, jj in nested_dict.items():
        
#        for set_arcs, rest in jj.items():
            
#           b = rest.values.x_star - sum(mir_val.mu_s*(inst.tanks[tan].h_star[ii[0]] - inst.tanks[tan].hmin[ii[0]]) for tan, mir_val in rest.items()) + sum(mir_val.mu_f*(inst.tanks[tan].h_star[ii[1]] - inst.tanks[tan].hmin[ii[1]]) for tan, mir_val in rest.items()) 
           
#           b_ = math.floor(b)
           
#           b_hat = b-b_
           
#           if b_hat > 0.0001:
               
#               milp.addConstr(s3_var[ii[0]]>= (1/b_hat)*((gp.quicksum(mir_val.mu_s*hvar_ts[j, ii[0]] for tan, mir_val in rest.items()))-gp.quicksum(mir_val.mu_f*hvar_ts[j, ii[1]] for tan, mir_val in rest.items())+b_+1))



def mir_cuts_implementation_n(milp, inst, D, s3_var, hvar):

    nested_dict = {}
    for key, value in inst._tank_mir_cuts.items():
    # Unpacking the original key
        tank, set_arcs, ts, tf = key
    # Creating the outer key based on ts and tf
        outer_key = (ts, tf)
    
    # Creating the middle layer if not exists
        if outer_key not in nested_dict:
            nested_dict[outer_key] = {}
    
    # Checking if set_arcs already exists in the middle dictionary
        if set_arcs not in nested_dict[outer_key]:
            nested_dict[outer_key][set_arcs] = {}
    
    # Assigning the value with tank as the most inner key
        nested_dict[outer_key][set_arcs][tank] = value





    hvar_ts = {}
    hvar_tf = {}
    
    hvar_ts_alternate = {}
    hvar_tf_alternate = {}
    for j, tank in inst.tanks.items():
        for t in inst.horizon():
            if t == 0:
                pass
            else:
                lbt = tank.head(tank.vinit) if t == 0 else D[j, t][0]
                ubt = tank.head(tank.vinit) if t == 0 else D[j, t][1]
                hvar_ts[j, t] = milp.addVar(lb=lbt-ubt, ub=0, name=f'hts({j},{t})')
                hvar_ts_alternate[j, t] = milp.addVar(lb=0, ub=ubt-lbt, name=f'hts({j},{t})')
                milp.addConstr(hvar_ts_alternate[j, t] == hvar[j, t] - D[j, t][0])
                milp.addConstr(hvar_ts[j, t] == hvar[j, t] - D[j, t][1])
                hvar_tf[j, t] = milp.addVar(lb=0, ub=ubt-lbt, name=f'htf({j},{t})')
                hvar_tf_alternate[j, t] =  milp.addVar(lb=lbt-ubt, ub=0, name=f'htf({j},{t})')
                milp.addConstr(hvar_tf[j, t] == hvar[j, t] - D[j, t][0])         
                milp.addConstr(hvar_tf_alternate[j, t] == hvar[j, t] - D[j, t][1])
                
            

    kkkk=0
    for ii, jj in nested_dict.items():
        
        for set_arcs, rest in jj.items():
#            
           for ij, jk in rest.items():
               bn = jk.x_star
#           b = bn - sum(mir_val.mu_s*(inst.tanks[tan].h_star[ii[0]] - inst.tanks[tan].hmin[ii[0]]) for tan, mir_val in rest.items()) + sum(mir_val.mu_f*(inst.tanks[tan].h_star[ii[1]] - inst.tanks[tan].hmin[ii[1]]) for tan, mir_val in rest.items()) 
           
           temp_h_hat = 0
           for tan, mir_val in rest.items():
               
               if mir_val.mu_s >= 0:
                   temp_h_hat += mir_val.mu_s*(-inst.tanks[tan].h_star[ii[0]] + inst.tanks[tan].hmax[ii[0]])
               elif mir_val.mu_s < 0:
                   temp_h_hat += mir_val.mu_s*(-inst.tanks[tan].h_star[ii[0]] + inst.tanks[tan].hmin[ii[0]])
    
                   
               if mir_val.mu_f >=0:
                   temp_h_hat += mir_val.mu_f*(-inst.tanks[tan].h_star[ii[1]] + inst.tanks[tan].hmax[ii[1]])
               elif mir_val.mu_f < 0:
                   temp_h_hat += mir_val.mu_f*(-inst.tanks[tan].h_star[ii[1]] + inst.tanks[tan].hmin[ii[1]])                  
               
           

 #@          b = bn + sum((mir_val.mu_s*(inst.tanks[tan].h_star[ii[0]] - inst.tanks[tan].hmax[ii[0]])) for tan, mir_val in rest.items()) + sum((mir_val.mu_f*(inst.tanks[tan].h_star[ii[1]] - inst.tanks[tan].hmin[ii[1]])) for tan, mir_val in rest.items()) 
           b = bn + temp_h_hat
           
#           b_alternate = bn + sum(abs(mir_val.mu_s*(inst.tanks[tan].h_star[ii[0]] - inst.tanks[tan].hmax[ii[0]])) for tan, mir_val in rest.items()) + sum(abs(mir_val.mu_f*(inst.tanks[tan].h_star[ii[1]] - inst.tanks[tan].hmax[ii[1]])) for tan, mir_val in rest.items()) 
#          


           
           b_ = math.floor(b)
           b_up = math.ceil(b)
           
#           b_alternate_ = math.floor(b_alternate)
           
           b_hat = b-b_
           
#           b_hat_alternate = b_alternate - b_alternate_
           
#@           if b_hat > 0.00001  :
               
#@                   milp.addConstr(s3_var[ii[0]]>= (-1/b_hat)*((gp.quicksum((mir_val.mu_s)*hvar_ts[j, ii[0]] for tan, mir_val in rest.items() if mir_val.mu_s<=0))+gp.quicksum((mir_val.mu_f)*hvar_tf[j, ii[1]] for tan, mir_val in rest.items() if mir_val.mu_f>=0))+b_+1)
                   
#@           if b_hat > 0.00001 and ii[0] !=3 :
               
#@                   milp.addConstr(s3_var[ii[0]]>= (-1/b_hat)*((gp.quicksum((mir_val.mu_s)*hvar_ts[j, ii[0]] for tan, mir_val in rest.items()))+gp.quicksum((mir_val.mu_f)*hvar_tf[j, ii[1]] for tan, mir_val in rest.items()))+b_+1)
                   
                   
 #                  milp.addConstr(s3_var[ii[0]]>= (1/b_hat_alternate)*((gp.quicksum(abs(mir_val.mu_s)*hvar_ts_alternate[j, ii[0]] for tan, mir_val in rest.items()))+gp.quicksum(abs(mir_val.mu_f)*hvar_tf_alternate[j, ii[1]] for tan, mir_val in rest.items()))+b_alternate_+1)
                   
# For items where mir_val.mu_s <= 0, sum (mir_val.mu_s * hvar_ts[j, ii[0]])
#@           sum_neg_mu_s = gp.quicksum(mir_val.mu_s * hvar_ts[j, ii[0]] for tan, mir_val in rest.items() if mir_val.mu_s <= 0)

# For items where mir_val.mu_f >= 0, sum (mir_val.mu_f * hvar_tf[j, ii[1]])
#@           sum_pos_mu_f = gp.quicksum(mir_val.mu_f * hvar_tf[j, ii[1]] for tan, mir_val in rest.items() if mir_val.mu_f >= 0)

# Full constraint
#@#@           constraint_expr = (-1 / b_hat) * (sum_neg_mu_s + sum_pos_mu_f) + b_ + 1
#@#@           milp.addConstr(s3_var[ii[0]] >= constraint_expr)
           seems_no_need = True
           
           for ta, mir_v in rest.items():
               if abs(mir_val.mu_s) <= 0.001 or abs(mir_val.mu_s) >=20:
                   seems_no_need = False
  #                 print("false now")
               elif abs(mir_val.mu_f) <= 0.001 or abs(mir_val.mu_f) >=20:
                   print("false now")
#                   seems_no_need = False
           if b_hat>= 0.001 and ii[0] and seems_no_need == True :
               print("*********some cuts*******")
               kkkk = kkkk +1
               print("######################")
               print(kkkk)
               milp.addConstr(s3_var[ii[0]] >= (1 / b_hat) * (
    gp.quicksum((mir_val.mu_s * hvar_ts_alternate[tan, ii[0]]) if mir_val.mu_s <= 0 else (mir_val.mu_s * hvar_ts[tan, ii[0]]) for tan, mir_val in rest.items()) +
    gp.quicksum((mir_val.mu_f * hvar_tf[tan, ii[1]]) if mir_val.mu_f <= 0 else (mir_val.mu_f * hvar_tf_alternate[tan, ii[1]]) for tan, mir_val in rest.items()))
 + b_up)

                   
                   