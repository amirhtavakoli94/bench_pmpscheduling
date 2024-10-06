# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 13:59:15 2024

@author: amirhossein.tavakoli
"""

# Removed unnecessary imports
import BT_one_st_with_diff
import bt_h_
import numpy as np


OA_GAP=0.01

def initial_bounds(instance):
    
    
    if instance.name == 'Simple_Network':
        z_flow= np.load('..//data//Simple_Network//Bound0_flow_arcs_fsd.npy',allow_pickle=True)
        zz = z_flow.tolist()
        c_head = np.load('..//data//Simple_Network//Bound0_head_arcs_fsd.npy',allow_pickle=True)
        cc = c_head.tolist()
    
    
    if instance.name == 'Richmond':
        z_flow = np.load('..//data//Richmond//Bound0_flow_arcs_ric.npy',allow_pickle=True)
        zz = z_flow.tolist()
        c_head = np.load('..//data//Richmond//Bound0_head_arcs_ric.npy',allow_pickle=True)
        cc = c_head.tolist()
        
        
    if instance.name == 'Vanzyl':
        z_flow = np.load('..//data//Vanzyl//Bound0_flow_arcs_van.npy',allow_pickle=True)
        zz = z_flow.tolist()
        c_head = np.load('..//data//Vanzyl//Bound0_head_arcs_van.npy',allow_pickle=True)
        cc = c_head.tolist()
    
    

    bt_tank_infl = []
    for K, k in instance.tanks.items():
        for t in range(0, len(list(instance.horizon()))):
            bt_tank_infl_min = float(sum(instance.arcs[a].qmin for a in instance.inarcs(K)) - sum(instance.arcs[a].qmax for a in instance.outarcs(K)))
            bt_tank_infl_max = float(sum(instance.arcs[a].qmax for a in instance.inarcs(K)) - sum(instance.arcs[a].qmin for a in instance.outarcs(K)))
            bt_tank_infl.append([K, t, bt_tank_infl_min, bt_tank_infl_max])
            bt_tank_infl_arr =  np.array(bt_tank_infl)
            
            k.update_qtbound_BT(t, bt_tank_infl_min, bt_tank_infl_max)
            
            
        
    a2 = np.array(list(map(int, bt_tank_infl_arr[:, 1])))
    x2 = list(zip(bt_tank_infl_arr[:, 0], a2[:]))
    bt_tank_infl_dic = dict([
        (key, [float(bt_tank_infl_arr[i][2]), float(bt_tank_infl_arr[i][3])]) for key, i in
        zip(x2, range(len(bt_tank_infl_arr)))])
    
    BT = {}    
    for K, k in instance.arcs.items():
        for t in range(0, len(list(instance.horizon()))):
            bt_arcs= dict([((K, t), [zz[K][0], zz[K][1]])])
            
            if k.control:
                k.update_qbounds_control_BT(t, zz[K][0], zz[K][1])
            else:
                k.update_qbound_BT(t, zz[K][0], zz[K][1])
            BT = {**BT, **bt_arcs}
    BT={**BT, **bt_tank_infl_dic}
        
        
    pump_off_dic = {}       
    for K, k in instance.arcs.items():
        if k.control:
            for t in range(0, len(list(instance.horizon()))):
                bt_arcs = dict([((K, t), [cc[K][0], cc[K][1]])])
                k.update_dhbounds_BT(t, cc[K][0], cc[K][1])
                pump_off_dic = {**pump_off_dic, **bt_arcs}


    tank_lev=[]
    for K, k in instance.tanks.items():
        for t in range(0, len(list(instance.horizon()))):
            if t==0:
                tank_lev_min = k.head(k.vinit)
                tank_lev_max = k.head(k.vinit)
            else:
                tank_lev_min = k.head(k.vmin)
                tank_lev_max = k.head(k.vmax)
                tank_lev.append([K, t, tank_lev_min, tank_lev_max])
                tank_lev_arr = np.array(tank_lev)
    
    a3 = np.array(list(map(int, tank_lev_arr[:, 1])))
    x3 = list(zip(tank_lev_arr[:, 0], a3[:]))
    tank_lev_dic= dict([
        (key, [float(tank_lev_arr[i][2]), float(tank_lev_arr[i][3])]) for key, i in
        zip(x3, range(len(tank_lev_arr)))])
    
    
    tank_lev_diff={}
    for t in range(0, len(list(instance.horizon()))):
        tank_lev_diff['T56', t] = [-250, 250]
        




    P0 = {}
    P1 = {}
    Z = BT
    C = pump_off_dic
    D = tank_lev_dic
    
    D = {**D, **tank_lev_diff}
    
    return Z, C, D, P0, P1


def build_and_optimize_model_tanks(instance, Z, C, D, P0, P1, K, t, model_type, accuracy, envelop, oagap, arcvals):
    model = BT_one_st_with_diff.build_common_model(instance, Z, C, D, P0, P1, K, K, t, model_type, 'partial_SOS', accuracy=accuracy, envelop=envelop, oagap=oagap)
    model.optimize()
    objective = model.getObjective()
    return objective.getValue()


def optimize_arc(instance, arc, time, optimization_type, parameters):
    """Optimizes arc for a given time and type (min or max) with specified parameters."""
    optimizer = BT_one_st_with_diff.build_common_model(instance, *parameters, arc[0], arc[1], time, optimization_type, 'full_SOS', accuracy=0.02, envelop=0.04, oagap=OA_GAP)
    optimizer.optimize()
    if optimizer.status == 2:
        objective_value = optimizer.getObjective().getValue()
        return objective_value
    else:
        return None

def update_arcs_flow_bounds(k, arc, time, min_val, max_val, control_flag, q_bounds):
    """Updates bounds, probabilities, and control flags based on optimization results."""
    if min_val is not None:
        if abs(min_val) <= 10e-4:
            min_val = 0
#@        updated_bounds = (max(min_val, q_bounds[0]), min(max_val, q_bounds[1]))
        
        updated_bounds = max(min_val, q_bounds[0]), min(max_val, q_bounds[1])
    else:
        min_val = max_val = q_bounds[0]  # Assumes min_val == max_val in this context
#@        updated_bounds = (q_bounds[0], q_bounds[0] + 0.0001)
        updated_bounds = q_bounds[0], q_bounds[0] + 10e-4

    # Control logic based on the flag
    if control_flag:
        k.update_qbounds_control_BT(time, *updated_bounds)
    else:
        k.update_qbound_BT(time, *updated_bounds)

    return {(arc, time): updated_bounds}





def variables_bounds(instance, Z, C, D, P0, P1, Tau, first_time = False, arcvals=None):
    OA_GAP = 0.01
    
    if first_time == True:
        Z, C, D, P0, P1 = initial_bounds(instance)

#    Tau = 0

    for tau in range(0, Tau):
        Z, P0 = update_arc_bounds(instance, Z, C, D, P0, P1)
        Z = update_tank_bounds(instance, Z, C, D, P0, P1)

        C, P1 = update_pump_bounds(instance, Z, C, D, P0, P1)
        D = update_hbounds(instance, Z, C, D, P0, P1)
        
        D = update_difference_tanks(instance, Z, C, D, P0, P1)
#@        D = update_tanks_diff(instance, Z, C, D, P0, P1)
    
    return Z, C, D, P0, P1

def update_arc_bounds(instance, Z, C, D, P0, P1):
    for (i, j), k in instance.arcs.items():
        arc = (i, j)
        for time in range(0, len(list(instance.horizon()))):
            if (arc, time) in P0:
                pass
            else:
                min_val = optimize_arc(instance, arc, time, 'ARC_MIN', [Z, C, D, P0, P1])
                if min_val is not None:
                    max_val = optimize_arc(instance, arc, time, 'ARC_MAX', [Z, C, D, P0, P1])
                    bounds_update = update_arcs_flow_bounds(k, arc, time, min_val, max_val, k.control, [k.qmin, k.qmax])
                else:

#@                    P0.update({(arc, time): 0.01})
                    k.update_status_off_BT(time)
                    P0[(arc, time)] = 10e-4
                    bounds_update = update_arcs_flow_bounds(k, arc, time, min_val, max_val, k.control, [k.qmin, k.qmax])
                    # Handle failed optimization or other logic as needed
                Z.update(bounds_update)
    return Z, P0

def update_tank_bounds(instance, Z, C, D, P0, P1):
    for K, k in instance.tanks.items():
        for t in range(0, len(list(instance.horizon()))):
            t_min = build_and_optimize_model_tanks(instance, Z, C, D, P0, P1, K, t, 'TANK_MIN', 0.02, 0.05, OA_GAP, None)
            t_max = build_and_optimize_model_tanks(instance, Z, C, D, P0, P1, K, t, 'TANK_MAX', 0.02, 0.05, OA_GAP, None)

            k.update_qtbound_BT(t, t_min, t_max)

            bt_t_infl = dict([((K, t), [t_min-10e-4, t_max+10e-4]) ])
            Z = {**Z, **bt_t_infl}
    return Z

def update_pump_bounds(instance, Z, C, D, P0, P1):
    for (mm, nn), k in instance.arcs.items():
        for t in range(0, len(list(instance.horizon()))):
#@        for t in range(0, 0):
            
            if k.control:
                if ((mm, nn), t) in P1:
                    pass
                else:
                    dh_off_mi = BT_one_st_with_diff.build_common_model(instance, Z, C, D, P0, P1, mm, nn, t, 'PUMP_OFF_MIN', 'full_SOS', accuracy=0.02, envelop=0.05, oagap=OA_GAP)
                    dh_off_ma = BT_one_st_with_diff.build_common_model(instance, Z, C, D, P0, P1, mm, nn, t, 'PUMP_OFF_MAX', 'full_SOS', accuracy=0.02, envelop=0.05, oagap=OA_GAP)

                    dh_off_mi.optimize()
                    dh_off_ma.optimize()

                    dh_p_mi_off = dh_off_mi.getObjective()
                    dh_p_ma_off = dh_off_ma.getObjective()

                    if dh_off_mi.status == 2:
                        d0_min_ = dh_p_mi_off.getValue()
                        d0_max_ = dh_p_ma_off.getValue()

                        k.update_dhbounds_BT(t, d0_min_, d0_max_)

                        bt_P_off = dict([(((mm, nn), t), [max(d0_min_-10e-4, k.dhmin), min(d0_max_+10e-4, k.dhmax)]) ])

                    else:
                        bt_P_off = dict([(((mm, nn), t), [k.dhmin, k.dhmax]) ])
                        P1[(mm, nn), t] = 1
                        k.update_status_on_BT(t)
                    C={**C, **bt_P_off}
    return C, P1

def update_hbounds(instance, Z, C, D, P0, P1):
    for K, k in instance.tanks.items():
        for t_ in range(1, len(list(instance.horizon()))):
            h_min_milp = bt_h_.build_model_BT_h(instance, Z, C, D, P0, P1, K, t_, 'MILP', 'oa_cuts', accuracy=0.02, envelop=0.5, Minim=True, two_h=False, oagap=OA_GAP)
            h_max_milp = bt_h_.build_model_BT_h(instance, Z, C, D, P0, P1, K, t_, 'MILP', 'oa_cuts', accuracy=0.02, envelop=0.5, Minim=False, two_h=False, oagap=OA_GAP)

            h_min_milp.optimize()
            h_max_milp.optimize()

            h_min = h_min_milp.ObjBound
            h_max = h_max_milp.ObjBound

            bt_h_milp = dict([((K, t_), [max(h_min, D[K, t_][0]), min(h_max, D[K, t_][1])]) ])

            k.update_hbound_BT(t_, h_min, h_max)
            D = {**D, **bt_h_milp}
    return D


def update_difference_tanks(instance, Z, C, D, P0, P1):
    for k, tanks in instance.tanks_couples.items():
        for t_ in range(1, len(list(instance.horizon()))):
            
            h_min_milp = bt_h_.build_model_BT_h(instance, Z, C, D, P0, P1, k, t_, 'MILP', 'oa_cuts', accuracy=0.02, envelop=0.5, Minim=True, two_h=False, oagap=OA_GAP)
            h_max_milp = bt_h_.build_model_BT_h(instance, Z, C, D, P0, P1, k, t_, 'MILP', 'oa_cuts', accuracy=0.02, envelop=0.5, Minim=False, two_h=False, oagap=OA_GAP)

            h_min_milp.optimize()
            h_max_milp.optimize()

            h_min = h_min_milp.ObjBound
            h_max = h_max_milp.ObjBound

            bt_h_milp = dict([((k, t_), [h_min, h_max]) ])
            
            tanks.update_dhbound_BT(t_, h_min, h_max)
            
            D = {**D, **bt_h_milp}
    return D
            
def update_tanks_diff(instance, Z, C, D, P0, P1):
    for K in instance.tanks_diff['tanks_diff']:
        for t_ in range(1, len(list(instance.horizon()))):
            
            h_min_milp = bt_h_.build_model_BT_h(instance, Z, C, D, P0, P1, K, t_, 'MILP', 'oa_cuts', accuracy=0.02, envelop=0.5, Minim=True, two_h=False, oagap=OA_GAP)
            h_max_milp = bt_h_.build_model_BT_h(instance, Z, C, D, P0, P1, K, t_, 'MILP', 'oa_cuts', accuracy=0.02, envelop=0.5, Minim=False, two_h=False, oagap=OA_GAP)

            h_min_milp.optimize()
            h_max_milp.optimize()

            h_min = h_min_milp.ObjBound
            h_max = h_max_milp.ObjBound

            bt_h_milp = dict([((K, t_), [h_min, h_max]) ])
            D = {**D, **bt_h_milp}
    return D
