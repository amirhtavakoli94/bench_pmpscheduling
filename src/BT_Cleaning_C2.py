# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 13:59:15 2024

@author: amirhossein.tavakoli
"""

import numpy as np
import cond_C2_Cleaned
import new_probing_cleanup
import new_cvx_last_cleaning_grad_cuts
import reduced_version_new_cvx_last_cleaning_grad_cuts


OA_GAP=0.01

def probing_arcs(instance, Z, C, D, P0, P1, Tau, oagap = OA_GAP):
    ### for various configurations, it finds conditional bounds ####
    """ for various configurations, it finds vonditional bounds
        Inputs: instance
        Z: a dict for variables flow bounds
        C: static flow bounds
        D: head and difference of coupled bounds
        P0: probed to zero
        P1: probed to one
        Tau: the number of iteration
        oagap: the gap between the curve and cutting planes
        
        Output:
            conditional bounds related to possible configurations (activity of the pumps and valves)
 
    """
    ###here is obviously wrong!!!!
    for com in instance.arcs_probing.items():
        conf = com[0][0]
        a = com[0][1]
        for t in instance.horizon():
            
            if (t, conf) in P0:
                for com in instance.arcs_probing.items():
                    if conf == com[0][0]:
#                    conf = com[0][0]
                        a = com[0][1]
                        com[1].update_status_off_BT(t)
                        bt_arc_dic= dict([((a, t, conf), [0 , 10e-4]) ])
                        com[1].update_qbound_BT(t, 0, 10e-4)
                    
                    pass
            else:
   #@         arc_min_cond = cond_C2_Cleaned.build_common_model(instance, Z, C, D, P0, P1, a, t, 'MIN', conf, 'full_SOS', accuracy = 0.01, envelop = 0.05, oagap = OA_GAP)
                arc_min_cond = cond_C2_Cleaned.build_common_model(instance, Z, C, D, P0, P1, a, t, 'MIN', 'full_SOS', 0.01, 0.05, OA_GAP, conf)
                arc_min_cond.optimize()
                arc_min = arc_min_cond.getObjective()
                if arc_min_cond.status == 2:
                    arc_min = arc_min.getValue()
#@                arc_max_cond = cond_C2_Cleaned.build_common_model(instance, Z, C, D, P0, P1, a, t, 'MAX', conf, 'full_SOS', accuracy = 0.01, envelop = 0.05, oagap = OA_GAP)
                    arc_max_cond = cond_C2_Cleaned.build_common_model(instance, Z, C, D, P0, P1, a, t, 'MAX', 'full_SOS', 0.01, 0.05, OA_GAP, conf)
                    arc_max_cond.optimize()
                    arc_max = arc_max_cond.getObjective()
                    arc_max = arc_max.getValue()
                
                    if abs(arc_min) <= 10e-4:
                        arc_min = 0
                
                    bt_arc_dic= dict([((a, t, conf), [arc_min , arc_max]) ])
                    com[1].update_qbound_BT(t, arc_min, arc_max)
                else:
                    P0[(t, conf)] = 10e-4
                    
                    P0[(a, t, conf)] = 10e-4
                
                    com[1].update_status_off_BT(t)
                
                    bt_arc_dic= dict([((a, t, conf), [0 , 10e-4]) ])
                    com[1].update_qbound_BT(t, 0, 10e-4)
            
                
            Z= {**Z, **bt_arc_dic}
        
    return Z

def probing_tanks_inflow(instance, Z, C, D, P0, P1, Tau, oagap):
    
    """ for various configurations, it finds vonditional bounds for tanks inflow variables
        Inputs: instance
        Z: a dict for variables flow bounds
        C: static flow bounds
        D: head and difference of coupled bounds
        P0: probed to zero
        P1: probed to one
        Tau: the number of iteration
        oagap: the gap between the curve and cutting planes
        
        Output:
            conditional bounds for tanks inflow qexpr related to possible configurations (activity of the pumps and valves)
 
    """
    for conf_tank in instance.tanks_inflow_probing.items():
        conf = conf_tank[0][0]
        tank = conf_tank[0][1]
        for t in instance.horizon():
#@            tank_min_cond = cond_C2_Cleaned.build_common_model(instance, Z, C, D, P0, P1, tank, t, 'MIN', conf, 'full_SOS', accuracy = 0.01, envelop = 0.05, oagap = OA_GAP)
            tank_min_cond = cond_C2_Cleaned.build_common_model(instance, Z, C, D, P0, P1, tank, t, 'MIN', 'full_SOS', 0.01, 0.05, OA_GAP, conf)
            tank_min_cond.optimize()
            tankinflow_min_probing = tank_min_cond.getObjective()
            if tank_min_cond.status == 2:
                tank_inf_min = tankinflow_min_probing.getValue()
#@                tank_max_cond = cond_C2_Cleaned.build_common_model(instance, Z, C, D, P0, P1, tank, t, 'MAX', conf, 'full_SOS', accuracy = 0.01, envelop = 0.05, oagap = OA_GAP)
                tank_max_cond = cond_C2_Cleaned.build_common_model(instance, Z, C, D, P0, P1, tank, t, 'MAX', 'full_SOS', 0.01, 0.05, OA_GAP, conf)
                tank_max_cond.optimize()
                tankinflow_max_probing = tank_max_cond.getObjective()
                tank_inf_max = tankinflow_max_probing.getValue()
                
                if abs(tank_inf_min) <= 10e-4:
                    tank_inf_min = 0
                
                bt_arc_dic = dict([((tank, conf, t), [tank_inf_min, tank_inf_max]) ])
                
                conf_tank[1].update_qbound_BT(t, tank_inf_min, tank_inf_max)
                
            else:
                
                
                P0[(tank, conf, t)] = 10e-4
                
                conf_tank[1].update_status_off_BT(t)
                
                bt_arc_dic = dict([((tank, conf, t), [0, 10e-4] )])
                
                conf_tank[1].update_qbound_BT(t, 0, 10e-4)                
                
                
            Z = {**Z, **bt_arc_dic}
        
    return Z


def probing_discretization(instance, Z, C, D, P0, P1, Tau, number_disc: dict, oagap):
    
    """ for various discretization of the tanks heads' differences, it finds vonditional bounds for arcs variables
        Inputs: instance
        Z: a dict for variables flow bounds
        C: static flow bounds
        D: head and difference of coupled bounds
        P0: probed to zero
        P1: probed to one
        Tau: the number of iterations
        number_disc = {tanks_couples: number of discretization}
        oagap: the gap between the curve and cutting planes
        
        Output:
            conditional bounds for related arcs flows related to possible discretization of the tanks heads
 
    """
    dh = instance.disc_tank_couples
    
    instance.set_number_disc(number_disc)
    for coupled_tank_arc, relevant_arc in instance.arcs_discretization_initialization.items():
        tanks_ckey = coupled_tank_arc[0]
        arcs_ckey = coupled_tank_arc[1]
        
        temp_dict = {}
        for t in instance.horizon():
            if t>=1:
                for n in range(0, number_disc[tanks_ckey]):                
                    temp_dict[tanks_ckey] = [dh[tanks_ckey][f'discrete_{n}', t][0], dh[tanks_ckey][f'discrete_{n}', t][1]]
                    arc_min_disc = cond_C2_Cleaned.build_common_model(instance, Z, C, D, P0, P1, arcs_ckey, t, 'MIN', 'full_SOS', 0.01, 0.05, OA_GAP, [], temp_dict, Discrete = True)
                    arc_min_disc.optimize()
                    if arc_min_disc.status == 2:
                        arc_min_dis = arc_min_disc.getObjective()
                        arc_min = arc_min_dis.getValue()
                        arc_max_disc = cond_C2_Cleaned.build_common_model(instance, Z, C, D, P0, P1, arcs_ckey, t, 'MAX', 'full_SOS', 0.01, 0.05, OA_GAP, [], temp_dict, Discrete = True)
                        arc_max_disc.optimize()
                        arc_max_dis = arc_max_disc.getObjective()
                        arc_max = arc_max_dis.getValue()
                    
                        if abs(arc_min) <= 10e-4:
                            arc_min = 0
                    
                        bt_arc_dic = dict([((arcs_ckey, f'discrete_{n}', t), [arc_min, arc_max]) ])
                        
                        relevant_arc.update_qbound_BT_disc(t, n, arc_min, arc_max)
                    
                
                    else:
                        
                        P0[(arcs_ckey, f'discrete_{n}', t)] = 10e-4

                        bt_arc_dic = dict([((arcs_ckey, f'discrete_{n}', t), [0, 10e-4] )])
                    
                        relevant_arc.update_qbound_BT_disc(t, n, 0, 10e-4)        
                        
                        relevant_arc.update_status_off_BT(t, n)
                
                
                    Z = {**Z, **bt_arc_dic}
        
    return Z                


def cardinality_vi(instance, Z, C, D, P0, P1, number_disc, oagap):
    
    """ for various discretization of the tanks heads' differences, it finds the min number of required number of ON pumps until and from each times tep
        Inputs: instance
        Z: a dict for variables flow bounds
        C: static flow bounds
        D: head and difference of coupled bounds
        P0: probed to zero
        P1: probed to one
        Tau: the number of iterations
        number_disc = {tanks_couples: number of discretization}
        oagap: the gap between the curve and cutting planes
        
        Output:
            it finds the min number of required number of ON pumps until and from each times while the relaxation can be different from the one deployed in branch and check
 
    """
    nm = {}
    
    for t in instance.horizon():

#@    for t in range(0, 1):
        min_numb = new_probing_cleanup.build_model(instance, Z, C, D, P0, P1, t, number_disc, False, True, oagap)
        min_numb.optimize()
        n1 = min_numb.getObjective()
        nn1 = n1.getValue()
        min_numb_reverse = new_probing_cleanup.build_model(instance, Z, C, D, P0, P1, t, number_disc, True, True, oagap)
        min_numb_reverse.optimize()
        n2 = min_numb_reverse.getObjective()
        nn2 = n2.getValue()
        nm[t] = [nn1, nn2]
    
    return nm

def mir_cuts(instance, Z, C, D, P0, P1, sn_b_all, oagap):
    
    k_h_=[]                           
    cv_p1 = new_cvx_last_cleaning_grad_cuts.build_model(instance, Z, C, D, P0, P1, sn_b_all, oagap)
    cv_p1=cv_p1.relax()
    cv_p1.optimize()
    for i in range(0,len(cv_p1.getVars())):
#%#for i in range(0,0):
            k_h_.append([cv_p1.VarName[i], cv_p1.X[i]])

            k_h_arr=np.array(k_h_)
            keydicts= k_h_arr[:, 0]
            bt_h_p= dict([
                (key, [float(k_h_arr[i][1])]) for key, i in zip(keydicts, range(len(k_h_arr)))])
            
    for t in instance.horizon():
        for j, tank in instance.tanks.items():
            tank.activate_h_star(bt_h_p[f'ht({j},{t})'][0], t)
            
    
    
            
            
    
                            
    bt_pih_prob= {}
    pi_n_ts_tf_Q_1= {}
    s_n_ts_tf_lp_Q_1= {}
                            
    s_n_ts_tf_lp_Q_1_sep= {} 
    
    
    bt_pih_prob= {}
    pi_n_ts_tf_Q_1= {}
    s_n_ts_tf_lp_Q_1= {}
                            
    s_n_ts_tf_lp_Q_1_sep= {}
    
    aa = instance.arcs
    pumps_considered = set()
    for (i, j), arc in instance.arcs.items():
            if arc.control:
                if aa[(i, j)].type == 'FSD':
                    pumps_considered.add((i, j))
    pumps_considered = sorted(pumps_considered)
                    
    for K, k in instance.tanks.items():
        dur=1
        for ts in range(1, len(list(instance.horizon()))-dur):
            
####    for ts in range(1, 1):
###    for ts in range(1, 1):
            tf= ts+dur
            
            instance.mircuts_generation(k, pumps_considered, ts, tf)
#            hs= bt_h_p[f'ht({K},{ts})'][0]
#            hf= bt_h_p[f'ht({K},{tf})'][0]
##            for t_ in range(16, 24):
            s_t = reduced_version_new_cvx_last_cleaning_grad_cuts.build_model(instance, Z, C, D, P0, P1, bt_h_p, ts, tf, oagap, arcvals=None)
            s_t = s_t.relax()
            s_t.optimize()
            sts_tf_=s_t.getObjective()
            sts_tf= sts_tf_.getValue()
            h_srzt=[]
        
        
            k_h_prob_1=[]
            for i in range(0,len(s_t.getVars())):
                k_h_prob_1.append([s_t.VarName[i], s_t.RC[i]])

                k_h_prob_1_arr=np.array(k_h_prob_1)
                keydicts= k_h_prob_1_arr[:, 0]
                bt_op_prob_1= dict([
                (key, [float(k_h_prob_1_arr[i][1])]) for key, i in zip(keydicts, range(len(k_h_prob_1_arr)))])
####            for j, tank in instance.tanks.items():
                
            bt_pi_h_prob_1= dict([( (K, tuple(pumps_considered), ts, tf), [ bt_op_prob_1[f'ht({K},{ts})'][0], bt_op_prob_1[f'ht({K},{tf})'][0]  ]) ])
            instance._tank_mir_cuts[(K, tuple(pumps_considered), ts, tf)].update_x(sts_tf)
            
            instance._tank_mir_cuts[(K, tuple(pumps_considered), ts, tf)].update_x(sts_tf)
            instance._tank_mir_cuts[(K, tuple(pumps_considered), ts, tf)].update_mu_s(bt_op_prob_1[f'ht({K},{ts})'][0])
            instance._tank_mir_cuts[(K, tuple(pumps_considered), ts, tf)].update_mu_f(bt_op_prob_1[f'ht({K},{tf})'][0])
            
            
            bt_pih_prob= {**bt_pih_prob,**bt_pi_h_prob_1}
            pi_n_ts_tf_Q_1={**pi_n_ts_tf_Q_1, **bt_pi_h_prob_1}
        

        
        

        
#                h_min= h_min_milp_.getValue()
#                h_max= h_max_milp_.getValue()
        
            bt_s_lp_Q_1= dict([( (K, tuple(pumps_considered), ts, tf), [sts_tf]) ])
        
            s_n_ts_tf_lp_Q_1={**s_n_ts_tf_lp_Q_1, **bt_s_lp_Q_1}
###        H_n_n={**H_n_n, **bt_h_sztT}
            pi_n_ts_tf_Q_1={**pi_n_ts_tf_Q_1, **bt_pi_h_prob_1}
###        pi_n1= {**pi_n1, **bt_pi_h1}

            SQ_n1= s_n_ts_tf_lp_Q_1
            piQ_1= pi_n_ts_tf_Q_1
    
    
#    build_model(inst: Instance, Z, C, D, P0, P1, sn_b_all, bt_h_p, ts, tf, oagap: float, arcvals=None)
    
    return bt_h_p, SQ_n1, piQ_1
    
        
                


                
    

