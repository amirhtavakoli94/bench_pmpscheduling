# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 16:33:34 2022

@author: amirhossein.tavakoli
"""

from instance import Instance
from datetime import datetime
import convexrelaxation as rel
import lpnlpbb as bb
import csv
import graphic
from hydraulics import HydraulicNetwork
from pathlib import Path
from stats import Stat
import os
import new_cvx_last_cleaning
import bound_tightening_C1_cop
import BT_Cleaning_C2
import time
import new_cvx_last_cleaning_grad_cuts
import pickle
import datetime as dt



OA_GAP = 1e-2
MIP_GAP = 1e-6

BENCH = {
    'FSD': {'ntk': 'Simple_Network', 'D0': 1, 'H0': '/01/2013 00:00'},
    'RIC': {'ntk': 'Richmond', 'D0': 21, 'H0': '/05/2013 07:00'},
    'ANY': {'ntk': 'Anytown', 'D0': 1, 'H0': '/01/2013 00:00'},
    'VAN': {'ntk': 'Vanzyl', 'D0': 1, 'H0': '/01/2011 08:00'},
#    'VAN': {'ntk': 'Vanzyl', 'D0': 21, 'H0': '/05/2013 07:00'}
}
PROFILE = {'s': 'Profile_5d_30m_smooth', 'n': 'demand_tank0_2011_lazy_bum_van_gaus'}
STEPLENGTH = {'12': 4, '24': 2, '48': 1}


# ex of instance id: "FSD s 24 3"
def makeinstance(instid: str):
    a = instid.split()
    assert len(a) == 4, f"wrong instance key {instid}"
    d = BENCH[a[0]]
    dbeg = f"{(d['D0'] + int(a[3]) - 1):02d}" + d['H0']
    dend = f"{(d['D0'] + int(a[3])):02d}" + d['H0']
    return Instance(d['ntk'], PROFILE[a[1]], dbeg, dend, STEPLENGTH[a[2]])

def makeinstance_mein_loop(instid: str, ite):
    a = instid.split()
    assert len(a) == 4, f"wrong instance key {instid}"
    d = BENCH[a[0]]
    dbeg= dt.datetime.strftime(dt.datetime.strptime('01/01/2011 08:00','%d/%m/%Y %H:%M')+(ite)*dt.timedelta(hours=24),'%d/%m/%Y %H:%M')
    dend= dt.datetime.strftime(dt.datetime.strptime('01/01/2011 08:00','%d/%m/%Y %H:%M')+(ite+1)*dt.timedelta(hours=24),'%d/%m/%Y %H:%M')
#    dbeg = f"{(d['D0'] + int(a[3]) - 1):02d}" + d['H0']
#    dend = f"{(d['D0'] + int(a[3])):02d}" + d['H0']
    return Instance(d['ntk'], PROFILE[a[1]], dbeg, dend, STEPLENGTH[a[2]])

FASTBENCH = [
    'VAN n 24 2',
##    'VAN s 12 2',
##    'VAN s 12 3',
##    'VAN s 12 4',
#    'RIC s 12 2',
#    'VAN s 24 1',
#    'VAN s 24 3',
#    'VAN s 24 4',
#    'VAN s 48 1',
#    'VAN s 48 2',
#    'VAN s 48 3',
#    'VAN s 48 4'
]

OUTDIR = Path("../output/")
defaultfilename = Path(OUTDIR, f'resall.csv')
SOLFILE = Path(OUTDIR, f'solutions.csv')


# RECORD (default: gurobi manages incumbent), FATHOM (cut feas int nodes) or CVX (MIP relaxation only)
# NOADJUST (default: no adjustment heuristic), ADJUST (run heur) or ADJUSTNCUT (cut with heur solutions)
MODES = {"solve": ['RECORD', 'FATHOM', 'CVX'],
         "adjust": ['NOADJUST', 'ADJUST', 'ADJUSTNCUT']}


def parsemode(modes):
    pm = {k: mk[0] for k, mk in MODES.items()}
    if modes is None:
        return pm
    elif type(modes) is str:
        modes = [modes]
    for k, mk in MODES.items():
        for mode in mk:
            if mode in modes:
                pm[k] = mode
                break
    return pm


def solve(instance, oagap, mipgap, drawsolution, stat, arcvals=None):
    print('***********************************************')
    print(instance.tostr_basic())
    print(instance.tostr_network())

    print("obbt: parse bounds")
    try:
        instance.parse_bounds()
    except UnicodeDecodeError as err:
        print(f'obbt bounds not read: {err}')
    # instance.print_all()

    print("create model")
    cvxmodel = rel.build_model(instance, oagap, arcvals=arcvals)
    # cvxmodel.write('convrel.lp')
    cvxmodel.params.MIPGap = mipgap
    cvxmodel.params.timeLimit = 5
    # cvxmodel.params.OutputFlag = 0
    cvxmodel.params.Threads = 1
    # cvxmodel.params.FeasibilityTol = 1e-5

    print("solve model")
    costreal, plan = bb.solveconvex(cvxmodel, instance, drawsolution=drawsolution) if stat.solveconvex() \
        else bb.lpnlpbb(cvxmodel, instance, stat.modes, drawsolution=drawsolution)

    stat.fill(cvxmodel, costreal)
    print('***********************************************')
    print(f"solution for {instance.tostr_basic()}")
    print(stat.tostr_basic())

    cvxmodel.terminate()
    return costreal, plan


def solveinstance(instid, oagap=OA_GAP, mipgap=MIP_GAP, modes=None, drawsolution=True, stat=None, file=defaultfilename):
    instance = makeinstance(instid)
    stat = Stat(parsemode(modes)) if stat is None else stat
    now = datetime.now().strftime("%y%m%d-%H%M")
    print(now)
    cost, plan = solve(instance, oagap, mipgap, drawsolution, stat)
    if cost:
        writeplan(instance, plan, f"{now}, {instid}, {cost},")
    fileexists = os.path.exists(file)
    f = open(file, 'a')
    if not fileexists:
        f.write(f"date, oagap, mipgap, mode, ntk T day, {stat.tocsv_title()}\n")
    f.write(f"{now}, {oagap}, {mipgap}, {stat.getsolvemode()}, {instid}, {stat.tocsv_basic()}\n")
    f.close()


def writeplan(instance, activity, preamb, solfile=SOLFILE):
    assert len(activity) == instance.nperiods() and len(activity[0]) == len(instance.arcs)
    plan = {a: [activity[t][a] for t in instance.horizon()] for a in instance.varcs}
    f = open(solfile, 'a')
    f.write(f"{preamb} {plan}\n")
    f.close()
    
def writevolume(instance, activity, preamb, solfile=SOLFILE):
#    assert len(activity) == instance.nperiods() and len(activity[0]) == len(instance.arcs)
    vol = {a: [activity[t][a] for t in range(0,len(instance.horizon())+1)] for a in instance.tanks}
    f = open(solfile, 'a')
    f.write(f"{preamb} {vol}\n")
    f.close()



def testsolution(instid, solfilename, oagap=OA_GAP, mipgap=MIP_GAP, modes='CVX', drawsolution=True):
    instance = makeinstance(instid)
    inactive = instance.parsesolution(solfilename)
    network = HydraulicNetwork(instance, feastol=mipgap)
    flow, hreal, volume, nbviolations = network.extended_period_analysis(inactive, stopatviolation=False)
    cost = sum(instance.eleccost(t) * sum(pump.power[0] + pump.power[1] * flow[t][a]
                                          for a, pump in instance.pumps.items() if a not in inactive[t])
               for t in instance.horizon())
    print(f'real plan cost (without draw cost) = {cost} with {nbviolations} violations')
    graphic.pumps(instance, flow)
    graphic.tanks(instance, flow, volume)

    stat = Stat(parsemode(modes))
    arcvals = {(a, t): 0 if a in inactive[t] else 1 for a in instance.varcs for t in instance.horizon()}
    solve(instance, oagap, mipgap, drawsolution, stat, arcvals=arcvals)


def testfullsolutions(instid, solfilename, oagap=OA_GAP, mipgap=MIP_GAP, modes='CVX', drawsolution=True):
    csvfile = open(solfilename)
    rows = csv.reader(csvfile, delimiter=',')
    data = [[float(x.strip()) for x in row] for row in rows]
    csvfile.close()

    print('************ TEST SOLUTIONS ***********************************')
    instance = makeinstance(instid)
    print(instance.tostr_basic())
    print(instance.tostr_network())

    print("obbt: parse bounds")
    try:
        instance.parse_bounds()
    except UnicodeDecodeError as err:
        print(f'obbt bounds not read: {err}')

    stat = Stat(parsemode(modes))
    print("create model")
    for i, d in enumerate(data):
        print(f"create model {i}")
        cvxmodel = rel.build_model(instance, oagap)
        rel.postsolution(cvxmodel, d)
        cvxmodel.params.MIPGap = mipgap
        cvxmodel.params.timeLimit = 1200
        print("solve model")
        costreal, plan = bb.solveconvex(cvxmodel, instance, drawsolution=drawsolution) if stat.solveconvex() \
            else bb.lpnlpbb(cvxmodel, instance, stat.modes, drawsolution=drawsolution)

        stat.fill(cvxmodel, costreal)
        print('***********************************************')
        print(f"solution for {instance.tostr_basic()}")
        print(stat.tostr_basic())

        cvxmodel.terminate()



def bound_t_check(instance, Z, C, D, P0, P1, oagap, first_time, arcvals= None):
    
    initial_time = time.time()
    tau= 5
    Z,C,D,P0,P1 = bound_tightening_C1_cop.variables_bounds(instance, Z, C, D, P0, P1, tau, first_time, arcvals= None)
    return Z, C, D, P0, P1, instance
        


def solve_BT(instance, Z, C, D, P0, P1, S_n_all, number_disc, oagap, mipgap, drawsolution, stat, arcvals=None):
    print('***********************************************')
    print(instance.tostr_basic())
    print(instance.tostr_network())

    print("obbt: parse bounds")

    print("create model")
    cvxmodel = new_cvx_last_cleaning_grad_cuts.build_model(instance, Z, C, D, P0, P1, S_n_all, oagap, arcvals=arcvals)
    cvxmodel.params.MIPGap = mipgap
    cvxmodel.params.timeLimit = 3600
    cvxmodel.params.Threads = 10
    print("solve model")
    br_ch_ini= time.time()
    costreal, plan, volu, lower_boun, solu = bb.solveconvex(cvxmodel, instance, drawsolution=drawsolution) if stat.solveconvex() \
        else bb.lpnlpbb(cvxmodel, instance, stat.modes, drawsolution=drawsolution)
        
        
        
####    def build_model(inst: Instance, Z, C, D, P0, P1, sn_b_all, number_disc, oagap: float, arcvals=None):
    
    br_ch_fin= time.time()
    time_br_ch= br_ch_fin-br_ch_ini
    return costreal, plan, volu, lower_boun, solu, time_br_ch, cvxmodel.NodeCount

def solveinstance_BT(instid, oagap=OA_GAP, mipgap=MIP_GAP, modes=None, drawsolution=True, stat=None, file=defaultfilename):
    instance = makeinstance(instid)
    stat = Stat(parsemode(modes)) if stat is None else stat
    now = datetime.now().strftime("%y%m%d-%H%M")
    print(now)
    Z, C, D, P0, P1, N_Cardi, dic_ca, s1p, piq, number_disc, instance = solveinstance_just_bound_check(instance, oagap)        
    cost, plan, volu, lower_bound, time_br_check = solve_BT(instance, Z, C, D, P0, P1, N_Cardi, number_disc, oagap, mipgap, drawsolution, stat)

    if cost:
        writeplan(instance, plan, f"{now}, {instid}, {cost},")
        writevolume(instance, volu, f"{now}, {instid}, {cost},")

    return cost, plan, volu, lower_bound, time_br_check, Z, C, D, P0, P1, N_Cardi, instance



def solvebench_BT(bench, oagap=OA_GAP, mipgap=MIP_GAP, modes=None, drawsolution=False):
    stat = Stat(parsemode(modes))
    now = datetime.now().strftime("%y%m%d-%H%M")
    resfilename = Path(OUTDIR, f'res{now}.csv')
    for i in bench:
        aa,bb,cc, low_b, t_br_check, z_tau, c_tau, d_tau, time_req = solveinstance_BT(i, oagap=oagap, mipgap=mipgap, drawsolution=drawsolution, stat=stat, file=resfilename)



def solveinstance_just_bound_check(instid, ii, oagap=OA_GAP, mipgap=MIP_GAP, modes=None, drawsolution=True, stat=None, file=defaultfilename):
#    instance = makeinstance(instid)
    instance = makeinstance_mein_loop(instid_n, ii)
    stat = Stat(parsemode(modes)) if stat is None else stat
    now = datetime.now().strftime("%y%m%d-%H%M")
    print(now)
    Z, C, D, P0, P1 = {}, {}, {}, {}, {}
#    first_time = True
#    Z, C, D, P0, P1, instance = bound_t_check(instance, Z, C, D, P0, P1, oagap, first_time, arcvals=None)    
#    Z = BT_Cleaning_C2.probing_arcs(instance, Z, C, D, P0, P1, 1, oagap = OA_GAP)
#    Z = BT_Cleaning_C2.probing_tanks_inflow(instance, Z, C, D, P0, P1, 1, oagap)
    # for vanzyl: coupled tank is {('t5', 't6'): number of discretization}
    number_disc = {('t5', 't6'): 5}
    instance.set_number_disc(number_disc)
#    Z = BT_Cleaning_C2.probing_discretization(instance, Z, C, D, P0, P1, 1, number_disc, oagap)
#@    N_Cardi = BT_Cleaning_C2.cardinality_vi(instance, Z, C, D, P0, P1, number_disc, oagap)
    
#1    Z, C, D, P0, P1, instance = bound_t_check(instance, Z, C, D, P0, P1,oagap, False, arcvals=None)  
    
    ###N_Cardi = BT_Cleaning_C2.cardinality_vi(instance, Z, C, D, P0, P1, number_disc, oagap)
    
#    with open('..//bounds//VAN s 12.pkl', 'rb') as file:
#        loaded_data = pickle.load(file)

    with open(f'..//bounds//50_instances//my_objects_van25_1_5discrete_50_instances_{ii-1}.pkl', 'rb') as file:
        loaded_data = pickle.load(file)
    
    

    Z = loaded_data['Z']
    C = loaded_data['C']
    D = loaded_data['D']
    P0 = loaded_data['P0']
    P1 = loaded_data['P1']
    N_Cardi = loaded_data['nnm']
    
    instance = loaded_data['instance']
    
    instance.set_mir()
    
#    Z, C, D, P0, P1, instance = bound_t_check(instance, Z, C, D, P0, P1, oagap, False, arcvals=None) 

#    Z = BT_Cleaning_C2.probing_arcs(instance, Z, C, D, P0, P1, 1, oagap = OA_GAP)
#    Z = BT_Cleaning_C2.probing_tanks_inflow(instance, Z, C, D, P0, P1, 1, oagap)
#    Z = BT_Cleaning_C2.probing_discretization(instance, Z, C, D, P0, P1, 1, number_disc, oagap)
    
###    N_Cardi = BT_Cleaning_C2.cardinality_vi(instance, Z, C, D, P0, P1, number_disc, oagap)
#@    N_Cardi = {}
    dic_ca, s1p, piq = BT_Cleaning_C2.mir_cuts(instance, Z, C, D, P0, P1, N_Cardi, oagap)

#@    dic_ca, s1p, piq = {}, {}, {}
    
 #@   s_t = new_cvx_last_cleaning_grad_cuts.build_model(instance, Z, C, D, P0, P1, N_Cardi, oagap)
#@    s_t = s_t.relax()
#@    s_t.optimize()
    return Z, C, D, P0, P1, N_Cardi, dic_ca, s1p, piq, number_disc, instance

            
instid='VAN n 24 2'
instance=makeinstance(instid)   

 
            
#@@@Z, C, D, P0, P1, Z_tau, C_tau, D_tau, dh_star, S_n_all_, tim_req= solveinstance_just_bound(instid, oagap=OA_GAP, mipgap=MIP_GAP, modes=None, drawsolution=False) 

###@@@@Z, C, D, P0, P1, nnm, bt_h, sip, piq, number_disc, instance = solveinstance_just_bound_check(instid, oagap=OA_GAP, mipgap=MIP_GAP, modes=None, drawsolution=False) 





def solveinstance_BT_mein(instid, oagap=OA_GAP, mipgap=MIP_GAP, modes=None, drawsolution=True, stat=None, file=defaultfilename):
    instance = makeinstance(instid)
    stat = Stat(parsemode(modes)) if stat is None else stat
    now = datetime.now().strftime("%y%m%d-%H%M")
    print(now)
    Z, C, D, P0, P1, N_Cardi, dic_ca, s1p, piq, number_disc, instance = solveinstance_just_bound_check(instance, oagap)        
    cost, plan, volu, lower_bound, solu, time_br_check = solve_BT(instance, Z, C, D, P0, P1, N_Cardi, number_disc, oagap, mipgap, drawsolution, stat)

    if cost:
        writeplan(instance, plan, f"{now}, {instid}, {cost},")
        writevolume(instance, volu, f"{now}, {instid}, {cost},")

    return cost, plan, volu, lower_bound, solu, time_br_check, Z, C, D, P0, P1, N_Cardi, instance


def solve_milp(instance, Z, C, D, P0, P1, S_n_all, number_disc, oagap, mipgap, drawsolution, stat, arcvals=None):
    print('***********************************************')
    print(instance.tostr_basic())
    print(instance.tostr_network())

    print("obbt: parse bounds")

    print("create model")
    cvxmodel = new_cvx_last_cleaning.build_model(instance, Z, C, D, P0, P1, S_n_all, number_disc, oagap, arcvals=arcvals)
    cvxmodel.params.MIPGap = mipgap
    cvxmodel.params.timeLimit = 3000
    cvxmodel.params.Threads = 10
    cv1=cvxmodel.relax()
    cv1.optimize()
#    cvxmodel.optimize()





def solve_BT_MIR(instance, Z, C, D, P0, P1, S_n_all, number_disc, oagap, mipgap, drawsolution, stat, arcvals=None):
    print('***********************************************')
    print(instance.tostr_basic())
    print(instance.tostr_network())

    print("obbt: parse bounds")

    print("create model")
    cvxmodel = new_cvx_last_cleaning_grad_cuts.build_model(instance, Z, C, D, P0, P1, nnm, 0.01)
    cvxmodel.params.MIPGap = mipgap
    cvxmodel.params.timeLimit = 7200
    cvxmodel.params.Threads = 10
    print("solve model")
    br_ch_ini= time.time()
    costreal, plan, volu, lower_boun, solu = bb.solveconvex(cvxmodel, instance, drawsolution=drawsolution) if stat.solveconvex() \
        else bb.lpnlpbb(cvxmodel, instance, stat.modes, drawsolution=drawsolution)
        
        
####    def build_model(inst: Instance, Z, C, D, P0, P1, sn_b_all, number_disc, oagap: float, arcvals=None):
    
    br_ch_fin= time.time()
    time_br_ch= br_ch_fin-br_ch_ini
    return costreal, plan, volu, lower_boun, solu, time_br_ch





FASTBENCH = [
    'VAN n 24 2',
##    'VAN s 12 2',
##    'VAN s 12 3',
##    'VAN s 12 4',
#    'RIC s 12 2',
#    'VAN s 24 1',
#    'VAN s 24 3',
#    'VAN s 24 4',
#    'VAN s 48 1',
#    'VAN s 48 2',
#    'VAN s 48 3',
#    'VAN s 48 4'
]

stat = Stat(parsemode(None))
####solve_milp(instance, Z, C, D, P0, P1, nnm, number_disc, 0.01, 10e-6, False, stat)
dic_macht={}
dict_fe={}
number_disc = 5
for ii in range(49, 50):

    instid_n = f'VAN n 24 {ii}'
#@    instance_n = makeinstance(instid_n)
#@    instance.update_tariff(instance_n.tariff)
    Z, C, D, P0, P1, nnm, bt_h, sip, piq, number_disc, instance = solveinstance_just_bound_check(instid_n,ii, oagap=OA_GAP, mipgap=MIP_GAP, modes=None, drawsolution=False) 
    cost, plan, volu, lower_bound, solu, time_br_check, numb_nodes = solve_BT(instance, Z, C, D, P0, P1, nnm, number_disc, 0.01, 10e-6, False, stat)
#    for itera in range(317, 364):

    now = datetime.now().strftime("%y%m%d-%H%M")
    resfilename = Path(OUTDIR, f'res{now}.csv')
        

        

        
    dic_macht[ii]= { 'height': volu, 'plan': plan, 'all solu': solu, 'lower_bound': lower_bound, 'time_br':time_br_check, 'number_node': numb_nodes, 'cost':cost }
    dict_fe[ii]= { 'height': volu, 'plan': plan, 'instance':instance, 'all solu': solu, 'lower_bound': lower_bound, 'time_br':time_br_check, 'number_node': numb_nodes, 'cost':cost}


    

#    with open(f'check_s_correc2011_24_van_dep_new_prob_noCard_MIR_disc_agg{ii}_seemsfine.pickle', 'wb') as f:
                    # write JSON data to file
#                    pickle.dump(dict_fe[ii], f)
                    
                    
            
                
#    with open(f'check_s_correc_macht2011_24_van_dep_new_prob_noCard_MIR_disc_agg{ii}_seemsfine.pickle', 'wb') as f:
                    # write JSON data to file
#                    pickle.dump(dic_macht[ii], f)

    


#@solve_BT_MIR(instance, Z, C, D, P0, P1, nnm, number_disc, 0.01, 10e-6, False, stat)
#@instid2 = 'VAN s 48 2'   
#@instance2 = makeinstance(instid2)  
#@instance.update_tariff(instance2.tariff)   
##cost, plan, volu, lower_bound, solu, time_br_check = solve_BT(instance, Z, C, D, P0, P1, nnm, number_disc, 0.01, 10e-6, False, stat)
#@solve_milp(instance, Z, C, D, P0, P1, nnm, number_disc, 0.01, 10e-6, False, stat)