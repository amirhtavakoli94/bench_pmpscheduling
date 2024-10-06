#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 11:07:46 2020

@author: Sophie Demassey, Gratien Bonvin
"""

import csv
import math

import pandas as pd
import datetime as dt
from pathlib import Path

TARIFF_COLNAME = 'elix'
TRUNCATION = 8


def myround(val: float) -> float:
    return round(val, TRUNCATION)


def myfloat(val: str) -> float:
    return myround(float(val))


def update_min(oldlb: float, newlb: float) -> float:
    """Update lower bound only if better: returns max(oldlb, newlb)."""
    if newlb <= oldlb:
        print(f'do not update min {oldlb:.3f} to {newlb:.3f}')
        return oldlb
    return newlb


def update_max(oldub: float, newub: float) -> float:
    """Update upper bound only if better: returns min(oldub, newub)."""
    if newub >= oldub:
        print(f'do not update max {oldub:.3f} to {newub:.3f}')
        return oldub
    return newub


class _Node:
    """Generic network node. coordinates X, Y, Z (in m)."""

    def __init__(self, id_, x, y, z):
        self.id = id_
        self.coord = {'x': x, 'y': y, 'z': z}

    def altitude(self):
        return self.coord['z']


class _Tank(_Node):
    """Network node of type cylindrical water tank.

    vmin    : minimum volume value (in m3)
    vmax    : maximum volume value (in m3)
    vinit   : initial volume       (in m3)
    surface : surface              (in m2)
    """

    def __init__(self, id_, x, y, z, vmin, vmax, vinit, surface, nperiods):
        _Node.__init__(self, id_, x, y, z)
        self.vmin = vmin
        self.vmax = vmax
        self.vinit = vinit
        self.surface = surface
        self.nperiods = nperiods
        
        self.qtmin = {t: -10e+8 for t in range(0, self.nperiods)}
        self.qtmax = {t: 10e+8 for t in range(0, self.nperiods)}
        
        self.hmin = {0: self.head(self.vinit)}
        self.hmax = {0: self.head(self.vinit)}
        
        self.hmin.update({t: self.head(self.vmin) for t in range(1, self.nperiods)})
        self.hmax.update({t: self.head(self.vmax) for t in range(1, self.nperiods)})
        
        self._h_star = None #initialization h_star as None, not activated yet
        
        self._mu_star = None
        
        self._h_minbase_star = None
        

    def activate_h_star(self, value, time_step):
        """ Activate (if not already) and update the h_star dictionary for a specific time step. """
        if self._h_star is None:
            self._h_star = {}  # Initialize as an empty dictionary if not already activated
        self._h_star[time_step] = value
        
    def activate_mu_star(self, value, time_step):
        """ Activate (if not already) and update the h_star dictionary for a specific time step. """
        if self._mu_star is None:
            self._mu_star = {}  # Initialize as an empty dictionary if not already activated
        self._mu_star[time_step] = value
        
    def activate_h_minbase_star(self):
        """ Activate (if not already) and update the h_star dictionary for a specific time step. """
        if self._h_minbase_star is None:
            self._h_minbase_star = {}  # Initialize as an empty dictionary if not already activated
        self._h_minbase_star = {t: self.hmin[t] for t in range(0, self.nperiods)}

    @property
    def h_star(self):
        """ Ensure h_star is accessed only after activation. """
        if self._h_star is None:
            raise AttributeError("h_star attribute is not activated yet.")
        return self._h_star
    
    @property
    def mu_star(self):
        """ Ensure mu_star is accessed only after activation. """
        if self._mu_star is None:
            raise AttributeError("mu_star attribute is not activated yet.")
        return self._mu_star 
    
    @property
    def h_minbase_star(self):
        if self._h_minbase_star is None:
            raise AttributeError("h_minbase_star attribute is not activated yet.")
        return self._h_minbase_star
    
        
    def head(self, volume):
        return self.altitude() + volume / self.surface

    def volume(self, head):
        return (head - self.altitude()) * self.surface
    
    def update_qtbound_BT(self, timestep, qtmin, qtmax):
        self.qtmin[timestep] = qtmin
        self.qtmax[timestep] = qtmax
        
    def update_hbound_BT(self, timestep, hmin, hmax):
        self.hmin[timestep] = hmin
        self.hmax[timestep] = hmax
        
class _TankPair:
    def __init__(self, tank1: _Tank, tank2: _Tank, nperiods):
        # Ensure tank1 and tank2 are instances of Tank
        if not isinstance(tank1, _Tank) or not isinstance(tank2, _Tank):
            raise ValueError("tank1 and tank2 must be instances of Tank")
        self.tank1 = tank1
        self.tank2 = tank2
        self.nperiods = nperiods
        
        self.dhmin = {t: -1000 for t in range(1, self.nperiods)}
        self.dhmax = {t: 1000 for t in range(1, self.nperiods)}

    def update_dhbound_BT(self, timestep, dhmin, dhmax):
        self.dhmin[timestep] = dhmin
        self.dhmax[timestep] = dhmax         
        



        


class _Junction(_Node):
    """Network node of type junction.

    dmean    : mean demand (in L/s)
    dprofile : demand pattern profile id
    """

    def __init__(self, id_, x, y, z, dmean, profileid):
        _Node.__init__(self, id_, x, y, z)
        self.dmean = dmean
        self.profileid = profileid
        self.dprofile = None
        self.demands = None

    def setprofile(self, profile):
        self.dprofile = profile[self.profileid]
        self.demands = [myround(self.dmean * p) for p in self.dprofile]

    def demand(self, t):
        return self.demands[t]


class _Reservoir(_Node):
    """Network node of type infinite reservoirs (sources).

    hprofile : head profile id
    drawmax  : maximal daily withdrawal (in m3)
    drawcost : withdrawal cost  (in euro/m3)
    """

    def __init__(self, id_, x, y, z, profileid, drawmax, drawcost):
        _Node.__init__(self, id_, x, y, z)
        assert drawmax == 'NO', "max source withdrawal not completetly supported"
        self.drawmax = None if (drawmax == 'NO') else drawmax
        self.drawcost = drawcost
        if drawcost != 0:
            raise "drawcost not yet supported"
        self.profileid = profileid
        self.hprofile = None
        self.heads = None

    def setprofile(self, profile):
        self.hprofile = profile[self.profileid]
        self.heads = [myround(self.altitude() * p) for p in self.hprofile]

    def head(self, t):
        return self.heads[t]


class _Arc:
    """Generic network arc.

    id      : identifier
    nodes   : '(i,j)' with i the start node id, and j the end node id
    qmin    : minimum flow value <= q(i,j) (in L/s)
    qmax    : maximum flow value >= q(i,j) (in L/s)
    hloss   : head loss polynomial function: h(i) - h(j) = sum_n hloss[n] q(i,j)^n (L/s -> m)
    control : is the arc controllable or not ? (valved pipe or pump)
    """

    def __init__(self, id_, nodes, qmin, qmax, hloss, nperiods):
        self.id = id_
        self.nodes = nodes
        self.qmin = qmin
        self.qmax = qmax
        self.qmin_t = {t: qmin for t in range(0, nperiods)}
        self.qmax_t = {t: qmax for t in range(0, nperiods)}
        self.hloss = hloss
        self.nperiods = nperiods
        self.control = False

    def abs_qmin(self):
        return self.qmin

    def abs_qmax(self):
        return self.qmax

    def update_qbounds(self, qmin, qmax):
        self.qmin = update_min(self.qmin, qmin)
        self.qmax = update_max(self.qmax, qmax)
    
    def update_qbound_BT(self, timestep, qmin, qmax):
        self.qmin_t[timestep] = qmin
        self.qmax_t[timestep] = qmax

    def hlossval(self, q):
        """Value of the quadratic head loss function at q."""
        return self.hloss[0] + self.hloss[1] * q + self.hloss[2] * q * abs(q)

    def hlossprimitiveval(self, q):
        """Value of the primitive of the quadratic head loss function at q."""
        return self.hloss[0]*q + self.hloss[1] * q * q / 2 + self.hloss[2] * q * q * abs(q) / 3

    def hlossinverseval(self, dh):
        """Value of the inverse of the quadratic head loss function at dh."""
        sgn = -1 if self.hloss[0] > dh else 1
        return sgn * (math.sqrt(self.hloss[1]*self.hloss[1] + 4 * self.hloss[2] * (dh - self.hloss[0]))
                      - self.hloss[1]) / (2 * self.hloss[2])

    def gval(self, q, dh):
        """Value of the duality function g at (q,dh)."""
        q2 = self.hlossinverseval(dh)
        return self.hlossprimitiveval(q) - self.hlossprimitiveval(q2) + dh*q2

    def hlosstan(self, q):
        """Tangent line of the head loss function at q: f(q) + f'(q)(x-q)."""
        return [self.hloss[0] - self.hloss[2] * q * abs(q),
                self.hloss[1] + 2 * self.hloss[2] * abs(q)]

    def hlosschord(self, q1, q2):
        """Line intersecting the head loss function at q1 and q2."""
        c0 = self.hlossval(q1)
        c1 = (c0 - self.hlossval(q2)) / (q1 - q2)
        return [c0-c1*q1, c1]

    def nonnull_flow_when_on(self):
        return False

    def __str__(self):
        return f'{self.id} [{self.qmin}, {self.qmax}] {self.hloss}'


class _ControllableArc(_Arc):
    """Controllable network arc: valved pipe or pump
    dhmin   : minimum head loss value when arc is off (valve open or pump off)
    dhmax   : maximum head loss value when arc is off (valve open or pump off)
    """

    def __init__(self, id_, nodes, qmin, qmax, hloss, dhmin, dhmax, nperiods):
        super().__init__(id_, nodes, qmin, qmax, hloss, nperiods)
        #_Arc.__init__(self, id_, nodes, qmin, qmax, hloss)
        self.dhmin = dhmin
        self.dhmax = dhmax
        self.dhmin_t = {t: dhmin for t in range(0, self.nperiods)}
        self.dhmax_t = {t: dhmax for t in range(0, self.nperiods)}
        self.is_fixed_on = {t: False for t in range(0, self.nperiods)}
        self.is_fixed_off = {t: False for t in range(0, self.nperiods)}
        self.qmin_if_on = {t: qmin for t in range(0, self.nperiods)}
        self.qmax_if_on = {t: qmax for t in range(0, self.nperiods)}
        self.control = True
        

    def abs_qmin(self):
        return min(0, self.qmin)

    def abs_qmax(self):
        return max(0, self.qmax)

    def update_dhbounds(self, dhmin, dhmax):
        self.dhmin = update_min(self.dhmin, dhmin)
        self.dhmax = update_max(self.dhmax, dhmax)
    
    def update_dhbounds_BT(self, timestep, dhmin, dhmax):
        self.dhmin_t[timestep] = dhmin
        self.dhmax_t[timestep] = dhmax
        
    def update_qbounds_control_BT(self, timestep, qmin_ifon, qmax_ifon):
        self.qmin_if_on[timestep] = qmin_ifon
        self.qmax_if_on[timestep] = qmax_ifon
        
    def update_status_on_BT(self, timestep):
        self.is_fixed_on[timestep] = True
    
    def update_status_off_BT(self, timestep):
        self.is_fixed_off[timestep] = True
        

    def __str__(self):
        return f'{self.id} [{self.qmin}, {self.qmax}] {self.hloss} [{self.dhmin}, {self.dhmax}]'


class _ValvedPipe(_ControllableArc):
    """Network arc of type pipe + valve.

    valve type     : 'GV' or 'PRV' or 'CV'
    """
    def __init__(self, id_, nodes, type_, dhmin, dhmax, qmin, qmax, nperiods):
        super().__init__(id_, nodes, qmin, qmax, None, dhmin, dhmax, nperiods)
        self.type = type_
        if type_ != 'GV':
            raise NotImplementedError('pressure reducing valves are not yet supported')
        self.valve = nodes
        self.pipe = None

    def __str__(self):
        return f'V{self.id} {self.type} [{self.qmin}, {self.qmax}] {self.hloss} [{self.dhmin}, {self.dhmax}]'

    def merge_pipe(self, pipe):
        print(f'merge valve {self.nodes} + pipe {pipe.nodes}')
        self.pipe = pipe.nodes
        assert self.nodes[0] == pipe.nodes[1], f'valve {self.nodes} + pipe {pipe.nodes}'
        auxnode = self.nodes[0]
        self.nodes = (pipe.nodes[0], self.nodes[1])
        self.hloss = pipe.hloss
        print(f'valve bounds = [{self.qmin}, {self.qmax}]')
        print(f'pipe bounds = [{pipe.qmin}, {pipe.qmax}]')
        self.update_qbounds(pipe.qmin, pipe.qmax)
        return auxnode


class _Pump(_ControllableArc):
    """Network arc of type pump.

    type    : 'FSD' or 'VSD'
    power   : power polynomial function: p = sum_n power[n]q(i,j)^n (L/s -> W)
    """

    def __init__(self, id_, nodes, hloss, power, qmin, qmax, dhmin, dhmax, type_, nperiods):
        super().__init__(id_, nodes, qmin, qmax, hloss, dhmin, dhmax, nperiods)
        self.type = type_
        self.power = power
        if type_ == 'VSD':
            raise NotImplementedError('variable speed pumps are not yet supported')

    def powerval(self, q):
        assert len(self.power) == 2
        return self.power[0] + self.power[1] * q if q else 0

    def nonnull_flow_when_on(self):
        return True

    def __str__(self):
        return f'K{self.id} [{self.qmin}, {self.qmax}] {self.hloss} [{self.dhmin}, {self.dhmax}] ' \
               f'{self.power} {self.type} '





class _TankMir:
    def __init__(self, tank1: _Tank, setarcs, ts, tf):
        # Ensure tank1 and tank2 are instances of Tank
        if not isinstance(tank1, _Tank):
            raise ValueError("tank1 must be instances of tanks")
        self.tank1 = tank1
        self.setarcs = setarcs
        self.ts = ts
        self.tf = tf
        
        
        self.mu_s =  None
        self.mu_f = None
        
        self.x_star =  None
        
        
    def update_x(self, x):
        self.x_star = x
    
    def update_mu_s(self, mu_s):
        self.mu_s = mu_s
    
    def update_mu_f(self, mu_f):
        self.mu_f = mu_f
        
#        self.qmin_prob = {t: -1000 for t in range(0, self.nperiods)}
#        self.qmax_prob = {t: 1000 for t in range(0, self.nperiods)}
        
#        self.is_fixed_off = {t: False for t in range(0, self.nperiods)}

#    def update_qbound_BT(self, timestep, qmin, qmax):
#        self.qmin_prob[timestep] = qmin
#        self.qmax_prob[timestep] = qmax  
        
#    def update_status_off_BT(self, timestep):
#        self.is_fixed_off[timestep] = True

#arc2, no need 
class _ArcPair:
    def __init__(self, arc1: _Arc, conf, nperiods):
        # Ensure tank1 and tank2 are instances of Tank
        if not isinstance(arc1, _Arc):
            raise ValueError("tank1 must be instances of arcs")
        self.arc1 = arc1
        self.conf = conf
        self.nperiods = nperiods
        
        self.qmin_prob = {t: -1000 for t in range(0, self.nperiods)}
        self.qmax_prob = {t: 1000 for t in range(0, self.nperiods)}
        
        self.is_fixed_off = {t: False for t in range(0, self.nperiods)}

    def update_qbound_BT(self, timestep, qmin, qmax):
        self.qmin_prob[timestep] = qmin
        self.qmax_prob[timestep] = qmax  
        
    def update_status_off_BT(self, timestep):
        self.is_fixed_off[timestep] = True
        
        

#inflow and outflow related to activities of the related pumps based on configurations
class _TankInflowPair:
    def __init__(self, tank: _Tank, conf, nperiods):
        # Ensure tank1 and tank2 are instances of Tank
        if not isinstance(tank, _Tank):
            raise ValueError("tank1 must be instances of arcs")
        self.tank = tank
        self.conf = conf
        self.nperiods = nperiods
        
        self.qexprmin_prob = {t: -10e6 for t in range(0, self.nperiods)}
        self.qexprmax_prob = {t: 10e6 for t in range(0, self.nperiods)}
        
        self.is_fixed_off = {t: False for t in range(0, self.nperiods)}
        


    def update_qbound_BT(self, timestep, qmin, qmax):
        self.qexprmin_prob[timestep] = qmin
        self.qexprmax_prob[timestep] = qmax      
    
    def update_status_off_BT(self, timestep):
        self.is_fixed_off[timestep] = True
        
        
#I'm not sure to mention the tanks_couples or not
class _ArcDisc:
    def __init__(self, arc1:_Arc, tanks_coupling: tuple, nperiods, numb_disc = 1):
        
        self.arc1 = arc1
        self.tanks_coupling = tanks_coupling
        self.nperiods = nperiods
        
        self.qtmin_disc = {(t, f'discrete_{n}'): -1000 for n in range(0, numb_disc) for t in range(1, self.nperiods)}
        self.qtmax_disc = {(t, f'discrete_{n}'): 1000 for n in range(0, numb_disc) for t in range(1, self.nperiods)}
        
        self.is_fixed_off = {(t, f'discrete_{n}'): False for n in range(0, numb_disc) for t in range(1, self.nperiods)}
        
        
    def update_qbound_BT_disc(self, timestep, disc_number, qmin, qmax):
        self.qtmin_disc[(timestep, f'discrete_{disc_number}')] = qmin
        self.qtmax_disc[(timestep, f'discrete_{disc_number}')] = qmax
        
    def update_status_off_BT(self, timestep, disc_number):
        self.is_fixed_off[(timestep, f'discrete_{disc_number}')] = True        
        


class Instance:
    """Instance of the Pump Scheduling Problem."""

    DATADIR = Path("../data/")
    BNDSDIR = Path("../bounds/")

    def __init__(self, name, profilename, starttime, endtime, aggregatesteps):
        self.name = name

        periods, profiles = self._parse_profiles(f'{profilename}.csv', starttime, endtime, aggregatesteps)
        self.periods = periods
        self.profiles = profiles
        self.tanks = self._parse_tanks('Reservoir.csv', self._parse_initvolumes('History_V_0.csv'))
        self.junctions = self._parse_junctions('Junction.csv')
        self.reservoirs = self._parse_reservoirs('Source.csv')
        self.pumps = self._parse_pumps('Pump.csv')
        self.fpipes = self._parse_pipes('Pipe.csv')
        self._valves = self._parse_valves('Valve_Set.csv')
        self.vpipes = self._merge_pipes_and_valves()
        self.farcs = self.fpipes
        self.varcs = {**self.pumps, **self.vpipes}
        self.arcs = {**self.varcs, **self.farcs}
        self.nodes = {**self.junctions, **self.tanks, **self.reservoirs}
        self.incidence = self._getincidence()
        self.tariff = profiles[TARIFF_COLNAME]
        self.tsduration = self._get_timestepduration(periods)
        self.dependencies = self._dependencies()
        self.symmetries = self._pump_symmetric()
        self.tanks_diff = self._tank_relevance()
        
        self._number_disc = None
        
        self.configurations_probing = self._all_possible_configurations()
        
        self.tanks_couples = self._tanks_coupling_initilization()
        
        self.arcs_probing = self._arcs_coupling_initilization()
        
        self.arcs_discretizing = self._discretizing_arcs_relevance()
        
#@        self.disc_tank_couples = self._discritization_tank_pairs
        
        self._disc_cache = None
        
 #@       self.arcs_discretization_initialization = self._arcs_discretization_initilization()
        
        self.tanks_inflow_probing = self._tanks_inflow_initilization()
        
        self._tanks_lp_h = None
        
        self._tank_mir_cuts = None
        

        
        



        for r in self.reservoirs.values():
            r.setprofile(profiles)
        for j in self.junctions.values():
            j.setprofile(profiles)
            
            
    def update_tariff(self, tariff_vec):
        self.tariff = tariff_vec 
            

    def set_mir(self):
        self._tank_mir_cuts = {}  # Initialize as an empty dictionary

    def mircuts_generation(self, tank, set_arcs, start_time, final_time):
        set_arcs_key = tuple(sorted(set_arcs))
        key = (tank.id, set_arcs_key, start_time, final_time)  # Create a unique key for each instance
        self._tank_mir_cuts[key] = _TankMir(tank, set_arcs, start_time, final_time)
            
#    def mircuts_generation(self, tank, set_arcs, start_time, final_time):
#        self._tank_mir_cuts = _TankMir(tank, set_arcs, start_time, final_time)
        
            
            
    def set_number_disc(self, number_disc):
        """Method to set number_disc."""
        self._number_disc = number_disc
        self.arcs_discretization_initialization = self._arcs_discretization_initilization()


    def nperiods(self):
        return len(self.periods) - 1

    def horizon(self):
        return range(self.nperiods())

    def tsinhours(self):
        return self.tsduration.total_seconds() / 3600  # in hour

    def eleccost(self, t):
        return self.tsinhours() * self.tariff[t] / 1000  # in euro/W

    def inarcs(self, node):
        return self.incidence[node, 'in']

    def outarcs(self, node):
        return self.incidence[node, 'out']

    def inflowmin(self, node):
        return (sum(self.arcs[a].abs_qmin() for a in self.inarcs(node))
                - sum(self.arcs[a].abs_qmax() for a in self.outarcs(node)))

    def inflowmax(self, node):
        return (sum(self.arcs[a].abs_qmax() for a in self.inarcs(node))
                - sum(self.arcs[a].abs_qmin() for a in self.outarcs(node)))

    def flowtoheight(self, tank):
        return self.tsduration.total_seconds() / tank.surface / 1000  # in m / (L / s)

    def flowtovolume(self):
        return self.tsduration.total_seconds() / 1000  # in m3` / (L / s)

    #  PARSERS

    def _parsecsv(self, filename):
        csvfile = open(Path(self.DATADIR, self.name, filename))
        rows = csv.reader(csvfile, delimiter=';')
        data = [[x.strip() for x in row] for row in rows]
        return data

    def _parse_initvolumes(self, filename):
        data = self._parsecsv(filename)
        return {A[0]: myfloat(A[1]) for A in data[1:]}

    def _parse_pumps(self, filename):
        data = self._parsecsv(filename)
        return dict({(A[1], A[2]): _Pump(A[0], (A[1], A[2]),
                                         [-myfloat(A[c]) for c in [5, 4, 3]],
                                         [myfloat(A[c]) for c in [7, 6]],
                                         myfloat(A[8]), myfloat(A[9]),
                                         -myfloat(A[10]), -myfloat(A[11]),
                                         A[12], self.nperiods()) for A in data[1:]})

    def _parse_valves(self, filename):
        data = self._parsecsv(filename)
        return dict({(A[1], A[2]): _ValvedPipe(A[0], (A[1], A[2]), A[3],
                                               myfloat(A[4]), myfloat(A[5]),
                                               myfloat(A[6]), myfloat(A[7]), self.nperiods()) for A in data[1:]})

    def _parse_pipes(self, filename):
        data = self._parsecsv(filename)
        return dict({(A[1], A[2]): _Arc(A[0], (A[1], A[2]), myfloat(A[5]), myfloat(A[6]),
                                        [0, myfloat(A[4]), myfloat(A[3])], self.nperiods()) for A in data[1:]})

    def _parse_junctions(self, filename):
        data = self._parsecsv(filename)
        return dict({A[0]: _Junction(A[0], myfloat(A[1]), myfloat(A[2]),
                                     myfloat(A[3]), myfloat(A[4]), A[5]) for A in data[1:]})

    def _parse_reservoirs(self, filename):
        data = self._parsecsv(filename)
        return dict({A[0]: _Reservoir(A[0], myfloat(A[1]), myfloat(A[2]), myfloat(A[3]),
                                      A[4], A[5], myfloat(A[6])) for A in data[1:]})

    def _parse_tanks(self, filename, initvolumes):
        data = self._parsecsv(filename)
        return dict({A[0]: _Tank(A[0], myfloat(A[1]), myfloat(A[2]), myfloat(A[3]),
                                 myfloat(A[4]), myfloat(A[5]), initvolumes[A[0]],
                                 myfloat(A[6]), self.nperiods()) for A in data[1:]})

    def _parse_profiles(self, filename, starttime, endtime, aggsteps):
        data = self._parsecsv(filename)
        i = 1
        while i < len(data) and data[i][0] != starttime:
            # print(f'{data[i][0]} == {starttime}')
            i += 1
        assert i < len(data), f'starting time {starttime} not found in {filename}'

        assert data[0][1] == TARIFF_COLNAME, f'2nd column of {filename} is electricity tariff'
        profilename = data[0][1:]
        profiles = {n: [] for n in profilename}
        periods = []
        while i < len(data) and data[i][0] != endtime:
            for j, n in enumerate(profilename):
                profiles[n].append(float(data[i][j + 1]))
            periods.append(dt.datetime.strptime(data[i][0], '%d/%m/%Y %H:%M'))
            i += aggsteps
        assert i < len(data), f'end time {endtime} not found in {filename}'
        periods.append(dt.datetime.strptime(endtime, '%d/%m/%Y %H:%M'))
        return periods, profiles

    def _parse_profiles_aggregate(self, filename, starttime, endtime, aggsteps):
        data = self._parsecsv(filename)
        i = 1
        while i < len(data) and data[i][0] != starttime:
            i += 1
        assert i < len(data), f'starting time {starttime} not found in {filename}'

        assert data[0][1] == TARIFF_COLNAME, f'2nd column of {filename} is electricity tariff'
        profilename = data[0][1:]
        profiles = {n: [] for n in profilename}
        periods = []
        sumagg = [0 for _ in profilename]
        cntagg = 0
        while i < len(data) and data[i][0] != endtime:
            i += 1
            if cntagg == aggsteps:
                cntagg = 0
                for j, s in enumerate(sumagg):
                    profiles[profilename[j]].append(s / aggsteps)
                    sumagg[j] = 0
                    periods.append(dt.datetime.strptime(data[i][0], '%d/%m/%Y %H:%M'))
            cntagg += 1
            for j, s in enumerate(sumagg):
                sumagg[j] += float(data[i][j + 1])
        assert i < len(data), f'{filename}: not found end {endtime}'
        return periods, profiles

    @staticmethod
    def _get_timestepduration(periods):
        duration = periods[1] - periods[0]
        for i in range(len(periods) - 1):
            assert duration == periods[i + 1] - periods[i]
        return duration

    # !!! inverse the dh bounds for pumps in the hdf file
    # !!! do not substract the error margin  when lb = 0 (for pumps especially !)
    def parse_bounds(self, filename=None):
        """Parse bounds in the hdf file."""
        file = Path(Instance.BNDSDIR, filename if filename else self.name)
        bounds = pd.read_hdf(file.with_suffix('.hdf'), encoding='latin1').to_dict()
        margin = 1e-6
        for i, b in bounds.items():
            a = (i[0][0].replace('Tank ', 'T'), i[0][1].replace('Tank ', 'T'))
            arc = self.arcs.get(a)
            if not arc:
                (a, arc) = [(na, narc) for na, narc in self.vpipes.items() if narc.valve == a or narc.pipe == a][0]
            if i[1] == 'flow':
                arc.update_qbounds(myround(b[0] - margin), myround(b[1] + margin))
            elif i[1] == 'head':
                if a in self.pumps:
                    arc.update_dhbounds(myround(-b[0] - margin), myround(-b[1] + margin))
                else:  # if a in self.valves:
                    arc.update_dhbounds(myround(b[0] - margin), myround(b[1] + margin))

    def _merge_pipes_and_valves(self):
        vpipes = {}
        for (iv, j), valve in self._valves.items():
            inpipe = [(i, jv) for (i, jv) in self.fpipes if jv == iv]
            assert len(inpipe) == 1, f'valve {(iv,j)} is not attached to exactly one pipe: {inpipe}'
            i = inpipe[0][0]
            pipe = self.fpipes.pop((i, iv))
            auxnode = valve.merge_pipe(pipe)
            print(f'valved pipe {(i, j)} = {valve.pipe} + {valve.valve}')
            assert self.junctions[auxnode].dmean == 0
            self.junctions.pop(auxnode)
            vpipes[(i, j)] = valve
        return vpipes

    def _getincidence(self):
        incidence = {}
        for node in self.nodes:
            incidence[node, 'in'] = set()
            incidence[node, 'out'] = set()
            for arc in self.arcs:
                if arc[1] == node:
                    incidence[node, 'in'].add(arc)
                elif arc[0] == node:
                    incidence[node, 'out'].add(arc)
        return incidence

    def pumps_without_sym(self):
        """Aggregate symmetric pumps as a fictional 'sym' pump."""
        uniquepumps = set(self.pumps.keys())
        symgroup = self.symmetries
        if symgroup:
            uniquepumps -= set(symgroup)
            uniquepumps.add('sym')
        return sorted(uniquepumps, key=str)

    def _pump_symmetric(self):
        """Return a list of symmetric pumps."""
        if self.name == 'Simple_Network':
            return [('R1', 'J2'), ('R2', 'J2'), ('R3', 'J2')]
        elif self.name == 'Anytown':
            return [('R1', 'J20'), ('R2', 'J20'), ('R3', 'J20')]
        elif self.name == 'Richmond':
            return [('196', '768'), ('209', '766')]
        elif self.name == 'Vanzyl':
            return [('n10', 'n11'), ('n12', 'n13')]
        elif self.name == 'SAUR':
            return [('Arguenon_IN_1', 'Arguenon_OUT'), ('Arguenon_IN_2', 'Arguenon_OUT'),
                    ('Arguenon_IN_3', 'Arguenon_OUT'), ('Arguenon_IN_4', 'Arguenon_OUT')]
        return []

    def _dependencies(self):
        """Return 4 types of control dependencies as a dict of lists of dependent pumps/valves."""
        dep = {'p1 => p0': set(), 'p0 xor p1': set(),
               'p0 = p1 xor p2': set(), 'p1 => not p0': set()}

        if self.name == 'Richmond':
            # dep['p1 => p0'].add((('196', '768'), ('209', '766'))) already in symmetry
            dep['p1 => p0'].add((('196', '768'), ('175', '186')))
            dep['p1 => p0'].add((('312', 'TD'), ('264', '112')))
            dep['p1 => p0'].add((('264', '112'), ('312', 'TD')))

            dep['p0 xor p1'].add((('201', '770'), ('196', '768')))
            dep['p0 xor p1'].add((('321', '312'), ('264', '112')))

            dep['p0 = p1 xor p2'].add((('196', '768'), ('164', '197'), ('175', '186')))
            
        elif self.name == 'Vanzyl':
            dep['p0 = p1 xor p2'].add((('n10', 'n11'), ('n361', 'n365'), ('n362', 'n364')))
        else:
            dep = None

        return dep
    
    
    
    #recently added for arcs probing
    
    def _probing_arcs_relevance(self):
        dep = {}
        if self.name == 'Vanzyl':

            dep = {
            'conf1': set(),
            'conf2': set(),
            'conf3': set(),
            'conf4': set(),
            'conf5': set(),

        }
            

            dep['conf1'].add((('n10','n11'),('n12','n13'),('n362','n364'),('n365','t6'),('n3','t5'),('n2','n3')))
            dep['conf2'].add((('n10','n11'),('n12','n13'),('n361','n365'),('n365','t6'),('n3','t5'),('n2','n3')))
            dep['conf3'].add((('n10','n11'),('n361','n365'),('n365','t6'),('n3','t5'),('n2','n3')))
            dep['conf4'].add((('n10','n11'),('n362','n364'),('n365','t6'),('n3','t5'),('n2','n3')))
            
#@            dep['conf5'].add((('n10','n11'),('n12','n13'),('n362','n364'),('n361','n365'),('n365','t6'),('n3','t5'),('n2','n3')))
            
            

            
        if self.name == 'Richmond':
            
            dep ={
                'conf1': set(),
                'conf2': set(),
                'conf3': set(),
                'conf4': set(),
                'conf5': set(),
                'conf6': set(),
                'conf7': set(),
                'conf8': set(),
                'conf9': set(),
                'conf10': set(),
                'conf11': set(),
                'conf12': set(),
                'conf13': set(),
                'conf14': set(),
                'conf15': set(),
                'conf16': set(),
                'conf17': set(),
                
                }
            
            dep['conf1'].add((('209','766'),('766','768'),('768','770'),('633','196'),('632','209'),('Bache_O','632'),('770','771'),('771','9'),('9','42'),('42','164'),('164','175'),('164','197'),('175','186'),('186','197'),('197','284')))
            dep['conf2'].add((('209','766'),('766','768'),('768','770'),('633','196'),('632','209'),('Bache_O','632'),('770','771'),('771','9'),('9','42'),('42','164'),('164','175'),('164','197'),('175','186'),('186','197'),('197','284')))
            dep['conf3'].add((('209','766'),('766','768'),('768','770'),('633','196'),('632','209'),('Bache_O','632'),('770','771'),('771','9'),('9','42'),('42','164'),('164','175'),('164','197'),('175','186'),('186','197'),('197','284')))
            dep['conf4'].add((('209','766'),('766','768'),('768','770'),('633','196'),('632','209'),('Bache_O','632'),('770','771'),('771','9'),('9','42'),('42','164'),('164','175'),('164','197'),('175','186'),('186','197'),('197','284')))
            dep['conf5'].add((('209','766'),('766','768'),('768','770'),('633','196'),('632','209'),('Bache_O','632'),('770','771'),('771','9'),('9','42'),('42','164'),('164','175'),('164','197'),('175','186'),('186','197'),('197','284')))
            
            dep['conf6'].add((('284','125'),('353','364'),('364','TB')))
            dep['conf7'].add((('284','125'),('353','364'),('364','TB')))
            
            # here can be strangely dangerous because you define auxiliary binary variables however you might not gain as much
            
            dep['conf8'].add((('TA','4'),('4','104')))
            dep['conf9'].add((('TA','4'),('4','104')))
            dep['conf10'].add((('TA','4'),('4','104')))
            dep['conf11'].add((('TA','4'),('4','104')))
            
            # if you really wanna define the auxilairy binary vars then you need to consider the using the auxialiary vars here as well.
            
            dep['conf12'].add((('104','264'),('112','312'),('TD','320'),('320','325'),('320','321')))
            dep['conf13'].add((('104','264'),('112','312'),('TD','320'),('320','325'),('320','321')))
            
            dep['conf14'].add((('321','701'),('701','729'),('729','Reservoir_E'),('Reservoir_E','745'),('745','753'),('753','TF')))
            dep['conf15'].add((('321','701'),('701','729'),('729','Reservoir_E'),('Reservoir_E','745'),('745','753'),('753','TF')))
            
            dep['conf16'].add((('637','TC'),('636','637'),('634','635'),('249','634'),('206','249'),('104','206')))
            dep['conf17'].add((('637','TC'),('636','637'),('634','635'),('249','634'),('206','249'),('104','206')))
            
            
            
            
            
            
            
        return dep
    
    

    def _probing_tanks_inflow_relevance(self):
        dep = {}
        if self.name == 'Vanzyl':

            dep = {
            'conf1': set(),
            'conf2': set(),
            'conf3': set(),
            'conf4': set(),
            'conf5': set(),

        }
            
            dep['conf1'].add(('t5','t6'))
            dep['conf2'].add(('t5','t6'))
            dep['conf3'].add(('t5','t6'))
            dep['conf4'].add(('t5','t6'))
            dep['conf5'].add(('t5','t6'))
            
        elif self.name == 'Richmond':
            dep = {
                'conf1': set(),
                'conf2': set(),
                'conf3': set(),
                'conf4': set(),
                'conf5': set(),
                'conf6': set(),
                'conf7': set(),
                'conf8': set(),
                'conf9': set(),
                'conf10': set(),
                'conf11': set(),
                'conf12': set(),
                'conf13': set(),
                'conf14': set(),
                'conf15': set(),
                'conf16': set(),
                'conf17': set(),
                
                }
            
#            dep['conf1'] = 'TA'
            
            dep['conf1'].add(('TA',))
            dep['conf2'].add(('TA',))
            dep['conf3'].add(('TA',))
            dep['conf4'].add(('TA',))
            dep['conf5'].add(('TA',))
            
            dep['conf6'].add(('TB',))
            dep['conf7'].add(('TB',))
            
            dep['conf12'].add(('TD',))
            dep['conf13'].add(('TD',))
            
            dep['conf14'].add(('TF',))
            dep['conf15'].add(('TF',))
            
            dep['conf16'].add(('TC',))
            dep['conf17'].add(('TC',))            
            
            
                
                
                
            
            
            
        return dep
    
    
    
    
    def _tank_relevance(self):
       dep = {'tanks_diff': set()}
       if self.name == 'Vanzyl':
           dep['tanks_diff'].add(('t5', 't6'))
           
       elif self.name == 'Richmond':
           
           dep['tanks_diff'].add(('TA', 'TB'))
           dep['tanks_diff'].add(('TA', 'TC'))
           dep['tanks_diff'].add(('TA', 'TD'))
           dep['tanks_diff'].add(('TD', 'TF'))
           
       return dep


   
    
   
    #recently added for arcs discretization
    
    def _discretizing_arcs_relevance(self):
        dep = {}
        if self.name == 'Vanzyl':
            for i in self.tanks_couples.keys():

                dep[i] = set()
                    
            
                dep[i].add((('n3','t5'),('t5','n5'),('t6','n6'),('n6','n5')))
            
        return dep
    
#@    def _tank_relevance(self):
#@       dep = {'tanks_diff': set()}
#@       if self.name == 'Vanzyl':
#@           dep['tanks_diff'].add(('t5', 't6'))
           
#@       return dep
   
   
#@    def update_tank_relevance(self, related_tanks, time_step, diff_min, diff_max):
        
        
    #recently added for tanks differences
    def initialization_tanks_diff(self):
        relevant_tanks = self._tank_relevance()['tanks_diff']
        for i in relevant_tanks:
            diff_min = -1000
            diff_max = 1000
        return {(relevant_tanks, t): [diff_min, diff_max] for t in self.horizon()}




    #recently added for arcs_discritization
#@    def _arcs_discretization_initilization(self, number_disc: dict):
#@        temp_dic = self._discretizing_arcs_relevance()
#@        temp_lis = []
#@        for i, j in temp_dic.items():
#@            for jj in j:
#@                for jjj in jj: 
#@                    temp_lis.append((i, jjj))
#@                    print(f'I wanna see{(i, jj)}')
#        for i in temp_dic['conf1']:
#            temp_lis.append(i)
        
#@        print("gosh that is not what I expected of temp lis:", temp_lis)
#@        print("Available keys in self.tanks:", self.arcs.keys())
#@        print("Accessed keys:", [(A[0], A[1]) for A in temp_lis])

            
#@        return dict({A: _ArcDisc(self.arcs[A[1]], A[0], len(self.periods)-1, number_disc[A[0]]) for A in temp_lis})
    
    
    def _arcs_discretization_initilization(self, number_disc: dict = None):
        # Use the stored number_disc if not passed
        if number_disc is None:
            if self._number_disc is None:
                raise ValueError("number_disc is not set.")
        number_disc = self._number_disc

        temp_dic = self._discretizing_arcs_relevance()
        temp_lis = []
        for i, j in temp_dic.items():
            for jj in j:
                for jjj in jj: 
                    temp_lis.append((i, jjj))
###                    print(f'I wanna see{(i, jj)}')
        
###        print("gosh that is not what I expected of temp lis:", temp_lis)
###        print("Available keys in self.tanks:", self.arcs.keys())
###        print("Accessed keys:", [(A[0], A[1]) for A in temp_lis])
            
        return dict({A: _ArcDisc(self.arcs[A[1]], A[0], len(self.periods)-1, number_disc[A[0]]) for A in temp_lis})


    

    def _tanks_inflow_initilization(self):
        temp_dic = self._probing_tanks_inflow_relevance()
        temp_lis = []
        for i, j in temp_dic.items():
            for jj in j:
                for jjj in jj: 
                    temp_lis.append((i, jjj))
###                    print(f'I wanna see{(i, jj)}')
#        for i in temp_dic['conf1']:
#            temp_lis.append(i)
        
        print("gosh that is not what I expected of temp lis:", temp_lis)
        print("Available keys in self.tanks:", self.arcs.keys())
        print("Accessed keys:", [(A[0], A[1]) for A in temp_lis])

            
        return dict({A: _TankInflowPair(self.tanks[A[1]], A[0], len(self.periods)-1) for A in temp_lis})

    
    #recently added for tanks coupling
    def _tanks_coupling_initilization(self):
        temp_dic = self._tank_relevance()
        temp_lis = []
        for i in temp_dic['tanks_diff']:
            temp_lis.append(i)
        
        print(temp_lis)
        print("Available keys in self.tanks:", self.tanks.keys())
        print("Accessed keys:", [(A[0], A[1]) for A in temp_lis])

            
        return dict({A: _TankPair(self.tanks[A[0]], self.tanks[A[1]], len(self.periods)-1) for A in temp_lis})
    
    
    #check here 
    
    #recently added for probing
    def _arcs_coupling_initilization(self):
        temp_dic = self._probing_arcs_relevance()
        temp_lis = []
        for i, j in temp_dic.items():
            for jj in j:
                for jjj in jj: 
                    temp_lis.append((i, jjj))
###                    print(f'I wanna see{(i, jj)}')
#        for i in temp_dic['conf1']:
#            temp_lis.append(i)
        
        print("gosh that is not what I expected of temp lis:", temp_lis)
        print("Available keys in self.tanks:", self.arcs.keys())
        print("Accessed keys:", [(A[0], A[1]) for A in temp_lis])

            
        return dict({A: _ArcPair(self.arcs[A[1]], A[0], len(self.periods)-1) for A in temp_lis})
    

        

    
##    def update_diff_h(self, relevant_tank, time_step, min_value, max_value):
        
        
        
        
    #recently added for probing
    def _all_possible_configurations(self):
        conf = {}
        if self.name == 'Vanzyl':

            conf['conf1'] = {'ON': set({('n12','n13'), ('n362','n364')}), 'OFF': set()}
            conf['conf2'] = {'ON': set({('n12','n13'), ('n361', 'n365')}), 'OFF': set()}
            conf['conf3'] = {'ON': set({('n10','n11'),('n361','n365')}), 'OFF': set({('n12','n13')})}
            conf['conf4'] = {'ON': set({('n10','n11'),('n362','n364')}), 'OFF': set({('n12','n13')})}
            
            conf['conf5'] = {'ON': set(), 'OFF': set({('n10','n11')})}
            
        if self.name == 'Richmond':
            
            conf['conf1'] = {'ON': set({('196','768'), ('164','197')}), 'OFF': set({('209','766')})}
            conf['conf2'] = {'ON': set({('196','768'), ('164','197')}), 'OFF': set({('209','766')})}
            conf['conf3'] = {'ON': set({('209','766'),('164','197')}), 'OFF': set()}
            conf['conf4'] = {'ON': set({('209','766'),('175','186')}), 'OFF': set()}            
            conf['conf5'] = {'ON': set({('201','770')}), 'OFF': set()}
            
            conf['conf6'] = {'ON': set({('125','353')}), 'OFF': set()}
            conf['conf7'] = {'ON': set(), 'OFF': set({('125','353')})}
            
            
            conf['conf8'] = {'ON': set({('264','112'), ('635','636')}), 'OFF': set()}
            conf['conf9'] = {'ON': set({('264','112')}), 'OFF': set({('635','636')})}
            conf['conf10'] = {'ON': set({('635','636')}), 'OFF': set({('264','112')})}
            conf['conf11'] = {'ON': set(), 'OFF': set({('635','636'),('264','112')})}
            

            conf['conf12'] = {'ON': set({('264','112')}), 'OFF': set()}
            conf['conf13'] = {'ON': set(), 'OFF': set({('264','112')})}
            

            conf['conf14'] = {'ON': set({('745','753')}), 'OFF': set()}
            conf['conf15'] = {'ON': set(), 'OFF': set({('745','753')})}            
            

            conf['conf16'] = {'ON': set({('635','636')}), 'OFF': set()}
            conf['conf17'] = {'ON': set(), 'OFF': set({('635','636')})}
            
        return conf
    
    
    ###To do it requires some good architecture
#@    def _categorical_arcs(self):
#@        cat = {}
#@        if self.name == 'Richmond':
#@            
#@            dep = {
#@                'categ1': set(),
#@                'categ2': set(),
 #@               'categ3': set(),
 #@               'categ4': set(),
 #@               'categ5': set(),
 #@               'categ6': set(),
 #@               'categ7': set(),
  #@              'categ8': set(),
  #@              'categ9': set(),
  #@              'categ10': set(),
  #@              'categ11': set(),
  #@              'categ12': set(),
  #@              'categ13': set(),
  #@              'categ14': set(),
  #@              'categ15': set(),
  #@              'categ16': set(),
  #@              'categ17': set(),
                
  #@              }
            
            
  #@          dep['categ1'].add((('','')))
  #@          dep['conf2'].add(('TA'))
  #@          cat['category1'] = {}
    
    
    #@to_do
#    def _probing_tanks(self):
        
        
    #recently added for surrogate function and other stuffs
#@    def _discritization_tank_pairs(self, number: dict):
#@        """discretization of the difference of the two connected tanks
#@        taking each pair of relevant tanks, the timestep, and the number of discretizations
        
#@        returning a nested disc{pair_tanks}{[nth discretization, timestep]} 
        
#@        """
#@        tanks_pairs = self.tanks_couples
#@        dh = {}
        
#@        for j, tanks in tanks_pairs.items():
#@                dh[j] = {}
#@                for timestep in range(1, len(self.periods)-1):
#@                    for i in range(0, number[j]):
#@                        step = (tanks.dhmax[timestep] - tanks.dhmin[timestep])/number[j]
#@                        count = i
#@                        dh[j][f'discrete_{count}', timestep] = [tanks.dhmin[timestep] + count*step, tanks.dhmin[timestep] + (count+1)*step] 
                
#@        return dh

    def _discritization_tank_pairs(self):
        """discretization of the difference of the two connected tanks
        taking each pair of relevant tanks, the timestep, and the number of discretizations
        
        returning a nested disc{pair_tanks}{[nth discretization, timestep]} 
        
        """
        tanks_pairs = self.tanks_couples
        dh = {}
        
        for j, tanks in tanks_pairs.items():
                dh[j] = {}
                for timestep in range(1, len(self.periods)-1):
                    for i in range(0, self._number_disc[j]):
                        step = (tanks.dhmax[timestep] - tanks.dhmin[timestep])/self._number_disc[j]
                        count = i
                        dh[j][f'discrete_{count}', timestep] = [tanks.dhmin[timestep] + count*step, tanks.dhmin[timestep] + (count+1)*step] 
                
        return dh
    
    @property
    def disc_tank_couples(self):
        """Lazy loading the discretization values only when needed."""
        if self._disc_cache is None:
            self._disc_cache = self._discritization_tank_pairs()
        return self._disc_cache

    def force_recompute(self):
        """Force a re-computation of the discretization values."""
        self._disc_cache = self._discritization_tank_pairs()
    
            

    def tostr_basic(self):
        return f'{self.name} {self.periods[0]} {self.horizon()}'

    def tostr_network(self):
        return f'{len(self.pumps)} pumps, {len(self.vpipes)} valved pipes, {len(self.fpipes)} fixed pipes, ' \
               f'{len(self.tanks)} tanks'

    def print_all(self):
        print(f'{self.tostr_basic()} {self.tostr_network()}')
        print(f'{len(self.arcs)} arcs:')
        for i, a in self.arcs.items():
            print(i)
            print(str(a))

    # def transcript_bounds(self, csvfile):
    #    """Parse bounds in the hdf file."""
    #     file = Path(Instance.BNDSDIR, self.name)
    #     pd.read_hdf(file.with_suffix('.hdf')).to_csv(csvfile)

    def parsesolution(self, filename):
        csvfile = open(filename)
        rows = csv.reader(csvfile, delimiter=',')
        data = [[x.strip() for x in row] for row in rows]
        csvfile.close()
        assert float(data[0][1]) == self.nperiods(), f"different horizons in {data[0]} and {self.tostr_basic()}"
        inactive = {t: set((A[0], A[1]) for A in data[1:] if A[t + 2] == '0') for t in self.horizon()}
        return inactive
    

            
        

