#!/usr/bin/env python3
'''
#todo enforcing the dimension of the system model matrices (A,B,C,D) to 2 before entering the MPC, there are 2 solutions:
    - self.B = np.array(..., ndmin=2)
    - self.B = np.atleast_2d(self.B)
#todo using of np.block to assemble matrices
#todo using of np.tipe for repetition of matrices
#todo optimization for speed:
    - np.array and np.matrix have copy=False option
'''

from mpc_dataclass import MpcParamData, SignalsData
from lti_mpc import LinearTimeIvariantMPC
from dummy_model import SisoLtiModel
import utils

import time
import math
import numpy as np
import scipy
# import matplotlib.pyplot as plt
import logging



def main_loop():
    print('\nStart MPC loop')
    t_sim = 10
    tic = time.perf_counter()
    for i in np.arange(start=0.0, stop=t_sim, step=Ts, dtype=np.single):
        t = time.monotonic()
        # signals.ref = math.sin(t)
        signals.ref = -5

        # compute optimal input
        u = lti_mpc.on_step(signals=signals)

        # update system state
        out = system.on_step(x=signals.x, u=u)

        print(f'u = {u}')
        print(f'x_dot = {out.x_dot}')
        print(f'y = {out.y}')

        # x(k+1) = x(k) + v(k) dt
        signals.x = signals.x + out.x_dot*Ts
        # signals.x = signals.x - out.x_dot*Ts
        # signals.x = out.x_dot*Ts - signals.x



    toc = time.perf_counter()
    print(f'  Execution time: {toc-tic:0.4}s')



if __name__ == '__main__':
    print('Start MPC control')

    # system
    T_sys = 0.1
    system = SisoLtiModel(Ts=T_sys)

    # system model
    sys_model = system.get_model

    # signals
    x0 = np.ndarray(shape=(1,1), dtype=float, buffer=np.array([0], copy=True))  # vo
    ref = np.array([0]) # A*sin(t)
    u0 = np.array([0])  # x
    signals = SignalsData(x=x0, ref=ref, u=u0)

    # mpc parameters
    nx = 1
    ny = 1
    nu = 1
    Ts = 0.1
    Np = 2
    Nc = 2
    Q = np.diag([1])
    R = np.diag([1])
    x_c_min = np.array([[-10]], ndmin=2)
    x_c_max = np.array([[10]], ndmin=2)
    u_c_min = np.array([[-2]], ndmin=2)
    u_c_max = np.array([[2]], ndmin=2)
    du_c_min = np.array([[-2]], ndmin=2)
    du_c_max = np.array([[2]], ndmin=2)
    mpc_params = MpcParamData(
        nx=nx, ny=ny, nu=nu, Ts=Ts, Np=Np, Nc=Nc, Q=Q, R=R,
        x_c_min=x_c_min, x_c_max=x_c_max,
        u_c_min=u_c_min, u_c_max=u_c_max,
        du_c_min=du_c_min, du_c_max=du_c_max)

    # mpc
    lti_mpc = LinearTimeIvariantMPC(system_model=sys_model, mpc_parameters=mpc_params)


    # data validation
    utils.check_system_dimensions(system_model=sys_model, mpc_params=mpc_params, signals=signals)


    # start mpc main loop
    main_loop()


    print('\nEnd MPC control')

