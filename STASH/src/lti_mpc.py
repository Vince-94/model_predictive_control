#!/usr/bin/env python3
'''
Linear Model Predictive Control
------------------------
1. System prediction
2. Optimization problem: Quadratic Program
'''
from mpc_dataclass import SystemModelData, MpcParamData, SignalsData
from system_discretization import SystemDiscretization
from system_prediction import SystemPrediction
from qp_optimization import QuadraticProgramOptimization

import numpy as np

import logging



class LinearTimeIvariantMPC():
    def __init__(self, system_model: SystemModelData, mpc_parameters: MpcParamData):
        self.nu = mpc_parameters.nu

        self.pred = SystemPrediction(sys_model=system_model, mpc_params=mpc_parameters)
        self.pred.prediction_invariant_matrices()
        self.pred_mtx = self.pred.get_prediction_matrices

        self.qp_opt = QuadraticProgramOptimization(
            nx=mpc_parameters.nx, nu=mpc_parameters.nu, Np=mpc_parameters.Np, Nc=mpc_parameters.Nc, prediction_matrices=self.pred_mtx)


    def on_step(self, signals: SignalsData) -> np.array:
        self.pred.prediction_varying_matrices(signals=signals)
        self.pred_mtx = self.pred.get_prediction_matrices
        self.U_star = self.qp_opt.qp_solve(x0=signals.x, prediction_matrices=self.pred_mtx)

        self.u_star = self.U_star[0:self.nu]

        return self.u_star


if __name__ == '__main__':
    print('Linear Model Predictive Control')
