#!/usr/bin/env python3
'''
System Prediction
-----------------
- state/output:
    A_bar = nx Np,nx | B_bar = nx Np,nu Nc
    C_bar = ny Np,nx | D_bar = ny Np,nu Nc
- contraints:
    X_c_min = nx Np,1 | X_c_min = nx Np,1
    U_c_min = nu Nc,1 | U_c_min = nu Nc,1
    dU_c_min = nu Nc,1 | dU_c_min = nu Nc,1
- weithg matrices
    Q_bar = ny Np,ny Np | R = nu Np,nu Np
- signals:
    REF = ny Np,1
    U0 = nu Nc,1
'''
from mpc_dataclass import SystemModelData, PredictionModelData, MpcParamData, SignalsData
from system_discretization import SystemDiscretization

import numpy as np
import logging


class SystemPrediction():
    def __init__(self, sys_model: SystemModelData, mpc_params: MpcParamData):

        # init mpc_params
        self.nx, self.ny, self.nu  = mpc_params.nx, mpc_params.ny, mpc_params.nu
        self.Np, self.Nc = mpc_params.Np, mpc_params.Nc
        self.Ts = mpc_params.Ts
        self.Q, self.R = mpc_params.Q, mpc_params.R
        self.x_c_min, self.x_c_max = mpc_params.x_c_min, mpc_params.x_c_max
        self.u_c_min, self.u_c_max = mpc_params.u_c_min, mpc_params.u_c_max
        self.du_c_min, self.du_c_max = mpc_params.du_c_min, mpc_params.du_c_max

        # system discretization
        self.sys_disc = SystemDiscretization(cont_sys_model=sys_model, Ts=self.Ts).euler_discretization()

        # init prediction matrices
        self.A_bar = np.zeros((self.nx*self.Np, self.nx))
        self.B_bar = np.zeros((self.nx*self.Np, self.nu*self.Nc))
        self.C_bar = np.zeros((self.ny*self.Np, self.nx))
        self.D_bar = np.zeros((self.ny*self.Np, self.nu*self.Nc))

        self.REF = np.zeros((self.ny*self.Np, 1))
        self.U0 = np.zeros((self.nu*self.Nc, 1))

        self.Q_bar = np.zeros((self.ny*self.Np, self.ny*self.Np))
        self.R_bar = np.zeros((self.nu*self.Nc, self.nu*self.Nc))

        self.X_c_min = np.zeros((self.nx*self.Np, 1))
        self.X_c_max = np.zeros((self.nx*self.Np, 1))
        self.U_c_min = np.zeros((self.nu*self.Nc, 1))
        self.U_c_max = np.zeros((self.nu*self.Nc, 1))
        self.dU_c_min = np.zeros((self.nu*self.Nc, 1))
        self.dU_c_max = np.zeros((self.nu*self.Nc, 1))

    def prediction_invariant_matrices(self) -> PredictionModelData:
        '''
        Compute time-invariant prediction matrices
        '''
        for i in range(0, self.Np):
            self.A_bar[i*self.nx:(i+1)*self.nx, 0:self.nx] = np.linalg.matrix_power(self.sys_disc.A, i)
            self.C_bar[i*self.ny:(i+1)*self.ny, 0:self.nx] = np.matmul(self.sys_disc.C, np.linalg.matrix_power(self.sys_disc.A, i))
            self.X_c_min[i*self.nx:(i+1)*self.nx, :] = self.x_c_min
            self.X_c_max[i*self.nx:(i+1)*self.nx, :] = self.x_c_max

            if i < self.Nc:
                self.U_c_min[i*self.nu:(i+1)*self.nu, :] = self.u_c_min
                self.U_c_max[i*self.nu:(i+1)*self.nu, :] = self.u_c_max
                self.dU_c_min[i*self.nu:(i+1)*self.nu, :] = self.du_c_min
                self.dU_c_max[i*self.nu:(i+1)*self.nu, :] = self.du_c_max

            for j in range(0, self.Np):

                if (j <= i and j < self.Nc):

                    self.A_bar_pow = np.linalg.matrix_power(self.sys_disc.A, i-j)
                    self.B_bar[i*self.nx:(i+1)*self.nx, j*self.nu:(j+1)*self.nu] = np.matmul(self.A_bar_pow, self.sys_disc.B)
                    self.D_bar[i*self.ny:(i+1)*self.ny, j*self.nu:(j+1)*self.nu] = np.matmul(np.matmul(self.sys_disc.C, self.A_bar_pow), self.sys_disc.B)

                    if i == j:
                        self.Q_bar[i*self.ny:(i+1)*self.ny, i*self.ny:(i+1)*self.ny] = self.Q
                        self.R_bar[i*self.nu:(i+1)*self.nu, i*self.nu:(i+1)*self.nu] = self.R

    def prediction_varying_matrices(self, signals: SignalsData) -> PredictionModelData:
        '''
        Compute time-varying prediction matrices
        '''
        for i in range(0, self.Np):
            self.REF[i*self.ny:(i+1)*self.ny, :] = signals.ref  #todo use np.tipe(A, rep)

            if i < self.Nc:
                self.U0[i*self.nu:(i+1)*self.nu, :] = signals.u #todo use np.tipe(A, rep)

    @property
    def get_prediction_matrices(self):
        '''
        Getter for PredictionModelData
        '''
        return PredictionModelData(
            A_bar=self.A_bar, B_bar=self.B_bar, C_bar=self.C_bar, D_bar=self.D_bar,
            Q_bar=self.Q_bar, R_bar=self.R_bar,
            REF=self.REF, U0=self.U0,
            X_c_min=self.X_c_min, X_c_max=self.X_c_max,
            U_c_min=self.U_c_min, U_c_max=self.U_c_max,
            dU_c_min=self.dU_c_min, dU_c_max=self.dU_c_max
        )


if __name__ == '__main__':
    print('System prediction')
