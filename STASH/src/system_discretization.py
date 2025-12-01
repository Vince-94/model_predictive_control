#!/usr/bin/env python3
'''
System model discretization:
    Ad = nx,nx | Bd = nx,nu | Cd = ny,nx | Dd = ny,nu
'''
from mpc_dataclass import SystemModelData

import numpy as np
import logging


class SystemDiscretization():
    def __init__(self, cont_sys_model: SystemModelData, Ts: float):
        self.A = cont_sys_model.A
        self.B = cont_sys_model.B
        self.C = cont_sys_model.C
        self.D = cont_sys_model.D
        self.nx = np.size(self.A, 0)
        self.Ts = Ts

    def euler_discretization(self) -> SystemModelData:
        '''
        Matrices discretization through 'Euler Approximation' method:
        Ad = I + A*Ts
        Bd = B*Ts
        Cd = C*Ts
        Dd = D*Ts
        '''
        self.Ad = np.eye(self.nx, self.nx) + self.A*self.Ts
        self.Bd = self.B*self.Ts
        self.Cd = self.C*self.Ts
        self.Dd = self.D*self.Ts

        return SystemModelData(A=self.Ad, B=self.Bd, C=self.Cd, D=self.Dd)


if __name__ == '__main__':
    print('System discretization')
