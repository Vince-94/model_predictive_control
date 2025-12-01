#!/usr/bin/env python3

import numpy as np

from mpc_dataclass import SignalsData, SystemModelData


class SisoLtiModel():
    '''
    v = a*t
    v = v0 + x*dt -> x_dot = x + u*dt
    y = v0 -> y = x
    '''
    def __init__(self, Ts):
        self.dt = Ts
        self.A = np.array([1], ndmin=2)
        self.B = np.array([self.dt], ndmin=2)
        self.C = np.array([1], ndmin=2)
        self.D = np.array([0], ndmin=2)

    @property
    def get_model(self):
        return SystemModelData(A=self.A, B=self.B, C=self.C, D=self.D)

    def on_step(self, x: np.array, u: np.array):
        x_dot = np.matmul(self.A, x) + np.matmul(self.B, u)
        y = np.matmul(self.C, x) + np.matmul(self.D, u)
        return SignalsData(x_dot=x_dot, y=y)
        # return x_dot, y


class MimoLtiModel():
    def __init__(self, Ts):
        self.dt = Ts




    @property
    def get_model(self):
        self.A = np.array([
            [0.22, 0.44],
            [0, 0.88]])
        self.B = np.array([
            [1],
            [1]])
        self.C = np.array([
            [1, 0]])
        self.D = np.array([
            [0]])
        return SystemModelData(A=self.A, B=self.B, C=self.C, D=self.D)


    def on_step(self, x: np.array, u: np.array):
        x_dot = np.matmul(self.A, x) + np.matmul(self.B, u)
        y = np.matmul(self.C, x) + np.matmul(self.D, u)
        return SignalsData(x_dot=x_dot, y=y)



class LinearTimeInvariantModel():

    def __init__(self) -> None:
        self.A = np.array([
            [0, 0, 1, 0],
            [0, 0, 0, 1],
            [0, 0, 0, 0],
            [0, 0, 0, 0]])
        self.B = np.array([
            [1, 0],
            [0, 1],
            [1, 0],
            [0, 1]])
        self.C = np.array([
            [1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, 0]])
        self.D = np.array([
            [0, 0],
            [0, 0],
            [0, 0]])

    @property
    def get_model(self):
        return SystemModelData(A=self.A, B=self.B, C=self.C, D=self.D)

    def on_step(self, x: np.array, u: np.array) -> SignalsData:
        x_dot = np.matmul(self.A, x) + np.matmul(self.B, u)
        y = np.matmul(self.C, x) + np.matmul(self.D, u)
        return SignalsData(x_dot=x_dot, y=y)


if __name__ == '__main__':
    print('Linear Model Predictive Control')
