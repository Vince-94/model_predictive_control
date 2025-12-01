#!/usr/bin/env python3
'''
Optimization problem (Quadratic Program)
    U_star = min_(u(k)) U(k)^T H U(k) + 2f^T U(k)
                s.t.| A_ineq U(k) <= b_ineq
                    | A_eq U(k) = b_eq
'''
from mpc_dataclass import PredictionModelData, QuadraticProblemData

import numpy as np
import scipy
import quadprog

import logging


class QuadraticProgramOptimization():
    def __init__(self, nx: float, nu: float, Np: float, Nc: float, prediction_matrices: PredictionModelData):
        self.nx = nx
        self.nu = nu
        self.Np = Np
        self.Nc = Nc
        self.A_bar = prediction_matrices.A_bar
        self.B_bar = prediction_matrices.B_bar
        self.C_bar = prediction_matrices.C_bar
        self.D_bar = prediction_matrices.D_bar
        self.Q_bar = prediction_matrices.Q_bar
        self.R_bar = prediction_matrices.R_bar
        self.X_c_min = prediction_matrices.X_c_min
        self.X_c_max = prediction_matrices.X_c_max
        self.U_c_min = prediction_matrices.U_c_min
        self.U_c_max = prediction_matrices.U_c_max
        self.dU_c_min = prediction_matrices.dU_c_min
        self.dU_c_max = prediction_matrices.dU_c_max

        # init optimization matrices
        self.H = np.zeros((self.nu*self.Nc, self.nu*self.Nc))
        self.f = np.zeros((1, self.nu*self.Nc))
        self.A_ineq = np.zeros((2*self.nx*self.Np+4*self.nu*self.Nc, self.nu*self.Nc))
        self.b_ineq = np.zeros(2*self.nx*self.Np+4*self.nu*self.Nc)
        self.U_star = np.zeros((self.nu*self.Nc))

    def qp_matrices(self, x0: np.array, prediction_matrices: PredictionModelData) -> QuadraticProblemData:
        '''
        cost_func = U(k)^T H U(k) + 2f^T U(k) + g
        ineq_constr = A_ineq U(k) <= b_ineq
        '''
        # H = (D_bar^T Q_bar D_bar + R_bar)
        self.H = np.matmul(np.matmul(self.D_bar.T, self.Q_bar), self.D_bar) + self.R_bar

        # f = (C_bar x(k) - REF(k))^T Q_bar D_bar
        self.f = np.matmul(np.matmul((np.matmul(self.C_bar, x0) - prediction_matrices.REF).T, self.Q_bar), self.D_bar)

        # A_ineq = [-B_bar, B_bar, -I_(q Nc), I_(q Nc), -I_(q Nc), I_(q Nc)]^T
        self.A_ineq = np.vstack((
            -self.B_bar,
            self.B_bar,
            -np.identity(self.nu*self.Nc),
            np.identity(self.nu*self.Nc),
            -np.identity(self.nu*self.Nc),
            np.identity(self.nu*self.Nc),
        ))

        # b_ineq = [-X_min+A_bar x(k), X_max+A_bar x(k), -U_min, U_max, -dU_min-U(k-1), dU_max+U(k-1)]^T
        self.b_ineq = np.vstack((
            -self.X_c_min + np.matmul(self.A_bar, x0),
            self.X_c_max - np.matmul(self.A_bar, x0),
            -self.U_c_min,
            self.U_c_max,
            -self.dU_c_min - prediction_matrices.U0,
            self.dU_c_max + prediction_matrices.U0,
        ))

        return QuadraticProblemData(H=self.H, f=self.f, A_ineq=self.A_ineq, b_ineq=self.b_ineq)

    def qp_solve(self, x0: np.array, prediction_matrices: PredictionModelData) -> np.array:
        '''
        U_star = min_(u(k)) U(k)^T H U(k) + 2f^T U(k)
                 s.t.| A_ineq U(k) <= b_ineq
                     | A_eq U(k) = b_eq
        '''
        qp_matrices = self.qp_matrices(x0=x0, prediction_matrices=prediction_matrices)
        U_star = self.qp_solver_quadprog(qp_matrices.H, qp_matrices.f, qp_matrices.A_ineq, qp_matrices.b_ineq)
        return U_star

    #todo in another class
    def qp_solver_quadprog(self, H: np.array, f: np.array, A_ineq: np.array, b_ineq: np.array, A_eq=None, b_eq=None) -> np.array:
        '''
        https://scaron.info/doc/qpsolvers/quadratic-programming.html
        U_star = min_(u(k)) U(k)^T H U(k) - 2f^T U(k)
                    s.t.| A_ineq U(k) >= b_ineq
                        | A_eq U(k) = b_eq

        quadprog.solve_qp(H, f, A_ineq, b_ineq, A_eq, b_eq, solver)
        solvers: ['cvxopt', 'ecos', 'gurobi', 'osqp', 'quadprog', 'scs']
        '''
        f = np.squeeze(f)
        b_ineq = np.squeeze(b_ineq)
        U_star, fun, xu, iters, lagr, iact = quadprog.solve_qp(H, -f, -A_ineq.T, -b_ineq)
        return U_star


    #todo in another class
    def qp_solver_scipy(self, x0: float, H: np.array, f: np.array, A_ineq: np.array, b_ineq: np.array, A_eq=None, b_eq=None, meq=0) -> np.array:
        '''
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html
        U_star = min_(u(k)) U(k)^T H U(k) - 2f^T U(k)
                    s.t.| A_ineq U(k) >= b_ineq
                        | A_eq U(k) = b_eq

        scipy.optimize.minimize(fun, x, method, constraints, tol, options

        )
        '''
        def cost_function(x):
            return np.dot(x, H).dot(x) - np.dot(f, x)

        constraints = [] #todo not needed
        constraints = [
            {
                'type': 'eq' if i < meq else 'ineq',
                'fun': lambda x, A_ineq=A_ineq, b=b_ineq, i=i: (np.dot(A_ineq.T, x) - b)[i]
            } for i in range(A_ineq.shape[1])
        ]

        U_star = scipy.optimize.minimize(
            fun=cost_function, x0=x0, method='SLSQP', constraints=constraints,
            tol=1e-10, options={'maxiter': 2000})

        return U_star



if __name__ == '__main__':
    print('Optimization Problem: Quadratic Program')
