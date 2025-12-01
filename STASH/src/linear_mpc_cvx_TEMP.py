#!/usr/bin/env python3
'''
'''
from mpc_dataclass import SystemModelData, PredictionModelData

import numpy as np
import scipy
import cvxpy

import logging





class LinearMPC(SystemPrediction):
    def __init__(self):

        # init optimization matrices
        self.H = np.zeros((self.nu*self.Nc, self.nu*self.Nc))
        self.f = np.zeros((1, self.nu*self.Nc))

        self.A_ineq = np.zeros((2*self.nx*self.Np+4*self.nu*self.Nc, self.nu*self.Nc))
        self.b_ineq = np.zeros((2*self.nx*self.Np+4*self.nu*self.Nc, 1))

        self.U_star = np.zeros((self.nu*self.Nc, 1))


    def qp_solver(self):
        # self.U_star = np.zeros((self.nu*self.Nc, 1))
        self.opt = {'disp':True}
        self.U_star = scipy.optimize.minimize(self.loss, self.x0, jac=self.jac, method='SLSQP', options=self.opt)

        self.u_star = self.U_star[0:self.nu]

        # sysyem evolution
        self.x = 0 # x = f(x0, u)

        # update past state and output
        self.x0 = self.x
        self.u0 = self.u_star


    def loss(self, x):
        return np.dot(x.T, np.dot(self.H, x)) + np.dot(self.f, x)


    def jac(self, x):
        return np.dot(x.T, self.H) + self.f



    # def qp_solver(self, x0):
    #     '''
    #     '''
    #     opt = {'disp':True}
    #     U_star = scipy.optimize.minimize(self.loss, x0, jac=self.jac, method='SLSQP', options=opt)

    #     u_star = U_star[0:self.nu]

    #     return u_star


    def qp_solve_cvxpy(self):
        '''
        U_k = min_(u(k)) U(k)^T H U(k) + 2f^T U(k) + g
        s.t.| A_ineq U(k) <= b_ineq
        '''

        # x = cp.Variable(n)
        # prob = cp.Problem(cp.Minimize((1/2)*cp.quad_form(x, P) + q.T @ x),
        #                 [G @ x <= h,
        #                 A @ x == b])
        # prob.solve()


        u = cvxpy.Variable(self.Nc)

        obj_func = cvxpy.Minimize(cost)
        qp_problem = cvxpy.Problem(obj_func, constr)

        if cvxpy.CVX == True:
            qp_problem.solve(verbose=cvxpy.verbose, solver=cvxpy.ECOS) # I find that ECOS is better please use it when solving QPs
        else:
            qp_problem.solve(verbose=cvxpy.verbose)


