# Author: Katiana Kontolati, Ph.D. candidate
# Department of Civil and Systems Engineering, Johns Hopkins University
# Last update: March 25, 2021

#######################################################################################################################
#######################################################################################################################
#                                    Grassmannian Diffusion Maps PCE surrogates                                       #
#######################################################################################################################
#######################################################################################################################

# Import all necessary libraries
from random import randrange
import numpy as np
import math
from mpl_toolkits.mplot3d import Axes3D
from UQpy.Distributions import Normal, Uniform, JointInd
from sklearn.model_selection import train_test_split
from UQpy.Surrogates import *
from DimensionReduction import Grassmann
from DimensionReduction import DiffusionMaps
from UQpy.SampleMethods import LHS
from scipy.integrate import odeint

from skopt.space import Real, Categorical, Integer
from skopt.searchcv import BayesSearchCV
# to install skopt run: $ pip install scikit-optimize

import matplotlib.pyplot as plt
from matplotlib import cm

from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import KMeans

import datafold.pcfold as pfold
from datafold.dynfold import GeometricHarmonicsInterpolator as GHI

import random
import itertools as it
import os, subprocess
from matplotlib.ticker import MaxNLocator
import time
import sys

sys.path.append('./data/')
from DiffusionEquation import diffusion
from electric_potential import function
from rand_cmap import rand_cmap
from LoktaVoltera import LV

#######################################################################################################################
#                                            Step 1: Create dataset                                                   #
#######################################################################################################################

# Provided 


#######################################################################################################################
#                                         Step 2: Grassmannian Diffusion Maps                                         #
#######################################################################################################################

class GDMaps:
    """
    Performs GDMaps for a given dataset.
    n_evecs must be greater then n_parsim

    """

    def __init__(self, data, n_evecs, n_parsim, p, verbose=False):
        self.data = data
        self.n_evecs = n_evecs
        self.n_parsim = n_parsim
        self.p = p
        self.verbose = verbose

    def get(self):
        Gr = Grassmann(distance_method=Grassmann.grassmann_distance, kernel_method=Grassmann.projection_kernel,
                       karcher_method=Grassmann.gradient_descent)
        Gr.manifold(p=self.p, samples=self.data)

        dfm = DiffusionMaps(alpha=0.5, n_evecs=self.n_evecs + 1, kernel_object=Gr, kernel_grassmann='prod')
        g, evals, evecs = dfm.mapping()

        # Parsimonious representation
        index, residuals = dfm.parsimonious(num_eigenvectors=self.n_evecs, visualization=False)

        coord = index[1:self.n_parsim + 1]

        g_k = g[:, coord]

        # g_k = g[:, 1:]  # without parsimonious
        # coord = np.arange(1, g_k.shape[1]+1)  # diffusion coordinates numbers

        print('Grassmann projection rank is: ', Gr.p)

        return g_k, coord, Gr, residuals, index, evals


def plot_diff_coord(x, data, coord, labels):
    """
    Plots the diffusion coordinates from the GDMaps.

    """
    plt.rcParams.update({'font.size': 24})

    nlabels = np.unique(labels).shape[0]
    cmap = rand_cmap(nlabels=nlabels, type='bright', first_color_black=False)

    comb1 = list(it.combinations(list(coord), 2))
    comb2 = list(it.combinations([i for i in range(coord.shape[0])], 2))

    if os.path.exists('figures'):
        command = ['rm', '-r', 'figures']
        subprocess.run(command)

    command = ['mkdir', 'figures']
    subprocess.run(command)

    for i in range(len(comb1)):
        plt.figure(figsize=(8, 6), constrained_layout=True)
        plt.scatter(data[:, comb2[i][0]], data[:, comb2[i][1]], s=30, c=labels, cmap=cmap)
        plt.xlabel(r'$\psi_{}$'.format(comb1[i][0]), fontsize=26)
        plt.ylabel(r'$\psi_{}$'.format(comb1[i][1]), fontsize=26)
        plt.grid(True)
        plt.savefig('figures/Psi_{},{}.png'.format(comb1[i][0], comb1[i][1]), bbox_inches='tight')

    # Plot first three plots
    if coord.shape[0] > 2:
        fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(22, 5), constrained_layout=True)
        for i in range(3):
            ax[i].scatter(data[:, comb2[i][0]], data[:, comb2[i][1]], s=30, c=labels, cmap=cmap)
            ax[i].set_xlabel(r'$\psi_{}$'.format(comb1[i][0]), fontsize=28)
            ax[i].set_ylabel(r'$\psi_{}$'.format(comb1[i][1]), fontsize=28)
            ax[i].grid('True')
            ax[i].ticklabel_format(style='sci', axis='both', scilimits=(0, 0))
            # plt.legend()
            # ax[i].set_title('Training realizations: {}'.format(trunc[i]))
        plt.savefig('figures/Diffusion-coord.png', bbox_inches='tight', dpi=300)

    fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)
    plt.scatter(x[:, 0], x[:, 1], c=labels, cmap=cmap)
    plt.xlabel(r'$x_1$', fontsize=22)
    plt.ylabel(r'$x_2$', fontsize=22)
    plt.title('Input parameters colored by \n the clusters on diffusion manifold')
    plt.savefig('figures/stochastic-inputs.png', bbox_inches='tight', dpi=300)


#######################################################################################################################
#                                                 Step 3: PCE surrogate                                               #
#######################################################################################################################

class PceModel:
    """
    Constructs a PCE surrogate on the Grassmannian diffusion manifold.

    """

    def __init__(self, x, g, dist_obj, max_degree, verbose=False):
        self.x = x
        self.g = g
        self.dist_obj = dist_obj
        self.max_degree = max_degree
        self.verbose = verbose

    def get(self):

        # Polynomial basis
        dim_in = self.x.shape[1]
        polys = Polynomials(dist_object=self.dist_obj, degree=self.max_degree)

        n_basis = math.factorial(self.max_degree + dim_in) / \
                  (math.factorial(self.max_degree) * math.factorial(dim_in))
        if self.verbose:
            print('Basis terms: ', int(n_basis))

        # Regression method
        reg = PolyChaosLstsq(poly_object=polys)
        # reg = PolyChaosLasso(poly_object=polys, learning_rate=0.001, iterations=1000, penalty=0.05)
        # reg = PolyChaosRidge(poly_object=polys, learning_rate=0.001, iterations=10000, penalty=0)

        pce = PCE(method=reg)

        x_train, x_test, \
        g_train, g_test = train_test_split(self.x, self.g, train_size=2 / 3, random_state=1)

        # Design matrix / conditioning
        D = polys.evaluate(self.x)
        cond_D = np.linalg.cond(D)
        if self.verbose:
            print('Condition number: ', cond_D)

        # Fit model
        pce.fit(x_train, g_train)

        error_val = ErrorEstimation(surr_object=pce).validation(x_test, g_test)

        if self.verbose:
            # Plot accuracy of PCE
            if os.path.exists('pce_accuracy'):
                command = ['rm', '-r', 'pce_accuracy']
                subprocess.run(command)

            command = ['mkdir', 'pce_accuracy']
            subprocess.run(command)

            print(g_test[0, :])
            print(pce.predict(x_test)[0, :])

            for i in range(5):
                r = random.randint(0, x_test.shape[0])
                plt.figure()
                plt.plot(g_test[r, :], 'b-o', label='true')
                plt.plot(pce.predict(x_test)[r, :], 'r-*', label='pce')
                plt.legend()
                plt.savefig('pce_accuracy/pce_{}.png'.format(i), bbox_inches='tight')

        return pce, error_val


#######################################################################################################################
#                                     Step 4: Adaptive clustering                                                     #
#######################################################################################################################

def AdaptClust(n_clust_max, Gr_object, data):
    """
    Adaptive clustering to find the optimal number of clusters of the data onto the diffusion manifold

    """
    Gr = Gr_object
    n_clust = 1
    n_clust_max = n_clust_max
    clusters, error, mat_all, ind_all, kmeans_models, L_all = [0], [], [], [], [], []

    while clusters[-1] < n_clust_max:

        if n_clust == 1:
            clusters.pop(-1)

        n_clust += 1

        # K-means clustering
        kmeans = KMeans(n_clusters=n_clust, random_state=0).fit(data)
        C, L = kmeans.cluster_centers_, kmeans.labels_

        # Get indices of clusters
        indices = []
        cluster_num = [i for i in range(n_clust)]
        for i in cluster_num:
            indices.append(np.where(L == i)[0])

        ind_all.append(indices)
        kmeans_models.append(kmeans)
        clusters.append(n_clust)
        L_all.append(L)

        n = np.array([ind_all[-1][i].shape[0] for i in range(len(ind_all[-1]))])
        if np.any(n < 5):
            ind_all.pop(-1)
            kmeans_models.pop(-1)
            clusters.pop(-1)
            L_all.pop(-1)
            print('A cluster of less than 5 points was detected. The algorithm stopped.')
            break

        # Compute psi and phi matrices
        acc, mat_ = [], []
        for k in range(n_clust):
            psi = np.array([Gr.psi[indices[k][i]] for i in range(len(indices[k]))])
            phi = np.array([Gr.phi[indices[k][i]] for i in range(len(indices[k]))])

            # Compute Karcher mean of points
            karcher_psi = Gr.karcher_mean(points_grassmann=psi, acc=False, tol=1e-3, maxiter=1000)
            karcher_phi = Gr.karcher_mean(points_grassmann=phi, acc=False, tol=1e-3, maxiter=1000)

            tan_psi = Gr.log_map(points_grassmann=psi, ref=karcher_psi)
            tan_phi = Gr.log_map(points_grassmann=phi, ref=karcher_phi)

            back_psi = Gr.exp_map(points_tangent=tan_psi, ref=karcher_psi)
            back_phi = Gr.exp_map(points_tangent=tan_phi, ref=karcher_phi)

            # Mean-squared error (MSE)
            mean_psi = np.mean(
                [(np.square(psi[i] - back_psi[i])).mean(axis=None) for i in range(indices[k].shape[0])])
            mean_phi = np.mean(
                [(np.square(phi[i] - back_phi[i])).mean(axis=None) for i in range(indices[k].shape[0])])

            mean_all = [mean_psi, mean_phi]
            acc.append(np.mean(mean_all))

            mat_.append([tan_psi, tan_phi, karcher_psi, karcher_phi])

        mat_all.append(mat_)
        error.append(np.mean(acc))

        # print('Error for {} clusters:'.format(n_clust), error[-1])

        if n_clust > 2:
            imp = (error[-2] - error[-1]) / error[-2]
            if imp < 0:
                error.pop(-1)
                clusters.pop(-1)
                mat_all.pop(-1)
                ind_all.pop(-1)
                kmeans_models.pop(-1)
                L_all.pop(-1)

    mat = mat_all[-1]
    indices = ind_all[-1]
    kmeans = kmeans_models[-1]
    L = L_all[-1]
    n_clust = len(mat)

    return mat, indices, kmeans, L, n_clust, error, clusters


#######################################################################################################################
#                                  Step 5: Geometric harmonics and PCE interpolators                                  #
#######################################################################################################################


def Interpolators(x, data, mat, indices, n_clust, Gr, joint):
    """
    Constructs a GHI model to find a map between SVD matrices of knn points and diffusion coordinates
     on the GDMaps manifold.

    """

    # Create GH and PCE models
    random_state = 1
    models_all = []

    for k in range(n_clust):

        models = []  # save GH for psi and phi for each cluster
        tan_psi, tan_phi, karcher_psi, karcher_phi = mat[k]

        # Convert lists to numpy arrays
        m1, m2 = np.array(tan_psi), np.array(tan_phi)
        l, m, n = m1.shape
        m1, m2 = m1.reshape(l, m * n), m2.reshape(l, m * n)

        M = [m1, m2]
        for i in range(2):
            p_train, p_test, \
            g_train, g_test = train_test_split(M[i], data[indices[k]],
                                               train_size=2 / 3, random_state=random_state)

            pcm = pfold.PCManifold(g_train)
            pcm.optimize_parameters(random_state=random_state)

            train_indices, test_indices = train_test_split(np.random.permutation(p_train.shape[0]),
                                                           train_size=2 / 3, test_size=1 / 3)

            # GH training (no Bayesian optimization)
            opt_epsilon = pcm.kernel.epsilon
            opt_cutoff = pcm.cut_off
            opt_n_eigenpairs = train_indices.shape[0] - 1

            # test the interpolation quality with PCManifold optimization
            optimal_GHI = GHI(pfold.GaussianKernel(epsilon=opt_epsilon),
                              n_eigenpairs=opt_n_eigenpairs,
                              dist_kwargs=dict(cut_off=opt_cutoff))

            optimal_GHI.fit(g_train[train_indices, :], p_train[train_indices, :])

            # Get error and residual
            # residual = optimal_GHI.score(g_train, p_train)
            # error = optimal_GHI.score(g_test, p_test)

            models.append(optimal_GHI)

        models_all.append(models)

    dims = (m, n)
    # Construct PCE models for sigmas
    sigmas = np.array(Gr.sigma)

    max_degree = 2
    polys = Polynomials(dist_object=joint, degree=max_degree)
    reg = PolyChaosLstsq(poly_object=polys)
    pce = PCE(method=reg)
    # Design matrix / conditioning
    D = polys.evaluate(x)
    cond_D = np.linalg.cond(D)
    # print('Condition number: ', cond_D)
    # Fit model
    pce.fit(x, sigmas)

    error_val2 = ErrorEstimation(surr_object=pce).validation(x, sigmas)
    print('Validation error of PCE of sigmas: ', error_val2)

    return models_all, pce, dims


#######################################################################################################################
#                                          Step 6: Out-of-sample predictions                                          #
#######################################################################################################################


def Prediction(x_pred, y_real, models_all, kmeans, mat, pce, pce_sigmas, Gr, dims):
    """
    Out-of-sample predictions by using the GH and pce models of the clusters on the diffusion manifold

    """

    # In case of one global pce
    # pce = pce[0]

    m, n = dims  # dimensions of points onto the Grassmann manifold
    y_recon, l2, r2, diff = [], [], [], []
    num = x_pred.shape[0]

    for k in range(num):

        # for one global PCE
        x_new = x_pred[k, :]  # new sample
        y_pce = pce.predict(x_new.reshape(1, -1))
        l = int(kmeans.predict(y_pce.reshape(1, -1)))  # predicted label

        # For multiple PCEs
        # print('iteration of predictions: {}'.format(k))
        # x_new = x_pred[k, :]  # new sample
        # get_index = knn.kneighbors(x_new.reshape(1, -1), return_distance=False)
        # get = int(kmeans.labels_[get_index])
        # print(get)
        # y_pce = pce[get].predict(x_new.reshape(1, -1))

        # l = int(kmeans.predict(y_pce.reshape(1, -1)))  # predicted label
        # print(l)
        # print('')

        # Return to ambient space
        psi_tan_point = models_all[l][0].predict(y_pce.reshape(1, -1))
        phi_tan_point = models_all[l][1].predict(y_pce.reshape(1, -1))
        sigma_grass = pce_sigmas.predict(x_new.reshape(1, -1))
        sigma_grass = np.diag(sigma_grass.flatten())

        # Project psi and phi back on Grassmann
        back_psi = np.squeeze(np.array(Gr.exp_map(points_tangent=[psi_tan_point.reshape(m, n)],
                                                  ref=mat[l][2])))
        back_phi = np.squeeze(np.array(Gr.exp_map(points_tangent=[phi_tan_point.reshape(m, n)],
                                                  ref=mat[l][3])))

        # Reconstruct sample with reverse SVD
        y_recon.append(np.dot(np.dot(back_psi, sigma_grass), back_phi.T))

        # Relative L2 error
        error = np.linalg.norm(y_real[k] - y_recon[k]) / np.linalg.norm(y_real[k])
        l2.append(error)

        # Coefficient of determination, aka R2 score
        mean_ref = np.mean(y_real[k])
        r2.append(1 - (np.sum((y_recon[k].ravel() - y_real[k].ravel()) ** 2))
                  / (np.sum((y_real[k].ravel() - mean_ref) ** 2)))

        # Relative error between reference and prediction
        diff.append(np.abs((y_real[k] - y_recon[k]) / y_real[k]))

    return y_recon, l2, r2, diff


if __name__ == '__main__':

    print('Nothing on main')