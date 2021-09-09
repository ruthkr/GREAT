#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import numpy as np
import os
import gpflow
import subprocess
import sys
import pandas as pd
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import concurrent.futures
from collections import defaultdict
import warnings

# suppress X shape warning generated in logdensities.py
logNightmare = logging.getLogger('gpflow.logdensities')
logNightmare.setLevel(logging.ERROR)

def get_lmls(input_dim, f_var, n_var, length, X, Y, limitsDict):
    # model construction + optimization + get optimal parameters, and start
    # and optimised lml

    # see https://gpflow.readthedocs.io/en/master/_modules/gpflow/kernels.html
    # for kernel sourcecode from which can get params

    # unpack limits dict to constrain the bounds on the search space for
    # hyperparam optimization. Can only set single global limit for all of
    # lengtscales.
    L_bounds = limitsDict['lengthscale']
    Vf_bounds = limitsDict['Vf']
    Vn_bounds = limitsDict['Vn']

    with gpflow.defer_build():
        # set search start values
        k = gpflow.kernels.Matern32(input_dim=input_dim, variance=f_var, lengthscales=length)
        m = gpflow.models.GPR(X, Y, kern=k)
        m.likelihood.variance = n_var

        # set limits on search space, same limits although diff dims can have diff. start values
        m.kern.lengthscales.transform = gpflow.transforms.Logistic(L_bounds[0], L_bounds[1])
        m.kern.variance.transform = gpflow.transforms.Logistic(Vf_bounds[0], Vf_bounds[1])
        # nb limit Vn possibilities to enforce that there is SOME noise - otherwise, seems to really overweight towards 0 noise, flexible functions.
        # enforce the same Vn limits for all gene combos.
        m.likelihood.variance.transform = gpflow.transforms.Logistic(Vn_bounds[0], Vn_bounds[1])

    try:
        m.compile()
        lml = m.compute_log_likelihood()
    except:
        warnings.warn('WARNING: model compilation failed')
        return -999999, f_var, n_var, length, -999999

    # optimization seems can fail when very near boundaries
    # as only seems to occur for some "random" starts, just have safely fail, and
    # calc continue.
    try:
        gpflow.train.ScipyOptimizer().minimize(m)
    except:
        warnings.warn('WARNING: model optimization failed!\nParams: L={}, Vf={}, Vn={}, restarting loop'.format(length, f_var, n_var))

    P = m.read_trainables()
    optVf = P['GPR/kern/variance']
    optVn = P['GPR/likelihood/variance']
    optL = P['GPR/kern/lengthscales']
    optLml = m.compute_log_likelihood()

    return lml, optVf, optVn, optL, optLml

def predict_expression(input_dim, f_var, n_var, length, xx, X, Y):
    # predict expression for plotting ribbon at optimal solutions
    # parameter arguments have already been optimised!

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        k = gpflow.kernels.Matern32(input_dim=input_dim, variance=f_var, lengthscales=length)
        m = gpflow.models.GPR(X, Y, kern=k)
        m.likelihood.variance = n_var
        lml = m.compute_log_likelihood()

    y_pred, var = m.predict_y(xx)

    return y_pred, var

def get_columns(sym, colnames):
    # given a symbol, return all the labels with the gene in
    keepcols = [col for col in colnames if sym in col]
    return keepcols

def get_predictors_and_target(parent_cols, curr_target, df):
    # gets actual values as arrays of correct dimensions for appropriate columns
    # from df master table
    for i, p in enumerate(parent_cols):
        curr_X = df[p].values
        curr_X = curr_X.reshape(len(curr_X), 1)
        if i == 0:
            X = curr_X
        else:
            X = np.concatenate((X, curr_X), axis=1)

    # X and Y have to have 2 dimensions
    Y = df[curr_target].values
    Y = Y.reshape(len(Y), 1)
    if num_parents == 1:
        X = X.reshape((len(X), 1))

    assert X.shape[1] == num_parents

    return(X, Y)

def unpack_synth_info(info_string):
    info = info_string.split('_')
    infoDict = {}

    for i in info:
        k,v = i.split('=')
        infoDict[k] = v
    return(infoDict)

def make_multidim_X(xx, num_parents):
    numObs = len(xx)

    if num_parents == 1:
        return(xx)

    for i in range(num_parents):
        curr_xx = np.linspace(0,1,numObs).reshape(numObs, 1)
        if i == 0:
            xx = curr_xx
        else:
            xx = np.concatenate((xx, curr_xx), axis=1)

    xx = xx.transpose()
    G = np.meshgrid(*xx)
    for i,j in enumerate(G):
        curr_xx = j.ravel().reshape(numObs**num_parents, 1)
        if i == 0:
            XX = curr_xx
        else:
            XX = np.concatenate((XX, curr_xx), axis=1)
    return(XX)

def do_the_work(X, Y, curr_parent, length_bounds, Vf_bounds, Vn_bounds, numStarts):

    numParents = X.shape[1]

    # estimate approximate starting hyperparameter values using the heuristics described
    # https://drafts.distill.pub/gp/
    var_f_est = np.var(Y)  # fine that only one, because only ever 1 target.

    length_raw = np.var(X, axis=0)  # get length estimate seperately for each predictor
    assert len(length_raw) == X.shape[1],'Screwed up length raw estimate shape!'

    # calculate boundaries for optimization of hyperparameters
    # dict of {length : ((lower_p1, lower_p2), (upper_p1, upper_p2))}
    l_b = ( (length_raw/length_bounds[1])-1e-3, (length_raw/length_bounds[0])+1e-3)
    Vf_b = ( (var_f_est/Vf_bounds[1])-1e-3, (var_f_est/Vf_bounds[0])+1e-3)
    # don't set hard boundary limits on Vn optimization

    # boundaries for start values for hyperparameter optimization search
    start_boundaries = {'lengthscale':l_b, 'Vf':Vf_b}

    # boundaries for optimization of hyperparameters
    opt_boundaries = {'lengthscale':( np.min(start_boundaries['lengthscale']),
                                      np.max(start_boundaries['lengthscale']) ),
                      'Vf':( np.min(start_boundaries['Vf']),
                              np.max(start_boundaries['Vf']) ),
                      'Vn':( np.min(Vn_bounds),
                              np.max(Vn_bounds))}

    # setup storage of param scan values
    var_f_start = [0]*numStarts
    var_n_start = [0]*numStarts
    lengthscale_start = [0]*numStarts
    lml_start = [0]*numStarts

    # scan for best hyperparam. values
    seenOptParams = set()  # store opt. param sets, so don't repeatedly report
    best_lmls = []  # -999999 # dummy lml value to start
    i = 0  # loop counter for reperting progress
    sol_num = 0  # the number of log marg likeli minima found

    # do numStarts random starts of the optimization
    #for i in range(1, numStarts):
    for i in range(numStarts):
                print('i:{}'.format(i))
                if i % 10 == 0:
                    print('on iter: ', str(i), '/', str(numStarts))

                # generates random values between the allowed lower and upper
                # bounds (output is same shape as bounds, and each element is
                # in the rannge of the pair in the same position).
                var_f = np.random.uniform(start_boundaries['Vf'][0], start_boundaries['Vf'][1]) # can only be one Vf, irrespective of how many inputs (Y axis has the same variability whichever X-axis look along!)
                # lenscale will be vector of length num parents, w/ start value
                # for each parent, BUT optimization limits are the most extreme
                # scalars.
                lenscale = np.random.uniform(start_boundaries['lengthscale'][0], start_boundaries['lengthscale'][1])
                # enforce that Vn start value must be within the boundaries of Vn search space
                #if (var_n < np.min(opt_boundaries['Vn']) or var_n < np.max(op_boundaries['Vn'])):
                var_n = np.random.uniform(np.min(opt_boundaries['Vn']), np.max(opt_boundaries['Vn'])) # start noise value for optimization

                # record search start parameter space point - used only if 1 parent
                # sanitise lenscale so plotting R script works if 1 parent
                if numParents == 1:
                    lenscale_record = float(lenscale)
                else:
                    lenscale_record = lenscale

                var_f_start[i] = var_f
                var_n_start[i] = var_n
                lengthscale_start[i] = lenscale_record

                # make model @ start+opt hyperparams + record log marginal likelihoods
                # AVOID MEMORY LEAK in GPFlow BY SPAWNING CHILD PROCESS
                start_soln = (np.around(var_f, 2), np.around(var_n, 2), np.around(lenscale, 2))
                with concurrent.futures.ProcessPoolExecutor(max_workers=1) as executor:
                    job = executor.submit(get_lmls, numParents, var_f, var_n, lenscale, X, Y, opt_boundaries)
                    while not job.done():
                        pass
                    lml, optVf, optVn, optL, optlml = job.result()

                # record lml at start point
                lml_start[i] = lml

                # need to be hashable type to add to compare to seenOptParams
                soln = (np.around(optVf, 2),
                        np.around(optVn, 2),
                        tuple(np.around(optL, 2)))

                # if find a new best optimal solution
                if not soln in seenOptParams:
                    sol_num += 1 # update count of unique minima found
                    # add current solution parameters to seen set
                    seenOptParams.add(soln)
                    # make prediction using optimised model for plots

                    # if only one parent make grid of values to make predictions for. Otherwise
                    # probably not too useful to make predictions + visualise.
                    xx = np.linspace(np.min(X), np.max(X), 100).reshape(100, 1)

                    if numParents == 1:
                        with concurrent.futures.ProcessPoolExecutor(max_workers=1) as executor:
                            pred_job = executor.submit(predict_expression, 1, optVf, optVn, optL, xx, X, Y)
                            while not pred_job.done():
                                pass
                            mean, var = pred_job.result()
                    else:
                        mean = np.repeat('NA', xx.shape[0])
                        var = np.repeat('NA', xx.shape[0])

                    # record optimised model results - only actually used if 1 parent only
                    # sanitise L hyperparam output
                    if numParents == 1:
                        L_out = float(soln[2][0])
                    else:
                        L_out = str(soln[2])  # doesn't matter that crap -> table never read by R unless 1 parent only

                    curr_preds = pd.DataFrame.from_dict(
                                    {'Log_marg_lik':np.repeat(optlml, xx.shape[0]),
                                     'Vf':np.repeat(soln[0], xx.shape[0]),
                                     'Vn':np.repeat(soln[1], xx.shape[0]),
                                     'L':np.repeat(L_out, xx.shape[0]),
                                     curr_parent:np.transpose(xx)[0],
                                     'y_pred':np.transpose(mean)[0],
                                     'pred_var':np.transpose(var)[0],
                                     'sol_num':np.repeat(sol_num, xx.shape[0]),
                                     'L_upper':np.repeat(round(opt_boundaries['lengthscale'][1], 3), xx.shape[0]),
                                     'L_lower':np.repeat(round(opt_boundaries['lengthscale'][0], 3), xx.shape[0]),
                                     'Vf_upper':np.repeat(round(opt_boundaries['Vf'][1], 3), xx.shape[0]),
                                     'Vf_lower':np.repeat(round(opt_boundaries['Vf'][0], 3), xx.shape[0]),
                                     'num_parents':np.repeat(numParents, xx.shape[0])
                                    })
                    # append to all other optimized predictions
                    if sol_num == 1:
                        all_preds = curr_preds
                    else:
                        all_preds = all_preds.append(curr_preds)
                    # record all local optima found
                    best_lmls.append(optlml)

    # record the starting grid search values
    param_scans = pd.DataFrame.from_dict(
                                {'var.f.start':var_f_start,
                                 'var.n.start':var_n_start,
                                 'lengthscale.start':lengthscale_start,
                                 'log.marg.likeli.start':lml_start,
                                })

    return (all_preds, param_scans, best_lmls)


#### MAIN :
path_to_normed_data = sys.argv[1]
curr_target = sys.argv[2]
parent_symbol = sys.argv[3] # "," seperated string of candidate parent genes to evaluate
rds = sys.argv[4]

print('**************parent_CSI.py**************')
print('arguments:')
print('path_to_normed_data={}'.format(path_to_normed_data))
print('curr_target={}'.format(curr_target))
print('parent_symbol={}'.format(parent_symbol))
print('rds={}'.format(rds))


# set hyperparameter scanning parameters:
numStarts = 10 # number of random starts for optimisation
# parameters to control upper and lower limits on search space for each hyperparam.
length_bounds = (0.2, 0.4) # bigger values allow smaller lengthscales. if allow too big (min limit=0.1), then can get non-invertible matrix error...
Vf_bounds = (0.5, 5) # bigger values allow smaller function variance
Vn_bounds = (0.01, 1) # the limits for allowed noise for calculating likelihood (forces there to be some noise!).

# the number of regulating genes to consider at once for each target gene
parent_cols = parent_symbol.split(',')
num_parents = len(parent_cols)

# id string will be used for out file names etc.
data_path='T={}_P={}'.format(curr_target, parent_symbol)

# read scaled expression data into panda df
df = pd.read_csv('{}/express_data.csv'.format(path_to_normed_data))
cols=list(df.columns.values)

# get full name of target + check that uniquely IDs
targets = [c for c in cols if curr_target in c]
assert len(targets)==1, "parent_CSI.py: target symbol doesn't uniquely identify a target gene expression column! {} in {}".format(curr_target, cols)
curr_target = targets[0]

# get Y (target) and X (predictor values), need to be array, second value of X
# should be the number of parents
X, Y = get_predictors_and_target(parent_cols, curr_target, df)

# record the true values for plotting latter compared to fit GPs.
measured_vals = df[parent_cols+[curr_target]+['sample_id']+['group']]

# for current parent and target, search GP hyperparameter space and evaluate
# log marginal likelihood.
# Returns the optimal solutions found from each random start,
# the solutions at the start points of the optimization,
# and a vector of just the best lmls found.
all_preds, param_scans, best_lmls = do_the_work(X,
                                                Y,
                                                parent_symbol.replace(',', '+'),
                                                length_bounds,
                                                Vf_bounds,
                                                Vn_bounds, numStarts)


# make directory where will write output files if doesn't exist
final_dir = './lml_files/'
if not os.path.exists(final_dir):
    try:
        os.makedirs(final_dir)
    except FileExistsError:
        pass


# write current target:parent all the diff optimised scores to file
# get the details of the synth data run
# if synthetic data then want to record info about params to gen data too.
outPath = '{}{}_best_lml.tsv'.format(final_dir, data_path)
print('writing best lmls to ...')
print(outPath)

with open(outPath, 'w') as F:
    for curr_best_lml in best_lmls:
        F.write('{}\t{}\t{}\t{}\n'.format(curr_target,
                                    parent_symbol.replace(',','+'),
                                    curr_best_lml,
                                    num_parents) )


#  cleanup intermediate files for current parent:target set generated by this script
print('Cleaning up tmp files...')
cmd = 'rm {}{}*'.format(path_to_normed_data, data_path)
print(cmd)
subprocess.call(cmd, shell=True)

print('************end parent_CSI.py************')
