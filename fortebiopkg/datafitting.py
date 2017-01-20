import numpy as np
import math
import scipy.optimize
import matplotlib.pyplot as pl
import os
import warnings
from io import BytesIO as BIO


def __init__(self):
    """ No parameters required """
    pass


def linearfit(baselinex, baseliney, onlylast20percent=False):
    yintercept = 0
    sumxy = 0
    sumx = 0
    sumy = 0
    sumxsq = 0
    baselinex.tolist()
    baseliney.tolist()
    if onlylast20percent:
        i = int(len(baselinex) / 5) * 2
    else:
        i=0
    count = len(baselinex)
    while i < len(baselinex):
        if not np.isfinite(baselinex[i]) or not np.isfinite(baseliney[i]):
            count -= 1
            i += 1
            continue
        sumxy += baselinex[i] * baseliney[i]
        sumx += baselinex[i]
        sumy += baseliney[i]
        sumxsq += baselinex[i] * baselinex[i]
        i += 1
    slope = ((count * sumxy) - (sumx * sumy)) / ((count * sumxsq) - (sumx * sumx))
    yintercept = baseliney[1] - (slope * baselinex[1])
    return [slope, yintercept]


def fit1to1modeldecon(experiment, associationstepnumber=1, dissociationstepnumber=2, globalfit=False, *args, **kwargs):
    forcedecon = False
    verbose = False
    if 'force_decon' in kwargs.keys():
        forcedecon = kwargs['force_decon']
    if 'verbose' in kwargs.keys():
        verbose = kwargs['verbose']

    if not globalfit:
        np.seterr(all="ignore")
        warnings.filterwarnings("ignore")
        count = 0
        totalexp = len(experiment)
        for exp in experiment:
            percentcounter = 100 * (count + 1) / totalexp
            if verbose:
                print "\r1 to 1 model (decon) fitting progress: " + str(percentcounter) + "%",

            if exp['step_status'][-1].upper() is 'ERROR' or \
                    (exp['signal2noise'] is not None and float(exp['signal2noise']) <= 2.5):
                exp['param_name'] = ['Rmax', 'kd', 'ka', 'Cpro', 'm', 'c', 'X0', 'kds']
                exp['fit_param'] = [1, 1, 1, 1, 0, 0, exp['x_data'][dissociationstepnumber][0], 1]
                exp['param_error'] = [1E-12, 1E-12, 1E-12, 1E-12, 1E-12, 1E-12, 1E-12, 1E-12]
                exp['fit_fxn'] = [0, 1]
                if exp['step_status'][-1].upper() is 'ERROR':
                    exp['comments'].append('Data not fit: Sensor Error')
                elif (exp['signal2noise'] is not None and float(exp['signal2noise']) <= 2.5):
                    exp['comments'].append('Data not fit: s/n below threshold')
                paramcheck([exp], associationstepnumber, dissociationstepnumber)
                count += 1
                continue

            if exp['molarconcentration'][associationstepnumber] is None or (
                        float(exp['molarconcentration'][associationstepnumber]) < 0):
                exp['molarconcentration'][associationstepnumber] = 100
            templist = exp['x_data']
            data_x1 = np.array(templist[associationstepnumber])
            data_x2 = np.array(templist[dissociationstepnumber])
            templist = exp['y_data']
            data_y1 = np.array(templist[associationstepnumber])
            data_y2 = np.array(templist[dissociationstepnumber])
            # Get initial parameters
            fit1to1model([exp], 1, 2, False, partial_fit=True, verbose=False)
            p_partial = exp['fit_param']
            p_best_obs = exp['fit_param']
            # initial parameters: 'Rmax', 'kd', 'ka', 'Cpro', 'm', 'c', 'X0', 'kds'
            MCO = (-0.01146 * np.exp(-1*(data_x2[-1] - data_x2[0]) / 77.5588)) + (
                -0.07192 * np.exp(-1*(data_x2[-1] - data_x2[0]) / 10.3))
            exp['param_error'] = [1E-12, 1E-12, 1E-12, 1E-12, 1E-12, 1E-12, 1E-12, 1E-12]

            if forcedecon:
                p_fix = [None, None, None, (float(exp['molarconcentration'][associationstepnumber]) * 1E-9),
                         MCO*1.5, None, float(data_x2[0]), 1]
                p_global = [data_y2.max() * 1.3, p_best_obs[1]/5, p_best_obs[2], (data_y2[0]-data_y2[-1])/2]
                p_best, pcov, infodict, errmsg, success = scipy.optimize.leastsq(residualscalc_analytical1to1model,
                        p_global, args=(data_x1, data_y1, data_x2, data_y2, p_fix),
                            full_output=1, epsfcn=0.0001, maxfev=20000)
                p_fix = [None, p_best[1], None, (float(exp['molarconcentration'][associationstepnumber]) * 1E-9),
                         None, None, float(data_x2[0]), 1]
                p_global = [data_y2.max() * 1.3, p_best[2], MCO*1.5, p_best[3]]
                p_best, pcov, infodict, errmsg, success = scipy.optimize.leastsq(residualscalc_analytical1to1model,
                    p_global, args=(data_x1, data_y1, data_x2, data_y2, p_fix),
                    full_output=1, epsfcn=0.0001, maxfev=20000)
                p_fix = [None, None, None, (float(exp['molarconcentration'][associationstepnumber]) * 1E-9),
                         None, None, float(data_x2[0]), 1]
                p_global = [data_y2.max() * 1.3, p_best_obs[1], p_best[1], p_best[2], p_best[3]]
                p_best, pcov, infodict, errmsg, success = scipy.optimize.leastsq(residualscalc_analytical1to1model,
                    p_global, args=(data_x1, data_y1, data_x2, data_y2, p_fix),
                                                                             full_output=1, epsfcn=0.0001, maxfev=20000)
                p_fix = [None, p_best[1], None, (float(exp['molarconcentration'][associationstepnumber]) * 1E-9),
                         p_best[3], p_best[4], float(data_x2[0]), None]
                p_global = [data_y2.max() * 1.3, p_best[2], .8]
                p_best, pcov, infodict, errmsg, success = scipy.optimize.leastsq(residualscalc_analytical1to1model,
                    p_global, args=(data_x1, data_y1, data_x2, data_y2, p_fix),
                    full_output=1, epsfcn=0.0001, maxfev=20000)
                exp['fit_param'] = assemblefitparam(p_best, p_fix).tolist()
                exp['flags'].append("Override: Decon-Fit enforced")
                exp['comments'].append("Decon-Fit used")
                exp['param_error'] = np.average([calcparamerror(p_best, p_fix, data_x1, data_y1, 0, pcov),
                                                 calcparamerror(p_best, p_fix, data_x2, data_y2, 1, pcov)],
                                                axis=0)
                count += 1
                continue

            # if full fit required do to low amplitude change
            if (np.abs(data_y2[0] - data_y2[-1])/np.abs(data_y1[0] - data_y1[-1]) < 0.2):
                p_best_obs[1] = 1/(2*data_x2.max())
                if (np.abs((data_y2[0]) - data_y2[len(data_y2)/2])) / np.abs(data_y1[-1] - data_y1[0]) <= 0.05:
                    p_fix = [None, (-1*np.log(0.95)/(data_x2[-1]-data_x2[0])), None,
                         (float(exp['molarconcentration'][associationstepnumber])*1E-9), 0.0, 0.0, float(data_x2[0]), 1]
                    p_global = [data_y2.max() * 1.3, p_best_obs[2]]
                    exp['comments'].append(
                        "5% rule invoked(kd fixed as: " + str(-1 * np.log(0.95) / (data_x2[-1] - data_x2[0])) + ")")
                else:
                    p_fix = [None, None, None,
                             (float(exp['molarconcentration'][associationstepnumber]) * 1E-9), 0.0, 0.0,
                             float(data_x2[0]), 1]
                    p_global = [data_y2.max() * 1.3, (-1*np.log(0.95)/(data_x2[-1]-data_x2[0])), p_best_obs[2]]
                exp['flags'].append("Avid or aggregation like binding detected!")
                p_best, pcov, infodict, errmsg, success = scipy.optimize.leastsq(residualscalc_analytical1to1model,
                                p_global, args=(data_x1, data_y1, data_x2, data_y2, p_fix),
                                full_output=1, epsfcn=0.0001, maxfev=20000)
                if p_best[1] < (-1 * np.log(0.95) / (data_x2[-1] - data_x2[0])):
                    p_fix = [None, (-1 * np.log(0.95) / (data_x2[-1] - data_x2[0])), None,
                             (float(exp['molarconcentration'][associationstepnumber]) * 1E-9),
                             0, 0, float(data_x2[0]), 1]
                    p_global = [data_y2.max() * 1.3, 2 / p_fix[3]]
                    p_best, pcov, infodict, errmsg, success = scipy.optimize.leastsq(residualscalc_analytical1to1model,
                        p_global, args=(data_x1, data_y1, data_x2, data_y2, p_fix),
                            full_output=1, epsfcn=0.0001, maxfev=20000)
                    exp['flags'].append(
                        "kd limit exceeded(slower than: " + str(-1 * np.log(0.95) / (data_x2[-1] - data_x2[0])) + ")")
                exp['fit_param'] = assemblefitparam(p_best, p_fix).tolist()
                exp['param_error'] = np.average([calcparamerror(p_best, p_fix, data_x1, data_y1, 0, pcov),
                                                 calcparamerror(p_best, p_fix, data_x2, data_y2, 1, pcov)], axis=0)
                exp['fit_fxn'] = [0, 1]
                exp['comments'].append("Full Fit used")
                count += 1
                continue
            # else if full fit required due to majority loss in amplitude
            elif (np.abs(data_y2[0] - data_y2[-1]))/np.abs(data_y1[-1] - data_y1[0]) >= 0.85 and data_y2[-1] >= 0:
                p_fix = [None, None, None, (float(exp['molarconcentration'][associationstepnumber]) * 1E-9),
                0.0, 0.0, float(data_x2[0]), None]
                p_global = [data_y2.max() * 1.3, p_best_obs[1], p_best_obs[2], 1]
                p_best, pcov, infodict, errmsg, success = scipy.optimize.leastsq(residualscalc_analytical1to1model,
                        p_global, args=(data_x1, data_y1, data_x2, data_y2, p_fix),
                        full_output=1, epsfcn=0.0001, maxfev=20000)
                modeloff = analytical1to1modeloff(data_x2, assemblefitparam(p_best, p_fix))
                if np.max(modeloff - data_y2) / (data_y2[0] - data_y1[0]) <= 0.05 or calcrsq(data_y2, modeloff) >= 0.98:
                    #if rsq or max residual is too high partial fit is required. FXN passes through if full fit is not acceptable
                    exp['comments'].append("Full Fit used")
                    exp['fit_param'] = assemblefitparam(p_best, p_fix).tolist()
                    exp['param_error'] = np.average([calcparamerror(p_best, p_fix, data_x1, data_y1, 0, pcov),
                                                     calcparamerror(p_best, p_fix, data_x2, data_y2, 1, pcov)], axis=0)
                    exp['fit_fxn'] = [0, 1]
                    count += 1
                    continue
            #decide between partial and lin aaprox
            # fit partial first
            p_fix = [None, None, None, (float(exp['molarconcentration'][associationstepnumber]) * 1E-9),
                     0.0, None, float(data_x2[0]), 1]
            p_global = [data_y2.max() * 1.3, p_best_obs[1], p_best_obs[2], data_y2[-1]]
            p_best, pcov, infodict, errmsg, success = scipy.optimize.leastsq(residualscalc_analytical1to1model,
                        p_global, args=(data_x1, data_y1, data_x2, data_y2, p_fix),
                        full_output=1, epsfcn=0.0001, maxfev=20000)
            exp['fit_param'] = assemblefitparam(p_best, p_fix).tolist()
            modeloff = analytical1to1modeloff(data_x2, exp['fit_param'])
            if np.max(modeloff-data_y2)/(data_y2[0]-data_y1[0]) >= 0.05 or calcrsq(data_y2, modeloff) <= 0.98:
                # fit partial + linear
                p_fix = [None, None, None, (float(exp['molarconcentration'][associationstepnumber]) * 1E-9),
                         MCO * 1.5, None, float(data_x2[0]), 1]
                p_global = [data_y2.max() * 1.3, p_best_obs[1] / 5, p_best_obs[2], (data_y2[0] - data_y2[-1]) / 2]
                p_best, pcov, infodict, errmsg, success = scipy.optimize.leastsq(residualscalc_analytical1to1model,
                    p_global, args=(data_x1, data_y1, data_x2, data_y2, p_fix),
                                                                             full_output=1, epsfcn=0.0001, maxfev=20000)
                p_fix = [None, p_best[1], None, (float(exp['molarconcentration'][associationstepnumber]) * 1E-9),
                     None, None, float(data_x2[0]), 1]
                p_global = [data_y2.max() * 1.3, p_best[2], MCO * 1.5, p_best[3]]
                p_best, pcov, infodict, errmsg, success = scipy.optimize.leastsq(residualscalc_analytical1to1model,
                    p_global, args=(data_x1, data_y1, data_x2, data_y2, p_fix),
                                                                             full_output=1, epsfcn=0.0001, maxfev=20000)
                p_fix = [None, None, None, (float(exp['molarconcentration'][associationstepnumber]) * 1E-9),
                     None, None, float(data_x2[0]), 1]
                p_global = [data_y2.max() * 1.3, p_best_obs[1], p_best[1], p_best[2], p_best[3]]
                p_best, pcov, infodict, errmsg, success = scipy.optimize.leastsq(residualscalc_analytical1to1model,
                    p_global, args=(data_x1, data_y1, data_x2, data_y2, p_fix),
                                                                             full_output=1, epsfcn=0.0001, maxfev=20000)
                p_fix = [None, p_best[1], None, (float(exp['molarconcentration'][associationstepnumber]) * 1E-9),
                     p_best[3], p_best[4], float(data_x2[0]), None]
                p_global = [data_y2.max() * 1.3, p_best[2], .8]
                p_best, pcov, infodict, errmsg, success = scipy.optimize.leastsq(residualscalc_analytical1to1model,
                    p_global, args=(data_x1, data_y1, data_x2, data_y2, p_fix),
                                                                             full_output=1, epsfcn=0.0001, maxfev=20000)
                modeloff2 = analytical1to1modeloff(data_x2, assemblefitparam(p_best, p_fix).tolist())
                if (MCO*.2) > p_fix[4] > MCO and calcrsq(data_y2, modeloff2) > calcrsq(data_y2, modeloff):
                    exp['fit_param'] = assemblefitparam(p_best, p_fix).tolist()
                    exp['comments'].append("Decon Fit used")
                    exp['param_error'] = np.average([calcparamerror(p_best, p_fix, data_x1, data_y1, 0, pcov),
                        calcparamerror(p_best, p_fix, data_x2, data_y2, 1, pcov)], axis=0)
                    exp['fit_fxn'] = [0, 1]
                else:
                    p_fix = [None, None, None, None, None, None, None, None]
                    p_best, pcov, infodict, errmsg, success = scipy.optimize.leastsq(residualscalc_analytical1to1model,
                        p_partial, args=(data_x1, data_y1, data_x2, data_y2, p_fix),
                        full_output=1, epsfcn=0.0001, maxfev=20000)
                    exp['fit_param'] = assemblefitparam(p_best, p_fix).tolist()
                    if p_best[4] < MCO:
                        exp['flags'].append("MCO cutoff violated!")
                    exp['comments'].append("Partial Fit used")
                    exp['param_error'] = np.average([calcparamerror(p_partial, p_fix, data_x1, data_y1, 0, pcov),
                        calcparamerror(p_partial, p_fix, data_x2, data_y2, 1, pcov)], axis=0)
                    exp['fit_fxn'] = [0, 1]
            else:
                p_fix = [None, None, None, None, None, None, None, None]
                p_best, pcov, infodict, errmsg, success = scipy.optimize.leastsq(residualscalc_analytical1to1model,
                    p_partial, args=(data_x1, data_y1, data_x2, data_y2, p_fix),
                    full_output=1, epsfcn=0.0001, maxfev=20000)
                exp['fit_param'] = assemblefitparam(p_best, p_fix).tolist()
                exp['comments'].append("Partial Fit used")
                exp['param_error'] = np.average([calcparamerror(p_partial, p_fix, data_x1, data_y1, 0, pcov),
                    calcparamerror(p_partial, p_fix, data_x2, data_y2, 1, pcov)], axis=0)
                exp['fit_fxn'] = [0, 1]
            count += 1
        np.seterr(all=None)
        warnings.filterwarnings("default")
        if verbose:
            print '\n'
        paramcheck(experiment, associationstepnumber, dissociationstepnumber)
    else:
        p_global = [2.8E-2, 9.77E4]
        for exp in experiment:
            exp['param_name'] = ['Rmax', 'kd', 'ka', 'Cpro', 'm', 'c', 'X0']
            p_global.extend([5, -0.00001, 3])
        p_best, iternum = scipy.optimize.leastsq(residualscalc_analytical1to1model_global,
                                                 p_global, args=(experiment, associationstepnumber,
                                                dissociationstepnumber))
        count = 0
        for exp in experiment:
            exp['fit_fxn'].append([0, 1])
            p_best2 = [p_best[count+2], p_best[0], p_best[1], float(
                exp['molarconcentration'][associationstepnumber])*1E-9, p_best[count+3], p_best[count+4], float(
                exp['x_data'][dissociationstepnumber][0])]
            exp['fit_param'] = assemblefitparam(p_best2, fix=[None]).tolist()
            count += 3
    return experiment


def fit1to1model(experiment, associationstepnumber=1, dissociationstepnumber=2, globalfit=False, *args, **kwargs):
    offrange = None
    onrange = None
    partial = True
    verbose = False
    if 'off_xrange' in kwargs.keys():
        offrange = kwargs['off_xrange']
    if 'on_xrange' in kwargs.keys():
        onrange = kwargs['on_xrange']
    if 'partial_fit' in kwargs.keys():
        partial = kwargs['partial_fit']
    if 'verbose' in kwargs.keys():
        verbose = kwargs['verbose']

    if not globalfit:
        np.seterr(all="ignore")
        warnings.filterwarnings("ignore")
        count = 0
        totalexp = len(experiment)
        for exp in experiment:
            percentcounter = 100 * (count+1) / totalexp
            if not partial:
                if verbose:
                    print "\r1 to 1 model (std) fitting progress: " + str(percentcounter) +"%",
            else:
                if verbose:
                    print "\r1 to 1 model (parital) fitting progress: " + str(percentcounter)+"%",

            if exp['step_status'][-1].upper() is 'ERROR' or \
                    (exp['signal2noise'] is not None and float(exp['signal2noise']) <= 2.5):
                exp['param_name'] = ['Rmax', 'kd', 'ka', 'Cpro', 'm', 'c', 'X0', 'kds']
                exp['fit_param'] = [1, 1, 1, 1, 0, 0, exp['x_data'][dissociationstepnumber][0], 1]
                exp['param_error'] = [1E-12, 1E-12, 1E-12, 1E-12, 1E-12, 1E-12, 1E-12, 1E-12]
                exp['fit_fxn'] = [0, 1]
                if exp['step_status'][-1].upper() is 'ERROR':
                    exp['flags'].append('Sensor Error!')
                elif  (exp['signal2noise'] is not None and float(exp['signal2noise']) <= 2.5):
                    exp['flags'].append('Signal to noise below 2.5-fold')
                paramcheck([exp], associationstepnumber, dissociationstepnumber)
                count += 1
                continue

            if exp['molarconcentration'][associationstepnumber] is None or (
                    float(exp['molarconcentration'][associationstepnumber]) < 0):
                exp['molarconcentration'][associationstepnumber] = 100
            templist = exp['x_data']
            data_x1 = np.array(templist[associationstepnumber])
            data_x2 = np.array(templist[dissociationstepnumber])
            templist = exp['y_data']
            data_y1 = np.array(templist[associationstepnumber])
            data_y2 = np.array(templist[dissociationstepnumber])
            p_best = fitkdobs([exp], dissociationstepnumber)[0]
            # if segment fits only, data needs to be masked
            if onrange is not None and range(len(onrange)) == 2:
                m = np.ma.masked_outside(data_x1, onrange[0], onrange[1])
                data_x1 = data_x1[~m.mask]
                data_y1 = data_y1[~m.mask]
            if offrange is not None and range(len(offrange)) == 2:
                m = np.ma.masked_outside(data_x2, offrange[0], offrange[1])
                data_x2 = data_x2[~m.mask]
                data_y2 = data_y2[~m.mask]
            # initial parameters: Rmax, kd, ka, Cpro, m, c, X0
            exp['param_name'] = ['Rmax', 'kd', 'ka', 'Cpro', 'm', 'c', 'X0', 'kds']
            kdguess = p_best[1]
            # fit full
            p_fix = [None, p_best[1], None, (float(exp['molarconcentration'][associationstepnumber]) * 1E-9),
                     0, 0, float(data_x2[0]), 1]
            p_global = [data_y2.max() * 1.3, 2 / p_fix[3]]
            p_best, pcov, infodict, errmsg, success = scipy.optimize.leastsq(residualscalc_analytical1to1model,
                p_global, args=(data_x1, data_y1, data_x2, data_y2, p_fix),
                full_output=1, epsfcn=0.0001, maxfev=20000)
            kaguess = p_best[1]
            if partial:
                p_global = [data_y2.max() * 1.3, kdguess, kaguess]
                p_fix = [None, None, None, (float(exp['molarconcentration'][associationstepnumber]) * 1E-9),
                         0, data_y2[-1], float(data_x2[0]), 1]
                p_best, pcov, infodict, errmsg, success = scipy.optimize.leastsq(residualscalc_analytical1to1model,
                    p_global, args=(data_x1, data_y1, data_x2, data_y2, p_fix),
                    full_output=1, epsfcn=0.0001, maxfev=20000)
                p_fix = [None, None, None, (float(exp['molarconcentration'][associationstepnumber]) * 1E-9),
                         0, None, float(data_x2[0]), 1]
                p_global = [p_best[0], p_best[1], p_best[2], float(data_y2[-1])]
            else:
                if (np.abs(data_y2[0] - data_y2[-1]) / np.abs(data_y1[0] - data_y1[-1]) < 0.2):
                    kdguess = 1 / (2 * data_x2.max())
                    if (np.abs((data_y2[0]) - data_y2[len(data_y2) / 2])) / np.abs(data_y1[-1] - data_y1[0]) <= 0.05:
                        p_fix = [None, (-1 * np.log(0.95) / (data_x2[-1] - data_x2[0])), None,
                                 (float(exp['molarconcentration'][associationstepnumber]) * 1E-9), 0.0, 0.0,
                                 float(data_x2[0]), 1]
                        p_global = [data_y2.max() * 1.3, kaguess]
                        exp['comments'].append(
                            "5% rule invoked(kd fixed as: " + str(-1 * np.log(0.95) / (data_x2[-1] - data_x2[0])) + ")")
                    else:
                        p_fix = [None, None, None,
                                 (float(exp['molarconcentration'][associationstepnumber]) * 1E-9), 0.0, 0.0,
                                 float(data_x2[0]), 1]
                        p_global = [data_y2.max() * 1.3, kdguess, kaguess]
                    exp['flags'].append("Avid or aggregation like binding detected!")
                else:
                    p_fix = [None, None, None, (float(exp['molarconcentration'][associationstepnumber]) * 1E-9),
                             0.0, 0.0, float(data_x2[0]), 1]
                    p_global = [data_y2.max() * 1.3, kdguess, kaguess]
            #fit final
            p_best, pcov, infodict, errmsg, success = scipy.optimize.leastsq(residualscalc_analytical1to1model,
                p_global, args=(data_x1, data_y1, data_x2, data_y2, p_fix),
                full_output=1, epsfcn=0.0001, maxfev=20000)
            if p_best[1] < (-1*np.log(0.95)/(data_x2[-1]-data_x2[0])):
                exp['flags'].append(
                    "kd limit exceeded(slower than: " + str(-1 * np.log(0.95) / (data_x2[-1] - data_x2[0])) + ")")
            exp['fit_fxn'] = [0, 1]
            exp['fit_param'] = assemblefitparam(p_best, p_fix).tolist()
            exp['param_error'] = np.average([calcparamerror(p_best, p_fix, data_x1, data_y1, 0, pcov),
                                   calcparamerror(p_best, p_fix, data_x2, data_y2, 1, pcov)], axis=0)
            count += 1
        np.seterr(all=None)
        warnings.filterwarnings("default")
        if verbose:
            print '\n'
    else:
        p_global = [2.8E-2, 9.77E4]
        for exp in experiment:
            exp['param_name'] = ['Rmax', 'kd', 'ka', 'Cpro', 'm', 'c', 'X0']
            p_global.extend([5, -0.00001, 3])
        p_best, iternum = scipy.optimize.leastsq(residualscalc_analytical1to1model_global,
                                                 p_global, args=(experiment, associationstepnumber,
                                                dissociationstepnumber))
        count = 0
        for exp in experiment:
            exp['fit_fxn'].append([0, 1])
            p_best2 = [p_best[count+2], p_best[0], p_best[1], float(
                exp['molarconcentration'][associationstepnumber])*1E-9, p_best[count+3], p_best[count+4], float(
                exp['x_data'][dissociationstepnumber][0])]
            exp['fit_param'] = assemblefitparam(p_best2, fix=[None]).tolist()
            count += 3
    paramcheck(experiment, associationstepnumber, dissociationstepnumber)
    return experiment


def approxexp(x, y, PositiveAmp=False):
    if PositiveAmp == False:
        return (linearfit(x, np.log(y), False)[0])* -1
    else:
        temp = (linearfit(np.log(x)*-1, y, False))[0]*-1
        return temp


def fitkdobs(experiment, dissociationstepnumber, trylinearapproximation=True):
    p_return = []
    for exp in experiment:
        templist = exp['x_data']
        data_x = np.array(templist[dissociationstepnumber])
        temp_x = np.array(data_x[0:int(len(data_x) / 3)])
        templist = exp['y_data']
        data_y = np.array(templist[dissociationstepnumber])
        temp_y = np.array(data_y[0:int(len(data_y) / 3)])
        # fit single exp first with fixed amplitude
        exp['param_name'] = ['Amp', 'kd', 'X0']
        p_1 = approxexp(data_x, data_y)
        if p_1 <= 0:
            p_1 = p_1 * -1
        fix = [data_y.max()-data_y.min(), None, data_x[0]]
        p_best, iternum = scipy.optimize.leastsq(residualscalc,
                                                 p_1, args=(temp_x, temp_y, 3, fix))
        slope, intercept = linearfit(data_x, data_y, True)
        MCO = (-0.01146 * np.exp((data_x[-1]-data_x[0])/77.5588)) - (0.07192 * np.exp((data_x[-1]-data_x[0])/10.3))
        if slope < MCO:
            trylinearapproximation = False
            exp['comments'].append("MCO cutoff violated")
        # fit end range with linear term solve for rate, hold intercept
        if p_best[0] < 0:
            p_best[0] = float(p_best[0])*-1

        if not trylinearapproximation:
            p_1 = assemblefitparam(p_best, fix).tolist()
            fix = [None, None, None]
            p_best, iternum = scipy.optimize.leastsq(residualscalc,
                                                     p_1, args=(temp_x, temp_y, 3, fix))
            exp['fit_param'] = assemblefitparam(p_best, fix).tolist()
            p_return.append(exp['fit_param'])
        else:
            exp['param_name'] = ['Amp', 'kd', 'm', 'c', 'X0']
            p_2 = [p_best[0]]
            fix = [data_y.max()-data_y.min(), None, slope, intercept, data_x[0]]
            p_best, iternum = scipy.optimize.leastsq(residualscalc,
                                                     p_2, args=(data_x, data_y, 2, fix))
            kdval = p_best[0]
            # fit full range with linear term hold rate solve intercept
            exp['param_name'] = ['Amp', 'kd', 'm', 'c', 'X0']
            p_3 = [intercept]
            fix = [data_y.max()-data_y.min(), kdval, slope, None, data_x[0]]
            p_best, iternum = scipy.optimize.leastsq(residualscalc,
                                                     p_3, args=(data_x, data_y, 2, fix))
            intercept = p_best[0]
            # fre all
            exp['param_name'] = ['Amp', 'kd', 'm', 'c', 'X0']
            p_3 = [data_y.max()-data_y.min(), kdval, slope, intercept]
            fix1 = [None, None, None, None, data_x[0]]
            p_best1, iternum1 = scipy.optimize.leastsq(residualscalc,
                                                     p_3, args=(data_x, data_y, 2, fix1))

            # fre all no linear term
            exp['param_name'] = ['Amp', 'kd', 'm', 'c', 'X0']
            p_3 = [data_y.max()-data_y.min(), kdval, data_y[-1]]
            fix2 = [None, None, 0, None, data_x[0]]
            p_best2, iternum2 = scipy.optimize.leastsq(residualscalc,
                                                     p_3, args=(data_x, data_y, 2, fix2))


            explin = (kdobs(data_x, assemblefitparam(p_best1, fix1)) - data_y)
            exponly = (kdobs(data_x, assemblefitparam(p_best2, fix2)) - data_y)

            if (sum(explin) <= sum(exponly)*0.9) and (1 - data_y[-1]/data_y[0] >= 0.2) and p_best1[0] > 0 and \
                            p_best2[0] > 0 and exp['signal2noise'] is not None and exp['signal2noise'] >= 2.5:
                exp['fit_fxn'] = [2]
                exp['fit_param'] = assemblefitparam(p_best1, fix1).tolist()
                p_return.append(p_best1)
            else:
                exp['fit_fxn'] = [3]
                templist = assemblefitparam(p_best2, fix2).tolist()
                exp['fit_param'] = [templist[0], templist[1], templist[4]]
                exp['param_name'] = ['Amp', 'kd', 'X0']
                p_return.append([templist[0], templist[1], templist[4]])
    return p_return


def fitkaobs(experiment, associationstepnumber):
    p_return = []
    for exp in experiment:
        templist = exp['x_data']
        data_x = np.array(templist[associationstepnumber])
        temp_x = np.array(data_x[0:int(len(data_x) / 3)])
        templist = exp['y_data']
        data_y = np.array(templist[associationstepnumber])
        temp_y = np.array(data_y[0:int(len(data_y) / 3)])
        fiftyp = templist[associationstepnumber][int(len(templist[associationstepnumber])/2)]
        fmin = np.array(exp['y_data'][associationstepnumber - 1]).min()
        fmax = data_y.max()
        fiftyp = ((fiftyp - fmin) / (fmax - fmin))
        # fit single exp first with fixed amplitude
        if fiftyp < 0.5:
            exp['param_name'] = ['Amp', 'kd', 'X0']
            p_1 = approxexp(data_x, data_y, True)
            fix = [fmax-fmin, None, data_x[0]]
            p_best, iternum = scipy.optimize.leastsq(residualscalc,
                                                     p_1, args=(temp_x, temp_y, 4, fix))
        else:
            exp['param_name'] = ['Amp', 'kd', 'X0']
            p_1 = approxexp(data_x, data_y, True)
            fix = [fmax-fmin, None, data_x[0]]
            p_best, iternum = scipy.optimize.leastsq(residualscalc,
                                                     p_1, args=(data_x, data_y, 4, fix))
        if p_best[0] < 0:
            p_best[0] = float(p_best[0])*-1

        p_1 = [fmax-fmin, p_best[0]]
        fix = [None, None, data_x[0]]
        p_best, iternum = scipy.optimize.leastsq(residualscalc,
                                                 p_1, args=(data_x, data_y, 4, fix))

        p_best = assemblefitparam(p_best, fix).tolist()
        exp['fit_param'] = p_best
        p_return.append(p_best)
    return p_return


def singleexpPos(x,p):
    # parameters: Amp, k, xo
    a, k, xo = p
    return (a * (1-np.exp(-1 * ((x - xo) * k))))


def singleexp(x,p):
    # parameters: Amp, k, xo
    a, k, xo = p
    return (a * np.exp(-1 * ((x - xo) * k)))


def kdobs(x, p):
    #parameters: Amp, kd, m, c
    a, kd, m, c, xo = p
    return (a*np.exp(-1*((x-xo)*kd)) + m*(x-xo) + c)


def analytical1to1modelon(x, p):
    # parameters: Rmax, kd, ka, Cpro, m, c, X0, kds
    rm, kd, ka, cp, m, c, xo, kds = p
    if cp==0:
        cp=0.001
    return (rm / (1 + (kd / (ka * cp)))) * (1-(np.exp((-1 * x * ((ka * cp) + kd)))))


def analytical1to1modeloff(x, p):
    #parameters: Rmax, kd, ka, Cpro, m, c, X0, kds
    rm, kd, ka, cp, m, c, xo, kds = p
    if cp == 0:
        cp = 0.001
    return kds*(((rm/(1+(kd/(ka*cp))))*(1-np.exp(-1*xo*(ka*cp+kd)))-c))*(np.exp(-1*kd*(x-xo))) + m*(x-xo) + c


def residualscalc_analytical1to1model(p, x1, y1, x2, y2, fix):
    # parameters: Rmax, kd, ka, Cpro, m, c, X0, kds
    #setting parameters for each dataset
    err1 = residualscalc(p, x1, y1, 0, fix) #calculate the goodness of fit for parameterset 1, y1
    err2 = residualscalc(p, x2, y2, 1, fix) #calculate the goodness of fit for parameterset 2, y2
    return np.concatenate((err1, err2)) #combines err1 and err2 to reconstruct p.


def residualscalc_analytical1to1model_global(p, experiment, associationstepnumber,
                                                dissociationstepnumber):
    # parameters: Rmax, kd, ka, Cpro, m, c, X0
    test = True
    count = 0
    for exp in experiment:
        p1 = p[count + 2], p[0], p[1], float(exp['molarconcentration'][associationstepnumber])*1E-9, p[count + 3], p[count + 4] #setting parameters for each dataset
        p2 = p[count + 2], p[0], p[1], float(exp['molarconcentration'][associationstepnumber])*1E-9, \
             p[count + 3], p[count + 4], float(exp['x_data'][dissociationstepnumber][0])
        err1 = residualscalc(p1, np.array(exp['x_data'][associationstepnumber]),
                             np.array(exp['y_data'][associationstepnumber]), 0) #calculate the goodness of fit for parameterset 1, y1
        err2 = residualscalc(p2, np.array(exp['x_data'][dissociationstepnumber]),
                             np.array(exp['y_data'][dissociationstepnumber]), 1) #calculate the goodness of fit for parameterset 2, y2
        if test:
            summederror = np.concatenate((err1, err2))
        else:
            summederror = np.concatenate(summederror, np.concatenate((err1, err2)))
        count += 3
    return summederror

def residualscalc(p, x, y, fxn, fix): # calculate residuals to determine if the parameters are improving the fit
    p = assemblefitparam(p, fix)
    if fxn == 0:
        # parameters: Rmax, kd, ka, Cpro, m, c, X0, kds
        penalty = 1
        if p[0] < y.max()*.7 or p[0] > y.max()*10:
            penalty += 500
        if p[1] == 0 or p[2] == 0:
            penalty += 5000000
        if p[1] <= 0 or p[1] > 1:
            penalty += 5000000
        if p[2] <= 1000:
            penalty += 5000000
        if p[2] <= 0:
            penalty += 5000000
        if p[1] / p[2] < 1E-11:
            penalty += 50000
        return np.power((analytical1to1modelon(x, p) - y), 2) * penalty
    elif fxn == 1:
        # parameters: Rmax, kd, ka, Cpro, m, c, X0, kds
        penalty = 1
        if p[1] == 0 or p[2] == 0:
            penalty += 5000000
        if p[1] <= 0 or p[1] > 1:
            penalty += 5000000
        if p[2] <= 1000:
            penalty += 5000000
        if p[2] <= 0:
            penalty += 5000000
        if p[5] > y.max():
            penalty += 5000000
        if p[1] / p[2] < 1E-11:
            penalty += 50000
        return np.power((analytical1to1modeloff(x, p) - y), 2) * penalty
    elif fxn == 2:
        penalty = 1
        if p[1] <= 0 or p[1] > 1:
            penalty = 50000
        if p[2] > 0:
            penalty += 50000
        return np.power((kdobs(x, p) - y), 2) * penalty
    elif fxn == 3:
        penalty = 1
        if p[1] <= 0 or p[1] > 1:
            penalty = 50000
        return np.power((singleexp(x, p) - y), 2) * penalty
    elif fxn == 4:
        return np.power((singleexpPos(x, p) - y), 2)
    else:
        return False


def calcparamerror(pval, pfix, x, y, fxn, pcov):
    p = assemblefitparam(pval, pfix)
    s_sq = 0
    if (len(y) > len(p)) and pcov is not None:
        if fxn == 0:
            s_sq = np.power((analytical1to1modelon(x, p) - y), 2).sum()/(len(y)-len(p))
        elif fxn == 1:
            s_sq = np.power((analytical1to1modeloff(x, p) - y), 2).sum()/(len(y)-len(p))
        elif fxn == 2:
            s_sq = np.power((kdobs(x, p) - y), 2).sum()/(len(y)-len(p))
        elif fxn == 3:
            s_sq = np.power((singleexp(x, p) - y), 2).sum()/(len(y)-len(p))
        elif fxn == 4:
            s_sq = np.power((singleexpPos(x, p) - y), 2).sum() / (len(y) - len(p))
        else:
            return False
        pcov = pcov * s_sq
    else:
        pcov = np.inf
    perror = []
    for i in range(len(pval)):
        try:
            perror.append(np.absolute(pcov[i][i]) ** 0.5)
        except:
            perror.append(0.00)
    return assembleparamerror(perror, pfix)


def plotresults(exp, stackdepth=0):
    plot = pl.figure()
    ax = plot.add_subplot(1, 1, 1)
    pl.xlim(exp['x_data'][1][0], exp['x_data'][2][-1])
    pl.scatter(np.array(exp['x_data'][1]), np.array(exp['y_data'][1]), color="black", s=0.5)
    pl.scatter(np.array(exp['x_data'][2]), np.array(exp['y_data'][2]), color="black", s=0.5)
    if stackdepth >= 1:
        pl.plot(np.array(exp['x_data'][1]), analytical1to1modelon(np.array(exp['x_data'][1]),
                                                                     exp['fit_param']), color="blue")
        pl.plot(np.array(exp['x_data'][2]), analytical1to1modeloff(np.array(exp['x_data'][2]),
                                                                     exp['fit_param']), color="red")
    pl.show()

    return


def saveplot(exp, outfile, stackdepth=0, returnBytesIO=False, associationstep=1, dissociationstep=2, *args, **kwargs):
    offrange = None
    onrange = None
    if 'off_xrange' in kwargs.keys()and kwargs['off_xrange'] is not None:
        offrange = kwargs['off_xrange']
    if 'on_xrange' in kwargs.keys() and kwargs['on_xrange'] is not None:
        onrange = kwargs['on_xrange']
    filename = ''
    if exp == []:
        return False
    if not returnBytesIO:
        fileext = os.path.splitext(outfile)[1]
        outpath, outfile = os.path.split(os.path.abspath(outfile))
        if fileext == '':
            fileext = '.png'
            outpath = os.path.join(outpath, outfile)
            outfile = 'plot' + fileext
        if not os.path.exists(outpath):
            try:
                os.mkdir(outpath)
            except OSError:
                return False
        outfile = outfile[:len(outfile) - len(fileext)]
        filename = os.path.join(outpath, outfile + fileext)
    plot = pl.figure()
    ax = plot.add_subplot(1, 1, 1)
    #pl.ylabel('Response (nm)', size=20)
    pl.ylabel('Fractional Response', size=20)
    pl.xlabel('Time (sec)', size=20)
    pl.grid()
    pl.gcf().subplots_adjust(bottom=0.15)
    pl.gcf().subplots_adjust(left=0.15)
    pl.xticks(np.arange(min(exp['x_data'][associationstep]), max(exp['x_data'][dissociationstep]), 100))
    pl.tick_params(axis='both', which='major', labelsize=20)
    pl.xlim(exp['x_data'][associationstep][0], exp['x_data'][dissociationstep][-1])
    pl.plot(np.array(exp['x_data'][associationstep]), np.array(exp['y_data'][associationstep]),
            linewidth=2, color="blue")
    pl.plot(np.array(exp['x_data'][dissociationstep]), np.array(exp['y_data'][dissociationstep]),
            linewidth=2, color="blue")
    if stackdepth >= 1:
        data_x1 = np.ma.array(exp['x_data'][associationstep])
        data_x2 = np.ma.array(exp['x_data'][dissociationstep])
        data_y1 = analytical1to1modelon(np.ma.array(exp['x_data'][associationstep]), exp['fit_param'])
        data_y2 = analytical1to1modeloff(np.ma.array(exp['x_data'][dissociationstep]), exp['fit_param'])
        if onrange is not None and len(offrange) == 2:
            m = np.ma.masked_outside(data_x1, onrange[0], onrange[1])
            data_x1 = data_x1[~m.mask]
            data_y1 = data_y1[~m.mask]
        pl.plot(data_x1, data_y1, linewidth=2, color="red")
        if offrange is not None and len(offrange) == 2:
            m = np.ma.masked_outside(data_x2, offrange[0], offrange[1])
            data_x2 = data_x2[~m.mask]
            data_y2 = data_y2[~m.mask]
        pl.plot(data_x2, data_y2, linewidth=2, color="red")
    if not returnBytesIO:
        pl.savefig(filename)
        pl.close()
        return
    else:
        BIOstream = BIO()
        pl.savefig(BIOstream, format='png')
        BIOstream.seek(0)
        byte_png = BIOstream.getvalue()
        pl.close()
        return byte_png

    return


def assemblefitparam(p, fix):
    if p is None:
        return np.array(fix)
    if fix is None:
        return np.array(p)
    param = p
    temp = []
    counter = 0
    for f in fix:
        if f is None:
            temp.append(param[counter])
            counter += 1
        else:
            temp.append(f)
    p = np.array(temp)
    return p


def assembleparamerror(p, fix):
    if p is None:
        return p
    temp = []
    counter = 0
    for f in fix:
        if f is None:
            temp.append(p[counter])
            counter += 1
        else:
            temp.append(0.00)
    p = np.array(temp)
    return p


def paramcheck(experiment, associationstepnumber=0, dissociationstepnumber=1):
    for exp in experiment:
        for x in range(len(exp['fit_param'])):
            if exp['fit_param'][x] is None or math.isnan(float(exp['fit_param'][x])):
                exp['fit_param'][x] = 1
                exp['param_error'][x] = 1E-12
            if x<=3 and exp['fit_param'][x] < 0 or \
                            exp['y_data'][associationstepnumber][-1] <= exp['y_data'][associationstepnumber][0]:
                exp['fit_param'] = [1, 1, 1, 1, 0, 0, exp['x_data'][dissociationstepnumber][0], 1]
                exp['param_error'] = [1E-12, 1E-12, 1E-12, 1E-12, 1E-12, 1E-12, 1E-12, 1E-12]


def calcrsq(rawy, modely):
    rsq = np.sum(np.power(modely - rawy, 2))
    avgerror = np.sum(np.power(rawy - np.average(rawy), 2))
    return 1 - (rsq / avgerror)


def calcxsq(rawy, modely):
    SD = np.average(np.abs(rawy - modely))
    return np.sum(((rawy - modely)**2)/SD)

