import numpy as np
import scipy.optimize
import matplotlib.pyplot as pl
import os


def __init__(self):
    """ No parameters required """
    pass


def linearfit(baselinex, baseliney, onlylast30percent=False):
    yintercept = 0
    sumxy = 0
    sumx = 0
    sumy = 0
    sumxsq = 0
    baselinex.tolist()
    baseliney.tolist()
    if onlylast30percent:
        i = int(len(baselinex) / 3) * 2
    else:
        i=0
    count = len(baselinex) - i
    while i < len(baselinex):
        sumxy += baselinex[i] * baseliney[i]
        sumx += baselinex[i]
        sumy += baseliney[i]
        sumxsq += baselinex[i] * baselinex[i]
        i += 1
    slope = ((count * sumxy) - (sumx * sumy)) / ((count * sumxsq) - (sumx * sumx))
    yintercept = baseliney[1] - (slope * baselinex[1])
    return [slope, yintercept]

def fit1to1model(experiment, associationstepnumber=1, dissociationstepnumber=2, globalfit=False):
    if not globalfit:
        count = 0
        for exp in experiment:
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
            print str(count) + ' kobs fit'
            #p_fix = concentraion, xo, kd, m, c
            if exp['fit_fxn'][0] == 3:
                p_best[2] = 0
            if (1 - data_y2[-1]/data_y2[0] >= 0.85) or (1 - data_y2[-1]/data_y2[0] < 0.2):
                p_best[3] = 0
            p_fix = [None, p_best[1], None, (float(exp['molarconcentration'][associationstepnumber])*1E-9),
                     None, None, float(data_x2[0]), None]
            p_global = [data_y2.max() * 1.3, 9.77E4, p_best[2], p_best[3], 0.5]
            if 1 - data_y2[-1]/data_y2[0] <= 0.2:
                p_fix[5] = 0.0
                p_fix[4] = 0.0
                p_fix[1] = (1/(2*data_x2.max()))
                exp['flags'].append("Avid or aggregation like binding detected")
                exp['flags'].append("kd limit exceeded(slower than: "+str(1/(2*data_x2.max()))+")")
            # initial parameters: Rmax, kd, ka, Cpro, m, c, X0
            exp['param_name'] =['Rmax', 'kd', 'ka', 'Cpro', 'm', 'c', 'X0', 'koff scalar']
            p_best, iternum = scipy.optimize.leastsq(residualscalc_analytical1to1model,
                                                p_global, args=(data_x1, data_y1, data_x2, data_y2, p_fix))
            print str(count) + ' global fit'
            exp['fit_fxn'] = [0, 1]
            exp['fit_param'] = assemblefitparam(p_best, p_fix).tolist()
            count += 1
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
            p_best1 = [p_best[count+2], p_best[0], p_best[1], float(
                exp['molarconcentration'][associationstepnumber])*1E-9, p_best[count+3], p_best[count+4]]
            p_best2 = [p_best[count+2], p_best[0], p_best[1], float(
                exp['molarconcentration'][associationstepnumber])*1E-9, p_best[count+3], p_best[count+4], float(
                exp['x_data'][dissociationstepnumber][0])]
            exp['fit_param'] = assemblefitparam(p_best2, fix=[None]).tolist()
            count += 3
    return experiment


def fitkdobs(experiment, dissociationstepnumber):
    p_return = []
    for exp in experiment:
        templist = exp['x_data']
        data_x = np.array(templist[dissociationstepnumber])
        temp_x = np.array(data_x[0:int(len(data_x) / 3)])
        templist = exp['y_data']
        data_y = np.array(templist[dissociationstepnumber])
        temp_y = np.array(data_y[0:int(len(data_y) / 3)])
        fiftyp = templist[dissociationstepnumber][int(len(templist[dissociationstepnumber])/2)]
        fmin = np.array(exp['y_data'][dissociationstepnumber - 1]).min()
        fmax = data_y.max()
        fiftyp = ((fiftyp - fmin) / (fmax - fmin))
        # fit single exp first with fixed amplitude
        if fiftyp < 0.5:
            exp['param_name'] = ['Amp', 'kd', 'X0']
            p_1 = [.0001]
            fix = [data_y.max()-data_y.min(), None, data_x[0]]
            p_best, iternum = scipy.optimize.leastsq(residualscalc,
                                                     p_1, args=(temp_x, temp_y, 3, fix))
            slope, intercept = linearfit(data_x, data_y, True)
        else:
            exp['param_name'] = ['Amp', 'kd', 'X0']
            p_1 = [.000001]
            fix = [data_y.max()-data_y.min(), None, data_x[0]]
            p_best, iternum = scipy.optimize.leastsq(residualscalc,
                                                     p_1, args=(data_x, data_y, 3, fix))
            slope, intercept = linearfit(data_x, data_y, True)
        # fit full range with linear term solve for rate, hold intercept
        if p_best[0] < 0:
            p_best[0] = float(p_best[0])*-1

        exp['param_name'] = ['Amp', 'kd', 'm', 'c', 'X0']
        p_2 = [p_best[0]]
        fix = [data_y.max()-data_y.min(), None, slope, intercept, data_x[0]]
        p_best, iternum = scipy.optimize.leastsq(residualscalc,
                                                 p_2, args=(data_x, data_y, 2, fix))
        # fit full range with linear term hold rate solve intercept
        kdval = p_best[0]
        exp['param_name'] = ['Amp', 'kd', 'm', 'c', 'X0']
        p_3 = [intercept]
        fix = [data_y.max()-data_y.min(), kdval, slope, None, data_x[0]]
        p_best, iternum = scipy.optimize.leastsq(residualscalc,
                                                 p_3, args=(data_x, data_y, 2, fix))
        # fre all
        kdval = p_best[0]
        exp['param_name'] = ['Amp', 'kd', 'm', 'c', 'X0']
        p_3 = [data_y.max()-data_y.min(), kdval, slope, p_best[0]]
        fix1 = [None, None, None, None, data_x[0]]
        p_best1, iternum1 = scipy.optimize.leastsq(residualscalc,
                                                 p_3, args=(data_x, data_y, 2, fix1))

        # fre all no linear term
        kdval = p_best[0]
        exp['param_name'] = ['Amp', 'kd', 'm', 'c', 'X0']
        p_3 = [data_y.max()-data_y.min(), kdval, 0, p_best[0]]
        fix2 = [None, None, 0, None, data_x[0]]
        p_best2, iternum2 = scipy.optimize.leastsq(residualscalc,
                                                 p_3, args=(data_x, data_y, 2, fix2))


        explin = (kdobs(data_x, assemblefitparam(p_best1, fix1)) - data_y)
        exponly = (kdobs(data_x, assemblefitparam(p_best2, fix2)) - data_y)

        if (sum(explin) <= sum(exponly)*0.9) and (1 - data_y[-1]/data_y[0] >= 0.2) and p_best1[0] > 0 and \
                        p_best2[0] > 0 and exp['signal2noise'] is not None and exp['signal2noise'] >= 5:
            exp['fit_fxn'] = [2]
            exp['fit_param'] = assemblefitparam(p_best1, fix1).tolist()
            exp['comments'] = ("Linear approximation used")
            p_return.append(p_best1)
        else:
            exp['fit_fxn'] = [3]
            exp['fit_param'] = assemblefitparam(p_best2[:3], fix2).tolist()
            exp['param_name'] = ['Amp', 'kd', 'X0']
            p_return.append(p_best2)
    return p_return

def singleexp(x,p):
    # parameters: Amp, kd, xo
    a, kd, xo = p
    return (a * np.exp(-1 * ((x - xo) * kd)))

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
    if cp==0:
        cp=0.001
    return kds*((rm*(1/(1+(kd/(ka*cp))))*(1-np.exp(-1*xo*(ka*cp+kd)))*(np.exp(-1*kd*(x-xo))))) + m*(x-xo) + c


def residualscalc_analytical1to1model(p, x1, y1, x2, y2, fix):
    # parameters: Rmax, kd, ka, Cpro, m, c, X0, , kds
    #setting parameters for each dataset
    err1 = residualscalc(p, x1, y1, 0, fix) #calculate the goodness of fit for parameterset 1, y1
    err2 = residualscalc(p, x2, y2, 1, fix) #calculate the goodness of fit for parameterset 2, y2
    return np.concatenate((err1, err2)) #combines err1 and err2 to reconstruct p.  In this way the only one [antigen] value is returned for the minimization, effectively linking the parameter.


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
        return np.power((analytical1to1modelon(x, p) - y), 2)
    elif fxn == 1:
        penalty = 1
        if p[0] < 0 or p[1] < 0 or p[2] < 0 or p[7] < 0:
            penalty = 100000
        return (analytical1to1modeloff(x, p) - y)*penalty
    elif fxn == 2:
        return (kdobs(x, p) - y)
    elif fxn == 3:
        return (singleexp(x, p) - y)
    else:
        return False


def plotresults(exp, stackdepth=0):
    plot = pl.figure()
    ax = plot.add_subplot(1, 1, 1)
    pl.xlim(exp['x_data'][1][0], exp['x_data'][2][-1])
    pl.scatter(np.array(exp['x_data'][1]), np.array(exp['y_data'][1]), color="black", s=0.5)
    pl.scatter(np.array(exp['x_data'][2]), np.array(exp['y_data'][2]), color="black", s=0.5)
    if stackdepth == 1:
        pl.plot(np.array(exp['x_data'][1]), analytical1to1modelon(np.array(exp['x_data'][1]),
                                                                     exp['fit_param']), color="blue")
        pl.plot(np.array(exp['x_data'][2]), analytical1to1modeloff(np.array(exp['x_data'][2]),
                                                                     exp['fit_param']), color="red")
    pl.show()

    return


def saveplot(exp, outfile, stackdepth=0):
    if exp == []:
        return False
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
    pl.ylabel('Response (nm)', size=20)
    pl.xlabel('Time (sec)', size=20)
    pl.grid()
    pl.gcf().subplots_adjust(bottom=0.15)
    pl.gcf().subplots_adjust(left=0.15)
    pl.xticks(np.arange(min(exp['x_data'][1]), max(exp['x_data'][2]), 100))
    pl.tick_params(axis='both', which='major', labelsize=20)
    pl.xlim(exp['x_data'][1][0], exp['x_data'][2][-1])
    pl.plot(np.array(exp['x_data'][1]), np.array(exp['y_data'][1]), linewidth=2, color="blue")
    pl.plot(np.array(exp['x_data'][2]), np.array(exp['y_data'][2]), linewidth=2, color="blue")
    if stackdepth == 1:
        pl.plot(np.array(exp['x_data'][1]), analytical1to1modelon(np.array(exp['x_data'][1]),
                                                                  exp['fit_param']), linewidth=2, color="#F8BB41")
        pl.plot(np.array(exp['x_data'][2]), analytical1to1modeloff(np.array(exp['x_data'][2]),
                                                                   exp['fit_param']), linewidth=2, color="red")
    pl.savefig(filename)
    pl.close()

    return


def assemblefitparam(p, fix):
    if p is None:
        return np.array(fix)
    param = p.tolist()
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
