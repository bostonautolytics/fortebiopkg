import datafitting
import numpy as np

def __init__(self):
    """ No parameters required """
    pass


def removeblanksensors(experiment, sensor_name='all', sensor_info=False, concentration=False,
                       molarconcentration=False, molecularweight=False, sampleid=False):
    for i in reversed(range(0, len(experiment))):
        if sensor_name != 'all' and experiment[i].get('sensor') != sensor_name:
            continue
        else:
            if sensor_info and experiment[i].get('sensor_info') is None:
                experiment.pop(i)
                continue
            elif sensor_info and (experiment[i].get('sensor_info') == "NONE" or experiment[i].get(
                                  'sensor_info') == '' or experiment[i].get('sensor_info') == "NULL"):
                experiment.pop(i)
                continue
            if not concentration and not molarconcentration and not molecularweight and not sampleid:
                continue
            else:
                keepflag = True #keep by default
                concentration_list = experiment[i].get('concentration')
                molarconcentration_list = experiment[i].get('molarconcentration')
                molecularweight_list = experiment[i].get('molecularweight')
                sampleid_list = experiment[i].get('sampleid')
                welltype_list = experiment[i].get('welltype')
                for j in range(0, len(experiment[i].get('step_name'))):
                    if concentration and concentration_list[j] is None and str(welltype_list[j]).upper() != 'BUFFER':
                        keepflag = False
                    elif concentration and concentration_list[j] is not None and str(
                                          welltype_list[j]).upper() != 'BUFFER' and (str(concentration_list[j]).upper()
                                          == 'NULL' and
                                          str(concentration_list[j]).upper() == '' and
                                          concentration_list[j] <= 0):
                        keepflag = False
                    if molarconcentration and molarconcentration_list[j] is None and str(
                            welltype_list[j]).upper() != 'BUFFER':
                        keepflag = False
                    elif molarconcentration and molarconcentration_list[j] is not None and str(
                                                            welltype_list[j]).upper() != 'BUFFER' and (str(
                                                            molarconcentration_list[j]).upper() != 'NULL' and
                                                            str(molarconcentration_list[j]).upper() != '' and
                                                            molarconcentration_list[j] <= 0):
                        keepflag = False
                    if molecularweight and molecularweight_list[j] is None and str(
                            welltype_list[j]).upper() != 'BUFFER':
                        keepflag = False
                    elif molecularweight and molecularweight_list[j] is not None and str(
                                                            welltype_list[j]).upper() != 'BUFFER' and (str(
                                                            molecularweight_list[j]).upper() != 'NULL' and
                                                            str(molecularweight_list[j]).upper() != '' and
                                                            molecularweight_list[j] <= 0):
                        keepflag = False
                    if sampleid and sampleid_list[j] is None and str(welltype_list[j]).upper() != 'BUFFER':
                        keepflag = False
                    elif sampleid and sampleid_list[j] is not None and str(
                                                            welltype_list[j]).upper() != 'BUFFER' and (str
                                                            (sampleid_list[j]).upper() == 'NULL' and
                                                            str(sampleid_list[j]).upper() == ''):
                        keepflag = False
                    if not keepflag:
                        experiment.pop(i)
                        break
    return experiment


def blankstepsubtraction(experiment, stepnumberofblank=0, stepnumberofassociation=1, step2correct='all'):
    for exp in experiment:
        xdata = exp['x_data']
        ydata = exp['y_data']
        baselinex = xdata[stepnumberofblank]
        baseliney = ydata[stepnumberofblank]
        #### get linear fit of last 30% of baseline ####
        slope, yintercept = datafitting.linearfit(baselinex, baseliney, True)
        #### Correct data steps with linear model ####
        for j in range(0, len(xdata)):
            if j == step2correct or step2correct == 'all':###correct all data
                for k in range(0, len(xdata[j])):
                    ydata[j][k] -= (slope * xdata[j][k]) + yintercept
        if step2correct == stepnumberofassociation or step2correct == 'all': ##if association phase was corrected, check
            #### Check association slope (last 30%) is not negative #####
            slopeA, yinterceptA = datafitting.linearfit(xdata[stepnumberofassociation],
                                  ydata[stepnumberofassociation], True)
            if slopeA < 0: ###if negative slope correct to 0
                slopeA = slopeA * -1
                yinterceptA = ydata[j][1] - (slopeA * xdata[j][1])
                for j in range(0, len(xdata)):
                    if j == step2correct or step2correct == 'all':  ###correct all data
                        for k in range(0, len(xdata[j])):
                            ydata[j][k] += (slopeA * xdata[j][k]) + yinterceptA
        #### Check to make sure re-correction was sufficient ####
        slopeA, yinterceptA = datafitting.linearfit(xdata[stepnumberofassociation],
                                                    ydata[stepnumberofassociation], True)
        slopeB, yinterceptB = datafitting.linearfit(xdata[stepnumberofassociation+1],
                                                    ydata[stepnumberofassociation+1], True)
        if slopeA > 0 and slopeB > 0 :
            exp['comments'].append('Sensogram is not self-background subtracted, reference subtraction required')
        else:
            exp['x_data'] = xdata
            exp['y_data'] = ydata
    return experiment


def blanksubtractionsequential(experiment, firstbaselinestep=0):
    for exp in experiment:
        xdata = exp['x_data']
        ydata = exp['y_data']
        step_name = exp['step_name']
        baselinex = []
        baseliney = []
        #### Correct data steps with linear model ####
        j = 0
        while j < len(xdata):
            if j >= firstbaselinestep:
                baselinex = xdata[j]
                baseliney = ydata[j]
                #### get linear fit of last 30% of baseline ####
                slope, yintercept = datafitting.linearfit(baselinex, baseliney, True)
            for k in range(0, len(xdata[j])):
                ydata[j][k] -= (slope * xdata[j][k]) + yintercept
            for k in range(0, len(xdata[j+1])):
                ydata[j+1][k] -= (slope * xdata[j+1][k]) + yintercept
            for k in range(0, len(xdata[j+2])):
                ydata[j+2][k] -= (slope * xdata[j+2][k]) + yintercept
            #### Check association slope (last 30%) is not negative #####
            slopeA, yinterceptA = datafitting.linearfit(xdata[j+1],
                                  ydata[j+1], True)
            if slopeA < 0: ###if negative slope correct to 0
                slopeA = slopeA * -1
                yinterceptA = ydata[j+1][1] - (slopeA * xdata[j+1][1])
                for k in range(0, len(xdata[j])):
                    ydata[j][k] += (slopeA * xdata[j][k]) + yinterceptA
                for k in range(0, len(xdata[j+1])):
                    ydata[j+2][k] += (slopeA * xdata[j+2][k]) + yinterceptA
                for k in range(0, len(xdata[j+2])):
                    ydata[j+2][k] += (slopeA * xdata[j+2][k]) + yinterceptA
            exp['x_data'] = xdata
            exp['y_data'] = ydata
            j += 3
    return experiment


def normalizetoassociation (experiment, associationstepnumber=1):
    for exp in experiment:
        maxval = exp['y_data'][associationstepnumber][-1]
        minval = exp['y_data'][associationstepnumber][0]
        for step in exp['y_data']:
            for i in range(len(step)):
                step[i] = (step[i] - minval)/(maxval - minval)
    return experiment


def subtractreferncesensor (referenceexp, subtractfromexp):
    for step in range(len(subtractfromexp['y_data'])):
        for i in range(len(subtractfromexp['y_data'][step])):
            subtractfromexp['y_data'][step][i] -= referenceexp['y_data'][step][i]
    return subtractfromexp


def interstepcorrection(experiment):
    for exp in experiment:
        xdata = exp['x_data']
        ydata = exp['y_data']
        for j in range(0, len(xdata)-1):
            if ydata[j][-1] > ydata[j + 1][0]:
                for k in range(0, len(xdata[j+1])):
                    ydata[j+1][k] += (ydata[j][-1] - ydata[j + 1][0])
            else:
                for k in range(0, len(xdata[j+1])):
                    ydata[j+1][k] -= (ydata[j+1][0] - ydata[j][-1])
        exp['y_data'] = ydata
    return experiment


def setstepatorigin(experiment, stepnumber):
    for exp in experiment:
        xdata = exp['x_data']
        ydata = exp['y_data']
        originx = xdata[stepnumber][0]
        originy = ydata[stepnumber][0]
        for j in range(0, len(xdata)):
            for k in range(0, len(xdata[j])):
                xdata[j][k] -= originx
                ydata[j][k] -= originy
    return experiment


def setassociationtoorigin(experiment, baselinestepnumber, associationstepnumber, dissociationstepnumber):
    for exp in experiment:
        originx = float(exp['x_data'][associationstepnumber][0])
        originy = float(exp['y_data'][associationstepnumber][0])
        exp['x_data'][baselinestepnumber] = np.array(exp['x_data'][baselinestepnumber]) - originx
        exp['y_data'][baselinestepnumber] = np.array(exp['y_data'][baselinestepnumber]) - originy
        exp['x_data'][associationstepnumber] = np.array(exp['x_data'][associationstepnumber]) - originx
        exp['y_data'][associationstepnumber] = np.array(exp['y_data'][associationstepnumber]) - originy
        exp['x_data'][dissociationstepnumber] = np.array(exp['x_data'][dissociationstepnumber]) - originx
        exp['y_data'][dissociationstepnumber] = np.array(exp['y_data'][dissociationstepnumber]) - originy
    return experiment


def getsignaltonoise(experiment, baselinestepnumber=0, associationstepnumber=1):
    for exp in experiment:
        ydata = exp['y_data']
        backgroundnoise = ydata[baselinestepnumber]
        backgroundnoise = backgroundnoise[(len(backgroundnoise)/5):]
        fullresponse = ydata[associationstepnumber]
        backgroundamplitude = max(backgroundnoise)-min(backgroundnoise)
        signalamplitude = max(fullresponse)-min(fullresponse)
        exp['signal2noise'] = signalamplitude/backgroundamplitude
    return experiment


def scalesensortomatchamplitude(experiment, stepnumber, datapointindex, experimenttobescaled, experimenttoscaleto):
    if len(experiment[experimenttobescaled]['x_data']) < stepnumber or stepnumber < 0:
        return False
    elif len(experiment[experimenttobescaled]['x_data'][stepnumber]) < datapointindex or datapointindex < 0:
        return False
    elif len(experiment) < experimenttobescaled or experimenttobescaled < 0:
        return False
    elif len(experiment) < experimenttoscaleto or experimenttoscaleto < 0:
        return False
    elif experimenttoscaleto == experimenttobescaled:
        return experiment
    targetamp = experiment[experimenttoscaleto]['y_data'][stepnumber][datapointindex]
    amptoscale = experiment[experimenttobescaled]['y_data'][stepnumber][datapointindex]

    for ystep in experiment[experimenttobescaled]['y_data']:
        for i in range(len(ystep)):
            ystep[i] = ystep[i]*(targetamp/amptoscale)

    return experiment


def splitsequentialbysteptype(experiment, sequence):
    new_experiment = []
    indexcontrol = []
    indexperexp = []
    seqlen = len(sequence)
    for exp in experiment:
        for i in range(len(exp['step_type'])):
                for j in range(seqlen):
                    if exp['step_type'][i] == sequence[j]:
                        if (j < seqlen-1 and exp['step_type'][i+1] == sequence[j+1]) or j == seqlen - 1:
                            indexperexp.extend([i])
        indexcontrol.append(indexperexp)
        indexperexp = []

    for i in range(len(experiment)):
        for j in range((len(indexcontrol[j])/seqlen)):
            SensorName = None
            SensorRole = None
            SensorType = None
            SensorInfo = None
            Xdata = []
            Ydata = []
            StepName = []
            ActualTime = []
            StepStatus = []
            StepType = []
            Concentration = []
            MolarConcentration = []
            SampleID = []
            WellType = []
            MW = []
            SN = None
            Com = []
            Flags = []
            FXN = []
            FitParam = []
            ParamName = []
            for k in range(len(sequence)):
                if experiment[i]['sensor'] is not None:
                    SensorName = experiment[i]['sensor']
                if experiment[i]['sensor_role'] is not None:
                    SensorRole = experiment[i]['sensor_role']
                if experiment[i]['sensor_type'] is not None:
                    SensorType = experiment[i]['sensor_type']
                if experiment[i]['sensor_info'] is not None:
                    SensorInfo = experiment[i]['sensor_info']
                if experiment[i]['x_data'] is not None:
                    Xdata.append(experiment[i]['x_data'][indexcontrol[i][k + (j * seqlen)]])
                if experiment[i]['y_data'] is not None:
                    Ydata.append(experiment[i]['y_data'][indexcontrol[i][k + (j * seqlen)]])
                if experiment[i]['step_name'] is not None:
                    StepName.append(experiment[i]['step_name'][indexcontrol[i][k + (j * seqlen)]])
                if experiment[i]['actual_time'] is not None:
                    ActualTime.append(experiment[i]['actual_time'][indexcontrol[i][k + (j * seqlen)]])
                if experiment[i]['step_status'] is not None:
                    StepStatus.append(experiment[i]['step_status'][indexcontrol[i][k + (j * seqlen)]])
                if experiment[i]['step_type'] is not None:
                    StepType.append(experiment[i]['step_type'][indexcontrol[i][k + (j * seqlen)]])
                if experiment[i]['concentration'] is not None:
                    Concentration.append(experiment[i]['concentration'][indexcontrol[i][k + (j * seqlen)]])
                if experiment[i]['molarconcentration'] is not None:
                    MolarConcentration.append(experiment[i]['molarconcentration'][indexcontrol[i][k + (j * seqlen)]])
                if experiment[i]['sampleid'] is not None:
                    SampleID.append(experiment[i]['sampleid'][indexcontrol[i][k + (j * seqlen)]])
                if experiment[i]['welltype'] is not None:
                    WellType.append(experiment[i]['welltype'][indexcontrol[i][k + (j * seqlen)]])
                if experiment[i]['molecularweight'] is not None:
                    MW.append(experiment[i]['molecularweight'][indexcontrol[i][k + (j * seqlen)]])
                if experiment[i]['signal2noise'] is not None:
                    SN = experiment[i]['signal2noise']
                try:
                    Com = experiment[i]['comments']
                except:
                    Com = []
                try:
                    Flags = experiment[i]['flags']
                except:
                    Flags = []
                try:
                    FXN = experiment[i]['fit_fxn'][indexcontrol[i][k + (j * seqlen)]]
                except:
                    FXN = []
                try:
                    FitParam = experiment[i]['fit_param'][indexcontrol[i][k + (j * seqlen)]]
                except:
                    FitParam = []
                try:
                    ParamName = experiment[i]['param_name'][indexcontrol[i][k + (j * seqlen)]]
                except:
                    ParamName = []

            new_experiment.append({'sensor': SensorName, 'sensor_role': SensorRole, 'sensor_type': SensorType,
                'sensor_info': SensorInfo, 'x_data': Xdata, 'y_data': Ydata, 'step_name': StepName,
                'actual_time': ActualTime, 'step_status': StepStatus, 'step_type': StepType,
                'concentration': Concentration, 'molarconcentration': MolarConcentration,
                'sampleid': SampleID, 'welltype': WellType, 'molecularweight': MW, 'signal2noise': SN,
                'comments': Com, 'flags': Flags, 'fit_fxn': FXN, 'fit_param': FitParam, 'param_name': ParamName})

    return new_experiment