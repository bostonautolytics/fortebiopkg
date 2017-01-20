import os
import errno
import xml.etree.ElementTree as xmlio
import base64
import array
import xlsxwriter as XL
import datafitting
import numpy as np
from io import BytesIO as BIO
import xlrd
import xlwt


def __init__(self):
    """ No parameters required """
    pass


def read(experimentdirectory, ignoreregenerationaandneutralization):
    if not os.path.exists(experimentdirectory):
        return False
    experiment = []

    for files in os.listdir(experimentdirectory):
        if files.endswith('.frd'):
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
            Flags =[]
            infile = xmlio.parse(os.path.join(experimentdirectory, files)).getroot()
            for expinfo in infile.findall('ExperimentInfo'):
                SensorName = expinfo.find('SensorName').text
                SensorType = expinfo.find('SensorType').text
                SensorRole = expinfo.find('SensorRole').text
                SensorInfo = expinfo.find('SensorInfo').text
            for kindata in infile.findall('KineticsData'):
                for stepdata in kindata.findall('Step'):
                    for commondata in stepdata.findall('CommonData'):
                        if ignoreregenerationaandneutralization:
                            WellType.append(commondata.find('WellType').text)
                            if WellType[-1] == 'REGENERATION' or WellType[-1] == 'NEUTRALIZATION':
                                WellType.pop()
                                continue
                            else:
                                WellType.pop()
                        Concentration.append(commondata.find('Concentration').text)
                        MolarConcentration.append(commondata.find('MolarConcentration').text)
                        SampleID.append(commondata.find('SampleID').text)
                        WellType.append(commondata.find('WellType').text)
                        MW.append(commondata.find('MolecularWeight').text)
                        Xdata.append(array.array('f', base64.decodestring(stepdata.find('AssayXData').text)))
                        Ydata.append(array.array('f', base64.decodestring(stepdata.find('AssayYData').text)))
                        StepName.append(stepdata.find('StepName').text)
                        ActualTime.append(stepdata.find('ActualTime').text)
                        StepStatus.append(stepdata.find('StepStatus').text)
                        StepType.append(stepdata.find('StepType').text)
            for status in StepStatus:
                if not status == 'OK':
                    Flags.append('Sensor:' + status)
                    break
            experiment.append({'sensor': SensorName, 'sensor_role': SensorRole, 'sensor_type': SensorType,
                               'sensor_info': SensorInfo, 'x_data': Xdata, 'y_data': Ydata, 'step_name': StepName,
                               'actual_time': ActualTime, 'step_status': StepStatus, 'step_type': StepType,
                               'concentration': Concentration, 'molarconcentration': MolarConcentration,
                               'sampleid': SampleID, 'welltype': WellType, 'molecularweight': MW, 'signal2noise': 0,
                               'comments': [], 'flags': Flags, 'fit_fxn': [], 'fit_param': [],
                               'param_name': []})
    return experiment


def readbysteptype(experimentdirectory, steps2read):
    if not os.path.exists(experimentdirectory):
        return False
    experiment = []

    for files in os.listdir(experimentdirectory):
        if files.endswith('.frd'):
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
            Flags = []
            SampleGroup = []
            StepLoc = []
            steptypeflag = 1
            infile = xmlio.parse(os.path.join(experimentdirectory, files)).getroot()
            for expinfo in infile.findall('ExperimentInfo'):
                SensorName = expinfo.find('SensorName').text
                SensorType = expinfo.find('SensorType').text
                SensorRole = expinfo.find('SensorRole').text
                SensorInfo = expinfo.find('SensorInfo').text
            for kindata in infile.findall('KineticsData'):
                for stepdata in kindata.findall('Step'):
                    for commondata in stepdata.findall('CommonData'):
                        temp = stepdata.find('StepType').text
                        if temp not in steps2read:
                            continue
                        Concentration.append(commondata.find('Concentration').text)
                        MolarConcentration.append(commondata.find('MolarConcentration').text)
                        SampleID.append(commondata.find('SampleID').text)
                        if commondata.find('SampleGroup') is not None:
                            SampleGroup.append(commondata.find('SampleGroup').text)
                        else:
                            SampleGroup.append(None)
                        WellType.append(commondata.find('WellType').text)
                        MW.append(commondata.find('MolecularWeight').text)
                        Xdata.append(array.array('f', base64.decodestring(stepdata.find('AssayXData').text)))
                        Ydata.append(array.array('f', base64.decodestring(stepdata.find('AssayYData').text)))
                        StepName.append(stepdata.find('StepName').text)
                        ActualTime.append(stepdata.find('ActualTime').text)
                        StepStatus.append(stepdata.find('StepStatus').text)
                        StepLoc.append(commondata.find('SampleRow').text + commondata.find('SampleLocation').text)
                        if steptypeflag % 2 == 0:
                            StepType.append(stepdata.find('StepType').text)
                            steptypeflag += 1
                        StepType.append(temp)
            for status in StepStatus:
                if not status == 'OK':
                    Flags.append('Sensor:' + status)
                    break
            experiment.append({'sensor': SensorName, 'sensor_role': SensorRole, 'sensor_type': SensorType,
                               'sensor_info': SensorInfo, 'x_data': Xdata, 'y_data': Ydata, 'step_name': StepName,
                               'actual_time': ActualTime, 'step_status': StepStatus, 'step_type': StepType,
                               'step_loc': StepLoc, 'sample_group': SampleGroup,
                               'concentration': Concentration, 'molarconcentration': MolarConcentration,
                               'sampleid': SampleID, 'welltype': WellType, 'molecularweight': MW, 'signal2noise': 0,
                               'comments': [], 'flags': [], 'fit_fxn': [], 'param_error': [],
                               'fit_param': [], 'param_name': []})
    return experiment


def experimentcompleted(experimentdirectory):
    if not os.path.exists(experimentdirectory):
        return False
    dirlist = dirs = [d for d in os.listdir(experimentdirectory) if os.path.isdir(os.path.join(experimentdirectory, d))]
    for each in dirlist:
        lockfiles = []
        fmffiles = []
        frdfiles = []
        if not each.endswith('.db'):
            for files in os.listdir(os.path.join(experimentdirectory, each)):
                if files.endswith('.frd'):
                    frdfiles.append(files)
                elif files.endswith('.fmf'):
                    fmffiles.append(files)
                elif files.endswith('.lock'):
                    lockfiles.append(files)
            if len(lockfiles) == 0 and len(fmffiles) == 1 and len(frdfiles) > 1:
                return os.path.join(experimentdirectory, each)
    return False


def exportallxy(experiment, outfile, overwrite=True, onefilepersensor=False, *args, **kwargs):
    offstep = 2
    onstep = 1
    if 'off_step' in kwargs.keys():
        offstep = kwargs['off_step']
    if 'on_step' in kwargs.keys():
        onstep = kwargs['on_step']
    outpath, outfile, fileext = filewritecheck(experiment, outfile)
    if not outpath:
        return
    for exp in experiment:
        sensor = exp['sensor']
        if outfile == '*':
            filename = os.path.join(outpath, sensor + "Results" + fileext)
        else:
            filename = os.path.join(outpath, outfile + '_' + sensor + fileext)
        if os.path.exists(filename) and overwrite:
            try:
                os.remove(filename)
            except:
                return False
    if onefilepersensor:
        alpha = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
        for i in range(len(alpha)):
            for j in range(12):
                ka = []
                kd = []
                conc = []
                kaerror = []
                kderror = []
                x_data = []
                y_data = []
                step_type = []
                actual_time = []
                model = []
                for exp in experiment:
                    sensor = str(alpha[i] + str(j+1))
                    if exp['sensor'] == sensor:
                        ka.append(exp['fit_param'][2])
                        kaerror.append(exp['param_error'][2])
                        kd.append(exp['fit_param'][1])
                        kderror.append(exp['param_error'][1])
                        conc.append(exp['fit_param'][3])
                        x_data.append(np.concatenate((exp['x_data'][onstep], exp['x_data'][offstep]), axis=0))
                        y_data.append(np.concatenate((exp['y_data'][onstep], exp['y_data'][offstep]), axis=0))
                        step_type.append(exp['step_type'])
                        actual_time.append(exp['actual_time'])
                        model.append(np.concatenate((datafitting.analytical1to1modelon(np.array(exp['x_data'][onstep]),
                            exp['fit_param']), datafitting.analytical1to1modeloff(np.array(exp['x_data'][offstep]),
                                                                                  exp['fit_param'])), axis=0))
                if len(ka) == 0:
                    continue
                if outfile == '*':
                    filename = os.path.join(outpath, sensor + "Results" + fileext)
                else:
                    filename = os.path.join(outpath, outfile + '_' + sensor + fileext)
                ofile = open(filename, "w")
                ofile.write('Conc1')
                for c in conc:
                    ofile.write("\t%s" % c)
                ofile.write('\nStart1')
                for t in actual_time:
                    ofile.write("\t%s" % t[0])
                ofile.write('\nStop1')
                for x in range(len(x_data)):
                    ofile.write("\t%s" % str(float(actual_time[x][0]) + float(x_data[x][-1])))
                ofile.write('\nKobs')
                for k in range(len(ka)):
                    ofile.write("\t%s" % str((float(ka[k]) * float(conc[k])) + float(kd[k])))
                ofile.write('\nError Kobs')
                for k in range(len(ka)):
                    ofile.write("\t%s" % str(((float(kaerror[k])/float(ka[k]))) + float(kderror[k])))
                ofile.write('\nKoff')
                for k in kd:
                    ofile.write("\t%s" % str(k))
                ofile.write('\nError Koff')
                for k in kderror:
                    ofile.write("\t%s" % str(k))
                ofile.write('\n')
                for k in range(len(ka)):
                    ofile.write('Time' + str(k+1) + '\tData' + str(k+1) + '\tFit' + str(k+1))
                    if k < len(ka)-1:
                        ofile.write('\t')
                    else:
                        ofile.write('\n')
                maxrows = 0
                for x in x_data:
                    if len(x) > maxrows:
                        maxrows = len(x)
                for m in range(maxrows):
                    for n in range(len(ka)):
                        try:
                            ofile.write(str(float(x_data[n][m]) + float(actual_time[n][0])) + '\t' + str(y_data[n][m]) +
                                    '\t' + str(float(model[n][m])))
                        except:
                            ofile.write("")
                        if n < len(ka)-1:
                            ofile.write('\t')
                        else:
                            ofile.write('\n')

                ofile.close()
    elif not onefilepersensor:
        if outfile == '*':
            outfile = "AllData"
        counter = -1
        filename = os.path.join(outpath, outfile + fileext)
        if os.path.exists(filename) and not overwrite:
            return False
        elif os.path.exists(filename):
            try:
                os.remove(filename)
            except OSError:
                return False
        x_data = experiment[0]['x_data']
        x_data = [item for sublist in x_data for item in sublist]
        xlength = len(x_data)
        allx = []
        ally = []
        sensor = []
        for exp in experiment:
            x_data = exp['x_data']
            x_data = [item for sublist in x_data for item in sublist]
            y_data = exp['y_data']
            y_data = [item for sublist in y_data for item in sublist]
            allx.append(x_data)
            ally.append(y_data)
            sensor.append(exp['sensor'])

        ofile = open(filename, "w")
        while counter < xlength:
            for i in range(len(allx)):
                if counter < 0:
                    ofile.write(str(sensor[i]) + "\t")
                    if i < len(allx)-1:
                        ofile.write('\t')
                else:
                    ofile.write(str(allx[i][counter]) + "\t" + str(ally[i][counter]))
                    if i < len(allx)-1:
                        ofile.write('\t')
            counter += 1
            ofile.write('\n')
        ofile.close()
    return True


def xlsxmodelreportnew(experiment, outfile):
    outpath, outfile, fileext = filewritecheck(experiment, outfile)
    if not outpath:
        return
    outfile = 'fit'
    fileext = '.xlsx'
    filename = os.path.join(outpath, outfile + fileext)
    wb = XL.Workbook(filename)
    ws = wb.add_worksheet("Result Table")
    ws.set_column(2, 2, 32)
    ws.set_row(0, 30)
    ws.write(0, 0, 'flags')
    ws.write(0, 1, 'comments')
    ws.write(0, 2, 'Plot Data and Model')
    ws.write(0, 3, 'Sensor')
    ws.write(0, 4, 'SampleID')
    ws.write(0, 5, 'Peptide')
    ws.write(0, 6, 'KD')
    ws.write(0, 7, 'Req')
    ws.write(0, 8, 'avg Rsq (ka)')
    ws.write(0, 9, 'signal to noise')
    ws.write(0, 10, 'Amplitude')
    for i in range(len(experiment[0]['param_name'])):
        ws.write(0, i+11, str(experiment[0]['param_name'][i]))
    ws.freeze_panes(1, 0)

    for i in range(len(experiment)):
        if experiment[i]['sensor_info'] is None:
            continue
        if str(experiment[i]['step_status'])[0].upper() == "ERROR":
            ws.write(i + 1, 1, "Sensor Error")
            continue
        figure = BIO(datafitting.saveplot(experiment[i], None, 1, True))
        ws.set_row(i+1, 132)
        ws.write(i + 1, 0, str(experiment[i]['flags']))
        ws.write(i + 1, 1, str(experiment[i]['comments']))
        ws.insert_image(i + 1, 2, 'figure.png', {'image_data': figure, 'x_scale': 0.3, 'y_scale': 0.3})
        ws.write(i + 1, 3, experiment[i]['sensor'])
        ws.write(i + 1, 4, experiment[i]['sampleid'][1])
        ws.write(i + 1, 5, '')
        eqkdtemp = float(experiment[i]['fit_param'][1])/float(experiment[i]['fit_param'][2])
        ws.write(i + 1, 6, eqkdtemp)
        ctemp = float(experiment[i]['fit_param'][3])
        rmtemp = float(experiment[i]['fit_param'][0])
        ws.write(i + 1, 7, (ctemp/(ctemp+eqkdtemp)*rmtemp))
        model = datafitting.analytical1to1modelon(np.array(experiment[i]['x_data'][1]), experiment[i]['fit_param'])
        rawy = np.array(np.array(experiment[i]['y_data'][1]))
        rsqtemp = np.sum(np.power(model - rawy, 2))
        rawy = np.sum(np.power(rawy - np.average(rawy), 2))
        ws.write(i + 1, 8, (1 - (rsqtemp/rawy)))
        ws.write(i + 1, 9, experiment[i]['signal2noise'])
        ws.write(i + 1, 10, experiment[i]['y_data'][1][-1] - experiment[i]['y_data'][1][0])
        for j in range(len(experiment[i]['fit_param'])):
            ws.write(i+1, j+11, experiment[i]['fit_param'][j])
    wb.close()
    for i in range(len(experiment)):
        if str(experiment[i]['step_status'])[0].upper() == "ERROR" or experiment[i]['sensor_info'] is None:
            ws.write(i + 1, 1, "Sensor Error")
            continue
    return


def xlsmodelreport(experiment, outfile, *args, **kwargs):
    offrange = None
    onrange = None
    figures = False
    xlsout = False
    rsq = "full"
    if 'off_xrange' in kwargs.keys():
        offrange = kwargs['off_xrange']
    if 'on_xrange' in kwargs.keys():
        onrange = kwargs['on_xrange']
    if 'r_sq' in kwargs.keys():
        rsq = kwargs['r_sq']
    if 'figures' in kwargs.keys():
        figures = kwargs['figures']

    outpath, outfile, fileext = filewritecheck(experiment, outfile)
    if not outpath:
        return
    if fileext == '.xls':
        fileext = '.xlsx'
        xlsout = True
    elif fileext == '.xlsx':
        xlsout = False
    filename = os.path.join(outpath, outfile + fileext)
    wb = XL.Workbook(filename)
    ws = wb.add_worksheet("Result Table")
    ws.set_column(33, 33, 32)
    ws.set_row(0, 30)
    ws.write(0, 0, 'Comments')
    ws.write(0, 1, 'Flags')
    ws.write(0, 2, 'Index')
    ws.write(0, 3, 'Color')
    ws.write(0, 4, 'Sensor Location')
    ws.write(0, 5, 'Sensor Type')
    ws.write(0, 6, 'Sensor Info')
    ws.write(0, 7, 'Replicate Group')
    ws.write(0, 8, 'Baseline Loc.')
    ws.write(0, 9, 'Assoc. (Sample) Loc.')
    ws.write(0, 10, 'Sample ID')
    ws.write(0, 11, 'Dissoc. Loc.')
    ws.write(0, 12, 'Loading Well Location')
    ws.write(0, 13, 'Loading Sample ID')
    ws.write(0, 14, 'Cycle')
    ws.write(0, 15, 'Conc. (nM)')
    ws.write(0, 16, 'Response')
    ws.write(0, 17, 'KD (M)')
    ws.write(0, 18, 'KD Error')
    ws.write(0, 19, 'kon(1/Ms)')
    ws.write(0, 20, 'kon Error')
    ws.write(0, 21, 'kdis(1/s)')
    ws.write(0, 22, 'kdis Error')
    ws.write(0, 23, 'Rmax')
    ws.write(0, 24, 'Rmax Error')
    ws.write(0, 25, 'kobs(1/s)')
    ws.write(0, 26, 'Req')
    ws.write(0, 27, 'Req/Rmax(%)')
    if rsq == "on":
        ws.write(0, 28, 'Assoc X^2')
        ws.write(0, 29, 'Assoc R^2')
    elif rsq == "off":
        ws.write(0, 28, 'Dissoc X^2')
        ws.write(0, 29, 'Dissoc R^2')
    else:
        ws.write(0, 28, 'Full X^2')
        ws.write(0, 29, 'Full R^2')
    ws.write(0, 30, 'SSG KD')
    ws.write(0, 31, 'SSG Rmax')
    ws.write(0, 32, 'SSG R^2')

    ws.freeze_panes(1, 0)
    skipped = 0
    for i in range(len(experiment)):
        if experiment[i]['sensor_info'] is None:
            skipped = skipped + 1
            continue
        if str(experiment[i]['step_status'])[0].upper() == "ERROR":
            ws.write(i + 1, 1, "Sensor Error")
            continue
        ws.write(i - skipped + 1, 0, str(experiment[i]['comments']))
        ws.write(i - skipped + 1, 1, str(experiment[i]['flags']))
        ws.write(i - skipped + 1, 2, int(i-skipped))
        ws.write(i - skipped + 1, 3, str('0'))
        ws.write(i - skipped + 1, 4, str(experiment[i]['sensor']))
        ws.write(i - skipped + 1, 5, str(experiment[i]['sensor_type']))
        ws.write(i - skipped + 1, 6, str(experiment[i]['sensor_info']))
        ws.write(i - skipped + 1, 7, str(experiment[i]['sample_group'][1]))
        for w in range(len(experiment[i]['step_type'])):
            if experiment[i]['step_type'][w] == 'BASELINE':
                ws.write(i - skipped + 1, 8, str(experiment[i]['step_loc'][w]))
            elif experiment[i]['step_type'][w] == 'ASSOC':
                ws.write(i - skipped + 1, 9, str(experiment[i]['step_loc'][w]))
                rawy = np.array(np.array(experiment[i]['y_data'][w]))
            elif experiment[i]['step_type'][w] == 'DISASSOC':
                ws.write(i - skipped + 1, 11, str(experiment[i]['step_loc'][w]))
                rawyoff = np.array(np.array(experiment[i]['y_data'][w]))
            elif experiment[i]['step_type'][w] == 'LOADING':
                ws.write(i - skipped + 1, 12, str(experiment[i]['step_loc'][w]))
                ws.write(i - skipped + 1, 13, str(experiment[i]['sampleid'][w]))

        if rsq == "full" or "on":
            rawyon = rawy
            modelon = datafitting.analytical1to1modelon(np.array(experiment[i]['x_data'][1]), experiment[i]['fit_param'])
            rsqtempon = np.sum(np.power(modelon - rawyon, 2))
            rawyon = np.sum(np.power(rawyon - np.average(rawyon), 2))
            if rsq == "on":
                rsqreport = rsqtempon
                rawyreport = rawyon
        if rsq == "full" or "off":
            modeloff = datafitting.analytical1to1modeloff(np.array(experiment[i]['x_data'][2]), experiment[i]['fit_param'])
            rsqtempoff = np.sum(np.power(modeloff - rawyoff, 2))
            rawyoff = np.sum(np.power(rawyoff - np.average(rawyoff), 2))
            if rsq == "off":
                rsqreport = rsqtempoff
                rawyreport = rawyoff
            else:
                rsqreport = (rsqtempoff+rsqtempon) / 2
                rawyreport = (rawyoff+rawyon) / 2

        ws.write(i - skipped + 1, 10, experiment[i]['sampleid'][1])
        ws.write(i - skipped + 1, 14, str(experiment[i]['fit_param'][4]) + " " + str(experiment[i]['fit_param'][5])+ " " + str(experiment[i]['fit_param'][7]))
        ws.write(i - skipped + 1, 15, float(experiment[i]['molarconcentration'][1]))
        ws.write(i - skipped + 1, 16, float(np.max(rawy) - np.min(rawy)))
        #Rmax, kd, ka, Cpro, m, c, X0, kds
        eqkdtemp = float(experiment[i]['fit_param'][1]) / float(experiment[i]['fit_param'][2])
        ctemp = float(experiment[i]['fit_param'][3])
        rmtemp = float(experiment[i]['fit_param'][0])
        reqtemp = (ctemp / (ctemp + eqkdtemp) * rmtemp)
        ws.write(i - skipped + 1, 17, eqkdtemp)
        ws.write(i - skipped + 1, 18, ((float(experiment[i]['param_error'][1])/float(experiment[i]['fit_param'][1]))**2+
                            (float(experiment[i]['param_error'][2])/float(experiment[i]['fit_param'][2]))**2))**0.5
        ws.write(i - skipped + 1, 19, float(experiment[i]['fit_param'][2]))
        ws.write(i - skipped + 1, 20, float(experiment[i]['param_error'][2]))
        ws.write(i - skipped + 1, 21, float(experiment[i]['fit_param'][1]))
        ws.write(i - skipped + 1, 22, float(experiment[i]['param_error'][1]))
        ws.write(i - skipped + 1, 23, rmtemp)
        #ws.write(i - skipped + 1, 23, (float(experiment[i]['fit_param'][2]) *
        #           (float(experiment[i]['molarconcentration'][1]) / 1E9)) + float(experiment[i]['fit_param'][1]))
        ws.write(i - skipped + 1, 24, np.sqrt(float(experiment[i]['param_error'][2]) *
                    (float(experiment[i]['molarconcentration'][1]) / 1E9)**2 + float(experiment[i]['param_error'][1])**2))
        ws.write(i - skipped + 1, 25, reqtemp)
        ws.write(i - skipped + 1, 26, (np.sqrt((float(experiment[i]['param_error'][1])/
                    float(experiment[i]['fit_param'][1]))**2 +
                    (float(experiment[i]['param_error'][2])/float(experiment[i]['fit_param'][2]))**2 +
                    (float(experiment[i]['param_error'][0])/float(experiment[i]['fit_param'][0]))**2)) / ctemp)
        ws.write(i - skipped + 1, 27, 100*reqtemp/rmtemp)
        ws.write(i - skipped + 1, 28, str('')) #association chi sq
        ws.write(i - skipped + 1, 29, (1 - (rsqreport / rawyreport)))
        ws.write(i - skipped + 1, 30, str(float(experiment[i]['signal2noise'])))
        ws.write(i - skipped + 1, 31, str(''))
        ws.write(i - skipped + 1, 32, str(experiment[i]['fit_param'][7]))
        if figures:
            pngfigure = BIO(datafitting.saveplot(experiment[i], None, 1, True, on_xrange=onrange, off_xrange=offrange))
            ws.set_row(i - skipped + 1, 132)
            ws.insert_image(i - skipped + 1, 33, 'figure.png', {'image_data': pngfigure,
                                                                'x_scale': 0.3, 'y_scale': 0.3})
    wb.close()
    if xlsout:
        convertxlsx2xls(filename)
        if not figures:
            os.remove(filename)
    return


def convertxlsx2xls(filename):
    in_xlsx = xlrd.open_workbook(filename)
    out_xls = xlwt.Workbook()
    sheet_names = in_xlsx.sheet_names()
    for sheet_index in range(0, len(sheet_names)):
        sheet = out_xls.add_sheet(sheet_names[sheet_index])
        for row in range(0, in_xlsx.sheet_by_index(sheet_index).nrows):
            for col in range(0, in_xlsx.sheet_by_index(sheet_index).ncols):
                sheet.write(row, col, in_xlsx.sheet_by_index(sheet_index).cell_value(row, col))
    out_xls.save(filename[:len(filename)-1])
    return


def filewritecheck(experiment, outfile):
    if experiment == []:
        return False
    fileext = os.path.splitext(outfile)[1]
    outpath, outfile = os.path.split(os.path.abspath(outfile))
    if fileext == '':
        fileext = '.dat'
        outpath = os.path.join(outpath, outfile)
        outfile = 'data' + fileext
    if not os.path.exists(outpath):
        try:
            os.makedirs(outpath)
        except OSError:
            return False
    outfile = outfile[:len(outfile) - len(fileext)]
    return outpath, outfile, fileext
