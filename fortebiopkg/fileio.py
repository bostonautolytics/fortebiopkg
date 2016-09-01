import os
import xml.etree.ElementTree as xmlio
import base64
import array
import xlsxwriter as XL
import datafitting
import numpy as np


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
                        WellType.append(commondata.find('WellType').text)
                        MW.append(commondata.find('MolecularWeight').text)
                        Xdata.append(array.array('f', base64.decodestring(stepdata.find('AssayXData').text)))
                        Ydata.append(array.array('f', base64.decodestring(stepdata.find('AssayYData').text)))
                        StepName.append(stepdata.find('StepName').text)
                        ActualTime.append(stepdata.find('ActualTime').text)
                        StepStatus.append(stepdata.find('StepStatus').text)
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
                               'concentration': Concentration, 'molarconcentration': MolarConcentration,
                               'sampleid': SampleID, 'welltype': WellType, 'molecularweight': MW, 'signal2noise': 0,
                               'comments': [], 'flags': [], 'fit_fxn': [],
                               'fit_param': [], 'param_name': []})
    return experiment


def experimentcompleted(experimentdirectory):
    if not os.path.exists(experimentdirectory):
        return False
    lockfiles = []
    frmfiles = []
    frdfiles = []
    for each in os.listdir(experimentdirectory):
        if not each.endswith('.db'):
            for files in os.listdir(os.path.join(experimentdirectory)):
                if files.endswith('.frd'):
                    frdfiles.append(files)
                elif files.endswith('.frm'):
                    frmfiles.append(files)
                elif files.endswith('.lock'):
                    lockfiles.append(files)
            if len(lockfiles) == 0 and len(frmfiles) == 1 and len(frdfiles) > 1:
                return os.path.join(experimentdirectory)
    return False


def exportallxy(experiment, outfile, overwrite=True, onefilepersensor=False):
    outpath, outfile, fileext = filewritecheck(experiment, outfile)
    if not outpath:
        return
    for exp in experiment:
        sensor = exp['sensor']
        filename = os.path.join(outpath, outfile + '_' + sensor + fileext)
        if os.path.exists(filename) and overwrite:
            try:
                os.remove(filename)
            except:
                return False
    if onefilepersensor:
        for exp in experiment:
            sensor = exp['sensor']
            filename = os.path.join(outpath, outfile + '_' + sensor + fileext)
            if os.path.exists(filename):
                success = False
                count = 1
                while not success:
                    count += 1
                    filename = os.path.join(outpath, outfile + '_' + sensor + '_' + str(count) +fileext)
                    if not os.path.exists(filename):
                        success = True
            ofile = open(filename, "w")
            x_data = exp['x_data']
            y_data = exp['y_data']

            for j in range(0, (len(x_data))):
                for k in range(0, len(x_data[j])):
                    ofile.write(str(x_data[j][k]) + "\t" + str(y_data[j][k]) + "\n")
            ofile.close()
    elif not onefilepersensor:
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


def xlsxmodelreport(experiment, outfile):
    outpath, outfile, fileext = filewritecheck(experiment, outfile)
    if not outpath:
        return
    outfile = 'fit'
    fileext = '.xlsx'
    filename = os.path.join(outpath, outfile + fileext)
    wb = XL.Workbook(filename)
    ws = wb.add_worksheet("Experiment")
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
        datafitting.saveplot(experiment[i], os.path.join(outpath + '\\plot'+str(i)+'.png'), 1)
        ws.set_row(i+1, 132)
        ws.write(i + 1, 0, str(experiment[i]['flags']))
        ws.write(i + 1, 1, str(experiment[i]['comments']))
        ws.insert_image(i + 1, 2, os.path.join(outpath + '\\plot' + str(i) + '.png'), {'x_scale': 0.3, 'y_scale': 0.3})
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
        os.remove(outpath + '\\plot'+str(i)+'.png')
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
            os.mkdir(outpath)
        except OSError:
            return False
    outfile = outfile[:len(outfile) - len(fileext)]
    return outpath, outfile, fileext
