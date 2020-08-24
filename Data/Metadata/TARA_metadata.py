from pandas.io.excel import read_excel
from pandas.io.parsers import read_csv
from os.path import basename, join
from datetime import datetime
from pandas.io.pickle import read_pickle
from pandas.core.reshape.concat import concat
from Data.config import Biodata, RawFastq, Mapping

def _fix_header(dataline):
    def getstatfromline(line):
        return '_min' if line.startswith('minimum') else '_25p' if line.startswith('lower quartile') \
            else '_median' if line.startswith('median') else '_75p' if line.startswith('upper quartile') \
            else '_max' if line.startswith('maximum') else ''
    dataline = dataline.apply(getstatfromline)
    dataline = dataline.reset_index()
    dataline['index'] = dataline['index'].apply(lambda x:x.split('.')[0] if 'Sample ID' not in x else x)
    return dataline.sum(1).values

def readTARAxls(datafile, headerline):
    df = read_excel(datafile, header=headerline)
    df.columns = _fix_header(df.iloc[2])
    df = df.drop(range(4))
    df = df.rename(columns = {'Sample ID':'TARA_SampleID', 'Sample ID.2':'SampleID'})\
            .drop(['Sample ID.1', 'PARAMETER'], axis=1).set_index('SampleID')
    return df.loc[df.index != 'none']
    
def getbasemetadata():
    df = read_excel(Biodata.TARA.SeqContextXL, header =16)
    df = df.drop(range(4)).rename(columns={'Analysis ID.1':'experiment_accession'})\
            .set_index('experiment_accession')
    sample_sheet = read_csv(RawFastq.TARA.SampleSheet, sep = '\t').set_index('experiment_accession')
    df = sample_sheet.join(df)
    df['RunID'] = df['run_accession']
    df['SampleID'] = df['secondary_sample_accession']
    df = df[df['submitted_ftp'].apply(lambda x: '.fastq.gz' in x)]
    df['Fastq_1'] = df['submitted_ftp'].apply(lambda x:join(RawFastq.TARA.FastqDir, 
                                                            basename(x.split(';')[0])))
    df['Fastq_2'] = df['submitted_ftp'].apply(lambda x:join(RawFastq.TARA.FastqDir,
                                                            basename(x.split(';')[0])))
    df['ICRABAM_1'] =  df['submitted_ftp'].apply(lambda x:join(Mapping.OM_RGC.MapDir,
                                    basename(x.split(';')[0].replace('.fastq.gz','.icra.bam'))))
    df['ICRABAM_2'] =  df['submitted_ftp'].apply(lambda x:join(Mapping.OM_RGC.MapDir,
                                    basename(x.split(';')[1].replace('.fastq.gz','.icra.bam'))))
    df = df.rename(columns = {'Sample ID':'TARA_SampleID','Analysis label':'Analysis_label',
                              'Environmental feature':'Environmental_feature',
                              'Size fraction, lower threshold':'Size_min',
                              'Size fraction, upper threshold':'Size_max',
                              'Depth, top/min':'Depth_min','Depth, bottom/max':'Depth_max',
                              'Event label':'Event_label', 'Station label':'Station_label'})
    df['Collection_datetime'] = df['Event_label'].dropna()\
                                    .apply(lambda x: datetime.strptime(x.split('_')[1],
                                            '%Y%m%dT%H%MZ'))
    stations = read_csv(Biodata.TARA.Stations, sep = '\t')
    df['Latitude'] = df['Station_label'].dropna()\
                        .apply(lambda x: stations[stations['Station'] == x]['Latitude'].values[0])
    df['Longitude'] = df['Station_label'].dropna()\
                        .apply(lambda x: stations[stations['Station'] == x]['Longitude'].values[0])
    df['Depth'] = df[['Depth_min','Depth_max']].astype(float).mean(1)
    df['Cruise_series'] = 'TARA'
    df = df[['SampleID','RunID','Cruise_series','Collection_datetime','Latitude','Longitude','Depth',
             'TARA_SampleID','Analysis_label','Event_label','Station_label',
             'Environmental_feature','Size_min','Size_max','Depth_min','Depth_max',
             'Fastq_1','Fastq_2','ICRABAM_1','ICRABAM_2']]
    df = df.reset_index().set_index(['SampleID','RunID'])
    df.to_pickle(Biodata.TARA.metadataDF)
    return df['experiment_accession'].values
    
def get_measurements_TARA():
    tara_md = read_pickle(Biodata.TARA.metadataDF)
    tara_ixs = tara_md.groupby(level=0).first().index
    carbchem = readTARAxls(Biodata.TARA.CarbChemXL, Biodata.TARA.XLHeaderLine)
            
    nutrient = readTARAxls(Biodata.TARA.NutrientsXL, Biodata.TARA.XLHeaderLine)
    hplc = readTARAxls(Biodata.TARA.HPLCXL, Biodata.TARA.XLHeaderLine)
    sensors = read_csv(Biodata.TARA.DepthSensorCSV, sep = '\t', 
                       header=Biodata.TARA.DepthSensorHeader)
    sensors = sensors.rename(columns = {sensors.columns[2]:'SampleID', 
                                        sensors.columns[0]:'TARA_SampleID'})\
                     .drop(sensors.columns[1], axis = 1).set_index('SampleID')
    sensors.columns = [c.split(' (')[0] + ('_min' if '(minimum' in c \
                                           else '_25p' if '(lower quartile' in c \
                                           else '_median' if '(median' in c \
                                           else '_75p' if '(upper quartile' in c \
                                           else '_max' if '(maximum' in c \
                                           else '' if 'OXYGEN' in c and c.endswith(')') \
                                           else '!@#$%' if '(calculated' in c or '(Calculated' in c \
                                           else '') for c in sensors.columns]
    sensors = sensors[[c for c in sensors.columns if '!@#$%' not in c]]
    concat([sensors.loc[tara_ixs], 
            nutrient.loc[tara_md.groupby(level=0).first().index].loc[tara_ixs].iloc[:,16:],
            carbchem.loc[tara_md.groupby(level=0).first().index].loc[tara_ixs].iloc[:,16:],
            hplc.loc[tara_md.groupby(level=0).first().index].loc[tara_ixs].iloc[:,16:]], axis=1)\
        .to_pickle(Biodata.TARA.SampleMeasurementsDF)
