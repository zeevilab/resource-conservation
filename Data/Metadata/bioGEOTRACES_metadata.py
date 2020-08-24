import requests
import lxml.html as lh
from pandas.core.frame import DataFrame
import numpy as np
import re
from pandas.core.series import Series
from os.path import join, basename
from pandas.io.parsers import read_csv
from pandas.io.pickle import read_pickle
from datetime import datetime, timedelta
from glob import glob
from xarray.backends.api import open_dataset
from Data.config import Mapping, Biodata
from pandas.core.reshape.concat import concat

def _get_latlon(x):
    m = re.search(r"(\d+\.?\d*)\s?([NS])\s(\d+\.?\d*)\s?([EW])", x)
    try:
        lat, ns, lon, ew = m.groups()
    except:
        print('jackie')
    return (-float(lat) if ns == 'S' else float(lat)), (-float(lon) if ew == 'W' else float(lon))

def extract(raw_fastq, biodata):
    '''
    Code based on https://towardsdatascience.com/web-scraping-html-tables-with-python-c9baba21059
    '''
    url = 'https://www.nature.com/articles/sdata2018176/tables/3'
    #Create a handle, page, to handle the contents of the website
    page = requests.get(url)
    #Store the contents of the website under doc
    doc = lh.fromstring(page.content)
    #Parse data that are stored between <tr>..</tr> of HTML
    tr_elements = doc.xpath('//tr')
    df = DataFrame([[t.text_content().strip() for t in tr_e] for tr_e in tr_elements])
    df = df.T.set_index(0).T
    df = df.rename(columns={'NCBI SRA Accession':'RunID'})\
                .set_index('RunID')
    sample_sheet = read_csv(raw_fastq.SampleSheet, sep = '\t')\
                    .rename(columns={'run_accession':'RunID',
                                     'secondary_sample_accession':'SampleID'})\
                    .set_index('RunID')
    df = sample_sheet.join(df)
    df = df[df['Sample name'].notnull()]
    df['Fastq_1'] = df['fastq_ftp'].apply(lambda x:join(raw_fastq.FastqDir, 
                                                            basename(x.split(';')[0])))
    df['Fastq_2'] = df['fastq_ftp'].apply(lambda x:join(raw_fastq.FastqDir,
                                                            basename(x.split(';')[0])))
    df['ICRABAM_1'] =  df['fastq_ftp'].apply(lambda x:join(Mapping.OM_RGC.MapDir,
                                    basename(x.split(';')[0].replace('.fastq.gz','.icra.bam'))))
    df['ICRABAM_2'] =  df['fastq_ftp'].apply(lambda x:join(Mapping.OM_RGC.MapDir,
                                    basename(x.split(';')[1].replace('.fastq.gz','.icra.bam'))))
    df['Collection Date'] = df['Collection Date'].astype(np.datetime64)
    df['Depth (m)'] = df['Depth (m)'].astype(float)
    df[['Latitude','Longitude']] = df['Latitude and Longitude'].apply(_get_latlon).apply(Series)
    df = df.rename(columns = {'Sample name':'bioG_SampleID','Bottle ID':'Event_label',
                              'Cruise Station':'Station_label', 'Cruise ID':'CruiseID',
                              'GEOTRACES section':'GEOTRACES_section', 
                              'Cruise series':'Cruise_series','Depth (m)':'Depth',
                              'Collection Date':'Collection_datetime'})
    df = df[['SampleID','Cruise_series','Collection_datetime','Latitude','Longitude','Depth',
             'bioG_SampleID','Event_label','Station_label', 'GEOTRACES_section','CruiseID',
             'Fastq_1','Fastq_2','ICRABAM_1','ICRABAM_2']]
    df = df.reset_index().set_index(['SampleID','RunID'])
    df.to_pickle(biodata.metadataDF)

def _getrowmd(disc_df, row, depthcol='DEPTH [m]', match_section=True, depth_tolerance=None):
    d_row_t = disc_df[abs(disc_df.Collection_datetime - row.Collection_datetime) < timedelta(days=1)]
    if match_section:
        d_row_t = d_row_t[d_row_t.Cruise == row.GEOTRACES_section]
    d_row = d_row_t[abs(d_row_t[depthcol] - row.Depth) < (1 \
                    if depth_tolerance is None else depth_tolerance)]
    if d_row.shape[0] == 0:
        if d_row_t.shape[0] > 0:
            d_row_t = d_row_t.loc[abs(d_row_t[depthcol] - row.Depth).idxmin()]
            if abs(d_row_t[depthcol] - row.Depth) / row.Depth < 0.05:
                d_row = d_row_t
    elif d_row.shape[0] > 1 and depth_tolerance is None:
        d_row = d_row.loc[abs(d_row.Collection_datetime - row.Collection_datetime).idxmin()]
        if len(d_row.shape)>1 and d_row.shape[0] > 1:
            print('jackie')
    if type(d_row) == DataFrame and d_row.shape[0] == 1:
        d_row = d_row.iloc[0]
    d_row['Orig_datetime'] = row.Collection_datetime
    d_row['Orig_Latitude'] = row.Latitude
    d_row['Orig_Longitude'] = row.Longitude
    return d_row
    
def get_measurements_GEOTRACES(depth_tolerance=None):
    sample_md = read_pickle(Biodata.bioGEOTRACES.metadataDF)
    disc_df = read_csv(Biodata.bioGEOTRACES.DiscreteSampleTXT, sep = '\t')
    disc_df = disc_df.rename(columns={'yyyy-mm-ddThh:mm:ss.sss':'Collection_datetime'})
    disc_df.Collection_datetime = disc_df.Collection_datetime\
                                         .apply(lambda x:datetime.strptime(x,'%Y-%m-%dT%H:%M:%S'))
    disc_df = disc_df[disc_df.Collection_datetime > datetime(2000,1,1)]
    ret = {}
    for nm, row in sample_md.iterrows():
        d_row = _getrowmd(disc_df, row, depth_tolerance=depth_tolerance)
        if depth_tolerance is not None:
            if type(d_row) == Series:
                ret[(*nm, d_row.name)] = d_row[[i for i in d_row.index \
                        if (not i.startswith('QV') and (not i.startswith('STANDARD_DEV')))]]
            else:
                #if more than one date, take closest
                datediffs = (d_row.Collection_datetime.astype(np.datetime64)-d_row.Orig_datetime).apply(np.abs)
                d_row = d_row.loc[datediffs[datediffs==datediffs.min()].index]
                for rownm, row in d_row.iterrows():
                    ret[(*nm, rownm)] = row[[i for i in row.index \
                        if (not i.startswith('QV') and (not i.startswith('STANDARD_DEV')))]]
        else:
            ret[nm] = d_row[[i for i in d_row.index \
                            if (not i.startswith('QV') and (not i.startswith('STANDARD_DEV')))]]
    if depth_tolerance is None:
        DataFrame(ret).dropna(how='all').T.groupby(level=0).first()\
            .to_pickle(Biodata.bioGEOTRACES.SampleMeasurementsDF)
    else:
        df = DataFrame(ret).dropna(how='all').T.reset_index()\
            .rename(columns={'level_0':'SampleID','level_2':'MeasurementID'})\
            .set_index(['SampleID','MeasurementID']).drop('level_1',axis=1)
        df.to_pickle(Biodata.bioGEOTRACES.SampleMeasurementsDF\
                     .replace('.df','.tol_{}.df'.format(depth_tolerance)))

def get_measurements_HOT():
    hotmd = read_pickle(Biodata.ALOHA_BATS.metadataDF)
    hotmd = hotmd[hotmd.Cruise_series == 'HOT']
    def read_ds(f):
        ds = open_dataset(f, decode_times=False)
        df = ds.to_dataframe().reset_index()
        if ds.time_coverage_start != ds.time_coverage_end:
            raise
        dsdt = datetime.strptime(ds.time_coverage_start, '%Y-%m-%dT%H:%M:%SZ')
        if abs(hotmd.Collection_datetime-dsdt).min() > timedelta(days=2):
            return None
        df['TIME'] = dsdt
        return df
    ctd = concat([read_ds(f) for f in glob(join(Biodata.ALOHA_BATS.ALOHACTD, '*.nc'))])
    ctd = ctd.rename(columns={'TIME':'Collection_datetime'}).reset_index()
    water = concat([read_ds(f) for f in glob(join(Biodata.ALOHA_BATS.ALOHAWater, '*.nc'))])
    water = water.rename(columns={'TIME':'Collection_datetime'}).reset_index()
    ret = {}
    for nm, row in hotmd.iterrows():
        watrow = _getrowmd(water, row, 'DEPTH', False)
        ctdrow = _getrowmd(ctd, row, 'DEPTH', False)
        watrow.index = [i.replace('Orig','BotOrig') for i in watrow.index]
        ret[nm[0]] = concat([ctdrow, watrow[11:]])
    df = DataFrame({k:v for k,v in ret.items() if type(v) == Series}).T
    df = df[[c for c in df.columns if not c.endswith('_QC')]]
    df.replace('nan',np.nan).drop('index', axis=1).to_pickle(Biodata.ALOHA_BATS.ALOHASampleMeasurementsDF)

def get_measurements_BATS():
    ctd = concat([read_csv(f, header=None, sep='\t') \
                  for f in glob(join(Biodata.ALOHA_BATS.BATSASCIIDir, '*ctd.txt'))],
                 sort = False).dropna().replace(-999, np.nan)
    ctd.columns = ['ID','Collection_datetime','Latitude','Longitude','Pressure [dbar]',
                   'Depth','Temperature [c]','Conductivity [S/m]','Salinity','Oxygen [umol/kg]',
                   'Beam Attenuation Coefficient [1/m]','Flourescence','PAR [uE/m2/s]']
    def todatetime(x):
        year = int(x)
        rem = x-year
        basedt = datetime(year,1,1)
        return basedt + timedelta(seconds=(basedt.replace(year=year+1)-basedt).total_seconds()*rem)
    ctd.Collection_datetime = ctd.Collection_datetime.apply(todatetime)
    bot = read_csv(Biodata.ALOHA_BATS.BATSBottle, header=Biodata.ALOHA_BATS.BATSBottleHeader, 
                   sep='\t').reset_index().replace(-999, np.nan)
    bot.columns = ['ID','yyyymmdd','Collection_datetime','time' ,'Latitude','Longitude','Depth',
                   'Temp','CTD_S','Sal1','Sig-th','O2_1','OxFix','Anom1','CO2','Alk',
                   'Nitrate+Nitrite','Nitrite','Phosphate','Silicate','POC','PON','TOC','TN','Bact',
                   'POP','TDP','SRP','BSi','LSi','Pro','Syn','Piceu','Naneu']
    bot.Collection_datetime = bot.Collection_datetime.apply(todatetime)
    batsmd = read_pickle(Biodata.ALOHA_BATS.metadataDF)
    batsmd = batsmd[batsmd.Cruise_series == 'BATS']
    ret = {}
    for nm, row in batsmd.iterrows():
        botrow = _getrowmd(bot, row, 'Depth', False)
        ctdrow = _getrowmd(ctd, row, 'Depth', False)
        botrow.index = [i.replace('Orig','BotOrig') for i in botrow.index]
        ret[nm[0]] = concat([ctdrow, botrow[10:]])
    DataFrame({k:v for k,v in ret.items() if type(v) == Series}).T\
        .to_pickle(Biodata.ALOHA_BATS.BATSSampleMeasurementsDF)
