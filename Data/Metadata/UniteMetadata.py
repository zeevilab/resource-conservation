from pandas.core.reshape.concat import concat
from pandas.io.pickle import read_pickle
import numpy as np
from Data.config import Biodata

DEPTH = 'Depth_m'
TEMPERATURE = 'Temperature_degC'
SALINITY = 'Salinity'
OXYGEN = 'Oxygen_umolkg-1'
PHOSPHATE = 'Phosphate_umolkg-1'
SILICATE = 'Silicate_umolkg-1'
NITRATE = 'Nitrate_umolkg-1'
NITRIRA = 'NO2_N03_umolkg-1'

def unite_sampledata():
    concat([read_pickle(Biodata.bioGEOTRACES.metadataDF),
            read_pickle(Biodata.ALOHA_BATS.metadataDF),
            read_pickle(Biodata.TARA.metadataDF)], 
            sort=False).to_pickle(Biodata.United.metadataDF)

def unite_measurements():
    biog = read_pickle(Biodata.bioGEOTRACES.SampleMeasurementsDF)\
            .rename(columns = {'DEPTH [m]':DEPTH, 'CTDTMP [deg C]':TEMPERATURE, 
                               'CTDSAL':SALINITY, 'CTDOXY [umol/kg]':OXYGEN,
                               'PHOSPHATE_D_CONC_BOTTLE [umol/kg]':PHOSPHATE,
                               'SILICATE_D_CONC_BOTTLE [umol/kg]':SILICATE,
                               'NITRATE_D_CONC_BOTTLE [umol/kg]':NITRATE,
                               'NO2+NO3_D_CONC_BOTTLE [umol/kg]':NITRIRA})
    tara = read_pickle(Biodata.TARA.SampleMeasurementsDF)\
            .rename(columns = {'Depth, nominal':DEPTH, 'Temp [°C]_median':TEMPERATURE, 
                               'Sal_median':SALINITY, 'OXYGEN [µmol/kg]':OXYGEN,
                               'Phosphate_median':PHOSPHATE,  #umol/l - https://doi.pangaea.de/10.1594/PANGAEA.839233
                               'Silicate_median':SILICATE, #umol/l
                               '[NO3]- [µmol/l]_median':NITRATE,
                               'Nitrate and Nitrite_median':NITRIRA}) #umol/l
    tara[SILICATE] = tara[SILICATE].astype(float)
    bats = read_pickle(Biodata.ALOHA_BATS.BATSSampleMeasurementsDF)\
            .rename(columns = {'Depth':DEPTH, 'Temperature [c]':TEMPERATURE, 
                               'Salinity':SALINITY, 'Oxygen [umol/kg]':OXYGEN,
                               'Phosphate':PHOSPHATE, #umol/kg
                               'Silicate':SILICATE, #umol/kg
                               'Nitrate+Nitrite':NITRIRA}) #umol/kg
    bats[[DEPTH,TEMPERATURE,SALINITY,OXYGEN,PHOSPHATE,SILICATE,NITRIRA]] = \
            bats[[DEPTH,TEMPERATURE,SALINITY,OXYGEN,PHOSPHATE,SILICATE,NITRIRA]].astype(float) 
    alha = read_pickle(Biodata.ALOHA_BATS.ALOHASampleMeasurementsDF)\
            .rename(columns = {'DEPTH':DEPTH, 'TEMP':TEMPERATURE, 
                               'PSAL':SALINITY, 'DOXY1':OXYGEN,
                               'PO41':PHOSPHATE, #umol/kg
                               'SILC1':SILICATE, #umol/kg
                               'NO31':NITRATE}) #umol/kg
    tara[DEPTH] = tara[DEPTH].apply(lambda x: float(x) if '-' not in str(x) else np.nan)
    tara[[PHOSPHATE, SILICATE, NITRATE, NITRIRA]] = tara[[PHOSPHATE, SILICATE, NITRATE, NITRIRA]]\
            .astype(float)\
            .truediv(1 + tara['Sigma-theta [kg/m**3]_median']/1000, axis=0) #umol/l --> umol/kg
    concat([tara,bats,alha,biog], sort=False)\
            [[DEPTH,TEMPERATURE,SALINITY,OXYGEN,PHOSPHATE,SILICATE,NITRATE,NITRIRA]]\
            .to_pickle(Biodata.United.SampleMeasurementsDF)
    