### This file is the GlyLine helper function file
## It contains functions that are more broadly applicable
import pandas as pd
import numpy as np

def GSoftCSVRead(filename,subset=None,index=['scan_id'],PSMBool=True):
    temp=pd.read_csv(filename,index_col=index)
    if PSMBool:
        temp=temp.loc[temp['is_best_match']]
        
    if subset!=None:
        temp_reduce=temp.loc[~temp.duplicated(subset=subset,keep='first')]
    else:
        temp_reduce=temp
    
    return temp_reduce

#get the times from df1 onto df2
def TimeMatch(df1,df2,labelt='scan_time'):
    time_dict=df1[labelt].to_dict()
    return [time_dict[i] for i in df2.index.values]

#trapezoidal sum of product ions
def trap_sum_ion(x,label1='peak_intensity',labelt='scan_time'):
    temp=np.trapz(x.sort_values(labelt)[label1],x.sort_values(labelt)[labelt])
    return pd.Series(temp, index=['AUC'])

#trapezoidal sum of precursor ions
def trap_sum_psm(x,label1='precursor_abundance',labelt='scan_time'):
    temp=np.trapz(x.sort_values(labelt)[label1],x.sort_values(labelt)[labelt])
    return pd.Series(temp, index=['AUC'])

def trap_sum(x,label1,labelt='scan_time'):
    temp=np.trapz(x.sort_values(labelt)[label1],x.sort_values(labelt)[labelt])
    return pd.Series(temp, index=['AUC'])

#determine fragment type    
def FragmentType(x):
    pepbool=x['fragment_name'].str.contains(pat='y[0-9]$|b[0-9]$|y[0-9][0-9]$|b[0-9][0-9]$|peptide$')
    stubbool=x['fragment_name'].str.contains(pat='y[0-9]\\+|b[0-9]\\+|y[0-9][0-9]\\+|b[0-9][0-9]\\+|peptide\\+')
    glybool=np.logical_not(np.logical_or(pepbool,stubbool))
    x['fragment_type']=[np.nan]*len(x['fragment_name'])
    x.loc[pepbool,'fragment_type']='Peptide'
    x.loc[stubbool,'fragment_type']='Stub'
    x.loc[glybool,'fragment_type']='Glycan'
    
def PreProdRatio(x,y):
    x['precursor_abundance']=x.index.map(y['precursor_abundance'].to_dict())
    x['abun_ratio']=x['precursor_abundance']/x['peak_intensity']
    x['log_abun_ratio']=np.log(x['precursor_abundance'])/np.log(x['peak_intensity'])

def AUC_PreProdCalc(PSMFile,IonFile,index=['scan_id'],labelt='scan_time',ProdGrouping=0):
    dfPSM=GSoftCSVRead(PSMFile,subset=['glycopeptide'],index=index)
    dfIon=GSoftCSVRead(IonFile,subset=['glycopeptide','fragment_name'],index=index,PSMBool=False)
    dfIon=dfIon.loc[dfIon.index.get_level_values(0).isin(dfPSM.index.values)]
    dfIon[labelt]=TimeMatch(dfPSM,dfIon)
    PreProdRatio(dfIon,dfPSM)
    PSMAUC=dfPSM.groupby(['glycopeptide']).apply(trap_sum,label1='precursor_abundance')
    FragmentType(dfIon)
    if ProdGrouping==0:
        #group by fragment name & glycopeptide
        ProdAUC=dfIon.groupby(['glycopeptide','fragment_name']).apply(trap_sum,label1='peak_intensity')
        ProdAUC['n_Obs']=dfIon.groupby(['glycopeptide','fragment_name']).size()
        ProdAUC['fragment_type']=ProdAUC.index.get_level_values(1).map(pd.Series(dfIon.fragment_type.values,index=dfIon.fragment_name).to_dict())
    elif ProdGrouping==1:
        #group by mean of fragment type & glycopeptide
        tempdfIon=pd.DataFrame(dfIon.groupby(['glycopeptide','scan_time','fragment_type']).apply(np.mean))
        ProdAUC=tempdfIon.groupby(['glycopeptide','fragment_type']).apply(trap_sum,label1='peak_intensity')
        ProdAUC['n_Obs']=dfIon.groupby(['glycopeptide','fragment_type']).size()    
    else:
        print('Please pick 0 or 1 for ProdGrouping for product ion name or type respectively.')
        quit()
    
    #add precursorAUC to eventual output dataframe
    ProdAUC['PrecursorAUC']=ProdAUC.index.get_level_values(0).map(PSMAUC['AUC'].to_dict())
    chargedict=dfIon.groupby(['glycopeptide','fragment_name'])['peak_charge'].apply(np.unique).to_dict()
    ProdAUC['prod_charge']=ProdAUC.index.map(chargedict)
    chargedict=dfPSM.groupby(['glycopeptide'])['charge'].apply(np.unique).to_dict()
    ProdAUC['precursor_charge']=ProdAUC.index.get_level_values(0).map(chargedict)
    return ProdAUC



def MassCalc(dfIon):
    mass1=-dfIon.groupby(['glycopeptide','fragment_name'])['mass_accuracy_ppm'].apply(np.unique)*dfIon.groupby(['glycopeptide','fragment_name'])['peak_mass'].apply(np.unique)/1000000+dfIon.groupby(['glycopeptide','fragment_name'])['peak_mass'].apply(np.unique)
    mass=[np.median(m) for m in mass1.values]
    gps=[g[0] for g in mass1.keys()]
    ions=[g[1] for g in mass1.keys()]
    return mass, gps, ions

def PeptideReducer(dfMeta):
    peps=[gpep.split("{")[0] for gpep in dfMeta['glycopeptide']]
    dfMeta['peptide']=peps
   
    