### This file is the GlyLine helper function file
## It contains functions that are more broadly applicable
import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression

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

def LogRelativizeIntensity(DF,ratio=.9):
    logint=np.log(DF['Intensity']+1)
    DF['LR_Intensity']=(logint-np.min(logint)*.9)/(np.max(logint)-np.min(logint)*.9)

#turn a vector into proportions
def SumRatio(vec):
    return list(vec/np.sum(vec))

#determine fragment type    
def FragmentType(x):
    pepbool=x['fragment_name'].str.contains(pat='y[0-9]$|b[0-9]$|y[0-9][0-9]$|b[0-9][0-9]$|peptide$')
    stubbool=x['fragment_name'].str.contains(pat='y[0-9]\\+|b[0-9]\\+|y[0-9][0-9]\\+|b[0-9][0-9]\\+|peptide\\+')
    glybool=np.logical_not(np.logical_or(pepbool,stubbool))
    x['fragment_type']=['str']*len(x['fragment_name'])
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

def BoundIndex(peakvalue,targetdf):
    ID=targetdf.index[np.where((targetdf['upper']>=peakvalue) & (targetdf['lower']<=peakvalue))].values
    return ID

def binsearch(arr,x):
    low = 0
    high = len(arr) - 1
    mid = 0
    while low <= high:
        mid = (high + low) // 2
        if arr[mid] < x:
            low = mid + 1
        elif arr[mid] > x:
            high = mid - 1
        else:
            return mid
    return -1

def ConfirmedGPIdxMaker(PSMMetaMaster,MS1Data,RunID=None,margin=18):
    if RunID is None:
        confirmed_gps=PSMMetaMaster['GPID','scan_time']
    else:
        confirmed_gps=PSMMetaMaster.loc[PSMMetaMaster['RunID']==RunID,['GPID','scan_time']]
    timebounds=pd.DataFrame(None,columns=['lower','upper'])
    timebounds['lower']=MS1Data['scan_time'].tolist()
    timebounds['upper']=MS1Data['scan_time'].tolist()[1:]+[MS1Data['scan_time'].max()+MS1Data['scan_time'].min()]
    timebounds.index=MS1Data.index
    tidx=[]
    for j in confirmed_gps['scan_time']:
        tidx+=[BoundIndex(j,timebounds)[0]]
    confirmed_gps['PrecursorIdx']=tidx
    confirmed_times=pd.DataFrame([[0,0]],columns=['GPID','PrecursorIdx'])
    for cgidx in confirmed_gps.index:
        tempdf=pd.DataFrame(None,columns=['GPID','PrecursorIdx'])
        t=confirmed_gps.loc[cgidx,'PrecursorIdx']
        tempdf['PrecursorIdx']=np.arange(t-margin,t+margin)
        tempdf['GPID']=confirmed_gps.loc[cgidx,'GPID']
        confirmed_times=pd.concat([confirmed_times,tempdf])
    confirmed_times=confirmed_times.drop(0).reset_index().drop('index',axis=1)
    return confirmed_times

def MassCalc(dfIon,decimals=4):
    mass1=-dfIon.groupby(['glycopeptide','fragment_name'])['mass_accuracy_ppm'].apply(np.unique)*dfIon.groupby(['glycopeptide','fragment_name'])['peak_mass'].apply(np.unique)/1000000+dfIon.groupby(['glycopeptide','fragment_name'])['peak_mass'].apply(np.unique)
    dfTemp=pd.DataFrame(None)
    dfTemp['neutralmass']=[np.round(np.median(m),decimals) for m in mass1.values]
    dfTemp['glycopeptide']=[g[0] for g in mass1.keys()]
    dfTemp['fragment_name']=[g[1] for g in mass1.keys()]
    return dfTemp

def MZfromNM(nm,z,Hmass=1.00797,pos=True):
    if pos==True:
        out=(nm+z*Hmass)/z
    else:
        out=((nm-z*Hmass)/z)
    return out

def PeptideReducer(dfMeta):
    peps=[gpep.split("{")[0] for gpep in dfMeta['glycopeptide']]
    dfMeta['peptide']=peps



def BasicImputation(dataframe,dataname='Intensity',refname='PrecursorIdx',RefMin=None,RefMax=None,imputetype=None,linearspan=3,basesignal=0):
    outdf=pd.DataFrame(None,columns=[dataname,refname])
    if RefMin==None:
        RefMin=dataframe[refname].min()-1
    if RefMax==None:
        RefMax=dataframe[refname].max()+2
    outdf[refname]=[r for r in range(RefMin,RefMax)]
    outdf=outdf.set_index(refname)
    outdf[refname]=outdf.index.tolist()
    temp=[basesignal]*outdf.shape[0]
    for idx,j in enumerate(outdf.index.tolist()):
        if j in dataframe[refname].tolist():
            temp[idx]=dataframe.loc[dataframe[refname]==j,dataname].tolist()[0]
    if imputetype!=None:
        temphold=[0]+temp.copy()+[0]
        tempf=temphold.copy()
        tempb=temphold.copy()
        for t in range(len(temp)):
            if temp[t]==basesignal:
                tempf[t+1]=(tempf[t]+tempf[t+2])/2
            if temp[-(t+1)]==basesignal:
                tempb[-(t+2)]=(tempb[-(t+1)]+tempb[-(t+3)])/2
        if imputetype=='Average':
            tempa=[(tempb[t]+tempf[t])/2 for t in range(len(temphold))]
            temp=tempa[1:(len(tempa)-1)]
        elif imputetype=='Higher':
            temph=[np.max([tempb[t],tempf[t]]) for t in range(len(temphold))]
            temp=temph[1:(len(temph)-1)]
        elif imputetype=='Lower':
            temph=[np.min([tempb[t],tempf[t]]) for t in range(len(temphold))]
            temp=temph[1:(len(temph)-1)]
        elif imputetype=='Linear':
            templin=temphold.copy()
            holdidx=[]
            for idx,val in enumerate(temphold):
                if val==basesignal:
                    holdidx+=[idx]
                elif (val!=basesignal) & (linearspan>=len(holdidx)>0):
                    startval=temphold[holdidx[0]-1]
                    endval=temphold[idx]
                    dx=(endval-startval)/(len(holdidx)+1)
                    for jdx in holdidx:
                        templin[jdx]=templin[jdx-1]+dx
                    holdidx=[]
                else:
                    for jdx in holdidx:
                        templin[jdx]=(tempb[jdx]+tempf[jdx])/2
                    holdidx=[]
            for jdx in holdidx:
                templin[jdx]=(tempb[jdx]+tempf[jdx])/2
            temp=templin[1:(len(templin)-1)]
    outdf[dataname]=temp
    return outdf

def LinearModelGenerator(curvedf,tgtpids,positive=True):
    tempdf=curvedf.iloc[np.argwhere(curvedf[['Intensity']+tgtpids].sum(axis=1)>0).flatten(),]
    y=tempdf['Intensity'].array
    X=np.array([[0]]*tempdf.shape[0])
    for tgt in tgtpids:
        X=np.concatenate((X,np.array(tempdf[tgt].tolist()).reshape([-1,1])),axis=1)
    m,n=X.shape
    X= X[np.arange(n) != np.array([0]*m)[:,None]].reshape(m,-1)
    model=LinearRegression(fit_intercept=False,positive=True).fit(X,y)
    return model, X, y  