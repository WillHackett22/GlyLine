#GlyLine Signal class
#these functions parse the extracted data
def CWTMaker(signalvector,scalemin=2,scalemax=31):
    cwtmatr=signal.cwt(signalvector,signal.ricker,np.arange(scalemin,scalemax))
    cwtpeak=signal.find_peaks_cwt(signalvector,np.arange(scalemin,scalemax))
    return cwtmatr, cwtpeak

def tempcurveMaker(scale,lowerbound,upperbound,wave,cwtmass,pid,labels,gplab):
    tempcurve=pd.DataFrame(None,columns=['scanid','peakid','PeakIntensity','AddID','GPID'])
    tempcurve['scanid']=list(np.arange(lowerbound,upperbound))
    tempcurve['PeakIntensity']=list(wave/np.max(wave)*cwtmass)
    tempcurve['peakid']=[pid]*tempcurve.shape[0]
    tempcurve['AddID']=[labels]*tempcurve.shape[0]
    tempcurve['GPID']=[gplab]*tempcurve.shape[0]
    return tempcurve
                
def EmptyIndexFinder(theoreticalcurves,observeddata,AddIDTarget,exceptionthreshmin=3):
    d={'PeakIntensity':'SummedIntensity','peakid':'peakids'}
    grpedgp=theoreticalcurves.loc[theoreticalcurves['AddID']==AddIDTarget].groupby(['scanid']).agg({'PeakIntensity':'sum','peakid':'unique'}).rename(columns=d)
    #find the empty intensity scanids, first get the ones that are totally missing
    emptyids=list(set(observeddata['scanid'])-set(grpedgp.index.tolist()))
    #find the ones without any peakids assoc, find closest scans with peaks
    #if within exception distance, include, if not add to exception list
    for j in emptyids:
        emptydist=abs(emptyids-grpedgp.index)
        if np.min(emptydist)<=exceptionthreshmin:
            grpedgp=pd.concat([grpedgp,pd.DataFrame([[0,grpedgp['peakids'].iloc[emptydist.argmin()].tolist()]],columns=['SummedIntensity','peakids'],index=[j])])
    #now take the ones that are seen but have a total intensity of 0
    emptyout=grpedgp.loc[observeddata['scanid']].loc[(grpedgp['SummedIntensity']==0)].index.tolist()
    return emptyout

def TheoreticalCurveAdjuster(pkdatadf,cwtmatr,labels,gplab,curvesout):
    pids=curvesout.loc[(curvesout['scanid'].isin(emptyids)) & (curvesout['GPID']==labels),'peakid']
    #look at peak with most missing associations
    pidtgt=pids.mode()[0]
    pkdata=pkdatadf.loc[pkdatadf['peakid']==pidtgt]
    tempcurve=pd.DataFrame(None,columns=['scanid','peakid','PeakIntensity','AddID','GPID'])
    #if both pk ends are missing then increase scale
    #otherwise, just shift the peak
    if ((pkdata['startscan'].iloc[0] in emptyids) & (pkdata['endscan'].iloc[0] in emptyids)):
        scale=pkdata['scale'].iloc[0]+1
        lowerbound=pkdata['apex'].iloc[0]-scale
        upperbound=pkdata['apex'].iloc[0]+scale+1
        cwtmass=cwtmatr[scale,lowerbound:upperbound].sum()
        wave=list(signal.ricker(scale*2+1,scale))
    else:
        scale=pkdata['scale'].iloc[0]
        wave=list(signal.ricker(scale*2,scale))
        if (pkdata['startscan'].iloc[0] in emptyids):
            lowerbound=pkdata['startscan'].iloc[0]
            upperbound=pkdata['endscan'].iloc[0]-1
        else:
            lowerbound=pkdata['startscan'].iloc[0]+1
            upperbound=pkdata['endscan'].iloc[0]
        cwtmass=cwtmatr[scale,lowerbound:upperbound].sum()                
    tempcurve=tempcurveMaker(scale,lowerbound,upperbound,wave,cwtmass,pid,labels,gplab)
    pkdatadf.loc[pkdatadf['peakid']==pidtgt]=[pidtgt,pkdata['AddID'].iloc[0],scale,lowerbound,upperbound,pkdata['apex'].iloc[0],cwtmass]
    curvesout=curvesout.loc[curvesout['peakid']!=pidtgt]
    curvesout=pd.concat([curvesout,tempcurve])
    

def TheoreticalCurveGenerator(cwtmatr,pkn,pk,ms1temp,curvesout,pkdatadf):
    pid="pid"+str(labels)+"_"+str(pkn)
    scale=cwtmatr[:,pk].argmax()+scalemin
    lowerbound=ms1temp.index.values.min()+pk-scale
    upperbound=ms1temp.index.values.min()+pk+scale+1
    cwtmass=cwtmatr[scale,lowerbound:upperbound].sum()
    wave=list(signal.ricker(scale*2+1,scale))
    tempcurve=tempcurveMaker(scale,lowerbound,upperbound,wave,cwtmass,pid,labels,gplab)
    curvesout=pd.concat([curvesout,tempcurve])
    pkdatadf.loc[len(pkdatadf.index)]=[pid,labels,scale,lowerbound,upperbound,pk,cwtmass]

def TheoreticalCurvesOfAddID(labels,dfi,curvesout,pkdatadf,PreIDMin,PreIDMax,imputetype=None,basesignal=0):
    gplab=labels #CHANGE LINES FOR AddID vs GPID changes
    ms1temp=Helper.BasicImputation(dfi,dataname='intensity',refname='scanid',RefMin=PreIDMin,RefMax=PreIDMax,imputetype=imputetype,basesignal=basesignal)
    ms1temp['intensity']=ms1temp['intensity']+0.000001 # replace with basesignal variable in helper.basicimputation
    cwtmatr, cwtpeak=CWTMaker(ms1temp['intensity'].tolist())
    if len(cwtpeak)>0:
        for pkn, pk in enumerate(cwtpeak):
            if np.max(cwtmatr[:,pk])>1:
                TheoreticalCurveGenerator(cwtmatr,pkn,pk,ms1temp,curvesout,pkdatadf)
    else:
        pid='pid'+str(labels)+'_0'
        posidx=ms1temp.loc[ms1temp['intensity']>0].index.tolist()
        tempcurve=pd.DataFrame(None,columns=['scanid','peakid','PeakIntensity','AddID','GPID'])
        tempcurve['scanid']=ms1temp['scanid'].loc[posidx].tolist()
        tempcurve['PeakIntensity']=ms1temp['intensity'].loc[posidx].tolist()
        tempcurve['peakid']=[pid]*tempcurve.shape[0]
        tempcurve['AddID']=[labels]*tempcurve.shape[0]
        tempcurve['GPID']=[gplab]*tempcurve.shape[0]
        curvesout=pd.concat([curvesout,tempcurve])
        pkdatadf.loc[len(pkdatadf.index)]=[pid,labels,np.nan,np.min(posidx),np.max(posidx),np.nan,tempcurve['PeakIntensity'].sum()]
    emptyids=EmptyIndexFinder(curvesout,dfi,labels)
    while len(emptyids)>0:
        TheoreticalCurveAdjuster(pkdatadf,cwtmatr,labels,gplab,curvesout)
        #repeat until empty values gone 
        emptyids=EmptyIndexFinder(curvesout,dfi,labels)

def TheoreticalCurveGenerationMain(observed,PreIDMin=None,PreIDMax=None,imputetype=None,basesignal=0):
    curvesout=pd.DataFrame(None,columns=['scanid','peakid','PeakIntensity','AddID','GPID'])
    pkdatadf=pd.DataFrame(None,columns=['peakid','AddID','scale','startscan','endscan','apex','peakmass'])
    for labels, dfi in observed.groupby("AddID"):
        TheoreticalCurvesOfAddID(labels,dfi,curvesout,pkdatadf,PreIDMin,PreIDMax,imputetype,basesignal)