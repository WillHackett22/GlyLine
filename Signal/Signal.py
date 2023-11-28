#GlyLine Signal class
#these functions parse the extracted data
import pandas as pd
import numpy as np
from scipy import signal
from GlyLine.Helper import Helper

class IonAllocationByTheoreticalCurve:
    def __init__(self,observed,observedion,scalemin=2,scalemax=31,exceptionthreshmin=3,PreIDMin=0,PreIDMax=None,imputetype=None,basesignal=0):
        self.observed=observed
        self.observedion=observedion
        self.scalemin=scalemin
        self.scalemax=scalemax
        self.PreIDMin=PreIDMin
        self.PreIDMax=PreIDMax
        self.imputetype=imputetype
        self.basesignal=basesignal
        self.exceptionthreshmin=exceptionthreshmin
        
    def main(self,IonGPAssoc,FragEffTable,UseFragEff=True):
        self.TheoreticalCurveGenerationMain()
        self.IonAllocationMain(FragEffTable,IonGPAssoc,UseFragEff)
        
    def TheoreticalCurveGenerationMain(self):
        self.curvesout=pd.DataFrame(None,columns=['scanid','peakid','PeakIntensity','AddID','GPID'])
        self.pkdatadf=pd.DataFrame(None,columns=['peakid','AddID','scale','startscan','endscan','apex','peakmass'])
        for labels, dfi in self.observed.groupby("AddID"):
            self.TheoreticalCurvesOfAddID(labels,dfi)
        self.pkdatadf['conversion']=[self.observed.loc[self.observed['AddID']==prow['AddID'],'intensity'].sum()/self.pkdatadf.loc[self.pkdatadf['AddID']==prow['AddID'],'peakmass'].sum() for index, prow in self.pkdatadf.iterrows()]
        self.curvesout['PeakIntensity']=[crow['PeakIntensity']*self.pkdatadf.loc[self.pkdatadf['peakid']==crow['peakid'],'conversion'].iloc[0] for index, crow in self.curvesout.iterrows()]    
            
    def CWTMaker(self,signalvector):
        self.cwtmatr=signal.cwt(signalvector,signal.ricker,np.arange(self.scalemin,self.scalemax))
        self.cwtpeak=signal.find_peaks_cwt(signalvector,np.arange(self.scalemin,self.scalemax))
    
    def tempcurveMaker(self,scale,lowerbound,upperbound,wave,cwtmass,pid,labels,gplab):
        tempcurve=pd.DataFrame(None,columns=['scanid','peakid','PeakIntensity','AddID','GPID'])
        tempcurve['scanid']=list(np.arange(lowerbound,upperbound))
        tempcurve['PeakIntensity']=list(wave/np.max(wave)*cwtmass)
        tempcurve['peakid']=[pid]*tempcurve.shape[0]
        tempcurve['AddID']=[labels]*tempcurve.shape[0]
        tempcurve['GPID']=[gplab]*tempcurve.shape[0]
        return tempcurve
                    
    def EmptyIndexFinder(self,dfi,AddIDTarget):
        d={'PeakIntensity':'SummedIntensity','peakid':'peakids'}
        grpedgp=self.curvesout.loc[self.curvesout['AddID']==AddIDTarget].groupby(['scanid']).agg({'PeakIntensity':'sum','peakid':'unique'}).rename(columns=d)
        #find the empty intensity scanids, first get the ones that are totally missing
        emptyids=list(set(dfi['scanid'])-set(grpedgp.index.tolist()))
        #find the ones without any peakids assoc, find closest scans with peaks
        #if within exception distance, include, if not add to exception list
        for j in emptyids:
            emptydist=abs(emptyids-grpedgp.index)
            if np.min(emptydist)<=self.exceptionthreshmin:
                grpedgp=pd.concat([grpedgp,pd.DataFrame([[0,grpedgp['peakids'].iloc[emptydist.argmin()].tolist()]],columns=['SummedIntensity','peakids'],index=[j])])
        #now take the ones that are seen but have a total intensity of 0
        emptyout=grpedgp.loc[dfi['scanid']].loc[(grpedgp['SummedIntensity']==0)].index.tolist()
        return emptyout
    
    def TheoreticalCurveAdjuster(self,labels,gplab,emptyids):
        pids=self.curvesout.loc[(self.curvesout['scanid'].isin(emptyids)) & (self.curvesout['GPID']==labels),'peakid']
        #look at peak with most missing associations
        pidtgt=pids.mode()[0]
        pkdata=self.pkdatadf.loc[self.pkdatadf['peakid']==pidtgt]
        tempcurve=pd.DataFrame(None,columns=['scanid','peakid','PeakIntensity','AddID','GPID'])
        #if both pk ends are missing then increase scale
        #otherwise, just shift the peak
        if ((pkdata['startscan'].iloc[0] in emptyids) & (pkdata['endscan'].iloc[0] in emptyids)):
            scale=pkdata['scale'].iloc[0]+1
            lowerbound=pkdata['apex'].iloc[0]-scale
            upperbound=pkdata['apex'].iloc[0]+scale+1
            cwtmass=self.cwtmatr[scale,lowerbound:upperbound].sum()
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
            cwtmass=self.cwtmatr[scale,lowerbound:upperbound].sum()                
        tempcurve=self.tempcurveMaker(scale,lowerbound,upperbound,wave,cwtmass,pkdata['peakid'].iloc[0],labels,gplab)
        self.pkdatadf.loc[self.pkdatadf['peakid']==pidtgt]=[pidtgt,pkdata['AddID'].iloc[0],scale,lowerbound,upperbound,pkdata['apex'].iloc[0],cwtmass]
        self.curvesout=self.curvesout.loc[self.curvesout['peakid']!=pidtgt]
        self.curvesout=pd.concat([self.curvesout,tempcurve])
        
    
    def TheoreticalCurveGenerator(self,labels,gplab,pkn,pk,ms1temp):
        pid="pid"+str(labels)+"_"+str(pkn)
        scale=self.cwtmatr[:,pk].argmax()+self.scalemin
        lowerbound=ms1temp.index.values.min()+pk-scale
        upperbound=ms1temp.index.values.min()+pk+scale+1
        cwtmass=self.cwtmatr[scale,lowerbound:upperbound].sum()
        wave=list(signal.ricker(scale*2+1,scale))
        tempcurve=self.tempcurveMaker(scale,lowerbound,upperbound,wave,cwtmass,pid,labels,gplab)
        self.curvesout=pd.concat([self.curvesout,tempcurve])
        self.pkdatadf.loc[len(self.pkdatadf.index)]=[pid,labels,scale,lowerbound,upperbound,pk,cwtmass]
    
    def TheoreticalCurvesOfAddID(self,labels,dfi):
        gplab=labels #CHANGE LINES FOR AddID vs GPID changes
        ms1temp=Helper.BasicImputation(dfi,dataname='intensity',refname='scanid',RefMin=self.PreIDMin,RefMax=self.PreIDMax,imputetype=self.imputetype,basesignal=self.basesignal)
        ms1temp['intensity']=ms1temp['intensity']+0.000001 # replace with basesignal variable in helper.basicimputation
        self.CWTMaker(ms1temp['intensity'].tolist())
        if len(self.cwtpeak)>0:
            for pkn, pk in enumerate(self.cwtpeak):
                if np.max(self.cwtmatr[:,pk])>1:
                    self.TheoreticalCurveGenerator(labels,gplab,pkn,pk,ms1temp)
        else:
            pid='pid'+str(labels)+'_0'
            posidx=ms1temp.loc[ms1temp['intensity']>0].index.tolist()
            tempcurve=pd.DataFrame(None,columns=['scanid','peakid','PeakIntensity','AddID','GPID'])
            tempcurve['scanid']=ms1temp['scanid'].loc[posidx].tolist()
            tempcurve['PeakIntensity']=ms1temp['intensity'].loc[posidx].tolist()
            tempcurve['peakid']=[pid]*tempcurve.shape[0]
            tempcurve['AddID']=[labels]*tempcurve.shape[0]
            tempcurve['GPID']=[gplab]*tempcurve.shape[0]
            self.curvesout=pd.concat([self.curvesout,tempcurve])
            self.pkdatadf.loc[len(self.pkdatadf.index)]=[pid,labels,np.nan,np.min(posidx),np.max(posidx),np.nan,tempcurve['PeakIntensity'].sum()]
        emptyids=self.EmptyIndexFinder(dfi,labels)
        while len(emptyids)>0:
            self.TheoreticalCurveAdjuster(labels,gplab,emptyids)
            #repeat until empty values gone 
            emptyids=self.EmptyIndexFinder(dfi,labels)
    
    def IonAllocationMain(self,FragEffTable,IonGPAssoc,UseFragEff=True):
        if UseFragEff!=True:
            FragEffTable['FragEff']=1
        self.adjion_int=pd.DataFrame(None,columns=['GPID','IonID','scanid','adj_intensity'])
        for scid, dfi in self.observedion.groupby(['scanid']):
            gplist=self.curvesout.loc[self.curvesout['scanid']==scid,['GPID','PeakIntensity']].copy()
            for ion in dfi['IonID']:
               self.IonAllocationForIndIonID(dfi,ion,scid,gplist,FragEffTable,IonGPAssoc)
    
    def IonAllocationForIndIonID(self,dfi,ion,scid,gplist,FragEffTable,IonGPAssoc):
        hits=gplist.loc[gplist['GPID'].isin(IonGPAssoc.loc[ion,'GPID'])].copy()
        hits['FragEfficiency']=FragEffTable.loc[(FragEffTable['GPID'].isin(hits['GPID']) & (FragEffTable['IonID']==ion)),'FragEff'].tolist()
        hits['FragRatio']=Helper.SumRatio(hits['FragEfficiency'])
        hits['Ratio']=Helper.SumRatio(hits['PeakIntensity']*hits['FragRatio'])
        tempint=pd.DataFrame(None,columns=['GPID','IonID','scanid','adj_intensity'])
        tempint['adj_intensity']=list(dfi.loc[dfi['IonID']==ion,'intensity'].tolist()*np.array(hits['Ratio']))
        tempint['GPID']=hits['GPID'].tolist()
        tempint['IonID']=ion
        tempint['scanid']=scid
        self.adjion_int=pd.concat([self.adjion_int,tempint],ignore_index=True)
        

