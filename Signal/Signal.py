#GlyLine Signal class
#these functions parse the extracted data

import scipy.signal as sig
import numpy as np
import pandas as pd
import tables as tb
import ms_deisotope.averagine as msavg
from GlyLine.Helper import Helper
from GlyLine.Trawler import Trawler
from GlyLine.Meta import Meta

# for getting which acq bins the GPs from MS1Adducted will fall into
class MZWindowFinder:
    def __init__(self,chargerange=[2,6],windowbounds=[800,1600],nbins=50):
        self.chargerange=chargerange
        self.bounds=windowbounds
        self.nbins=nbins
        
    def main(self,NeutralMassDF,AcqWindows=None):
        if AcqWindows is None:
            AcqWindows=self.TheoreticalAcqWindowDFMaker()
            #add in the non acquiring bounds for weirdness
            if AcqWindows.shape[0]==self.nbins:
                AcqWindows.loc[AcqWindows.index.max()+1]=[AcqWindows.loc[AcqWindows.index.max(),'upper'],AcqWindows.loc[AcqWindows.index.max(),'upper']*3]
                AcqWindows.loc[AcqWindows.index.max()+1]=[0,AcqWindows.loc[AcqWindows.index.min(),'lower']]
        self.AcqWindows=AcqWindows
        ChargeDF=self.ChargeDFMaker(NeutralMassDF)
        self.MZDF=self.MZDFMaker(ChargeDF)
    
    def TheoreticalAcqWindowDFMaker(self):
        AcqWindows=pd.DataFrame(None,columns=['lower','upper'])
        stepsize=(self.bounds[1]-self.bounds[0])/self.nbins
        AcqWindows['lower']=np.arange(self.bounds[0],self.bounds[1],step=stepsize)
        AcqWindows['upper']=np.arange(self.bounds[0]+stepsize,self.bounds[1]+stepsize,step=stepsize)
        AcqWindows.loc[AcqWindows.index.max()+1]=[0,self.bounds[0]]
        AcqWindows.loc[AcqWindows.index.max()+1]=[self.bounds[1],self.bounds[1]*2]
        return AcqWindows
        
    def ChargeDFMaker(self,NeutralMassDF):
        ChargeDF=NeutralMassDF.copy()
        for z in range(self.chargerange[0],self.chargerange[1]):
            ChargeDF[z]=Helper.MZfromNM(ChargeDF['neutralmass'], z)
        ChargeDF=ChargeDF.reset_index()
        return ChargeDF
        
    def MZDFMaker(self,ChargeDF):
        MZDF=pd.DataFrame([[0,np.nan,0]],columns=['AddID','mz','charge'])
        for j in range(self.chargerange[0],self.chargerange[1]):
            tempdf=pd.DataFrame(None,columns=['AddID','GPID','mz','charge'])
            tempdf['AddID']=ChargeDF['AddID'].tolist()
            tempdf['GPID']=ChargeDF['GPID'].tolist()
            tempdf['mz']=ChargeDF[j].tolist()
            tempdf['charge']=j
            MZDF=pd.concat([MZDF,tempdf])
        MZDF=MZDF.dropna()
        MZDF['AcqIdx']=[Helper.BoundIndex(amz,self.AcqWindows).tolist()[0] for amz in MZDF['mz']]
        MZDF=MZDF.reset_index().drop('index',axis=1)
        isoclusts=msavg.AveragineCache(msavg.glycopeptide)
        tempdf=pd.DataFrame([[0,np.nan,0,0,0]],columns=['AddID','mz','charge','GPID','AcqIdx'])
        for g in MZDF.index:
            tempiso=isoclusts.isotopic_cluster(MZDF.loc[g,'mz'],MZDF.loc[g,'charge'],truncate_after=0.99)
            tempmzs=[pk.mz for pk in tempiso.peaklist]
            if any(self.AcqWindows.loc[MZDF.loc[g,'AcqIdx'],'lower']>tempmzs):
                tempdf.loc[tempdf.index.max()+1]=MZDF.loc[g]
                tempdf.loc[tempdf.index.max(),'AcqIdx']=MZDF.loc[g,'AcqIdx']-1
            if any(self.AcqWindows.loc[MZDF.loc[g,'AcqIdx'],'upper']<tempmzs):
                tempdf.loc[tempdf.index.max()+1]=MZDF.loc[g]
                tempdf.loc[tempdf.index.max(),'AcqIdx']=MZDF.loc[g,'AcqIdx']+1
        tempdf=tempdf.dropna()
        MZDF=pd.concat([MZDF,tempdf])
        MZDF=MZDF.reset_index().drop('index',axis=1)
        return MZDF
    
# for getting the confirmed indices a product ion will be considered at
class ConfirmedIndexes:
    def __init__(self,PSMMetaMaster,IndexData,GPIonAssoc,MZDF,RunID=None,margin=18):
        self.RunID=RunID
        self.margin=margin
        self.ConfirmedGPIdxMaker(PSMMetaMaster,IndexData,GPIonAssoc,MZDF)
        self.ConfirmedGPInAllWindows(MZDF)
        
    def ConfirmedGPIdxMaker(self,PSMMetaMaster,IndexData,GPIonAssoc,MZDF):
        if self.RunID is None:
            confirmed_gps=PSMMetaMaster.loc[PSMMetaMaster['GPID'].isin(GPIonAssoc.index),['AddID','GPID','scan_time','scan_id','charge']].copy() #add AddID once updated
            timebounds=pd.DataFrame(None,columns=['lower','upper'])
            timebounds['lower']=IndexData.MS1Data['scan_time'].tolist()
            timebounds['upper']=IndexData.MS1Data['scan_time'].tolist()[1:]+[IndexData.MS1Data['scan_time'].max()+IndexData.MS1Data['scan_time'].min()]
            timebounds.index=IndexData.MS1Data.index
            tidx=[]
            aidx=[]
            for j in confirmed_gps.index:
                tidx+=[Helper.BoundIndex(confirmed_gps.loc[j,'scan_time'],timebounds)[0]]
                aidx+=[MZDF.loc[(MZDF['GPID']==confirmed_gps.loc[j,'GPID'])&(MZDF['charge']==confirmed_gps.loc[j,'charge']),'AcqIdx'].iloc[0]]
        else:
            confirmed_gps=PSMMetaMaster.loc[(PSMMetaMaster['RunID']==self.RunID)&(PSMMetaMaster['GPID'].isin(GPIonAssoc.index)),['AddID','GPID','scan_time','scan_id']]
            tidx=[]#^add AddID once updated
            aidx=[]
            for j in confirmed_gps['scan_id']:
                tidx+=[IndexData.MS2Data.loc[IndexData.MS2Data['scan_id']==j.split('.')[0],'PrecursorIdx'].iloc[0]]
                aidx+=[IndexData.MS2Data.loc[IndexData.MS2Data['scan_id']==j.split('.')[0],'AcqIdx'].iloc[0]]
        confirmed_gps['PrecursorIdx']=tidx
        confirmed_gps['AcqIdx']=aidx
        confirmed_times=pd.DataFrame([[np.nan,0,0,0]],columns=['AddID','GPID','PrecursorIdx','AcqIdx'])
        for cgidx in confirmed_gps.index:
            tempdf=pd.DataFrame(None,columns=['AddID','GPID','PrecursorIdx'])
            t=confirmed_gps.loc[cgidx,'PrecursorIdx']
            tempdf['PrecursorIdx']=np.arange(t-self.margin,t+self.margin)
            tempdf['AddID']=confirmed_gps.loc[cgidx,'AddID']
            tempdf['GPID']=confirmed_gps.loc[cgidx,'GPID']
            tempdf['AcqIdx']=confirmed_gps.loc[cgidx,'AcqIdx']
            confirmed_times=pd.concat([confirmed_times,tempdf])
        self.confirmed_times=confirmed_times.dropna().reset_index().drop('index',axis=1)
        self.confirmed_gps=confirmed_gps
    
    def ConfirmedGPInAllWindows(self,MZDF):
        submz=MZDF.loc[MZDF['AddID'].isin(self.confirmed_gps['AddID'].unique())] 
        self.allwindowadds=submz.groupby(['AcqIdx'])['AddID'].apply(lambda x: list(set(x)))
        self.confwindowadds=self.confirmed_gps.groupby(['AcqIdx'])['AddID'].apply(lambda x: list(set(x)))
        
    def AcqIdxInit(self,AcqIdx,GPIonAssoc,AddGPAssoc):
        #now use this to get confirmed times of gpid in a window
        if AcqIdx in self.confwindowadds.index:
            self.CGPinAcq=self.confwindowadds.loc[AcqIdx]
            self.CGPfromIoninAcqDF=pd.DataFrame([[0,0,np.nan]],columns=['IonID','GPID','AddID'])
            for addx in self.CGPinAcq:
                gpdx=AddGPAssoc.loc[addx,'GPID']
                tempdf=pd.DataFrame(None)
                tempdf['IonID']=GPIonAssoc.loc[gpdx,'IonID']
                tempdf['GPID']=gpdx
                tempdf['AddID']=addx
                self.CGPfromIoninAcqDF=pd.concat([ self.CGPfromIoninAcqDF,tempdf])
            self.CGPfromIoninAcqDF=self.CGPfromIoninAcqDF.dropna().reset_index().drop('index',axis=1)
        else:
            self.CGPinAcq=[]
            self.CGPfromIoninAcqDF=None
        self.UGPinAcq=self.allwindowadds.loc[AcqIdx]
        self.UGPfromIoninAcqDF=pd.DataFrame([[0,0,np.nan]],columns=['IonID','GPID','AddID'])        
        for addx in self.UGPinAcq:
            gpdx=AddGPAssoc.loc[addx,'GPID']
            tempdf=pd.DataFrame(None)
            tempdf['IonID']=GPIonAssoc.loc[gpdx,'IonID']
            tempdf['GPID']=gpdx
            tempdf['AddID']=addx
            self.UGPfromIoninAcqDF=pd.concat([ self.UGPfromIoninAcqDF,tempdf])
        self.UGPfromIoninAcqDF=self.UGPfromIoninAcqDF.dropna().reset_index().drop('index',axis=1)
        self.IonGroupsofUGP=self.UGPfromIoninAcqDF.groupby('IonID')['GPID'].apply(lambda x: list(set(x)))
        if AcqIdx in self.confwindowadds.index:
            self.IonGroupsofCGP=self.IonGroupsofUGP.loc[self.IonGroupsofUGP.index.isin(self.CGPfromIoninAcqDF['IonID'])]
        else:
            self.IonGroupsofCGP=None
        
    def GetFinalConfirmedIndexForIonPrecursor(self,IonID):
        #get valid times for the product to be observed
        validgptimes=self.confirmed_times.loc[self.confirmed_times['GPID'].isin(self.IonGroupsofCGP.loc[IonID])]
        validIontimes=np.sort(validgptimes['PrecursorIdx'].unique()).tolist()
        return validIontimes,validgptimes

class ProductIonSignal:
    def __init__(self,MS2ObsDF,IonID,ConfGPIdxObj,AcqIdx,scalebounds=[1,6],referencebounds=[0,1200],
                 imputetype=None,basesignalrat=.5,linearspan=3,peak_dist_thresh=None,lmscorethreshold=.7,InitWidthMin=4,
                 signalminrat=.05,ObsIdxOpt='Confirmed_Strict',confirmed_gp_margin=3,forcelogint=True,intensityname='LR_Intensity'):
        self.scalemin=scalebounds[0]
        self.scalemax=scalebounds[1]
        self.RefMin=referencebounds[0]
        self.RefMax=referencebounds[1]
        self.imputetype=imputetype
        self.linearspan=linearspan
        self.forcelogint=forcelogint
        self.IonID=IonID
        self.AcqIdx=AcqIdx
        self.DF=MS2ObsDF.loc[MS2ObsDF['IonID']==self.IonID].copy()
        if ('LR_Intensity' not in self.DF.columns)&(intensityname=='LR_Intensity'):
            Helper.LogRelativizeIntensity(self.DF)
        self.intensityname=intensityname
        self.basesignal=self.DF[intensityname].min()*basesignalrat
        self.signalmin=self.DF[intensityname].max()*signalminrat
        self.ConfGPIdxObj=ConfGPIdxObj
        self.ObsIdxOpt=ObsIdxOpt
        self.margin=confirmed_gp_margin*self.scalemax
        self.lmscorethreshold=lmscorethreshold
        self.InitWidthMin=InitWidthMin
        if peak_dist_thresh==None:
            self.exceptionthreshmin=self.scalemax
        else:
            self.exceptionthreshmin=peak_dist_thresh
    
    def Main(self,IonGPAssoc,MS2IdxData):
        #get the confirmed observed indexes and subset DF to them
        self.ObservedIndexes(IonGPAssoc)
        if self.DFsub.shape[0]>0:
            self.SingletHandling()
            #if there is all data singlets, skip all but parameterization, otherwise proceed
            if self.DFsub['singlet'].all():
                allsingles=True
                self.InitialCurveParameterization(allsingles)
                outpks=self.potpks
                outpks['Coef']=1
                outpks['score']=1
                outpks['ApexWeight']=1
                pkref=pd.DataFrame([['str',0,np.nan,0]],columns=['pid','IonID','AddID','ApexDist'])
                d_scantime=MS2IdxData.loc[(MS2IdxData['AcqIdx']==self.AcqIdx)&(MS2IdxData['PrecursorIdx'].isin([1,2])),'scan_time'].diff().tolist()[1]
                for idx in outpks.index:
                    outpks.loc[idx,'AUC']=list(self.DFsub.loc[self.DFsub['PrecursorIdx']==self.potpks.loc[idx,'apex'],'Intensity']*d_scantime/2)[0]
                    hits=self.ConfGPIdxObj.confirmed_gps.loc[(self.ConfGPIdxObj.confirmed_gps['PrecursorIdx']==outpks.loc[idx,'apex'])&(self.ConfGPIdxObj.confirmed_gps['AcqIdx']==self.AcqIdx),'AddID'].tolist()
                    if len(hits)>0:
                        for h in hits:
                            pkref.loc[pkref.index.max()+1]=[outpks.loc[idx,'pid'],self.IonID,h,0]
            else:
                imputedf=Helper.BasicImputation(self.DFsub.loc[self.DFsub['singlet']==False],RefMin=self.RefMin,RefMax=self.RefMax,
                                                imputetype=self.imputetype,basesignal=self.basesignal,
                                                linearspan=self.linearspan,dataname=self.intensityname)
                self.CWTMaker(imputedf[self.intensityname])
                allsingles=False
                self.InitialCurveParameterization(allsingles)
                self.InitialCurveGeneration()
                emptyids,curvesout=self.EmptyIndexFinder(self.exceptionthreshmin)
                self.ContinueCheck=True
                self.ignorepid=[]
                whilecount=0
                wcmax=10
                while (len(emptyids)>0) & (self.ContinueCheck) & (whilecount<wcmax):
                    self.EmptyIndexAdjuster(emptyids,curvesout)
                    #repeat until empty values gone 
                    emptyids,curvesout=self.EmptyIndexFinder(self.exceptionthreshmin)
                    whilecount+=1
                #if there are still missing values after 10 iterations
                #make a minimum scale wavelet at each remaining emptyid
                if (whilecount==wcmax) & (len(emptyids)>0):
                    self.EmergencyAdjust(emptyids)
                lm,X,y=Helper.LinearModelGenerator(self.curvedf, self.potpks.loc[self.potpks['InModel'],'pid'].tolist())
                while any(lm.coef_<=self.signalmin):
                    self.RemovePeaksBelowThresh(lm.coef_)
                    lm,X,y=Helper.LinearModelGenerator(self.curvedf, self.potpks.loc[self.potpks['InModel'],'pid'].tolist())
                currentscore=lm.score(X,y)
                if currentscore<self.lmscorethreshold:
                    currentscore=self.ScoreBelowThresholdProcess()
                lm,X,y=Helper.LinearModelGenerator(self.curvedf, self.potpks.loc[self.potpks['InModel'],'pid'].tolist())
                while any(lm.coef_<=self.signalmin):
                    self.RemovePeaksBelowThresh(lm.coef_)
                    lm,X,y=Helper.LinearModelGenerator(self.curvedf, self.potpks.loc[self.potpks['InModel'],'pid'].tolist())
                outpks,pkref=self.CheckModelEndProcessing(currentscore,IonGPAssoc,MS2IdxData)
        else:
            outpks=None
            pkref=None
        return outpks,pkref
    
    def ObservedIndexes(self,IonGPAssoc):
        if self.ObsIdxOpt=='Observed':
            self.ObsIdx=self.DF['PrecursorIdx'].tolist()
        elif self.ObsIdxOpt=='Confirmed':
            self.ObsIdx=self.ConfGPIdxObj.GetFinalConfirmedIndexForIonPrecursor(self.IonID)[0]
        elif self.ObsIdxOpt=='Confirmed_Strict':
            cgp=list(set(IonGPAssoc.loc[self.IonID,'GPID'])&set(self.ConfGPIdxObj.confirmed_gps.loc[self.ConfGPIdxObj.confirmed_gps['AcqIdx']==self.AcqIdx,'GPID'])) 
            self.ObsIdx=self.ConfGPIdxObj.confirmed_times.loc[(self.ConfGPIdxObj.confirmed_times['AcqIdx']==self.AcqIdx)&(self.ConfGPIdxObj.confirmed_times['GPID'].isin(cgp)),'PrecursorIdx'].tolist()
        else:
            self.ObsIdx=np.arange(self.RefMin,self.RefMax+1)
        obsidxes=list(set(self.DF['PrecursorIdx'])&set(self.ObsIdx))
        self.conf_obsidxes=list(np.sort(obsidxes))
        self.DFsub=self.DF.loc[self.DF['PrecursorIdx'].isin(self.conf_obsidxes)].copy()
        
    def SingletHandling(self):
        preidxdiff=self.DFsub['PrecursorIdx'].diff().tolist()[1:]
        self.DFsub['priordiff']=[self.scalemax+1]+preidxdiff
        self.DFsub['postdiff']=preidxdiff+[self.scalemax+1]
        self.DFsub['singlet']=(self.DFsub['priordiff']>self.scalemax)&(self.DFsub['postdiff']>self.scalemax)
        self.conf_obsidxes=list(set(self.conf_obsidxes)-set(self.DFsub.loc[self.DFsub['singlet'],'PrecursorIdx'].tolist()))
        
    def CWTMaker(self,signalvector):
        self.cwtmatr=sig.cwt(signalvector,sig.ricker,np.arange(self.scalemin,self.scalemax))
        
    def InitialCurveParameterization(self,AllSingles,RemoveScale1=False):
        potpks=pd.DataFrame([[0,0,np.nan,0,0,'NA','NA','NA',0,False,False,0.0]],columns=['apex','scale','wavelen','startscan','endscan','pid','originpid','LastTransform','LTDirection','InSolution','InModel','Coef'])
        self.pidcounter=1
        for singidx in self.DFsub.loc[self.DFsub['singlet']==True].index.tolist():
            pid='pid_'+str(self.IonID)+'_'+str(self.AcqIdx)+'_'+str(self.pidcounter)
            j=self.DFsub.loc[singidx,'PrecursorIdx']
            pot=[j,1,3,j-1,j+1+1,pid,pid,'start',0,True,False,1]
            potpks.loc[len(potpks.index)+1]=pot
            self.pidcounter+=1
        if AllSingles==False:
            scopes=pd.DataFrame(np.array([[0.0]*(self.scalemax+1-(self.scalemin-1))]*self.cwtmatr.shape[1]),
                                columns=list(np.arange(self.scalemin-1,self.scalemax+1)),index=np.arange(self.RefMin,self.RefMax))
            keep=[]
            for j in range(self.cwtmatr.shape[0]):
                hits=list(sig.argrelextrema(self.cwtmatr[j,:], np.greater)[0])
                keep=list(set(keep)|set(hits))
                scopes[j+self.scalemin].iloc[hits]=self.cwtmatr[j,hits].tolist()
            keep.sort()
            scopes=scopes.iloc[keep]
            for j in scopes.index:
                hits=sig.argrelextrema(np.array(scopes.loc[j].tolist()), np.greater)[0].tolist()
                if len(hits)>0:
                    for h in hits:
                        if scopes.loc[j,h]>=self.signalmin:
                            scale=scopes.columns[h]
                            inobs=set(list(range(j-scale,j+scale+2)))&set(self.ObsIdx)
                            if len(inobs)>0:
                                pid='pid_'+str(self.IonID)+'_'+str(self.AcqIdx)+'_'+str(self.pidcounter)
                                pot=[j,scale,scale*2+1,j-scale,j+scale+1,pid,pid,'start',0,False,True,0]
                                potpks.loc[len(potpks.index)+1]=pot
                                self.pidcounter+=1
        self.potpks=potpks.dropna()
        if RemoveScale1==True:
            if any(potpks['scale']==1):
                self.ScaleMinRemover()
    
    def ScaleMinRemover(self):
        TargetObj=Trawler.IonTargetList()
        targetdf=pd.DataFrame(None,columns=['upper','lower'])
        targetdf['upper']=self.potpks.loc[self.potpks['scale']!=1,'endscan'].tolist()
        targetdf['lower']=self.potpks.loc[self.potpks['scale']!=1,'startscan'].tolist()
        for j in self.potpks.loc[self.potpks['scale']==1].index:
            if TargetObj.BoundBool(self.potpks.loc[j,'apex'],targetdf):
                self.potpks.loc[j,'InModel']=False
    
    def InitialCurveGeneration(self,tgtpids=None):
        curvedf=pd.DataFrame(None,columns=['PrecursorIdx','Intensity'])
        curvedf['PrecursorIdx']=list(np.arange(self.RefMin-1,self.RefMax+1))
        curvedf['Intensity']=0.0
        curvedf=curvedf.set_index(['PrecursorIdx'])
        curvedf.loc[self.conf_obsidxes,'Intensity']=self.basesignal
        curvedf.loc[self.DFsub.loc[self.DFsub['singlet']==False,'PrecursorIdx'].tolist(),'Intensity']=self.DFsub.loc[self.DFsub['singlet']==False,self.intensityname].tolist()
        self.curvedf=curvedf
        self.CurveAdder()
        
    def CurveAdder(self,colbool='InModel'):
        newcurves=list(set(self.potpks.loc[self.potpks[colbool]==True,'pid'])-set(self.curvedf.columns))
        if len(newcurves)>0:
            idxlist=self.potpks.loc[self.potpks['pid'].isin(newcurves)].index.tolist()
            newpiddf=pd.DataFrame(0.0,columns=self.potpks.loc[idxlist,'pid'].tolist(),index=self.curvedf.index)
            for rix in idxlist:
                pid=self.potpks.loc[rix,'pid']
                tempidx=list(np.arange(self.potpks.loc[rix,'apex']-self.potpks.loc[rix,'scale'],
                                       self.potpks.loc[rix,'apex']+self.potpks.loc[rix,'wavelen']-self.potpks.loc[rix,'scale']))
                newpiddf.loc[tempidx,pid]=sig.ricker(self.potpks.loc[rix,'wavelen'],self.potpks.loc[rix,'scale'])
            self.curvedf=pd.concat([self.curvedf,newpiddf],axis=1)
    
    def EmptyIndexFinder(self,apexdist):
        curvesout=pd.DataFrame([[0,np.nan,'str']],columns=['PrecursorIdx','PeakIntensity','pid'])
        pidtgts=self.potpks.loc[self.potpks['InModel']==True,'pid'].tolist()
        for pid in pidtgts:
            tempcurves=pd.DataFrame(None,columns=['PrecursorIdx','PeakIntensity','pid'])
            tempcurves['PrecursorIdx']=np.arange(self.potpks.loc[self.potpks['pid']==pid,'startscan'].iloc[0],
                                                 self.potpks.loc[self.potpks['pid']==pid,'endscan'].iloc[0])
            tempcurves['pid']=pid
            tempcurves['PeakIntensity']=self.curvedf.loc[tempcurves['PrecursorIdx'],pid].tolist()
            curvesout=pd.concat([curvesout,tempcurves])
        curvesout=curvesout.dropna()
        d={'PeakIntensity':'SummedIntensity','pid':'pids'}
        grpedgp=curvesout.groupby(['PrecursorIdx']).agg({'PeakIntensity':'sum','pid':lambda x: list(x)}).rename(columns=d)
        grpedgp=grpedgp.loc[grpedgp['SummedIntensity']>0]
        #find the empty intensity PrecursorIdx, first get the ones that are totally missing
        emptyids=list(set(self.conf_obsidxes)-set(grpedgp.index.tolist()))
        #find the ones without any peakids assoc, find closest scans with peaks
        #if within exception distance, include, if not generate new peak
        if len(emptyids)>0:
            for j in emptyids:
                emptydist=abs(j-self.potpks.loc[self.potpks['InModel']==True,'apex'])
                if np.min(emptydist)<=apexdist:
                    mindistilocs=np.where(emptydist==emptydist.min())[0].tolist()
                    grpedgp=pd.concat([grpedgp,pd.DataFrame([[0,self.potpks.loc[self.potpks['InModel']==True,'pid'].iloc[mindistilocs].tolist()]],columns=['SummedIntensity','pids'],index=[j])])
                    for k in grpedgp['pids'].loc[j]:
                        curvesout=pd.concat([curvesout,pd.DataFrame([[j,0,k]],columns=['PrecursorIdx','PeakIntensity','pid'])])
                else:
                    #add all the emptyids outside of exception range to potpks
                    if any((self.potpks['apex']==j) & (self.potpks['scale']==1) & (self.potpks['wavelen']==3)):
                        self.potpks.loc[(self.potpks['apex']==j) & (self.potpks['scale']==1) & (self.potpks['wavelen']==3),'InModel']=True
                    else:
                        pid='pid_'+str(self.IonID)+'_'+str(self.AcqIdx)+'_'+str(self.pidcounter)
                        self.potpks.loc[self.potpks.index.max()+1]=[j,1,3,j-1,j+1,pid,pid,'start',0,False,True,0]
                        self.pidcounter+=1
                    self.CurveAdder()
                    grpedgp=pd.concat([grpedgp,pd.DataFrame([[0.9,[self.potpks.loc[(self.potpks['apex']==j) & (self.potpks['scale']==1) & (self.potpks['wavelen']==3),'pid']]]],columns=['SummedIntensity','pids'],index=[j])])
                    curvesout=pd.concat([curvesout,pd.DataFrame([[j,0.9,self.potpks.loc[(self.potpks['apex']==j) & (self.potpks['scale']==1) & (self.potpks['wavelen']==3),'pid'].iloc[0]]],columns=['PrecursorIdx','PeakIntensity','pid'])])
        #now take the ones that are now seen but have a total intensity of 0
        emptyout=grpedgp.loc[self.conf_obsidxes].loc[(grpedgp['SummedIntensity']==0)].index.tolist()
        return emptyout, curvesout
    
    def EmptyIndexAdjuster(self,emptyids,curvesout):
        pids=curvesout.loc[curvesout['PrecursorIdx'].isin(emptyids),'pid']
        #look at peak with most missing associations
        if len(self.ignorepid)>0:
            pids=pids.loc[~pids.isin(self.ignorepid)]
        if pids.shape[0]==0:
            self.ContinueCheck=False
        else:
            self.ContinueCheck=True
            pidtgt=pids.mode()
            for pt in pidtgt:
                ptidx=self.potpks.loc[self.potpks['pid']==pt].index
                pkdata=self.potpks.loc[ptidx]
                #if both pk ends are missing then increase scale
                #if one is missing shift the scale
                #if the peak is beyond range increase scale
                skip=False #skip adding peak if maxscale is reached
                if ((pkdata['startscan'].iloc[0] in emptyids) | (pkdata['endscan'].iloc[0] in emptyids)):
                    #all cases where emptiness is in tails
                    if ((pkdata['startscan'].iloc[0] in emptyids) & (pkdata['endscan'].iloc[0] in emptyids)):
                        #where emptiness in both tails
                        pkmaxloc=pkdata['scale'].iloc[0]
                        if (pkmaxloc<self.scalemax):
                            scale=pkmaxloc+1
                            apex=pkdata['apex'].iloc[0]
                            lowerbound=apex-scale
                            upperbound=apex+scale+1
                            wavelen=scale*2+1
                            transform='scale'
                            ltdirection=1
                            newpid=pt+'_scUP'
                        else:
                            #where scale is maxed out
                            self.ignorepid+=[pt]
                            skip=True
                    else:
                        #where emptiness is single tailed
                        pkmaxloc=pkdata['scale'].iloc[0]-1
                        scale=pkdata['scale'].iloc[0]
                        wavelen=scale*2
                        transform='shoulder'
                        if (pkdata['startscan'].iloc[0] in emptyids):
                            #lower tailed
                            apex=pkdata['apex'].iloc[0]
                            newpid=pt+'_sh0'
                            ltdirection=0
                        else:
                            #upper tailed
                            apex=pkdata['apex'].iloc[0]+1
                            newpid=pt+'_sh1'
                            ltdirection=1
                        lowerbound=apex-scale
                        upperbound=apex+scale
                else:
                    #all cases where emptiness is beyond tails
                    missingpres=curvesout.loc[(curvesout['pid']==pt)&(curvesout['PrecursorIdx'].isin(emptyids)),'PrecursorIdx'].tolist()
                    if (pkdata['LastTransform'].iloc[0]=='scale'):
                        pkmaxloc=pkdata['scale'].iloc[0]-1
                        scale=pkdata['scale'].iloc[0]
                        wavelen=scale*2
                        transform='shoulder'
                        if any(missingpres<pkdata['startscan'].iloc[0]):
                            #lower tailed
                            apex=pkdata['apex'].iloc[0]
                            newpid=pt+'_sh0'
                            ltdirection=0
                        else:
                            #upper tailed
                            apex=pkdata['apex'].iloc[0]+1
                            newpid=pt+'_sh1'
                            ltdirection=1
                        lowerbound=apex-scale
                        upperbound=apex+scale
                    else:
                        #either is beyond both tails or just did a shoulder adjustment
                        pkmaxloc=pkdata['scale'].iloc[0]
                        if (pkmaxloc<self.scalemax):
                            scale=pkmaxloc+1
                            if all(missingpres<pkdata['startscan'].iloc[0]):
                                apex=pkdata['apex'].iloc[0]-1
                                newpid=pt+'_scUP-ap'
                            else:
                                apex=pkdata['apex'].iloc[0]
                                newpid=pt+'_scUP'
                            lowerbound=apex-scale
                            upperbound=apex+scale+1
                            wavelen=scale*2+1
                            transform='scale'
                            ltdirection=1
                        else:
                            #where scale is maxed out
                            self.ignorepid+=[pt]
                            skip=True
                if skip==False:
                    self.potpks.loc[ptidx,'InModel']=False
                    self.curvedf=self.curvedf.drop(pt,axis=1)
                    if any((self.potpks['apex']==apex) & (self.potpks['scale']==scale) & (self.potpks['wavelen']==wavelen)):
                        tempidx=self.potpks.loc[(self.potpks['apex']==apex) & (self.potpks['scale']==scale) & (self.potpks['wavelen']==wavelen)].index
                        self.potpks.loc[tempidx,'InModel']=True
                    else:
                        self.potpks.loc[self.potpks.index.max()+1]=[apex,scale,wavelen,lowerbound,upperbound,newpid,pt,transform,ltdirection,False,True,0]
                    self.CurveAdder()
                    
    def EmergencyAdjust(self,emptyids):
        for e in emptyids:
            newpid='prod_'+str(self.IonID)+'_'+str(self.AcqIdx)+'_'+str(self.pidcounter)
            self.pidcounter+=1
            self.potpks.loc[self.potpks.index.max()+1]=[e,1,3,e-1,e+2,newpid,newpid,'start',0,False,True,0]
        self.CurveAdder()
        
    def RemovePeaksBelowThresh(self,coef,report=False):
        reportl=[]
        for idx,pdx in enumerate(self.potpks.loc[self.potpks['InModel']].index):
            if coef[idx]<=self.signalmin:
                self.potpks.loc[pdx,'InModel']=False
                self.potpks.loc[pdx,'InSolution']=False
                self.curvedf=self.curvedf.drop(self.potpks.loc[pdx,'pid'],axis=1)
                reportl+=[self.potpks.loc[pdx,'pid']]
        if report==True:
            return reportl
    
    def ScoreBelowThresholdProcess(self):
        newscore=self.MultiCurveApexRemoval()
        if newscore<self.lmscorethreshold:
            newscore=self.CheckRemovalofApexMissingPIDs(newscore)
            if newscore<self.lmscorethreshold:
                currentcoef,X,y=Helper.LinearModelGenerator(self.curvedf, self.potpks.loc[self.potpks['InModel'],'pid'].tolist())
                newscore=self.AdjacentCurveRemoval(currentcoef.coef_,newscore)
                # if newscore<self.lmscorethreshold:
                #     lm,X,y=Helper.LinearModelGenerator(self.curvedf, self.potpks.loc[self.potpks['InModel'],'pid'].tolist())
                #     PEDF,RSSDF=self.ResidualFinder(lm,X,y)
        return newscore
    
    def MultiCurveApexRemoval(self):
        currentmodel=self.potpks.loc[self.potpks['InModel']]
        if any(currentmodel.duplicated('apex',keep=False)):
            apexrpts=currentmodel.loc[currentmodel.duplicated('apex'),'apex'].tolist()
            for ap in apexrpts:
                rptpids=currentmodel.loc[currentmodel['apex']==ap].sort_values('scale',ascending=False)
                maxpid=currentmodel.loc[(currentmodel['apex']==ap)&(currentmodel['scale'].max()==currentmodel['scale']),'pid']
                for rdx in rptpids.loc[~rptpids['pid'].isin(maxpid)].index:
                    self.potpks.loc[self.potpks['pid']==rptpids.loc[rdx,'pid'],'InModel']=False
                    self.curvedf=self.curvedf.drop(rptpids.loc[rdx,'pid'],axis=1)
        lmtemp,Xtemp,ytemp=Helper.LinearModelGenerator(self.curvedf, self.potpks.loc[self.potpks['InModel'],'pid'].tolist())
        newscore=lmtemp.score(Xtemp,ytemp)
        return newscore
    
    def CheckRemovalofApexMissingPIDs(self,currentscore):
        relscore=currentscore
        currentmodel=self.potpks.loc[self.potpks['InModel']]
        MissingApexes=list(set(self.potpks.loc[self.potpks['InModel'],'apex'])&set(self.curvedf.loc[self.curvedf['Intensity']<=self.basesignal]))
        if len(MissingApexes)>0:
            for MA in MissingApexes:
                mapids=self.potpks.loc[self.potpks['apex']==MA,'pid']
                jscore=pd.DataFrame(None,columns=['pid','score'])
                jscore['pid']=mapids
                jscore['score']=currentscore
                currentbestmodel=mapids
                for mapid in mapids:
                    pkdata=self.potpks.loc[self.potpks['pid']==mapid].iloc[0]
                    self.potpks.loc[self.potpks['pid']==mapid,'InModel']=False
                    self.curvedf=self.curvedf.drop(mapid,axis=1)
                    if pkdata['scale']>1:
                        scale=pkdata['scale']
                        newpids=[mapid+'_sp0',mapid+'_sp1']
                        if scale%2==0:
                            newscale=scale/2
                            newcurve1=[pkdata['apex']-newscale,newscale,newscale*2,pkdata['startscan'],pkdata['apex']-1,newpids[0],mapid,'split',-1,False,True] #add in 0 if after coef
                            if pkdata['wavelen']%2:
                                newcurve2=[pkdata['apex']+newscale,newscale,newscale*2+1,pkdata['apex'],newscale*2+1+pkdata['apex'],newpids[1],mapid,'split',1,False,True]
                            else:
                                newcurve2=[pkdata['apex']+newscale+1,newscale,newscale*2,pkdata['apex']+1,pkdata['apex']+1+newscale*2,newpids[1],mapid,'split',1,False,True]
                        else:
                            if pkdata['wavelen']%2==1:
                                newwavelen=((pkdata['wavelen']-1)/2)+2
                                newscale=(newwavelen-1)/2
                                newcurve1=[pkdata['apex']-newscale,newscale,newwavelen,pkdata['startscan'],pkdata['apex'],newpids[0],mapid,'split',-1,False,True] #add in 0 if after coef
                                newcurve2=[pkdata['apex']+newscale,newscale,newwavelen,pkdata['apex'],pkdata['endscan'],newpids[1],mapid,'split',1,False,True]
                            else:
                                newwavelen=((pkdata['wavelen'])/2)+2
                                newscale1=(newwavelen-1)/2
                                newscale2=newscale1-1
                                newcurve1=[pkdata['apex']-newscale1,newscale1,newwavelen,pkdata['startscan'],pkdata['apex'],newpids[0],mapid,'split',-1,False,True] #add in 0 if after coef
                                newcurve2=[pkdata['apex']+newscale2+1,newscale2,newscale2*2,pkdata['apex']+1,pkdata['endscan'],newpids[1],mapid,'split',1,False,True]
                        if 'Coef' in self.potpks.columns:
                            newcurve1+=[0]
                            newcurve2+=[0]
                        self.potpks.loc[self.potpks.index.max()+1]=newcurve1
                        self.potpks.loc[self.potpks.index.max()+1]=newcurve2
                    else:
                        newpids=[mapid+'_ap']
                        newcurve1==[pkdata['apex']-1,1,3,pkdata['apex']-2,pkdata['apex'],newpids[0],mapid,'shift',-1,False,True]
                        if 'Coef' in self.potpks.columns:
                            newcurve1+=[0]
                        
                        self.potpks.loc[self.potpks.index.max()+1]=newcurve1
                    self.CurveAdder()
                    lmtemp,Xtemp,ytemp=Helper.LinearModelGenerator(self.curvedf, self.potpks.loc[self.potpks['InModel'],'pid'].tolist())
                    newscore=lmtemp.score(Xtemp,ytemp)
                    if newscore<relscore:
                        self.potpks.loc[self.potpks['pid'].isin(newpids),'InModel']=False
                        self.curvedf=self.curvedf.drop(newpids,axis=1)
                        self.potpks.loc[self.potpks['pid']==mapid,'InModel']=True
                        self.CurveAdder()
                    else:
                        relscore=newscore
        return relscore
    
    def AdjacentCurveRemoval(self,currentcoef,currentscore):
        relscore=currentscore
        AdjDF=pd.DataFrame(None)
        currentmodel=self.potpks.loc[self.potpks['InModel'],'pid']
        AdjDF['pid']=currentmodel.tolist()
        AdjDF['apex']=self.potpks.loc[self.potpks['InModel'],'apex'].tolist()
        AdjDF['scale']=self.potpks.loc[self.potpks['InModel'],'scale'].tolist()
        AdjDF['Coef']=currentcoef.tolist()
        AdjDF=AdjDF.sort_values('apex')
        diffs=abs(AdjDF['apex'].diff()).dropna().values
        AdjDF['ApexDist']=diffs.tolist()+[0]
        idxhold=np.where(diffs<=1)[0].tolist()
        for idx in np.where(diffs<=1)[0].tolist():
            if idx in idxhold:
                sigrat=self.DFsub.loc[self.DFsub['PrecursorIdx'].isin(AdjDF['apex'].iloc[[idx,idx+1]]),'LR_Intensity'].mean()/self.DFsub['LR_Intensity'].max()
                if sigrat<=.05:
                    firstdrop=AdjDF.loc[AdjDF['scale'].iloc[[idx,idx+1]]==AdjDF['scale'].min(),'pid'].tolist()[0]
                    self.potpks.loc[self.potpks['pid']==firstdrop,'InModel']=False
                    self.curvedf=self.curvedf.drop(firstdrop,axis=1)
                    idxhold=list(set(idxhold)-set([idx]))
                else:
                    coefrat=abs(AdjDF['Coef'].iloc[idx]-AdjDF['Coef'].iloc[idx+1])/AdjDF['Coef'].iloc[[idx,idx+1]].max()
                    if ((coefrat<0.05)&(sigrat<.5)):
                        firstdrop=AdjDF.loc[AdjDF['Coef']==AdjDF['Coef'].iloc[[idx,idx+1]].min(),'pid'].tolist()[0]
                        self.potpks.loc[self.potpks['pid']==firstdrop,'InModel']=False
                        self.curvedf=self.curvedf.drop(firstdrop,axis=1)
                        lmtemp,Xtemp,ytemp=Helper.LinearModelGenerator(self.curvedf, self.potpks.loc[self.potpks['InModel'],'pid'].tolist())
                        newscore=lmtemp.score(Xtemp,ytemp)
                        if newscore<relscore:
                            self.potpks.loc[self.potpks['pid']==firstdrop,'InModel']=True
                            self.CurveAdder()
                        else:
                            relscore=newscore
                            idxhold=list(set(idxhold)-set([idx]))
        return relscore
    
    def CheckModelEndProcessing(self,priorlmscore,IonGPAssoc,MS2IdxData):
        emptyids,curvesout=self.EmptyIndexFinder(2)
        #keep model for later reference
        priormodel=self.potpks.loc[self.potpks['InModel'],'pid']
        relscore=priorlmscore
        #get pertinent pid to test adjustments
        cgp=list(set(IonGPAssoc.loc[self.IonID,'GPID'])&set(self.ConfGPIdxObj.confirmed_gps.loc[self.ConfGPIdxObj.confirmed_gps['AcqIdx']==self.AcqIdx,'GPID']))
        cidx=self.ConfGPIdxObj.confirmed_gps.loc[(self.ConfGPIdxObj.confirmed_gps['AcqIdx']==self.AcqIdx)&(self.ConfGPIdxObj.confirmed_gps['GPID'].isin(cgp)),'PrecursorIdx'].tolist()
        relpid=curvesout.loc[curvesout['PrecursorIdx'].isin(cidx),'pid'].unique().tolist()
        fwdscore=0
        bwdscore=0
        #get model optimized from one direction
        if len(relpid)==0:
            self.potpks.loc[self.potpks.index.max()+1]=[0,np.nan,0,0,0,'NA','NA','start',0,False,True,0]
            outpks=self.potpks.loc[self.potpks['InSolution']].copy()
            pkref=pd.DataFrame([['NA',0,np.nan,0]],columns=['pid','IonID','AddID','ApexDist'])
        else:
            if len(relpid)==1:
                fwdscore=self.ValidCurveNeighborCheck(relpid[0],relscore)
                bwdscore=fwdscore
                fwdmodel=self.potpks.loc[self.potpks['InModel'],'pid'].tolist()
                bwdmodel=fwdmodel
            else:
                for r in relpid:
                    fwdscore=self.ValidCurveNeighborCheck(r,relscore)
                fwdmodel=self.potpks.loc[self.potpks['InModel'],'pid'].tolist()
                #reset parameters
                relscore=priorlmscore
                self.curvedf=self.curvedf.drop(self.curvedf.columns.tolist()[1:],axis=1)
                self.potpks['InModel']=False
                self.potpks.loc[self.potpks['pid'].isin(priormodel),'InModel']=True
                self.CurveAdder()
                #get model from opposite direction
                relpid.reverse()
                for r in relpid:
                    bwdscore=self.ValidCurveNeighborCheck(r,relscore)
                bwdmodel=self.potpks.loc[self.potpks['InModel'],'pid'].tolist()
            #find best model
            maxscore=np.max([priorlmscore,fwdscore,bwdscore])
            #reset parameters
            self.potpks['InModel']=False
            self.curvedf=self.curvedf.drop(self.curvedf.columns.tolist()[1:],axis=1)
            if maxscore>0:
                if (fwdscore==maxscore):
                    self.potpks.loc[self.potpks['pid'].isin(fwdmodel),'InSolution']=True
                    self.potpks.loc[self.potpks['pid'].isin(fwdmodel),'InModel']=True
                elif (bwdscore==maxscore):
                    self.potpks.loc[self.potpks['pid'].isin(bwdmodel),'InSolution']=True
                    self.potpks.loc[self.potpks['pid'].isin(bwdmodel),'InModel']=True
                else:
                    self.potpks.loc[self.potpks['pid'].isin(priormodel),'InSolution']=True
                    self.potpks.loc[self.potpks['pid'].isin(priormodel),'InModel']=True
                self.CurveAdder()
                #get best model
                lm,X,y=Helper.LinearModelGenerator(self.curvedf, self.potpks.loc[self.potpks['InModel'],'pid'].tolist())
                while any(lm.coef_<=self.signalmin):
                    self.RemovePeaksBelowThresh(lm.coef_)
                    lm,X,y=Helper.LinearModelGenerator(self.curvedf, self.potpks.loc[self.potpks['InModel'],'pid'].tolist())
                #record coefficients
                coefdict={pk:lm.coef_[ix] for ix,pk in enumerate(self.potpks.loc[self.potpks['InModel'],'pid'].tolist())}
                coefdict.update({'Intensity':1})
                self.potpks['Coef']=0.0
                for col in self.curvedf.columns:
                    self.potpks.loc[self.potpks['pid']==col,'Coef']=coefdict[col]
                for ap in self.DFsub.loc[self.DFsub['singlet'],'PrecursorIdx']:
                    self.potpks.loc[(self.potpks['apex']==ap)&(self.potpks['scale']==1),'Coef']=1
                #produce output of final pid solution for ion in window
                outpks=self.potpks.loc[self.potpks['InSolution']].copy()
                outpks['score']=maxscore
                #produce output of pkreference for AddID correlation
                pkref=pd.DataFrame([['NA',0,np.nan,0]],columns=['pid','IonID','AddID','ApexDist'])
                #get the confirmed AddID in the window
                gpOfj_ionid=self.ConfGPIdxObj.CGPfromIoninAcqDF.loc[self.ConfGPIdxObj.CGPfromIoninAcqDF['IonID']==self.IonID,'AddID'].tolist()
                for odx in outpks.index:
                    #get the addid and related precursor for each pid
                    if outpks.loc[odx,'wavelen']%2==0:
                        oadds=self.ConfGPIdxObj.confirmed_gps.loc[(outpks.loc[odx,'endscan']>=self.ConfGPIdxObj.confirmed_gps['PrecursorIdx'])&(self.ConfGPIdxObj.confirmed_gps['PrecursorIdx']>=outpks.loc[odx,'startscan'])&(self.ConfGPIdxObj.confirmed_gps['AcqIdx']==self.AcqIdx)&(self.ConfGPIdxObj.confirmed_gps['AddID'].isin(gpOfj_ionid)),['AddID','PrecursorIdx']]
                    else:
                        oadds=self.ConfGPIdxObj.confirmed_gps.loc[((outpks.loc[odx,'endscan']-1)>self.ConfGPIdxObj.confirmed_gps['PrecursorIdx'])&(self.ConfGPIdxObj.confirmed_gps['PrecursorIdx']>outpks.loc[odx,'startscan'])&(self.ConfGPIdxObj.confirmed_gps['AcqIdx']==self.AcqIdx)&(self.ConfGPIdxObj.confirmed_gps['AddID'].isin(gpOfj_ionid)),['AddID','PrecursorIdx']]
                    if len(oadds)>0:
                        for o in oadds.index:
                            #get distance to apex of cgp
                            pkref.loc[pkref.index.max()+1]=[outpks.loc[odx,'pid'],self.IonID,oadds.loc[o,'AddID'],outpks.loc[odx,'apex']-oadds.loc[o,'PrecursorIdx']]
                    else:
                        #if none record 0 event
                        pkref.loc[pkref.index.max()+1]=[outpks.loc[odx,'pid'],0,0,0]
                pkref=pkref.dropna().reset_index(drop=True)
                self.CurveAdder(colbool='InSolution')
                self.AUC_WeightedbyPID(MS2IdxData,outpks)
            else:
                outpks=None
                pkref=None
        return outpks,pkref
    
    def AUC_WeightedbyPID(self,MS2IdxData,outpks,aucmethod='exp_wt_log'):
        outpks['AUC']=[0.0]*outpks.shape[0]
        if aucmethod=='both':
            outpks['AUCalt']=[0.0]*outpks.shape[0]
        outpks['ApexWeight']=[1.0]*outpks.shape[0]
        if hasattr(self,'curvedf'):
            WeightDF=self.curvedf.copy()
            WeightDF=WeightDF.drop(WeightDF.index[~WeightDF.index.isin(MS2IdxData['PrecursorIdx'])].tolist())
            WeightDF=WeightDF.loc[WeightDF['Intensity']>self.basesignal].drop('Intensity',axis=1)
            for col in outpks.loc[outpks['InSolution'],'pid']:
                WeightDF[col]=np.exp(WeightDF[col]*outpks.loc[outpks['pid']==col,'Coef'].iloc[0])-1
            rowtotals=WeightDF.sum(axis=1)
            WeightDF=WeightDF.loc[WeightDF.sum(axis=1)>0]
            ExpIntDF=pd.DataFrame(0.0,columns=['ExpInt'],index=WeightDF.index)
            ExpIntDF.loc[self.DFsub.loc[self.DFsub['PrecursorIdx'].isin(WeightDF.index),'PrecursorIdx'].tolist(),'ExpInt']=self.DFsub.loc[self.DFsub['PrecursorIdx'].isin(WeightDF.index),'Intensity'].tolist()
            # if self.forcelogint==True:
            #     LogIntDF.loc[self.DFsub.loc[self.DFsub['PrecursorIdx'].isin(WeightDF.index),'PrecursorIdx'].tolist(),'LogInt']=np.log(self.DFsub.loc[self.DFsub['PrecursorIdx'].isin(WeightDF.index),'Intensity']).tolist()
            # else:
            #     LogIntDF.loc[self.DFsub.loc[self.DFsub['PrecursorIdx'].isin(WeightDF.index),'PrecursorIdx'].tolist(),'LogInt']=self.DFsub.loc[self.DFsub['PrecursorIdx'].isin(WeightDF.index),'Intensity'].tolist()
            scan_times=[MS2IdxData.loc[(MS2IdxData['AcqIdx']==self.AcqIdx)&(MS2IdxData['PrecursorIdx']==wdix),'scan_time'].iloc[0] for wdix in WeightDF.index]
            WeightDF['scan_time']=scan_times
            for col in outpks.loc[outpks['InSolution'],'pid']:
                WeightDF[col]=WeightDF[col]/rowtotals
                apextemp=outpks.loc[outpks['pid']==col,'apex'].iloc[0]
                SubWeight=WeightDF.loc[WeightDF[col]>0]
                if apextemp in WeightDF.index:
                    outpks.loc[outpks['pid']==col,'ApexWeight']=WeightDF.loc[apextemp,col]
                else:
                    if outpks.loc[outpks['pid']==col,'wavelen'].iloc[0]%2==0:
                        if apextemp-1 in WeightDF.index:
                            apextemp=apextemp-1
                            outpks.loc[outpks['pid']==col,'ApexWeight']=WeightDF.loc[apextemp,col]
                        else:
                            outpks.loc[outpks['pid']==col,'ApexWeight']=0.0
                    else:
                        outpks.loc[outpks['pid']==col,'ApexWeight']=0.0
                if SubWeight.shape[0]:
                    apextemp=SubWeight.index[0]
                    outpks.loc[outpks['pid']==col,'ApexWeight']=SubWeight.loc[apextemp,col]
                if aucmethod=='exp_wt_log':
                    auctemp=np.trapz(ExpIntDF.loc[SubWeight.index,'ExpInt']*SubWeight[col].tolist(),SubWeight['scan_time'].tolist())
                    outpks.loc[outpks['pid']==col,'AUC']=auctemp
                elif aucmethod=='both':
                    outpks.loc[outpks['pid']==col,'AUC']=np.trapz(np.log(ExpIntDF.loc[SubWeight.index,'ExpInt']*SubWeight[col].tolist()),SubWeight['scan_time'].tolist())
                    outpks.loc[outpks['pid']==col,'AUCalt']=np.trapz(ExpIntDF.loc[SubWeight.index,'ExpInt']*SubWeight[col].tolist(),SubWeight['scan_time'].tolist())
                else:
                    outpks.loc[outpks['pid']==col,'AUC']=np.trapz(np.log(ExpIntDF.loc[SubWeight.index,'ExpInt']*SubWeight[col].tolist()),SubWeight['scan_time'].tolist())
                #handle scale 1
                if ((outpks.loc[(outpks['pid']==col),'scale'].iloc[0]==1)&(outpks.loc[(outpks['pid']==col),'wavelen'].iloc[0]==3))|(SubWeight.shape[0]==1):
                    #if a singlet
                    dt=MS2IdxData.loc[MS2IdxData['PrecursorIdx']==2,'scan_time'].iloc[0]-MS2IdxData.loc[MS2IdxData['PrecursorIdx']==1,'scan_time'].iloc[0]
                    if apextemp not in WeightDF.index:
                        outpks.loc[outpks['pid']==col,'ApexWeight']=1
                    if aucmethod=='exp_wt_log':
                        outpks.loc[outpks['pid']==col,'AUC']=self.DFsub.loc[(self.DFsub['PrecursorIdx']==apextemp)&(self.DFsub['IonID']==self.IonID),'Intensity'].iloc[0]*outpks.loc[outpks['pid']==col,'ApexWeight'].iloc[0]*dt/2
                    elif aucmethod=='both':
                        outpks.loc[outpks['pid']==col,'AUCalt']=np.exp(np.log(self.DFsub.loc[(self.DFsub['PrecursorIdx']==apextemp)&(self.DFsub['IonID']==self.IonID),'Intensity'].iloc[0])*outpks.loc[outpks['pid']==col,'ApexWeight'].iloc[0])*dt/2
                        outpks.loc[outpks['pid']==col,'AUC']=np.log(self.DFsub.loc[(self.DFsub['PrecursorIdx']==apextemp)&(self.DFsub['IonID']==self.IonID),'Intensity'].iloc[0])*outpks.loc[outpks['pid']==col,'ApexWeight'].iloc[0]*dt/2
                    else:
                        outpks.loc[outpks['pid']==col,'AUC']=np.log(self.DFsub.loc[(self.DFsub['PrecursorIdx']==apextemp)&(self.DFsub['IonID']==self.IonID),'Intensity'].iloc[0])*outpks.loc[outpks['pid']==col,'ApexWeight'].iloc[0]*dt/2
                        
    def ValidCurveNeighborCheck(self,r,relscore):
        pkdata=self.potpks.loc[self.potpks['pid']==r].iloc[0]
        self.potpks.loc[self.potpks['pid']==r,'InModel']=False
        self.curvedf=self.curvedf.drop(r,axis=1)
        jdf=pd.DataFrame(None)
        jpids=[r+'_scDN',r+'_scUP',r+'_sh0',r+'_sh1']
        if len(set(jpids)&set(self.potpks['pid']))<4:
            jdf['pid']=jpids
            jdf['LastTransform']=['scale','scale','shoulder','shoulder']
            jdf['LTDirection']=[-1,1,0,1]
            jdf['scale']=[pkdata['scale']-1,pkdata['scale']+1,pkdata['scale'],pkdata['scale']]
            jdf['startscan']=[pkdata['startscan']+1,pkdata['startscan']-1,pkdata['startscan'],pkdata['startscan']+1]
            jdf['endscan']=[pkdata['endscan']-1,pkdata['endscan']+1,pkdata['endscan']-1,pkdata['endscan']]
            jdf['apex']=[pkdata['apex']]*3+[pkdata['apex']+1]
            jdf['wavelen']=[pkdata['wavelen']-2,pkdata['wavelen']+2,pkdata['wavelen']-1,pkdata['wavelen']-1]
            jdf['originpid']=r
            jdf['InSolution']=False
            jdf['InModel']=False
            jdf['Coef']=0
            if len(set(jpids)&set(self.potpks['pid']))>0:
                self.potpks=self.potpks.drop(self.potpks.loc[self.potpks['pid'].isin(jpids)].index)
            self.potpks=pd.concat([self.potpks,jdf]).reset_index(drop=True)
        else:
            jdf=self.potpks.loc[self.potpks['pid'].isin(jpids)]
        jscore=pd.DataFrame(None,columns=['pid','score'])
        jscore['pid']=[r]+jpids
        jscore['score']=-1.0
        jscore.loc[jscore['pid']==r,'score']=relscore
        for jp in jpids:
            if (self.scalemax>=jdf.loc[jdf['pid']==jp,'scale'].iloc[0])&(self.scalemax>=jdf.loc[jdf['pid']==jp,'scale'].iloc[0]>=self.scalemin):
                self.potpks.loc[self.potpks['pid']==jp,'InModel']=True
                self.CurveAdder()
                lmtemp,Xtemp,ytemp=Helper.LinearModelGenerator(self.curvedf, self.potpks.loc[self.potpks['InModel'],'pid'].tolist())
                # if any(lmtemp.coef_<=self.signalmin):
                #     revert=self.RemovePeaksBelowThresh(lmtemp.coef_,report=True)
                #     lmtemp,Xtemp,ytemp=Helper.LinearModelGenerator(self.curvedf, self.potpks.loc[self.potpks['InModel'],'pid'].tolist())
                jscore.loc[jscore['pid']==jp,'score']=lmtemp.score(Xtemp,ytemp)
                self.potpks.loc[self.potpks['pid']==jp,'InModel']=False                    
                self.curvedf=self.curvedf.drop(jp,axis=1)
                # if len(revert)>0:
                #     self.potpks.loc[self.potpks['pid'].isin(revert),'InModel']=True
            else:
                jscore.loc[jscore['pid']==jp,'score']=jscore['score'].min()-1
        relscore=jscore['score'].max()
        bestpid=jscore.loc[jscore['score']==relscore,'pid'].iloc[0]
        self.potpks.loc[self.potpks['pid']==bestpid,'InModel']=True
        self.CurveAdder()
        return relscore
    
    def ResidualFinder(self,lm,X,y,ApexWeights=True):
        RSSDF=pd.DataFrame(None)
        RSSDF['Error']=(lm.predict(X)-y)**2
        RSSDF['PrecursorIdx']=self.curvedf.loc[self.curvedf.sum(axis=1)>0].index
        RSSDF=RSSDF.set_index('PrecursorIdx')
        RSSDF['ErrRatio']=np.inf
        for ex in RSSDF.index:
            if self.curvedf.loc[ex,'Intensity']>0:
                RSSDF.loc[ex,'ErrRatio']=RSSDF.loc[ex,'Error']/self.curvedf.loc[ex,'Intensity']
            else:
                RSSDF.loc[ex,'ErrRatio']=RSSDF.loc[ex,'Error']/self.basesignal
        PeakErrorDF=pd.DataFrame(None,columns=['pid','TotalErr','ApexErr','ErrRatSum','ApexWeight'])
        PeakErrorDF['pid']=self.potpks.loc[self.potpks['InModel'],'pid'].tolist()
        for pid in PeakErrorDF['pid']:
            preidx=self.curvedf.loc[self.curvedf[pid]>0].index.tolist()
            PeakErrorDF.loc[PeakErrorDF['pid']==pid,'TotalErr']=RSSDF.loc[preidx,'Error'].sum()
            PeakErrorDF.loc[PeakErrorDF['pid']==pid,'ErrRatSum']=RSSDF.loc[preidx,'ErrRatio'].sum()
            PeakErrorDF.loc[PeakErrorDF['pid']==pid,'ApexErr']=RSSDF.loc[self.potpks.loc[self.potpks['pid']==pid,'apex'].iloc[0],'Error']
        PeakErrorDF=PeakErrorDF.set_index('pid')
        pidtgts=self.curvedf.columns.tolist()[1:]
        PeakErrorDF['ApexRatio']=0.0
        for px in range(len(pidtgts)):
            denom=self.curvedf.loc[self.potpks.loc[self.potpks['pid']==pidtgts[px],'apex'],'Intensity'].iloc[0]
            if denom>0:
                apexratio=self.curvedf.loc[self.potpks.loc[self.potpks['pid']==pidtgts[px],'apex'],pidtgts[px]].iloc[0]*lm.coef_[px]/denom
            else:
                apexratio=self.curvedf.loc[self.potpks.loc[self.potpks['pid']==pidtgts[px],'apex'],pidtgts[px]].iloc[0]*lm.coef_[px]/self.basesignal
            PeakErrorDF.loc[pidtgts[px],'ApexRatio']=apexratio
        if ApexWeights==True:
            Xco=X
            for idx,pid in enumerate(PeakErrorDF.index.tolist()):
                Xco[:,idx]=Xco[:,idx]*lm.coef_[idx]
            Weights=Xco/(Xco.sum(axis=1)[:,None]+.000001)
            for idx,pid in enumerate(PeakErrorDF.index.tolist()):
                apexix=np.where(RSSDF.index==self.potpks.loc[self.potpks['pid']==pid,'apex'].iloc[0])[0][0]
                PeakErrorDF.loc[pid,'ApexWeight']=Weights[apexix,idx]
        return PeakErrorDF, RSSDF
                    
class SignalWorkflow:
    def __init__(self,hdf5file,RunID=None,scalebounds=[1,6],referencebounds=[0,1200],
                 imputetype=None,basesignalrat=.5,linearspan=3,peak_dist_thresh=None,lmscorethreshold=.7,
                 signalminrat=.05,ObsIdxOpt='Confirmed_Strict',confirmed_gp_margin=3,logopt=False,
                 AcqWindowType='Theoretical',chargerange=[2,6],windowbounds=[800,1600],nbins=50,
                 startfromacq=0,endatacq=50,InitWidthMin=4):
        self.hdf5file=hdf5file
        self.RunID=RunID
        self.scalemin=scalebounds[0]
        self.scalemax=scalebounds[1]
        self.RefMin=referencebounds[0]
        self.RefMax=referencebounds[1]
        self.lmscorethreshold=lmscorethreshold
        self.imputetype=imputetype
        self.linearspan=linearspan
        self.basesignalrat=basesignalrat
        self.signalminrat=signalminrat
        self.ObsIdxOpt=ObsIdxOpt
        self.validacqrng=np.arange(startfromacq,endatacq)
        self.confirmed_gp_margin=confirmed_gp_margin
        self.margin=confirmed_gp_margin*self.scalemax
        self.peak_dist_thresh=peak_dist_thresh
        self.InitWidthMin=InitWidthMin
        if self.peak_dist_thresh==None:
            self.exceptionthreshmin=self.scalemax
        else:
            self.exceptionthreshmin=peak_dist_thresh
        self.AcqWindowType=AcqWindowType
        self.chargerange=chargerange
        self.windowbounds=windowbounds
        self.nbins=nbins
        self.logopt=logopt
        
    def main(self):
        #first get the meta data needed
        self.GrabNecessaryMetaData()
        #then set up the mz window assoc
        MZWFObj=MZWindowFinder(self.chargerange,self.windowbounds,self.nbins)
        if self.AcqWindowType=='Theoretical':
            MZWFObj.main(self.psmobj.dfMS1Adducted)
        elif self.AcqWindowType=='Observed':
            MZWFObj.main(self.psmobj.dfMS1Adducted,self.msdt.AcqWindows)
        #then set up the confirmed indices object
        ConfGPIdxObj=ConfirmedIndexes(PSMMetaMaster=self.psmobj.dfPSMMetaMaster,GPIonAssoc=self.gpiobj.dfGPIonAssoc,IndexData=self.msdt,MZDF=MZWFObj.MZDF,RunID=self.RunID,margin=self.margin)
        #then get the MS2 data connection from file
        self.OpenHDF5File()
        self.MasterCurveDF=pd.DataFrame(None)
        #next do signal processing for each acquisition window
        acqrange=list(set(self.validacqrng)&set(ConfGPIdxObj.confwindowadds.index))
        for acq_j in acqrange:
            self.CurrentAcq=acq_j
            ConfGPIdxObj.AcqIdxInit(AcqIdx=acq_j,GPIonAssoc=self.gpiobj.dfGPIonAssoc,AddGPAssoc=self.psmobj.dfMS1Adducted)
            AcqDF=self.GetAcqWindowProductIonData(AcqIdx=acq_j)
            AcqjIonIDs=set(AcqDF['IonID'].unique().tolist())
            IonIDGrped=ConfGPIdxObj.IonGroupsofCGP.reset_index()
            IonIDGrped['N_ConfGP']=[len(x) for x in IonIDGrped['GPID']]
            IonIDGrped=IonIDGrped.sort_values('N_ConfGP').set_index('IonID')
            self.iondex=list(AcqjIonIDs&set(IonIDGrped.index))
            AcqPeakMaster=pd.DataFrame([[0,0,np.nan,0,0,'NA','NA','NA',0,False,False,0]],columns=['apex','scale','wavelen','startscan','endscan','pid','originpid','LastTransform','LTDirection','InSolution','InModel','Coef'])
            PeakRefMaster=pd.DataFrame([['NA',np.nan,0]],columns=['pid','AddID','ApexDist'])
            for j_ionid in self.iondex:
                self.CurrentIon=j_ionid
                IonSigObj=ProductIonSignal(AcqDF, j_ionid, ConfGPIdxObj,acq_j,
                                            scalebounds=[self.scalemin,self.scalemax],referencebounds=[self.RefMin,self.RefMax],lmscorethreshold=self.lmscorethreshold,
                                            imputetype=self.imputetype,basesignalrat=self.basesignalrat,linearspan=self.linearspan,peak_dist_thresh=self.peak_dist_thresh,
                                            signalminrat=self.signalminrat,ObsIdxOpt=self.ObsIdxOpt,confirmed_gp_margin=self.confirmed_gp_margin)
                outpks,pkref=IonSigObj.Main(IonGPAssoc=self.gpiobj.dfIonGPAssoc,MS2IdxData=self.msdt.MS2Data)
                if outpks is not None:
                    if len(outpks)>0:
                        AcqPeakMaster=pd.concat([AcqPeakMaster,outpks],ignore_index=True)
                if pkref is not None:
                    if len(pkref)>0:
                        PeakRefMaster=pd.concat([PeakRefMaster,pkref],ignore_index=True)
            PeakRefMaster=PeakRefMaster.dropna()
            AcqPeakMaster=AcqPeakMaster.dropna()
            #handle reference docs
            if (PeakRefMaster.shape[0]>0)&(AcqPeakMaster.shape[0]>0):
                PIDtoAddID,PeakRefRed=self.PIDtoAddIDGenerator(AcqPeakMaster,PeakRefMaster,ConfGPIdxObj,acq_j)
                #Add fragment types
                AcqPeakMaster['FragType']=self.gpiobj.dfIonMetaUnique.loc[AcqPeakMaster['IonID'],'fragment_type'].tolist()
                PIDtoAddID['FragType']=self.gpiobj.dfIonMetaUnique.loc[PIDtoAddID['IonID'],'fragment_type'].tolist()
                #initialize apex info
                InitialApexDF=self.EstablishApexLocations(PeakRefRed,PIDtoAddID,ConfGPIdxObj,acq_j)
                #get new solution table from new apex
                NewPIDtoAddID=self.NewPIDtoAddIDGenerator(PeakRefRed,AcqPeakMaster,InitialApexDF)
                #reduce overlaps in solution
                ReducedPIDtoAddID,tempdf=self.ReduceOverlapsInNewPID(NewPIDtoAddID)
                #replace lost duplicates
                ReducedPIDtoAddID,AddWithLostIons=self.ReplaceLostDuplicatePIDinReduced(ReducedPIDtoAddID,tempdf,InitialApexDF)
                # last additions phase
                ReducedPIDtoAddID=self.FinalCurveAdditions(ReducedPIDtoAddID,AddWithLostIons,NewPIDtoAddID)
                #add helpful info to output solution table
                AcqPeakMaster=AcqPeakMaster.set_index('pid')
                ReducedPIDtoAddID['scale']=AcqPeakMaster.loc[ReducedPIDtoAddID['pid'],'scale'].tolist()
                ReducedPIDtoAddID['Coef']=AcqPeakMaster.loc[ReducedPIDtoAddID['pid'],'Coef'].tolist()
                ReducedPIDtoAddID['AcqIdx']=acq_j
                self.MasterCurveDF=pd.concat([self.MasterCurveDF,ReducedPIDtoAddID],ignore_index=True)
        self.HDF5Data.close()
        delattr(self,'HDF5Data')
        self.MasterCurveDF.to_hdf(self.hdf5file,key=self.RunID+'_IonCurves')
            
    def GrabNecessaryMetaData(self):
        self.msdt=Trawler.IndexedMSInfo(runidentifier=self.RunID, h5file=self.hdf5file)
        #get index info
        self.msdt.ReadFromHDF5() 
        #get psm and adduct info
        self.psmobj=Meta.PSMMetaDataTable(self.hdf5file)
        self.psmobj.ListReader()
        #get unique ion info
        self.gpiobj=Meta.GPIonAssociation(self.hdf5file)
        self.gpiobj.ListReader()
        #get associations
        self.gpiobj.AssocReader()
    
    def OpenHDF5File(self):
        self.HDF5Data=tb.open_file(self.hdf5file,mode='r')
        group=getattr(self.HDF5Data.root,self.RunID)
        self.MS2OverallData=group.MS2
    
    def GetAcqWindowProductIonData(self,AcqIdx):
        prodidxsofacq=self.msdt.MS2Data.loc[self.msdt.MS2Data['AcqIdx']==AcqIdx,'SrcProductIdx'].unique().tolist()
        ms2dfobs=pd.DataFrame(None)
        ms2dfobs['Intensity']=[x['Intensity'] for x in self.MS2OverallData.iterrows() if x['ProductIdx'] in prodidxsofacq]
        ms2dfobs['IonID']=[x['IonID'] for x in self.MS2OverallData.iterrows() if x['ProductIdx'] in prodidxsofacq]
        ms2dfobs['PrecursorIdx']=[x['PrecursorIdx'] for x in self.MS2OverallData.iterrows() if x['ProductIdx'] in prodidxsofacq]
        ms2dfobs['Charge']=[x['Charge'] for x in self.MS2OverallData.iterrows() if x['ProductIdx'] in prodidxsofacq]
        ms2dfobs['scan_time']=self.msdt.MS1Data.loc[ms2dfobs['PrecursorIdx'],'scan_time'].tolist()
        ms2dfobs=ms2dfobs.groupby(['PrecursorIdx','IonID']).agg({'Intensity':'sum','scan_time':'mean','Charge':'max'}).reset_index()
        Helper.LogRelativizeIntensity(ms2dfobs,self.basesignalrat)
        return ms2dfobs

    def PIDtoAddIDGenerator(self,AcqPeakMaster,PeakRefMaster,ConfGPIdxObj,acq_j):
        AcqPeakMaster['IonID']=[int(AcqPeakMaster.loc[x,'pid'].split('_')[1]) for x in AcqPeakMaster.index]
        PeakRefMaster['AbsApexD']=abs(PeakRefMaster['ApexDist'])
        PeakRefRed=PeakRefMaster.loc[PeakRefMaster['IonID']!=0].copy()
        PeakRefRed['apex']=0
        for idx in PeakRefRed.index.tolist():
            pidtgt=PeakRefRed.loc[idx,'pid']
            PeakRefRed.loc[idx,'apex']=AcqPeakMaster.loc[AcqPeakMaster['pid']==pidtgt,'apex'].iloc[0]
        #make initial solution table
        PIDtoAddID=pd.DataFrame([['str',np.nan,0,0,0.0,0.0,0,0]],columns=['pid','apex','wavelen','ApexDist','AUC','score','AddID','IonID'])
        for adx in ConfGPIdxObj.confwindowadds.loc[acq_j]:
            PRM=PeakRefRed.loc[PeakRefRed['AddID']==adx]
            for lbl,ion in PRM.groupby('IonID'):
                minval=ion['AbsApexD'].min()
                phx=ion.loc[ion['AbsApexD']==minval].index.tolist()
                for ph in phx:
                    pidtgt=ion.loc[ph,'pid']
                    piddata=AcqPeakMaster.loc[AcqPeakMaster['pid']==pidtgt].iloc[0]
                    PIDtoAddID.loc[PIDtoAddID.index.max()+1]=[pidtgt,piddata['apex'],piddata['wavelen'],minval,piddata['AUC'],piddata['score'],adx,lbl]
        PIDtoAddID=PIDtoAddID.dropna()
        for pidrpt,rpt in PIDtoAddID.groupby('pid'):
            if rpt['AddID'].nunique()>1:
                for adx in rpt['AddID'].unique():
                    red=PeakRefRed.loc[(PeakRefRed['AddID']==adx)&(PeakRefRed['IonID']==rpt['IonID'].iloc[0])]
                    minval=red['AbsApexD'].min()
                    red=red.loc[red['AbsApexD']!=minval]
                    if red.shape[0]>1:
                        minval=red['AbsApexD'].min()
                        phx=red.loc[red['AbsApexD']==minval].index.tolist()
                        for ph in phx:
                            pidtgt=red.loc[ph,'pid']
                            piddata=AcqPeakMaster.loc[AcqPeakMaster['pid']==pidtgt].iloc[0]
                            PIDtoAddID.loc[PIDtoAddID.index.max()+1]=[pidtgt,piddata['apex'],piddata['wavelen'],minval,piddata['AUC'],piddata['score'],adx,rpt['IonID'].iloc[0]]
        return PIDtoAddID, PeakRefRed
    
    def EstablishApexLocations(self,PeakRefRed,PIDtoAddID,ConfGPIdxObj,acq_j):
        InitialApexDF=pd.DataFrame([[0,np.nan,0]],columns=['AddID','apex','width'])
        # for a given confirmed addid
        for addid_j in PeakRefRed['AddID'].unique():
            gp_j=self.psmobj.dfMS1Adducted.loc[addid_j,'GPID']
            interferingadds=list(set(PIDtoAddID.loc[PIDtoAddID['pid'].isin(PIDtoAddID.loc[PIDtoAddID['AddID']==addid_j,'pid']),'AddID'].unique())-set([addid_j]))
            interferinggps=self.psmobj.dfMS1Adducted.loc[interferingadds,'GPID'].tolist()
            distinctobs=set(PeakRefRed.loc[PeakRefRed['AddID']==addid_j,'IonID'])-set(PeakRefRed.loc[PeakRefRed['AddID']!=addid_j,'IonID'])
            if len(interferinggps)==0:
                possiblydistinctproducts=set(self.gpiobj.dfGPIonAssoc.loc[self.gpiobj.dfGPIonAssoc['GPID']==gp_j,'IonID'].sum())
            else:
                possiblydistinctproducts=set(self.gpiobj.dfGPIonAssoc.loc[self.gpiobj.dfGPIonAssoc['GPID']==gp_j,'IonID'].sum())-set(self.gpiobj.dfGPIonAssoc.loc[self.gpiobj.dfGPIonAssoc['GPID'].isin(interferinggps),'IonID'].sum())
            distinctobsprod=list(possiblydistinctproducts&distinctobs)
            #get the distinctly observed products
            if len(distinctobsprod)>0:
                tentativeapexpid=PIDtoAddID.loc[(PIDtoAddID['IonID'].isin(distinctobsprod))]
                #get prod curves in score threshold
                validscore=tentativeapexpid.loc[(tentativeapexpid['score']>=self.lmscorethreshold)&(tentativeapexpid['score']<1)]
                if validscore.shape[0]==1:
                    initapex=validscore['apex'].iloc[0]
                    initwidth=validscore['wavelen'].iloc[0]
                elif validscore.shape[0]>1:
                    #if multiple get closest to observed
                    closest=validscore.loc[validscore['ApexDist']==validscore['ApexDist'].min()]
                    if closest.shape[0]==1:
                        initapex=closest['apex'].iloc[0]
                        initwidth=closest['wavelen'].iloc[0]
                    elif closest.shape[0]>1:
                        #if multiple closest get best score with biggest model
                        bestlm=closest.loc[closest['score']==closest['score'].max()]
                        initapex=bestlm.loc[bestlm['AUC']==bestlm['AUC'].max(),'apex'].iloc[0]
                        initwidth=bestlm.loc[bestlm['AUC']==bestlm['AUC'].max(),'wavelen'].iloc[0]
                else:
                    #if there wasnt a curve that passed that didnt have score=1, include score=1 
                    validscore=tentativeapexpid.loc[(tentativeapexpid['score']>=self.lmscorethreshold)]
                    if validscore.shape[0]==1:
                        initapex=validscore['apex'].iloc[0]
                        initwidth=self.InitWidthMin
                    elif validscore.shape[0]>1:
                        closest=validscore.loc[validscore['ApexDist']==validscore['ApexDist'].min()]
                        initwidth=validscore['apex'].max()-validscore['apex'].min()+1
                        if closest.shape[0]==1:
                            initapex=closest['apex'].iloc[0]
                        elif closest.shape[0]>1:
                            initapex=closest['apex'].mean() 
            elif len(list(distinctobs))>0:
                #if no distinctly observed distinct product ions with good enough scores
                #look at just distinctly observed ions
                tentativeapexpid=PIDtoAddID.loc[(PIDtoAddID['IonID'].isin(list(distinctobs)))]
                validscore=tentativeapexpid.loc[(tentativeapexpid['score']>=self.lmscorethreshold)&(tentativeapexpid['score']<1)]
                if validscore.shape[0]==1:
                    initapex=validscore['apex'].iloc[0]
                    initwidth=validscore['wavelen'].iloc[0]
                elif validscore.shape[0]>1:
                    #if multiple get closest to observed
                    closest=validscore.loc[validscore['ApexDist']==validscore['ApexDist'].min()]
                    if closest.shape[0]==1:
                        initapex=closest['apex'].iloc[0]
                        initwidth=closest['wavelen'].iloc[0]
                    elif closest.shape[0]>1:
                        #if multiple closest get best score with biggest model
                        bestlm=closest.loc[closest['score']==closest['score'].max()]
                        initapex=bestlm.loc[bestlm['AUC']==bestlm['AUC'].max(),'apex'].iloc[0]
                        initwidth=bestlm.loc[bestlm['AUC']==bestlm['AUC'].max(),'wavelen'].iloc[0]
                else:
                    #if there wasnt a curve that passed that didnt have score=1, check score=1 
                    validscore=tentativeapexpid.loc[(tentativeapexpid['score']==1)]
                    if validscore.shape[0]==1:
                        initapex=validscore['apex'].iloc[0]
                        initwidth=self.InitWidthMin
                    elif validscore.shape[0]>1:
                        closest=validscore.loc[validscore['ApexDist']==validscore['ApexDist'].min()]
                        if closest.shape[0]==1:
                            initapex=closest['apex'].iloc[0]
                            initwidth=validscore['apex'].max()-validscore['apex'].min()+1
                        elif closest.shape[0]>1:
                            initapex=closest['apex'].mean()
                            initwidth=validscore['apex'].max()-validscore['apex'].min()+1
                    else:
                        initapex=ConfGPIdxObj.confirmed_gps.loc[(ConfGPIdxObj.confirmed_gps['AddID']==addid_j)&(ConfGPIdxObj.confirmed_gps['AcqIdx']==acq_j),'PrecursorIdx'].mean()
                        initwidth=self.InitWidthMin
            if initwidth<self.InitWidthMin:
                initwidth=self.InitWidthMin
            InitialApexDF.loc[InitialApexDF.index.max()+1]=[addid_j,initapex,initwidth]
        InitialApexDF=InitialApexDF.dropna()
        return InitialApexDF
    
    def NewPIDtoAddIDGenerator(self,PeakRefRed,AcqPeakMaster,InitialApexDF):
        NewPIDtoAddID=pd.DataFrame([['str',np.nan,0,0,0.0,0.0,0,0,0]],columns=['pid','apex','wavelen','ApexDist','AUC','score','AddID','IonID','FragType'])
        for addid_j in PeakRefRed['AddID'].unique():
            gp_j=self.psmobj.dfMS1Adducted.loc[addid_j,'GPID']
            ions=self.gpiobj.dfGPIonAssoc.loc[gp_j,'IonID']
            subset=AcqPeakMaster.loc[AcqPeakMaster['IonID'].isin(ions),['pid','apex','wavelen','AUC','score','IonID','FragType']]
            subset['ApexDist']=abs(subset['apex']-InitialApexDF.loc[InitialApexDF['AddID']==addid_j,'apex'].iloc[0])
            subset['AddID']=addid_j
            rng=np.floor(subset.loc[subset['ApexDist']==0,'wavelen'].max()/2)
            subset=subset.loc[(subset['ApexDist']<=rng)&(subset['AUC']>0)&(subset['score']>self.lmscorethreshold)]
            InitialApexDF.loc[InitialApexDF['AddID']==addid_j,'width']=subset.loc[subset['ApexDist']==0,'wavelen'].max()
            NewPIDtoAddID=pd.concat([NewPIDtoAddID,subset],ignore_index=True)
        NewPIDtoAddID=NewPIDtoAddID.dropna()
        return NewPIDtoAddID
    
    def ReduceOverlapsInNewPID(self,NewPIDtoAddID):
        ReducedPIDtoAddID=pd.DataFrame([['str',np.nan,0,0,0.0,0.0,0,0,0]],columns=['pid','apex','wavelen','ApexDist','AUC','score','AddID','IonID','FragType'])
        for iontgt in NewPIDtoAddID['IonID'].unique():
            ion=NewPIDtoAddID.loc[NewPIDtoAddID['IonID']==iontgt].copy()
            ion['N_rpt_pid']=[sum(ion['pid']==pidtgt) for pidtgt in ion['pid']]
            ion=ion.sort_values(['ApexDist','N_rpt_pid'])
            missingadds=ion['AddID'].unique()
            for apexd in ion['ApexDist'].unique():
                addhits=[]
                sub1=ion.loc[ion['ApexDist']==apexd].sort_values('N_rpt_pid')
                for pidtgt in sub1['pid'].unique():
                    red=sub1.loc[(sub1['pid']==pidtgt)&(sub1['AddID'].isin(missingadds))]
                    if red.shape[0]>0:
                        for ax in red.index:
                            ReducedPIDtoAddID.loc[ReducedPIDtoAddID.index.max()+1]=red.loc[ax]
                            addhits+=[red.loc[ax,'AddID']]
                missingadds=list(set(missingadds)-set(addhits))
        ReducedPIDtoAddID=ReducedPIDtoAddID.dropna().reset_index(drop=True)
        #get rid of duplicates
        tempdf=ReducedPIDtoAddID[ReducedPIDtoAddID['pid'].duplicated(keep=False)]
        for ax in ReducedPIDtoAddID['AddID']:
            sub1=ReducedPIDtoAddID.loc[ReducedPIDtoAddID['AddID']==ax]
            for ix in sub1['IonID']:
                ion=sub1.loc[sub1['IonID']==ix]
                allionidx=set(ReducedPIDtoAddID.loc[(ReducedPIDtoAddID['AddID']==ax)&(ReducedPIDtoAddID['IonID']==ix)].index.tolist())
                if ion.shape[0]==1:
                    keep=set(ReducedPIDtoAddID.loc[(ReducedPIDtoAddID['AddID']==ax)&(ReducedPIDtoAddID['IonID']==ix)&(ReducedPIDtoAddID['pid']==ion['pid'].iloc[0])].index.tolist())
                else:
                    ionsub=ion.loc[ion['wavelen']==ion['wavelen'].max()]
                    if ionsub.shape[0]>1:
                        keep=set(ReducedPIDtoAddID.loc[(ReducedPIDtoAddID['AUC']==ionsub['AUC'].max())&(ReducedPIDtoAddID['AddID']==ax)&(ReducedPIDtoAddID['IonID']==ix)].index.tolist())
                    else:
                        keep=set(ReducedPIDtoAddID.loc[(ReducedPIDtoAddID['AddID']==ax)&(ReducedPIDtoAddID['IonID']==ix)&(ReducedPIDtoAddID['pid']==ionsub['pid'].iloc[0])].index.tolist())
                drop=list(allionidx-keep)
                ReducedPIDtoAddID=ReducedPIDtoAddID.drop(drop)
        ReducedPIDtoAddID.drop_duplicates('pid',keep=False,inplace=True)
        return ReducedPIDtoAddID, tempdf
    
    def ReplaceLostDuplicatePIDinReduced(self,ReducedPIDtoAddID,tempdf,InitialApexDF):
        AddWithLostIons=pd.DataFrame([[np.nan,np.nan]],columns=['AddID','IonID'])
        for pidtgt in tempdf['pid'].unique():
            sub1=tempdf.loc[tempdf['pid']==pidtgt]
            #if there's only one apex then go immediately into width determination
            if sub1['ApexDist'].nunique()==1:
                ApexSub=InitialApexDF.loc[InitialApexDF['AddID'].isin(sub1['AddID'])]
                sub3=ApexSub.loc[ApexSub['width']==ApexSub['width'].max()]
                if sub3.shape[0]==1:
                    hitadd=sub3['AddID'].iloc[0]
                    missadd=list(set(sub1['AddID'])-set([hitadd]))
                    ReducedPIDtoAddID.loc[ReducedPIDtoAddID.index.max()+1]=sub1.loc[sub1['AddID']==hitadd].iloc[0]
                else:
                    psmsub=self.psmobj.dfPSMMetaMaster.loc[(self.psmobj.dfPSMMetaMaster['AddID'].isin(sub3['AddID']))&(self.psmobj.dfPSMMetaMaster['RunID']==self.RunID)]
                    hitadd=psmsub.loc[psmsub['precursor_intensity']==psmsub['precursor_intensity'].max(),'AddID'].iloc[0]
                    missadd=list(set(sub1['AddID'])-set([hitadd]))
                    ReducedPIDtoAddID.loc[ReducedPIDtoAddID.index.max()+1]=sub1.loc[sub1['AddID']==hitadd].iloc[0]
            else:
                minval=sub1['ApexDist'].min()
                sub2=sub1.loc[sub1['ApexDist']==minval]
                if sub2.shape[0]==1:
                    hitadd=sub2['AddID'].iloc[0]
                    missadd=list(set(sub1['AddID'])-set([hitadd]))
                    ReducedPIDtoAddID.loc[ReducedPIDtoAddID.index.max()+1]=sub1.loc[sub1['AddID']==hitadd].iloc[0]
                else:
                    ApexSub=InitialApexDF.loc[InitialApexDF['AddID'].isin(sub2['AddID'])]
                    sub3=ApexSub.loc[ApexSub['width']==ApexSub['width'].max()]
                    if sub3.shape[0]==1:
                        hitadd=sub3['AddID'].iloc[0]
                        missadd=list(set(sub1['AddID'])-set([hitadd]))
                        ReducedPIDtoAddID.loc[ReducedPIDtoAddID.index.max()+1]=sub1.loc[sub1['AddID']==hitadd].iloc[0]
                    else:
                        psmsub=self.psmobj.dfPSMMetaMaster.loc[(self.psmobj.dfPSMMetaMaster['AddID'].isin(sub3['AddID']))&(self.psmobj.dfPSMMetaMaster['RunID']==self.RunID)]
                        hitadd=psmsub.loc[psmsub['precursor_intensity']==psmsub['precursor_intensity'].max(),'AddID'].iloc[0]
                        missadd=list(set(sub1['AddID'])-set([hitadd]))
                        ReducedPIDtoAddID.loc[ReducedPIDtoAddID.index.max()+1]=sub1.loc[sub1['AddID']==hitadd].iloc[0]
            for m in missadd:
                AddWithLostIons.loc[AddWithLostIons.index.max()+1]=[m,sub1['IonID'].iloc[0]]
        AddWithLostIons=AddWithLostIons.dropna()
        return ReducedPIDtoAddID, AddWithLostIons
    
    def FinalCurveAdditions(self,ReducedPIDtoAddID,AddWithLostIons,NewPIDtoAddID):
        for iontgt in AddWithLostIons['IonID'].unique():
            addtgt=AddWithLostIons.loc[AddWithLostIons['IonID']==iontgt,'AddID'].tolist()
            for ax in addtgt:
                if ReducedPIDtoAddID.loc[(ReducedPIDtoAddID['AddID']==ax)&(ReducedPIDtoAddID['IonID']==iontgt)].shape[0]>0:
                    AddWithLostIons=AddWithLostIons.drop(AddWithLostIons.loc[(AddWithLostIons['AddID']==ax)&(AddWithLostIons['IonID']==iontgt)].index)
            addtgt=AddWithLostIons.loc[AddWithLostIons['IonID']==iontgt,'AddID'].tolist()
            pidinuse=set(ReducedPIDtoAddID.loc[ReducedPIDtoAddID['IonID']==iontgt,'pid'])
            pidposs=set(NewPIDtoAddID.loc[(NewPIDtoAddID['IonID']==iontgt)&(NewPIDtoAddID['AddID'].isin(addtgt)),'pid'])
            pidrem=list(pidposs-pidinuse)
            if len(pidrem)>0:
                pot=NewPIDtoAddID.loc[(NewPIDtoAddID['pid'].isin(pidrem))&(NewPIDtoAddID['AddID'].isin(addtgt))]
                while (pot.shape[0]>0)&(len(addtgt)>0):
                    ptx=pot.loc[pot['ApexDist']==pot['ApexDist'].min()].index.tolist()[0]
                    ReducedPIDtoAddID.loc[ReducedPIDtoAddID.index.max()+1]=pot.loc[ptx]
                    addtgt=list(set(addtgt)-set([pot.loc[ptx,'AddID']]))
                    pot=pot.drop(pot.loc[pot['AddID']==pot.loc[ptx,'AddID']].index)
            else:
                for ax in addtgt:
                    AddWithLostIons=AddWithLostIons.drop(AddWithLostIons.loc[(AddWithLostIons['AddID']==ax)&(AddWithLostIons['IonID']==iontgt)].index)
        return ReducedPIDtoAddID
    
    def ReadMS1(HDF5Connection,RunID,tgtcol=None,tgtvals=None,graball=False):
        groupinit=getattr(HDF5Connection.root,RunID)
        tbgroup=groupinit.MS1
        ms1data=pd.DataFrame(None)
        if graball==False:
            ms1data['AddID']=[x['AddID'] for x in tbgroup.iterrows() if x[tgtcol] in tgtvals]
            ms1data['Intensity']=[x['Intensity'] for x in tbgroup.iterrows() if x[tgtcol] in tgtvals]
            ms1data['PrecursorIdx']=[x['PrecursorIdx'] for x in tbgroup.iterrows() if x[tgtcol] in tgtvals]
            ms1data['Charge']=[x['Charge'] for x in tbgroup.iterrows() if x[tgtcol] in tgtvals]
        else:
            ms1data['AddID']=[x['AddID'] for x in tbgroup.iterrows()]
            ms1data['Intensity']=[x['Intensity'] for x in tbgroup.iterrows()]
            ms1data['PrecursorIdx']=[x['PrecursorIdx'] for x in tbgroup.iterrows()]
            ms1data['Charge']=[x['Charge'] for x in tbgroup.iterrows()]    
        return ms1data
    
def TheoreticalCurveGenerator(scale,apex,magnitude=1,wavelen=None,noiselvl=0,refpid='NA',refion=0,dt=0.15,exp=True):
    if wavelen is None:
        wavelen=scale*2+1
    outdf=pd.DataFrame(None,columns=['Intensity','pid','PrecursorIdx','scan_time','IonID'])
    tgtrng=np.arange(apex-scale,apex-scale+wavelen)
    outdf['PrecursorIdx']=tgtrng
    outdf['scan_time']=outdf['PrecursorIdx']*dt
    outdf['pid']=refpid
    outdf['IonID']=refion
    outdf['Intensity']=0.0
    wave=sig.ricker(wavelen,scale)*magnitude
    if noiselvl>0:
        adj=np.random.normal(0,scale=noiselvl,size=len(wave))
        wave=wave+wave*adj
    if exp==True:
        outdf.loc[outdf['PrecursorIdx']==tgtrng,'Intensity']=np.exp(wave)-1
    else:
        outdf.loc[outdf['PrecursorIdx']==tgtrng,'Intensity']=wave
    if any(wave<0):
        outdf.loc[outdf['Intensity']<0,'Intensity']=0
    return outdf

class FakeConfObj:
    def __init__(self,refids,apexes,refgps=0,acqidxs=0):
        self.refids=refids,
        self.apexes=apexes
        self.refgps=refgps
        self.acqidxs=acqidxs
        
    def main(self):
        self.FakeConfGPs()
        self.FakeIonFrom()
        
    def FakeConfGPs(self):
        confirmed_gps=pd.DataFrame(None,columns=['PrecursorIdx','AcqIdx','AddID','GPID'])
        confirmed_gps['PrecursorIdx']=self.apexes
        confirmed_gps['AddID']=self.refgps
        confirmed_gps['GPID']=self.refgps
        confirmed_gps['AcqIdx']=self.acqidxs
        self.confirmed_gps=confirmed_gps
        
    def FakeIonFrom(self):
        CGPfromIoninAcqDF=pd.DataFrame(None,columns=['IonID','GPID','AddID'])
        CGPfromIoninAcqDF['IonID']=self.refids
        CGPfromIoninAcqDF['GPID']=self.refgps
        CGPfromIoninAcqDF['AddID']=self.refgps
        self.CGPfromIoninAcqDF=CGPfromIoninAcqDF

class FakeIonGPAssoc:
    def __init__(self,refids,refgps):
        self.refids=refids
        self.refgps=refgps
        
    def main(self):
        dfIonGPAssocPre=pd.DataFrame(None,columns=['IonID','GPID'])
        dfIonGPAssocPre['IonID']=self.refids
        dfIonGPAssocPre['GPID']=self.refgps
        dfIonGPAssoc=dfIonGPAssocPre.groupby('IonID')['GPID'].apply(list).reset_index()
        self.dfIonGPAssoc=dfIonGPAssoc
        
class FakeMS2Index:
    def __init__(self,refrange,dt,acqidxs):
        self.refrange=refrange
        self.dt=dt
        self.acqidxs=acqidxs
        
    def main(self):
        MS2Data=pd.DataFrame([[0,np.nan,np.nan]],columns=['PrecursorIdx','AcqIdx','scan_time'])
        ddt=self.dt/len(self.acqidxs)
        for j in self.refrange:
            for idx,k in enumerate(self.acqidxs):
                MS2Data.loc[MS2Data.index.max()+1]=[j,k,j*self.dt+idx*ddt]
        self.MS2Data=MS2Data.dropna()
        
# class IonAllocationByTheoreticalCurve:
#     def __init__(self,observed,observedion,scalemin=2,scalemax=31,exceptionthreshmin=3,PreIDMin=0,PreIDMax=None,imputetype=None,basesignal=0):
#         self.observed=observed
#         self.observedion=observedion
#         self.scalemin=scalemin
#         self.scalemax=scalemax
#         self.PreIDMin=PreIDMin
#         self.PreIDMax=PreIDMax
#         self.imputetype=imputetype
#         self.basesignal=basesignal
#         self.exceptionthreshmin=exceptionthreshmin
        
#     def main(self,IonGPAssoc,FragEffTable,UseFragEff=True):
#         self.TheoreticalCurveGenerationMain()
#         self.IonAllocationMain(FragEffTable,IonGPAssoc,UseFragEff)
        
#     def TheoreticalCurveGenerationMain(self):
#         self.curvesout=pd.DataFrame([[0,'pidNA_0',0,0,0]],columns=['scanid','peakid','PeakIntensity','AddID','GPID'])
#         self.pkdatadf=pd.DataFrame(None,columns=['peakid','AddID','scale','startscan','endscan','apex','peakmass'])
#         for labels, dfi in self.observed.groupby("AddID"):
#             self.TheoreticalCurvesOfAddID(labels,dfi)
#         self.curvesout=self.curvesout.loc[self.curvesout['peakid']!='pidNA_0']
#         self.pkdatadf['conversion']=[self.observed.loc[self.observed['AddID']==prow['AddID'],'intensity'].sum()/self.pkdatadf.loc[self.pkdatadf['AddID']==prow['AddID'],'peakmass'].sum() for index, prow in self.pkdatadf.iterrows()]
#         self.curvesout['PeakIntensity']=[crow['PeakIntensity']*self.pkdatadf.loc[self.pkdatadf['peakid']==crow['peakid'],'conversion'].iloc[0] for index, crow in self.curvesout.iterrows()]    
            
#     def CWTMaker(self,signalvector):
#         self.cwtmatr=sig.cwt(signalvector,sig.ricker,np.arange(self.scalemin,self.scalemax))
#         self.cwtpeak=sig.find_peaks_cwt(signalvector,np.arange(self.scalemin,self.scalemax))
    
#     def tempcurveMaker(self,scale,lowerbound,upperbound,wave,cwtmass,pid,labels,gplab):
#         tempcurve=pd.DataFrame(None,columns=['scanid','peakid','PeakIntensity','AddID','GPID'])
#         tempcurve['scanid']=list(np.arange(lowerbound,upperbound))
#         tempcurve['PeakIntensity']=list(wave/np.max(wave)*cwtmass)
#         tempcurve['peakid']=[pid]*tempcurve.shape[0]
#         tempcurve['AddID']=[labels]*tempcurve.shape[0]
#         tempcurve['GPID']=[gplab]*tempcurve.shape[0]
#         return tempcurve
                    
#     def EmptyIndexFinder(self,dfi,AddIDTarget):
#         d={'PeakIntensity':'SummedIntensity','peakid':'peakids'}
#         grpedgp=self.curvesout.loc[self.curvesout['AddID']==AddIDTarget].groupby(['scanid']).agg({'PeakIntensity':'sum','peakid':'unique'}).rename(columns=d)
#         #find the empty intensity scanids, first get the ones that are totally missing
#         emptyids=list(set(dfi['scanid'])-set(grpedgp.index.tolist()))
#         #find the ones without any peakids assoc, find closest scans with peaks
#         #if within exception distance, include, if not add to exception list
#         for j in emptyids:
#             emptydist=abs(emptyids-grpedgp.index)
#             if np.min(emptydist)<=self.exceptionthreshmin:
#                 grpedgp=pd.concat([grpedgp,pd.DataFrame([[0,grpedgp['peakids'].iloc[emptydist.argmin()].tolist()]],columns=['SummedIntensity','peakids'],index=[j])])
#         #now take the ones that are seen but have a total intensity of 0
#         emptyout=grpedgp.loc[dfi['scanid']].loc[(grpedgp['SummedIntensity']==0)].index.tolist()
#         return emptyout
    
#     def TheoreticalCurveAdjuster(self,labels,gplab,emptyids,ignorepid=[]):
#         pids=self.curvesout.loc[(self.curvesout['scanid'].isin(emptyids)) & (self.curvesout['AddID']==labels),'peakid']
#         #look at peak with most missing associations
#         if len(ignorepid)>0:
#             pids=pids.loc[~pids.isin(ignorepid)]
#         if pids.shape[0]==0:
#             self.ContinueCheck=False
#         else:
#             self.ContinueCheck=True
#             pidtgt=pids.mode()[0]
#             pkdata=self.pkdatadf.loc[self.pkdatadf['peakid']==pidtgt]
#             tempcurve=pd.DataFrame(None,columns=['scanid','peakid','PeakIntensity','AddID','GPID'])
#             #if both pk ends are missing then increase scale
#             #otherwise, just shift the peak
#             if ((pkdata['startscan'].iloc[0] in emptyids) & (pkdata['endscan'].iloc[0] in emptyids)):
#                 pkmaxloc=pkdata['scale'].iloc[0]
#                 if (pkmaxloc<self.cwtmatr.shape[0]):
#                     scale=pkmaxloc+1
#                     lowerbound=pkdata['apex'].iloc[0]-scale
#                     upperbound=pkdata['apex'].iloc[0]+scale+1
#                     cwtmass=self.cwtmatr[pkmaxloc,lowerbound:(upperbound+1)].sum()
#                     wave=list(sig.ricker(scale*2+1,scale))
#                 else:
#                     ignorepid+=[pidtgt]
#             else:
#                 pkmaxloc=pkdata['scale'].iloc[0]-1
#                 scale=pkdata['scale'].iloc[0]
#                 wave=list(sig.ricker(scale*2,scale))
#                 if (pkdata['startscan'].iloc[0] in emptyids):
#                     lowerbound=pkdata['startscan'].iloc[0]
#                     upperbound=pkdata['endscan'].iloc[0]-1
#                 else:
#                     lowerbound=pkdata['startscan'].iloc[0]+1
#                     upperbound=pkdata['endscan'].iloc[0]
#                 cwtmass=self.cwtmatr[pkmaxloc,lowerbound:(upperbound+1)].sum()                
#             tempcurve=self.tempcurveMaker(scale,lowerbound,upperbound,wave,cwtmass,pkdata['peakid'].iloc[0],labels,gplab)
#             self.pkdatadf.loc[self.pkdatadf['peakid']==pidtgt]=[pidtgt,pkdata['AddID'].iloc[0],scale,lowerbound,upperbound,pkdata['apex'].iloc[0],cwtmass]
#             self.curvesout=self.curvesout.loc[self.curvesout['peakid']!=pidtgt]
#             self.curvesout=pd.concat([self.curvesout,tempcurve])
            
        
    
#     def TheoreticalCurveGenerator(self,labels,gplab,pkn,pk,ms1temp):
#         pid="pid"+str(labels)+"_"+str(pkn)
#         pkmaxloc=self.cwtmatr[:,pk].argmax()
#         scale=pkmaxloc+self.scalemin
#         lowerbound=ms1temp.index.values.min()+pk-scale
#         upperbound=ms1temp.index.values.min()+pk+scale+1
#         cwtmass=self.cwtmatr[pkmaxloc,lowerbound:upperbound+1].sum()
#         wave=list(sig.ricker(scale*2+1,scale))
#         tempcurve=self.tempcurveMaker(scale,lowerbound,upperbound,wave,cwtmass,pid,labels,gplab)
#         self.curvesout=pd.concat([self.curvesout,tempcurve])
#         self.pkdatadf.loc[len(self.pkdatadf.index)]=[pid,labels,scale,lowerbound,upperbound,pk,cwtmass]
    
#     def TheoreticalCurvesOfAddID(self,labels,dfi):
#         gplab=labels #CHANGE LINES FOR AddID vs GPID changes
#         ms1temp=Helper.BasicImputation(dfi,dataname='intensity',refname='scanid',RefMin=self.PreIDMin,RefMax=self.PreIDMax,imputetype=self.imputetype,basesignal=self.basesignal)
#         ms1temp['intensity']=ms1temp['intensity']+0.000001 # replace with basesignal variable in helper.basicimputation
#         self.CWTMaker(ms1temp['intensity'].tolist())
#         if len(self.cwtpeak)>0:
#             for pkn, pk in enumerate(self.cwtpeak):
#                 if np.max(self.cwtmatr[:,pk])>self.basesignal+1:
#                     self.TheoreticalCurveGenerator(labels,gplab,pkn,pk,ms1temp)
#         else:
#             pid='pid'+str(labels)+'_0'
#             posidx=ms1temp.loc[ms1temp['intensity']>0].index.tolist()
#             tempcurve=pd.DataFrame(None,columns=['scanid','peakid','PeakIntensity','AddID','GPID'])
#             tempcurve['scanid']=ms1temp['scanid'].loc[posidx].tolist()
#             tempcurve['PeakIntensity']=ms1temp['intensity'].loc[posidx].tolist()
#             tempcurve['peakid']=[pid]*tempcurve.shape[0]
#             tempcurve['AddID']=[labels]*tempcurve.shape[0]
#             tempcurve['GPID']=[gplab]*tempcurve.shape[0]
#             self.curvesout=pd.concat([self.curvesout,tempcurve])
#             self.pkdatadf.loc[len(self.pkdatadf.index)]=[pid,labels,np.nan,np.min(posidx),np.max(posidx),np.nan,tempcurve['PeakIntensity'].sum()]
#         emptyids=self.EmptyIndexFinder(dfi,labels)
#         self.ContinueCheck=True
#         while (len(emptyids)>0) & (self.ContinueCheck):
#             self.TheoreticalCurveAdjuster(labels,gplab,emptyids)
#             #repeat until empty values gone 
#             emptyids=self.EmptyIndexFinder(dfi,labels)
    
#     def IonAllocationMain(self,FragEffTable,IonGPAssoc,UseFragEff=True):
#         if UseFragEff!=True:
#             FragEffTable['FragEff']=1
#         self.adjion_int=pd.DataFrame([[0,0,0,np.nan,'NA']],columns=['GPID','IonID','scanid','adj_intensity','peakid'])
#         for scid, dfi in self.observedion.groupby(['scanid']):
#             gplist=self.curvesout.loc[self.curvesout['scanid']==scid,['GPID','PeakIntensity','peakid']].copy()
#             for ion in dfi['IonID']:
#                self.IonAllocationForIndIonID(dfi,ion,scid,gplist,FragEffTable,IonGPAssoc)
#         self.adjion_int=self.adjion_int.dropna()
    
#     def IonAllocationForIndIonID(self,dfi,ion,scid,gplist,FragEffTable,IonGPAssoc):
#         hits=gplist.loc[gplist['GPID'].isin(IonGPAssoc.loc[ion,'GPID'])].copy()
#         hits['FragEfficiency']=[FragEffTable.loc[(FragEffTable['GPID']==gp) & (FragEffTable['IonID']==ion), 'FragEff'].iloc[0] for gp in hits['GPID'].tolist()]
#         hits['FragRatio']=Helper.SumRatio(hits['FragEfficiency'])
#         hits['Ratio']=Helper.SumRatio(hits['PeakIntensity']*hits['FragRatio'])
#         tempint=pd.DataFrame(None,columns=['GPID','IonID','scanid','adj_intensity','peakid'])
#         tempint['adj_intensity']=list(dfi.loc[dfi['IonID']==ion,'intensity'].tolist()*np.array(hits['Ratio']))
#         tempint['GPID']=hits['GPID'].tolist()
#         tempint['peakid']=hits['peakid'].tolist()
#         tempint['IonID']=ion
#         tempint['scanid']=scid[0]
#         self.adjion_int=pd.concat([self.adjion_int,tempint],ignore_index=True)