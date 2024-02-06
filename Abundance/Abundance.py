### GlyLine/Abundance.py 
#this file contains functions for iterating through data sets

from GlyLine.Helper import Helper as GLH
from GlyLine.Meta import Meta
from GlyLine.Trawler import Trawler
import pandas as pd
import numpy as np
import tables as tb



#this finds the ions related to a single glycopeptide
class AbundanceCalc:
    def trap_sum_DF(x,label1,labelt='Time'):
        temp=np.trapz(x.sort_values(labelt)[label1],x.sort_values(labelt)[labelt])
        return pd.Series(temp, index=['AUC'])
   
    def trap_sum_arrays(x,y):
        temp=np.trapz(y[x.argsort()],x[x.argsort()])
        return temp
    
    def trap_sum_list(self,x,y):
        temp=self.trap_sum_arrays(np.array(x),np.array(y))
        return temp
    
class Associator:
    def __init__(self,hdf5file):
        self.h5file=hdf5file
        self.dfGlyIonAssoc=pd.read_hdf(self.h5file,key='GlyIonAssoc',mode='a')
        
    def Gatherer(self):
        return

            
#Ion XIC maker
class IonXIC(AbundanceCalc):
    def __init__(self,IonID,GPID,MS2=None,decon=0,overlap=0,runID=None,h5file=None):
        self.IonID=IonID
        self.decon=decon
        self.overlap=overlap
        self.strtarg="(IonID=="+str(self.IonID)+") & (Decon>="+str(self.decon)+") & (Overlap=="+str(self.overlap)+")"
        if GPID!=None:
            self.GPID=GPID
            self.strtarg=self.strtarg+"& (GPID=="+str(self.GPID)+")"
            self.ionkey=str(self.GPID)+'_'+str(self.IonID)+'_'+str(self.overlap)
        else:
            self.ionkey='None_'+str(self.IonID)+'_'+str(self.overlap)
        if runID!=None:
            self.runID=runID
            self.strtarg=self.strtarg+"& (RunID==" +str(self.runID)+")"
            self.ionkey=str(self.runID)+'_'+self.ionkey
        else:
            self.ionkey='NoRunID_'+self.ionkey
        try:
            self.h5file=h5file
            self.h5xic=h5file.replace('.h5','_XIC.h5')
            self.h5connection=tb.open_file(self.h5file,mode='a')
            self.MS2=self.h5connection.root.MS2.MS2
            self.MS1=self.h5connection.root.MS1.MS1
        except:
            print("Please provide the hdf5file with .h5 extension")
            
    def IonGrab(self):
        intensity=[y['Intensity'] for y in self.MS2.where(self.strtarg)]
        time=[x['Time'] for x in self.MS2.where(self.strtarg)]
        preID=[z['PrecursorID'] for z in self.MS2.where(self.strtarg)]
        return time, intensity, preID
    
    def IonXICGen(self):
        time, inte, preid=self.IonGrab()
        self.XIC=pd.DataFrame(columns=['Time','Intensity','PrecursorID'])
        self.XIC['Time']=time
        self.XIC['Intensity']=inte
        self.XIC['PrecursorID']=preid
        if self.XIC.shape[0]>0:
            self.XIC.to_hdf(self.h5xic, self.ionkey,mode='a')
            
    
#GP XIC maker
class GPXIC(IonXIC):
    def __init__(self,GPID,MS1,decon=0,overlap=0,runID=None,h5file=None):
        super().__init__(GPID=GPID,IonID=None,decon=decon,overlap=overlap,runID=runID,h5file=h5file)
        self.strtarg="(GPID=="+str(self.GPID)+") & (Decon>="+str(self.decon)+") & (Overlap=="+str(self.overlap)+")"
        self.h5xic=h5file.replace('.h5','_XIC.h5')
        if runID!=None:
            self.runID=runID
            self.strtarg=self.strtarg+"& (RunID==" +str(self.runID)+")"
            self.gpkey=str(self.runID)+'_'+str(self.GPID)+'_GP_'+str(self.overlap)
        else:
            self.gpkey='NoRunID_'+str(self.GPID)+'_GP_'+str(self.overlap)
        
    def GPGrab(self):
        intensity=[y['Intensity'] for y in self.MS1.where(self.strtarg)]
        time=[x['Time'] for x in self.MS1.where(self.strtarg)]
        preID=[z['PrecursorID'] for z in self.MS1.where(self.strtarg)]
        return time, intensity, preID
    
    def GPXICGen(self):
        time, inte, preid=self.GPGrab()
        self.XIC=pd.DataFrame(columns=['Time','Intensity','PrecursorID'])
        self.XIC['Time']=time
        self.XIC['Intensity']=inte
        self.XIC['PrecursorID']=preid
        self.XIC.to_hdf(self.h5xic, self.gpkey,mode='w')    




#DO NOT USE THIS ONE IT IS TOO MANY ONLY MAKE TEMPORARY XICs
#make every XIC within an MS run
class XICwithinRun(GPXIC):
    def __init__(self,h5file,deconms1=0,deconms2=0,overlap1=0,overlap2=0,runID=None):
        self.runID=runID
        self.deconms1=deconms1
        self.deconms2=deconms2
        self.overlap1=overlap1
        self.overlap2=overlap2
        try:
            self.h5file=h5file
            self.h5xic=h5file.replace('.h5','_XIC.h5')
            self.h5connection=tb.open_file(self.h5file,mode='a')
            self.MS2=self.h5connection.root.MS2.MS2
            self.MS1=self.h5connection.root.MS1.MS1
        except:
            print("Please provide the hdf5file")
        
    def GPAssoc(self):
        self.dfGlyIonAssoc=pd.read_hdf(self.h5file,key='GlyIonAssoc',mode='a')
    def IonMeta(self):
        self.dfIonMeta=pd.read_hdf(self.h5file,key='IonMetaData',mode='a')
    def PreMeta(self):
        self.dfPreMeta=pd.read_hdf(self.h5file,key='PrecursorMetaData',mode='a')
    
    def IonXICMaker(self):
        self.IonMeta()
        IonTargets=self.dfIonMeta['IonID'].to_list()
        IonAUC=pd.DataFrame(None,index=IonTargets,columns=['AUC'])  
        for i, ix in enumerate(IonTargets):
            IonObj=IonXIC(ix,GPID=None,MS2=self.MS2,decon=self.deconms2,overlap=self.overlap2,runID=self.runID,h5file=self.h5file)
            IonObj.IonXICGen()
            if IonObj.XIC.shape[0]>1:
                IonAUC['AUC'].iloc[i]=AbundanceCalc.trap_sum_DF(IonObj.XIC,'Intensity')
                # assocGPtemp=[0]*IonObj.XIC.shape[0]
                # for j, jx in enumerate(IonObj.XIC['PrecursorID']):
                #     tempval=[x['GPID'] for x in self.MS1.where('(PrecursorID=='+str(jx)+')')]
                #     if len(tempval)==0:
                #         assocGPtemp[j]=0
                #     else:
                #         assocGPtemp[j]=tempval
                # GPIonAUC=pd.DataFrame(None,index=assocGPtemp,columns=['AUC'])
                # for jx in assocGPtemp:
                #     GPIonObj=IonXIC(ix,GPID=jx,MS2=self.MS2,decon=self.deconms2,overlap=self.overlap2,runID=self.runID,h5file=self.h5file)
                #     GPIonObj.IonXICGen()
                #     if GPIonObj.XIC.shape[0]>1:
                #         GPIonAUC['AUC'].iloc[i]=AbundanceCalc.trap_sum_DF(GPIonObj.XIC,'Intensity')
                # GPIonAUC.to_hdf(self.h5xic,key='Ion'+str(ix)+'_XIC_byGP_AUCs_of_'+str(self.runID),mode='a')
            else:
                IonAUC['AUC'].iloc[i]=0
        IonAUC.to_hdf(self.h5xic,key='Ion_WholeXIC_AUCs_of_'+str(self.runID),mode='a')
        
        
    def GPXICMaker(self):
        self.GPAssoc()
        GPTargets=self.dfGlyIonAssoc['GPID']
        GPAUC=pd.DataFrame(None,index=GPTargets,columns=['AUC'])
        #GP_Pre=pd.DataFrame(None,columns=['GPID','PrecursorID'])
        for i, ix in enumerate(GPTargets):
            GPObj=GPXIC(ix,MS1=self.MS1,decon=self.deconms1,overlap=self.overlap1,runID=self.runID,h5file=self.h5file)
            GPObj.GPXICGen()
            if GPObj.XIC.shape[0]>1:
                GPAUC['AUC'].iloc[i]=AbundanceCalc.trap_sum_DF(GPObj.XIC,'Intensity')[0]
            else:
                GPAUC['AUC'].iloc[i]=0
        GPAUC.to_hdf(self.h5xic,key='GP_WholeXIC_AUCs_of_'+str(self.runID),mode='a')
        #GPAUC['PrecursorID'].groupby('GPID')
        
    def main(self):
        self.GPXICMaker()
        self.IonXICMaker()
 
class XIC_all_runs(XICwithinRun):
    def __init__(self,h5file,deconms1=0,deconms2=0,overlap1=0,overlap2=0,runIDs=None):
        super().__init__(h5file,deconms1,deconms2,overlap1,overlap2)
        self.runIDs=runIDs
        
    def Loop(self):
        for rid in self.runIDs:
            self.runID=rid
            self.main()
            
    
class XICAssociator:
    def __init__(self,GPID,h5file):
        self.GPID=GPID
        try:
            self.h5file=h5file
        except:
            print("Please provide the hdf5file")
    
    def GPAssoc(self):
        self.dfGlyIonAssoc=pd.read_hdf(self.h5file,key='GlyIonAssoc',mode='a')
    def IonMeta(self):
        self.dfIonMeta=pd.read_hdf(self.h5file,key='IonMetaData',mode='a')
    def PreMeta(self):
        self.dfPreMeta=pd.read_hdf(self.h5file,key='PrecursorMetaData',mode='a')
        
    def GP_Ion_Links(self):
        self.GPAssoc()
        self.dfGlyIonAssoc.loc[self.dfGlyIonAssoc['GPID']==self.GPID]['IonID']
        
        
def AddIDStatistics(AddIDPeakTable,AddToGP,GPIonAssoc):
    AddIDTable=pd.DataFrame([[np.nan,0.0,0.0,0.0]],columns=['AddID','GPID','AcqIdx','apex'])
    for adx in AddIDPeakTable['AddID'].unique():
        for jdx in AddIDPeakTable.loc[AddIDPeakTable['AddID']==adx,'AcqIdx'].unique():
            if AddIDPeakTable.loc[(AddIDPeakTable['AddID']==adx)&(AddIDPeakTable['AcqIdx']==jdx)&(AddIDPeakTable['ApexDist']==0),'apex'].shape[0]>1:
                AddIDTable.loc[AddIDTable.index.max()+1]=[adx,AddToGP.loc[adx,'GPID'],jdx,AddIDPeakTable.loc[(AddIDPeakTable['AddID']==adx)&(AddIDPeakTable['AcqIdx']==jdx)&(AddIDPeakTable['ApexDist']==0),'apex'].iloc[0]]
            else:
                sub=AddIDPeakTable.loc[(AddIDPeakTable['AddID']==adx)&(AddIDPeakTable['AcqIdx']==jdx)]
                AddIDTable.loc[AddIDTable.index.max()+1]=[adx,AddToGP.loc[adx,'GPID'],jdx,sub.loc[sub['ApexDist']==sub['ApexDist'].min(),'apex'].iloc[0]]
    AddIDTable=AddIDTable.dropna().copy()
    AddIDTable['Wavelen']=0.0
    AddIDTable['TotalSignal']=0.0
    AddIDTable['UniqueSignal']=0.0
    AddIDTable['PotUniSignal']=0.0
    AddIDTable['N_IonsTotal']=0.0
    AddIDTable['N_IonsUni']=0.0
    AddIDTable['N_IonsPotUni']=0.0
    AddIDTable['RunID']=AddIDPeakTable['RunID'].iloc[0]
    for adx in AddIDTable['AddID'].index:
        tempapex=AddIDTable.loc[adx,'apex']
        tacq=AddIDTable.loc[adx,'AcqIdx']
        AddIDTable.loc[adx,'Wavelen']=AddIDPeakTable.loc[(AddIDPeakTable['AddID']==AddIDTable.loc[adx,'AddID'])&(AddIDPeakTable['apex']==tempapex),'wavelen'].max()
        potcompadds=list(set(AddIDTable.loc[(AddIDTable['AcqIdx']==tacq)&(AddIDTable['apex']>=(tempapex-12))&(AddIDTable['apex']<=(tempapex+12)),'AddID'].unique())-set([AddIDTable.loc[adx,'AddID']]))
        if len(potcompadds)>0:
            competingions=list(set(GPIonAssoc.loc[AddToGP.loc[potcompadds,'GPID'],'IonID'].sum()))
        else:
            competingions=list(set())
        sub1=AddIDPeakTable.loc[(AddIDPeakTable['AddID']==AddIDTable.loc[adx,'AddID'])&(AddIDPeakTable['AcqIdx']==tacq)]
        AddIDTable.loc[adx,'N_IonsTotal']=sub1['IonID'].nunique()
        AddIDTable.loc[adx,'TotalSignal']=sub1['AUC'].sum()
        defuniqueions=list(set(sub1['IonID'].tolist())-set(competingions))
        AddIDTable.loc[adx,'N_IonsUni']=len(defuniqueions)
        if len(defuniqueions)>0:
            AddIDTable.loc[adx,'UniqueSignal']=sub1.loc[sub1['IonID'].isin(defuniqueions),'AUC'].sum()
        seenions=set(AddIDPeakTable.loc[(AddIDPeakTable['AddID'].isin(potcompadds))&(AddIDPeakTable['AcqIdx']==tacq),'IonID'].tolist())
        possuniions=list(set(sub1['IonID'].tolist())-seenions)
        AddIDTable.loc[adx,'N_IonsPotUni']=len(possuniions)
        if len(possuniions)>0:
            AddIDTable.loc[adx,'PotUniSignal']=sub1.loc[sub1['IonID'].isin(possuniions),'AUC'].sum()
    AddIDTable['PercDefUnique']=np.log(AddIDTable['UniqueSignal']+1)/np.log(AddIDTable['TotalSignal']+1)
    AddIDTable['PercPotUnique']=np.log(AddIDTable['PotUniSignal']+1)/np.log(AddIDTable['TotalSignal']+1)
    return AddIDTable

def DifferentSummations(CurveDF,grptgt='AddID',sumtgt='AUC',summethod='sum'):
    TotalSum=CurveDF.groupby(grptgt).agg({sumtgt:summethod}).reset_index()
    TotalSum.rename({sumtgt:sumtgt+'_Broad'},inplace=True,axis=1)
    CloseSum=CurveDF.loc[CurveDF['ApexDist']<=1].groupby(grptgt).agg({sumtgt:summethod}).reset_index()
    CloseSum.rename({sumtgt:sumtgt+'_Close'},inplace=True,axis=1)
    ExactSum=CurveDF.loc[CurveDF['ApexDist']==0].groupby(grptgt).agg({sumtgt:summethod}).reset_index()
    ExactSum.rename({sumtgt:sumtgt+'_Exact'},inplace=True,axis=1)
    OneOffSum=CurveDF.loc[CurveDF['ApexDist']==1].groupby(grptgt).agg({sumtgt:summethod}).reset_index()
    OneOffSum.rename({sumtgt:sumtgt+'_OneOff'},inplace=True,axis=1)
    outdf=pd.merge(ExactSum,CloseSum,on=grptgt)
    outdf=pd.merge(outdf,TotalSum,on=grptgt)
    outdf=pd.merge(outdf,OneOffSum,on=grptgt)
    return outdf
    
    
    
