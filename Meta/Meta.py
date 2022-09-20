##GlyLine Meta Data Classes GlyLine.Meta
import pandas as pd
import numpy as np
from GlyLine.Helper import Helper

class DataTable:
    def __init__(self,hdf5file,key=None):
        self.h5file=hdf5file
        self.key=key

    
class IonMetaDataTable(DataTable):
#class IonMetaDataTable:
    # def __init__(self,hdf5file,filelist=None,index=['scan_id'],key='IonMetaData'):
    #     self.h5file=hdf5file
    #     self.files=filelist
    #     self.index=index
    #     self.key=key
    def __init__(self,hdf5file,filelist=None,index=['scan_id'],key='IonMetaData'):
        super().__init__(hdf5file,key)
        self.files=filelist
        self.index=index
        self.key=key
        
        
    def IonMetadata(self,dfIon):
        dfIonMeta=pd.DataFrame(data=None,columns=['neutralmass','glycopeptide','fragment_name','fragment_type'])
        mass, gps, ions=Helper.MassCalc(dfIon)
        dfIonMeta['glycopeptide']=gps
        dfIonMeta['fragment_name']=ions
        dfIonMeta['neutralmass']=mass
        Helper.FragmentType(dfIonMeta)
        return dfIonMeta
    
    def GlycanSimplifier(self,dfMeta):
        temp=dfMeta.where(dfMeta['fragment_type']=='Glycan').groupby(['fragment_name'])['IonID'].min().reset_index()
        for i, frag in enumerate(temp['fragment_name']):
            dfMeta.loc[dfMeta['fragment_name']==frag,'IonID']=temp['IonID'][i]
    
    def IonMasterListMaker(self):
        dfIonMetaMasterTemp=pd.DataFrame(data=None,columns=['neutralmass','glycopeptide','fragment_name','fragment_type'])
        for f in self.files:
            dfIonTemp=Helper.GSoftCSVRead(f,subset=['glycopeptide','fragment_name'],index=self.index,PSMBool=False)
            MetaTemp=self.IonMetadata(dfIonTemp)
            dfIonMetaMasterTemp=pd.concat([dfIonMetaMasterTemp,MetaTemp],ignore_index=True)
        
        ix=[i for i,x in enumerate(dfIonMetaMasterTemp.duplicated(['glycopeptide','fragment_name'])) if x==False]
        dfIonMetaMaster=dfIonMetaMasterTemp.loc[ix,]
        dfIonMetaMaster['IonID']=range(dfIonMetaMaster.shape[0])
        dfIonMetaMaster['IonID']=dfIonMetaMaster['IonID']
        self.GlycanSimplifier(dfIonMetaMaster)
        self.dfIonMetaMaster = dfIonMetaMaster
        
    def ListWriter(self):
        self.dfIonMetaMaster.to_hdf(self.h5file,key=self.key,mode='a')
    
    def ListReader(self):
        self.dfIonMetaMaster=pd.read_hdf(self.h5file, key=self.key)


        
### This is the class for the glycopeptide to ion association table
class GPIonAssociation(IonMetaDataTable):
    def __init__(self,hdf5file,filelistspec=None,GPIkey='GlyIonAssoc'):
        super().__init__(hdf5file,filelistspec)
        self.GPIkey=GPIkey
    
    def GP_Ion_Data(self):
        if hasattr(self,'dfIonMetaMaster')==False:
            self.ListReader()
        self.dfGPIonAssoc=self.dfIonMetaMaster.groupby(['glycopeptide'])['IonID'].apply(np.unique).reset_index()
        self.dfGPIonAssoc['GPID']=range(self.dfGPIonAssoc.shape[0])
        self.dfGPIonAssoc['GPID']=self.dfGPIonAssoc['GPID']+1
        
    def GPIonWriter(self):
        self.dfGPIonAssoc.to_hdf(self.h5file,key=self.GPIkey,mode='a')
    
    def GPIonReader(self):
        self.dfGPIonAssoc=pd.read_hdf(self.h5file,key=self.GPIkey)
        

#this is to gather the info on precursors
#which scans are confirmed in which samples, the times and intensities associated therein
#the meta data are stored like this
## RunID | protein | Glycopeptide | scan_id | time | neutralmass | mass_accuracy | charge | ms2_score | q_value
class PSMMetaDataTable(DataTable):
    def __init__(self,hdf5file,filelist=None,runidentifier=None,index=['scan_id'],PSMkey='PrecursorMetaData',GPIkey='GlyIonAssoc'):
        super().__init__(hdf5file)
        self.files=filelist
        self.runIDs=runidentifier
        self.index=index
        self.key=PSMkey
        self.GPIkey=GPIkey
        
        
    def PSMMetadata(self,dfPSM,runID):
        dfPSMMeta=pd.DataFrame(data=None,
                               columns=['RunID','protein_name','glycopeptide','scan_id',
                                        'scan_time','neutralmass','mass_accuracy','charge',
                                        'ms2_score','q_value'])
        dfPSMMeta['RunID']=[runID]*dfPSM.shape[0]
        dfPSMMeta['protein_name']=dfPSM['protein_name'].values
        dfPSMMeta['glycopeptide']=dfPSM['glycopeptide'].values
        dfPSMMeta['scan_id']=dfPSM.index.values
        dfPSMMeta['scan_time']=dfPSM['scan_time'].values
        dfPSMMeta['neutralmass']=dfPSM['neutral_mass'].values
        dfPSMMeta['mass_accuracy']=dfPSM['mass_accuracy'].values
        dfPSMMeta['charge']=dfPSM['charge'].values
        dfPSMMeta['ms2_score']=dfPSM['ms2_score'].values
        dfPSMMeta['q_value']=dfPSM['q_value'].values
        return dfPSMMeta
    
    def PSMMasterListMaker(self):
        dfPSMMetaMasterTemp=pd.DataFrame(data=None,
                               columns=['RunID','protein_name','glycopeptide','scan_id',
                                        'scan_time','neutralmass','mass_accuracy','charge',
                                        'ms2_score','q_value'])
        for f in range(len(self.files)):
            dfPSMTemp=Helper.GSoftCSVRead(self.files[f],subset=['glycopeptide'],index=self.index,PSMBool=True)
            MetaTemp=self.PSMMetadata(dfPSMTemp,self.runIDs[f])
            dfPSMMetaMasterTemp=pd.concat([dfPSMMetaMasterTemp,MetaTemp],ignore_index=True)
        
        ix=[i for i,x in enumerate(dfPSMMetaMasterTemp.duplicated(['glycopeptide','RunID'])) if x==False]
        dfPSMMetaMaster=dfPSMMetaMasterTemp.loc[ix,]
        if (hasattr(self,'dfGPIonAssoc')==False):
            self.ListReader(self.GPIkey)
        dl=dfPSMMetaMaster.shape[0]
        tempvec=[0]*dl
        for j in range(dl):
            tempval=self.dfGPIonAssoc.loc[self.dfGPIonAssoc['glycopeptide']==dfPSMMetaMaster['glycopeptide'].iloc[j],'GPID'].values
            if tempval.size>0:
                tempvec[j]=tempval
            else:
                tempvec[j]=np.max(self.dfGPIonAssoc['GPID'].values)+1
        dfPSMMetaMaster['GPID']=tempvec
        self.dfPSMMetaMaster = dfPSMMetaMaster
        
    def ListWriter(self):
        self.dfPSMMetaMaster.to_hdf(self.h5file,key=self.key,mode='a')
    
    def ListReader(self,key1=None):
        if key1==None:
            self.dfPSMMetaMaster=pd.read_hdf(self.h5file, self.key)
        else:
            self.dfGPIonAssoc=pd.read_hdf(self.h5file, key=key1)
        
