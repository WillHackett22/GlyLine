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
    
    def IonMasterListMaker(self):
        dfIonMetaMasterTemp=pd.DataFrame(data=None,columns=['neutralmass','glycopeptide','peptide','fragment_name','fragment_type'])
        for f in self.files:
            dfIonTemp=Helper.GSoftCSVRead(f,subset=['glycopeptide','fragment_name'],index=self.index,PSMBool=False)
            MetaTemp=self.IonMetadata(dfIonTemp)
            dfIonMetaMasterTemp=pd.concat([dfIonMetaMasterTemp,MetaTemp],ignore_index=True)
        
        #reduce to unique ions by glycopeptide
        ix=[i for i,x in enumerate(dfIonMetaMasterTemp.duplicated(['glycopeptide','fragment_name'])) if x==False]
        dfIonMetaMasterTemp=dfIonMetaMasterTemp.loc[ix,]
        # reduce to glycans
        GlyMeta=dfIonMetaMasterTemp.loc[dfIonMetaMasterTemp['fragment_type']=='Glycan']
        #produce IDs for unique glycnas
        GlyGroups=GlyMeta.groupby(['fragment_name'])['neutralmass'].apply(list).reset_index()
        GlyGroups['IonID']=range(GlyGroups.shape[0])
        GlyGroups=GlyGroups.drop(['neutralmass'],axis=1)
        #attach ids back to glycans by glycopep
        GlyMetaIDed=pd.merge(GlyGroups,GlyMeta,on='fragment_name')
        # reduce to peptide
        PepMeta=dfIonMetaMasterTemp.loc[dfIonMetaMasterTemp['fragment_type']=='Peptide']
        # give peptide ions ids
        PepGroups=PepMeta.groupby(['peptide','fragment_name'])['neutralmass'].apply(list).reset_index()
        PepGroups['IonID']=range(PepGroups.shape[0])
        PepGroups=PepGroups.drop(['neutralmass'],axis=1)
        # make pep ions start after gly
        PepGroups['IonID']=PepGroups['IonID']+GlyGroups.shape[0]
        # attach back to peptides by glycopep
        PepMetaIDed=pd.merge(PepGroups,PepMeta,on=['peptide','fragment_name'])
        #isolate stub
        StubMeta=dfIonMetaMasterTemp.loc[dfIonMetaMasterTemp['fragment_type']=='Stub']
        #stub ion ids
        StubGroups=StubMeta.groupby(['peptide','fragment_name'])['neutralmass'].apply(list).reset_index()
        StubGroups['IonID']=range(StubGroups.shape[0])
        StubGroups=StubGroups.drop(['neutralmass'],axis=1)
        # start stub ion ids after pep and gly
        StubGroups['IonID']=StubGroups['IonID']+GlyGroups.shape[0]+PepGroups.shape[0]
        # reattach stubs
        StubMetaIDed=pd.merge(StubGroups,StubMeta,on=['peptide','fragment_name'])
        # bring it all together
        dfIonMetaMaster=pd.concat([StubMetaIDed,PepMetaIDed,GlyMetaIDed])
        #getting the unique values for target list
        unidx=[i for i,x in enumerate(dfIonMetaMaster.duplicated(['IonID'])) if x==False]
        dfIonMetaUnique=dfIonMetaMaster.iloc[unidx,].reset_index()
        dfIonMetaUnique.loc[dfIonMetaUnique['fragment_type']=="Glycan",'peptide']='Oxonium'
        dfIonMetaUnique=dfIonMetaUnique.drop(['glycopeptide','index'],axis=1)
        dfIonMetaUnique=dfIonMetaUnique.set_index('IonID')
        self.dfIonMetaMaster = dfIonMetaMaster
        self.dfIonMetaUnique = dfIonMetaUnique
        
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
        dfGPIonAssoc=self.dfIonMetaMaster.groupby(['glycopeptide'])['IonID'].apply(list).reset_index()
        dfGPIonAssoc['GPID']=range(dfGPIonAssoc.shape[0])
        dfGPIonAssoc.set_index('GPID')
        self.dfGPIonAssoc=dfGPIonAssoc
        # let's get the inverse association table
        dfGPIonAssocMod=dfGPIonAssoc.drop(['IonID'],axis=1)
        dfTemp=pd.merge(self.dfIonMetaMaster,dfGPIonAssocMod,on='glycopeptide')
        dfIonGPAssoc=dfTemp.groupby(['IonID'])['GPID'].apply(list).reset_index()
        self.dfIonGPAssoc=dfIonGPAssoc
        
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
            tempval=self.dfGPIonAssoc.loc[self.dfGPIonAssoc['glycopeptide']==dfPSMMetaMaster['glycopeptide'].iloc[j],'GPID'].values.tolist()
            if len(tempval)==1:
                tempvec[j]=tempval[0]
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
        
