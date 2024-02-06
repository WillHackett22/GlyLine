##GlyLine Meta Data Classes GlyLine.Meta
import pandas as pd
import numpy as np
import re
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
    def __init__(self,hdf5file,filelistspec=None,index=['scan_id'],key='IonMetaData'):
        super().__init__(hdf5file,key)
        self.specfiles=filelistspec
        self.index=index
        self.key=key
        self.unikey=key+'Unique'
        
        
    def IonMetadata(self,dfIon):
        dfIonMeta=pd.DataFrame(data=None,columns=['neutralmass','glycopeptide','fragment_name','fragment_type'])
        dfIonMeta=Helper.MassCalc(dfIon)
        Helper.FragmentType(dfIonMeta)
        Helper.PeptideReducer(dfIonMeta)
        return dfIonMeta
    
    def IonMasterListMaker(self):
        dfIonMetaMasterTemp=pd.DataFrame(data=[[np.nan,'str','str','str','str']],columns=['neutralmass','glycopeptide','peptide','fragment_name','fragment_type'])
        for f in self.specfiles:
            dfIonTemp=Helper.GSoftCSVRead(f,subset=['glycopeptide','fragment_name'],index=self.index,PSMBool=False)
            MetaTemp=self.IonMetadata(dfIonTemp)
            dfIonMetaMasterTemp=pd.concat([dfIonMetaMasterTemp,MetaTemp],ignore_index=True)
        dfIonMetaMasterTemp=dfIonMetaMasterTemp.dropna()
        #reduce to unique ions by glycopeptide
        ix=[i for i,x in enumerate(dfIonMetaMasterTemp.duplicated(['glycopeptide','fragment_name'])) if x==False]
        dfIonMetaMasterTemp=dfIonMetaMasterTemp.iloc[ix]
        # reduce to glycans
        GlyMeta=dfIonMetaMasterTemp.loc[dfIonMetaMasterTemp['fragment_type']=='Glycan']
        #produce IDs for unique glycans
        GlyGroups=GlyMeta.groupby(['fragment_name'])['neutralmass'].apply(list).reset_index()
        GlyGroups['IonID']=range(GlyGroups.shape[0])
        GlyGroups=GlyGroups.drop(['neutralmass'],axis=1)
        #attach ids back to glycans by glycopep
        GlyMetaIDed=pd.merge(GlyGroups,GlyMeta,on='fragment_name')
        GlyMetaIDed['new_fragment_name']=GlyMetaIDed['fragment_name']
        # reduce to peptide backbone units
        PepMeta=dfIonMetaMasterTemp.loc[dfIonMetaMasterTemp['fragment_type']=='Peptide'].copy()
        PepMeta['new_fragment_name']='str'
        for idx in PepMeta.index:
            backbone=re.sub("\(.*?\)","",PepMeta.loc[idx,'peptide'])
            if PepMeta.loc[idx,'fragment_name']=='peptide':
                newfragname=backbone
            else:
                fraglength=list(map(int,re.findall(r'\d+', PepMeta.loc[idx,'fragment_name'])))[0]
                if 'b' in PepMeta.loc[idx,'fragment_name']:
                    newfragname=backbone[:fraglength]
                else:
                    newfragname=backbone[-fraglength:]
            PepMeta.loc[idx,'new_fragment_name']=newfragname
        # give peptide ions ids
        PepGroups=PepMeta.groupby(['neutralmass','new_fragment_name'])['fragment_name'].apply(list).reset_index()
        PepGroups['IonID']=range(PepGroups.shape[0])
        PepGroups=PepGroups.drop(['fragment_name'],axis=1)
        # make pep ions start after gly
        PepGroups['IonID']=PepGroups['IonID']+GlyGroups.shape[0]
        # attach back to peptides by glycopep
        PepMetaIDed=pd.merge(PepGroups,PepMeta,on=['neutralmass','new_fragment_name'])
        #stub ion ids
        StubMeta=dfIonMetaMasterTemp.loc[dfIonMetaMasterTemp['fragment_type']=='Stub'].copy()        
        StubMeta['new_fragment_name']='str'
        for idx in StubMeta.index:
            backbone=re.sub("\(.*?\)","",StubMeta.loc[idx,'peptide'])
            fragname=StubMeta.loc[idx,'fragment_name'].split('+')
            if fragname[0]=='peptide':
                newfragname=backbone+'+'+fragname[1]
            else:
                fraglength=list(map(int,re.findall(r'\d+', fragname[0])))[0]
                if 'b' in fragname[0]:
                    newfragname=backbone[:fraglength]+'+'+fragname[1]
                else:
                    newfragname=backbone[-fraglength:]+'+'+fragname[1]
            StubMeta.loc[idx,'new_fragment_name']=newfragname
        StubGroups=StubMeta.groupby(['neutralmass','new_fragment_name'])['fragment_name'].apply(list).reset_index()
        StubGroups['IonID']=range(StubGroups.shape[0])
        StubGroups=StubGroups.drop(['fragment_name'],axis=1)
        # start stub ion ids after pep and gly
        StubGroups['IonID']=StubGroups['IonID']+GlyGroups.shape[0]+PepGroups.shape[0]
        # reattach stubs
        StubMetaIDed=pd.merge(StubGroups,StubMeta,on=['neutralmass','new_fragment_name'])
        # bring it all together
        dfIonMetaMaster=pd.concat([StubMetaIDed,PepMetaIDed,GlyMetaIDed])
        dfIonMetaMaster=dfIonMetaMaster.reset_index()
        dfIonMetaMaster=dfIonMetaMaster.drop(['index'],axis=1)
        #getting the unique values for target list
        unidx=[i for i,x in enumerate(dfIonMetaMaster.duplicated(['IonID'])) if x==False]
        dfIonMetaUnique=dfIonMetaMaster.iloc[unidx,].reset_index()
        dfIonMetaUnique.loc[dfIonMetaUnique['fragment_type']=="Glycan",'peptide']='Oxonium'
        dfIonMetaUnique=dfIonMetaUnique.set_index('IonID')
        dfIonMetaUnique=dfIonMetaUnique.drop(['glycopeptide','index'],axis=1)
        self.dfIonMetaMaster = dfIonMetaMaster
        self.dfIonMetaUnique = dfIonMetaUnique
        
    def ListWriter(self):
        self.dfIonMetaMaster.to_hdf(self.h5file,key=self.key,mode='a')
        self.dfIonMetaUnique.to_hdf(self.h5file,key=self.unikey,mode='a')
    
    def ListReader(self):
        self.dfIonMetaMaster=pd.read_hdf(self.h5file, key=self.key)
        self.dfIonMetaUnique=pd.read_hdf(self.h5file,key=self.unikey)
        

class OverlappingInformation:
    def __init__(self,gpiobj,coreglys=None,stubcoreglys=None):
        if coreglys is None:
            self.coreglys=['HexNAc','HexNAcHexNAc','HexHexNAcHexNAc','HexHexHexNAcHexNAc','HexHexHexHexNAcHexNAc','Hex','HexHex','HexHexHex']
        else:
            self.coreglys=coreglys
        if stubcoreglys is None:
            self.stubcoreglys=['HexNAc','HexNAc1','HexNAc2','2HexNAc','Hex1HexNAc2','Hex2HexNAc2','Hex3HexNAc2']
        else:
            self.stubcoreglys=stubcoreglys
        self.dfIon=gpiobj.dfIonMetaMaster
        self.dfIonUnique=gpiobj.dfIonMetaUnique
        self.dfGPIonAssoc=gpiobj.dfGPIonAssoc
    
    def main(self):
        #get groups of oxonium ions
        glycangroups=self.GlycanIonGrouper()
        #check which of those are in the universal glycans
        presentcore=[cgly for cgly in self.coreglys if cgly in glycangroups.index.values]
        #get all the peptide backbones
        unipeps=self.dfIon['peptide'].unique().tolist()
        #get associations between backbones
        pepassocy=pd.DataFrame(False,columns=[unipeps],index=unipeps)
        pepassocb=pd.DataFrame(False,columns=[unipeps],index=unipeps)
        for pidx in unipeps:
            pepassocy.loc[pidx]=[pidx.endswith(pep) for pep in pepassocy.columns.tolist()]
            pepassocb.loc[pidx]=[pidx.startswith(pep) for pep in pepassocb.columns.tolist()]
        Y_Assoc={}
        B_Assoc={}
        for pidx in unipeps:
            tempy=list(set(pepassocy.index[np.where(pepassocy[pidx])[0].tolist()].tolist())|set(pepassocy.columns[pepassocy.loc[pidx]].get_level_values(0).tolist()))
            tempb=list(set(pepassocb.index[np.where(pepassocb[pidx])[0].tolist()].tolist())|set(pepassocb.columns[pepassocb.loc[pidx]].get_level_values(0).tolist()))
            Y_Assoc.update({pidx:tempy})
            B_Assoc.update({pidx:tempb})
            
        for pep in unipeps:
            #for each peptide get pep fragments and core stub ions
            unipepset, unipepfrags=self.UniquePeptideGrouper(pep,Y_Assoc,B_Assoc)
            unistubset,unistubfrags=self.UniqueStubGrouper(pep,Y_Assoc,B_Assoc)
            #get all gps for a pep backbone
            subgps=self.dfIon.loc[self.dfIon['peptide']==pep,'glycopeptide'].unique().tolist()
            for gp in subgps:
                # for each gp get difference between existing ions and implied ions
                existpepset=self.ExistingSubset(gp,'Peptide')
                missingfrags=list(unipepset-existpepset)
                existstubset=self.ExistingSubset(gp,'Stub')
                missingfrags=missingfrags+list(unistubset-existstubset)
                gfs=self.dfIon.loc[(self.dfIon['glycopeptide']==gp) & (self.dfIon['fragment_type']=='Glycan'),'new_fragment_name'].unique().tolist()    
                existglyset=self.ExistingSubset(gp,'Glycan')
                gfsimp=np.unique([re.sub("-.*","",s) for s in gfs]+presentcore).tolist()
                # add fuc and neu5ac to pertinent gps
                if 'Fuc' in gp:
                    gfsimp=gfsimp+['Fuc']
                if 'Neu5Ac' in gp:
                    gfsimp=gfsimp+['Neu5Ac']
                gidlist=[]
                for g in gfsimp:
                    gidlist=gidlist+np.unique(glycangroups.loc[g].tolist()[0]).tolist()
                missingfrags=missingfrags+list(set(gidlist)-existglyset)
                #update dfIon aka dfIonMetaMaster
                if len(missingfrags)>0:
                    self.IonFragmentAdder(missingfrags,gp,pep)
        self.GP_Ion_Data()
       
    def GlycanIonGrouper(self):
        #this function collects oxonium ions and groups them by base glycan composition (eg ignores the -h2o etc)
        glycanfrags=self.dfIon.loc[self.dfIon['fragment_type']=='Glycan','fragment_name'].tolist()
        glycantranslator=pd.DataFrame(None,columns=['fragment_name','simplified','IonID'])
        glycantranslator['fragment_name']=glycanfrags
        glycantranslator['simplified']=[re.sub("-.*","",s) for s in glycanfrags]
        glycantranslator['IonID']=self.dfIon.loc[self.dfIon['fragment_type']=='Glycan','IonID'].tolist()
        glycangroups=glycantranslator.groupby('simplified')['IonID'].apply(lambda x: list(set(x))).reset_index()
        glycangroups=glycangroups.set_index('simplified')
        return glycangroups
    
    def IonFragmentAdder(self,missingfrags,gp,pep):
        #this function creates a dataframe of missingfrags, replaces the pertinent info, and adds it onto the ms2 unique ion list
        tempdf=self.dfIonUnique.loc[missingfrags].copy()
        tempdf['IonID']=tempdf.index.values.tolist()
        tempdf['glycopeptide']=gp
        tempdf['peptide']=pep
        self.dfIon=pd.concat([self.dfIon,tempdf])
        self.dfIon=self.dfIon.reset_index()
        self.dfIon=self.dfIon.drop(['index'],axis=1)
            
    def UniquePeptideGrouper(self,pep,Y_Assoc,B_Assoc):
        #get unique peptide frags by peptide backbone, give to all members of peptide backbone
        unipepids_y=self.dfIon.loc[(self.dfIon['fragment_name'].str[0]=='y')&(self.dfIon['peptide'].isin(Y_Assoc[pep])) & (self.dfIon['fragment_type']=='Peptide'),'IonID'].unique().tolist()
        unipepids_b=self.dfIon.loc[(self.dfIon['fragment_name'].str[0]=='b')&(self.dfIon['peptide'].isin(B_Assoc[pep])) & (self.dfIon['fragment_type']=='Peptide'),'IonID'].unique().tolist()
        unipepids=unipepids_y+unipepids_b
        unipepfrags=self.dfIonUnique.loc[unipepids,'new_fragment_name'].tolist()
        return set(unipepids), unipepfrags

    def UniqueStubGrouper(self,pep,Y_Assoc,B_Assoc):
        # get stub ions part of core, gve to all stub ions
        stubids_y=self.dfIon.loc[(self.dfIon['fragment_name'].str[0]=='y')&(self.dfIon['fragment_type']=='Stub') & (self.dfIon['peptide'].isin(Y_Assoc[pep])),'IonID'].unique().tolist()
        stubids_b=self.dfIon.loc[(self.dfIon['fragment_name'].str[0]=='b')&(self.dfIon['fragment_type']=='Stub') & (self.dfIon['peptide'].isin(B_Assoc[pep])),'IonID'].unique().tolist()
        stubids=stubids_y+stubids_b
        stubfrags_y=self.dfIon.loc[(self.dfIon['fragment_name'].str[0]=='y')&(self.dfIon['fragment_type']=='Stub') & (self.dfIon['peptide'].isin(Y_Assoc[pep])),'new_fragment_name'].unique().tolist()
        stubfrags_b=self.dfIon.loc[(self.dfIon['fragment_name'].str[0]=='b')&(self.dfIon['fragment_type']=='Stub') & (self.dfIon['peptide'].isin(B_Assoc[pep])),'new_fragment_name'].unique().tolist()
        stubfrags=stubfrags_y+stubfrags_b
        keepidx=[i for i,s in enumerate(stubfrags) if re.sub(".*\\+","",s) in self.stubcoreglys]
        stubfragid=[stubids[i] for i in keepidx]
        unistubfrags=[stubfrags[i] for i in keepidx]
        return set(stubfragid), np.unique(unistubfrags).tolist()
    
    def ExistingSubset(self,gp,targetiontype):
        # get existing subset of an ion type for a glycoeptide
        return set(self.dfIon.loc[(self.dfIon['glycopeptide']==gp) & (self.dfIon['fragment_type']==targetiontype),'IonID'])
    
    def GP_Ion_Data(self):
        gpiddict={self.dfGPIonAssoc.loc[u,'glycopeptide']:u for u in self.dfGPIonAssoc.index.tolist()}
        dfGPIonAssocNew=self.dfIon.groupby(['glycopeptide'])['IonID'].apply(list).reset_index()
        dfGPIonAssocNew['GPID']=[gpiddict[gp] for gp in dfGPIonAssocNew['glycopeptide'].tolist()]
        dfGPIonAssocNew=dfGPIonAssocNew.set_index('GPID')
        dfGPIonAssocNew['GPID']=dfGPIonAssocNew.index.values.tolist()
        self.dfGPIonAssoc=dfGPIonAssocNew
        # let's get the inverse association table
        dfGPIonAssocMod=dfGPIonAssocNew.drop(['IonID'],axis=1)
        dfTemp=pd.merge(self.dfIon,dfGPIonAssocMod,on='glycopeptide')
        dfIonGPAssocNew=dfTemp.groupby(['IonID'])['GPID'].apply(list).reset_index()
        self.dfIonGPAssoc=dfIonGPAssocNew
    
### This is the class for the glycopeptide to ion association table and viceversa
class GPIonAssociation(IonMetaDataTable):
    def __init__(self,hdf5file,filelistspec=None,GPIkey='GlyIonAssoc',IGPkey='IonGlyAssoc'):
        super().__init__(hdf5file,filelistspec)
        self.GPIkey=GPIkey
        self.IGPkey=IGPkey
        
    def GP_Ion_Data(self):
        if hasattr(self,'dfIonMetaMaster')==False:
            self.ListReader()
        dfGPIonAssoc=self.dfIonMetaMaster.groupby(['glycopeptide'])['IonID'].apply(list).reset_index()
        dfGPIonAssoc['GPID']=[r+1 for r in range(dfGPIonAssoc.shape[0])]
        #add 1 so that GPID 0 is miscellaneous
        dfGPIonAssoc=dfGPIonAssoc.set_index('GPID')
        self.dfGPIonAssoc=dfGPIonAssoc
        dfGPIonAssoc['GPID']=dfGPIonAssoc.index.values.tolist()
        # let's get the inverse association table
        dfGPIonAssocMod=dfGPIonAssoc.drop(['IonID'],axis=1)
        dfTemp=pd.merge(self.dfIonMetaMaster,dfGPIonAssocMod,on='glycopeptide')
        dfIonGPAssoc=dfTemp.groupby(['IonID'])['GPID'].apply(list).reset_index()
        self.dfIonGPAssoc=dfIonGPAssoc
    
    def IonGPWriter(self):
        self.dfIonGPAssoc.to_hdf(self.h5file,key=self.IGPkey,mode='a')
    
    def GPIonWriter(self):
        self.dfGPIonAssoc.to_hdf(self.h5file,key=self.GPIkey,mode='a')
        
    def AssocWriter(self):
        self.IonGPWriter()
        self.GPIonWriter()
    
    def GPIonReader(self):
        self.dfGPIonAssoc=pd.read_hdf(self.h5file,key=self.GPIkey)
    
    def IonGPReader(self):
        self.dfIonGPAssoc=pd.read_hdf(self.h5file,key=self.IGPkey)
        
    def AssocReader(self):
        self.GPIonReader()
        self.IonGPReader()
    
class MetaIonDataMainFunction(DataTable):
    def __init__(self,hdf5file,filelistspec=None,runidentifier=None,index=['scan_id'],key='IonMetaData',
                 coreglys=None,stubcoreglys=None,GPIkey='GlyIonAssoc',IGPkey='IonGlyAssoc'):
        super().__init__(hdf5file,key)
        self.specfiles=filelistspec
        if runidentifier==None:
            self.runIDs=['r'+str(x) for x in range(len(self.specfiles))]
        self.index=index
        self.key=key
        self.unikey=key+'Unique'
        self.coreglys=coreglys
        self.stubcoreglys=stubcoreglys
        self.GPIkey=GPIkey
        self.IGPkey=IGPkey
        
    def MainGenerate(self):
        self.gpiobj=GPIonAssociation(self.h5file,self.specfiles)
        self.gpiobj.IonMasterListMaker()
        self.gpiobj.ListWriter()
        self.gpiobj.GP_Ion_Data()
        self.oi=OverlappingInformation(self.gpiobj,self.coreglys,self.stubcoreglys)
        self.oi.main()
        self.gpiobj.dfGPIonAssoc=self.oi.dfGPIonAssoc
        self.gpiobj.dfIonGPAssoc=self.oi.dfIonGPAssoc
        self.gpiobj.AssocWriter()
        
    def MainRead(self):
        self.gpiobj=GPIonAssociation(self.h5file,self.specfiles)
        self.gpiobj.ListReader()
        self.gpiobj.AssocReader()
        
#this is to gather the info on precursors
#which scans are confirmed in which samples, the times and intensities associated therein
#the meta data are stored like this
## RunID | protein | Glycopeptide | scan_id | time | neutralmass | mass_accuracy | charge | ms2_score | q_value
class PSMMetaDataTable(DataTable):
    def __init__(self,hdf5file,filelist=None,runidentifier=None,index=['scan_id'],PSMkey='PrecursorMetaData',GPIkey='GlyIonAssoc',adductdict=None,massdict=None,adductmultiplicity=1):
        super().__init__(hdf5file)
        self.files=filelist
        if (runidentifier==None) & (filelist is not None):
            self.runIDs=['r'+str(x) for x in range(len(self.files))]
        else:
            self.runIDs=runidentifier
        self.index=index
        self.key=PSMkey
        self.unikey=PSMkey+'Unique'
        self.addkey=PSMkey+'UniqueAdducted'
        self.GPIkey=GPIkey
        #adductdict has format: adductdict={'adductname':{'H':1,'N':1,'C':1,'O':1,'Na':1,'S':1,'K':1}} 
        # where adductname has chemical formula HNCONaSK
        #eg: adductdict={'Ammonium':{'H':4,'N':1},'sodium':{'Na':1},'waterloss':{'H':-2,'O':-1}}
        self.adductdict=adductdict
        if massdict is None:
            massdict={'H':1.00797,'O':15.9994,'N':14.0067,'C':12,'Na':22.98977,'K':39.0983,'S':32.05}
        self.massdict=massdict
        self.adductmultiplicity=adductmultiplicity
        
    def PSMMetadata(self,dfPSM,runID):
        dfPSMMeta=pd.DataFrame(data=None,
                               columns=['RunID','protein_name','glycopeptide','adduction','scan_id',
                                        'scan_time','neutralmass','mass_accuracy','charge',
                                        'ms2_score','q_value','precursor_intensity'])
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
        dfPSMMeta['adduction']=dfPSM['mass_shift_name'].values
        dfPSMMeta['precursor_intensity']=dfPSM['precursor_abundance'].values
        return dfPSMMeta
    
    def PSMMasterListMaker(self):
        dfPSMMetaMasterTemp=pd.DataFrame(data=[['str','str','str','str','str',np.nan,np.nan,np.nan,0,np.nan,np.nan,np.nan]],
                               columns=['RunID','protein_name','glycopeptide','adduction','scan_id',
                                        'scan_time','neutralmass','mass_accuracy','charge',
                                        'ms2_score','q_value','precursor_intensity'])
        for f in range(len(self.files)):
            dfPSMTemp=Helper.GSoftCSVRead(self.files[f],subset=['glycopeptide'],index=self.index,PSMBool=True)
            dfPSMTemp=dfPSMTemp.loc[dfPSMTemp['is_best_match']==True,]
            MetaTemp=self.PSMMetadata(dfPSMTemp,self.runIDs[f])
            dfPSMMetaMasterTemp=pd.concat([dfPSMMetaMasterTemp,MetaTemp],ignore_index=True)
        dfPSMMetaMasterTemp=dfPSMMetaMasterTemp.dropna()

        #ix=[i for i,x in enumerate(dfPSMMetaMasterTemp.duplicated(['glycopeptide','RunID'])) if x==False]
        #dfPSMMetaMaster=dfPSMMetaMasterTemp.loc[ix,]
        dfPSMMetaMasterTemp['obs_neutralmass']=dfPSMMetaMasterTemp['neutralmass'].tolist()
        dfPSMMetaMasterTemp['neutralmass']=dfPSMMetaMasterTemp['neutralmass']+dfPSMMetaMasterTemp['neutralmass']*dfPSMMetaMasterTemp['mass_accuracy']
        dfPSMMetaMaster=dfPSMMetaMasterTemp.loc[dfPSMMetaMasterTemp['adduction']=='Unmodified'].copy()
        if self.adductdict!=None:
            dfAdducted=dfPSMMetaMasterTemp.loc[dfPSMMetaMasterTemp['adduction'].isin(list(self.adductdict.keys()))].copy()
            AdductsSeen=dfPSMMetaMasterTemp['adduction'].unique()
            ValidTgts=list(set(self.adductdict.keys())&set(AdductsSeen))
            for Add in ValidTgts:
                mshift=0
                for comp in self.adductdict[Add]:
                    mshift+=self.massdict[comp]*self.adductdict[Add][comp]
                dfAdducted.loc[dfAdducted['adduction']==Add,'neutralmass']=dfAdducted.loc[dfAdducted['adduction']==Add,'neutralmass']-mshift
            dfPSMMetaMaster=pd.concat([dfPSMMetaMaster,dfAdducted],ignore_index=True)
        if (hasattr(self,'dfGPIonAssoc')==False):
            self.ListReader(self.GPIkey)
        dl=dfPSMMetaMaster.shape[0]
        tempvec=[0]*dl
        for j in range(dl):
            tempval=self.dfGPIonAssoc.loc[self.dfGPIonAssoc['glycopeptide']==dfPSMMetaMaster['glycopeptide'].iloc[j],'GPID'].values.tolist()
            if len(tempval)==1:
                tempvec[j]=tempval[0]
            else:
                newgpid=np.max(self.dfGPIonAssoc['GPID'].values)+1
                tempvec[j]=newgpid
                self.dfGPIonAssoc.loc[self.dfGPIonAssoc.index.max()+1]=[dfPSMMetaMaster['glycopeptide'].iloc[j],[],newgpid]
        dfPSMMetaMaster['GPID']=tempvec
        self.dfPSMMetaMaster = dfPSMMetaMaster
        unidx=[i for i,x in enumerate(dfPSMMetaMaster.duplicated(['GPID'])) if x==False ]
        self.dfMS1Unique=pd.DataFrame(None,columns=['glycopeptide','neutralmass','GPID'])
        self.dfMS1Unique['glycopeptide']=dfPSMMetaMaster['glycopeptide'].iloc[unidx].tolist()
        self.dfMS1Unique['neutralmass']=dfPSMMetaMaster['neutralmass'].iloc[unidx].tolist()
        self.dfMS1Unique['neutralmass']=[np.round(nm,5) for nm in self.dfMS1Unique['neutralmass']]
        self.dfMS1Unique['GPID']=dfPSMMetaMaster['GPID'].iloc[unidx].tolist()
        self.dfMS1Unique=self.dfMS1Unique.set_index('GPID')
        self.dfPSMMetaMaster['AddID']=self.dfPSMMetaMaster['GPID']
        if self.adductdict is not None:
            self.dfMS1Adducted=self.AdductionListMaker()
            for pdx in self.dfPSMMetaMaster.index:
                self.dfPSMMetaMaster.loc[pdx,'AddID']=self.dfMS1Adducted.loc[(self.dfMS1Adducted['GPID']==self.dfPSMMetaMaster.loc[pdx,'GPID'])&(self.dfMS1Adducted['adducts']==self.dfPSMMetaMaster.loc[pdx,'adduction'])].index.values[0]
        else:
            self.dfMS1Adducted=self.dfMS1Unique.copy()
            self.dfMS1Adducted['GPID']=self.dfMS1Adducted.index.tolist()
            self.dfMS1Adducted['AddID']=self.dfMS1Adducted.index.tolist()
            self.dfMS1Adducted=self.dfMS1Adducted.set_index('AddID')
            self.dfMS1Adducted['adducts']=['Unmodified']*self.dfMS1Adducted.shape[0]
            
        
    def AdductAdder(self,adduct,basedf):
        tempdf=pd.DataFrame(None,columns=['GPID','neutralmass','adducts','AddID'])
        if basedf['adducts'].iloc[0]=='Unmodified':
            tempdf['adducts']=[adduct]*basedf.shape[0]
        else:
            tempdf['adducts']=basedf['adducts']+','+adduct
        tempdf['GPID']=basedf['GPID'].tolist()
        mshift=0
        for j in self.adductdict[adduct]:
            mshift+=self.adductdict[adduct][j]*self.massdict[j]
        tempdf['neutralmass']=[tmass+mshift for tmass in basedf['neutralmass'].tolist()]
        return tempdf
    
    def AdductionListMaker(self):
        targetsbase=self.dfMS1Unique
        adducteddf=pd.DataFrame(None,columns=['GPID','neutralmass','adducts','AddID'])
        adducteddf['GPID']=targetsbase.index.tolist()
        adducteddf['neutralmass']=targetsbase['neutralmass'].tolist()
        adducteddf['adducts']=['Unmodified']*targetsbase.shape[0]
        basedf=adducteddf.copy()
        adducts=list(self.adductdict)
        for idx,ad in enumerate(adducts):
            tempdf=self.AdductAdder(ad,basedf)
            adducteddf=pd.concat([adducteddf,tempdf])
            if self.adductmultiplicity>=2:
                combo=list(range(idx,len(adducts)))
                if len(combo)>0:
                    for c in combo:
                        adc=adducts[c]
                        subdf=self.AdductAdder(adc,tempdf)
                        adducteddf=pd.concat([adducteddf,subdf])
                        if self.adductmultiplicity>2:
                           combotri=list(range(c,len(adducts)))
                           if len(combotri)>0:
                               for t in combotri:
                                   adc=adducts[t]
                                   subsubdf=self.AdductAdder(adc,subdf)
                                   adducteddf=pd.concat([adducteddf,subsubdf])
        adducteddf['AddID']=[r+1 for r in range(adducteddf.shape[0])]
        adducteddf=adducteddf.set_index(['AddID'])
        return adducteddf
        
    def ListWriter(self):
        self.dfPSMMetaMaster.to_hdf(self.h5file,key=self.key,mode='a')
        self.dfMS1Unique.to_hdf(self.h5file,key=self.unikey,mode='a')
        self.dfMS1Adducted.to_hdf(self.h5file,key=self.addkey,mode='a')
    
    def ListReader(self,key1=None):
        if key1==None:
            self.dfPSMMetaMaster=pd.read_hdf(self.h5file, self.key)
            self.dfMS1Unique=pd.read_hdf(self.h5file, self.unikey)
            self.dfMS1Adducted=pd.read_hdf(self.h5file, self.addkey)
        else:
            self.dfGPIonAssoc=pd.read_hdf(self.h5file, key=key1)