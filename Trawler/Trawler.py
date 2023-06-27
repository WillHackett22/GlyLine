#GlyLine Trawler Class
#these functions are for going through mzML data
#thanks to Josh Klein (@mobiusklein) for ms_deisotope!

import pandas as pd
import numpy as np
import tables as tb
import json
import ms_deisotope.output.mzml as msdomzml

#let's iterate through things safely
def safeNext(source):
    try:
        result = next(source)
    except KeyError:
        result = next(source)
    return result

def safeIter(source):
    try:
        while True:
            yield safeNext(source)
    except StopIteration:
        pass

def PPMBounds(array,threshold=20):
    upper=array+array*threshold/1000000
    lower=array-array*threshold/1000000
    return lower, upper

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

#this class is to help speed up our target matching by filtering neutral masses with binsearch first
class ReduceByRounding:
    def __init__(self,targetdf):
        self.targetdf=targetdf
        self.colid=targetdf.index.name
    
    def RoundedDFGenerator(self):
        roundeddf=pd.DataFrame(None,columns=['rounded',self.colid])
        roundeddf['rounded']=round(self.targetdf['upper']).tolist()+round(self.targetdf['lower']).tolist()
        roundeddf[self.colid]=self.targetdf.index.values.tolist()*2
        runix=[i for i,x in enumerate(roundeddf.duplicated(['rounded',self.colid])) if x==False]
        roundeddfout=roundeddf.iloc[runix,]
        self.roundeddf=roundeddfout
    
    def RoundedGrouper(self):
        groupeddf=self.roundeddf.groupby(['rounded'])
        self.groupeddf=groupeddf

    def RoundedUniqueNSort(self):
        uniround=self.roundeddf['rounded'].unique()
        uniroundsort=uniround[uniround.argsort()].tolist()
        self.uniroundsort=uniroundsort
        
    def RoundedInitialize(self):
        self.RoundedDFGenerator()
        self.RoundedGrouper()
        self.RoundedUniqueNSort()
        
    def RoundedTargetIDs(self,peakvalue):
        roundpeak=round(peakvalue)
        roundpeakloc=binsearch(self.uniroundsort, roundpeak)
        if roundpeakloc!=-1:
            rval=self.uniroundsort[roundpeakloc]
            temptargetids=np.unique(self.groupeddf.get_group(rval)[self.colid].tolist()).tolist()
        else:
            temptargetids=[]
        return temptargetids
    
    def ReducedTarget(self,temptargetids):
        temptargetdf=self.targetdf.loc[temptargetids]
        return temptargetdf

#this class makes the target list, msppm is the ppm
class IonTargetList:
    def __init__(self,msppm=[10,20]):
        self.msppm=msppm
        
    def Maker(self,df,mslvl):
        l, u =PPMBounds(df['neutralmass'],self.msppm[mslvl])
        dfTargets=pd.DataFrame(data=None,columns=['upper','lower'])
        dfTargets['upper']=u
        dfTargets['lower']=l
        dfTargets.index=df.index
        return dfTargets
        
    def MS2Maker(self,df,mslvl=1):
        dfMS2Targets=self.Maker(df,mslvl)
        self.dfMS2Objects=self.RoughMaker(dfMS2Targets)
    
    def MS1Maker(self,df,mslvl=0):
        dfMS1Targets=self.Maker(df,mslvl)
        self.dfMS1Objects=self.RoughMaker(dfMS1Targets)
    
    def RoughMaker(self,targetdf):
        roundedobj=ReduceByRounding(targetdf)
        roundedobj.RoundedInitialize()
        return roundedobj
    
    def BoundBool(self,peakvalue,targetdf):
        return any((targetdf['upper']>=peakvalue) & (targetdf['lower']<=peakvalue))
        
    def BoundIndex(self,peakvalue,targetdf):
        ID=targetdf.index[np.where((targetdf['upper']>=peakvalue) & (targetdf['lower']<=peakvalue))].values
        return ID
        
    def BoundMS2Bool(self,peakvalue):
         return any((self.dfMS2Targets['upper']>=peakvalue) & (self.dfMS2Targets['lower']<=peakvalue))
            
    def BoundMS2Index(self,peakvalue):
        ID=self.dfMS2Targets.index[np.where((self.dfMS2Targets['upper']>=peakvalue) & (self.dfMS2Targets['lower']<=peakvalue))].values
        return ID
    
    def BoundMS1Bool(self,peakvalue):
         return any((self.dfMS1Targets['upper']>=peakvalue) & (self.dfMS1Targets['lower']<=peakvalue))
            
    def BoundMS1Index(self,peakvalue):
        ID=self.dfMS1Targets.index[np.where((self.dfMS1Targets['upper']>=peakvalue) & (self.dfMS1Targets['lower']<=peakvalue))].values
        return ID
            
    
#class of functions to be used in trawler to get the indexed data from json file
class IndexedMSInfo:
    def __init__(self,jsonfile,runID,h5file,key='IndexInfo',verbose=True):
        self.jsonfile=jsonfile
        self.verbose=verbose
        self.runID=runID
        self.h5file=h5file
        self.idxkey=key+'_'+runID
        f=open(jsonfile)
        self.jdata=json.load(f)
        
    def MS1Info(self):
        ms1idxdf=pd.DataFrame(None,columns=['PrecursorIdx','scan_time','scan_id'])
        ms1idxdf['PrecursorIdx']=list(range(len(self.jdata['ms1_ids'])))
        scan_times=[]
        scan_ids=[]
        for jv in self.jdata['ms1_ids']:
            scan_ids=scan_ids+[jv[0]]
            scan_times=scan_times+[jv[1]['scan_time']]
        ms1idxdf['scan_time']=scan_times
        ms1idxdf['scan_id']=scan_ids
        ms1idxdf=ms1idxdf.set_index('PrecursorIdx')
        if self.verbose:
            self.MS1Data=ms1idxdf
        ms1idxdf.to_hdf(self.h5file,key=self.idxkey+'_MS1',mode='a')
    
    def MS2Info(self,targetlist):
        ms2idxdf=pd.DataFrame(None,columns=['ProductIdx','scan_time','scan_id','neutralmass','coisolation','Pre_scan_id'])
        ms2idxdf['ProductIdx']=list(range(len(self.jdata['msn_ids'])))
        scan_times=[]
        scan_ids=[]
        pre_scan_ids=[]
        nms=[]
        gpids=[]
        coi=[]
        for jv in self.jdata['msn_ids']:
            scan_ids=scan_ids+[jv[0]]
            scan_times=scan_times+[jv[1]['scan_time']]
            pre_scan_ids=pre_scan_ids+[jv[1]['precursor_scan_id']]
            if len(jv[1]['coisolation'])>0:
                coi=coi+[True]
            else:
                coi=coi+[False]
            nmtemp=jv[1]['neutral_mass']
            nms=nms+[nmtemp]
            tempids=targetlist.dfMS1Objects.RoundedTargetIDs(nmtemp)
            if len(tempids)>0:
                temptarg=targetlist.dfMS1Objects.ReducedTarget(tempids)
                hits=targetlist.BoundIndex(nmtemp,temptarg).tolist()
            else:
                hits=[0]
            if len(hits)==0:
                hits=[0]
            gpids=gpids+[hits]
        ms2idxdf['scan_time']=scan_times
        ms2idxdf['scan_id']=scan_ids
        ms2idxdf['Pre_scan_id']=pre_scan_ids
        ms2idxdf['neutralmass']=nms
        ms2idxdf=ms2idxdf.set_index('ProductIdx')
        ms2idxdf['PrecursorIdx']=[None]*ms2idxdf.shape[0]
        ms2idxdf['GPID']=gpids
        ms2idxdf['coisolation']=coi
        for j in self.MS1Data.index.tolist():
            ms2idxdf.loc[ms2idxdf['Pre_scan_id']==self.MS1Data['scan_id'].loc[j],'PrecursorIdx']=j
        if self.verbose:
            self.MS2Data=ms2idxdf
        ms2idxdf.to_hdf(self.h5file,key=self.idxkey+'_MS2',mode='a')
        
    def main(self,targetlist,force=False):
        if force:
            self.MS1Info()
            self.MS2Info(targetlist)
        else:
            try:
                self.MS1Data=pd.read_hdf(self.h5file,self.idxkey+'_MS1')
            except:
                self.MS1Info()
            try:
                self.MS2Data=pd.read_hdf(self.h5file,self.idxkey+'_MS2')
            except:
                self.MS2Info(targetlist)
        
#describe the data structure that trawler writes to
# DDA: runid | ionid | overlap | time | neutralmass | charge | intensity | score | precursorIdx | productIdx
# DIA: runid | ionid | overlap | time | neutralmass | charge | intensity | score | precursorIdx | productIdx | acq_range
#new data recorded to account for indexeddata
class ProductPeakData(tb.IsDescription):
    RunID = tb.StringCol(16)
    IonID = tb.Int32Col()
    Overlap = tb.Int8Col()
    NeutralMass = tb.Float64Col()
    Charge = tb.Int8Col()
    Intensity = tb.Float64Col()
    Decon = tb.Float32Col()
    ProductIdx = tb.Int32Col()
    
#precursor data structure
#new data recorded to account for indexeddata
class PrecursorPeakData(tb.IsDescription):
    NeutralMass = tb.Float64Col()
    Charge = tb.Int8Col()
    Intensity = tb.Float64Col()
    Decon = tb.Float32Col()
    PrecursorIdx = tb.Int32Col()
    GPID = tb.Int32Col()
    Overlap = tb.Int8Col()    
    
# ms1_deconvolution_args = {
#             "scorer": ms_deisotope.scoring.PenalizedMSDeconVFitter(score_threshold, isotopic_strictness),
#             "max_missed_peaks": missed_peaks,
#             "averagine": averagine,
#             "truncate_after": workflow.SampleConsumer.MS1_ISOTOPIC_PATTERN_WIDTH,
#             "ignore_below": workflow.SampleConsumer.MS1_IGNORE_BELOW,
#             "deconvoluter_type": ms1_deconvoluter_type,
#             "use_quick_charge": True
#         }
# default_glycresoft_ms1_deconvolution_args = {
#             "scorer": ms_deisotope.scoring.PenalizedMSDeconVFitter(20, isotopic_strictness),
#             "max_missed_peaks": 3,
#             "averagine": ms_deisotope.averagine.glycan,
#             "truncate_after": workflow.SampleConsumer.MS1_ISOTOPIC_PATTERN_WIDTH,
#             "ignore_below": workflow.SampleConsumer.MS1_IGNORE_BELOW,
#             "deconvoluter_type": ms1_deconvoluter_type,
#             "use_quick_charge": True
#         }
# preferred_glycresoft_ms1_deconvolution_args = {
#             "scorer": ms_deisotope.scoring.PenalizedMSDeconVFitter(score_threshold, isotopic_strictness),
#             "max_missed_peaks": missed_peaks,
#             "averagine": averagine,
#             "truncate_after": workflow.SampleConsumer.MS1_ISOTOPIC_PATTERN_WIDTH,
#             "ignore_below": workflow.SampleConsumer.MS1_IGNORE_BELOW,
#             "deconvoluter_type": ms1_deconvoluter_type,
#             "use_quick_charge": True
#         }



#this class will go through 1 mzML file
class Trawler:
    def __init__(self,mzML,hdf5file,runidentifier=None,ms1key='MS1',ms2key='MS2',
                 title=None,ms1title=None,ms2title=None,h5index=False,start_from=None,end_at=None,ms1_deconargs=None,collect_allMS1=True):
        self.scan_object=msdomzml.ProcessedMzMLDeserializer(mzML)
        self.scan_source=mzML
        self.h5file=hdf5file
        self.ms1key=ms1key
        self.ms2key=ms2key
        self.title=title
        if ms1title==None:
            self.ms1title=ms1key
        else:
            self.ms1title=ms1title
        if ms2title==None:
            self.ms2title=ms2key
        else:
            self.ms2title=ms2title
        self.h5index=h5index
        self.startscan=start_from
        self.endscan=end_at
        self.runID=runidentifier
        self.PrecursorIdx=0
        self.ProductIdx=0
        self.ms1counter=0
        self.ms2counter=0
        self.collect_allMS1=collect_allMS1
    
    def main(self,dfIonMetaObject,dfPSMMetaObject,msppm=[10,20]):
        self.IteratorGen()
        self.targetlist=IonTargetList(msppm)
        self.targetlist.MS2Maker(dfIonMetaObject)
        self.targetlist.MS1Maker(dfPSMMetaObject)
        self.h5connection=tb.open_file(self.h5file,mode='a',title=self.title)
        self.PrecursorTbMake()
        self.ProductTbMake()
        self.Trawling()
        self.h5connection.close()
            
    def IteratorGen(self):
        if self.startscan!=None:
            self.iter=safeIter(self.scan_object.start_from_scan(self.startscan))
        else:
            self.iter=safeIter(self.scan_object)
    
    def PrecursorTbMake(self):
        try:
            self.group1=self.h5connection.create_group('/','MS1','MS1 Data')
            self.ms1table=self.h5connection.create_table(self.group1,self.ms1key,PrecursorPeakData,title=self.ms1title)
        except:
            self.group1=self.h5connection.root.MS1
            self.ms1table=self.group1.MS1
    
    def ProductTbMake(self):
        try:
            self.group2=self.h5connection.create_group('/','MS2','MS2 Data')
            self.ms2table=self.h5connection.create_table(self.group2,self.ms2key,ProductPeakData,title=self.ms2title)
        except:
            self.group2=self.h5connection.root.MS2
            self.ms2table=self.group2.MS2
    
    def MS2RowCollect(self,prod,peak,hit,hitct):
        ms2row=self.ms2table.row
        ms2row['ProductIdx']=self.ProductIdx
        ms2row['Overlap']=hitct
        ms2row['NeutralMass']=peak.neutral_mass
        ms2row['Charge']=peak.charge
        ms2row['Intensity']=peak.intensity.real
        ms2row['Decon']=peak.score
        ms2row['IonID']=hit
        ms2row.append()
    
    def MS1RowCollect(self,scan,peak,hit,hitct):
        ms1row=self.ms1table.row
        ms1row['PrecursorIdx']=self.PrecursorIdx 
        ms1row['Overlap']=hitct
        ms1row['NeutralMass']=peak.neutral_mass
        ms1row['Charge']=peak.charge
        ms1row['Intensity']=peak.intensity.real
        ms1row['Decon']=peak.score
        ms1row['GPID']=hit
        ms1row.append()
    
    def CheckMSTargetsSub(self,roundedobj,peakmass):
        tempids=roundedobj.RoundedTargetIDs(peakmass)
        if len(tempids)>0:
            temptarg=roundedobj.ReducedTarget(tempids)
            hits=self.targetlist.BoundIndex(peakmass,temptarg)
        else:
            hits=[]
        return hits
    
    def CheckMS1Targets(self,scan,peak):
        hits=self.CheckMSTargetsSub(self.targetlist.dfMS1Objects,peak.neutral_mass)
        if len(hits)>0:
            for hit in hits:
                self.MS1RowCollect(scan,peak,hit,len(hits))
                self.ms1counter+=1
    
    def CheckMS2Targets(self,prod,peak,prehits):
        hits=self.CheckMSTargetsSub(self.targetlist.dfMS2Objects,peak.neutral_mass)
        if len(hits)>0:
            for hit in hits:
                self.MS2RowCollect(prod,peak,hit,len(hits))
                self.ms2counter+=1
    
    def ProductScanChecker(self,prod):
        if prod.precursor_information.orphan | prod.precursor_information.defaulted:
            prehits=[0]
        else:
            prehits=self.CheckMSTargetsSub(self.targetlist.dfMS1Objects,prod.precursor_information.neutral_mass)
            if len(prehits)==0:
                prehits=[0]
        for peak in prod.deconvoluted_peak_set:
            self.CheckMS2Targets(prod,peak,prehits)
        self.ProductIdx+=1
        
    def Trawling(self):
        for scan in self.iter:
            self.Scooper(scan)
            self.PrecursorIdx+=1
        self.ms1table.flush()
        self.ms2table.flush()
                
    def Scooper(self,scan):
        if self.collect_allMS1: 
            for peak in scan.precursor.deconvoluted_peak_set:    
                self.CheckMS1Targets(scan,peak)                
        for prod in scan.products:
            self.ProductScanChecker(prod)
        if self.ms1counter>=200:
            self.ms1table.flush()
            self.ms1counter=0
        if self.ms2counter>=200:
            self.ms2table.flush()
            self.ms2counter=0