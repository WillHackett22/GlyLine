#GlyLine Trawler Class
#these functions are for going through mzML data
#thanks to Josh Klein (@mobiusklein) for ms_deisotope!

import pandas as pd
import numpy as np
import tables as tb
import ms_deisotope
import ms_deisotope.output.mzml as msdomzml
from GlyLine.Meta import Meta



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

#this class makes the target list, msppm is the ppm
class IonTargetList:
    def __init__(self,msppm=[10,20]):
        self.msppm=msppm
        
    def Maker(self,df,colid,mslvl):
        ix=[i for i,x in enumerate(df[colid].duplicated()) if x==False]
        temp=df.iloc[ix,]
        l, u =PPMBounds(temp['neutralmass'],self.msppm[mslvl])
        dfTargets=pd.DataFrame(data=None,columns=['upper','lower'])
        dfTargets['upper']=u
        dfTargets['lower']=l
        dfTargets.index=temp[colid]
        return dfTargets
        
    def MS2Maker(self,df,colid='IonID',mslvl=1):
        self.dfMS2Targets=self.Maker(df,colid,mslvl)
    
    def MS1Maker(self,df,colid='GPID',mslvl=0):
        self.dfMS1Targets=self.Maker(df,colid,mslvl)
        
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
            
#describe the data structure that trawler writes to
# DDA: runid | ionid | overlap | time | neutralmass | charge | intensity | score | precursorID
# DIA: runid | ionid | overlap | time | neutralmass | charge | intensity | score | precursorID | acq_range
class ProductPeakData(tb.IsDescription):
    RunID = tb.StringCol(16)
    IonID = tb.Int32Col()
    Overlap = tb.Int8Col()
    Time = tb.Float32Col()
    NeutralMass = tb.Float64Col()
    Charge = tb.Int8Col()
    Intensity = tb.Float64Col()
    Decon = tb.Float32Col()
    PrecursorID = tb.Int32Col()
    


#precursor data structure
# runid | time | neutralmass | charge | intensity | decon | precursorid
#precursorid= #MS1scan in run _ #acquisition in scanbunch: #MS1*100+#acq

class PrecursorPeakData(tb.IsDescription):
    RunID = tb.StringCol(16)
    Time = tb.Float32Col()
    NeutralMass = tb.Float64Col()
    Charge = tb.Int8Col()
    Intensity = tb.Float64Col()
    Decon = tb.Float32Col()
    PrecursorID = tb.Int32Col()
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
    def __init__(self,mzML,hdf5file,runidentifier=None,ms1key='MS1',ms2key='MS2',title=None,ms1title=None,ms2title=None,h5index=False,start_from=None,end_at=None,ms1_deconargs=None,collect_allMS1=True):
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
        self.PrecursorID=0
        self.collect_allMS1=collect_allMS1
        self.h5connection=tb.open_file(self.h5file,mode='a',title=self.title)
    
    def main(self,dfIonMetaMaster,dfPSMMetaMaster,msppm=[10,20]):
        self.IteratorGen()
        self.targetlist=IonTargetList(msppm)
        self.targetlist.MS2Maker(dfIonMetaMaster)
        self.targetlist.MS1Maker(dfPSMMetaMaster)
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
    
    def Trawling(self):
        for scan in self.iter:
            self.Scooper(scan)
            self.ms1table.flush()
            self.ms2table.flush()
    
    def LittleScoop(self,peak,acq):
        tempidx=self.targetlist.BoundMS2Index(peak.neutral_mass)
        ms2row=self.ms2table.row
        for iidx in range(len(tempidx)):
            ms2row['Time']=acq.scan_time.real
            ms2row['RunID']=self.runID
            ms2row['PrecursorID']=self.PrecursorID+self.prodN
            ms2row['Overlap']=sum([len(tempidx)>1])
            ms2row['IonID']=tempidx[iidx]
            ms2row['NeutralMass']=peak.neutral_mass
            ms2row['Charge']=peak.charge
            ms2row['Intensity']=peak.intensity
            ms2row['Decon']=peak.score
            ms2row.append()
    
    def BigScoop(self,acq,pre,time,mz=False):
        if pre!=None:
            tempidx=self.targetlist.BoundMS1Index(pre.neutral_mass)
            if self.collect_allMS1:
                if len(tempidx)==0:
                    ms1row=self.ms1table.row
                    ms1row['GPID']=0
                    ms1row['Overlap']=sum([len(tempidx)>1])
                    ms1row['RunID']=self.runID
                    ms1row['Time']=time
                    ms1row['NeutralMass']=pre.neutral_mass
                    ms1row['Charge']=pre.charge
                    ms1row['Intensity']=pre.intensity
                    if acq.precursor_information.defaulted | acq.precursor_information.orphan:
                        ms1row['Decon']=0
                    else:
                        ms1row['Decon']=pre.score
                    ms1row['PrecursorID']=self.PrecursorID+self.prodN
                    ms1row.append()
                else:  
                    for iidx in range(len(tempidx)):
                        ms1row=self.ms1table.row
                        ms1row['GPID']=tempidx[iidx]
                        ms1row['Overlap']=sum([len(tempidx)>1])
                        ms1row['RunID']=self.runID
                        ms1row['Time']=time
                        ms1row['NeutralMass']=pre.neutral_mass
                        ms1row['Charge']=pre.charge
                        ms1row['Intensity']=pre.intensity
                        if acq.precursor_information.defaulted | acq.precursor_information.orphan:
                            ms1row['Decon']=0
                        else:
                            ms1row['Decon']=pre.score
                        ms1row['PrecursorID']=self.PrecursorID+self.prodN
                        ms1row.append()


class TrawlerDDA(Trawler):
    def __init__(self,mzML,hdf5file,runidentifier=None,ms1key='MS1',ms2key='MS2',title=None,ms1title=None,ms2title=None,h5index=False,start_from=None,end_at=None,ms1_deconargs=None):
        super().__init__(mzML,hdf5file,runidentifier,ms1key,ms2key,title,ms1title,ms2title,h5index,start_from,end_at,ms1_deconargs)
            
    def Scooper(self,scan):
        self.PrecursorID+=1000
        self.prodN=0
        if len(scan.products)!=0:
            for acq in scan.products:  
                self.prodN+=1
                premass=acq.precursor_information.extracted_neutral_mass
                pre=scan.precursor.deconvoluted_peak_set.has_peak(premass)
                self.BigScoop(acq,pre,scan.precursor.scan_time.real)
                for peak in acq.deconvoluted_peak_set:
                    if self.targetlist.BoundMS2Bool(peak.neutral_mass):
                        self.LittleScoop(peak,acq)
            

### DIA Functionality below this line

class TrawlerDIA(Trawler):
    def __init__(self,mzML,hdf5file,runidentifier=None,ms1key='MS1',ms2key='MS2',title=None,ms1title=None,ms2title=None,h5index=False,start_from=None,end_at=None,ms1_deconargs=None):
        super().__init__(mzML, hdf5file,runidentifier,ms1key,ms2key,title,ms1title,ms2title,h5index,start_from,end_at,ms1_deconargs)
            
    def Scooper(self,scan):
        self.PrecursorID+=1000
        self.prodN=0
        plist=np.array([s.mz for s in scan.precursor.deconvoluted_peak_set])
        for acq in scan.products:
            self.prodN+=1
            if any((acq.isolation_window.upper_bound > plist) & (plist > acq.isolation_window.lower_bound)):
                hold=plist[(acq.isolation_window.upper_bound > plist) & (plist > acq.isolation_window.lower_bound)]
                for pretemp in hold:
                    pre=scan.precursor.deconvoluted_peak_set.has_peak(pretemp,use_mz=True)
                    self.BigScoop(pre,scan.precursor.scan_time.real,mz=True)
                for peak in acq.deconvoluted_peak_set:
                    if self.targetlist.BoundMS2Bool(peak.neutral_mass):
                        self.LittleScoop(peak,acq)
            
        
        