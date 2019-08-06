from __future__ import print_function

import os
import sys
import time
import re
from io import StringIO
import base64

import requests

import pandas as pd
import numpy as np

__all__ = ["ProcessAnalysisData", "ProcessDataTableText", "RetrieveStudiesAnalysisAndResultsData", "SetupUIFDataForStudiesAnalysisAndResults, SetupCSVDownloadLink"]

def RetrieveStudiesAnalysisAndResultsData(StudyIDs, MWBaseURL):
    """Retrieve analysis and results data for a study or studies."""
    
    StudiesResultsData = {}
    
    StudyIDs = re.sub("[ ]+", " ", StudyIDs)
    
    for StudyID in StudyIDs.split(" "):
        MWDataURL = MWBaseURL + "/study/study_id/" + StudyID + "/analysis/"
        
        print("Initiating request: %s" % MWDataURL)
        Response = requests.get(MWDataURL)
        if Response.status_code != 200:
            print("Request failed: status_code: %d" % Response.status_code)
            continue
        
        AnalysisData = Response.json()
        
        print("Processing analysis data...")
        ProcessAnalysisData(StudiesResultsData, AnalysisData)
    
    for StudyID in StudiesResultsData:
        for AnalysisID in StudiesResultsData[StudyID]:
            print("\nRetrieving datatable for analysis ID, %s, in study ID, %s..." % (AnalysisID, StudyID))
            
            MWDataURL = MWBaseURL + "/study/analysis_id/" + AnalysisID + "/datatable"
            
            print("Initiating request: %s" % MWDataURL)
            Response = requests.get(MWDataURL)
            if Response.status_code != 200:
                print("***Error: Request failed: status_code: %d" % Response.status_code)
                continue
            
            print("Processing datatable text...")
            ResultsDataTable, ClassNamesToNumsMap = ProcessDataTableText(Response.text, AddClassNum = True)
            StudiesResultsData[StudyID][AnalysisID]["class_names_to_nums"] = ClassNamesToNumsMap
            
            print("Setting up Pandas dataframe...")
            RESULTSDATATABLE = StringIO(ResultsDataTable)
            StudiesResultsData[StudyID][AnalysisID]["data_frame"] = pd.read_csv(RESULTSDATATABLE, sep="\t", index_col = "Samples")
    
    return StudiesResultsData

def ProcessAnalysisData(StudiesResultsData, AnalysisData):
    """Process analysis data retrieved in JSON format for a study or set of studies"""
    
    if "study_id" in AnalysisData:
        # Turn single study with single analysis data set into dictionary
        # with multiple studies/analysis data set...
        AnalysisData = {"1" : AnalysisData}
    
    for DataSetNum in AnalysisData:
        StudyID = AnalysisData[DataSetNum]["study_id"]
        AnalysisID = AnalysisData[DataSetNum]["analysis_id"]
        
        # Intialize data...
        if StudyID not in StudiesResultsData:
            StudiesResultsData[StudyID] = {}
            
        
        StudiesResultsData[StudyID][AnalysisID] = {}
        
        # Track data...
        for DataType in AnalysisData[DataSetNum]:    
            DataValue = AnalysisData[DataSetNum][DataType]    
            if re.match("^(study_id|analysis_id)$", DataType, re.I):
                continue
            
            StudiesResultsData[StudyID][AnalysisID][DataType] = DataValue

def ProcessDataTableText(DataTableText, AddClassNum = True):
    """Process datatable retrieved retrieves in text format for a specific analysis ID"""
    
    DataLines = []
    
    TextLines = DataTableText.split("\n")
    
    # Process data labels...
    LineWords = TextLines[0].split("\t")
    
    DataLabels = []
    DataLabels.append(LineWords[0])
    DataLabels.append(LineWords[1])
    if AddClassNum:
        DataLabels.append("ClassNum")
    
    for Index in range(2, len(LineWords)):
        Name = LineWords[Index]
        DataLabels.append(Name)
    
    DataLines.append("\t".join(DataLabels))
    
    # Process data...
    ClassNamesMap = {}
    ClassNum = 0
    for Index in range(1, len(TextLines)):
        LineWords = TextLines[Index].split("\t")
        
        if len(LineWords) <= 2:
            continue
        
        # Handle sample ID and class name...
        DataLine = []
        DataLine.append(LineWords[0])
        DataLine.append(LineWords[1])
        
        if AddClassNum:
            ClassName = LineWords[1]
            if ClassName not in ClassNamesMap:
                ClassNum += 1
                ClassNamesMap[ClassName] = ClassNum
            DataLine.append("%s" % ClassNamesMap[ClassName])
            
        for Index in range(2, len(LineWords)):
            DataLine.append(LineWords[Index])
        
        DataLines.append("\t".join(DataLine))
    
    return ("\n".join(DataLines), ClassNamesMap)

def SetupUIFDataForStudiesAnalysisAndResults(StudiesResultsData, MinClassCount = None):
    """Retrieve data for setting up UIF from analysis and results data for a study or studies."""
    
    StudiesUIFData = {}
    StudiesUIFData["StudyIDs"] = []
    StudiesUIFData["AnalysisIDs"] = {}
    StudiesUIFData["MetaboliteIDs"] = {}
    StudiesUIFData["ClassIDs"] = {}
    
    for StudyID in StudiesResultsData:
        NewStudy = True
        
        for AnalysisID in StudiesResultsData[StudyID]:
            if "data_frame" not in StudiesResultsData[StudyID][AnalysisID]:
                print("***Warning: Excluding study ID, %s, analysis ID, %s, from further analysis: No named metabolities data available..." % (StudyID, AnalysisID))
                continue
            
            ResultsDataFrame = StudiesResultsData[StudyID][AnalysisID]["data_frame"]
            ColumnNames = list(ResultsDataFrame.columns.values)
            if len(ColumnNames) <= 3:
                print("***Warning: Excluding study ID, %s, analysis ID, %s, from further analysis: No named metabolities data available..." % (StudyID, AnalysisID))
                continue
            
            if MinClassCount is not None:
                ClassCount = len(StudiesResultsData[StudyID][AnalysisID]["class_names_to_nums"])
                if ClassCount < MinClassCount:
                    print("***Warning: Excluding study ID, %s, analysis ID, %s, from further analysis: Contains less than %s classes..." % (StudyID, AnalysisID, MinClassCount))
                    continue
            
            if NewStudy:
                NewStudy = False
                StudiesUIFData["StudyIDs"].append(StudyID)
                StudiesUIFData["AnalysisIDs"][StudyID] = []
                StudiesUIFData["MetaboliteIDs"][StudyID] = {}
                StudiesUIFData["ClassIDs"][StudyID] = {}
            
            StudiesUIFData["AnalysisIDs"][StudyID].append(AnalysisID)
            StudiesUIFData["MetaboliteIDs"][StudyID][AnalysisID] = []
            StudiesUIFData["ClassIDs"][StudyID][AnalysisID] = []
            
            StudiesUIFData["MetaboliteIDs"][StudyID][AnalysisID].extend(ColumnNames[3:])
            StudiesUIFData["ClassIDs"][StudyID][AnalysisID].extend(StudiesResultsData[StudyID][AnalysisID]["class_names_to_nums"])
    
    if len(StudiesUIFData["StudyIDs"]) == 0:
        print("***Warning: No studies available for further analysis...")

    return StudiesUIFData

#
# Reference:
# https://stackoverflow.com/questions/31893930/download-csv-from-an-ipython-notebook
#
def SetupCSVDownloadLink(DataFrame, Title = "Download CSV file", CSVFilename = "DataFrameDownload.csv"):  
    """Setup a HTML link for downloading a dataframe as a CSV file."""
    
    CSVData = DataFrame.to_csv()
    Base64EncodedData = base64.b64encode(CSVData.encode()).decode()
    HTMLText = '<a download="%s" href="data:text/csv;base64,%s" target="_blank">%s</a>' % (CSVFilename, Base64EncodedData, Title)
    
    return HTMLText
