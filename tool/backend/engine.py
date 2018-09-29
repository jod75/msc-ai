#!/usr/bin/env python

##################################################################################
# engine.py
from pyspark.sql.types import *
from pyspark.sql.functions import col, desc, concat, lit
from pyspark.sql import SQLContext

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem import MACCSkeys

from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast import NCBIXML

from moleculehelper import *
from chemblhelper import ChEMBLHelper
from pythonhelper import *
from hdfshelper import HDFSHelper

import os.path
import logging
import json
import subprocess
import shutil
import shlex
import re

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class ICS5200Engine(object):
    ################################################################################################################################################################
    # Private methods
    def __setup(self):
        """ Performs start-up routines:
            1. Loads data from ChEMBL MySQL and store it in a TSV file on Hadoop            
            2. Create DataFrames:
                a. Ligands:
                    - Ligands = unique list of Ligands                    
                b. Proteins:
                    - Proteins = unique list of Proteins
                c. Bindings:
                    - bindingsData = all bindings
            3. Create local blast database
            4. Summarize data             
        """

        # Step 1 - Get raw data
        chemblhelper = ChEMBLHelper()

        hdfsServer = "http://hadoop1:50070" # hdfs path
        localHome = "/home/hduser/Lab"
        hdfsHome = "/user/hduser/ics5200"
        datasetCount = 500000 # dataset count of bindings from ChEMBL

        self.sparkFastaFile = "/home/hduser/Lab/chembl" + \
            str(datasetCount) + \
            ".fasta" # file created by spark task
        self.localFasta = "/home/hduser/Lab/proteinbank/chembl" + str(datasetCount) + ".fasta"
        self.blastOutFile = "/home/hduser/Lab/blastResult.xml"

        self.databankFilename = chemblhelper.createDataBank(hdfsServer,
                                                            localHome,
                                                            hdfsHome,
                                                            datasetCount)

        # Step 2 - load data
        self.__loadData()

        # Step 3 - Create Blast DB
        self.__createLocalBlastDb()

        # Step 4 - Summarize data
        #self.__summarizeData()
        #==============================================
        #>Proteins:
        #>Total Proteins: 1875
        #>Test Bank Proteins: 1875
        #>Test Bank Bindings: 428266
        #>Ligands
        #>Total Ligands: 210090
        #>Test Bank Ligands: 210090
        #>Test Bank Bindings: 428266
        #==============================================
        
        PythonHelper.writeToJupyterConsole(">Engine ICS5200 is ready.")

    def __summarizeData(self):
        totalProteins = self.proteins.select("comp_id").distinct().count()
        totalLigands = self.ligands.select("mol_reg_no").distinct().count()
        bindings = self.bindingsData.count()
        PythonHelper.writeToJupyterConsole("==============================================")
        PythonHelper.writeToJupyterConsole(">Total Bindings: " + str(bindings))        
        PythonHelper.writeToJupyterConsole(">Total Proteins: " + str(totalProteins))                        
        PythonHelper.writeToJupyterConsole(">Total Ligands: " + str(totalLigands))                
        PythonHelper.writeToJupyterConsole("==============================================")

    def __loadData(self):
        """  Load raw data and create data sets
        """
        # dataRDD will contain all data imported from MySQL
        dataRDD = self.sc.textFile(self.databankFilename).map(lambda line: line.split('\t'))

        # convert each line (currently an list) to a tuple.  This makes it easier to manipulate
        # the data, especially to convert to DataFrames
        dataRDD = dataRDD.map(lambda l: tuple(l))

        ### DataFrames
        # Schemas
        ligandsSchema = StructType([StructField("mol_reg_no", IntegerType(), False),
                                    StructField("canonical_smiles", StringType(), False),
                                    StructField("mol_pref_name", StringType(), True)])

        proteinsSchema = StructType([StructField("comp_id", LongType(), False),
                                     StructField("prot_accession", StringType(), True),
                                     StructField("sequence", StringType(), False),
                                     StructField("prot_pref_name", StringType(), True),
                                     StructField("prot_short_name", StringType(), True),
                                     StructField("prot_description", StringType(), True)])

        bindingsSchema = StructType([StructField("row_id", LongType(), False),
                                     StructField("assay_id", LongType(), False),
                                     StructField("molregno", LongType(), False),
                                     StructField("std_relation", StringType(), True),
                                     StructField("std_value", DecimalType(), True),
                                     StructField("std_units", StringType(), True),
                                     StructField("std_type", StringType(), True),
                                     StructField("pchembl_value", DecimalType(), True),
                                     StructField("component_id", LongType(), False),
                                     StructField("canonical_smiles", StringType(), False)])

        # Unique Lists
        ligandsRDD = dataRDD.map(lambda t: (long(t[2]), LigandUtils.getCanonicalSmiles(str(t[11])), str(t[12]))).distinct()
        proteinsRDD = dataRDD.map(lambda t: (long(t[8]), str(t[9]), str(t[10]), str(t[13]), str(t[14]), str(t[15]))) \
                             .distinct()
        self.ligands = self.sqlContext.createDataFrame(ligandsRDD, ligandsSchema)
        self.proteins = self.sqlContext.createDataFrame(proteinsRDD, proteinsSchema)

        # Data sets
        self.bindingsData = self.sqlContext.createDataFrame(dataRDD.map(
            lambda t: (long(t[0].rstrip(".0")),
                       long(t[1]),
                       long(t[2]),
                       str(t[3]),
                       PythonHelper.getDecimal(t[4]),
                       str(t[5]),
                       str(t[6]),
                       PythonHelper.getDecimal(t[7]),
                       long(t[8]),
                       LigandUtils.getCanonicalSmiles(str(t[11])))), bindingsSchema)

    def __createLocalBlastDb(self):
        """ Create BLAST local db
        """

        # clean up before saving file to disk
        shutil.rmtree(self.sparkFastaFile, ignore_errors=True)

        # manipulate raw data rdd and create FASTA file
        # using "prot_pref_name" may result in duplicate values as it is not unique
        blastDb = self.bindingsData \
                    .join(self.proteins,
                          self.bindingsData.component_id == self.proteins.comp_id) \
                    .select("prot_accession", "sequence") \
                    .distinct() \
                    .rdd \
                    .map(lambda t: str(">ebl|" + t[0] + "|\r" + "\r".join(re.findall(".{1,80}",t[1]))))
        
        dump = blastDb.collect()
        myfile = open(self.localFasta, "w")
        for item in dump:
            myfile.write("%s\r" % item)            
        
        # create local blast db
        process = subprocess.Popen(shlex.split(
            "makeblastdb -in %s -parse_seqids -dbtype prot" % os.path.basename(self.localFasta)),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,  
            cwd=os.path.dirname(self.localFasta),
            shell=False)

        out, err = process.communicate()

        PythonHelper.writeToJupyterConsole(">Engine init makeblastdb out: " + out)
        PythonHelper.writeToJupyterConsole(">Engine init makeblastdb err: " + err)

    def __doLigandExperiment(self, querySmiles, molHelper=MoleculeHelper, fingerprintFunction=None, similarityFunction=None, similarityThreshold=0.5):
        """ Runs a ligand experiment

            Args:
                querySmiles: SMILES string

            Returns:
                List of known bindings
        """
        
        queryLigand = dict()
        queryLigand.update({0: querySmiles})
        queryRDD = self.sc.parallelize(queryLigand).map(lambda k:(k, molHelper(queryLigand[k], fingerprintFunction, similarityFunction)))
        ligandsRDD = self.ligands.select(self.ligands.mol_reg_no, self.ligands.canonical_smiles) \
                                              .distinct() \
                                              .rdd.map(lambda (k, v):(k, molHelper(v, fingerprintFunction, similarityFunction)))
        simRDD = ligandsRDD.cartesian(queryRDD) \
                           .map(lambda ((k1,v1),(k2,v2)): (k1, float(v1.similarity(v2)))) \
                           .filter(lambda (k1, v): v >= similarityThreshold)

        #queryMol = molHelper(querySmiles, fingerprintFunction, similarityFunction)        
        
        #uniqueLigandsRdd = self.ligands.select(self.ligands.mol_reg_no, self.ligands.canonical_smiles) \
        #                                      .distinct() \
        #                                      .rdd.map(lambda (k, v):(k, molHelper(v, fingerprintFunction, similarityFunction)))
        #uniqueFingerprints = uniqueLigandsRdd.map(lambda (k, v): (k, v.getFingerprint())).filter(lambda (k, v): v != None)
        #uniqueLigandsMap = uniqueFingerprints.collectAsMap()
        #sims = DataStructs.BulkTanimotoSimilarity(queryMol.getFingerprint(), uniqueLigandsMap.values())  
        #keys = uniqueLigandsMap.keys()
        #similarLigands = {} # ligand, tc
        #for x in range(len(sims)):
        #    if sims[x] >= similarityThreshold and keys[x] not in similarLigands:
        #        similarLigands[keys[x]] = sims[x]         
                
        #simRDD = self.sc.parallelize(similarLigands)        
        

        simSchema = StructType([StructField("molregno", LongType(), False),                                
                                StructField("similarity", FloatType(), False)])
        sim = self.sqlContext.createDataFrame(simRDD, simSchema)
        
        # get protein details
        simProteinsAndBindings = self.bindingsData.join(self.proteins.select("prot_description", "prot_accession", "comp_id"), 
                    self.bindingsData.component_id == self.proteins.comp_id)

        return sim.join(simProteinsAndBindings, simProteinsAndBindings.molregno == sim.molregno) \
                    .join(self.ligands.select("mol_reg_no", "mol_pref_name"),  sim.molregno == self.ligands.mol_reg_no) \
                    .orderBy(desc("similarity")).collect()

    ################################################################################################################################################################
    # Public methods  *** Ligands ***    
    
    def getSmiles(self, molRegNo):
        """ Get SMILES representation of a molecule.

            Args:
                molRegNo: long, molRegNo in ChEMBL                

            Returns:
                string: SMILES representation.
        """

        # there should be only one entry and one field
        ligandsSmilesList = self.ligands.filter(col("mol_reg_no") == long(molRegNo)).select("canonical_smiles").collect()
        if len(ligandsSmilesList) == 0:
            PythonHelper.writeToJupyterConsole(">Engine getSmiles no molecules with molregno: " + str(molRegNo))
            
        return ligandsSmilesList[0][0]       
    

    def doLigandExperiment(self, molRegNo, molHelper=MoleculeHelper, fingerprintFunction=None, similarityFunction=None, similarityThreshold=0.5):
        """ Runs a ligand experiment

            Args:
                molRegNo: long

            Returns:
                List of known bindings
        """

        querySmiles = self.getSmiles(molRegNo)
        return self.__doLigandExperiment(querySmiles, molHelper, fingerprintFunction, similarityFunction, similarityThreshold) 

    def doLigandExperimentFromSmiles(self, querySmiles, molHelper=MoleculeHelper, fingerprintFunction=None, similarityFunction=None, similarityThreshold=0.5):
        """ Runs a ligand experiment

            Args:
                querySmiles: SMILES string

            Returns:
                List of known bindings
        """

        PythonHelper.writeToJupyterConsole(">doLigandExperimentFromSmiles: " + querySmiles)
        return self.__doLigandExperiment(querySmiles, molHelper, fingerprintFunction, similarityFunction, similarityThreshold)  
    
    def isLigandInChEMBL(self, querySmiles):
        """ Checks whether the input querySmiles exists in ChEMBL data (or our sample)            

            Args:
                querySmiles: SMILES string

            Returns:
                long[]: If the ligand exists in our sample, then the function returns an array of molRegNo, else it returns an empty array.
        """

        querySmiles = LigandUtils.getCanonicalSmiles(querySmiles)
        molRegNoList = self.ligands.filter(col("canonical_smiles") == querySmiles).select("mol_reg_no").collect()
        if len(molRegNoList) == 0:
            PythonHelper.writeToJupyterConsole(">Engine isLigandInChEMBL no molecules with given canonical SMILES: " + querySmiles)
            
        return molRegNoList

    def getLigandBindings(self, querySmiles):
        """ Find bindings

            Args:
                querySmiles: SMILES string

            Returns:
                List of known bindings
        """
        
        querySmiles = LigandUtils.getCanonicalSmiles(querySmiles)
        PythonHelper.writeToJupyterConsole(">getLigandBindings: " + querySmiles)
        
        compList = self.bindingsData.filter(col("canonical_smiles") == querySmiles).select("component_id")
        result = compList.join(self.proteins.select("prot_accession", "comp_id"), 
                    compList.component_id == self.proteins.comp_id).select("comp_id", "prot_accession").distinct().orderBy("comp_id")
        return result.collect()
    
    ################################################################################################################################################################
    # Public methods  *** Proteins ***

    def doProteinExperiment(self, compId):    
        query_seq = self.proteins.filter(col("comp_id") == compId).select("sequence").collect()[0][0]
        blastp_cline = NcbiblastxCommandline(cmd = "blastp",
                                     db = self.localFasta,
                                     evalue = 0.01,
                                     outfmt = 5,
                                     remote = False,
                                     out = self.blastOutFile)

        (out, err) = blastp_cline(stdin = query_seq)
        PythonHelper.writeToJupyterConsole(">Engine doProteinExperiment out: " + out)
        PythonHelper.writeToJupyterConsole(">Engine doProteinExperiment err: " + err)
        
        # create matching proteins list.
        # each list entry is a tuple in the format:
        #  (Accession, Hit, High-Scoring-Pair expect value, hsp)        
        result_handle = open(self.blastOutFile)
        blast_record = NCBIXML.read(result_handle)
        matches = []
        for alignment in blast_record.alignments:   
            i = 1
            for hsp in alignment.hsps:
                matches.append((alignment.accession, i, hsp.expect, str(hsp)))
                i = i + 1 
        
        simRDD = self.sc.parallelize(matches)       
        simSchema = StructType([StructField("accession", StringType(), False),
                                StructField("hit", IntegerType(), False),
                                StructField("similarity", FloatType(), False),                                
                                StructField("hsp", StringType(), False)])
        simDF = self.sqlContext.createDataFrame(simRDD, simSchema)
        
        simProteins = simDF.join(self.proteins.select("prot_short_name", "prot_accession", "comp_id"), simDF.accession == self.proteins.prot_accession)
        fullKnownProteinDF = simProteins.join(self.bindingsData, simProteins.comp_id == self.bindingsData.component_id)
        return fullKnownProteinDF.select([c for c in fullKnownProteinDF.columns if c not in {"row_id", "assay_id", "prot_accession", "comp_id"}]).distinct().collect()
        
    def isProteinInChEMBL(self, sequence):        
        return self.proteins.filter(col("sequence") == sequence).select("comp_id", "prot_accession").distinct().orderBy("comp_id").collect()
        
    ################################################################################################################################################################
    # Constructor
    def __init__(self, sc, dataset_path):
        """Init the recommendation engine given a Spark context and a dataset path
        """

        logger.info("Starting up the ICS5200 Engine: ")

        self.sc = sc
        self.sqlContext = SQLContext(sc)

        # Load ratings data for later use
        logger.info("Loading data...")

        self.__setup()

