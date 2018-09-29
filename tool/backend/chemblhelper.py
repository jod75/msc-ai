
# coding: utf-8

# # ChEMBL Helper methods
# 
# The puspose of this notebook is to provide helper classes and methods related to ChEMBL data extraction.

# In[1]:

import mysql.connector
from mysql.connector import errorcode
from pythonhelper import *
from hdfshelper import HDFSHelper
import os.path


# In[ ]:

class ChEMBLHelper:
    """ This helper class provides methods to get data from ChEMBL MySQL database.
    """
    
    __config = None
    
    def __init__(self, config = None):
        """ Creates an instance of the helper class using predefined configuration if [config] is None.
        
            Args:
                config: a MySQL connection configuration.
        """
        
        if config == None:
            self.__config = {
              'user': 'joseph',
              'password': 'password',
              'host': '192.168.151.11',
              'database': 'chembl_22',
              'raise_on_warnings': False,
            }
        else:
            self.__config = config  
            
    def createDataBank(self, hdfsServer, localLocation, hdfsLocation, datasetCount = 100):
        # get small dataset from ChEMBL
        datasetTSVFilename = "sample" + str(datasetCount) + ".tsv"
        datasetTSVPath = os.path.join(localLocation, datasetTSVFilename)
        hdfsDatasetFilename = os.path.join(hdfsLocation, datasetTSVFilename)

        # check if the dataset file exists in hdfs and if it does not, then load data from ChEMBL database
        if not HDFSHelper.fileExists(hdfsServer, hdfsDatasetFilename):
            # get data from ChEMBL and stores the file to local dfs
            self.saveBindingsToTSV(datasetTSVPath, datasetCount) 
            # upload data to hdfs so that it is accessbile from all cluster worker nodes
            HDFSHelper.putFile(hdfsServer, datasetTSVPath, hdfsDatasetFilename)
            
        return hdfsDatasetFilename
        
    def getMolecules(self, limit):
        """ Gets a number (limit) of molecule related data and stores it in a dictionary.  The dictionary key is the MOLREGNO.
            
            Args:
                limit: the number of molecules to return. If limit is less than 1, it returns all data.
                
            Returns:
                dict() in the format:
                    <molregno>, <canonical_smiles>
        """
        molecules = dict([])

        try:
            cnx = mysql.connector.connect(**(self.__config))

            cursor = cnx.cursor()
            
            
            limitClause = ""
            if limit > 1:
                limitClause = " LIMIT " + str(limit)

            query = ("SELECT DISTINCT cs.molregno, cs.canonical_smiles " +
                   "FROM compound_structures cs, activities act, " +
                   "     assays asy " +
                   "WHERE cs.molregno = act.molregno AND" +
                   "     act.assay_id = asy.assay_id AND" +
                   "     asy.assay_tax_id = 9606 " + limitClause + ";")

            cursor.execute(query)

            for (molregno, canonical_smiles) in cursor:
                molecules.update({molregno: canonical_smiles})

        except mysql.connector.Error as err:
            if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
                print("Something is wrong with your user name or password")
            elif err.errno == errorcode.ER_BAD_DB_ERROR:
                print("Database does not exist")
            else:
                print(err)
        else:
            cnx.close()
            
        return molecules
    
    def getBindings(self, limit):
        """ Get molecule-protein bindings from ChEMBL. This is an expensive operation as multiple tables need to be
            joined together and filtered.  Only those bindings related to Homo Sapiens are returned.  These are identified by
            taxanomy id 9606.  To improve performance, a table called [binding] is created.  This is a heap table with all the
            data required for this study.  This table is created the first time this method is called.
            Data from [binding] data is retreived in batches of 50000.  This value is chosen large enough to reduce the 
            number of calls hitting the database but small enough for MySQL Python Connector not to give timeouts. This method
            is chosen over increasing the timeout value so as that it is not bound to the specifications of the current 
            infrastructure.            
            
            Args:
                limit: the number of bindings to return. Limit < 1 returns all data.  Note that this may take some time 
                to complete.
               
            Returns:
                a dict() with all binding details in the following format:
                    <row_id>, <assay_id, molrgeno, standard_relation, standard_value, 
                               standard_units, standard_type, pchembl_value, componend_id, 
                               accession, sequence, canonical_smiles>                           
        """
        # create a custom table with all our needed data
        # this is required to page data
        self.__createChEMBLBindingTable()
        
        bindings = dict([])

        try:
            cnx = mysql.connector.connect(**(self.__config))

            cursor = cnx.cursor()
            
            pagesize = 50000

            limitClause = " LIMIT " + str(limit)
            if limit < 1:                
                limitClause = " LIMIT " + str(pagesize) + " OFFSET %(pageoffset)s " 
                
            query = ("SELECT row_id, "+
                     "     assay_id, molregno, " +
                     "     standard_relation, CAST(standard_value as DECIMAL(38,30)), " +
                     "     standard_units, standard_type, " +
                     "     CAST(pchembl_value as DECIMAL(38,30)), component_id, " +
                     "     accession, sequence, canonical_smiles "+
                     "FROM binding " +                     
                     "ORDER BY row_id " +
                     limitClause + ";")            
            
            i = 0            
            PythonHelper.writeToJupyterConsole(">getBindings - pagesize: " + str(pagesize))            
            while True:                     
                PythonHelper.writeToJupyterConsole(">getBindings - Count: " + str(i) + ", Offset: " + str(i*pagesize))
                cursor.execute(query, {"pageoffset":i * pagesize})                
                rowsRead = False # need to track if there were any rows in resultset due to non buffered cursor mode                 
                for (row_id, assay_id, molregno, std_relation, std_value, std_units, std_type,
                     pchembl_value, component_id, accession, sequence, canonical_smiles) in cursor:
                    bindings.update({row_id: [assay_id, molregno, std_relation, std_value, std_units, 
                                        std_type, pchembl_value, component_id, accession, 
                                        sequence, canonical_smiles]})                    
                    if not rowsRead:   #checking a flag must be faster than setting it <- confirm this?!
                        rowsRead = True                
                i = i + 1
                if not rowsRead or limit > 0:                    
                    break;                 

        except mysql.connector.Error as err:
            if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
                print("Something is wrong with your user name or password")
            elif err.errno == errorcode.ER_BAD_DB_ERROR:
                print("Database does not exist")
            else:
                print(err)
        else:
            cnx.close()
            
        return bindings
    
    def saveBindingsToTSV(self, filename, limit):
        """ Store binding data in a Tab Seperated Values file.
            Get molecule-protein bindings from ChEMBL. This is an expensive operation as multiple tables need to be
            joined together and filtered.  Only those bindings related to Homo Sapiens are returned.  These are identified by
            taxanomy id 9606.  To improve performance, a table called [binding] is created.  This is a heap table with all the
            data required for this study.  This table is created the first time this method is called.
            Data from [binding] data is retreived in batches of 50000.  This value is chosen large enough to reduce the 
            number of calls hitting the database but small enough for MySQL Python Connector not to give timeouts. This method
            is chosen over increasing the timeout value so as that it is not bound to the specifications of the current 
            infrastructure.            
            
            Args:
                filename: the tsv filename where to save the data.  The columns are in this order:
                    <row_id, assay_id, molrgeno, standard_relation, standard_value, 
                     standard_units, standard_type, pchembl_value, componend_id, 
                     accession, sequence, canonical_smiles, mol_pref_name, prot_pref_name, 
                     prot_short_name, prot_description>
                limit: the number of bindings to return.
        """
        # create a custom table with all our needed data
        # this is required to page data
        self.__createChEMBLBindingTable()                

        try:
            cnx = mysql.connector.connect(**(self.__config))

            cursor = cnx.cursor()
            
            pagesize = 50000

            limitClause = " LIMIT " + str(limit)
            if limit < 1:                
                limitClause = " LIMIT " + str(pagesize) + " OFFSET %(pageoffset)s " 
                PythonHelper.writeToJupyterConsole(">storeBindingsToTSV - pagesize: " + str(pagesize))
                
            query = ("SELECT row_id, "+
                     "     assay_id, molregno, " +
                     "     standard_relation, CAST(standard_value as DECIMAL(38,30)), " +
                     "     standard_units, standard_type, " +
                     "     CAST(pchembl_value as DECIMAL(38,30)), component_id, " +
                     "     accession, sequence, canonical_smiles, mol_pref_name, prot_pref_name, short_name, prot_description "+
                     "FROM binding " +                     
                     "ORDER BY row_id " +
                     limitClause + ";")            
            
            i = 0                        
            
            PythonHelper.writeToJupyterConsole(">storeBindingsToTSV - " + filename)
            with open(filename, 'w') as tsv:
                while True:  
                    if limit < 1:
                        PythonHelper.writeToJupyterConsole(">storeBindingsToTSV - Count: " + 
                                                           str(i) + ", Offset: " + str(i*pagesize))
                    cursor.execute(query, {"pageoffset":i * pagesize})                
                    rowsRead = False # need to track if there were any rows in resultset due to non buffered cursor mode                 
                    for (row_id, assay_id, molregno, std_relation, std_value, std_units, std_type,
                         pchembl_value, component_id, accession, sequence, canonical_smiles,
                         mol_pref_name, prot_pref_name, prot_short_name, prot_description) in cursor:
                        # as per os.linesp - use \n for terminator on all platforms
                        tsv.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % 
                                  (row_id, assay_id, molregno, std_relation, std_value, std_units,
                                   std_type, pchembl_value, component_id, accession, 
                                   sequence, canonical_smiles, 
                                   mol_pref_name, prot_pref_name, prot_short_name, prot_description))
                        if not rowsRead:   #checking a flag must be faster than setting it <- confirm this?!
                            rowsRead = True                
                    i = i + 1
                    if not rowsRead or limit > 0:                    
                        break;                 

        except mysql.connector.Error as err:
            if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
                print("Something is wrong with your user name or password")
                raise
            elif err.errno == errorcode.ER_BAD_DB_ERROR:
                print("Database does not exist")
                raise
            else:
                print(err)
                raise
        else:
            cnx.close()            
    
    def __createChEMBLBindingTable(self):
        try:
            cnx = mysql.connector.connect(**(self.__config))

            cursor = cnx.cursor()
            
            
            mysqlCreateBindingTable = (            
                "CREATE TABLE IF NOT EXISTS binding AS " +
                "SELECT @rownum := @rownum + 1 AS row_id, " +
                "     act.assay_id, act.molregno, " +
                "     act.standard_relation, act.standard_value, " +
                "     act.standard_units, act.standard_type, " +
                "     act.pchembl_value, tc.component_id, " +
                "     cs.accession, cs.sequence, cs.description AS prot_description, cps.canonical_smiles, "+
                "     md.pref_name AS mol_pref_name, pc.pref_name AS prot_pref_name," +
                "     pc.short_name " +
                "FROM (select @rownum := 0) r, activities act, assays asy,  " +
                "     target_components tc, component_sequences cs, " +
                "     compound_structures cps, molecule_dictionary md, " +
                "     component_class cc, protein_classification pc " +
                "WHERE act.assay_id = asy.assay_id AND " +
                "      asy.tid = tc.tid AND tc.component_id = cs.component_id AND "+ 
                "      cs.component_id = cc.component_id AND " +
                "      cc.protein_class_id = pc.protein_class_id AND " +
                "      asy.assay_tax_id = 9606 AND act.molregno = md.molregno AND " +
                "      md.molregno = cps.molregno AND " +                
                "      (act.pchembl_value IS NULL OR act.pchembl_value >= 4.75) AND " +
                "      act.standard_value > 0 AND " +
                "      ( " +
                "        (act.standard_type = 'Inhibition' and lower(act.standard_units) = '%' and act.standard_value >= 70 and act.standard_value <= 100 and act.standard_relation in ('>', '=', '>=', '~')) " +
                "        OR (act.standard_type = 'IC50' and lower(act.standard_units) = 'nm' and act.standard_value <= 25000 and act.standard_relation in ('<', '=', '<=', '~')) " +
                "        OR (act.standard_type = 'IC50' and lower(act.standard_units) in ('ug.ml-1', 'ug ml-1') and act.standard_value <= 8 and act.standard_relation in ('<', '=', '<=', '~')) " +
                "        OR (act.standard_type = 'IC90' and lower(act.standard_units) = 'nm' and act.standard_value <= 100000 and act.standard_relation in ('<', '=', '<=', '~')) " +
                "        OR (act.standard_type = 'IC90' and lower(act.standard_units) in ('ug.ml-1', 'ug ml-1') and act.standard_value <= 32 and act.standard_relation in ('<', '=', '<=', '~')) " +
                "        OR (act.standard_type = 'MIC' and lower(act.standard_units) = 'nm' and act.standard_value <= 25000 and act.standard_relation in ('<', '=', '<=', '~'))     " +
                "        OR (act.standard_type = 'MIC' and lower(act.standard_units) in ('ug.ml-1', 'ug ml-1') and act.standard_value <= 8 and act.standard_relation in ('<', '=', '<=', '~')) " +
                "        OR (act.standard_type = 'MIC50' and lower(act.standard_units) = 'nm' and act.standard_value <= 25000 and act.standard_relation in ('<', '=', '<=', '~'))     " +
                "        OR (act.standard_type = 'MIC50' and lower(act.standard_units) in ('ug.ml-1', 'ug ml-1') and act.standard_value <= 8 and act.standard_relation in ('<', '=', '<=', '~')) " +
                "        OR (act.standard_type = 'MIC95' and lower(act.standard_units) = 'um' and act.standard_value <= 100 and act.standard_relation in ('<', '=', '<=', '~'))     " +
                "        OR (act.standard_type = 'MIC95' and lower(act.standard_units) in ('ug.ml-1', 'ug ml-1') and act.standard_value <= 32 and act.standard_relation in ('<', '=', '<=', '~')) " +
                "        OR (act.standard_type = 'EC50' and lower(act.standard_units) = 'nm' and act.standard_value <= 25000 and act.standard_relation in ('<', '=', '<=', '~')) " +
                "        OR (act.standard_type = 'EC50' and lower(act.standard_units) in ('ug.ml-1', 'ug ml-1') and act.standard_value <= 8 and act.standard_relation in ('<', '=', '<=', '~')) " +
                "        OR (act.standard_type = 'ED50' and lower(act.standard_units) = 'nm' and act.standard_value <= 25000 and act.standard_relation in ('<', '=', '<=', '~')) " +
                "        OR (act.standard_type = 'ED50' and lower(act.standard_units) = 'mg.kg-1' and act.standard_value <= 250 and act.standard_relation in ('<', '=', '<=', '~')) " +
                "        OR (act.standard_type = 'MBC' and lower(act.standard_units) in ('ug.ml-1', 'ug ml-1') and act.standard_value <= 8 and act.standard_relation in ('<', '=', '<=', '~')) " +
                "        OR (act.standard_type = 'MBC90' and lower(act.standard_units) in ('ug.ml-1', 'ug ml-1') and act.standard_value <= 32 and act.standard_relation in ('<', '=', '<=', '~')) " +
                "        OR (act.standard_type = 'GI' and lower(act.standard_units) = '%'  and act.standard_value >= 70 and act.standard_value <= 100 and act.standard_relation in ('>', '=', '>=', '~')) " +
                "        OR (act.standard_type = 'GI' and lower(act.standard_units) = 'nm' and act.standard_value <= 25000 and act.standard_relation in ('<', '<=', '=', '~')) " +
                "        OR (act.standard_type = 'GI' and lower(act.standard_units) = 'um' and act.standard_value <= 25 and act.standard_relation in ('<', '<=', '=', '~')) " +
                "        OR (act.standard_type = 'GI50' and lower(act.standard_units) in ('ug.ml-1', 'ug ml-1') and act.standard_value <= 8 and act.standard_relation in ('<', '<=', '=', '~')) " +
                "        OR (act.standard_type = 'GI50' and lower(act.standard_units) = 'nm' and act.standard_value <= 25000 and act.standard_relation in ('<', '<=', '=', '~')) " +
                "      ) AND " +
                "      md.molecule_type = 'Small molecule' AND " +
                "      md.polymer_flag = 0;"
            )
            
            cursor.execute(mysqlCreateBindingTable)            

        except mysql.connector.Error as err:
            if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
                print("Something is wrong with your user name or password")
                raise
            elif err.errno == errorcode.ER_BAD_DB_ERROR:
                print("Database does not exist")
                raise
            else:
                print(err)
                raise
        else:
            cnx.close()


# In[ ]:



