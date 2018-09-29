# add required imports
from pyspark.sql.types import *
from pyspark.sql import Row
from graphframes import *
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem import MACCSkeys

# import data
rawdata = sc.textFile('/FileStore/tables/chembl.tsv')
data = rawdata.filter(lambda x: x != '').map(lambda x: x.split('\t'))

# build bindings dataframe
bindingsRdd = data.map(lambda x: Row(ligand=x[0], protein=x[2]))
bindings = bindingsRdd.toDF()
print('Number of active bindings: %s' % bindings.count())

# build ligands sim matrix

# get fingerprint for a given smiles
def createFingerprint(smiles):
    try:
        m = Chem.MolFromSmiles(smiles)
        if m == None:
            return None
        else:
            return MACCSkeys.GenMACCSKeys(m)
    except:
        return None

# compute similarity between two fingerprints
def simFunction(fp1, fp2):
    return FingerprintSimilarity(fp1, fp2,) # default is Tanimoto

# create a ligands/fp rdd
ligands = data.map(lambda x: Row(ligand=x[0], fp=createFingerprint(x[1]))).distinct()
print('Number of unique ligands: %s' % ligands.count())

# create matrix DF
sim = ligands.cartesian(ligands).filter(lambda x: x[0][1] < x[1][1]).map(lambda x: Row(src=x[0][1], dst=x[1][1], sim=simFunction(x[0][0], x[1][0]))).toDF()

# create similarity graph
g = GraphFrame(data.map(lambda x: Row(id=x[0])).distinct().toDF(), sim)