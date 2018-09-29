# add required imports
from pyspark.sql.types import *
from pyspark.sql import Row
from graphframes import *
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from e3fp.config.params import default_params
from e3fp.fingerprint.fprint import Fingerprint, CountFingerprint
from e3fp.fingerprint.generate import *


# import data
rawdata = sc.textFile('/FileStore/tables/chembl.tsv')
data = rawdata.filter(lambda x: x != '').map(lambda x: x.split('\t'))

# build bindings dataframe
bindingsRdd = data.map(lambda x: Row(ligand=x[0], protein=x[2]))
bindings = bindingsRdd.toDF()
print('Number of active bindings: %s' % bindings.count())

# create molecules dictionary so that we can get molecule object by ligand name
molDict = {}
suppl = Chem.SDMolSupplier('/FileStore/tables/chembl.sdf', removeHs=False)
for mol in suppl:
    d.update({mol.GetProp('_Name').lower() : mol})

# build ligands sim matrix

# get fingerprint for a given molecule
def createFingerprint(chemblMolId):
    try:        
        mol = molDict[chemblMolId]
        fps = fprints_dict_from_mol(mol)
        ligandName = mol.GetProp('_Name').lower()
        # taken from https://e3fp.readthedocs.io/en/latest/overview.html
        fp = fps[5][0]
        fp = fp.fold().to_rdkit()
        return fp
    except:
        return None

# compute similarity between two fingerprints
def simFunction(fp1, fp2):
  return FingerprintSimilarity(fp1, fp2,) # default is Tanimoto

# create a ligands/fp rdd
ligands = data.map(lambda x: Row(ligand=x[0], fp=createFingerprint(x[0]))).distinct()
print('Number of unique ligands: %s' % ligands.count())

# create matrix DF
sim = ligands.cartesian(ligands).filter(lambda x: x[0][1] < x[1][1]).map(lambda x: Row(src=x[0][1], dst=x[1][1], sim=simFunction(x[0][0], x[1][0]))).toDF()

# create similarity graph
g = GraphFrame(data.map(lambda x: Row(id=x[0])).distinct().toDF(), sim)


