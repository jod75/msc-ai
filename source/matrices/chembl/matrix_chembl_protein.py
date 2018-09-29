# Steps:
# 
# 1. Download fasta file (from ChEMBL)
# 2. Create a local BLAST database with the file downloaded in (1)
# 3. Distribute local BLAST database to all nodes in cluster
# 4. For each protein sequence, caluclate similarities and store them in a DataFrame
# 5. Create a matrix from the DataFrame obtained in (4)

# ## Preperation

# imports
from pyspark import SparkContext, SparkConf
from hdfs import InsecureClient
import datetime
import subprocess
import shlex
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


# initialize Spark
if not 'sc' in globals():
    conf = SparkConf().setMaster('yarn')
    sc = SparkContext(conf=conf)
    
sc.setLogLevel("INFO")


# global variables
distFilename = "/home/hduser/Lab/blastdistfiles.sh"
fastaFile = "/home/hduser/Lab/prosim/proteins.fasta"
blastScript = "/home/hduser/proteins/blast.sh"


# ## Create local BLAST database
print("Creating local BLAST database")
process = subprocess.Popen(shlex.split("makeblastdb -in %s -parse_seqids -dbtype prot" % os.path.basename(fastaFile)),                           
                           stdout = subprocess.PIPE,
                           stderr = subprocess.PIPE,  
                           cwd = os.path.dirname(fastaFile),
                           shell = False)

out, err = process.communicate()
print("Output: %s" % out)
print("Error: %s" % err)


# ## Distribute BLAST database and helper Shell scripts
# run blastdistfiles.sh to distribute local BLAST db and blast.sh
print("Distributing BLAST database to all worker nodes")
process = subprocess.Popen(distFilename)
out, err = process.communicate()
print("Output: %s" % out)
print("Error: %s" % err)


# ## Run BLAST for each FASTA sequence
x = sc.parallelize(SeqIO.parse(fastaFile, "fasta"))
y = x.map(lambda q: ">%s\r\n%s" % (q.description, str(q.seq)))
result = y.pipe(blastScript)

# result.saveAsTextFile("/user/hduser/ics5200/proteinsresults")

nodesDf = x.map(lambda x: Row(id=x[0], sequence=x[1])).toDF()
edgesDf = y.map(lambda x: Row(src=x[0], dst=x[1], score=float(x[2]), evalue=float(x[3]))).toDF()

gsimGraph = GraphFrame(nodesDf, edgesDf)
