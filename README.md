# Discovery of Medicinal Molecules based on Similarity Networks

Source code for MSc AI distinction dissertation submitted at the University of Malta: "Discovery of Medicinal Molecules based on Similarity Networks"

At a molecular level, ligands work by binding to a therapeutic target of interest, usually a protein, to solicit a response or block its function.  Many modern drug discovery projects start with a computational search for these "active" molecules, called Virtual Screening (VS).  One of the approaches in these computational searches is Ligand-Based Virtual Screening (LBVS), based on the Similarity Property Principle (SPP) that states that molecules with similar chemical composition should exhibit similar behaviour. This similarity search is performed over a library with millions of molecules. But, how would 'similarity' be defined? Many ligand and protein similarity metrics exist, and the success of any particular technique depends on the target of interest. Big data techniques are adequate for such problems involving complex computations and huge data sets. In this dissertation we present a large-scale (LBVS) approach leveraging the parallelism and scalability of Apache Spark.  It uses GraphFrames, a package for Apache Spark which provides graph abstraction, to build multiple similarity networks of proteins and ligands using different molecule representations.  To the author's knowledge, this is the first work that uses GraphFrames to tackle VS problems.  Using ligand and protein networks and bridging them via known bindings in ChEMBL repository, the aim of this study is to be able to discover new putative ligand binders for a given protein target, and to find other proteins which are likely to interact with a given ligand. A web portal for visualising results was built. We evaluate with success our approach by comparing its enrichment and recall with similar tools, using DUD-E data set. A final test, which was carried out in collaboration with the Department of Pharmacy at the same university, we correctly rejected the farnesyl pyrophosphate synthase (FPPS) and suggested oestrogen receptor (ER) as a putative protein targets for maltanedienol that were also confirmed via _in vitro_ testing.

### Repo structure
<pre>
|
+---maltanedienol
|       maltanedienol_chembl_maccs_tc_0_8_results.csv  : VS results
|       maltanedienol_smiles.txt          : SMILES representation for maltanedienol
|       
+---source
|   +---evaluation                        : Evaluation python scripts
|   |   +---chembl
|   |   \---dud-e
|   |           
|   \---matrices                          : Pyspark scripts to create matrices
|       +---chembl
|       \---dud-e
|               
\---tool : web tool
    +---backend                           : These files must be installed on Spark's worker nodes
    |                                     :  Update chemblhelper.py with your ChEMBL MySql server details
    |                                     :   and engine.py to reflect your folder structure and hadoop urls
    |                                     :  Run start-web-server.sh bash script to start backend server
    \---www                               : These file must be installed in a web server, e.g. Apache 
                                          :  Update www/js/ics5200.js with your backend hadoop url 
</pre>
