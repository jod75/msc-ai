# import our matrix code
import matrix_chembl_protein

# get similar proteins within the given e-value threshold for the given query protein
def getSimProteins(evalue, queryProtein):
  t1 = g.find("(a)-[e]->(b)").filter("e.evalue <= " + str(evalue) + " and a.id != b.id and e.src = '" + str(queryProtein) + "'").select('e.dst', 'e.evalue')
  t2 = g.find("(a)-[e]->(b)").filter("e.evalue <= " + str(evalue) + " and a.id != b.id and e.dst = '" + str(queryProtein) + "'").select(sqlf.col('e.src').alias('dst'), 'e.evalue')
  t3 = t1.union(t2).join(bindings, sqlf.col('protein') == sqlf.col('dst')).select('dst', 'ligand', 'evalue').orderBy(sqlf.desc('evalue'))
  t3.createOrReplaceTempView('t3')
  return spark.sql('select row_number() over (order by "evalue") as num, * from t3') # rank, dst protein, ligand, evalue score

# get enrichment factor
def getEnrichment(x, evalue, total , top):
  tLigand = x[0]
  qProtein = x[1]
  result = getSimProteins(evalue, qProtein)
  result.createOrReplaceTempView('result')  
  bindings.createOrReplaceTempView('bindings')
  actives = spark.sql('select count(*) as c from bindings where ligand = "'+tLigand+'"').collect()[0]['c']
  found = spark.sql('select count(*) as c from result, bindings where result.ligand = "'+tLigand+'" and result.dst = bindings.protein').collect()[0]['c']
  ya =  found * 1.0 / top
  return (qProtein, tLigand, evalue, ya, ya / (actives * 1.0 / total)) 

# import chembl data
# ==================
# import data
rawdata = sc.textFile('/FileStore/tables/chembl.tsv')
data = rawdata.filter(lambda x: x != '').map(lambda x: x.split('\t'))

# build bindings dataframe
bindingsRdd = data.map(lambda x: Row(ligand=x[0], protein=x[2]))
bindings = bindingsRdd.toDF()
print('Number of active bindings: %s' % bindings.count())
  
# EF1 computation
# ===============

total = nodes.count()
top = int(round(total * 0.1, 0))
print('Total number of proteins in data set: %s' % (total))
print('Top 1 percent: %s' % top)

for x in bindings.collect():
    print(getEnrichment(x, 0.0001, total, top))