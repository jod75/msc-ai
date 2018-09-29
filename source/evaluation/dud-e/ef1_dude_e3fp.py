# import our matrix code
import matrix_dude_e3fp

# get similar molecules within the given Tc threshold for the given query molecule
def getSimMolecules(tc, queryMolecule):
  t1 = g.find("(a)-[e]->(b)").filter("e.sim >= " + str(tc) + " and a.id != b.id and e.src = '" + queryMolecule + "'").select('e.dst', 'e.sim')
  t2 = g.find("(a)-[e]->(b)").filter("e.sim >= " + str(tc) + " and a.id != b.id and e.dst = '" + queryMolecule + "'").select(sqlf.col('e.src').alias('dst'), 'e.sim')
  t3 = t1.union(t2).join(bindings, sqlf.col('ligand') == sqlf.col('dst')).select('dst', 'protein', 'sim').orderBy(sqlf.desc('sim'))
  t3.createOrReplaceTempView('t3')
  return spark.sql('select row_number() over (order by "sim") as num, * from t3') # rank, dst ligand, target protein, sim score

# get enrichment factor
def getEnrichment(x, tc, total , top):
  qLigand = x[0]
  tProtein = x[1]
  result = getSimMolecules(tc, qLigand)
  result.createOrReplaceTempView('result')  
  bindings.createOrReplaceTempView('bindings')
  actives = spark.sql('select count(*) as c from bindings where protein = "'+tProtein+'"').collect()[0]['c']
  found = spark.sql('select count(*) as c from result, bindings where result.protein = "'+tProtein+'" and result.dst = bindings.ligand').collect()[0]['c']
  ya =  found * 1.0 / top
  return (qLigand, tProtein, tc, ya, ya / (actives * 1.0 / total)) 

# EF1 computation
# ===============

total = ligands.count()
top = int(round(total * 0.1, 0))
print('Total number of ligands in data set: %s' % (total))
print('Top 1 percent: %s' % top)

for tc in range(1, 10):
    for x in bindings.collect():
        print(getEnrichment(x, tc / 10.0, total, top))