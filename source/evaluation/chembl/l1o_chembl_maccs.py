# import our matrix code
import matrix_chembl_maccs

# get similar molecules within the given Tc threshold for the given query molecule
def getSimMolecules(tc, queryMolecule):
  t1 = g.find("(a)-[e]->(b)").filter("e.sim >= " + str(tc) + " and a.id != b.id and e.src = '" + queryMolecule + "'").select('e.dst', 'e.sim')
  t2 = g.find("(a)-[e]->(b)").filter("e.sim >= " + str(tc) + " and a.id != b.id and e.dst = '" + queryMolecule + "'").select(sqlf.col('e.src').alias('dst'), 'e.sim')
  t3 = t1.union(t2).join(bindings, sqlf.col('ligand') == sqlf.col('dst')).select('dst', 'protein', 'sim').orderBy(sqlf.desc('sim'))
  t3.createOrReplaceTempView('t3')
  return spark.sql('select row_number() over (order by "sim") as num, * from t3') # rank, dst ligand, target protein, sim score

# get the rank at which the hidden binding is uncovered
def getUncoverBindingRank(x, tc):
  qLigand = x[0]
  tProtein = x[1]
  result = getSimMolecules(tc, qLigand)
  result.createOrReplaceTempView('result')  
  rank = spark.sql('select num from result where protein = "'+tProtein+'" limit 1').collect()
  if len(rank) > 0 :
    rank = rank[0]['num']
  else:
    rank = -1 # not found
  return (qLigand, tProtein, rank) 


# leave-one-out test
# ==================

for tc in range(1, 10):
    for x in bindings:
        print(getUncoverBindingRank(x, tc / 10.0))
        
        
