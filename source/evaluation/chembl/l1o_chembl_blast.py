# import our matrix code
import matrix_chembl_ecfp4

# get similar proteins within the given e-value threshold for the given query protein
def getSimProteins(evalue, queryProtein):
  t1 = g.find("(a)-[e]->(b)").filter("e.evalue <= " + str(evalue) + " and a.id != b.id and e.src = '" + str(queryProtein) + "'").select('e.dst', 'e.evalue')
  t2 = g.find("(a)-[e]->(b)").filter("e.evalue <= " + str(evalue) + " and a.id != b.id and e.dst = '" + str(queryProtein) + "'").select(sqlf.col('e.src').alias('dst'), 'e.evalue')
  t3 = t1.union(t2).join(bindings, sqlf.col('protein') == sqlf.col('dst')).select('dst', 'ligand', 'evalue').orderBy(sqlf.desc('evalue'))
  t3.createOrReplaceTempView('t3')
  return spark.sql('select row_number() over (order by "evalue") as num, * from t3') # rank, dst protein, ligand, evalue score

# get the rank at which the hidden binding is uncovered
def getUncoverBindingRank(x, tc):
  tLigand = x[0]
  qProtein = x[1]
  result = getSimProteins(evalue, qProtein)
  result.createOrReplaceTempView('result')  
  rank = spark.sql('select num from result where ligand = "'+tLigand+'" limit 1').collect()
  if len(rank) > 0 :
    rank = rank[0]['num']
  else:
    rank = -1 # not found
  return (qProtein, tLigand, rank) 


# leave-one-out test
# ==================

for x in bindings:
    print(getUncoverBindingRank(x, 0.0001))
        
        
