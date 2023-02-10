import sys
from neo4j import GraphDatabase
from tqdm import tqdm

def namestr(obj, namespace):
    return [name for name in namespace if namespace[name] is obj]

cond=sys.argv[1]

#index
RNA='CREATE INDEX index_uno IF NOT EXISTS FOR (n:mRNA) ON (n.GENEXP, n.HYPO, n.HYPER)'
CH3='CREATE INDEX index_due IF NOT EXISTS FOR (n:CH3) ON (n.GENEXP, n.HYPO, n.HYPER)'
MP='CREATE INDEX index_tre IF NOT EXISTS FOR (n:MP) ON (n.GENEXP, n.HYPO, n.HYPER)'

index=(RNA,CH3,MP)


#WHERE='MATCH p=(g:mRNA{CHR:"20", GENEXP_vM:"UP_REG"})-[*..2]->(g1:CH3) WHERE g1.HYPO_vM IS NOT NULL AND  all(n IN nodes(p)[1..-1] WHERE n.GENEXP_vM = "UP_REG") RETURN DISTINCT nodes(p),relationships(p)'

#lineX='MATCH (g:mRNA{GENEXP_vM:"UP_REG"}),(g1:CH3{SYMBOL:"IWS1"}), p=(g1)-[*..3]->(g) RETURN DISTINCT nodes(p),relationships(p)'
line1='MATCH (g:mRNA{GENEXP:"UP_REG"}),(g1:CH3{SYMBOL:"IWS1"}), p=((g1)-[*..2]->(g)) RETURN DISTINCT nodes(p),relationships(p)'
line2='MATCH (g:mRNA{GENEXP:"UP_REG"}), (g1:CH3{SYMBOL:"FBXO31"}), p=((g1)-[*..2]->(g)) WHERE all(n IN nodes(p)[1..-1] WHERE n.GENEXP IS NOT NULL) RETURN DISTINCT nodes(p),relationships(p)'
line3='MATCH (g:mRNA{isCancerG:"YES"}), (g1:CH3), p=((g1)-[*..2]->(g)) WHERE g1.HYPO IS NOT NULL AND all(n IN nodes(p)[1..-1] WHERE n.GENEXP IS NOT NULL) RETURN DISTINCT nodes(p),relationships(p)'


query=(line1, line2, line3)

driver=GraphDatabase.driver("bolt://localhost:7687",auth=("neo4j", "antonio"), encrypted=False) 
session=driver.session()

for i in index:
    result = session.run(i)
    print(result)

for q in query:
    extract= namestr(q,globals())
    name= extract[0]
    message="Query %s on %s" %(name, cond)
    print(message)
    #driver = GraphDatabase.driver("bolt://localhost:7687", auth=("neo4j", "antonio"), encrypted=False)
    #session = driver.session()
    result = session.run(q)
    relationships = list()
    nodes = list()
    print("Recording edges and node...")
    for record in tqdm(iter(result)):
        record_data = record.data()
        for path in record_data["relationships(p)"]:
            relationships.append(path[0]["SEQ"] + "," + path[1] + "," + path[2]["SEQ"])
        for subrecord in record_data["nodes(p)"]:
            nodes.append(subrecord["SEQ"] + "," + subrecord["ENTREZID"] + "," + subrecord["SYMBOL"])

    name_r='%s_%s_edges.csv' %(name, cond)
    relationships_csv = open(name_r, "w")
    relationships_csv.write("SRG,TYPE,TRG\n")
    print("Create csv file for edges")
    for record in relationships:
        relationships_csv.write(record + "\n")

    relationships_csv.close()

#nodes_csv = open("WHERE_gain_nodes.csv", "w")
    name_n='%s_%s_node.csv' %(name, cond)
    nodes_csv = open(name_n, "w")
#nodes_csv.write("ENTREZ,SYMBOL,CH3_STATE_vM,CH3_STATE_GvD,CH3_POSITION,logFC_vM,GENEXP_vM,logFC_GvD,GENEXP_GvD,PROTEXP,MUTATION,TRANFAC,MAP\n")
    nodes_csv.write("SEQ,ENTREZID,SYMBOL\n")
    print("Create csv file for nodes")
    for record in nodes:
        nodes_csv.write(record + "\n")

    nodes_csv.close()

#print(record_data["relationships(p)"])

#nodes.append(subrecord["ENTREZ"] + "," + subrecord["SYMBOL"] + "," + subrecord["CH3_STATE_vM"] + "," + subrecord["CH3_STATE_GvD"]
#+ "," + subrecord["CH3_POSITION"]+ "," + subrecord["logFC_vM"] + "," + subrecord["GENEXP_vM"] + "," + subrecord["logFC_GvD"]
#+ "," + subrecord["GENEXP_GvD"]+ "," + subrecord["PROTEXP"] + "," + subrecord["MUTATION"] + "," + subrecord["TRANFAC"] + "," + subrecord["MAP"])
