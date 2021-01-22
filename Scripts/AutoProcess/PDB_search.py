import os
import urllib2
import requests
import datetime
from datetime import timedelta
import xml.etree.ElementTree as ET

today = datetime.date.today()
last_week = today - timedelta(days=7)
today = str(today)
last_week = str(last_week)
#only for first use
last_week = "1900-01-01"
#last_week = "2019-06-11"
#today = "2019-06-20"

url = 'http://www.rcsb.org/pdb/rest/search'
queryText = """
<orgPdbCompositeQuery version="1.0">
 <queryRefinement>
  <queryRefinementLevel>0</queryRefinementLevel>
  <orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.ExpTypeQuery</queryType>
    <description>Experimental Method is X-RAY</description>
    <mvStructure.expMethod.value>X-RAY</mvStructure.expMethod.value>
  </orgPdbQuery>
 </queryRefinement>
 <queryRefinement>
  <queryRefinementLevel>1</queryRefinementLevel>
  <conjunctionType>and</conjunctionType>
  <orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.ResolutionQuery</queryType>
    <description>Resolution is 3.0 or less</description>
    <refine.ls_d_res_high.comparator>between</refine.ls_d_res_high.comparator>
    <refine.ls_d_res_high.max>3.0</refine.ls_d_res_high.max>
  </orgPdbQuery>
 </queryRefinement>
 <queryRefinement>
  <queryRefinementLevel>2</queryRefinementLevel>
  <conjunctionType>and</conjunctionType>
  <orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.ChemFormulaWeightQuery</queryType>
    <description>Molecular Weight (Chemical component): Molecular weight min is 300</description>
    <molecularWeightComparator>between</molecularWeightComparator>
    <molecularWeightMin>300</molecularWeightMin>
  </orgPdbQuery>
 </queryRefinement>
 <queryRefinement>
  <queryRefinementLevel>3</queryRefinementLevel>
  <conjunctionType>and</conjunctionType>
  <orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.ChainTypeQuery</queryType>
    <description>Chain Type: there is a Protein chain but not any DNA or RNA or Hybrid</description>
    <containsProtein>Y</containsProtein>
    <containsDna>N</containsDna>
    <containsRna>N</containsRna>
    <containsHybrid>N</containsHybrid>
  </orgPdbQuery>
 </queryRefinement>
 <queryRefinement>
  <queryRefinementLevel>4</queryRefinementLevel>
  <conjunctionType>and</conjunctionType>
  <orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.ReleaseDateQuery</queryType>
    <description>Released on """ + last_week + """ and later</description>
    <queryId>null</queryId>
    <pdbx_audit_revision_history.revision_date.comparator>between</pdbx_audit_revision_history.revision_date.comparator>
    <pdbx_audit_revision_history.revision_date.min>""" + last_week + """</pdbx_audit_revision_history.revision_date.min>
    <pdbx_audit_revision_history.ordinal.comparator>=</pdbx_audit_revision_history.ordinal.comparator>
    <pdbx_audit_revision_history.ordinal.value>1</pdbx_audit_revision_history.ordinal.value>
  </orgPdbQuery>
 </queryRefinement>
</orgPdbCompositeQuery>
"""

#print "querying PDB...\n"

req = urllib2.Request(url, data=queryText)
f = urllib2.urlopen(req)
result = f.read()
if result:
    print "Found number of PDB entries:", result.count('\n')
    result = result.split('\n')
    smile_ligands = []
    pdbs = result[:-1]
    for pdb in pdbs:
        link = "http://www.rcsb.org/pdb/rest/ligandInfo?structureId=" + pdb
        r = requests.get(link, stream=True)
        xml = ""
        for line in r.iter_lines():
            xml = xml + line + '\n'
        root = ET.fromstring(xml)
        for ligand in root.iter('ligand'):
            if float(ligand.get('molecularWeight')) >= 300:
                smiles = ligand.find('smiles').text
                if 'Fe' in smiles or 'Hg' in smiles:
                    continue
                ligID = ligand.get('chemicalID')
                new_ligand = smiles + '\t' + ligID + '\n'
                if not new_ligand in smile_ligands:
                    smile_ligands.append(new_ligand)
    os.mkdir(today)
    os.chdir(today)
    with open('ligands.smi', 'w') as f:
        for ligand in smile_ligands:
            f.write(ligand)
    with open('pdb.list', 'w') as f:
        for pdb in pdbs:
            f.write(pdb + '\n')
    os.chdir('../')

print today

#else:
#    print "Failed to retrieve results" 
