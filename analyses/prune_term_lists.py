from pyhpo import Ontology
_=Ontology()

import pandas as pd

phenotype_data = pd.read_csv('input.csv')
'''
"Terms" column should have list of HPO terms in the form 'Term1; Term2; Term3'
'''
phenotype_data

prenatal_terms = ['HP:0012188','HP:0008071','HP:0009800','HP:0030244','HP:0100622','HP:0011438','HP:0100603','HP:0011436','HP:0001511','HP:0001562','HP:0001561','HP:0001558','HP:0010519','HP:0001787','HP:0001622','HP:0001518','HP:0001520','HP:0003561','HP:0003517','HP:0011451','HP:0004488','HP:0002643','HP:0006579','HP:0002033','HP:0001998','HP:0040187',
                  'HP:0011410', 'HP:0030364', 'HP:0001787', 'HP:0030369','HP:0001788', 'HP:0100603','HP:0001194'
]
print(len(prenatal_terms))
for term in prenatal_terms:
    print(Ontology.get_hpo_object(term))

n=0
removed_terms = []
no_prenatal = []
changed_status=[]
separator='; '
for i, row in phenotype_data.iterrows():
    terms = row['Terms'].split('; ')
    prenatal_assigned_terms = list(set(terms) & set(prenatal_terms))
    if len(prenatal_assigned_terms) != 0:
        n+=1
        changed_status.append('Changed')
        removed_terms += prenatal_assigned_terms
        terms_no_prenatal = separator.join(list(set(terms) - set(prenatal_terms)))
        
    else:
        terms_no_prenatal = separator.join(terms[:])
        changed_status.append('Unchanged')
    no_prenatal.append(terms_no_prenatal)

print(n, 'patients had prenatal terms removed')

print(len(set(removed_terms)))
for term in set(removed_terms):
    print(Ontology.get_hpo_object(term))

phenotype_data['Pruned_Terms'] = no_prenatal
phenotype_data['Status'] = changed_status


## save file
phenotype_data.to_csv('pruned_terms.csv', index=None)
