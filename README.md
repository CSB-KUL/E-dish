# E-dish
A Repository for Functional Analyses of E coli Genomes on Desk 

# Usage
To map mutations to functional regions

$python edish.py

(Optional)

Enrichments analyses files can be mapped mutation files for summary of all the mutations which 
overlap the functional regions and their GO enrichments. (PS: for our analyses, we created enrichments
with string db (http://string-db.org/))

import edish

from edish import go

go("enrichment.txt", "ecoli_id.txt","mutation.txt")

# Contributor
Ahmed Arslan
