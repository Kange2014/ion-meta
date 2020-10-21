#!/usr/bin/env python

from ete3 import NCBITaxa
import os

os.system("centrifuge-download -o taxonomy taxonomy") 
os.system("touch taxdb.sqlite")
ncbi_taxa = NCBITaxa("taxdb.sqlite")
os.system("mv taxdb.sqlite* taxonomy")

os.system("mkdir krona")
#os.system("cd krona") # not really works
os.chdir("./krona")
os.system("ktUpdateTaxonomy.sh --only-build --preserve ./")

#os.system('centrifuge-download -o library -P 8 -m -d "archaea,bacteria,viral,fungi,protozoa" refseq > seqid2taxid.map')


