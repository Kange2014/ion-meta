cur_date=`date +%Y-%m-%d`
mkdir $cur_date
centrifuge-download -o ${cur_date} -P 10 -m -d "bacteria,viral,archaea" refseq > ${cur_date}/seqid2taxid.map
# for fungi or protozoa, download all refseq genomes because of no complete genomes for many pathogens
centrifuge-download -o ${cur_date} -P 10 -m -d "fungi,protozoa" -a 'Any' refseq >> ${cur_date}/seqid2taxid.map

