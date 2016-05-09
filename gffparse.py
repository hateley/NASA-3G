#pulls transcript and gene info out of bac gff
#this is useful because Ensembl Bac doesn't support biomart
#just change the in_file

from BCBio import GFF

in_file = "Smar.gff"
out_file = open("parsed.txv", 'w')

in_handle = open(in_file)
for line in in_handle:
    if line.startswith('Chromosome'):
        parts = line.split('\t')
        feature = parts[2]
        if feature == "transcript":
            names = parts[8].split(';')
            tran = names[0].split(':')[1]
            gene = names[1].split(':')[1]
            check = names[2].split('=')[0]
            ID = names[2].split('=')[1]
            if check != 'Name':
                out_file.write(tran+'\t'+gene+'\t'+gene+'\n')
            else:
                out_file.write(tran+'\t'+gene+'\t'+ID+'\n')

in_handle.close()
out_file.close()
