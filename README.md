**Encoding example**

Create a compressed tar archive:
```
tar -b1 -czvf info_to_code.tar.gz ./info_to_code/
```

Zero-padding to make the input a multiple of 512bytes
```
truncate -s2116608 ./info_to_code.tar.gz
```

Or download the original archive:
```
wget http://files.teamerlich.org/dna_fountain/dna-fountain-input-files.tar.gz
```

Actual encoding of data as DNA (output is a FASTA file):
```
python encode.py \
--file_in info_to_code.tar.gz \
--size 32 \
-m 3 \
--gc 0.05 \
--rs 2 \
--delta 0.001 \
--c_dist 0.025 \
--out info_to_code.tar.gz.dna \
--stop 72000
```

Add annealing sites:
```
cat info_to_code.tar.gz.dna | grep -v '>' |\
awk '{print "GTTCAGAGTTCTACAGTCCGACGATC"$0"TGGAATTCTCGGGTGCCAAGG"}' \
> info_to_code.tar.gz.dna_order
```

Output file is ready to order synthetic DNA.

**Decoding example**

Convert BCL to FASTQ using picard (https://github.com/broadinstitute/picard):
```
for i in {1101..1119} {2101..2119}; do
mkdir ~/Downloads/fountaincode/seq_data3/$i/;
done

for i in {1101..1119} {2101..2119}; do
java -jar ~/Downloads/picard-tools-2.5.0/picard.jar \
IlluminaBasecallsToFastq \
BASECALLS_DIR=./raw/19854859/Data/Intensities/BaseCalls/ \
LANE=1 \
OUTPUT_PREFIX=./seq_data3/$i/ \
RUN_BARCODE=19854859 \
MACHINE_NAME=M00911 \
READ_STRUCTURE=151T6M151T \
FIRST_TILE=$i \
TILE_LIMIT=1 \
FLOWCELL_BARCODE=AR4JF;
done
```

(Sequencing data are available at www.ebi.ac.uk/ena/data/view/PRJEB19305 and www.ebi.ac.uk/ena/data/view/PRJEB19307)

Read stitching using PEAR (Zhang J et al., Bioinformatics, 2014).
This step takes the 150nt reads and places them together to get back the full oligo.
```
for i in {1101..1119} {2101..2119}; do
pear -f ./$i.1.fastq -r ./$i.2.fastq -o $i.all.fastq;
done
```

Retain only fragments with 152nt (the original oligo size):
```
awk '(NR%4==2 && length($0)==152){print $0}' *.all.fastq.assembled.fastq > all.fastq.good
```

Sort to prioritize highly abundant reads:
```
sort -S4G all.fastq.good | uniq -c > all.fastq.good.sorted
gsed -r 's/^\s+//' all.fastq.good.sorted |\
sort -r -n -k1 -S4G > all.fastq.good.sorted.quantity
```

Exclude column 1 specifying the number of times a read was seen and exclude reads with N:
```
cut -f2 -d' ' all.fastq.good.sorted.quantity |\
grep -v 'N' > all.fastq.good.sorted.seq
# Decoding:
python ~/Downloads/fountaincode/receiver.py \
-f ./seq_data3/all.fastq.good.sorted.seq \
--header_size 4 \
--rs 2 \
--delta 0.001 \
--c_dist 0.025 \
-n 67088 \
-m 3 \
--gc 0.05 \
--max_hamming 0 \
--out decoder.out.bin
```

checksum verification:
```
md5 decoder.out.bin
```
expected output is 8651e90d3a013178b816b63fdbb94b9b
```
md5 info_to_code.tar.gz
```
expected output is 8651e90d3a013178b816b63fdbb94b9b
