# Meloidogyne_Genomes

Quality check, *de novo* assembly and annotation of Meloidogyne species genomes

## Preliminary information
Species         | Strain | Location    | Insert size | Reads       | Size (bp)       | Exp. coverage |
| :---          | :---:  |  :---:      |  :---:      | :---:       |  :---:          |          ---: |
|*M. javanica*  | VW4    | California  | 300         | 62,075,861  | 7,635,287,416   | 100x          |
|*M. javanica*  | VW4    | California  | 500         | 123,728,247 | 15,168,259,765  | 200x          |
|*M. javanica*  | VW5    | California  | 350         | 193,072,088 | 24,134,011,000  | 320x          |
|*M. javanica*  | L57    |Morocco      |             | 32,669,417  | 4,083,677,125   | 27x           |
|*M. javanica*c | L15    | Thailand    | -           | 29,324,182  | 3,665,522,750   | 24x           |
|*M. javanica*c | L17    | Burkina Faso| -           | 31,332,441  | 3,916,555,125   | 26x           |
|*M. incognita* | W1     | ?           | 350         | 38,260,145  | 4,782,518,125   | 63x           |
|*M. incognita* | W1     | ?           | 550         | 30,290,198  | 3,786,274,750   | 50x           |
|*M. incognita* | VW6    | ?           | 350         | 28,840,610  | 3,605,076,250   | 48x           |
|*M. incognita* | VW6    | ?           | 550         | 25,746,808  | 3,218,351,000   | 42x           |
|*M. incognita* | HarC   | ?           | 350         | 26,844,521  | 3,355,565,125   | 44x           |
|*M. incognita* | HarC   | ?           | 550         | 35,340,761  | 4,417,595,125   | 58x           |
|*M. incognita* | 557R   | ?           | 550         | 62,745,198  | 7,843,149,750   | 104x          |
|*M. incognita* | L9     | Ivory Coast |             | 19,009,603  | 2,376,200,375   |  18x          |
|*M. incognita* | L19    | French West Indies | -    | 33,486,356  | 4,185,794,500   | 28x           |
|*M. incognita* | L27    | USA         | -           | 35,218,809  | 4,402,351,125   | 29x           |
|*M. incognita* | A14    | Lybia       | -           | 20,025,193  | 2,503,149,125   | 17x           |
|*M. arenaria*  | HarA   | California	 | 350         | 49,813,878  | 6,226,734,750   | 41x           |
|*M. arenaria*  | HarA   | California  | 550         | 46,643,017  | 9,656,831,750   | 64x           |
|*M. arenaria*  | L28    | French West Indies |  -   | 16,744,391  | 2,093,048,875   | 14x           |
|*M. arenaria*  | L32    | French West Indies | -    | 14,159,397  | 176,9924,625    | 11x           |
|*M. konaensis* | -      | ?                  | 350  | 55,647,288  | 6,955,911,000   | 46x           |
|*M. konaensis* | -      | ?                  | 550  | 52,813,037  | 6,601,629,625   | 44x           |
|*M. enterolobii*| L30   | ?                  | 350  | 143,672,079 |17,959,009,875	 | 120x          |
|*M. enterolobii*| L30   | ?                  | 550  | 100,032,455 | 12,504,056,875  | 83x           |
|*M. floridensis*|       | Florida            | | | | |
## Quality Control
### Fastqc
```
/fastqc_v0.11.2/fastqc file.fastq.gz
```
### Preqc (2014.03.05)
```
# Filter and quality trim reads
sga preprocess --pe-mode 1 file_R1.fastq.gz file_R2.fastq.gz > file_ilv.fastq

#Build index
sga index -a ropebwt --no-reverse -t 8 file_ilv.fastq

#Pre-assembly quality check
sga preqc -t 8 file_ilv.fastq > file.preqc

#Create a report
source /virt_env/python/matplotlib_1.4/bin/activate
sga-preqc-report.py file.preqc
```
## Filtering of reads (for protocol see https://blobtools.readme.io/)
### Preliminary assembly with clc (5.0.0), mapping bowtie (2-2.2.9) and sort with samtools (1.3.1)
```
#Assembly with clc
clc_assembler \
-o assembly_clc.fasta \
-q raw_R1.fastq.gz raw_R2.fastq.gz

#Mapping with bowtie
bowtie2 \
-x assembly_clc.fasta \
-1 raw_R1.fastq.gz \
-2 raw_R2.fastq.gz \
-p 16 \
| samtools-1.3.1/samtools view - -b > assembly_map.bowtie
bam

#Sort with samtools
samtools sort \
-n assembly_map.bowtie.bam \
-@ 16 \
-o assembly_map.bowtie.sorted
```
###Blobplots
```
#Blast (2.2.30+) the assembly against the nt database
blastn -task megablast 
-query assembly_clc.fasta \
-db nt \
-evalue 1e-5 \
-outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
-out assembly.25cul1.1e25.megablast.out \
-max_target_seqs 25 \
-culling_limit 2 \
-num_threads 8

#Create blobplot
/blobtools-light/makeblobs.py \
-a assembly_clc.fasta \
-bam assembly_map.bowtie.sorted.bam \
-blast assembly.25cul1.1e25.megablast.out \
-o assembly -taxdb /blast_db/

#Plot (matplotlib is needed)
source /exports/virt_env/python/matplotlib/bin/activate

/blobtools-light/plotblobs.py assembly.blobplot.txt
```
### Filter out reads
```
#Index the assembly with samtools
samtools faidx assembly_clc.fasta

#Find bacteria contigs in the assembly (e.g. Proteobacteria)
grep "tax=Proteobacteria" assembly.25cul1.1e25_lib2.megablast.blobplot.txt |  cut -f1 > assembly_proteobacteria.contigs

#Remove the the contigs from the indexed assembly
#grep -v -f assembly_proteobacteria.contigs assembly_clc.fasta.fai > assembly_clc_noprot.fasta.fai

#Get non-bacteria contigs with both reads mapped
samtools view \
-t assembly_clc_noprot.fasta.fai  \
-bS -F12 assembly_map_lib1.bowtie.sorted.bam \
> Sp_lib1.m_m.bam

#Pull the reads that mapped to non-bacteria contigs
/samtools-1.3.1/samtools view \
-uf64 Sp_lib1.m_m.bam \
| /samtools-1.3.1/samtools bam2fq - \
| gzip > Sp_lib1.reads1.fq.gz

```

## Trimming with Skewer (0.2.2-linux-x86_64)
```
skewer-0.2.2-linux-x86_64 \
-n \
-y adapters.fa \
-q 20 \
-l 25 \
-m pe \
o Sp_lib \
-t 32 \
Sp_lib1.reads1.fq.gz \
Sp_lib1.reads2.fq.gz

```
## Assembly with Platanus (1.2.4)
```
#Contigs
platanus assemble \
-t 48 \
-m 400 \
-o Sp \
-f \
Sp_lib1-trimmed-pair1.fastq \
Sp_lib1-trimmed-pair2.fastq \
Sp_lib2-trimmed-pair1.fastq \
Sp_lib2-trimmed-pair2.fastq

#Scaffolds
platanus scaffold \
-t 48 \
-b Sp_contigBubble.fa \
-c Sp_contig.fa \
-o Sp \
-IP1  Sp_lib1-trimmed-pair1.fastq Sp_lib1-trimmed-pair2.fastq \
-IP2  Sp_lib2-trimmed-pair1.fastq Sp_lib2-trimmed-pair2.fastq

#Gap closing
platanus gap_close \
-t 48 \
-c Sp_scaffold.fa \
-IP1  Sp_lib1-trimmed-pair1.fastq Sp_lib1-trimmed-pair2.fastq \
-IP2  Sp_lib2-trimmed-pair1.fastq Sp_lib2-trimmed-pair2.fastq

```
## Quality check of Assembly
###CEGMA
```
cegma --genome --threads 32 --output SP_assembly.platanus.fa
```
###BUSCO
```
python3 BUSCO.py \
-c 16 \
-i SP_assembly.platanus.fa \
-o sp.busco \
-m geno \
-l busco/busco-v2.0/nematoda_odb9/ \
-t sp_busco
```
###Map reads back to assembly with bwa
```
 bwa mem \
 Sp_assembly.platanus.fa \
 raw_R1.fastq raw_R2.fastq \
 | samtools view -buS -  | samtools sort - Sp.bwa.sorted
```

## Annotation
###Repeatmodeler 
```
RepeatModeler/BuildDatabase -name Sp \
Sp_assembly.platanus.fa
-engine ncbi \
```
###Repeatmasker(v4-0-3)
```
#Create a library with nematoda.repeatmasker from repeatmask, repeatmodeler output (consensi.fa.classified) and Nematoda.lib from Amir (ref?)

cat ematoda.repeatmasker consensi.fa.classified Nematoda.lib > nematoda_LibTE.repeats

#Run cdhit-est (4.6.4) with high cutoff to avoid redundancy
cd-hit-est -i \
nematoda_LibTE.repeats \
-o nematoda_LibTE.repeats.cdhitest \
-c 0.95 \
-d 0 \
-T 16 \
-M 3000

#Run repeatmasker
RepeatMasker \
Sp_assembly.platanus.fa \
-lib nematoda_LibTE.repeats.cdhitest \
-dir .
```

###SNAP (v2013-11-29)
```
maker/maker-2.31/bin/cegma2zff \
/cegma_output/Sp.cegma.gff \
Sp_assembly.platanus.fa #this creates files genome.ann, genome.dna

fathom genome.ann genome.dna \
-categorize 1000 \

fathom -export 1000 -plus uni.ann uni.dna

mkdir parameters

cd parameters

forge ../export.ann ../export.dna

cd ../

hmm-assembler.pl \
Sp_assembly.platanus.fa \
parameters > Sp.cegmasnap.hmm
```

###Genemark (es-2.3)
```
gm_es.pl  \
--BP OFF \
-max_nnn 500 \
-min_contig 10000 \
Sp_assembly.platanus.fa.masked #output of repeatmask

#the output of genemark is  mod.es in the mod folder as symbolik link
```
###Maker (2.311)
```
mpiexec -n 48 maker -fix_nucleotides
cd *.fa.maker.output/
gff3_merge -d *.fasta_master_datastore_index.log
```
###Augustus (3.0.1)
```
maker/bin/maker2zff Sp_assembly.platanus.fa.masked.all.gff

snap/zff2gff3.pl genome.ann | perl -plne 's/\t(\S+)$/\t\.\t$1/' > augustus.train.gff3

#/exports/software/snap/snap-2013-11-29/zff2gff3.pl genome.ann | perl -plne 's/\t(\S+)$/\t\.\t$1/' > augustus.train.gff3

autoAug.pl --species=Sp \
--genome=Sp_assembly.platanus.fa.masked \
--trainingset=augustus.train.gff3 \
--useexisting -v \
--singleCPU \
--optrounds=3 > augustus.out.txt
```

## Genome stats

| | *M. javanica* VW4 illumina | *M. javanica* VW4 pacbio | *M. incognita* W1 | *M. arenaria* HarA| *M. enterolobii* L30| *M. floridensis* |
| :--- | :---:   |    :---:       | :---:          |  :---:        |  :---:           |             ---: |
Scaffolds |  34,394   |   5,527        | 33,735        |  46,509          | 46,090  | |
Genome Span | 142,608,877 | 212,245,641 | 122,043,328  | 163,770,989 | 162,361,678 | |
Longest scaffold | 223,460 | 510,271 |  248,829 | 163,224 | 94,967 | |
N50 | 14,133 | 44,602 | 16,498 | 10,504 |  9,280 | |
GC | 29.6 | 29.3 | 0.299 |  29.5 |  29.7 | |
Mapped reads | 98.82% | 97.16% | 99.1% | 98.8% | 90.19% | |
18S | yes | yes | yes | yes | yes | |
mitochondria | yes | yes | yes | yes | yes | |
CEGMA completness | 90.32% | 93.95% | 82.66% | 91.13% | 81.45% | |
CEGMA average | 2.51 | 3.17 | 2.38 | 2.76 | 2.61 | |
Predicted Genes | 26.917 | 29.413 | 24.714 | 30.308 | 31.051 | |
Functional Annotated | 17.659 | | 15.938 | 20.813  |  | |
