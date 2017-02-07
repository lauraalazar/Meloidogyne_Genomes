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
### Preqc
```
# Filter and quality trim reads
/sga-2014.03.05/bin/sga preprocess --pe-mode 1 file_R1.fastq.gz file_R2.fastq.gz > file_ilv.fastq

#Build index
/sga-2014.03.05/bin/sga index -a ropebwt --no-reverse -t 8 file_ilv.fastq

#Pre-assembly quality check
/sga-2014.03.05/bin/sga preqc -t 8 file_ilv.fastq > file.preqc

#Create a report
source /virt_env/python/matplotlib_1.4/bin/activate
/sga-2014.03.05/src/bin/sga-preqc-report.py file.preqc
```
## Filtering of reads (for the description of the protocol see https://blobtools.readme.io/)
### Preliminary assembly and mapping with clc
```
#Assembly with clc
module load /exports/modules/clc/5.0.0
clc_assembler -o assembly_clc.fasta -q file_R1.fastq.gz file_R2.fastq.gz

#Mapping with bowtie
/bwa index -a bwtsw assembly_clc.fasta 

/bwa mem 
-t 16 
-v assembly_clc.fasta \
file_R1.fastq.gz \
file_R2.fastq.gz | \
/samtools-1.1/samtools view -buS - | \
/samtools-1.1/samtools sort - assembly_map.bowtie.sorted.bam

```
###Blobplots
```
#Blast the assembly against the nt database
/blast/2.2.30+/blastn -task megablast 
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
## Trimming
```

```
## Assembly
```
#Contigs
/platanus/1.2.4/platanus assemble \
-t 48 \
-m 400 \
-o Sp \
-f \
Sp_lib1-trimmed-pair1.fastq \
Sp_lib1-trimmed-pair2.fastq \
Sp_lib2-trimmed-pair1.fastq \
Sp_lib2-trimmed-pair2.fastq

#Scaffolds
/platanus/1.2.4/platanus scaffold \
-t 48 \
-b Sp_contigBubble.fa \
-c Sp_contig.fa \
-o Sp \
-IP1  Sp_lib1-trimmed-pair1.fastq Sp_lib1-trimmed-pair2.fastq \
-IP2  Sp_lib2-trimmed-pair1.fastq Sp_lib2-trimmed-pair2.fastq

#Gap closing
/platanus/1.2.4/platanus gap_close \
-t 48 \
-c Sp_scaffold.fa \
-IP1  Sp_lib1-trimmed-pair1.fastq Sp_lib1-trimmed-pair2.fastq \
-IP2  Sp_lib2-trimmed-pair1.fastq Sp_lib2-trimmed-pair2.fastq

```

## Quality check of Assembly

## Annotation

## Final information
