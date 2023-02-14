import os
import re
from snakemake.utils import R

configfile: "config.yaml"

#Location of raw files
raw_seq_folder = config["raw_dir"]
qc_seq_folder = config["qc_dir"]

#Cleanup file names and make sample list
file_list=[f for f in sorted(os.listdir(raw_seq_folder)) if (str(f))[-9:] == ".fastq.gz"]
read_dirs = ["R1","R2"]
samples = []
for file in file_list:
    dir_and_file = raw_seq_folder + "/" + file
    if re.search("_S[0-9]_R1_[0-9]*.fastq.gz",file) and re.sub("_S[0-9]_R1_[0-9]*.fastq.gz", "",file) not in samples:
      samples.append(re.sub("_S[0-9]_R1_[0-9]*.fastq.gz", "",file))
    if re.search("_R1.fastq.gz",file) and re.sub("_R1.fastq.gz", "",file) not in samples:
      samples.append(re.sub("_R1.fastq.gz", "",file))
    if re.search("_R1_",dir_and_file):
       new_file = re.sub("_S[0-9]_R1_[0-9]*.fastq.gz", "_R1.fastq.gz",dir_and_file)
       os.rename(dir_and_file,new_file)
    if re.search("_R2_",dir_and_file):
       new_file = re.sub("_S[0-9]_R2_[0-9]*.fastq.gz", "_R2.fastq.gz",dir_and_file)
       os.rename(dir_and_file,new_file)


localrules: all, clean, foam_concat, make_rulegraph

rule all:
    input:
        #'data/assembly_metaspades/{sample}'
        expand("{raw_folder}/fastqc_reports/{sample}_{read_dir}_fastqc.html", raw_folder=raw_seq_folder, sample=samples, read_dir=read_dirs),
        expand("{qc_folder}/fastqc_reports/{sample}_{read_dir}_fastqc.html", qc_folder=qc_seq_folder, sample=samples, read_dir=read_dirs),
        #expand("data/assembly_metaspades/{sample}",sample=samples)
        #expand("data/shogun/input_fastqs/{sample}_{read_dir}.fasta", sample=samples, read_dir=read_dirs)
        "data/shogun/combined_samples.fasta",
        "results/sequence_quality/",
        expand("data/assembly_megahit/{sample}", sample=samples),
        expand("data/assembly_metaspades/{sample}", sample=samples),
        expand("data/prodigal/{sample}",sample=samples),
        #expand("data/FOAM/{sample}_{dir}_sixframe_summary.txt",sample=samples,dir=read_dirs),
        "data/shogun/bt2/",
        "data/shogun/burst/",
        expand("data/kraken2/{sample}_report.txt", sample = samples),
        expand("data/kraken2/{sample}_brackenReport.txt", sample = samples),
        expand("data/kraken2_euk/{sample}_brackenReport.txt", sample = samples),
        expand("data/kraken2/{sample}_brackenMpa.txt", sample = samples),
        expand("data/humann3/{sample}/", sample = samples),
        expand("data/FOAM/{sample}_combined.rds", sample = samples)

rule paper_outputs:
    input:
        expand("{raw_folder}/fastqc_reports/S3_Fallen_{read_dir}_fastqc.html", raw_folder=raw_seq_folder, read_dir=read_dirs),
        expand("{qc_folder}/fastqc_reports/S3_Fallen_{read_dir}_fastqc.html", qc_folder=qc_seq_folder, read_dir=read_dirs),
        "results/sequence_quality/",
        "data/kraken2/S3_Fallen_report.txt",
        "data/kraken2/S3_Fallen_brackenReport.txt",
        "data/kraken2_euk/S3_Fallen_brackenReport.txt",
        "data/kraken2/S3_Fallen_brackenMpa.txt",
        "data/humann3/S3_Fallen/",
        "data/assembly_megahit/S3_Fallen",
        "data/assembly_metaspades/S3_Fallen/contigs.fasta",
        "data/prodigal/S3_Fallen",
        "data/BLAST/ectoine_genes_BLAST.blastn",
        "data/contig_coverage.tsv",
        "rulegraph.png",
        "data/human_genome_coverage.tsv",
        "data/mmseqs/concreteMetaG_report"


rule tax_profile:
    input:
        expand("data/kraken2/{sample}_report.txt", sample = samples),
        expand("data/kraken2/{sample}_brackenReport.txt", sample = samples),
        expand("data/kraken2_euk/{sample}_brackenReport.txt", sample = samples),
        expand("data/kraken2/{sample}_brackenMpa.txt", sample = samples)

rule assemble:
    input:
        expand("data/assembly_metaspades/{sample}", sample = samples),
        expand("data/assembly_megahit/{sample}", sample = samples),
        expand("data/assembly_megahit/{sample}/final.contigs.fa",sample = samples)

rule make_rulegraph:
    output:
        "rulegraph.pdf",
        "rulegraph.png"
    shell:
        """
        snakemake paper_outputs --rulegraph --dry-run | dot -Tpdf > rulegraph.pdf
        snakemake paper_outputs --rulegraph --dry-run | dot -Tpng > rulegraph.png
        """

rule fastqc_raw:
    input:
        #fastqs=expand("{raw_folder}/{sample}_{read_dir}.fastq.gz",raw_folder = raw_seq_folder, sample=samples, read_dir=read_dirs)
        fastqs = "raw_data/{file}.fastq.gz"
    output:
        report = "raw_data/fastqc_reports/{file}_fastqc.html"
        #expand("{raw_folder}/fastqc_reports/{sample}_{read_dir}_fastqc.html", raw_folder=raw_seq_folder, sample=samples, read_dir=read_dirs)
        #"raw_data/fastqc_reports/{fastqc_report}.html"
    conda:
          "code/master_env.yaml"
    shell:
        """
        mkdir -p raw_data/fastqc_reports
        fastqc -o raw_data/fastqc_reports -t {threads} {input.fastqs}
        """

rule trim_galore:
    input:
        fwd_fastq="raw_data/{sample}_R1.fastq.gz",
        rev_fastq="raw_data/{sample}_R2.fastq.gz"
    output:
        f_fq="data/qc_sequence_files/{sample}_R1.fastq.gz",
        r_fq="data/qc_sequence_files/{sample}_R2.fastq.gz"
    conda:
          "code/master_env.yaml"
    shell:
        """
        mkdir -p data/qc_sequence_files
        trim_galore --cores {threads} --paired -o data/qc_sequence_files --basename {wildcards.sample} {input.fwd_fastq} {input.rev_fastq}
        mv data/qc_sequence_files/{wildcards.sample}_val_1.fq.gz {output.f_fq}
        mv data/qc_sequence_files/{wildcards.sample}_val_2.fq.gz {output.r_fq}
        """

#This hasn't been tested/run on the data as of 5/19/21 but shoould run and ideally should be included.
# rule phix_removal:
#     input:
#         phix="data/reference/phiX/phiX.fasta",
#         f_fq="data/qc_sequence_files/{sample}_val_1.fq.gz",
#         r_fq="data/qc_sequence_files/{sample}_val_2.fq.gz"
#     output:
#         f_fq="data/qc_sequence_files/{sample}_R1.fastq.gz",
#         r_fq="data/qc_sequence_files/{sample}_R2.fastq.gz"
#     conda: "code/bwa_env.yaml"
#     shell:
#         """
#         mkdir -p data/phix_removal
#
#         #align qc'd reads to phix genome with bwa mem
#         bwa mem {input.phix} {input.f_fq} {input.r_fq} -t {resources.cpus} > data/phix_removal/{wildcards.sample}.sam
#
#         #convert sam to bam
#         samtools view -bS data/phix_removal/{wildcards.sample}.sam > data/phix_removal/{wildcards.sample}.bam
#
#         #Extract unmapped reads
#         samtools view -b -f 12 -F 256 data/phix_removal/{wildcards.sample}.bam > data/phix_removal/{wildcards.sample}_unmapped.bam
#
#         #Sort unmapped reads with 2 threads
#         samtools sort -n -m 5G -@ 2 data/phix_removal/{wildcards.sample}_unmapped.bam -o data/phix_removal/{wildcards.sample}_unmapped_sorted.bam
#
#         samtools fastq -@ {resources.cpus} data/phix_removal/{wildcards.sample}_unmapped_sorted.bam \
#             -1 {output.f_fq} \
#             -2 {output.r_fq} \
#             -0 /dev/null -s /dev/null -n
#
#         #cleanup tmp files
#         rm data/phix_removal/{wildcards.sample}*
#         """
#


rule fastqc_trimmed:
    input:
        fastqs="data/qc_sequence_files/{sample}_{read_dir}.fastq.gz"
    output:
        report = "data/qc_sequence_files/fastqc_reports/{sample}_{read_dir}_fastqc.html"
    conda:
          "code/master_env.yaml"
    shell:
        """
        mkdir -p data/qc_sequence_files/fastqc_reports
        fastqc -o data/qc_sequence_files/fastqc_reports -t {threads} {input.fastqs}
        """

rule multiqc:
    input:
        reports = expand("{qc_folder}/fastqc_reports/{sample}_{read_dir}_fastqc.html", qc_folder=qc_seq_folder, sample=samples, read_dir=read_dirs),
        raw_fastqc = expand("raw_data/fastqc_reports/{sample}_{read_dir}_fastqc.html", sample=samples, read_dir=read_dirs)
    output: directory("results/sequence_quality")
    conda: "code/multiqc_env.yaml"
    shell:
        """
        multiqc -d "raw_data/fastqc_reports" "data/qc_sequence_files/fastqc_reports" -o {output}
        """

rule megahit_assemble:
    input:
      f_reads = "data/qc_sequence_files/{sample}_R1.fastq.gz",
      r_reads = "data/qc_sequence_files/{sample}_R2.fastq.gz"
    output:
        results_dir = directory("data/assembly_megahit/{sample}"),
        contigs = "data/assembly_megahit/{sample}/final.contigs.fa"
    conda: "code/master_env.yaml"
    resources: cpus=16, mem_mb=200000, time_min=1440, mem_bytes = 20000000000
    shell:
        """
        rm -r data/assembly_megahit/{wildcards.sample} #megahit doesn't like the folder to already exist, but good to let snakemake have directory as output

        megahit --presets meta-sensitive -m {resources.mem_bytes} -1 {input.f_reads} -2 {input.r_reads} -o {output.results_dir}
        """

rule metaspades_assemble:
    input:
          f_reads = "data/qc_sequence_files/{sample}_R1.fastq.gz",
          r_reads = "data/qc_sequence_files/{sample}_R2.fastq.gz"
    output:
        directory("data/assembly_metaspades/{sample}"),
        "data/assembly_metaspades/{sample}/contigs.fasta"
    conda: "code/master_env.yaml"
    resources: cpus=16, mem_mb=1000000, time_min=7200, mem_gb = 500
    shell:
          """
          metaspades.py -t {resources.cpus} -m {resources.mem_gb} -1 {input.f_reads} -2 {input.r_reads} -o {output}
          """

##Alternative to SHI7 using previously QCd reads:
rule format_for_shogun:
    input:
        samples=expand("data/qc_sequence_files/{sample}_{read_dir}.fastq.gz", sample=samples, read_dir=read_dirs)
    output:
        "data/shogun/combined_samples.fasta"
    conda: "code/master_env.yaml"
    shell:
        """
        for i in {input};
        do
        s=$i
        #s={input.samples}
        s1=${{s##*/}}
        filename="${{s1%.fastq.gz}}"

        bioawk -c fastx -v fname="$filename" \'{{print ">"fname"_"++i" "$name"\\n"$seq}}\' $i >> {output}
        done
        #cat data/shogun/input_fastqs/*.fasta >> data/shogun/combined_seqs.fasta
        """


rule shogun_bt2: ###Need to get database setup
    input:
        seqs="data/shogun/combined_samples.fasta",
        db="/work/akiledal/WebOfLife/shogun/"
    output: directory("data/shogun/bt2")
    conda: "code/master_env.yaml"
    resources: cpus=60, mem_mb=100000, time_min=1440, mem_gb = 100
    shell:
        """
        shogun pipeline --no-function -i {input.seqs} -d {input.db} -o {output} -a bowtie2 --threads {resources.cpus}
        """

rule shogun_burst: ###Need to get database setup
    input:
        seqs="data/shogun/combined_samples.fasta",
        db="/work/akiledal/WebOfLife/shogun/"
    output: directory("data/shogun/burst")
    conda: "code/master_env.yaml"
    resources: cpus=60, mem_mb=800000, time_min=1440, mem_gb = 800
    shell:
        """
        shogun pipeline --no-function -i {input.seqs} -d {input.db} -o {output} -a burst --threads {resources.cpus}
        """

rule kraken2_gtdb: ##Run kraken2
    input:
        f_seq = "data/qc_sequence_files/{sample}_R1.fastq.gz",
        r_seq = "data/qc_sequence_files/{sample}_R2.fastq.gz"
    params:
        db = "/work/akiledal/kraken_gtdb"
    output:
        report = "data/kraken2/{sample}_report.txt",
        out = "data/kraken2/{sample}_out.txt",
        bracken = "data/kraken2/{sample}_bracken.txt",
        bracken_report = "data/kraken2/{sample}_brackenReport.txt",
        bracken_mpa = "data/kraken2/{sample}_brackenMpa.txt",
        unclass_f = "data/kraken2/{sample}_unclassified_1.fasta",
        unlcass_r = "data/kraken2/{sample}_unclassified_2.fasta"
    conda: "code/master_env.yaml"
    resources: cpus=8, mem_mb=500000, time_min=1440, mem_gb = 500
    shell:
        """
        kraken2 \
            --threads {resources.cpus} \
            --classified-out data/kraken2/{wildcards.sample}_classified#.fasta \
            --unclassified-out data/kraken2/{wildcards.sample}_unclassified#.fasta \
            --report {output.report} \
            --output {output.out} \
            --db {params.db} \
            --paired {input.f_seq} {input.r_seq}

        bracken -d {params.db} -i {output.report} -o {output.bracken} -w {output.bracken_report}

        ./code/kreport2mpa.py -r {output.bracken_report} -o {output.bracken_mpa} --percentages
        """

rule kraken2_refseq: ##Run kraken2
    input:
        f_seq = "data/kraken2/{sample}_unclassified_1.fasta",
        r_seq = "data/kraken2/{sample}_unclassified_2.fasta"
    params:
        db = "/work/akiledal/kraken"
    output:
        report = "data/kraken2_euk/{sample}_report.txt",
        out = "data/kraken2_euk/{sample}_out.txt",
        bracken = "data/kraken2_euk/{sample}_bracken.txt",
        bracken_report = "data/kraken2_euk/{sample}_brackenReport.txt"
    conda: "code/master_env.yaml"
    resources: cpus=8, mem_mb=500000, time_min=1440, mem_gb = 500
    shell:
        """
        kraken2 \
            --threads {resources.cpus} \
            --classified-out data/kraken2_euk/{wildcards.sample}_classified#.fasta \
            --unclassified-out data/kraken2_euk/{wildcards.sample}_unclassified#.fasta \
            --report {output.report} \
            --output {output.out} \
            --db {params.db} \
            --paired {input.f_seq} {input.r_seq}

        bracken -d {params.db} -i {output.report} -o {output.bracken} -w {output.bracken_report}
        """

rule humann_db_download:
    output:
        "/work/akiledal/humann3/all_genes_annot.1.bt2l",
        "/work/akiledal/humann3/all_genes_annot.2.bt2l",
        "/work/akiledal/humann3/all_genes_annot.3.bt2l",
        "/work/akiledal/humann3/all_genes_annot.4.bt2l",
        "/work/akiledal/humann3/all_genes_annot.rev.1.bt2l",
        "/work/akiledal/humann3/all_genes_annot.rev.2.bt2l",
        "/work/akiledal/humann3/genome_reps_filt_annot.faa.gz",
        "/work/akiledal/humann3/genome_reps_filt_annot.fna.gz",
        "/work/akiledal/humann3/genome_reps_filt_annot.tsv.gz",
        "/work/akiledal/humann3/protein_database/uniref90_201901.dmnd"
    resources: time_min=1440,
    shell:
        """
        cd /work/akiledal/humann3/
        wget -nv -r -nH --cut-dirs=6 -nc ftp://ftp.tue.mpg.de/ebio/projects/struo2/GTDB_release95/humann3/uniref90/
        """

rule humann:
    input:
        NUC_DB = "/work/akiledal/humann3/genome_reps_filt_annot.fna.gz",
        PROT_DB = "/work/akiledal/humann3/protein_database/uniref90_201901.dmnd",
        NUC_fol = "/work/akiledal/humann3/",
        PROT_fol = "/work/akiledal/humann3/protein_database/",
        f_seq = "data/qc_sequence_files/{sample}_R1.fastq.gz",
        r_seq = "data/qc_sequence_files/{sample}_R2.fastq.gz",
        bracken_mpa = "data/kraken2/{sample}_brackenMpa.txt"
    output:
        humann_output = directory("data/humann3/{sample}")
    params:
        mem_use = "maximum"
    conda: "code/master_env.yaml"
    resources: cpus=63, mem_mb=800000, time_min=1440, mem_gb = 800
    shell:
        """

        gunzip -c {input.f_seq} > data/qc_sequence_files/{wildcards.sample}_combined.fastq
        gunzip -c {input.r_seq} >> data/qc_sequence_files/{wildcards.sample}_combined.fastq

        #Copy databases to scratch
        #pre this: --nucleotide-database {input.NUC_fol} \
        #          --protein-database {input.PROT_fol} \
        #          --input data/qc_sequence_files/{wildcards.sample}_combined.fastq.gz \

        mkdir -p {config[scratch]} || exit $?
        cp -r {input.NUC_fol}*  {config[scratch]}db/  || exit $?

        cp data/qc_sequence_files/{wildcards.sample}_combined.fastq {config[scratch]}
        rm data/qc_sequence_files/{wildcards.sample}_combined.fastq

        cd {config[scratch]}

        humann3 --bypass-nucleotide-index \
            --threads {resources.cpus} \
            --memory-use {params.mem_use} \
            --taxonomic-profile {input.bracken_mpa} \
            --nucleotide-database ./db/ \
            --protein-database ./db/protein_database/ \
            --output-basename {wildcards.sample}_humann \
            --input {wildcards.sample}_combined.fastq \
            --output {wildcards.sample}/

        mkdir -p {output.humann_output}
        cp {wildcards.sample}/* {output.humann_output}/
        rm -rf {config[scratch]} || exit $?

        """


#Predict genes with prodigal
rule prodigal:
    input:
        contigs="data/assembly_metaspades/{sample}/contigs.fasta"
    output: directory("data/prodigal/{sample}")
    conda: "code/master_env.yaml"
    resources: cpus=1, mem_mb=50000, time_min=1440
    shell:
        """
        mkdir {output}
        prodigal -i {input.contigs} -o {output}/genes.gbk -a {output}/proteins.faa -p meta
        """

#FOAM hmm functional profiling
##redirected hmmsearch .hm output to null dev because it was not needed and large, could restore with: -o data/FOAM/{wildcards.sample}_{wildcards.dir}_sixframe.hm \
#doesn't scale well beyond ~6 cpus. Could be sped up in the future by splitting the fastqs into chunks, running hmmsearch on each in parallel and the re-merging

#splits = [x for x in range(20)] ##didn't work because pyfasta makes 01 instead of 1

splits = ["00", "01", "02" , "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19"]

rule fastq_split:
    input: "raw_data/{sample}_{dir}.fastq.gz"
    output:
        gunzip = "raw_data/{sample}_{dir}.fastq",
        sixframe = "data/FOAM/six_frame_{sample}-{dir}.fasta",
        split_six = expand("data/FOAM/split/six_frame_{sample}-{dir}.{num}.fasta", num = splits, allow_missing=True)
    params:
        nsplit = 20
    resources: cpus=1, mem_mb=100000, time_min=240
    shell:
        """
        #transeq doesn't seem to support zipped fastqs
        gunzip -k {input}

        echo "Files unzipped"

        #export qc reads in all six translation frames for hmmsearch comparison against FOAM database
        transeq {output.gunzip} -outseq {output.sixframe} -frame=6

        echo "Six frame translations made"

        mkdir -p data/FOAM/split

        pyfasta split -n {params.nsplit} {output.sixframe}

        echo "Files split"

        mv  data/FOAM/six_frame_{wildcards.sample}-{wildcards.dir}.*.fasta data/FOAM/split/
        """
    #
    # """
    # #transeq doesn't seem to support zipped fastqs
    # gunzip -k {input}
    #
    # #export qc reads in all six translation frames for hmmsearch comparison against FOAM database
    # transeq {output.gunzip} -outseq {output.sixframe} -frame=6
    #
    # #rm raw_data/{wildcards.sample}_{wildcards.dir}.fastq    #remove temp file
    #
    # mkdir -p data/FOAM/split
    # #split -d -l 40000 --additional-suffix .fasta data/FOAM/six_frame_{wildcards.sample}_{wildcards.dir}.fasta data/FOAM/split/six_frame_{wildcards.sample}-{wildcards.dir}_split
    #
    # pyfasta split -n {params.nsplit} {output.sixframe}
    #
    # mv data/FOAM/six_frame_*.*.fasta data/FOAM/split/
    # """



rule foam_raw:
    input:
        fasta= "data/FOAM/split/six_frame_{sample}-{dir}.{num}.fasta"
    output:
        counts="data/FOAM/{sample}-{dir}.{num}_sixframe_summary.txt"
    conda: "code/master_env.yaml"
    resources: cpus=4, mem_mb=25000, time_min=28800
    shell:
        """
        fp={input.fasta}
        fp2="${{fp%.fasta}}"
        num=${{fp2: -2}}

        #find hits against FOAM database using hmmsearch
        hmmsearch \
        	--cpu {resources.cpus} \
        	--tblout data/FOAM/{wildcards.sample}_{wildcards.dir}-${{num}}_sixframe.txt \
        	--domT 14 \
        	data/reference/FOAM-hmm_rel1_converted.hmm \
        	{input.fasta} >/dev/null

        #Make data more easily parseable and remove unneeded columns
        tail -n+4 data/FOAM/{wildcards.sample}_{wildcards.dir}-${{num}}_sixframe.txt | sed 's/ * / /g' | cut -f 1-10 -d " " > {output.counts}

        rm {input.fasta}
        #rm data/FOAM/{wildcards.sample}_{wildcards.dir}_${{num}}_sixframe.txt

        #cat data/FOAM/{wildcards.sample}_{wildcards.dir}_sixframe_summary.txt >> data/FOAM/{wildcards.sample}_{wildcards.dir}_sixframe_summary.txt

        #rm data/FOAM/{wildcards.sample}_{wildcards.dir}_${{num}}_sixframe_summary.txt
        """

rule foam_concat:
    input:
        sixframe_simp = expand("data/FOAM/{sample}-{dir}.{num}_sixframe_summary.txt", num = splits, allow_missing=True)
    output: "data/FOAM/{sample}-{dir}_sixframe_summary.txt"
    params:
        rm_file = expand("data/FOAM/{sample}-{dir}.{num}_sixframe.txt", num = splits, allow_missing=True)
    conda: "code/master_env.yaml"
    resources: cpus=1, mem_mb=16000, time_min=30
    shell:
        """
        #remove unsimplified split hmmer results
        for file in {params.rm_file}
        do
            rm -f $file
        done

        #concatenate split simplified hmmer/FOAM results
        for file in {input.sixframe_simp}
        do
            cat $file >> {output}
        done

        #delete split files
        for file in {input.sixframe_simp}
        do
            rm -f $file
        done
        """


# rule foam_raw:
#     input:
#         fastq="raw_data/{sample}_{dir}.fastq.gz"
#     output:
#         counts="data/FOAM/{sample}_{dir}_sixframe_summary.txt"
#     conda: "code/master_env.yaml"
#     resources: cpus=4, mem_mb=25000, time_min=28800
#     shell:
#         """
#         #transeq doesn't seem to support zipped fastqs
#         gunzip -k {input.fastq}
#
#         #export qc reads in all six translation frames for hmmsearch comparison against FOAM database
#         transeq raw_data/{wildcards.sample}_{wildcards.dir}.fastq -outseq data/FOAM/six_frame_{wildcards.sample}_{wildcards.dir}.fasta -frame=6
#
#         rm raw_data/{wildcards.sample}_{wildcards.dir}.fastq    #remove temp file
#
#         #find hits against FOAM database using hmmsearch
#         hmmsearch \
#         	--cpu {resources.cpus} \
#         	--tblout data/FOAM/{wildcards.sample}_{wildcards.dir}_sixframe.txt \
#         	--domT 14 \
#         	data/reference/FOAM-hmm_rel1_converted.hmm \
#         	data/FOAM/six_frame_{wildcards.sample}_{wildcards.dir}.fasta >/dev/null
#
#         #Make data more easily parseable and remove unneeded columns
#         tail -n+4 data/FOAM/{wildcards.sample}_{wildcards.dir}_sixframe.txt | sed 's/ * / /g' | cut -f 1-10 -d " " > {output.counts}
#
#         rm data/FOAM/six_frame_{wildcards.sample}_{wildcards.dir}.fasta
#
#         """

rule foam_process:
    input:
        foam_f = "data/FOAM/{sample}-R1_sixframe_summary.txt",
        foam_r = "data/FOAM/{sample}-R2_sixframe_summary.txt"
    output:
        combined="data/FOAM/{sample}_combined.rds"
    resources: cpus=6, mem_mb=200000, time_min=28800
    shell:
        """
        R --vanilla << 'RSCRIPT'
        library(data.table)
        library(vroom)
        library(tidyverse)

        names <- c(
        		"target.name",
        		"target.accession",
        		"query.name",
        		"query.accession",
        		"full.sequence.E.value",
        		"full.sequence.score",
        		"full.sequence.bias",
        		"best.1.domain.E.value",
        		"best.1.domain.score",
        		"best.1.domain.bias")

        foam_f <- vroom("{input.foam_f}",col_names = names,num_threads = 8) %>%
         filter(best.1.domain.score >= 14) %>%  #Min value used in the FOAM HMMerBestHit.py script
         mutate(target.name = str_remove(target.name, "_[0-9]$")) %>%
         group_by(target.name) %>%
         top_n(1,full.sequence.E.value) %>%
         separate_rows(query.name,sep = ",") %>%
         ungroup() %>% group_by(query.name) %>%
         summarise(meta_abund = n()) %>%
         rename(id = "query.name") %>%
         mutate(sample = "S3_fallen",
                id = str_remove(id,"KO:"))


        foam_r <- vroom("{input.foam_r}",col_names = names,num_threads = 8) %>%
         filter(best.1.domain.score >= 14) %>%  #Min value used in the FOAM HMMerBestHit.py script
         mutate(target.name = str_remove(target.name, "_[0-9]$")) %>%
         group_by(target.name) %>%
         top_n(1,full.sequence.E.value) %>%
         separate_rows(query.name,sep = ",") %>%
         ungroup() %>% group_by(query.name) %>%
         summarise(meta_abund = n()) %>%
         rename(id = "query.name") %>%
         mutate(sample = "S3_fallen",
                id = str_remove(id,"KO:"))


        foam_combined <- bind_rows(foam_f,foam_r) %>%
            group_by(sample, id) %>%
            mutate(meta_abund = sum(meta_abund))

        write_rds(foam_combined,"{output.combined}")

        'RSCRIPT'
        """


# Run BLASTS for ectoine synthesis and degradation

rule make_nuc_blastdb:
    input:
        db = "data/assembly_metaspades/S3_Fallen/contigs.fasta.gz",
    output:
        db_index = "data/assembly_metaspades/S3_Fallen/contigs.fasta.nin"
    log:
        "logs/BLAST/makeblastdb_contigs_nuc.log"
    resources: cpus=1, mem_mb=5000, time_min=120
    shell:
        """
        gunzip -c {input.db} > data/assembly_metaspades/S3_Fallen/contigs.fasta
        
        makeblastdb -in data/assembly_metaspades/S3_Fallen/contigs.fasta -dbtype nucl -logfile logs/BLAST/makeblastdb_contigs_nuc.dblog
        """


rule blast_nuc:
    input:
        blast_db = "data/assembly_metaspades/S3_Fallen/contigs.fasta",
        blast_db_index = rules.make_nuc_blastdb.output.db_index,
        genes = "data/reference/ectoine_metabolism_fasta.fasta"
    output:
        blast_res = "data/BLAST/ectoine_genes_BLAST.blastn"
    log:
        "logs/BLAST/ectoine_genes_nuc.log"
    resources: cpus=8, mem_mb=5000, time_min=120
    shell:
        """
        tblastn -query {input.genes} \
            -db {input.blast_db} \
            -out {output.blast_res} \
            -outfmt '6 std qcovs stitle qseq sseq' \
            -num_threads {resources.cpus}
        """


rule kraken2_gtdb_contigs: ##Run kraken2
    input:
        f_seq = "data/assembly_metaspades/S3_Fallen/contigs.fasta"
    params:
        db = "/work/akiledal/kraken2_db/gtdb_r202"
    output:
        report = "data/kraken2_contigs/gtdb/report.txt",
        out = "data/kraken2_contigs/gtdb/out.txt",
        bracken = "data/kraken2_contigs/gtdb/bracken.txt",
        bracken_report = "data/kraken2_contigs/gtdb/brackenReport.txt",
        bracken_mpa = "data/kraken2_contigs/gtdb/brackenMpa.txt",
        unclass_f = "data/kraken2_contigs/gtdb/unclassified.fasta"
    conda: "code/kraken.yaml"
    resources: cpus=8, mem_mb=500000, time_min=1440, mem_gb = 500
    shell:
        """
        kraken2 \
            --threads {resources.cpus} \
            --classified-out data/kraken2_contigs/gtdb/classified.fasta \
            --unclassified-out data/kraken2_contigs/gtdb/unclassified.fasta \
            --report {output.report} \
            --output {output.out} \
            --db {params.db} \
            {input.f_seq}

        bracken -d {params.db} -i {output.report} -o {output.bracken} -w {output.bracken_report}

        ./code/kreport2mpa.py -r {output.bracken_report} -o {output.bracken_mpa} --percentages
        """

rule kraken2_refseq_contigs: ##Run kraken2
    input:
        f_seq = "data/kraken2_contigs/gtdb/unclassified.fasta"
    params:
        db = "/work/akiledal/kraken2_db/refseq"
    output:
        report = "data/kraken2_contigs/refseq/report.txt",
        out = "data/kraken2_contigs/refseq/out.txt",
        bracken = "data/kraken2_contigs/refseq/bracken.txt",
        bracken_report = "data/kraken2_contigs/refseq/brackenReport.txt"
    conda: "code/kraken.yaml"
    resources: cpus=8, mem_mb=500000, time_min=1440, mem_gb = 500
    shell:
        """
        kraken2 \
            --threads {resources.cpus} \
            --classified-out data/kraken2_contigs/refseq/classified.fasta \
            --unclassified-out data/kraken2_contigs/refseq/unclassified.fasta \
            --report {output.report} \
            --output {output.out} \
            --db {params.db} \
            {input.f_seq}

        bracken -d {params.db} -i {output.report} -o {output.bracken} -w {output.bracken_report}
        """


rule translate_kraken: ##Run kraken2
    input:
        script = "code/translateKraken2.py",
        refseq = "data/kraken2_contigs/refseq/out.txt",
        gtdb = "data/kraken2_contigs/gtdb/out.txt"
    params:
        db_gtdb = "/work/akiledal/kraken2_db/gtdb_r202",
        db_refseq = "/work/akiledal/kraken2_db/gtdb_r202"
    output:
        refseq = "data/kraken2_contigs/refseq/translated_out.txt",
        gtdb = "data/kraken2_contigs/gtdb/translated_out.txt"
    conda: "code/kraken.yaml"
    resources: cpus=1, mem_mb=50000, time_min=1440, mem_gb = 50
    shell:
        """
        python {input.script} \
            --krakenout {input.refseq} \
            --translatedout {output.refseq} \
            --taxdatadir {params.db_refseq}/taxonomy \
            --taxaprefixes true
        
        python {input.script} \
            --krakenout {input.gtdb} \
            --translatedout {output.gtdb} \
            --taxdatadir {params.db_gtdb}/taxonomy \
            --taxaprefixes true
        """


rule contig_coverage: 
    input:
        contigs = "data/assembly_metaspades/S3_Fallen/contigs.fasta",
        fwd_reads = "data/qc_sequence_files/S3_Fallen_R1.fastq.gz",
        rev_reads = "data/qc_sequence_files/S3_Fallen_R2.fastq.gz"
    output:
        coverage = "data/contig_coverage.tsv"
    conda: "code/kraken.yaml"
    resources: cpus=16, mem_mb=50000, time_min=1440, mem_gb = 50
    shell:
        """
        coverm contig \
            -1 {input.fwd_reads} \
            -2 {input.rev_reads} \
            -r {input.contigs} \
            --min-read-percent-identity 95 \
            --methods mean rpkm tpm length covered_bases \
            -o {output.coverage} \
            -t {resources.cpus}
        """

rule uniref_lca_taxonomy:
    input:
        fwd_reads = "data/qc_sequence_files/S3_Fallen_R1.fastq.gz",
        rev_reads = "data/qc_sequence_files/S3_Fallen_R2.fastq.gz"
    output:
        report = "data/mmseqs/concreteMetaG_report"
    conda:  "code/mmseqs_env.yaml"
    params:
        unirefDB = "/home/akiledal/work_akiledal/mmseqs_unirefdb/mmseqs2/uniref100",
        out_prefix = "data/mmseqs/concreteMetaG",
        tmp_dir = "/home/akiledal/work_akiledal/mmseqs_unirefdb/tmp/mmseqs2"
    resources:
        mem_mb = 500000, cpus=64, time_min=20000
    shell:
        """
        mkdir -p data/mmseqs

        mmseqs \
            easy-taxonomy \
            {input.fwd_reads} {input.rev_reads} \
            {params.unirefDB} \
            ./{params.out_prefix} \
            {params.tmp_dir} \
            --lca-mode 3 \
            --tax-lineage 1 \
            --threads {resources.cpus} \
            --split-memory-limit 450G
        """

rule get_human_genome:
    output: "data/reference/human_genome.fasta.gz"
    resources: time_min=5000
    shell:
        """
        wget -O {output} http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.p13.genome.fa.gz
        """

rule human_genome_coverage:
    input:
        human_genome = rules.get_human_genome.output,
        fwd_reads = "data/qc_sequence_files/S3_Fallen_R1.fastq.gz",
        rev_reads = "data/qc_sequence_files/S3_Fallen_R2.fastq.gz"
    params: 
        bam_dir = "data/map_to_human_genome"
    output: "data/human_genome_coverage.tsv"
    conda: "coverm"
    resources: cpus=64, mem_mb=250000, time_min=10000
    shell:
        """
        # Currently coverm doesn't support gziped fasta genomes
        zcat {input.human_genome} > /tmp/human_genome.fasta

        coverm genome \
            -t {resources.cpus} \
            -m relative_abundance mean covered_bases variance length count rpkm tpm \
            --contig-end-exclusion 0 \
            --min-covered-fraction 0 \
            -1 {input.fwd_reads} \
            -2 {input.rev_reads} \
            --genome-fasta-files /tmp/human_genome.fasta \
            --bam-file-cache-directory {params.bam_dir} \
            -o {output}

        rm /tmp/human_genome.fasta
        """


### DELETE ALL RESULTS ###
### ONLY USE TO RE-BUILD ALL RESULTS ###

rule clean:
    shell:
        """
        rm -rf data/assembly_metaspades/*
        rm -rf raw_data/fastqc_reports
        rm -rf data/qc_sequence_files
        rm -rf data/shi7_out
        rm -rf results/sequence_quality
        rm -f data/shogun/combined_samples.fasta
        rm -rf data/prodigal
        rm -rf data/FOAM
        """
