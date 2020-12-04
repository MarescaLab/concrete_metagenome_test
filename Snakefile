import os
import re

#Location of raw files
raw_seq_folder = "raw_data"
qc_seq_folder = "data/qc_sequence_files"

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

localrules: all, clean

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
        "data/shogun/bt2/",
        "data/shogun/burst/"

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
        reports = expand("{qc_folder}/fastqc_reports/{sample}_{read_dir}_fastqc.html", qc_folder=qc_seq_folder, sample=samples, read_dir=read_dirs)
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
        directory("data/assembly_megahit/{sample}")
    conda: "code/master_env.yaml"
    resources: cpus=8, mem_mb=100000, time_min=1440, mem_bytes = 100000000000
    shell:
        """
        megahit -m {resources.mem_bytes} -1 {input.f_reads} -2 {input.r_reads} -o {output}
        """

rule metaspades_assemble:
    input:
          f_reads = "data/qc_sequence_files/{sample}_R1.fastq.gz",
          r_reads = "data/qc_sequence_files/{sample}_R2.fastq.gz"
    output:
        directory("data/assembly_metaspades/{sample}")
    conda: "code/master_env.yaml"
    resources: cpus=8, mem_mb=100000, time_min=2160, mem_gb = 100
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

# rule shi7:
#     input: "raw_data/{sample}_{direction}.fastq.gz"
#     output:
#         "data/shi7_out/combined_seqs.fna"
#     conda:
#         "code/master_env.yaml"
#     shell:
#         """
#         mkdir data/shi7_out/input_files/
#         cp {input} data/shi7_out/input_files
#         gunzip data/shi7_out/input_files/*.fastq.gz
#
#         shi7 -i data/shi7_out/input_files/ -o data/shi7_out --adaptor Nextera --flash False --filter_length 50 -m 0
#
#         rm -rf data/shi7_out/input_files/
#         """

rule shogun_bt2: ###Need to get database setup
    input:
        seqs="data/shogun/combined_samples.fasta",
        db="/work/akiledal/shogun/"
    output: directory("data/shogun/bt2")
    conda: "code/master_env.yaml"
    resources: cpus=12, mem_mb=100000, time_min=1440, mem_gb = 100
    shell:
        """
        shogun pipeline -i {input.seqs} -d {input.db} -o {output} -a bowtie2 --threads {resources.cpus}
        """
rule shogun_burst: ###Need to get database setup
    input:
        seqs="data/shogun/combined_samples.fasta",
        db="/work/akiledal/shogun/"
    output: directory("data/shogun/burst")
    conda: "code/master_env.yaml"
    resources: cpus=48, mem_mb=800000, time_min=1440, mem_gb = 800
    shell:
        """
        shogun pipeline -i {input.seqs} -d {input.db} -o {output} -a burst --threads {resources.cpus}
        """


rule clean:
    shell:
        """
        rm -rf data/assembly_metaspades/*
        rm -rf raw_data/fastqc_reports
        rm -rf data/qc_sequence_files
        rm -rf data/shi7_out
        rm -rf results/sequence_quality
        rm -f data/shogun/combined_samples.fasta
        """
