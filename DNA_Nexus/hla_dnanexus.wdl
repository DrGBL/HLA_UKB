version 1.0

#need to change docker image in the hla_calling task runtime section

workflow hla_calling_wf {
    input {
        File cram_file
        File ref_genome         #https://github.com/broadinstitute/gatk/raw/master/src/test/resources/large/Homo_sapiens_assembly38.fasta.gz
        File ref_genome_index   #https://github.com/broadinstitute/gatk/raw/master/src/test/resources/large/Homo_sapiens_assembly38.fasta.gz.fai
        File ref_genome_dict    #https://github.com/broadinstitute/gatk/raw/master/src/test/resources/large/Homo_sapiens_assembly38.dict
        File ref_genome_gzi     #https://github.com/broadinstitute/gatk/raw/master/src/test/resources/large/Homo_sapiens_assembly38.fasta.gz.gzi
        Int min_read_length
    }

    call fatsq_conversion { 
        input: cram_file = cram_file, ref_genome = ref_genome, ref_genome_index = ref_genome_index, ref_genome_dict = ref_genome_dict, ref_genome_gzi = ref_genome_gzi
    }

    call hla_calling {
        input: fastq1 = fatsq_conversion.fastq1, fastq2 = fatsq_conversion.fastq2, min_read_length = min_read_length
    }

    output {
        File hla_out = hla_calling.out_hla
    }
  
}


task fatsq_conversion {
    
    input {
        File cram_file
        File ref_genome
        File ref_genome_index
        File ref_genome_dict
        File ref_genome_gzi
    }

    command <<<
        gatk SamToFastq -I ~{cram_file} -F sample.1.fastq.gz -F2 sample.2.fastq.gz -R ~{ref_genome}
    >>>

    runtime {
        docker: "broadinstitute/gatk"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }

    output {
        File fastq1 = "sample.1.fastq.gz"
        File fastq2 = "sample.2.fastq.gz"
    }

}

task hla_calling {
    input {
        File fastq1
        File fastq2
        Int min_read_length

    }

    command <<<
        hlahd.sh -t 8 -m ~{min_read_length} -f /HLA/hlahd.1.4.0/freq_data/ ~{fastq1} ~{fastq2}  /HLA/hlahd.1.4.0/HLA_gene.split.3.32.0.txt /HLA/hlahd.1.4.0/dictionary/ test_call .
        tar -zcf  hla_res.tar.gz test_call/result/
    >>>

    runtime {
        docker: "" #need to add the required HLA-HD docker image
        dx_instance_type: "mem3_ssd1_v2_x8" #adjust if want more memore and RAM
    }

    output {
        File out_hla = "hla_res.tar.gz"
    }

}
