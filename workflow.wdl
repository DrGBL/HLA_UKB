version 1.0

#these are files that are required for my workflow. They should be on your dna nexus folder
#https://github.com/broadinstitute/gatk/raw/master/src/test/resources/large/Homo_sapiens_assembly38.fasta.gz
#https://github.com/broadinstitute/gatk/raw/master/src/test/resources/large/Homo_sapiens_assembly38.fasta.gz.fai
#https://github.com/broadinstitute/gatk/raw/master/src/test/resources/large/Homo_sapiens_assembly38.fasta.gz.gzi
#https://github.com/broadinstitute/gatk/raw/master/src/test/resources/large/Homo_sapiens_assembly38.dict



workflow hla_calling_wf {
    input {
        File cram_file
        File ref_genome
        File ref_genome_index
        File ref_genome_dict
        File ref_genome_gzi
        Int cores_hla
    }

    call fatsq_conversion { 
        input: cram_file = cram_file, ref_genome = ref_genome, ref_genome_index = ref_genome_index, ref_genome_dict = ref_genome_dict, ref_genome_gzi = ref_genome_gzi
    }

    call hla_calling {
        input: fastq1 = fatsq_conversion.fastq1, fastq2 = fatsq_conversion.fastq2, cores_hla = cores_hla
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
        Int cores_hla
    }

    command <<<
        hlahd.sh -t ~{cores_hla} -m 50 -f /HLA/hlahd.1.4.0/freq_data/ ~{fastq1} ~{fastq2}  /HLA/hlahd.1.4.0/HLA_gene.split.3.32.0.txt /HLA/hlahd.1.4.0/dictionary/ test_call .
        tar -zcf  hla_res.tar.gz test_call/result/
    >>>

    runtime {
        #if you want to choose the specific instances for the calling.
        #dx_instance_type: "mem3_ssd1_v2_x8"
        #dx_instance_type: "mem3_ssd1_v2_x32"
    }

    output {
        File out_hla = "hla_res.tar.gz"
    }

}
