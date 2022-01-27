version 1.0

import "tasks/task_consensus.wdl" as consensus

workflow DeNovo_build_bacterial_consensus {

  input {
    String    id
    File      read1
    File      read2
    File      reference_seq
  }

  call consensus.bwa_pe_de_novo{
    input:
      id=id,
      read1_trim=read1,
      read2_trim=read1,
      reference_seq=reference_seq
  }

  call consensus.sam_to_bam {
    input:
      id=id,
      samfile=bwa_pe_de_novo.samfile
  }

  call consensus.bcftools_consensus {
    input:
      id=id,
      reference_seq=reference_seq,
      sorted_bam=sam_to_bam.sorted_bam
  }

  call consensus.consensus_qc {
    input:
      id=id,
      reference_seq=reference_seq,
      assembly_fasta=bcftools_consensus.consensus_seq
  }


  output {
    File    bwa_pe_de_novo_software=bwa_pe_de_novo.image_software.txt
    File    de_novo_bam=sam_to_bam.bamfile
    File    de_novo_sorted_bam=sam_to_bam.sorted_bam
    File    de_novo_indexed_bam=sam_to_bam.indexed_bam
    File    de_novo_sam_to_bam_software=sam_to_bam.image_software
    File    de_novo_consensus_seq=bcftools_consensus.consensus_seq
    File    de_novo_consensus_variants=bcftools_consensus.consensus_variants
    File    de_novo_bcftools_consensus_software=bcftools_consensus.image_software
    File    de_novo_consensus_qc_software=consensus_qc.image_software

    Int  de_novo_consensus_number_N=consensus_qc.consensus_number_N
    Int  de_novo_consensus_number_ATCG=consensus_qc.consensus_number_ATCG
    Int  de_novo_consensus_number_Degenerate=consensus_qc.consensus_number_Degenerate
    Int  de_novo_consensus_number_Total=consensus_qc.consensus_number_Total
    Int  de_novo_consensus_percent_reference_coverage=consensus_qc.percent_reference_coverage




  }
}
