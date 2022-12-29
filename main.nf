#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process get_images {
  stageInMode 'symlink'
  stageOutMode 'move'

  script:
    """
    if [[ "${params.run_type}" == "r2d2" ]] || [[ "${params.run_type}" == "raven" ]] ; 
      then
        cd ${params.image_folder}
        if [[ ! -f samtools-1.16.1.sif ]] ;
          then
            singularity pull samtools-1.16.1.sif docker://index.docker.io/mpgagebioinformatics/samtools:1.16.1
        fi
        if [[ ! -f bwa-0.7.17.sif ]] ;
          then
            singularity pull samtools-1.16.1.sif docker://index.docker.io/mpgagebioinformatics/bwa:0.7.17
        fi
        
    fi
    if [[ "${params.run_type}" == "local" ]] ; 
      then
        docker pull mpgagebioinformatics/samtools:1.16.1
        docker pull mpgagebioinformatics/bwa:0.7.17
    fi
    """

}

workflow images {
  main:
    get_images()
}
