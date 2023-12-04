#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process get_images {
  stageInMode 'symlink'
  stageOutMode 'move'

  script:
    """

    if [[ "${params.containers}" == "singularity" ]] ; 
    
      then

        cd ${params.image_folder}

        if [[ ! -f samtools-1.16.1.sif ]] ;
          then
            singularity pull samtools-1.16.1.sif docker://index.docker.io/mpgagebioinformatics/samtools:1.16.1
        fi

        if [[ ! -f bwa-0.7.17.sif ]] ;
          then
            singularity pull bwa-0.7.17.sif docker://index.docker.io/mpgagebioinformatics/bwa:0.7.17
        fi
        
    fi


    if [[ "${params.containers}" == "docker" ]] ;

      then

        docker pull mpgagebioinformatics/samtools:1.16.1
        docker pull mpgagebioinformatics/bwa:0.7.17

    fi

    """

}


process bwa_indexer {
  stageInMode 'symlink'
  stageOutMode 'move'

  when:
    ( ! file("${params.genomes}/${params.organism}/${params.release}/toplevel_bwa/index.fa").exists() ) 
  
  script:
    """

    target_folder=/genomes/${params.organism}/${params.release}/toplevel_bwa

    if [[ ! -e \$target_folder ]] ; then mkdir -p \$target_folder ; fi

    cd \$target_folder

    ##ln -s ${params.genomes}/${params.organism}/${params.release}/${params.organism}.${params.release}.fa index.fa

    ln -s ../${params.organism}.${params.release}.fa index.fa
    
    bwa index -a bwtsw -p index.fa index.fa
    """
}

process mapping {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    tuple val(pair_id), path(fastq)

  output:
    val pair_id

  when:
    ( ! file("${params.project_folder}/bwa_output/${pair_id}.sam").exists() ) 
  
  script:
  def single = fastq instanceof Path

  if ( single ) {
    """
      mkdir -p /workdir/bwa_output
      cd /raw_data
      
      echo ${pair_id}
      bwa mem -t ${task.cpus} -M "${params.genomes}/${params.organism}/${params.release}/toplevel_bwa/index.fa" ${pair_id}.READ_1.fastq.gz > /workdir/bwa_output/${pair_id}.sam
  
    """
  } 
  else { 
    """
      mkdir -p /workdir/bwa_output
      cd /raw_data
      
      echo ${pair_id}
      bwa mem -t ${task.cpus} -M "${params.genomes}/${params.organism}/${params.release}/toplevel_bwa/index.fa" ${pair_id}.READ_1.fastq.gz ${pair_id}.READ_2.fastq.gz > /workdir/bwa_output/${pair_id}.sam 
    """
  }

}

process flagstat {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val pair_id
    tuple val(pair_id), path(fastq)

  when:
    ( ! file("${params.project_folder}/bwa_output/${pair_id}.sorted.bam.bai").exists() ) 


  script:
    """
    cd /workdir/bwa_output/

    samtools view -bS ${pair_id}.sam > ${pair_id}.bam
    samtools flagstat ${pair_id}.bam > ${pair_id}.bam.stat 
    samtools sort -@ ${task.cpus} -o ${pair_id}.sorted.bam ${pair_id}.bam
    samtools index ${pair_id}.sorted.bam

    """
}

workflow images {
  main:
    get_images()
}

workflow index {
  main:
    bwa_indexer()
}

workflow map_reads {
  main:
    read_files=Channel.fromFilePairs( "${params.bwa_raw_data}/*.READ_{1,2}.fastq.gz", size: -1 )
    mapping( read_files )
    flagstat( mapping.out.collect(), read_files )
}





