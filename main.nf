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
            singularity pull bwa-0.7.17.sif docker://index.docker.io/mpgagebioinformatics/bwa:0.7.17
        fi
        
    fi
    if [[ "${params.run_type}" == "local" ]] ; 
      then
        docker pull mpgagebioinformatics/samtools:1.16.1
        docker pull mpgagebioinformatics/bwa:0.7.17
    fi
    """

}

process genome_collector {
  stageInMode 'symlink'
  stageOutMode 'move'

  output:
    val "finished", emit: get_genome_status

  when:
    ( ! file("${params.genomes}/${params.organism}/${params.release}/${params.organism}.${params.release}.genome").exists() )
  
  script:
    """
    target_folder=/genomes/${params.organism}/${params.release}/

    if [[ ! -e \$target_folder ]] ; then mkdir -p \$target_folder ; fi

    cd \$target_folder

    if [[ ! -e ${params.organism}.${params.release}.gtf ]] ; 
      then
        curl -#O ${params.url_gtf} && gtf=`basename ${params.url_gtf}` || gtf=`curl -#l ${params.url_gtf} | grep "gtf" | grep -v abinitio` && curl -#O ${params.url_gtf}/\$gtf
        if [[ "\$gtf" == *".gz" ]] ; then unpigz -p ${task.cpus} \$gtf ; gtf=\${gtf%.gz} ; fi
        mv \$gtf ${params.organism}.${params.release}.gtf
        grep -v -i 'biotype "rRNA' ${params.organism}.${params.release}.gtf | grep -v -i "Mt_rRNA" | grep -v -i srrna > ${params.organism}.${params.release}.no.rRNA.gtf
    fi

    if [[ ! -e ${params.organism}.${params.release}.fa ]] ; 
      then
        curl -#O ${params.url_dna} && dna=\$(basename ${params.url_dna} ) || dna=""
        if [[ ! -f \$dna ]] ;
          then 
            dna=\$(curl -#l ${params.url_dna} | grep .dna.toplevel.fa.gz)
            curl -#O ${params.url_dna}/\$dna
        fi
        if [[ "\$dna" == *".gz" ]] ; then unpigz \$dna ; dna=\${dna%.gz} ; fi
        mv \$dna ${params.organism}.${params.release}.fa
    fi

    if [[ ! -e ${params.organism}.${params.release}.genome ]] ;
      then
        samtools faidx ${params.organism}.${params.release}.fa
        awk '{print \$1"\t"\$2}' ${params.organism}.${params.release}.fa.fai > ${params.organism}.${params.release}.genome
    fi

    """
}

process indexer {
  stageInMode 'symlink'
  stageOutMode 'move'

  when:
    ( ! file("${params.genomes}/${params.organism}/${params.release}/toplevel_bwa/index.fa").exists() ) 
  
  script:
    """

    target_folder=/genomes/${params.organism}/${params.release}/toplevel_bwa

    if [[ ! -e \$target_folder ]] ; then mkdir -p \$target_folder ; fi

    cd \$target_folder

    ln -s ${params.genomes}/${params.organism}/${params.release}/${params.organism}.${params.release}.fa index.fa

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
    samtools sort -@ 10 -o ${pair_id}.sorted.bam ${pair_id}.bam
    samtools index ${pair_id}.sorted.bam

    """
}

workflow images {
  main:
    get_images()
}

workflow get_genome {
  main:
    genome_collector()
}

workflow index {
  main:
    indexer()
}

workflow map_reads {
  main:
    // Channel
    //   .fromFilePairs( "${params.kallisto_raw_data}*.READ_{1,2}.fastq.gz", size: -1 )
    //   .ifEmpty { error "Cannot find any reads matching: ${params.kallisto_raw_data}*.READ_{1,2}.fastq.gz" }
    //   .set { read_files } 
    read_files=Channel.fromFilePairs( "${params.bwa_raw_data}/*.READ_{1,2}.fastq.gz", size: -1 ) 
    mapping( read_files )
    flagstat( mapping.out.collect(), read_files )
}





