process readTrimming {
    /**
    * Trims paired fastq using trim_galore (https://github.com/FelixKrueger/TrimGalore)
    * @input tuple(sampleName, path(forward), path(reverse))
    * @output trimgalore_out tuple(sampleName, path("*_val_1.fq.gz"), path("*_val_2.fq.gz"))
    */

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: '*_val_{1,2}.fq.gz', mode: 'copy'

    cpus 2

    input:
    tuple(sampleName, path(forward), path(reverse))

    output:
    tuple(sampleName, path("*_val_1.fq.gz"), path("*_val_2.fq.gz")) optional true

    script:
    """
    if [[ ! -s ${forward} ]]; then
      #exit 0
      touch ${sampleName}_val_1.fq.gz ${sampleName}_val_2.fq.gz
    else
      trim_galore --paired $forward $reverse
    fi
    """
}

process indexReference {
    /**
    * Indexes reference fasta file in the scheme repo using bwa.
    */

    tag { ref }

    input:
        path(ref)

    output:
        tuple path('ref.fa'), path('ref.fa.*')

    script:
        """
        ln -s ${ref} ref.fa
        bwa index ref.fa
        """
}

process readMapping {
    /**
    * Maps trimmed paired fastq using BWA (http://bio-bwa.sourceforge.net/)
    * Uses samtools to convert to BAM, sort and index sorted BAM (http://www.htslib.org/doc/samtools.html)
    * @input 
    * @output 
    */

    tag { sampleName }

    label 'largecpu'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.sorted.bam", mode: 'copy'
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.sorted.bam.bai", mode: 'copy'

    input:
        tuple sampleName, path(forward), path(reverse), path(ref), path("*"), path(posc)

    output:
        tuple sampleName, path("${sampleName}.sorted.bam"), emit: bam
        path "${sampleName}.sorted.bam.bai"
        path "posc/${sampleName}.tsv", emit: posc

    script:
      """
      if [[ ! -s ${forward} ]]; then
        mkdir posc
        touch ${sampleName}.sorted.bam ${sampleName}.sorted.bam.bai posc/${sampleName}.tsv
      else
        bwa mem -t ${task.cpus} ${ref} ${forward} ${reverse} | \
        samtools sort -o ${sampleName}.sorted.bam
        samtools index ${sampleName}.sorted.bam

        mkdir posc
        bwa index ${posc}
        bwa mem -a -t ${task.cpus} ${posc} ${forward} ${reverse} | \
        samtools sort -o ${sampleName}.posc.bam
        samtools index ${sampleName}.posc.bam
        bedtools genomecov -ibam "${sampleName}.posc.bam" -d > "${sampleName}.depth.tsv"
        for ref in \$(cut -f 1 "${sampleName}.depth.tsv" | sort -u); do cov=\$(bc <<< "scale=3; (\$(awk -v ref=\$ref 'BEGIN {bp=0} \$1==ref && \$3>=10 {bp+=\$2} {print bp}' "${sampleName}.depth.tsv" | tail -1)/\$(awk -v ref=\$ref 'BEGIN {bp=0} \$1==ref {bp+=\$2} {print bp}' "${sampleName}.depth.tsv" | tail -1))"); echo -e "\${ref}\\t\${cov}" >> "posc/${sampleName}.tsv"; done;
      fi
      """
}

process trimPrimerSequences {

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.mapped.bam", mode: 'copy'
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.mapped.bam.bai", mode: 'copy'
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.mapped.primertrimmed.sorted.bam", mode: 'copy'
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.mapped.primertrimmed.sorted.bam.bai", mode: 'copy'

    input:
    tuple sampleName, path(bam), path(bedfile)

    output:
    path "${sampleName}.mapped.bam.bai"
    path "${sampleName}.mapped.primertrimmed.sorted.bam.bai"
    tuple sampleName, path("${sampleName}.mapped.bam"), emit: mapped
    tuple sampleName, path("${sampleName}.mapped.primertrimmed.sorted.bam" ), emit: ptrim

    script:
    if (params.allowNoprimer){
        ivarCmd = "ivar trim -e"
    } else {
        ivarCmd = "ivar trim"
    }
   
    if ( params.cleanBamHeader )
        """
        if [[ ! -s ${bam} ]]; then
          touch ${sampleName}.mapped.bam ${sampleName}.mapped.bam.bai ${sampleName}.mapped.primertrimmed.sorted.bam ${sampleName}.mapped.primertrimmed.sorted.bam.bai
        else
          samtools reheader --no-PG  -c 'sed "s/${sampleName}/sample/g"' ${bam} | \
          samtools view -F4 -o sample.mapped.bam

          mv sample.mapped.bam ${sampleName}.mapped.bam
        
          samtools index ${sampleName}.mapped.bam

          ${ivarCmd} -i ${sampleName}.mapped.bam -b ${bedfile} -m ${params.illuminaKeepLen} -q ${params.illuminaQualThreshold} -p ivar.out

          samtools reheader --no-PG  -c 'sed "s/${sampleName}/sample/g"' ivar.out.bam | \
          samtools sort -o sample.mapped.primertrimmed.sorted.bam

          mv sample.mapped.primertrimmed.sorted.bam ${sampleName}.mapped.primertrimmed.sorted.bam
          samtools index ${sampleName}.mapped.primertrimmed.sorted.bam
        fi
        """

    else
        """
        if [[ ! -s ${bam} ]]; then
          touch ${sampleName}.mapped.bam ${sampleName}.mapped.bam.bai ${sampleName}.mapped.primertrimmed.sorted.bam ${sampleName}.mapped.primertrimmed.sorted.bam.bai
        else
          samtools view -F4 -o ${sampleName}.mapped.bam ${bam}
          samtools index ${sampleName}.mapped.bam

          ${ivarCmd} -i ${sampleName}.mapped.bam -b ${bedfile} -m ${params.illuminaKeepLen} -q ${params.illuminaQualThreshold} -p ivar.out
          samtools sort -o ${sampleName}.mapped.primertrimmed.sorted.bam ivar.out.bam
          samtools index ${sampleName}.mapped.primertrimmed.sorted.bam
        fi
        """
}

process callVariants {

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.variants.tsv", mode: 'copy'

    input:
    tuple(sampleName, path(bam), path(ref))

    output:
    tuple sampleName, path("${sampleName}.variants.tsv"), emit: variants
    path "${sampleName}.count.tsv", emit: count
    path "${sampleName}.variants.tsv", emit: contaminant
    tuple sampleName, path("${sampleName}.variants.tsv"), path("${ref}"), path("${sampleName}.bed"), emit: contaminated

    script:
        """
        if [[ ! -s ${bam} ]]; then
          touch ${sampleName}.variants.tsv ${sampleName}.count.tsv ${sampleName}.bed
        else
          samtools mpileup -A -d 0 --reference ${ref} -B -Q 0 ${bam} |\
          ivar variants -r ${ref} -m ${params.ivarMinDepth} -p ${sampleName}.variants -q ${params.ivarMinVariantQuality} -t ${params.ivarMinFreqThreshold}
          refName=\$(grep ">" ${ref} | tr -d ">")
          echo -e "\${refName}\\t10-20%\\t\$(awk 'BEGIN {bp=0} \$11>=0.1 && \$11<0.2 {bp+=1} {print bp}' "${sampleName}.variants.tsv" | tail -1)\\t${sampleName}" > "${sampleName}.count.tsv"
          echo -e "\${refName}\\t20-50%\\t\$(awk 'BEGIN {bp=0} \$11>=0.2 && \$11<0.5 {bp+=1} {print bp}' "${sampleName}.variants.tsv" | tail -1)\\t${sampleName}" >> "${sampleName}.count.tsv"
          bedtools genomecov -ibam "${bam}" -bga > "${sampleName}.bed"
        fi
        """
}

process makeConsensus {

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.primertrimmed.consensus.fa", mode: 'copy'

    input:
        tuple(sampleName, path(bam), path(ref))

    output:
        tuple(sampleName, path("${sampleName}.primertrimmed.consensus.fa"))

    script:
        """
        if [[ ! -s ${bam} ]]; then
          echo ">${sampleName}" > ${sampleName}.primertrimmed.consensus.fa
          seq \$(grep -v ">" ${ref} | tr -d '\n' | wc -m) | sed "c N" | tr -d '\n' | sed -e '\$a\\' >> ${sampleName}.primertrimmed.consensus.fa
        else
          samtools mpileup -aa -A -B -d ${params.mpileupDepth} -Q0 ${bam} | \
          ivar consensus -t ${params.ivarFreqThreshold} -m ${params.ivarMinDepth} \
          -n N -p ${sampleName}.primertrimmed.consensus
        
          sed -i "s/>.*/>${sampleName}/" ${sampleName}.primertrimmed.consensus.fa
        fi
        """
}

process cramToFastq {
    /**
    * Converts CRAM to fastq (http://bio-bwa.sourceforge.net/)
    * Uses samtools to convert to CRAM, to FastQ (http://www.htslib.org/doc/samtools.html)
    * @input
    * @output
    */

    input:
        tuple sampleName, file(cram)

    output:
        tuple sampleName, path("${sampleName}_1.fastq.gz"), path("${sampleName}_2.fastq.gz")

    script:
        """
        samtools collate -u ${cram} -o tmp.bam
        samtools fastq -1 ${sampleName}_1.fastq.gz -2 ${sampleName}_2.fastq.gz tmp.bam
        rm tmp.bam
        """
}

process reportAllConsensus {
    /**
    * Concatenate consensus of all samples
    */

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "all_consensus.fasta", mode: 'copy'
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "all_trimcons.fasta", mode: 'copy'

    input:
        path("*.fa")

    output:
        path "all_consensus.fasta", emit: cons
        path "all_trimcons.fasta", emit: trimcons

    script:
      """
        cat *.fa > all_consensus.fasta
        seqtk trimfq -b 150 -e 150 all_consensus.fasta > all_trimcons.fasta
        awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%50==0){file=sprintf("consensus%d.fa",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < all_consensus.fasta
      """
}

process reportSampleCoverage {
    /**
    * Report coverage for each sample
    */

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.cov", mode: 'copy'

    input:
        tuple(sampleName, path(bam), path(ref))

    output:
        path("${sampleName}.cov")

    script:
      """
        if [[ ! -s ${bam} ]]; then
          refName=\$(grep ">" ${ref})
          refName="\${refName##>}"
          length=\$(grep -v ">" ${ref} | tr -d '\n' | wc -m)
          for i in \$(seq 1 \${length}); do echo -e "\${refName}\t\${i}\t0\t${sampleName}" >> ${sampleName}.cov; done;
        else
          samtools index ${bam}
          bedtools genomecov -ibam ${bam} -d -split | sed "s/\$/\t${sampleName}/" > ${sampleName}.cov
        fi
      """
}

process reportAllCoverage {
    /**
    * Report coverage summary of the run
    */

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "all_cov.cov", mode: 'copy'

    input:
        path("*.cov")

    output:
        path("all_cov.cov")

    script:
      """
        cat *.cov > all_cov.cov
      """
}

process reportAllCounts {
    /**
    * Report counts summary of the run
    */

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "all_counts.tsv", mode: 'copy'

    input:
        path("*.tsv")

    output:
        path("all_counts.tsv")

    script:
      """
        cat *.tsv > all_counts.tsv
      """
}


process seekContaminant {

    tag { sampleName }

    when:
      var.size() > 0

    input:
      tuple(sampleName, path(var), path(ref), path(bed), path(pool))
      file "vcf/*"
    
    output:
      tuple(path("raw_${poolname}_${sampleName}.tsv"), path("cov_${poolname}_${sampleName}.tsv"), path("common_${poolname}_${sampleName}.tsv"), path("expected_${poolname}_${sampleName}.tsv"))
    
    script:
    poolname = (pool =~ /([0-9a-zA-Z_\-]+)(.+)/)[0][1]
      """
      for file in vcf/*.tsv; do
      bn=\$(basename "\${file}")
      if [[ \${bn} != ${var} && -s \${file} ]]; then
          python2.7 ${params.scripts}/compare_vcf.py -r ${ref} -c \${file} -v ${var} -b ${bed} -d 100 -f 0.05 -m raw,cov,common,expected -R ${poolname}.bed
      fi
      done
      """
}

process mergeContaminant {

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "contamination_${mode}_${poolname}.tsv", mode: 'copy'

    input:
      file("*")
      tuple(path(pool), mode)
    
    output:
      file "contamination_${mode}_${poolname}.tsv"
    
    script:
    poolname = (pool =~ /([0-9a-zA-Z_\-]+)(.+)/)[0][1]
      """
      mkdir -p ${mode}/${poolname}/
      for file in ${mode}_${poolname}_*.tsv; do
        mv \${file} "${mode}/${poolname}/\${file##${mode}_${poolname}_}"
      done
      python2.7 ${params.scripts}/merge_tables.py ${mode}/${poolname}/*.tsv > contamination_${mode}_${poolname}.tsv
      """
}

process mergePosControls {

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "posc.tsv", mode: 'copy'

    input:
      file 'posc/*'
    
    output:
      file "posc.tsv"
    
    script:
      """
      python2.7 ${params.scripts}/merge_tables.py posc/*.tsv > posc.tsv
      """
}

process makeSummary {
    /**
    * Make final summary of the run
    */

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "summary.csv", mode: 'copy'

    input:
        tuple(path(consensus), path(trimcons), path(coverage), path(counts))

    output:
        path("summary.csv")

    script:
      """
        Rscript ${params.scripts}/AI_analysis.R \$PWD/
      """
}
