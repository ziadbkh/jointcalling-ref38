#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.splitJaabs = false
params.batchSize = 250

batch_size = params.batchSize
sample_name_map = file(params.sample_name_map)
callset_name = params.callset_name
num_gvcfs = sample_name_map.readLines().size()

gatk_path = "/g/data/np30/users/nn6960/gatk/gatk-4.1.8.0"

unpadded_intervals_file = file("/g/data/np30/users/nn6960/ref38/hg38/v0/hg38.even.handcurated.20k.intervals")
num_of_original_intervals = unpadded_intervals_file.readLines().size()

process DynamicallyCombineIntervals {
    module 'singularity'

    input:
    path intervals
    val merge_count

    output:
    path "out.intervals", emit: output_intervals

    script:
    """
#!/usr/bin/env python
def parse_interval(interval):
    colon_split = interval.split(":")
    chromosome = colon_split[0]
    dash_split = colon_split[1].split("-")
    start = int(dash_split[0])
    end = int(dash_split[1])
    return chromosome, start, end

def add_interval(chr, start, end):
    lines_to_write.append(chr + ":" + str(start) + "-" + str(end))
    return chr, start, end

count = 0
chain_count = ${merge_count}
l_chr, l_start, l_end = "", 0, 0
lines_to_write = []
with open("${intervals}") as f:
    with open("out.intervals", "w") as f1:
        for line in f.readlines():
            # initialization
            if count == 0:
                w_chr, w_start, w_end = parse_interval(line)
                count = 1
                continue
            # reached number to combine, so spit out and start over
            if count == chain_count:
                l_char, l_start, l_end = add_interval(w_chr, w_start, w_end)
                w_chr, w_start, w_end = parse_interval(line)
                count = 1
                continue

            c_chr, c_start, c_end = parse_interval(line)
            # if adjacent keep the chain going
            if c_chr == w_chr and c_start == w_end + 1:
                w_end = c_end
                count += 1
                continue
            # not adjacent, end here and start a new chain
            else:
                l_char, l_start, l_end = add_interval(w_chr, w_start, w_end)
                w_chr, w_start, w_end = parse_interval(line)
                count = 1
        if l_char != w_chr or l_start != w_start or l_end != w_end:
            add_interval(w_chr, w_start, w_end)
        f1.writelines("\\n".join(lines_to_write))
    """
}

process getIntervalLength {
    module 'singularity'

    input:
    val size

    output:
    path 'range'

    script:
    """
    string="${size}"
    length=\${#string}
    for i in \$(seq 0 \$(($size-1))); do
        ilen=\${#i}
        for z in \$(seq \$ilen \$((\$length-1))); do
            echo -n "0"
            done
        echo \$i
    done > range
    """
}

process CombineGenotypeSitesOnlyGVCF {
    module 'singularity'

    input:
    path sample_name_map
    val unpadded_intervals
    each idx

    output:
    tuple val(idx), val(interval), path("${output_vcf_filename}"), path("${output_vcf_filename}.tbi"), emit: GenotypeGVCFOutput
    tuple val(idx), val(interval), path("${variant_filtered_vcf_filename}"), path("${variant_filtered_vcf_filename}.tbi"), emit: VariantFilterationOutput
    tuple path("${sites_only_vcf_filename}"), path("${sites_only_vcf_filename}.tbi"), emit: MakeSitesOnlyOutput

    script:
    interval = unpadded_intervals.get(idx.toInteger())
    output_vcf_filename = "${idx}.output.vcf.gz"
    excess_het_threshold = 54.69
    variant_filtered_vcf_filename = "${callset_name}.${idx}.variant_filtered.vcf.gz"
    sites_only_vcf_filename = "${callset_name}.${idx}.variant_filtered.sites_only.vcf.gz"
    workspace_dir_name = "genomicsdb.${idx}"
    """
    export TILEDB_DISABLE_FILE_LOCKING=1
    ${gatk_path}/gatk --java-options "-Xmx80g -Xms80g" \
        GenomicsDBImport \
        --genomicsdb-workspace-path ${workspace_dir_name} \
        --batch-size ${batch_size} \
        -L ${interval} \
        --sample-name-map ${sample_name_map} \
        --reader-threads 5 \
        --genomicsdb-shared-posixfs-optimizations \
        -ip 500

    ${gatk_path}/gatk --java-options "-Xmx80g -Xms80g" \
        GenotypeGVCFs \
        -R \${ref_fasta} \
        -O ${output_vcf_filename} \
        -D \${dbsnp_vcf} \
        -G StandardAnnotation \
        --only-output-calls-starting-in-intervals \
        -V gendb://${workspace_dir_name} \
        -L ${interval}

    ${gatk_path}/gatk --java-options "-Xmx80g -Xms80g" \
        VariantFiltration \
        --filter-expression "ExcessHet > ${excess_het_threshold}" \
        --filter-name ExcessHet \
        -O ${variant_filtered_vcf_filename} \
        -V ${output_vcf_filename}

    java -Xmx80g -Xms80g -jar \${picard_path} \
        MakeSitesOnlyVcf \
        INPUT=${variant_filtered_vcf_filename} \
        OUTPUT=${sites_only_vcf_filename}

    rm -rf ${workspace_dir_name} && touch ${workspace_dir_name}
    """
}

process ImportGVCFs {
    module 'singularity'

    input:
    path sample_name_map
    val unpadded_intervals
    each idx

    output:
    tuple val(idx), val(interval), path("${idx}.combined.g.vcf.gz"), path("${idx}.combined.g.vcf.gz.tbi")

    script:
    interval = unpadded_intervals.get(idx.toInteger())
    """
    CombineGVCFs() {
        inputMap=\$1
        outputFile=\$2
        ${gatk_path}/gatk CombineGVCFs --java-options "-Xmx80g -Xms80g" \
            -R \${ref_fasta} \
            -O "\${outputFile}" \
            -ip 500 \
            -L ${interval} \
            \$(echo \$(cat \${inputMap} | sed 's/^/ -V /g'))
    }

    oneGVCF=1
    i=0
    currentInput=${sample_name_map}
    while [[ \$oneGVCF -eq 1 ]]; do
        split -d -a1 -l${batch_size} --additional-suffix=.toCombine.list \${currentInput} "\${i}."
        currentList=\$(ls \${i}.*.toCombine.list)
        for file in \$currentList; do
            outputName="tmp.\$(echo \${file} | sed 's/.toCombine.list//').${idx}.combined.g.vcf.gz"
            echo \$outputName >> \${i}.${sample_name_map}
            CombineGVCFs \${file} \${outputName}
        done
        currentInput=\${i}.${sample_name_map}
        noGVCFs=\$(cat \$currentInput | wc -l)
        if [[ \$noGVCFs -eq 1 ]]; then
            oneGVCF=0
            mv \$(cat \$currentInput) ${idx}.combined.g.vcf.gz
            mv \$(cat \$currentInput).tbi ${idx}.combined.g.vcf.gz.tbi
            rm -f tmp*.vcf.gz*
            rm -f *toCombine.list
            rm -f *.${sample_name_map}
        fi
        i=\$((\$i+1))
    done
    """
}

process GenotypeGVCFs {
    input:
    tuple val(idx), val(interval), path(gvcfs), path(tbi)

    output:
    tuple val(idx), val(interval), path('*.output.vcf.gz'), path('*.output.vcf.gz.tbi')

    script:
    output_vcf_filename = "${idx}.output.vcf.gz"
    """
    ${gatk_path}/gatk --java-options "-Xmx5g -Xms5g" \
        GenotypeGVCFs \
        -R \${ref_fasta} \
        -O ${output_vcf_filename} \
        -D \${dbsnp_vcf} \
        -G StandardAnnotation \
        --only-output-calls-starting-in-intervals \
        --use-new-qual-calculator \
        -V \$(echo $gvcfs | sed 's/ / -V /g') \
        -L ${interval}
    """
}

process VariantFilteration {
    input:
    tuple val(idx), val(interval), path(vcf), path(tbi)

    output:
    tuple val(idx), val(interval), path('*variant_filtered.vcf.gz'), path('*variant_filtered.vcf.gz.tbi')

    script:
    excess_het_threshold = 54.69
    variant_filtered_vcf_filename = "${callset_name}.${idx}.variant_filtered.vcf.gz"
    """
    ${gatk_path}/gatk --java-options "-Xmx3g -Xms3g" \
        VariantFiltration \
        --filter-expression "ExcessHet > ${excess_het_threshold}" \
        --filter-name ExcessHet \
        -O ${variant_filtered_vcf_filename} \
        -V ${vcf}
    """
}

process MakeSitesOnlyVcf {
    input:
    tuple val(idx), val(interval), path(vcf), path(tbi)

    output:
    tuple path('*variant_filtered.sites_only.vcf.gz'), path('*variant_filtered.sites_only.vcf.gz.tbi'), emit: MakeSitesOnlyOuput

    script:
    sites_only_vcf_filename = "${callset_name}.${idx}.variant_filtered.sites_only.vcf.gz"
    """
    java -Xmx3g -Xms3g -jar \${picard_path} \
        MakeSitesOnlyVcf \
        INPUT=${vcf} \
        OUTPUT=${sites_only_vcf_filename}
    """
}

process GatherSiteOnlyVCFs {
    module 'singularity'
    module '/g/data/np30/users/nn6960/centos6.10/marcow/tabix/gcc-4.4.6/0.2.5'

    input:
    path(vcfandtbis)

    output:
    tuple path("${output_vcf_name}"), path("${output_vcf_name}.tbi")

    script:
    output_vcf_name = "${callset_name}.sites_only.vcf.gz"
    """
    export PATH="/g/data/np30/users/nn6960/tabix/tabix:$PATH"
    ls -1 *vcf.gz > inputs.list

    # --ignore-safety-checks makes a big performance difference so we include it in our invocation.
    # This argument disables expensive checks that the file headers contain the same set of
    # genotyped samples and that files are in order by position of first record.
    ${gatk_path}/gatk --java-options "-Xmx6g -Xms6g" \
        GatherVcfsCloud \
        --ignore-safety-checks \
        --gather-type BLOCK \
        --input inputs.list \
        --output ${output_vcf_name}

    tabix ${output_vcf_name}
    """
}

process IndelsVariantRecalibrator {
    module 'singularity'

    input:
    tuple path(vcf), path(tbi)

    output:
    tuple path("${recalibration_filename}"), path("${recalibration_filename}.idx"), path("${tranches_filename}")

    script:
    recalibration_filename = "${callset_name}.indels.recal"
    tranches_filename = "${callset_name}.indels.tranches"
    """
    ls \${mills_resource_vcf}
    ${gatk_path}/gatk --java-options "-Xmx24g -Xms24g" \
        VariantRecalibrator \
        -V ${vcf} \
        -O ${recalibration_filename} \
        --tranches-file ${tranches_filename} \
        --trust-all-polymorphic \
        -tranche \$(echo \${indel_recalibration_tranche_values} | sed 's/,/ -tranche /g') \
        -an \$(echo \${indel_recalibration_annotation_values} | sed 's/,/ -an /g') \
        -mode INDEL \
        --max-gaussians 4 \
        -resource:mills,known=false,training=true,truth=true,prior=12 \${mills_resource_vcf} \
        -resource:axiomPoly,known=false,training=true,truth=false,prior=10 \${axiomPoly_resource_vcf} \
        -resource:dbsnp,known=true,training=false,truth=false,prior=2 \${dbsnp_vcf}
    """
}

process SNPsVariantRecalibratorClassic {
    module 'singularity'

    input:
    tuple path(vcf), path(tbi)

    output:
    tuple path("${recalibration_filename}"), path("${recalibration_filename}.idx"), path("${tranches_filename}")

    script:
    recalibration_filename = "${callset_name}.snps.recal"
    tranches_filename = "${callset_name}.snps.tranches"
    """
    ${gatk_path}/gatk --java-options "-Xmx20g -Xms20g" \
        VariantRecalibrator \
        -V ${vcf} \
        -O ${recalibration_filename} \
        --tranches-file ${tranches_filename} \
        --trust-all-polymorphic \
        -tranche \$(echo \${snp_recalibration_tranche_values} | sed 's/,/ -tranche /g') \
        -an \$(echo \${snp_recalibration_annotation_values} | sed 's/,/ -an /g') \
        -mode SNP \
        --max-gaussians 6 \
        -resource:hapmap,known=false,training=true,truth=true,prior=15 \${hapmap_resource_vcf} \
        -resource:omni,known=false,training=true,truth=true,prior=12 \${omni_resource_vcf} \
        -resource:1000G,known=false,training=true,truth=false,prior=10 \${one_thousand_genomes_resource_vcf} \
        -resource:dbsnp,known=true,training=false,truth=false,prior=7 \${dbsnp_vcf}
    """
}


process ApplyRecalibration {
    module 'singularity'

    input:
    tuple val(idx), val(interval), path(vcf), path(tbi)
    tuple path(snps_recalibration), path(snps_ecalibration_index), path(snps_tranches)
    tuple path(indels_recalibration), path(indels_ecalibration_index), path(indels_tranches)

    output:
    tuple path("${recalibrated_vcf_filename}"), path("${recalibrated_vcf_filename}.tbi")

    script:
    recalibrated_vcf_filename = "${callset_name}.filtered.${idx}.vcf.gz"
    """
    ${gatk_path}/gatk --java-options "-Xmx5g -Xms5g" \
        ApplyVQSR \
        -O tmp.indel.recalibrated.vcf \
        -V ${vcf} \
        --recal-file ${indels_recalibration} \
        --tranches-file ${indels_tranches} \
        --truth-sensitivity-filter-level \${indel_filter_level} \
        --create-output-variant-index true \
        -mode INDEL

    ${gatk_path}/gatk --java-options "-Xmx5g -Xms5g" \
        ApplyVQSR \
        -O ${recalibrated_vcf_filename} \
        -V tmp.indel.recalibrated.vcf \
        --recal-file ${snps_recalibration} \
        --tranches-file ${snps_tranches} \
        --truth-sensitivity-filter-level \${snp_filter_level} \
        --create-output-variant-index true \
        -mode SNP
    """
}

process FinalGather {

    publishDir 'output/', mode: "move"

    module 'singularity'

    input:
    path(vcfandtbis)

    output:
    tuple path("${output_vcf_name}"), path("${output_vcf_name}.tbi"), emit: FinalGatherOutput

    script:
    output_vcf_name = "${callset_name}.vcf.gz"
    """
    export PATH="/g/data/np30/users/nn6960/tabix/tabix:$PATH"
    ls -1 *vcf.gz > inputs.list

    # --ignore-safety-checks makes a big performance difference so we include it in our invocation.
    # This argument disables expensive checks that the file headers contain the same set of
    # genotyped samples and that files are in order by position of first record.
    ${gatk_path}/gatk --java-options "-Xmx6g -Xms6g" \
        GatherVcfsCloud \
        --ignore-safety-checks \
        --gather-type BLOCK \
        --input inputs.list \
        --output ${output_vcf_name}

    tabix ${output_vcf_name}
    """
}

workflow {
    // get merge count
    possible_merge_count = Math.round(Math.floor(num_of_original_intervals / num_gvcfs / 2.5))
    merge_count = 1
    if (possible_merge_count > 1) {
        merge_count = possible_merge_count
    }

    // run python script for dynamically combine intervals
    DynamicallyCombineIntervals(unpadded_intervals_file, merge_count)

    // get range from the dynamically combine interval size
    range = getIntervalLength(DynamicallyCombineIntervals.out.output_intervals.readLines().size())

    idRange = range.splitText() { it -> it.trim() }

    // Collect MakeSitesOnlyOutput with two different methods based on params.splitJaabs
    if (! params.splitJaabs) {
        CombineGenotypeSitesOnlyGVCF(sample_name_map, DynamicallyCombineIntervals.out.output_intervals.readLines(), idRange.collect())
        CombineGenotypeSitesOnlyGVCF.out.VariantFilterationOutput.set { ApplyRecalibrationVCF }
        CombineGenotypeSitesOnlyGVCF.out.MakeSitesOnlyOutput.set { MakeSitesOnlyOuput }
    } else {
        import_gvcf_out = ImportGVCFs(sample_name_map, DynamicallyCombineIntervals.out.output_intervals.readLines(), idRange.collect())

        import_gvcf_out.groupTuple(by:[0,1]).set { GenotypeGVCFInput }
        genotype_gvcf_output = GenotypeGVCFs(GenotypeGVCFInput)

        genotype_gvcf_output.set { VariantFilterationInput }
        variant_filteration_output = VariantFilteration(VariantFilterationInput)

        variant_filteration_output.set { ApplyRecalibrationVCF }

        MakeSitesOnlyVcf(variant_filteration_output)
        MakeSitesOnlyVcf.out.MakeSitesOnlyOuput.set { MakeSitesOnlyOuput }
    }

    // run Gather Vcfs file
    gather_site_only_output = GatherSiteOnlyVCFs(MakeSitesOnlyOuput.collect())

    // run recalibrate indel variant
    indel_output = IndelsVariantRecalibrator(gather_site_only_output)

    // run recalibrate SNPs variant
    snp_output = SNPsVariantRecalibratorClassic(gather_site_only_output)

    // Apply VQSR from the recalibrated SNPs and INDELs
    apply_recal_output = ApplyRecalibration(ApplyRecalibrationVCF, snp_output, indel_output)

    // run final gather on the recalibrated output
    FinalGather(apply_recal_output.collect())
}
