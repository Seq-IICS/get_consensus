import groovy.json.JsonSlurper

def jsonSlurper = new JsonSlurper()

params.debug = false
params.quality = false

def configFile = new File("${params.path}/config.json")
String configJSON = configFile.text
def myConfig = jsonSlurper.parseText(configJSON)

timestamp = workflow.start
treshold = 0.5
depth = 600000
min_map_q = 20

if (myConfig.treshold){
    treshold = myConfig.treshold
}
if (myConfig.min_map_q) {
    min_map_q = myConfig.min_map_q
}
if (myConfig.depth) {
    depth = myConfig.depth
}

process nanoplotWrapper {
    input:
        val myConfig
        path dataPath
    output:
        stdout
    script:
    """
    nanoplot_wrapper.py \
        -p ${dataPath}/ \
        -s ${myConfig.sample} \
        -t ${timestamp.format("dd-MM-yyyy_HH_mm_ss")} \
        -d ${params.debug} 
    """
}

process getResults {
    input:
        tuple val(x);
        path dataPath;
    output:
        stdout
    script:
    """
    myStringArray=("${x}")
    myArray=\$(echo \$myStringArray | tr -d "[],''" )
    echo \$myArray
    outputDir=${dataPath}/output/results_${timestamp.format("dd-MM-yyyy_HH_mm_ss")}/cats
    cp ${dataPath}/${myConfig.reference} \$outputDir
    cd \$outputDir

    for i in \${myArray[@]}; do
        echo "Processing \$i"
        minimap2 -a ${myConfig.reference} \${i}.gz > \${i}.sam
        samtools view \${i}.sam > \${i}.bam
        samtools sort \${i}.sam > \${i}_sorted.bam
        samtools index \${i}_sorted.bam
        #samtools mpileup -uf ${myConfig.reference} \${i}_sorted.bam | bcftools call -c | vcfutils.pl vcf2fq > \${i}_cns.fastq
        
        if [[ "${params.type}" == "l" ]];
        then
            echo "consensus long"
            bcftools mpileup -f ${myConfig.reference} \${i}_sorted.bam | bcftools call -c | vcfutils.pl vcf2fq > \${i}_cns.fastq
            seqtk seq -aQ64 -q 20 -n N \${i}_cns.fastq > \${i}_dirty.fasta
            tail -n +2 "\${i}_dirty.fasta" > \${i}_cns.fasta
            ## ELIMINA LA PRIMERA LINEA Y LE AGREGA LA CORRECTA
            sed -i "1i >\${i}" \${i}_cns.fasta
            rm \${i}_dirty.fasta
        elif [[ "${params.type}" == "c" ]];
        then
            echo "consensus"
            samtools mpileup -d ${depth} -A -Q 0 \${i}_sorted.bam | ivar consensus -p \${i}_cns.fasta -q ${min_map_q}  -t ${treshold} -n N
        fi

        if [[ "${params.quality}" == "true" ]];
        then
            echo "quality"
            echo "### CRAMINO ###" > \${i}.crami
            cramino \${i}_sorted.bam >> \${i}.crami
            echo "\n### SAMTOOLS DEPTH ####\n" >> \${i}.crami
            samtools depth \${i}_sorted.bam >>\${i}.crami
            echo "\n### SAMTOOLS COVERAGE ####\n" >> \${i}.crami
            samtools coverage \${i}_sorted.bam >> \${i}.crami
        fi
    done
    """
}

workflow {
    dataPath = Channel.fromPath(params.path);
    lstFiles = nanoplotWrapper(myConfig, dataPath);
    getResults(lstFiles, dataPath) | view;
}
