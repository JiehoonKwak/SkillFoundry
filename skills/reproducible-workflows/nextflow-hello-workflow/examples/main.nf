nextflow.enable.dsl = 2

params.outdir = "${workflow.projectDir}/results"

Channel
    .of('Hello', 'Bonjour', 'Hola', 'Ciao')
    .set { greetings }

process sayHello {
    publishDir params.outdir, mode: 'copy'

    input:
    val greeting

    output:
    path "${greeting}.txt"

    script:
    """
    echo '${greeting} world!' > ${greeting}.txt
    """
}

workflow {
    sayHello(greetings)
}
