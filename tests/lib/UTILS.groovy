// Helper functions for pipeline tests
class UTILS {

    public static def get_assertion = { Map args ->
        // Mandatory, as we always need an outdir
        def outdir = args.outdir

        // Use this args to run the test with stub
        // It will disable all assertions but versions and stable_name
        def stub = args.stub

        // Will print the summary instead of the md5sum for vcf files
        def no_vcf_md5sum = args.no_vcf_md5sum

        // stable_name: All files + folders in ${outdir}/ with a stable name
        def stable_name = getAllFilesFromDir(outdir, relative: true, includeDir: true, ignore: ['pipeline_info/*.{html,json,txt}'])
        // stable_content: All files in ${outdir}/ with stable content
        def stable_content = getAllFilesFromDir(outdir, ignoreFile: 'tests/.nftignore')
        // vcf_files: All vcf files
        def vcf_files = getAllFilesFromDir(outdir, include: ['**/*.vcf{,.gz}'])

        def assertion = []

        assertion.add(removeFromYamlMap("${outdir}/pipeline_info/nf_core_variantprioritization_software_mqc_versions.yml", "Workflow"))
        assertion.add(stable_name)

        if (!stub) {
            assertion.add(stable_content.isEmpty() ? 'No stable content' : stable_content)
            if (no_vcf_md5sum) {
                assertion.add(vcf_files.isEmpty() ? 'No VCF files' : vcf_files.collect { file -> [ file.getName(), path(file.toString()).vcf.summary ] })
            } else {
                assertion.add(vcf_files.isEmpty() ? 'No VCF files' : vcf_files.collect { file -> file.getName() + ":md5," + path(file.toString()).vcf.variantsMD5 })
            }
        }

        return assertion
    }

    public static def get_test = { scenario ->
        // This function returns a closure that will be used to run the test and the assertion
        // It will create tags or options based on the scenario

        return {

            if (scenario.stub) {
                options "-stub"
            }
            // If a tag is provided, add it to the test
            if (scenario.tag) {
                tag scenario.tag
            }
            when {
                params {
                    // Mandatory, as we always need an outdir
                    outdir = "${outputDir}"
                    // Apply scenario-specific params
                    scenario.params.each { key, value ->
                        delegate."$key" = value
                    }
                }
            }

            then {
                // Assert failure
                if (scenario.failure) {
                    // Early failure, so we don't pollute console with massive diffs
                    assert workflow.failed
                    // Check stdout if specified
                    if (scenario.stdout) {
                        assertAll(
                            { assert workflow.stdout.toString().contains(scenario.stdout) }
                        )
                    }
                    // Check stderr if specified
                    if (scenario.stderr) {
                        { assert snapshot(
                            workflow.stderr.toString().replaceAll(/\x1B\[[0-9;]*m/, '').replaceAll(/^\[/, '').replaceAll(/\]$/, '').replaceAll(/, /, ',').split(",").findAll { !it.matches(/.*Nextflow [0-9]+\.[0-9]+\.[0-9]+ is available.*/) }[scenario.stderr]
                        ).match() }
                    }
                // Assert success
                } else {
                    // Early failure, so we don't pollute console with massive diffs
                    assert workflow.success
                    assertAll(
                        { assert snapshot(
                            // Number of successful tasks
                            workflow.trace.succeeded().size(),
                            // All assertions based on the scenario
                            *UTILS.get_assertion(no_vcf_md5sum: scenario.no_vcf_md5sum, outdir: params.outdir, stub: scenario.stub)
                        ).match() }
                    )
                    // Check stdout if specified
                    if (scenario.stdout) {
                        assert workflow.stdout.toString().contains(scenario.stdout)
                    }
                }
            }
        }
    }
}