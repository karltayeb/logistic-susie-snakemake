rule render_ukbb_sim_report:
    input:
        pip_summary = 'results/ukbb_geno/chr1_100794065_101983457/pip_summary.rds',
        cs_summary = 'results/ukbb_geno/chr1_100794065_101983457/cs_summary.rds' 
    output: 'docs/ukbb_simulation_report.html'
    run:
        shell("Rscript -e \"rmarkdown::render('workflow/notebooks/ukbb_simulation_report.Rmd', output_dir = 'docs')\"")
