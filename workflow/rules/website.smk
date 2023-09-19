# note render rule does not look for inputs!
# make sure all the necessary inputs are complete
# mimics workflowr report generating by committing rmd and html at once.
rule render_notebook:
    params: report = '{report}'
    output: 'docs/{report}.html'
    run:
        shell("Rscript -e \"rmarkdown::render('workflow/notebooks/{params.report}.Rmd', output_dir = 'docs')\"")
        shell("git commit workflow/notebooks/{params.report}.Rmd docs/{params.report}.html -m 'rendered {params.report}'")

REPORTS = [
    'index',
    'gsea_simulation_report',
    'ukbb_simulation_report'
]

rule render_website:
    input: 
        expand('docs/{report}.html', report = REPORTS)
        
# rule render_ukbb_sim_report:
#     input:
#         pip_summary = 'results/ukbb_geno/chr1_100794065_101983457/pip_summary.rds',
#         cs_summary = 'results/ukbb_geno/chr1_100794065_101983457/cs_summary.rds',
#         rmd = 'workflow/notebooks/ukbb_simulation_report.Rmd' 
#     output: 'docs/ukbb_simulation_report.html'
#     run:
#         shell("Rscript -e \"rmarkdown::render('workflow/notebooks/ukbb_simulation_report.Rmd', output_dir = 'docs')\"")

# rule render_gsea_sim_report:
#     input:
#         pip_summary = 'results/gsea/X_bin_strong/linear_susie_default_cs_summary.rds',
#         cs_summary = 'results/gsea/X_bin_strong/linear_susie_default_pip_summary.rds',
#         # rmd = 'workflow/notebooks/gsea_simulation_report.Rmd'
#     output: 'docs/gsea_simulation_report.html'
#     run:
#         shell("Rscript -e \"rmarkdown::render('workflow/notebooks/gsea_simulation_report.Rmd', output_dir = 'docs')\"")
