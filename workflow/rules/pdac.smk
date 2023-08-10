import glob

extract_filename = lambda path: path.split('/')[-1].split('.')[0]
pdac_factor_paths = glob.glob('resources/yusha_sc_tumor/pdac/factor*')
PDAC_FACTORS = [extract_filename(path) for path in pdac_factor_paths]

rule load_pdac_data:
    output: 'results/pdac/pdac_data.rds'
    shell: "Rscript workflow/scripts/pdac/load_pdac_data.R"

rule load_pdac_msigdb_X:
    input: rules.load_pdac_data.output,
    output: 'results/pdac/pdac_msigdb_X.rds',
    shell: "Rscript workflow/scripts/pdac/load_pdac_msigdb_X.R"

rule pdac_fit_and_sim_lasso:
    input: 
        rules.load_pdac_data.output,
        rules.load_pdac_msigdb_X.output
    output:
        'results/pdac/pdac_lasso_fits.rds',
        'results/pdac/pdac_lasso_sims.rds'
    shell: "Rscript workflow/scripts/pdac/fit_and_sim_lasso.R"

rule pdac_r_to_py:
    input:
        rules.load_pdac_data.output,
        rules.load_pdac_msigdb_X.output,
        rules.pdac_fit_and_sim_lasso.output
    output:
        'results/pdac/pdac_data_py.pkl'
    shell: "Rscript workflow/scripts/pdac/pdac_r_to_py.R"

rule pdac_gibss_fit:
    input: rules.pdac_r_to_py.output
    output: 
        fit = 'results/pdac/fits/{factor}_gibss.pkl'
    script: '../scripts/pdac/pdac_fit.py'

rule pdac_gibss_sim:
    input: rules.pdac_r_to_py.output
    output: 
        simfit = 'results/pdac/sims/{factor}_{sim_id}_gibss.pkl'
    script: '../scripts/pdac/pdac_sim.py'

rule pdac_fit_and_sim:
    input: 
        expand("results/pdac/sims/{{factor}}_{sim_id}_gibss.pkl", sim_id=range(20)),
        expand("results/pdac/fits/{{factor}}_gibss.pkl", sim_id=range(20))
    output: 'results/pdac/pdac_fit_and_sim_{factor}.txt'
    shell: "ran logistic SuSiE via GIBSS for PDAC factors, and lasso simulations > {output}"

rule pdac_summarize_gibss:
    input: expand("results/pdac/sims/{factor}_{sim_id}_gibss.pkl", factor=PDAC_FACTORS, sim_id=range(20)),
    output: "results/pdac/sim_summary.pkl"
    script: "scripts/pdac_sim_summary.py"

rule pdac_lasso_lasso:
    input: "results/pdac/pdac_lasso_sims.rds"
    output: "results/pdac/sim_lasso_fits/{factor}_lasso_lasso.rds"
    script: "../scripts/pdac/pdac_lasso_lasso.R"

rule pdac_lasso_lasso2:
    input: expand("results/pdac/sim_lasso_fits/{factor}_lasso_lasso.rds", factor=PDAC_FACTORS)
    output: 'results/pdac/pdac_lasso_lasso.txt'
    shell: "fit lasso to lasso simulations > {output}"

rule pdac_all:
    input: 
        expand("results/pdac/sims/{factor}_{sim_id}_gibss.pkl", factor=PDAC_FACTORS, sim_id=range(20)),
        expand("results/pdac/fits/{factor}_gibss.pkl", factor=PDAC_FACTORS, sim_id=range(20)),
        expand("results/pdac/sim_lasso_fits/{factor}_lasso_lasso.rds", factor=PDAC_FACTORS)
    output: 'results/pdac/pdac_all.txt'
    shell: "ran logistic SuSiE via GIBSS for PDAC factors, and lasso simulations > {output}"
