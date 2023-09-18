# rules for simulating and analyzing GSEA examples 
configfile: "config/config.yaml"

rule msigdb_X:
    output: 'results/gsea/msigdb/X.rds'

rule sparse_X:
    params:
        params = lambda w: config['sparse_X'].get(w.sparse_x_id)
    output: 
        'results/gsea/{sparse_x_id}/X.rds',
        'results/gsea/{sparse_x_id}/ld.rds'
    script: "../scripts/gsea/simulate_sparse_X.R"

rule simulate_gene_list:
    input: 
        genotype = 'results/gsea/{x_id}/X.rds'
    output: 
        rds = 'results/gsea/{x_id}/{sim_id}/sim.rds'
    script: "../scripts/gsea/simulate_gene_list.R"

rule gsea_fit_gibss_labf:
    input: 
        X = 'results/gsea/{x_id}/X.rds',
        sim = 'results/gsea/{x_id}/{sim_id}/sim.rds'
    params:
        fit_params = lambda w: config['logistic_ibss_labf_params'].get(w.param_id, {})
    output: 
        fit = 'results/gsea/{x_id}/{sim_id}/gibss_labf_{param_id}.rds'
    script: "../scripts/fitting/fit_gibss_labf.R"

rule gsea_fit_gibss_abf:
    input: 
        X = 'results/gsea/{x_id}/X.rds',
        sim = 'results/gsea/{x_id}/{sim_id}/sim.rds'
    params:
        fit_params = lambda w: config['logistic_ibss_abf_params'].get(w.param_id, {})
    output: 
        fit = 'results/gsea/{x_id}/{sim_id}/gibss_abf_{param_id}.rds'
    script: "../scripts/fitting/fit_gibss_abf.R"

rule gsea_fit_linear_susie:
    input: 
        X = 'results/gsea/{x_id}/X.rds',
        sim = 'results/gsea/{x_id}/{sim_id}/sim.rds'
    params:
        fit_params = lambda w: config['linear_susie_params'].get(w.param_id, {})
    output: 
        fit = 'results/gsea/{x_id}/{sim_id}/linear_susie_{param_id}.rds'
    script: "../scripts/fitting/fit_linear_susie.R"


import pandas as pd
gsea_manifest = pd.read_csv('config/gsea_sim/gsea_sim_manifest.tsv', sep='\t')
METHODS = ['gibss_abf_default', 'gibss_labf_default', 'gibss_abf_L1', 'gibss_labf_L1', 'linear_susie_default', 'linear_susie_L1'] 
gsea_fit_manifest = pd.merge(gsea_manifest, pd.DataFrame({'method': METHODS}), how='cross')
rule run_gsea_simulation:
    input:
        expand([
            'results/gsea/{x_id}/{sim_id}/{method}.rds'
        ], zip, 
        x_id = gsea_fit_manifest.X, 
        sim_id = gsea_fit_manifest.sim_id,
        method=gsea_fit_manifest.method)
    output: 'results/gsea/run.txt'
    shell: "echo 'ran simulations' > {output}"

rule summarize_gsea_simulation:
    input:
        manifest = 'config/gsea_sim/gsea_sim_manifest.tsv',
        ld = 'results/gsea/{x_id}/ld.rds',
        sims = expand('results/gsea/{x_id}/{sim_id}/sim.rds', zip, 
            x_id = gsea_manifest.X, 
            sim_id = gsea_manifest.sim_id),
        fits = expand('results/gsea/{x_id}/{sim_id}/{{model}}.rds', zip, 
            x_id = gsea_manifest.X, 
            sim_id = gsea_manifest.sim_id)
    output: 
        cs_summary = 'results/gsea/{x_id}/{model}_cs_summary.rds',
        pip_summary = 'results/gsea/{x_id}/{model}_pip_summary.rds'
    script: "../scripts/fitting/summarize.R"

rule summarize_gsea_simultions:
    input: expand('results/gsea/{{x_id}}/{model}_cs_summary.rds', model=METHODS)
    output: 'results/gsea/{x_id}/summarized.txt'
    shell: "echo 'ran rule summarize_gsea_simulations' > {output}"