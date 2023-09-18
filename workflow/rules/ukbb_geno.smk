configfile: "config/config.yaml"

wildcard_constraints:
    param_id="[^_]+"

rule make_genotype_matrix:
    output: 
        'results/ukbb_geno/{region}/genotype.pkl',
        'results/ukbb_geno/{region}/genotype.rds'
    script: "../scripts/ukbb_geno/make_genotype_matrix_py.R"

rule make_ld_matrix:
    input: 'results/ukbb_geno/{region}/genotype.rds'
    output: 'results/ukbb_geno/{region}/ld.rds'
    script: "../scripts/ukbb_geno/make_ld.R"

rule simulate_case_control:
    input: 
        genotype = 'results/ukbb_geno/{region}/genotype.rds'
        #'config/ukbb_sim/ukbb_sim_manifest.tsv'
    output: 
        rds = 'results/ukbb_geno/{region}/{sim_id}/sim.rds',
        pkl = 'results/ukbb_geno/{region}/{sim_id}/sim.pkl'
    script: "../scripts/ukbb_geno/simulate_case_control.R"
    
rule ukbb_geno_fit_ser:
    input: 
        genotype='results/ukbb_geno/{region}/genotype.pkl',
        sim='results/ukbb_geno/{region}/{sim_id}/sim.pkl'
    params:
        fit_params = lambda w: config['logistic_ser_params'].get(w.param_id, {})
    output:
        'results/ukbb_geno/{region}/{sim_id}/ser_{param_id}.pkl'
    script: "../scripts/ukbb_geno/ukbb_geno_fit_ser.py" 

rule ukbb_geno_fit_susie_abf:
    input: 
        genotype='results/ukbb_geno/{region}/genotype.pkl',
        sim='results/ukbb_geno/{region}/{sim_id}/sim.pkl'
    params:
        fit_params = lambda w: config['logistic_ibss_abf_params'].get(w.param_id, {})
    output:
        'results/ukbb_geno/{region}/{sim_id}/gibss_abf_{param_id}.pkl'
    script: "../scripts/ukbb_geno/ukbb_geno_fit_susie_abf.py"

rule ukbb_geno_fit_susie_labf:
    input: 
        genotype='results/ukbb_geno/{region}/genotype.pkl',
        sim='results/ukbb_geno/{region}/{sim_id}/sim.pkl'
    params:
        fit_params = lambda w: config['logistic_ibss_labf_params'].get(w.param_id, {})
    output:
        'results/ukbb_geno/{region}/{sim_id}/gibss_labf_{param_id}.pkl'
    script: "../scripts/ukbb_geno/ukbb_geno_fit_susie_labf.py" 

rule ukbb_geno_fit_rss:
    input:
        ld_path = 'results/ukbb_geno/{region}/ld.rds',
        ser_path = 'results/ukbb_geno/{region}/{sim_id}/ser_default.pkl',
        sim_path = 'results/ukbb_geno/{region}/{sim_id}/sim.rds'
    params:
        fit_params = lambda w: config['logistic_rss_params'].get(w.param_id, {})
    output:
        'results/ukbb_geno/{region}/{sim_id}/rss_{param_id}.rds'
    script:
        '../scripts/ukbb_geno/ukbb_geno_fit_susie_rss.R'

rule ukbb_geno_fit_linear_susie:
    input:
        genotype_path = 'results/ukbb_geno/{region}/genotype.rds',
        sim_path = 'results/ukbb_geno/{region}/{sim_id}/sim.rds'
    params:
        fit_params = lambda w: config['linear_susie_params'].get(w.param_id, {})
    output:
        'results/ukbb_geno/{region}/{sim_id}/linear_susie_{param_id}.rds'
    script:
        '../scripts/ukbb_geno/ukbb_geno_fit_linear_susie.R'


import pandas as pd
ukbb_sim_manifest = pd.read_csv('config/ukbb_sim/ukbb_main_manifest.tsv', sep='\t')
rule run_ukbb_sim:
    input:
        expand([
            'results/ukbb_geno/{region}/{sim_id}/gibss_labf_default.pkl',
            'results/ukbb_geno/{region}/{sim_id}/gibss_abf_default.pkl',
            'results/ukbb_geno/{region}/{sim_id}/ser_default.pkl', 
            'results/ukbb_geno/{region}/{sim_id}/rss_L1-eb.rds',
            'results/ukbb_geno/{region}/{sim_id}/rss_L5-eb.rds',
            'results/ukbb_geno/{region}/{sim_id}/linear_susie_default.rds',
        ], zip, 
        region = ukbb_sim_manifest.region, 
        sim_id = ukbb_sim_manifest.sim_id)
    output: 'results/ukbb_geno/fit.txt'
    shell: "echo 'ran simulations' > {output}"

rule simulate_case_control2:
    input:
        expand([
            'results/ukbb_geno/{region}/{sim_id}/sim.rds',
        ], zip, 
        region = ukbb_sim_manifest.region, 
        sim_id = ukbb_sim_manifest.sim_id)
    output: 'results/ukbb_geno/sim.txt'
    shell: "echo 'generated simulations' > {output}"

rule summarize_simulations:
    input: 
        rules.run_ukbb_sim.input,
        ld = 'results/ukbb_geno/{region}/ld.rds'
    output: 
        pip_summary = 'results/ukbb_geno/{region}/pip_summary.rds',
        cs_summary = 'results/ukbb_geno/{region}/cs_summary.rds'
    script: "../scripts/ukbb_geno/ukbb_pip_and_cs_summary.R"


ukbb_ser_manifest = pd.read_csv('config/ukbb_sim/ukbb_ser_manifest.tsv', sep='\t')
rule run_ukbb_ser_sim:
    input:
        expand([
            'results/ukbb_geno/{region}/{sim_id}/ser_default.pkl', 
            'results/ukbb_geno/{region}/{sim_id}/rss_ukbb-L1-default.rds',
            'results/ukbb_geno/{region}/{sim_id}/rss_ukbb-L1-eb.rds',
            'results/ukbb_geno/{region}/{sim_id}/linear_susie_L1.rds',
        ], zip, 
        region = ukbb_ser_manifest.region, 
        sim_id = ukbb_ser_manifest.sim_id)
    output: 'results/ukbb_geno/fit_ser_sim.txt'
    shell: "echo 'ran simulations' > {output}"


