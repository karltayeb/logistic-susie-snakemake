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
    input: 'results/ukbb_geno/{region}/genotype.rds'
    output: 
        'results/ukbb_geno/{region}/sim.rds',
        'results/ukbb_geno/{region}/sim.pkl'
    script: "../scripts/ukbb_geno/simulate_case_control.R"
    
rule ukbb_geno_simulate_ser:
    input: 
        'results/ukbb_geno/{region}/genotype.pkl'
    output:
        'results/ukbb_geno/{region}/ser_sims.pkl'
    script: "../scripts/ukbb_geno/ukbb_geno_ser_sim.py"

rule ukbb_geno_fit_ser:
    input: 
        'results/ukbb_geno/{region}/genotype.pkl',
        'results/ukbb_geno/{region}/sim.pkl'
    output:
        'results/ukbb_geno/{region}/{sim_id}/ser.pkl'
    script: "../scripts/ukbb_geno/ukbb_geno_fit_ser.py" 

rule ukbb_geno_fit_susie_abf:
    input: 
        'results/ukbb_geno/{region}/genotype.pkl',
        'results/ukbb_geno/{region}/sim.pkl'
    output:
        'results/ukbb_geno/{region}/{sim_id}/gibss_abf.pkl'
    script: "../scripts/ukbb_geno/ukbb_geno_fit_susie_abf.py"

rule ukbb_geno_fit_susie_labf:
    input: 
        'results/ukbb_geno/{region}/genotype.pkl',
        'results/ukbb_geno/{region}/sim.pkl'
    output:
        'results/ukbb_geno/{region}/{sim_id}/gibss.pkl'
    script: "../scripts/ukbb_geno/ukbb_geno_fit_susie_labf.py" 

rule ukbb_geno_fit_rss_estimate_prior_variance:
    input:
        ld_path = 'results/ukbb_geno/{region}/ld.rds',
        ser_path = 'results/ukbb_geno/{region}/{sim_id}/ser.pkl'
    params:
        n = 50000,
        estimate_prior_variance = 'true'
    output:
        'results/ukbb_geno/{region}/{sim_id}/rss_eb.rds'
    script:
        '../scripts/ukbb_geno/ukbb_geno_fit_susie_rss.R'

rule ukbb_geno_fit_rss:
    input:
        ld_path = 'results/ukbb_geno/{region}/ld.rds',
        ser_path = 'results/ukbb_geno/{region}/{sim_id}/ser.pkl'
    params:
        n = 50000,
        estimate_prior_variance = 'false'
    output:
        'results/ukbb_geno/{region}/{sim_id}/rss.rds'
    script:
        '../scripts/ukbb_geno/ukbb_geno_fit_susie_rss.R'

rule ukbb_geno_fit_linear_susie:
    input:
        genotype_path = 'results/ukbb_geno/{region}/genotype.rds',
        sim_path = 'results/ukbb_geno/{region}/sim.rds'
    params:
        estimate_prior_variance = 'false'
    output:
        'results/ukbb_geno/{region}/{sim_id}/linear_susie.rds'
    script:
        '../scripts/ukbb_geno/ukbb_geno_fit_linear_susie.R'

import pandas as pd
ukbb_sim_manifest = pd.read_csv('config/ukbb_sim/ukbb_sim_manifest.tsv', sep='\t')
rule run_ukbb_sim:
    input:
        expand([
            'results/ukbb_geno/{region}/{sim_id}/gibss.pkl',
            'results/ukbb_geno/{region}/{sim_id}/gibss_abf.pkl',
            'results/ukbb_geno/{region}/{sim_id}/ser.pkl', 
            'results/ukbb_geno/{region}/{sim_id}/rss.rds',
            'results/ukbb_geno/{region}/{sim_id}/rss_eb.rds'
        ], zip, 
        region = ukbb_sim_manifest.region, 
        sim_id = ukbb_sim_manifest.sim_id)
    output: 'results/ukbb_geno/fit.txt'
    shell: "echo 'ran simulations' > {output}"