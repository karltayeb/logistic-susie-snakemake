import re
import pandas as pd

SIMS = pd.read_csv('config/rss_comparison/parameter_table.csv')\
    .set_index('sim_id')\
    .to_dict(orient='index')

unique_pfiles = list(set([SIMS[x]['pfile'] for x in SIMS]))
region2pfile = {pfile.split('/')[-1]: pfile for pfile in unique_pfiles}

def get_sim_params(sim_id):
    return SIMS.get(sim_id)

def get_pgen_path(wildcards):
    pfile = get_sim_params(wildcards.sim_id).get('pfile')
    return re.sub(r'\s+', '', f"{pfile}.pgen")

def get_ld_path(wildcards):
    region = get_sim_params(wildcards.sim_id)['pfile'].split('/')[-1]
    path = f"results/rss_comparison/ld/{region}.unphased.vcor1"
    return re.sub(r'\s+', '', path)


rule simulate_case_control_from_pgen:
    output:
        rds="results/rss_comparison/simulation/{sim_id}/{sim_id}.sim.rds",
        pheno="results/rss_comparison/simulation/{sim_id}/{sim_id}.sim.pheno"
    params:
        sim_params = lambda wildcards, input: get_sim_params(wildcards.sim_id)
    script:
        "../scripts/rss_comparison/simulate_case_control.R"


rule run_association_study_via_plink:
    input:
        pgen = get_pgen_path,
        pheno="results/rss_comparison/simulation/{sim_id}/{sim_id}.sim.pheno"
    params:
        pfile = lambda wildcards, input:  input.pgen[:-5]
    output:
        linear = 'results/rss_comparison/simulation/{sim_id}/{sim_id}.Continuous.glm.linear',
        logistic = 'results/rss_comparison/simulation/{sim_id}/{sim_id}.Binary.glm.logistic.hybrid'
    shell:
        """
        plink2 --pfile {params.pfile} \
          --pheno {input.pheno} \
          --glm allow-no-covars \
          --threads 12
          --out results/rss_comparison/simulation/{wildcards.sim_id}/{wildcards.sim_id}
        """

rule ld_matrix:
    params:
        pfile = lambda wildcards, input: region2pfile.get(wildcards.region) 
    output:
        ld_matrix="results/rss_comparison/ld/{region}.unphased.vcor1"
    shell:
        """
        plink2 --pfile {params.pfile} \
          --r-unphased square \
          --threads 12 \
          --out results/rss_comparison/ld/{wildcards.region}
        """

rule ld_matrix_all:
    output:
      expand('results/rss_comparison/ld/{region}.unphased/vcor1', region=list(region2pfile.keys()))
    shell:
        """
        touch results/asthma/ld/get_all_ld.txt
        """

rule fit_rss_linear:
    input:
        sumstat_path = "results/rss_comparison/simulation/{sim_id}/{sim_id}.Continuous.glm.linear",
        ld_path = get_ld_path
    output:
        rss = "results/rss_comparison/simulation/{sim_id}/{sim_id}.linear_rss_fit.rds"
    script:
        "../scripts/rss_comparison/fit_rss.R"


rule fit_rss_logistic:
    input:
        sumstat_path = "results/rss_comparison/simulation/{sim_id}/{sim_id}.Binary.glm.logistic.hybrid",
        ld_path = get_ld_path
    output:
        rss = "results/rss_comparison/simulation/{sim_id}/{sim_id}.logistic_rss_fit.rds"
    script:
        "../scripts/rss_comparison/fit_rss.R"


rule fit_all:
    input:
        expand('results/rss_comparison/simulation/{sim_id}/{sim_id}.logistic_rss_fit.rds', sim_id = list(SIMS.keys()))
    output:
        "results/rss_comparison/simulation/fit.txt"
    shell:
        """
        touch {output}
        """

rule pip_calibration:
    input:
        logistic_fits = expand('results/rss_comparison/simulation/{sim_id}/{sim_id}.logistic_rss_fit.rds', sim_id = list(SIMS.keys())),
        sims = expand('results/rss_comparison/simulation/{sim_id}/{sim_id}.sim.rds', sim_id = list(SIMS.keys()))
    output:
        pip = "results/rss_comparison/simulation/pip_calibration.pdf",
        cs = "results/rss_comparison/simulation/cs_coverage.pdf"
    script:
        "../scripts/rss_comparison/pip_calibration.R"
