import re
import pandas as pd

SIMS = pd.read_csv('config/rss_comparison/parameter_table.csv')\
    .set_index('sim_id')\
    .to_dict(orient='index')

def get_sim_params(sim_id):
    return SIMS.get(sim_id)

def get_pgen_path(wildcards):
    return re.sub(r'\s+', '', f"{get_sim_params(wildcards.sim_id).get('pfile')}.pgen")

def get_pfile_prefix(wildcards):
    return re.sub(r'\s+', '', f"resources/ukb-geno/{wildcards.region}")


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
        prefix = lambda wildcards, input:  input.pgen[:-5]
    output:
        linear = 'results/rss_comparison/simulation/{sim_id}/{sim_id}.Continuous.glm.linear',
        logistic = 'results/rss_comparison/simulation/{sim_id}/{sim_id}.Binary.glm.logistic.hybrid'
    shell:
        """
        plink2 --pfile {params.prefix} --pheno {input.pheno} --glm allow-no-covars --out results/rss_comparison/simulation/{wildcards.sim_id}/{wildcards.sim_id}
        """

rule ld_matrix:
    params:
        pfile = get_pfile_prefix
    output:
        ld_matrix="results/rss_comparison/ld/{region}.unphased.vcor1"
    shell:
        """
        plink2 --pfile {params.pfile} --r-unphased square yes-really --out results/rss_comparison/ld/{wildcards.region}
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
    sumstat_path = "results/rss_comparison/simulation/{sim_id}/{sim_id}.Continuous.glm.logistic.hybrid",
    ld_path = get_ld_path
  output:
    rss = "results/rss_comparison/simulation/{sim_id}/{sim_id}.logistic_rss_fit.rds"
  script:
    "../scripts/rss_comparison/fit_rss.R"


prefix = "/project2/mstephens/yuxin/ukb-bloodcells/genotypes/bloodcells_chr1.100794065.101983457"
"resources/ukb-geno/bloodcells_chr14.100926769.101431881.pgen"
"snakemake -c1 results/rss_comparison/simulation/chr14.100926769.101431881/Continuous.glm.linear"
"snakemake -c1 results/rss_comparison/simulation/chr14.100926769.101431881/chr14.100926769.101431881.ld"
"plink2 --pfile resources/ukb-geno/bloodcells_ch14.100926769.101431881 --r-unphased square yes-really"
# results/rss_comparison/simulation/chr14.100926769.101431881/f43ak.linear_rss_fit.rds
