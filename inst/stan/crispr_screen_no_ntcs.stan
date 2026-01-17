data {
    int<lower=0> n_guides; // number of guides
    int<lower=0> n_genes; // number of genes
    int<lower=0> n_samples; // number of samples
    vector[n_samples] sf; // normalization constant/offset
    array[n_guides, n_samples] int y; // observed data
    vector[n_samples] x; // design matrix (0 for control, 1 for treatment)
    array[n_guides] int guide_to_gene; // guide to gene mapping or group assignments
    vector<lower=0, upper=1>[n_guides] is_ntc; // non-targeting control indicator
    vector[n_guides] beta0_hat; // log-normalized counts for each guide in control samples (if n control samples = 1, use all samples)
}

parameters {
    vector[n_guides] beta0; // guide-level base expression
    vector[n_genes] z_mu; // gene-level expression z-parametrization
    real<lower=0> sigma; // gene-level expression global shrinkage
    vector<lower=0>[n_genes] lambda; // gene-level expression local shrinkage
    real<lower=0> tau; // global targeting guide standard deviation
    vector[n_guides] z_beta1; // guide-level expression z-parametrization
    real<lower=0> sigma_phi; // guide-wise dispersion standard deviation
    vector[n_guides] z_phi; // guide-wise dispersion z-parametrization
    real<lower=0> sigma_beta0; // guide-level base expression
    real<lower=0> a; // guide-wise dispersion coefficient
    real<lower=0> b; // guide-wise dispersion coefficient
    vector[n_guides] z_eff;
    real<lower=0> sigma_eff;
}


transformed parameters {
    // guide efficiency / weights
    vector[n_guides] eff;
    { // local variables
        vector[n_guides] eff_raw = exp(sigma_eff .* z_eff);
        vector[n_genes] sum_eff = rep_vector(0, n_genes);
        vector[n_genes] count_eff = rep_vector(0, n_genes);

        for (g in 1:n_guides) {
            sum_eff[guide_to_gene[g]] += eff_raw[g];
            count_eff[guide_to_gene[g]] += 1;
        }

        // constrain such that mean of eff within each gene = 1
        vector[n_genes] mean_eff = sum_eff ./ count_eff;
        // compute eff
        for (g in 1:n_guides) {
            eff[g] = eff_raw[g] / mean_eff[guide_to_gene[g]];
        }
    }

    vector[n_genes] mu = (sigma .* lambda) .* z_mu; // gene-level expression
    vector[n_guides] beta1; // guide-level expression
    for (g in 1:n_guides) {
        beta1[g] = (mu[guide_to_gene[g]] * eff[g]) + tau .* z_beta1[g];
    }

    vector[n_guides] log_trend = log(a + b / exp(beta0_hat)); // guide-wise dispersion trend function (ref. DESeq2)
    vector<lower=0>[n_guides] alpha = exp(log_trend + sigma_phi .* z_phi); // variation around trend function (ref. DESeq2)
    vector<lower=0>[n_guides] phi = inv(alpha);
}

model {
    z_eff ~ normal(0, 1);
    sigma_eff ~ normal(0, 1);

    // guide-level base expression
    beta0 ~ normal(mean(beta0_hat), sigma_beta0);
    sigma_beta0 ~ normal(0, 1);

    // gene-level expression
    sigma ~ normal(0, 1);
    lambda ~ normal(0, 1);
    z_mu ~ normal(0, 1);

    // guide-level expression
    tau ~ normal(0, 1);
    z_beta1 ~ normal(0, 1);

    // guide-wise dispersion standard deviation
    sigma_phi ~ normal(0, 1);
    z_phi ~ normal(0, 1);

    // likelihood
    for (n in 1:n_samples) {
        y[,n] ~  neg_binomial_2_log(log(sf[n]) + beta0 + (beta1 .* x[n]), phi);
    }
}
