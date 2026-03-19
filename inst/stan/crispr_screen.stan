data {
    int<lower=0> n_guides;                          // number of guides
    int<lower=0> n_genes;                           // number of genes
    int<lower=0> n_samples;                         // number of samples
    vector[n_samples] sf;                           // size factors (log-scale offset)
    array[n_guides, n_samples] int y;               // observed counts
    vector[n_samples] x;                            // design vector (0 = control, 1 = treatment)
    array[n_guides] int guide_to_gene;              // guide-to-gene mapping
    vector<lower=0, upper=1>[n_guides] is_ntc;      // non-targeting control indicator
    vector[n_guides] beta0_hat;                     // log-normalized control counts per guide
}

parameters {
    // base expression
    real<lower=0> sigma_beta0;                      // guide-level base expression SD
    vector[n_guides] beta0;                         // guide-level base expression

    // gene-level treatment effects (horseshoe)
    real<lower=0> sigma;                            // global shrinkage
    vector<lower=0>[n_genes] lambda;                // local shrinkage
    vector[n_genes] z_mu;                           // gene-level effects (non-centered)

    // guide-level treatment effects
    real mu_ntc;                                    // NTC guide mean effect
    real<lower=0> tau;                              // targeting guide effect SD
    real<lower=0> tau_ntc;                          // NTC guide effect SD
    vector[n_guides] z_beta1;                       // guide-level effects (non-centered)

    // guide efficiency
    vector[n_guides] z_eff;                         // guide efficiency (non-centered)
    real<lower=0> sigma_eff;                        // guide efficiency SD

    // dispersion
    real<lower=0> a;                                // dispersion trend intercept
    real<lower=0> b;                                // dispersion trend coefficient
    real<lower=0> sigma_phi;                        // dispersion SD around trend
    vector[n_guides] z_phi;                         // dispersion (non-centered)
}

transformed parameters {
    // guide efficiency, constrained so mean within each gene = 1
    vector[n_guides] eff;
    {
        vector[n_guides] eff_raw = exp(sigma_eff * z_eff);
        vector[n_genes] sum_eff   = rep_vector(0, n_genes);
        vector[n_genes] count_eff = rep_vector(0, n_genes);

        for (g in 1:n_guides) {
            if (is_ntc[g] == 0) {
                sum_eff[guide_to_gene[g]]   += eff_raw[g];
                count_eff[guide_to_gene[g]] += 1;
            }
        }

        vector[n_genes] mean_eff = sum_eff ./ count_eff;

        for (g in 1:n_guides) {
            eff[g] = is_ntc[g] == 0
                ? eff_raw[g] / mean_eff[guide_to_gene[g]]
                : 1.0;
        }
    }

    // gene-level effects
    vector[n_genes] mu = (sigma * lambda) .* z_mu;

    // guide-level effects
    vector[n_guides] beta1;
    for (g in 1:n_guides) {
        beta1[g] = is_ntc[g] == 1
            ? mu_ntc + tau_ntc * z_beta1[g]
            : mu[guide_to_gene[g]] * eff[g] + tau * z_beta1[g];
    }

    // dispersion: trend + deviation (DESeq2-style)
    vector[n_guides] log_trend = log(a + b * exp(-beta0_hat));
    vector<lower=0>[n_guides] alpha = exp(log_trend + sigma_phi * z_phi);
    vector<lower=0>[n_guides] phi = inv(alpha);
}

model {
    // base expression
    beta0 ~ normal(mean(beta0_hat), sigma_beta0);
    sigma_beta0 ~ normal(0, 1);

    // gene-level effects (horseshoe)
    sigma ~ normal(0, 1);
    lambda ~ normal(0, 1);
    z_mu ~ normal(0, 1);

    // guide-level effects
    mu_ntc ~ normal(0, 1);
    tau  ~ normal(0, 1);
    tau_ntc ~ normal(0, 1);
    z_beta1 ~ normal(0, 1);

    // guide efficiency
    sigma_eff ~ normal(0, 1);
    z_eff ~ normal(0, 1);

    // dispersion
    sigma_phi ~ normal(0, 1);
    z_phi ~ normal(0, 1);

    // likelihood
    for (n in 1:n_samples) {
        y[, n] ~ neg_binomial_2_log(log(sf[n]) + beta0 + beta1 * x[n], phi);
    }
}
