data {
    int<lower=0> n_guides; // number of guides
    int<lower=0> n_genes; // number of genes
    int<lower=0> n_samples; // number of samples
    vector[n_samples] sf; // normalization constant/offset
    array[n_guides, n_samples] int y; // observed data
    vector[n_samples] x; // design matrix
    array[n_guides] int guide_to_gene; // guide to gene mapping or group assignments
    vector<lower=0, upper=1>[n_guides] is_ntc; // non-targeting control indicator
    vector[n_guides] mu_g; // log normalized counts for each guide
}

parameters {
    // parameters
    vector[n_genes] mu; // gene-level expression
    vector[n_guides] beta0; // guide-level base expression

    // parameter variances
    real<lower=0> sigma; // gene-level expression standard deviation
    real<lower=0> tau; // global targeting guide standard deviation
    real<lower=0> tau_ntc; // global non-targeting control guide standard deviation
    vector[n_guides] z_beta1; // guide-level expression z-parametrization
    real<lower=0> sigma_phi; // guide-wise dispersion standard deviaion
    vector[n_guides] z_phi; // guide-wise dispersion z-parametrization
    real<lower=0> sigma_beta0; // guide-level base expression standard deviation

    // dispersion parabolic function coefficients
    real<lower=0> a; // guide-wise dispersion coefficient
    real<lower=0> b; // guide-wise dispersion coefficient
}


transformed parameters {
    // guide-level expression
    vector[n_guides] beta1 = mu[guide_to_gene] + tau .* z_beta1;

    // guide-wise dispersoin
    vector[n_guides] log_trend = log(a + b / exp(mu_g));
    vector[n_guides] log_phi = log_trend + sigma_phi .* z_phi;
    vector<lower=0>[n_guides] phi = exp(log_phi);
}

model {
    // guide-level base expression
    beta0 ~ normal(0, sigma_beta0);
    sigma_beta0 ~ normal(0, 1);

    // gene-level expression
    sigma ~ normal(0, 1);
    mu ~ cauchy(0, sigma);

    // guide-level expression
    tau ~ normal(0, 1);
    tau_ntc ~ normal(0, 1);
    z_beta1 ~ normal(0, 1);

    // guide-wise dispersion variance
    sigma_phi ~ normal(0, 1);
    z_phi ~ normal(0, 1);

    // likelihood
    for (n in 1:n_samples) {
        y[,n] ~  neg_binomial_2_log(log(sf[n]) + beta0 + (beta1 .* x[n]), 1 ./ phi);
    }
}
