Online POT Model
================
David Baptista

<!--
Other ways of making a link clickable:
[Efficient quantile tracking using an oracle](https://doi.org/10.1007/s10489-022-03489-1)
<https://doi.org/10.1007/s10489-022-03489-1>
&#10;message = FALSE hides only the output.
echo    = FALSE hides only the code but shows the output.
include = FALSE runs the chunk and hides both code and output.
&#10;**This is bold text**
*This is italic text*
***This is bold and italic text***
-->

<!-- comment -->

The main purpose of this document is to replicate both synthetic
experiments presented in Hammer et al. ([2023](#ref-rue2023oracle)). The
second part aims to use the findings from these experiments to
investigate how the methodology can be adapted for the continuous
tracking of extreme quantiles.

# Support Functions

Deterministic Multiplicative Incremental Quantile Estimator (DUMIQE),
based on Yazidi and Hammer ([2019](#ref-DUMIQE_algorithm)). This article
provides a more detailed explanation of the initialization step of the
DUMIQE algorithm. In the experiments presented in Hammer et al.
([2023](#ref-rue2023oracle)), there is no need to account for negative
quantiles or the possibility of the algorithm producing quantile
estimates equal to zero. Therefore, this document presents a simplified
version of the DUMIQE algorithm.

``` r
D1 <- function (q_prob, q_hat) {
    q_prob * q_hat
}
D2 <- function (q_prob, q_hat) {
    (1 - q_prob) * q_hat
}

# Incremental Algorithm (returns quantile estimates at iteration t for each lambda_vec value)
DUMIQE <- function (t, x_obs, q_hat_previous, lambda_vec, q_prob, Q_0) {
    if (t == 0) { # Initialization step (Q_0)
        # Q_0 chosen to be "outside" the data stream distribution
        rep(Q_0, times = length(lambda_vec)) # L was still not defined
    } else {
        ifelse(
            x_obs[t] >= q_hat_previous,                           # Logical vector
            # RECALL: between 2 vectors, * and + performs element-wise multiplication and addition, respectively
            q_hat_previous + lambda_vec * D1(q_prob, q_hat_previous), # If TRUE
            q_hat_previous - lambda_vec * D2(q_prob, q_hat_previous)  # If FALSE
        )
    }
}
```

``` r
# Quantile Loss function on iteration t 
QL <- function (x_obs_t, q_hat_t, q_prob) {
    ifelse(
        x_obs_t >= q_hat_t,                 # Logical vector
        q_prob * (x_obs_t - q_hat_t), # If TRUE
        (1 - q_prob) * (q_hat_t - x_obs_t)        # If FALSE
    )
}
```

# Figure 2

Replicating Figure 2 (synthetic) experiment.

## Data stream generation

Change-point at time t = T. The number of observations in each data
stream is $2\text{T}$. Given the streams are changing at every time
step, the same is true for the underlying quantiles.

``` r
tau1 <- 500
tau2 <- 10^4 # 1e4 also works?
T    <- 10^4
N    <- 2*T

# t %% (2 * T) is always lower than 2 * T
tau_fun <- function(t) {
    ifelse(t %% (2 * T) < T, tau1, tau2)
}

# Normally distributed data stream
miu <- 8; b <- 2; sigma <- 1
# Chi-Square distributed data stream
nu <- 6

generate_data <- function (family, N, q_prob) {
    
    iter <- 0:(N-1)

    if (family == "normal"){
        mean_t    <- miu + b * sin(2 * pi * iter / tau_fun(iter))
        stream    <- rnorm(n = N, mean = mean_t, sd = sigma)
        quantiles <- qnorm(q_prob, mean = mean_t, sd = sigma)
    }  
    if (family == "chisq"){
        df_t      <- nu + b * sin(2 * pi * iter / tau_fun(iter))
        stream    <- rchisq(n = N, df = df_t)
        quantiles <- qchisq(q_prob, df = df_t)
    }
    # Return
    list(
        stream = stream,
        quantiles = quantiles
    )
}
```

``` r
q_fig2 <- 0.7
q_fig3 <- c(0.5, 0.7, 0.9)
```

``` r
normal_stream_fig2 <- generate_data("normal", N, q_fig2)$stream
```

## *Oracle* Approach

``` r
possible_lambdas <- exp(seq (from = -7, to = 0, by = 0.05))
L <- length(possible_lambdas) # number of estimates per iteration 
gamma_fig2 <- 1 - (0.01)^(1/1000) # 0.005
```

``` r
oracle_fun <- function (N, data_stream, gamma, q_prob, lambdas_list = possible_lambdas, Q_initial) {

    q_hat   <- matrix(0, nrow = N, ncol = L)
    QL_hat  <- matrix(0, nrow = N, ncol = L)
    EQL_hat <- matrix(0, nrow = N, ncol = L)

    oracle_quantiles <- numeric(N) # Oracle Quantile Estimates
    oracle_QL        <- numeric(N)
    oracle_EQL       <- numeric(N)
    oracle_lambdas   <- numeric(N)

    for (t in 0:(N-1)) {

        # Quantiles estimates on iteration t for each value in lambdas_list 
        q_hat[t + 1, ] <- DUMIQE(
            t, 
            data_stream, 
            q_hat[max(t,1), ], # max is due to the initialization step
            lambdas_list,
            q_prob,
            Q_0 = Q_initial
        )

        QL_hat[t + 1, ] <- QL(data_stream[t + 1], q_hat[t + 1, ], q_prob)

        # Expected Quantile Loss (EQL) function on iteration t. It is based on the exponentially weighted moving average.
        # It was previously defined as a function, but the node stack overflow error was thrown given it was defined recursively (so, e.g., it was evaluating all the previous steps creating very deep stack call)
        if (t == 0L) {
            EQL_hat[t + 1, ] <- QL_hat[t + 1, ] # Initialization
        } else {
            EQL_hat[t + 1, ] <- (1 - gamma) * EQL_hat[t, ] + gamma * QL_hat[t + 1, ]
        }

        index <- which.min(EQL_hat[t + 1, ])
        # The Oracle approach stores the lambda value with the minimum estimated EQL
        oracle_lambdas [t + 1]  <- lambdas_list[index]
        oracle_quantiles[t + 1] <- q_hat[t + 1, index]
        oracle_QL[t + 1]        <- QL_hat[t + 1, index]
        oracle_EQL[t + 1]       <- EQL_hat[t + 1, index]
    }

    # Return a list
    list(
        lambdas   = oracle_lambdas,
        quantiles = oracle_quantiles,
        QL        = oracle_QL,
        EQL       = oracle_EQL
    )   
}
```

``` r
oracle_fig2 <- oracle_fun(2*T, normal_stream_fig2, gamma_fig2, q_fig2, Q_initial = 14)
```

## Diagnostics Plots

Tracking the q = 0.7 quantile using the Oracle approach.

<!-- 
In the upper left panel, we see that quantile loss quickly decreases and then stabilizes. The initial value of the quantile tracking algorithm is outside the data stream distribution. Therefore, in the initial phase, the quantile loss is high and a high value of $\lambda$ (update size) is used for the tracking algorithm to quickly get within the data stream distribution. When the quantile tracking algorithm is within the data stream distribution, a smaller $\lambda$ is optimal and the value of $\lambda$ is reduced. The expected quantile loss is reduced accordingly.
 -->

``` r
par(mfrow = c(1,1))
#| label: Estimated EQL by iteration graph
plot(
    x = 0:(N-1),
    y = oracle_fig2$EQL, # oracle_fig2$EQL, oracle_fig2$QL
    type = "l", # Line Plot
    ylim = c(0.00,1),
    xlab = "Iteration",
    ylab = expression("Estimated Expected Quantile Loss (" * hat(EQL) * ")")
)
```

![](3rd_session_files/figure-gfm/Estimated%20EQL%20by%20iteration%20graph-1.png)<!-- -->

How can the theoretically optimal value of $\lambda$, denoted by
$\lambda_{\text{Theoretical}}$, be computed? <!-- 
Optimal (oracle) value of $\lambda$ for the data stream dynamics before and after the concept change at iteration 10 000. We see that a larger $\lambda$ is needed to efficiently track rapid variations. We further see that the Oracle approach uses $\lambda$ close to the optimal theoretical $\lambda$. 
-->

``` r
# ALERT: They're not computed this way!
theoretical_lambdas <- c(rep(mean(oracle_fig2$lambdas[1:T]), times = T), rep(mean(oracle_fig2$lambdas[(T+1):(2*T)]), times = T))

plot(
    x = 0:(N-1),
    y = oracle_fig2$lambdas,
    type = "l", # Line Plot
    ylim = c(0,0.10),
    xlab = "Iteration",
    ylab = expression(lambda)
)

lines(
    x = 0:(N-1),
    y = theoretical_lambdas,
    col = "green"
)

legend(
    "topright",
    legend = c(expression(lambda[Oracle]), expression(lambda[Theoretical])),
    col = c("black", "green"),
    lty = c(1, 1)
)
```

![](3rd_session_files/figure-gfm/Recursive%20Lambda%20Updating%20Graph-1.png)<!-- -->

``` r
plot(
    x = 0:(N-1),
    y = normal_stream_fig2,
    type = "p", # Line Plot
    pch = 16, cex = 0.3,
    col = "grey",     
    ylim = c(0,14),
    xlab = "Iteration",
    ylab = expression("Quantile at " ~ q == 0.7)
)

lines(
    x = 0:(N-1),
    y = oracle_fig2$quantiles, # Estimated Quantiles...
    col = "black",
    lwd = 1
)

legend(
    "topright",
    legend = c(expression(hat(Q)[t * "," * q](lambda[t])), expression(x[t])),
    col = c("black", "grey"),
    lty = c(1, NA),
    pch = c(NA, 16)
)
```

![](3rd_session_files/figure-gfm/Quantile%20Tracking%20Plot-1.png)<!-- -->

# Streaming POT

Recall from Carvalho et al. ([2026](#ref-deCarvalho2026_handbook)) that
the maximum of an i.i.d. sample from a Student’s *t*-distribution with
$\nu > 0$ degrees of freedom belongs to the *max-domain of attraction*
$\mathcal{D}_{1/\nu}$ (Fréchet – Type II).

<!-- 
O Hill estimator só é usado com Type II distributions. In those cases, $\xi$ is also known as the tail index... 
# f(n/k) is known for the student-t, right? just use its density function, right? De Haan & Ferreira Book
# t_data <- rt(N, df = 1)
# t_data <- rcauchy(N, location = 0, scale = 1)
-->

``` r
# 5 years of observations
N    <- 5 * 365 
time <- 0:(N-1) 

# (Hyper-) Parameters to be fitted...
beta0  <- 1; beta1  <- 1; beta2 <- 1
alpha0 <- 1; alpha1 <- 0.5

# Underlying (distribution) parameters 
mu_gp    <- beta0 + beta1 * cos((2 * pi * time / 365) - beta2)
sigma_gp <- exp (alpha0 + alpha1 * time)
shape_gp <- 1

library(evd)
gp_data <- rgpd(N, loc = mu_gp, scale = sigma_gp, shape = shape_gp)
```

``` r
# It would be interesting to, at least, "find" some possible values for the threshold...
```

Several quantiles are tracked using the methodology proposed by Hammer
et al. ([2023](#ref-rue2023oracle)), providing a set of candidate
thresholds for subsequent analysis of the exceedances.

``` r
u_prob_track <- c(0.7) #, 
#0.75, 
#0.8, 
#0.85,
# 0.9) 
```

``` r
possible_lambdas <- exp(seq (from = -7, to = 0, by = 0.05))
L <- length(possible_lambdas) # number of estimates per iteration 
gamma <- 1 - (0.01)^(1/1000) # close to 0.005
```

When the underlying process is non-stationary, what is the value of
tracking a threshold without continuously assessing its suitability for
accurately modelling the exceedances in an online setting?

``` r
oracle_thresholds <- vector("list", length(u_prob_track))

for (i in seq_along(u_prob_track)) {
  oracle_thresholds[[i]] <- oracle_fun(N, gp_data, gamma, u_prob_track[i], Q_initial = 20)$quantiles
}
```

``` r
y_u <- vector("list", length(u_prob_track))

for (i in seq_along(u_prob_track)) {

    y_u[[i]] <- vector("list", length(time)) # List of lists

    for (j in time) {
        
        y_u[[i]][[j + 1]] <- ifelse(
            gp_data > oracle_thresholds[[i]][[j + 1]],
            gp_data - oracle_thresholds[[i]][[j + 1]],
            0
        ) 

    }

}
```

I could have used instead functions from the ‘ismev’ R package.

``` r
model_MLEs <- vector("list", length(u_prob_track))
y_u_hat    <- vector("list", length(u_prob_track))

for (i in seq_along(u_prob_track)) {

    u_vec <- oracle_thresholds[[i]]
    model_MLEs[[i]] <- vector("list", length(time))
    y_u_hat[[i]]    <- vector("list", length(time))

    for (j in time) { # time starts at j = 0;
        
        fit <- fpot(
            gp_data[1:max(365,(j + 1))], # max: to start with a decent n_obs
            threshold = u_vec[max(365,(j + 1))], # for stability purposes
            model = "gpd", 
            # "pp": non-homogeneous Point (Poisson) Process
            # "gpd": non-homogeneous Poisson Process
            cmax = FALSE,
            std.err = FALSE
        )
        model_MLEs[[i]][[j + 1]] <- fit$estimate
        y_u_hat[[i]][[j + 1]]    <- fitted(fit)
    }
}
```

<!-- 
It makes sense to me that the first thing to do is to fit a GP model and find a first level of a suitable threshold, I guess. Only then can I use stochastic approximation (online stochastic gradient descent algorithm) to continuously update $\xi$ and $\sigma_u$, right? I should have an initial estimate for them. Or should I use a rolling window Hill estimator (for $\xi$)?
-->

<!-- 
See what is the numerical algorithm behind the 'gpd.fit' function from the 'ismev' R package.
-->

<!-- 
gpd.fitrange and mrl.plot commands
-->

<!-- 
For stationary models four plots are produced; a probability plot, a quantile plot, a return level plot and a histogram of data with fitted density. For non-stationary models two plots are produced; a residual probability plot and a residual quantile plot.
-->

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-deCarvalho2026_handbook" class="csl-entry">

Carvalho, Miguel de, Raphaël Huser, Philippe Naveau, and Brian J. Reich.
2026. *Handbook of Statistics of Extremes*. Chapman & Hall/CRC.

</div>

<div id="ref-rue2023oracle" class="csl-entry">

Hammer, Hugo L., Anis Yazidi, Michael A. Riegler, and Håvard Rue. 2023.
“Efficient Quantile Tracking Using an Oracle.” *Applied Intelligence*
53: 289–300. <https://doi.org/10.1007/s10489-022-03489-1>.

</div>

<div id="ref-DUMIQE_algorithm" class="csl-entry">

Yazidi, Anis, and Hugo Lewi Hammer. 2019. “Multiplicative Update Methods
for Incremental Quantile Estimation.” *IEEE Transactions on Cybernetics*
49 (3): 746–56. <https://doi.org/10.1109/TCYB.2017.2779140>.

</div>

</div>
