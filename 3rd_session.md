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
    return(q_prob * q_hat)
}
D2 <- function (q_prob, q_hat) {
    return((1 - q_prob) * q_hat)
}

# Incremental Algorithm (returns quantile estimates at iteration t for each lambda_vec value)
DUMIQE <- function (t, x_obs, q_hat_previous, lambda_vec, q_prob) {
    if (t == 0) { # Initialization step (Q_0)
        Q_min <- 14 # chosen to be "outside" the data stream distribution
        return(rep (Q_min, times = length(lambda_vec))) # L was still not defined
    } else {
        return(ifelse (
            x_obs[t + 1] >= q_hat_previous,                           # Logical vector
            # RECALL: between 2 vectors * and + performs element-wise multiplication and addition, respectively
            q_hat_previous + lambda_vec * D1(q_prob, q_hat_previous), # If TRUE
            q_hat_previous - lambda_vec * D2(q_prob, q_hat_previous)  # If FALSE
        )) 
    }
}
```

Quantile Loss function on iteration t

``` r
QL <- function (x_obs_t, q_hat_t, q_prob) {
    return(ifelse(
        x_obs_t >= q_hat_t,                 # Logical vector
        (1 - q_prob) * (x_obs_t - q_hat_t), # If TRUE
        q_prob * (q_hat_t - x_obs_t)        # If FALSE
    ))
}
```

# Figure 2

Replicating Figure 2 (synthetic) experiment.

## Data stream generation

Change-point at time t = T. The number of observations in each data
stream is $2\text{T}$.

``` r
tau1  <- 500
tau2  <- 10^4
T     <- 10^4
tau_n <- function (t) {
    return(tau1 * as.integer(t < T) + tau2 * as.integer(t >= T & t < 2*T))
}
iter <- 0:(2*T - 1)

# Normally distributed data stream
miu    <- 8
b      <- 2
sigma  <- 1
mean_t <- miu + b * sin(2 * pi * iter / tau_n(iter))
normal_stream <- rnorm(
    n = 2*T,
    mean = mean_t,
    sd = sigma
)
```

``` r
q_fig2 <- 0.7
```

## *Oracle* Approach

``` r
possible_lambdas <- exp(seq (from = -7, to = 0, by = 0.05))
L <- length(possible_lambdas) # number of estimates per iteration 

gamma_Oracle <- 1 - (0.01)^(1/1000) # 0.005
```

``` r
q_hat   <- matrix(0, nrow = 2*T, ncol = L)
QL_hat  <- matrix(0, nrow = 2*T, ncol = L)
EQL_hat <- matrix(0, nrow = 2*T, ncol = L)

oracle_lambdas   <- numeric(2*T)
oracle_quantiles <- numeric(2*T) # Oracle Quantile Estimates
oracle_QL        <- numeric(2*T)
oracle_EQL       <- numeric(2*T)

for (t in iter) {
    # Quantiles estimates on iteration t for each value in possible_lambdas 
    q_hat[t + 1, ] <- DUMIQE(
        t, 
        normal_stream, 
        q_hat[max(t,1), ], # max is due to the initialization step
        possible_lambdas,
        q_fig2
    )

    QL_hat[t + 1, ] <- QL(normal_stream[t + 1], q_hat[t + 1, ], q_fig2)

    # Expected Quantile Loss (EQL) function on iteration t. It is based on the exponentially weighted moving average.
    # It was previously defined as a function, but the node stack overflow error was thrown given it was defined recursively (so, e.g., it was evaluating all the previous steps creating very deep stack call)
    if (t == 0) {
        EQL_hat[t + 1, ] <- QL_hat[t + 1, ] # Initialization
    } else {
        EQL_hat[t + 1, ] <-
            (1 - gamma_Oracle) * EQL_hat[t, ] +
            gamma_Oracle * QL_hat[t + 1, ]
    }

    index <- which.min(EQL_hat[t + 1, ])
    # The Oracle approach stores the lambda value with the minimum estimated EQL
    oracle_lambdas [t + 1]  <- possible_lambdas[index]
    oracle_quantiles[t + 1] <- q_hat[t + 1, index]
    oracle_QL[t + 1]        <- QL_hat[t + 1, index]
    oracle_EQL[t + 1]       <- EQL_hat[t + 1, index]
}
```

## Diagnostics Plots

Tracking the q = 0.7 quantile using the Oracle approach.

<!-- 
In the upper left panel, we see that quantile loss quickly decreases and then stabilizes. The initial value of the quantile tracking algorithm is outside the data stream distribution. Therefore, in the initial phase, the quantile loss is high and a high value of $\lambda$ (update size) is used for the tracking algorithm to quickly get within the data stream distribution. When the quantile tracking algorithm is within the data stream distribution, a smaller $\lambda$ is optimal and the value of $\lambda$ is reduced. The expected quantile loss is reduced accordingly.
 -->

``` r
plot(
    x = iter,
    y = oracle_EQL, # oracle_EQL, oracle_QL
    type = "l", # Line Plot
    ylim = c(0.00,1),
    xlab = "Iteration",
    ylab = "Expected Quantile Loss (EQL)"
)
```

![](3rd_session_files/figure-gfm/Estimated%20EQL%20by%20iteration%20graph-1.png)<!-- -->

How can the theoretically optimal value of $\lambda$, denoted by
$\lambda_{\text{Theoretical}}$, be computed? <!-- 
Optimal (oracle) value of $\lambda$ for the data stream dynamics before and after the concept change at iteration 10 000. We see that a larger $\lambda$ is needed to efficiently track rapid variations. We further see that the Oracle approach uses $\lambda$ close to the optimal theoretical $\lambda$. 
-->

``` r
plot(
    x = iter,
    y = oracle_lambdas,
    type = "l", # Line Plot
    ylim = c(0,0.50),
    xlab = "Iteration",
    ylab = expression(lambda)
)

# ALERT: They're not computed this way!
theoretical_lambdas <- c(rep(mean(oracle_lambdas[1:T]), times = T), rep(mean(oracle_lambdas[(T+1):(2*T)]), times = T))

lines(
    x = iter,
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

![](3rd_session_files/figure-gfm/Recursive%20Lambda%20Updating%20graph-1.png)<!-- -->

``` r
plot(
    x = iter,
    y = normal_stream,
    type = "p", # Line Plot
    pch = 16, cex = 0.3,
    col = "grey",     
    ylim = c(0,14),
    xlab = "Iteration",
    ylab = expression("Quantile at " ~ q == 0.7)
)

lines(
    x = iter,
    y = oracle_quantiles, # Estimated Quantiles...
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

<!-- # Figure 3 -->

<!-- ## _HIL_ Approach -->

# Streaming POT

Recall from Carvalho et al. ([2026](#ref-deCarvalho2026_handbook)) that
the maximum of an i.i.d. sample from a Student’s *t*-distribution with
$\nu > 0$ degrees of freedom belongs to the *max-domain of attraction*
$\mathcal{D}_{1/\nu}$ (Fréchet – Type II).
<!-- O Hill estimator só é usado com Type II distributions. In those cases, $\xi$ is also known as the tail index... -->

``` r
N <- 1000
t_data <- rt(N, df = 1)
# t_data <- rcauchy(N, location = 0, scale = 1)

# There's also the rgpd command from the 'evd' package to generate directly data from a GP model...
# f(n/k) is known for the student-t, right? just use its density function, right?
```

<!-- 
It makes sense to me that the first thing to do is to fit a GP model and find a first level of a suitable threshold, I guess. Only then can I use stochastic approximation (online stochastic gradient descent algorithm) to continuously update $\xi$ and $\sigma_u$, right? I should have an initial estimate for them. Or should I use a rolling window Hill estimator (for $\xi$)?
&#10;Instead of using Rue's method to keep track of the quantiles, one can apply the quantile(df = df, prob = ) function to every entry of the u_track vector to see what's the current sample quantile estimate...
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
