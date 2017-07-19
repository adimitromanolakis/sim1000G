



# Modified version of haplosim2 function from hapsim package.
# Used internally

haplosim2 = function (n, hap, which.snp = NULL, seed = NULL, force.polym = TRUE,
                      summary = TRUE)
{
    if (length(seed) > 0)
        set.seed(seed)
    nsubset <- length(which.snp)
    if (nsubset > 0) {
        which.snp <- sort(unique(which.snp))
        if ((which.snp[1] < 1) || (which.snp[nsubset] > length(hap$freqs)))
            stop("which.snp does not contain valid indeces")
        indexset <- which.snp
    } else { indexset <- seq(length(hap$freqs)) }

    nloci <- length(indexset)
    quants <- qnorm(hap$freqs[indexset])
    y <- matrix(0, nrow = n, ncol = nloci)
    colnames(y) <- names(hap$freqs[indexset])
    A <- mvrnorm(n, mu = rep(0, nloci), Sigma = hap$cov[indexset,
                                                        indexset], tol = 1e-06, empirical = FALSE)
    if (n == 1)
        A <- array(A, dim = c(1, nloci))
    for (i in 1:n) {
        for (j in 1:nloci) {
            if (A[i, j] <= quants[j])
                y[i, j] <- 0
            else y[i, j] <- 1
        }
    }
    y.freqs <- allelefreqs(y)

    y.freqs <- y.freqs$freqs
    if (summary) {
        y.cor <- cor(y)
        y.div <- divlocus(y)
        mse.freqs <- mse(hap$freqs[indexset], y.freqs)
        hap.cor <- hap$cor[indexset, indexset]
        tmat1 <- hap.cor[upper.tri(hap.cor)]
        tmat2 <- y.cor[upper.tri(y.cor)]
        mse.cor <- mse(tmat1, tmat2)
        return(list(data = y, freqs = y.freqs, cor = y.cor, div = y.div,
                    mse.freqs = mse.freqs, mse.cor = mse.cor))
    }
    else return(list(data = y, freqs = y.freqs))
}




