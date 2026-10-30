normalizeQuantiles <-
function (A, ties = TRUE) 
{
    n <- dim(A)
    if (is.null(n)) 
        return(A)
    if (n[2] == 1) 
        return(A)
    O <- S <- array(, n)
    nobs <- rep(n[1], n[2])
    i <- (0:(n[1] - 1))/(n[1] - 1)
    for (j in 1:n[2]) {
        Si <- sort(A[, j], method = "quick", index.return = TRUE)
        nobsj <- length(Si$x)
        if (nobsj < n[1]) {
            nobs[j] <- nobsj
            isna <- is.na(A[, j])
            S[, j] <- approx((0:(nobsj - 1))/(nobsj - 1), Si$x, 
                i, ties = list("ordered", mean))$y
            O[!isna, j] <- ((1:n[1])[!isna])[Si$ix]
        }
        else {
            S[, j] <- Si$x
            O[, j] <- Si$ix
        }
    }
    m <- rowMeans(S)
    for (j in 1:n[2]) {
        if (ties) 
            r <- rank(A[, j])
        if (nobs[j] < n[1]) {
            isna <- is.na(A[, j])
            if (ties) 
                A[!isna, j] <- approx(i, m, (r[!isna] - 1)/(nobs[j] - 
                  1), ties = list("ordered", mean))$y
            else A[O[!isna, j], j] <- approx(i, m, (0:(nobs[j] - 
                1))/(nobs[j] - 1), ties = list("ordered", mean))$y
        }
        else {
            if (ties) 
                A[, j] <- approx(i, m, (r - 1)/(n[1] - 1), ties = list("ordered", 
                  mean))$y
            else A[O[, j], j] <- m
        }
    }
    A
}
