

mahalanobisDistance <- function(x, center, cov) {
  x <- sweep(x, 2L, center)
  cov <- solve(cov)
  setNames(rowSums(x %*% cov * x), rownames(x))
}

multivariateOutliers <- function (X, quantile = 0.975, plot = TRUE) 
  {
  library('robustbase')
    X.mcd <- covMcd(X)
    md = sqrt(mahalanobisDistance(as.matrix(X), apply(X, 2, mean), cov(as.matrix(X))))
    rd = sqrt(mahalanobisDistance(as.matrix(X), X.mcd$center, X.mcd$cov))
    cutoff <- sqrt(qchisq(quantile, ncol(X)))
    if (plot) {
      par(mfrow = c(1, 2))
      plot(1:length(md), md, xlab = "Index of object", ylim = c(0, 
                                                                max(md, rd)), ylab = "Classical Mahalanobis distance")
      abline(h = sqrt(qchisq(quantile, ncol(X))), lty = 2)
      plot(1:length(rd), rd, xlab = "Index of object", ylim = c(0, 
                                                                max(md, rd)), ylab = "Robust Mahalanobis distance")
      abline(h = cutoff, lty = 2)
    }
    list(md = md, rd = rd, cutoff = cutoff)
}

adjustedQuantileOutliers <- function (x, delta = qchisq(0.975, df = ncol(x)), quan = 1/2, alpha = 0.05) 
{
  covr <- covMcd(x, alpha = quan)
  dist <- mahalanobisDistance(x, center = covr$center, cov = covr$cov)
  s <- sort(dist, index = TRUE)
  z <- x
  # par(mfrow = c(1, 2))

  # Plot of cumulative probability of individuals regarding mahalanobis distance
  plot(s$x, (1:length(dist))/length(dist), col = "green", xlab = "Ordered squared robust distance", ylab = "Cumulative probability", type = "n")
  # Introduce index as point
  text(s$x, (1:length(dist))/length(dist), as.character(s$ix), col = "green", cex = 0.8)
  
  # Introduce pink line according chi-square cumualtive probability (0.01 steps)
  t <- seq(0, max(dist), by = 0.01)
  lines(t, pchisq(t, df = ncol(x)), col = "pink")
  
  # Showing vertical line regarding the mahalanobis distance according the introduced quantile (e.g 97.5%)
  abline(v = delta, col = "blue")
  text(x = delta, y = 0.4, paste(100 * (pchisq(delta, df = ncol(x))), "% Quantile", sep = ""), col = "blue", pos = 2, srt = 90, cex = 0.8)
  
  # Showing vertical line regarding the mahalanobis distance according the adjusted quantile (e.g 97.5%, calculated by the 
  # Adaptative rightweighted estimator) for the outlier threshold
  xarw <- arw(x, covr$center, covr$cov, alpha = alpha)
  if (xarw$cn < Inf) {
    abline(v = xarw$cn, col = "cadetblue1")
    text(x = xarw$cn, y = 0.4, "Adjusted threshold", col = "cadetblue1", pos = 4, srt = 90, cex = 0.8)
  }
  
  if (ncol(x) > 2) {
    p <- princomp(x, covmat = covr)
    # Take the value of the score (distance) for the first 2 principal components.
    z <- p$scores[, 1:2]
  }
  
  plot(z, col = "green", type = "n", main = "Outliers (adjusted threshold)")
  if (xarw$cn < Inf) {
    # Showing individuals treated as outliers in the two first principal components. (Disance > adaptive outlier threshold)
    text(z[dist > xarw$cn, 1], z[dist > xarw$cn, 2], dimnames(as.data.frame(x)[dist > xarw$cn, ])[[1]], col = "red", cex = 0.8)
  }
  # Showing individuals in the two first principal components. (Disance <= adaptive outlier threshold)
  text(z[dist <= xarw$cn, 1], z[dist <= xarw$cn, 2], dimnames(as.data.frame(x)[dist <= xarw$cn, ])[[1]], col = "green", cex = 0.8)
  # List of all individuals. The outliers are computed by comparing the mahalanobis distance with the maximum between the adjusted outlier threshold 
  # and the chi square quantile function with 9 degrees of freedom (9 dimensions of the dataset).
  o <- (dist > max(xarw$cn, qchisq(0.975, dim(x)[2])))
  subset(o, o == TRUE)
}