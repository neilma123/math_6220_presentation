
```{r}
#| context: setup

library(magick)

# ===================================================================
# fBM SIMULATION METHODS IMPLEMENTATION
# ===================================================================
# This section implements all 6 methods for simulating Fractional Brownian Motion:
# 
# 1. HOSKING ALGORITHM: Sequential conditioning on past values
#    - Time: O(N²), Space: O(N²)
#    - Best for: Small samples (N < 5,000)
#    - Exact method capturing full covariance
#
# 2. CHOLESKY DECOMPOSITION: Direct matrix factorization
#    - Time: O(N³), Space: O(N²)
#    - Best for: Small samples (N < 1,000)
#    - Simplest and most intuitive exact method
#
# 3. FFT (CIRCULANT EMBEDDING): Circulant matrix via FFT
#    - Time: O(N log N), Space: O(N)
#    - Best for: Large samples (N > 5,000)
#    - Industry standard; near-exact with boundary corrections
#
# 4. CORRELATED RANDOM WALK: Aggregation-based approximation
#    - Time: O(N^{2-2H}), Space: O(N)
#    - Best for: Theory, error control, streaming
#    - Approximates covariance asymptotically
#
# 5. RANDOM MIDPOINT DISPLACEMENT: Bisection with local conditioning
#    - Time: O(N), Space: O(1) worst-case
#    - Best for: Streaming/unknown horizon
#    - Exact conditioning on local neighborhood
#
# 6. SPECTRAL METHOD: Frequency-domain approximation via FFT
#    - Time: O(N log N), Space: O(N)
#    - Best for: Very large N, speed critical
#    - ~4x faster than FFT with Paxson's spectral density approx
# ===================================================================


generate_fbm_hosking <- function(n, H) {
  # Compute autocovariance
  gamma <- function(k) {
    if (k == 0) return(1)
    0.5 * (abs(k + 1)^(2*H) - 2*abs(k)^(2*H) + abs(k - 1)^(2*H))
  }
  
  Z <- numeric(n)
  Z[1] <- rnorm(1)
  
  if (n == 1) return(cumsum(Z))
  
  # Build covariance matrix incrementally
  cov_matrix <- matrix(0, n, n)
  cov_matrix[1, 1] <- 1
  
  for (i in 1:(n-1)) {
    # Build i x i covariance matrix
    for (j in 1:i) {
      for (k in 1:i) {
        cov_matrix[j, k] <- gamma(abs(j - k))
      }
    }
    
    # Conditional mean and variance
    c_vec <- sapply(1:i, function(k) gamma(k + 1))
    
    cov_inv <- tryCatch(
      solve(cov_matrix[1:i, 1:i]),
      error = function(e) solve(cov_matrix[1:i, 1:i] + diag(1e-8, i))
    )
    
    mu <- as.numeric(crossprod(c_vec, cov_inv %*% Z[1:i]))
    sigma2 <- 1 - as.numeric(crossprod(c_vec, cov_inv %*% c_vec))
    sigma2 <- max(sigma2, 1e-8)
    
    Z[i + 1] <- mu + sqrt(sigma2) * rnorm(1)
  }
  
  cumsum(Z)
}

# ===================================================================
# METHOD 2: CHOLESKY DECOMPOSITION
# ===================================================================
generate_fbm_cholesky <- function(n, H) {
  # Create covariance matrix
  cov_matrix <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      cov_matrix[i, j] <- 0.5 * (abs(i - j + 1)^(2*H) - 
                                  2*abs(i - j)^(2*H) + 
                                  abs(i - j - 1)^(2*H))
    }
  }
  
  # Cholesky decomposition
  L <- tryCatch(
    chol(cov_matrix),
    error = function(e) {
      chol(cov_matrix + diag(1e-6, n))
    }
  )
  
  # Generate standard normal increments
  Z <- rnorm(n)
  
  # Generate FBM
  fbm <- crossprod(L, Z)
  cumsum(as.vector(fbm))
}

# ===================================================================
# METHOD 3: FFT (CIRCULANT EMBEDDING METHOD)
# ===================================================================
generate_fbm_fft <- function(n, H) {
  gamma <- function(k) {
    if (k == 0) return(1)
    0.5 * (abs(k + 1)^(2*H) - 2*abs(k)^(2*H) + abs(k - 1)^(2*H))
  }
  
  # Extend to power of 2
  M <- 2^ceiling(log2(2*n))
  
  # Build circulant covariance
  circ_cov <- c(sapply(0:(n-1), gamma), rep(0, M - n))
  circ_cov[(M-n+1):M] <- sapply((n-1):1, gamma)
  
  # FFT of first row
  lambda <- Re(fft(circ_cov))
  
  # Ensure positive (adjust if needed)
  if (any(lambda < -1e-10)) {
    M <- 2 * M
    circ_cov <- c(sapply(0:(n-1), gamma), rep(0, M - n))
    circ_cov[(M-n+1):M] <- sapply((n-1):1, gamma)
    lambda <- Re(fft(circ_cov))
  }
  
  # Take square root
  sqrt_lambda <- sqrt(pmax(lambda, 0))
  
  # Generate samples
  U <- rnorm(M)
  U_fft <- fft(U)
  Z_fft <- sqrt_lambda * U_fft
  Z <- Re(fft(Z_fft, inverse = TRUE)) / M
  
  cumsum(Z[1:n])
}

# ===================================================================
# METHOD 4: CORRELATED RANDOM WALK (Approximate)
# ===================================================================
generate_fbm_crw <- function(n, H, M = NULL) {
  if (is.null(M)) {
    # Heuristic for aggregation level
    M <- max(2, floor(sqrt(n / (1 - H))))
  }
  
  # Generate M correlated sequences
  X <- matrix(0, n, M)
  
  for (j in 1:M) {
    # Generate uniform random variables
    U <- runif(n)
    
    # Transform to get correlated increments
    mu_vec <- (1 - U)^(1 / (2 - 2*H))
    
    # Apply sign randomness
    signs <- sample(c(-1, 1), n, replace = TRUE)
    X[, j] <- cumsum(signs * mu_vec)
  }
  
  # Average and scale
  c_H <- sqrt(pi * H / sin(pi * H)) * 2^H
  fbm_approx <- c_H * rowMeans(X)
  
  fbm_approx
}

# ===================================================================
# METHOD 5: RANDOM MIDPOINT DISPLACEMENT (RMD)
# ===================================================================
generate_fbm_rmd <- function(n, H, m = 10, n_neighbors = 10) {
  gamma <- function(k) {
    if (k == 0) return(1)
    0.5 * (abs(k + 1)^(2*H) - 2*abs(k)^(2*H) + abs(k - 1)^(2*H))
  }
  
  # Bisection-based approach
  result <- numeric(n)
  result[1] <- 0
  
  if (n == 1) return(result)
  
  result[n] <- rnorm(1)
  
  # Recursive bisection
  bisect_fill <- function(left, right, result) {
    if (right - left < 2) return(result)
    
    mid <- (left + right) %/% 2
    
    # Determine neighbors
    n_left <- min(n_neighbors, mid - left)
    n_right <- min(n_neighbors, right - mid)
    
    left_idx <- max(left, mid - n_left)
    right_idx <- min(right, mid + n_right)
    
    neighbor_idx <- c(left_idx:(mid-1), (mid+1):right_idx)
    neighbor_vals <- result[neighbor_idx]
    
    if (length(neighbor_idx) > 0) {
      neighbor_distances <- abs(neighbor_idx - mid)
      
      # Build covariance
      cov_mat <- matrix(0, length(neighbor_idx) + 1, length(neighbor_idx) + 1)
      for (i in 1:(length(neighbor_idx)+1)) {
        for (j in 1:(length(neighbor_idx)+1)) {
          if (i == length(neighbor_idx) + 1 & j == length(neighbor_idx) + 1) {
            cov_mat[i, j] <- 1
          } else if (i == length(neighbor_idx) + 1 | j == length(neighbor_idx) + 1) {
            idx <- if (i == length(neighbor_idx) + 1) j else i
            cov_mat[i, j] <- gamma(abs(neighbor_idx[idx] - mid))
          } else {
            cov_mat[i, j] <- gamma(abs(neighbor_idx[i] - neighbor_idx[j]))
          }
        }
      }
      
      # Condition on neighbors
      cov_inv <- tryCatch(
        solve(cov_mat[-nrow(cov_mat), -ncol(cov_mat)]),
        error = function(e) diag(1e-8, nrow(cov_mat) - 1)
      )
      
      mu <- as.numeric(cov_mat[nrow(cov_mat), -ncol(cov_mat)] %*% cov_inv %*% neighbor_vals)
      sigma2 <- 1 - as.numeric(cov_mat[nrow(cov_mat), -ncol(cov_mat)] %*% cov_inv %*% 
                               cov_mat[-nrow(cov_mat), ncol(cov_mat)])
      sigma2 <- max(sigma2, 1e-8)
      
      result[mid] <- mu + sqrt(sigma2) * rnorm(1)
    } else {
      result[mid] <- rnorm(1)
    }
    
    result <- bisect_fill(left, mid, result)
    result <- bisect_fill(mid, right, result)
    result
  }
  
  result <- bisect_fill(1, n, result)
  cumsum(result)
}

# ===================================================================
# METHOD 6: SPECTRAL METHOD (Fast FFT approximation)
# ===================================================================
generate_fbm_spectral <- function(n, H) {
  # Extended to power of 2
  N <- 2^ceiling(log2(n))
  
  # Spectral density approximation for fGn
  spectral_density <- function(lambda, H) {
    if (lambda == 0) return(1)
    # Paxson's approximation
    term1 <- 2 * sin(pi * H) * gamma(2*H + 1) * (1 - cos(lambda)) * abs(lambda)^(-2*H - 1)
    term2 <- abs(lambda)^(1 - 2*H) * (1 - cos(lambda))
    max(term1 + term2, 0)
  }
  
  # Compute spectral density at Fourier frequencies
  freqs <- 2 * pi * (0:(N-1)) / N
  spec <- sapply(freqs, function(f) spectral_density(f, H))
  spec <- pmax(spec, 0)
  
  # Generate coefficients
  a <- complex(N)
  for (k in 1:(N-1)) {
    a[k] <- sqrt(spec[k] / N) * (rnorm(1) + 1i * rnorm(1)) / sqrt(2)
  }
  a[1] <- sqrt(spec[1] / N) * rnorm(1)
  
  # Enforce symmetry for real output
  for (k in 2:(N/2)) {
    a[N - k + 2] <- Conj(a[k])
  }
  
  # Inverse FFT
  Z <- Re(fft(a, inverse = TRUE)) / sqrt(N/2)
  
  cumsum(Z[1:n])
}

# Wrapper function to select method
generate_fbm <- function(n, H, method = "cholesky") {
  method <- tolower(method)
  
  switch(method,
    "hosking" = generate_fbm_hosking(n, H),
    "cholesky" = generate_fbm_cholesky(n, H),
    "fft" = generate_fbm_fft(n, H),
    "crw" = generate_fbm_crw(n, H),
    "rmd" = generate_fbm_rmd(n, H),
    "spectral" = generate_fbm_spectral(n, H),
    generate_fbm_cholesky(n, H)  # Default to Cholesky
  )
}

# Helper to compute frame indices used by both animations
calc_frame_indices <- function(n_steps) {
  frame_skip <- max(1, floor(n_steps / 100))
  sort(unique(c(seq(1, n_steps, by = frame_skip), n_steps)))
}

# Function to create animation frames
create_animation_frames <- function(path, H, n_steps, method = "cholesky", frame_indices = NULL) {
  frames <- list()
  
  if (is.null(frame_indices)) {
    frame_indices <- calc_frame_indices(n_steps)
  } else {
    frame_indices <- sort(unique(frame_indices))
  }
  
  for (step in frame_indices) {
    temp_file <- tempfile(fileext = ".png")
    
    png(temp_file, width = 700, height = 500, res = 100)
    
    plot(
      1:n_steps,
      path,
      type = "n",
      lwd = 2.5,
      col = "#2E86AB",
      xlab = "Time",
      ylab = "Position",
      main = paste0("fBM Simulation: H = ", round(H, 2), 
                    " | ", method, " | Step: ", step, "/", n_steps),
      cex.main = 1.5,
      cex.lab = 1.2
    )
    
    grid(nx = 20, ny = 20, col = "lightgray", lty = 2)
    
    if (step > 0) {
      lines(1:step, path[1:step], lwd = 2.5, col = "#2E86AB")
      points(step, path[step], col = "#E74C3C", cex = 1.5, pch = 16)
    }
    
    if (H > 0.55) {
      interpretation <- "Persistent (Trending)"
      text_col <- "#A23B72"
    } else if (H < 0.45) {
      interpretation <- "Mean-Reverting"
      text_col <- "#F18F01"
    } else {
      interpretation <- "Standard BM"
      text_col <- "#C73E1D"
    }
    
    mtext(interpretation, side = 3, line = -1.5, col = text_col, font = 2, cex = 1.2)
    
    dev.off()
    frames[[length(frames) + 1]] <- image_read(temp_file)
    unlink(temp_file)
  }
  
  animation <- do.call(c, frames)
  animation <- image_animate(animation, fps = 10, loop = 1)
  
  return(animation)
}

estimate_zero_crossings <- function(path) {
  dx <- diff(path)
  if (length(dx) < 2) return(NA_real_)
  dx <- dx - mean(dx)
  signed <- sign(dx)
  signed[signed == 0] <- NA
  crossings <- sum(signed[-1] * signed[-length(signed)] < 0, na.rm = TRUE)
  denom <- max(1, length(dx) - 1)
  p_star <- crossings / denom
  p_star <- pmin(pmax(p_star, 1e-6), 1 - 1e-6)
  1 + log2(cos(pi * p_star / 2))
}

estimate_quadratic_variation <- function(path, max_lag = 10) {
  n <- length(path)
  if (n < 5) return(NA_real_)
  max_h <- min(max_lag, floor((n - 1) / 2))
  if (max_h < 2) return(NA_real_)
  h_vals <- seq_len(max_h)
  V <- sapply(h_vals, function(h) {
    diffs <- path[(h + 1):n] - path[1:(n - h)]
    mean(diffs^2)
  })
  positive_v <- V > 0
  if (!any(positive_v)) return(NA_real_)
  fit <- lm(log(V[positive_v]) ~ log(h_vals[positive_v]))
  slope <- coef(fit)[2]
  as.numeric(slope / 2)
}

estimate_rs <- function(path) {
  n <- length(path)
  if (n < 2) return(NA_real_)
  mean_val <- mean(path)
  Y <- cumsum(path - mean_val)
  R <- max(Y) - min(Y)
  S <- sd(path)
  if (S <= 0 || R <= 0) return(NA_real_)
  log(R / S) / log(n)
}

estimate_whittle <- function(path,
                             H_grid = seq(0.1, 0.9, by = 0.02),
                             m_frac = 0.1) {
  dx <- diff(path)
  N <- length(dx)
  if (N < 64) return(NA_real_)
  dx <- dx - mean(dx)
  fft_vals <- fft(dx)
  periodogram <- (Mod(fft_vals)^2) / N
  max_j <- floor(N / 2)
  m_raw <- floor(N^0.6)
  m <- max(5, floor(m_frac * max_j), m_raw)
  m <- min(m, max_j)
  if (m < 5) return(NA_real_)
  j <- seq_len(m)
  lambda <- 2 * pi * j / N
  I_vals <- Re(periodogram[j + 1])
  scores <- sapply(H_grid, function(H) {
    shape <- lambda^(-(2 * H - 1))
    shape <- pmax(shape, 1e-8)
    sigma2_hat <- mean(I_vals / shape)
    sigma2_hat <- max(sigma2_hat, 1e-8)
    f_hat <- sigma2_hat * shape
    sum(log(f_hat) + I_vals / f_hat)
  })
  H_grid[which.min(scores)]
}

# --- PATCHED compute_segment_posteriors with logging ---
compute_segment_posteriors <- function(path, L, a = 0.5) {
  dx <- diff(path)
  n <- length(dx)
  message(sprintf("[post] L=%d, n_dx=%d", L, n))
  if (n < L * 2) {
    message(sprintf("[post] L=%d -> n_dx < 2L (%d < %d) -> NULL", L, n, 2 * L))
    return(NULL)
  }
  N_seg <- floor(n / L)
  if (N_seg < 2) {
    message(sprintf("[post] L=%d -> N_seg=%d < 2 -> NULL", L, N_seg))
    return(NULL)
  }
  message(sprintf("[post] L=%d -> N_seg=%d, L*N_seg=%d", L, N_seg, L * N_seg))
  dx <- dx[1:(L * N_seg)]
  segments <- split(dx, rep(seq_len(L), each = N_seg))
  out <- lapply(seq_along(segments), function(k) {
    seg <- segments[[k]]
    s <- sign(seg)
    s[s == 0] <- NA
    crossings <- sum(s[-1] * s[-length(s)] < 0, na.rm = TRUE)
    N_pairs <- max(1, length(seg) - 1)
    alpha_post <- crossings + 1 - a
    beta_post  <- (N_pairs - crossings) + 1 - a
    eps <- 1e-3
    alpha_post <- max(alpha_post, eps)
    beta_post  <- max(beta_post, eps)
    alpha_post <- min(alpha_post, 100)
    beta_post  <- min(beta_post, 100)
    if (k == 1) {
      message(sprintf(
        "[post] L=%d seg1: crossings=%d, N_pairs=%d, alpha=%.4f, beta=%.4f",
        L, crossings, N_pairs, alpha_post, beta_post
      ))
    }
    list(
      alpha = alpha_post,
      beta = beta_post,
      N_pairs = N_pairs
    )
  })
  out
}

segment_H_summary <- function(alpha, beta, n_draws = 2000, alpha_ci = 0.05) {
  p_draws <- rbeta(n_draws, alpha, beta)
  p_draws <- pmin(pmax(p_draws, 1e-6), 1 - 1e-6)
  H_draws <- 1 + log2(cos(pi * p_draws / 2))
  H_mean <- mean(H_draws)
  H_low <- quantile(H_draws, probs = alpha_ci / 2)
  H_high <- quantile(H_draws, probs = 1 - alpha_ci / 2)
  c(H_mean = H_mean, H_low = H_low, H_high = H_high)
}

# --- PATCHED compute_segment_H_trajectory ---
compute_segment_H_trajectory <- function(path, L, a = 0.5, n_draws = 2000, alpha_ci = 0.05) {
  dx <- diff(path)
  n <- length(dx)
  if (n < L * 2) return(NULL)
  N_seg <- floor(n / L)
  if (N_seg < 2) return(NULL)
  idx <- seq_len(L * N_seg)
  dx <- dx[idx]
  segments <- split(dx, rep(seq_len(L), each = N_seg))
  res <- data.frame(
    segment = integer(L),
    mid_idx = numeric(L),
    H_mean = numeric(L),
    H_low = numeric(L),
    H_high = numeric(L)
  )
  for (k in seq_len(L)) {
    seg <- segments[[k]]
    s <- sign(seg)
    s[s == 0] <- NA
    crossings <- sum(s[-1] * s[-length(s)] < 0, na.rm = TRUE)
    N_pairs <- max(1, length(seg) - 1)
    alpha_post <- crossings + 1 - a
    beta_post  <- (N_pairs - crossings) + 1 - a
    eps <- 1e-3
    alpha_post <- max(alpha_post, eps)
    beta_post  <- max(beta_post, eps)
    alpha_post <- min(alpha_post, 100)
    beta_post  <- min(beta_post, 100)
    stats <- segment_H_summary(alpha_post, beta_post, n_draws = n_draws, alpha_ci = alpha_ci)
    res$segment[k] <- k
    res$mid_idx[k] <- 1 + (k - 0.5) * N_seg
    res$H_mean[k] <- stats["H_mean"]
    res$H_low[k] <- stats["H_low"]
    res$H_high[k] <- stats["H_high"]
  }
  res
}

# ==== numerically stable f_L + robust normalization with logging ====
compute_f_L <- function(post_list,
                        grid_p = seq(1e-4, 1 - 1e-4, length.out = 600)) {
  if (length(post_list) == 0) {
    message("[fL] post_list empty -> NULL")
    return(NULL)
  }
  mat <- sapply(post_list, function(pars) {
    logd <- dbeta(grid_p, pars$alpha, pars$beta, log = TRUE)
    logd_shift <- logd - max(logd, na.rm = TRUE)
    d <- exp(logd_shift)
    d[!is.finite(d)] <- 0
    d
  })
  f_L <- rowMeans(mat)
  f_L[!is.finite(f_L)] <- 0
  f_min <- min(f_L, na.rm = TRUE)
  f_max <- max(f_L, na.rm = TRUE)
  f_sum <- sum(f_L, na.rm = TRUE)
  message(sprintf("[fL] range=[%.4e, %.4e], sum(f)=%.4e", f_min, f_max, f_sum))
  data.frame(p = grid_p, f = f_L)
}

# --- PATCHED UNIMODALITY TEST (robust to NAs, treat ambiguity as TRUE) ---
is_unimodal <- function(f) {
  if (length(f) < 5) return(TRUE)

  # Smooth and remove non-finite values
  fs <- stats::filter(f, rep(1/5, 5), sides = 2)
  fs <- as.numeric(fs)
  fs <- fs[is.finite(fs)]
  if (length(fs) < 5) return(TRUE)

  # Derivative and remove non-finite values
  d <- diff(fs)
  d <- d[is.finite(d)]
  if (length(d) < 3) return(TRUE)

  d[abs(d) < 1e-6] <- 0
  signs <- sign(d)
  signs <- signs[signs != 0]

  if (length(signs) < 3) return(TRUE)

  res <- all(diff(signs) <= 0)
  if (is.na(res)) TRUE else res
}

# --- PATCHED select_optimal_L with logging and robust unimodal_flags ---
select_optimal_L <- function(path,
                             candidate_Ls = c(2, 4, 8, 16),
                             a = 0.5,
                             grid_p = seq(1e-4, 1 - 1e-4, length.out = 600)) {
  dx <- diff(path)
  message(sprintf(
    "[select_L] path_len=%d, dx_len=%d, candidate_Ls=%s",
    length(path), length(dx), paste(candidate_Ls, collapse = ",")
  ))
  results <- list()
  for (L in sort(unique(candidate_Ls))) {
    message(sprintf("[select_L] Trying L=%d", L))
    post <- tryCatch(
      compute_segment_posteriors(path, L, a),
      error = function(e) {
        message(sprintf("[select_L] ERROR in compute_segment_posteriors(L=%d): %s", L, e$message))
        NULL
      }
    )
    if (is.null(post)) {
      message(sprintf("[select_L] L=%d -> post=NULL (skipping)", L))
      next
    }
    message(sprintf("[select_L] L=%d -> got %d segments", L, length(post)))
    fL_df <- tryCatch(
      compute_f_L(post, grid_p),
      error = function(e) {
        message(sprintf("[select_L] ERROR in compute_f_L(L=%d): %s", L, e$message))
        NULL
      }
    )
    if (is.null(fL_df)) {
      message(sprintf("[select_L] L=%d -> fL_df=NULL (skipping)", L))
      next
    }
    f_min <- min(fL_df$f, na.rm = TRUE)
    f_max <- max(fL_df$f, na.rm = TRUE)
    f_sum <- sum(fL_df$f, na.rm = TRUE)
    message(sprintf("[select_L] L=%d -> fL range=[%.4e, %.4e], sum(f)=%.4e", L, f_min, f_max, f_sum))
    uni <- is_unimodal(fL_df$f)
    message(sprintf("[select_L] L=%d -> unimodal=%s", L, uni))
    results[[as.character(L)]] <- list(
      L = L,
      post = post,
      fL = fL_df,
      unimodal = uni
    )
  }
  if (length(results) == 0) {
    message("[select_L] No valid L found -> returning NULL")
    return(NULL)
  }
  Ls <- as.numeric(names(results))
  unimodal_flags_raw <- sapply(results, function(entry) entry$unimodal)
  # Treat NA as TRUE (only reject when clearly multimodal)
  unimodal_flags <- ifelse(is.na(unimodal_flags_raw), TRUE, unimodal_flags_raw)
  message(sprintf(
    "[select_L] Finished: candidate Ls=%s, unimodal_flags=%s",
    paste(Ls, collapse = ","),
    paste(unimodal_flags, collapse = ",")
  ))
  eligible <- names(results)[unimodal_flags & Ls > 1]
  if (length(eligible) > 0) {
    chosen <- results[[eligible[1]]]
    message(sprintf("[select_L] Choosing smallest unimodal L>1: L_star=%s", eligible[1]))
  } else {
    chosen <- results[[1]]
    message(sprintf("[select_L] No unimodal L>1; falling back to first L=%d", chosen$L))
  }
  chosen$L_star <- chosen$L
  chosen
}

normalize_density <- function(fL_df) {
  dp <- diff(fL_df$p)[1]
  density <- fL_df$f
  density[!is.finite(density)] <- 0          # remove NaN/Inf
  density <- pmax(density, 0)
  total <- sum(density * dp)
  if (!is.finite(total) || total <= 0) {
    density <- rep(1 / nrow(fL_df), nrow(fL_df))
    total <- sum(density * dp)
  }
  density / total
}

# --- PATCHED credible_p_interval: NA-safe, robust, clamps/interpolates ---
credible_p_interval <- function(fL_df, alpha = 0.05) {
  density <- normalize_density(fL_df)
  dp <- diff(fL_df$p)[1]

  # Basic sanity checks
  if (!is.finite(dp) || dp <= 0) {
    return(c(p_min = NA_real_, p_max = NA_real_))
  }

  cdf_raw <- cumsum(density * dp)

  if (any(!is.finite(cdf_raw))) {
    return(c(p_min = NA_real_, p_max = NA_real_))
  }

  # Clamp CDF to [0,1]
  cdf_raw <- pmin(pmax(cdf_raw, 0), 1)

  # Enforce weakly increasing CDF
  eps <- 1e-10
  for (i in 2:length(cdf_raw)) {
    if (cdf_raw[i] < cdf_raw[i - 1]) {
      cdf_raw[i] <- cdf_raw[i - 1]  # or + eps for strict
    }
  }

  keep <- !duplicated(cdf_raw)
  cdf <- cdf_raw[keep]
  p   <- fL_df$p[keep]

  # Degenerate case: CDF is essentially flat (bad, fallback to mean!)
  if (length(cdf) < 2 || diff(range(cdf)) < 1e-8) {
    # Fall back to a point mass at the posterior mean of p
    p_mean <- sum(p * density[keep] * dp) / sum(density[keep] * dp)
    return(c(p_min = p_mean, p_max = p_mean))
  }

  # Safe to interpolate for quantiles
  p_min <- approx(cdf, p, xout = alpha / 2,      rule = 2)$y
  p_max <- approx(cdf, p, xout = 1 - alpha / 2,  rule = 2)$y

  c(p_min = as.numeric(p_min), p_max = as.numeric(p_max))
}

# --- PATCHED map_p_to_H_interval: NA-safe, clamps boundaries, robust ---
map_p_to_H_interval <- function(p_min, p_max) {
  if (!is.finite(p_min) || !is.finite(p_max)) {
    return(c(H_min = NA_real_, H_max = NA_real_))
  }
  # Clamp p away from 0/1 (where cos(pi p/2) = 0) to avoid log(-)
  eps <- 1e-6
  p_min <- min(max(p_min, eps), 1 - eps)
  p_max <- min(max(p_max, eps), 1 - eps)

  H_min <- 1 + log2(cos(pi * p_max / 2))
  H_max <- 1 + log2(cos(pi * p_min / 2))

  # If something still explodes, mark as NA
  if (!is.finite(H_min) || !is.finite(H_max)) {
    H_min <- NA_real_
    H_max <- NA_real_
  }

  c(H_min = H_min, H_max = H_max)
}

EH_and_varH <- function(fL_df) {
  density <- normalize_density(fL_df)
  dp <- diff(fL_df$p)[1]
  g <- 1 + log2(cos(pi * fL_df$p / 2))
  g[!is.finite(g)] <- 0
  density[!is.finite(density)] <- 0
  EH  <- sum(g * density * dp)
  EH2 <- sum(g^2 * density * dp)
  if (!is.finite(EH))  EH  <- NA_real_
  if (!is.finite(EH2)) EH2 <- NA_real_
  varH <- max(EH2 - EH^2, 0, na.rm = TRUE)
  c(EH = EH, varH = varH)
}

# --- PATCHED compute_segmentation_summary: NA-robust, fallback if CI fails ---
compute_segmentation_summary <- function(path,
                                         candidate_Ls = c(2, 4, 8, 16),
                                         a = 0.5,
                                         alpha_ci = 0.05,
                                         grid_p = seq(1e-4, 1 - 1e-4, length.out = 600)) {
  dx <- diff(path)
  message(sprintf("[SegSummary] path_len=%d, dx_len=%d", length(path), length(dx)))
  if (length(dx) < 4) {
    message("[SegSummary] dx too short (<4) -> returning NULL")
    return(NULL)
  }
  optimal <- select_optimal_L(path, candidate_Ls = candidate_Ls, a = a, grid_p = grid_p)
  if (is.null(optimal)) {
    message("[SegSummary] select_optimal_L returned NULL")
    return(NULL)
  }
  message(sprintf("[SegSummary] Chosen L_star=%d", optimal$L_star))
  H_stats <- EH_and_varH(optimal$fL)   # returns c(EH, varH)
  H_mean  <- H_stats["EH"]
  # If we can't even get a mean, abort this frame
  if (!is.finite(H_mean)) {
    message("[SegSummary] EH is not finite -> returning NULL for this frame")
    return(NULL)
  }
  credible   <- credible_p_interval(optimal$fL, alpha = alpha_ci)
  H_interval <- map_p_to_H_interval(credible["p_min"], credible["p_max"])
  # If CI failed, use a degenerate interval at the mean
  if (!is.finite(H_interval["H_min"]) || !is.finite(H_interval["H_max"])) {
    H_interval["H_min"] <- H_mean
    H_interval["H_max"] <- H_mean
  }
  message(sprintf(
    "[SegSummary] Global H: mean=%.4f, low=%.4f, high=%.4f",
    H_mean, H_interval["H_min"], H_interval["H_max"]
  ))
  segments <- compute_segment_H_trajectory(
    path,
    optimal$L_star,
    a = a,
    n_draws = 1200,
    alpha_ci = alpha_ci
  )
  list(
    L_star = optimal$L_star,
    fL = optimal$fL,
    global = c(
      H_mean = H_mean,
      H_low  = H_interval["H_min"],
      H_high = H_interval["H_max"]
    ),
    segments = segments
  )
}

compute_estimator_trajectory <- function(path, steps) {
  steps <- sort(unique(steps))
  df <- data.frame(
    step = steps,
    zero = NA_real_,
    quad = NA_real_,
    rs = NA_real_,
    whittle = NA_real_
  )
  for (i in seq_along(steps)) {
    sample_path <- path[1:steps[i]]
    df$zero[i] <- estimate_zero_crossings(sample_path)
    df$quad[i] <- estimate_quadratic_variation(sample_path)
    df$rs[i] <- estimate_rs(sample_path)
    df$whittle[i] <- estimate_whittle(sample_path)
  }
  df
}

# --------------------------------------------------------------------
# Improved create_estimator_animation: Bayesian line and CI always visible
# (patched: robustly recompute H_mean/H_low/H_high from fL)
# --------------------------------------------------------------------
create_estimator_animation <- function(est_df,
                                      H_true,
                                      path,
                                      frame_indices,
                                      n_steps = max(est_df$step),
                                      candidate_Ls = c(2, 4, 8, 16)) {
  frames <- list()
  method_colors <- c(
    zero    = "#1f78b4",
    quad    = "#33a02c",
    rs      = "#e31a1c",
    whittle = "#ff7f00",
    bayes   = "#6a3d9a"
  )
  method_labels <- c(
    zero    = "Zero-crossings",
    quad    = "Quadratic variation",
    rs      = "Rescaled Range (R/S)",
    whittle = "Whittle",
    bayes   = "Bayesian segmentation"
  )
  data_methods <- c("zero", "quad", "rs", "whittle")
  steps <- est_df$step
  n_frames <- length(steps)
  bayes_mean <- rep(NA_real_, n_frames)
  bayes_low  <- rep(NA_real_, n_frames)
  bayes_high <- rep(NA_real_, n_frames)
  for (i in seq_len(n_frames)) {
    step <- steps[i]
    if (step < 50) {
      message(sprintf("[Bayes] Step %d ignored, path too short for segmentation summary", step))
      next
    }

    seg_summary <- compute_segmentation_summary(
      path[1:step],
      candidate_Ls = candidate_Ls
    )

    if (is.null(seg_summary) || is.null(seg_summary$fL)) {
      message(sprintf("[Bayes] No usable segmentation summary at frame %d (step=%d)", i, step))
      next
    }

    # Recompute global H from the returned fL to avoid any issues with seg_summary$global
    H_stats   <- EH_and_varH(seg_summary$fL)   # c(EH, varH)
    H_mean    <- as.numeric(H_stats["EH"])

    if (!is.finite(H_mean)) {
      message(sprintf("[Bayes] Recomputed EH not finite at frame %d (step=%d)", i, step))
      next
    }

    cred_p    <- credible_p_interval(seg_summary$fL, alpha = 0.05)
    H_int     <- map_p_to_H_interval(cred_p["p_min"], cred_p["p_max"])

    H_low  <- as.numeric(H_int["H_min"])
    H_high <- as.numeric(H_int["H_max"])

    # If CI failed, collapse to a point at the mean
    if (!is.finite(H_low) || !is.finite(H_high)) {
      H_low  <- H_mean
      H_high <- H_mean
    }

    message(sprintf(
      "[Bayes] Global (recomputed): H_mean=%.4f, H_low=%.4f, H_high=%.4f at frame %d (step=%d)",
      H_mean, H_low, H_high, i, step
    ))

    bayes_mean[i] <- H_mean
    bayes_low[i]  <- H_low
    bayes_high[i] <- H_high
  }
  message("[Bayes] Summary stats over all frames:")
  message(sprintf("  > Number of frames: %d", n_frames))
  message(sprintf("  > Non-NA Bayes means: %d", sum(!is.na(bayes_mean))))
  message(sprintf("  > Bayes means (first 10, na.rm=TRUE): %s", 
                  paste(formatC(na.omit(bayes_mean[1:min(length(bayes_mean),10)]), digits=4), collapse=", ")))
  combined_vals <- c(
    unlist(est_df[, data_methods], use.names = FALSE),
    bayes_mean, bayes_low, bayes_high,
    H_true
  )
  combined_vals <- combined_vals[is.finite(combined_vals)]
  if (length(combined_vals) == 0) {
    ylim <- c(0.1, 0.9)
    message("[Bayes] Y-limits defaulted to (0.1, 0.9) because no valid data")
  } else {
    y_min  <- min(combined_vals, na.rm = TRUE)
    y_max  <- max(combined_vals, na.rm = TRUE)
    margin <- max(0.05, 0.15 * (y_max - y_min))
    ylim   <- c(max(-1, y_min - margin), min(1.5, y_max + margin))
    message(sprintf("[Bayes] Y-limits set to: c(%.4f, %.4f)", ylim[1], ylim[2]))
  }
  for (i in seq_len(n_frames)) {
    step    <- est_df$step[i]
    current <- est_df[1:i, ]
    temp_file <- tempfile(fileext = ".png")
    png(temp_file, width = 700, height = 500, res = 100)
    par(mar = c(4, 4, 3, 1))
    plot(
      current$step,
      current$zero,
      type = "n",
      xlim = c(1, n_steps),
      ylim = ylim,
      xlab = "Sample size",
      ylab = "Estimated H",
      main = paste0("Estimators approaching H = ", round(H_true, 2)),
      sub  = paste0("Step: ", step, " / ", n_steps),
      cex.main = 1.4,
      cex.lab  = 1.2
    )
    grid(nx = 10, ny = 10, col = "lightgray", lty = 3)
    abline(h = H_true, col = "#4d4d4d", lty = 2, lwd = 1.5)
    for (method in data_methods) {
      values <- current[[method]]
      valid  <- !is.na(values)
      if (any(valid)) {
        lines(
          current$step[valid],
          values[valid],
          col = method_colors[[method]],
          lwd = 2
        )
        last_idx <- tail(which(valid), 1)
        points(
          current$step[last_idx],
          values[last_idx],
          col = method_colors[[method]],
          pch = 16,
          cex = 1.2
        )
      }
    }
    valid_bayes <- !is.na(bayes_mean[1:i])
    message(sprintf("[Bayes-Plot] Frame %d/%d: valid_bayes sum=%d", i, n_frames, sum(valid_bayes)))
    if (sum(valid_bayes) > 0) {
      message(sprintf("[Bayes-Plot] Drawing line for steps: %s", paste(steps[valid_bayes], collapse=", ")))
      message(sprintf("[Bayes-Plot] Drawing mean values: %s", paste(formatC(bayes_mean[valid_bayes], digits=4), collapse=", ")))
      lines(
        steps[valid_bayes],
        bayes_mean[valid_bayes],
        col = method_colors[["bayes"]],
        lwd = 2,
        lty = 2
      )
      last_idx <- tail(which(valid_bayes), 1)
      message(sprintf("[Bayes-Plot] Drawing last bayes point at step %d, mean %.4f, low=%.4f, high=%.4f",
                      steps[last_idx], bayes_mean[last_idx], bayes_low[last_idx], bayes_high[last_idx]))
      points(
        steps[last_idx],
        bayes_mean[last_idx],
        col = method_colors[["bayes"]],
        pch = 17,
        cex = 1.2
      )
      if (!is.na(bayes_low[last_idx]) && !is.na(bayes_high[last_idx])) {
        arrows(
          x0 = steps[last_idx],
          y0 = bayes_low[last_idx],
          x1 = steps[last_idx],
          y1 = bayes_high[last_idx],
          code = 3,
          angle = 90,
          length = 0.05,
          col = method_colors[["bayes"]],
          lwd = 1.1
        )
      } else {
        message(sprintf("[Bayes-Plot] CI missing for last_idx=%d: bayes_low=%.4f, bayes_high=%.4f",
                        last_idx, bayes_low[last_idx], bayes_high[last_idx]))
      }
    } else {
      message(sprintf("[Bayes-Plot] No bayesian estimate available to plot for frames 1:%d", i))
    }
    legend(
      "topright",
      legend = method_labels,
      col    = method_colors,
      lwd    = 2,
      lty    = c(1, 1, 1, 1, 2),
      pch    = c(16, 16, 16, 16, 17),
      bty    = "n",
      cex    = 0.8
    )
    usr      <- par("usr")
    label_y  <- usr[4] - 0.05 * diff(usr[3:4])
    y_step   <- 0.05 * diff(usr[3:4])
    label_x  <- usr[1] + 0.02 * diff(usr[1:2])
    idx_lab  <- 0
    for (method in data_methods) {
      last_val <- tail(na.omit(current[[method]]), 1)
      if (length(last_val) == 0) next
      text(
        x = label_x,
        y = label_y - idx_lab * y_step,
        labels = paste0(
          method_labels[[method]], ": ",
          formatC(last_val, digits = 2, format = "f")
        ),
        col = method_colors[[method]],
        adj = c(0, 1),
        cex = 0.85
      )
      idx_lab <- idx_lab + 1
    }
    dev.off()
    frames[[length(frames) + 1]] <- image_read(temp_file)
    unlink(temp_file)
  }
  animation <- do.call(c, frames)
  image_animate(animation, fps = 10, loop = 1)
}
```
```{r}
#| context: server

anim_state <- reactiveValues(
  animation = NULL,
  H = 0.5,
  n_steps = 500,
  method = "cholesky"
)
est_state <- reactiveValues(
  animation = NULL
)

observeEvent(input$resample, {
  H <- input$hurst
  n <- input$n_steps
  method <- input$method
  
  frame_indices <- calc_frame_indices(n)

  path <- generate_fbm(n, H, method)
  anim_state$animation <- create_animation_frames(path, H, n, method, frame_indices)
  anim_state$H <- H
  anim_state$n_steps <- n
  anim_state$method <- method

  est_state$animation <- NULL
  estimator_df <- compute_estimator_trajectory(path, frame_indices)
  est_state$animation <- create_estimator_animation(
    estimator_df,
    H,
    path,
    frame_indices,
    n_steps = n
  )
})

output$simulation_plot <- renderImage({
  if (is.null(anim_state$animation)) {
    temp_file <- tempfile(fileext = ".png")
    png(temp_file, width = 700, height = 500)
    plot(0, type = "n", axes = FALSE, xlab = "", ylab = "")
    text(0.5, 0.5, "wait", 
         cex = 1.5, col = "gray60")
    dev.off()
    
    list(src = temp_file, contentType = "image/png", width = "100%", height = "450px")
  } else {
    temp_file <- tempfile(fileext = ".gif")
    image_write(anim_state$animation, path = temp_file)
    
    list(src = temp_file, contentType = "image/gif", width = "100%", height = "450px", deleteFile = TRUE)
  }
}, deleteFile = FALSE)

output$estimator_plot <- renderImage({
  if (is.null(est_state$animation)) {
    temp_file <- tempfile(fileext = ".png")
    png(temp_file, width = 700, height = 500)
    plot(0, type = "n", axes = FALSE, xlab = "", ylab = "")
    text(0.5, 0.5, "wait", 
         cex = 1.3, col = "gray50")
    dev.off()
    
    list(src = temp_file, contentType = "image/png", width = "100%", height = "450px")
  } else {
    temp_file <- tempfile(fileext = ".gif")
    image_write(est_state$animation, path = temp_file)
    
    list(src = temp_file, contentType = "image/gif", width = "100%", height = "450px", deleteFile = TRUE)
  }
}, deleteFile = FALSE)
```
::: {style="display: flex; gap: 20px; align-items: flex-start;"}

::: {style="flex: 0 0 25%; background-color: #f0f0f0; padding: 15px; border-radius: 5px;"}

```{r}
# Add CSS styling for compact controls
tags$style(HTML("
  .compact-controls .form-group {
    margin-bottom: 12px;
  }
  .compact-controls .form-group label {
    font-size: 13px;
    margin-bottom: 4px;
  }
  .compact-controls select {
    font-size: 12px;
    padding: 4px 6px;
  }
  .compact-controls input[type='slider'] {
    font-size: 12px;
  }
  .compact-controls .shiny-input-container {
    font-size: 12px;
  }
"))

tags$div(class = "compact-controls",
  selectInput(
    "method",
    strong("Simulation Method:"),
    choices = list(
      "Cholesky (Exact, O(N³))" = "cholesky",
      "Hosking (Exact, O(N²))" = "hosking",
      "FFT (Near-exact, O(N log N))" = "fft",
      "Spectral (Approx, O(N log N))" = "spectral",
      "Corr. Random Walk (Approx, fast)" = "crw",
      "Random Midpoint Displacement (Exact*)" = "rmd"
    ),
    selected = "cholesky"
  ),
  
sliderInput(
  "hurst",
  strong("Hurst Exponent (H):"),
  min = 0.1,
  max = 0.9,
  value = 0.5,
  step = 0.05
  ),

numericInput(
  "n_steps",
  strong("Time Steps:"),
  value = 500,
  min = 100,
  max = 2000,
  step = 100
  ),

actionButton("resample", "Generate Sample", class = "btn-primary btn-block", width = "100%")
)
```

:::

::: {style="flex: 1; min-width: 0;"}

```{r}
plotOutput("simulation_plot", height = "320px")
```


```{r}
imageOutput("estimator_plot", height = "320px")
```
</div>

:::

:::

---
