
source("useful_functions.R")

# ------------------------------------------
# Generate results
# ------------------------------------------

experiment = 0 # 0 implies skip experiments and go straight to plotting.
plot_first = TRUE # If TRUE, then either Figure 3 or 4 is produced. If FALSE, Figure 5 is produced.
normalize = FALSE # If TRUE, then Figure 4 is produced instead of Figure 3. This variable is only used if plot_first is TRUE.

if(experiment == 1){
  # extent is an output of extent_profile, and proportion is a vector of entries in [0,1].
  get_extremal_range = function(extent, proportion = 1){
    res = rep(Inf, length(proportion))
    start_ind = 1
    for (pix in sort(proportion, decreasing = TRUE, index.return = TRUE)$ix){
      for (i in start_ind:length(extent$y)){
        if (extent$y[i] < proportion[pix]){
          res[pix] = extent$x[i]
          break
        }
      }
      start_ind = i
    }
    return(res)
  }
  
  generator_path = "gaussian_generator.RData"
  if (file.exists(generator_path)){
    load(generator_path)
    print("Generator loaded!")
  } else {
    x = y = seq(-0.5, 0.5, length.out = 121)
    time2generateRF(x,y)
    cov_model = covModel("Matern", nu = 2.5)
    gen = randomFieldGenerator(x, y, cov_model)
    save(gen, file = generator_path)
  }
  
  
  grid = expand.grid(gen$x, gen$y)
  names(grid) = c("x", "y")
  origin = which(grid$x == 0 & grid$y == 0)
  # ps = seq(0.1,0.9,0.1)
  # us = qnorm(ps)
  # us = 1:10
  us = seq(2,6,0.5)
  
  
  
  # Make a bunch of realizations and check which ones have an exceedance at 0.
  N_fields = 5000
  results = matrix(rep(NA, length(us)*N_fields), nrow = N_fields)
  # ep = extent_profile(chi2RF(gen, thresh = 10), 9, max_dist = 0.2)
  # eps = lapply(us, function(u) return(matrix(rep(NA, length(ep$x)*N_fields), nrow = N_fields)))
  for (i in 1:N_fields){
    if (i %% 10 == 0) print(paste(i, "/", N_fields))
    # RF = studentRF(gen)
    # RF = generateRF(gen)
    # RF = chi2RF(gen)
    RF = mixtureRF(gen, alpha = 2)
    value_at_origin = RF$Z[origin]
    for (j in 1:length(us)){
      u = us[j]
      if (u < value_at_origin){
        ep = extent_profile(RF, u, max_dist = 0.2)
        er = get_extremal_range(ep, 1)
        results[i,j] = er
      } else {
        break
      }
    }
  }
  save(results, us, file = "matrix_of_extremal_ranges_mixture_2to6.RData")
  # x = ep$x
  # save(eps, x, us, file = "matrix_of_extent_profiles_mixture.RData")
} else if(experiment == 2){
  # extent is an output of extent_profile, and proportion is a vector of entries in [0,1].

  generator_path = "gaussian_generator_conditional.RData"
  if (file.exists(generator_path)){
    load(generator_path)
    print("Generator loaded!")
  } else {
    x = y = seq(-0.5, 0.5, length.out = 121)
    time2generateRF(x,y)
    cov_model = covModel("Matern", nu = 2.5)
    gen = randomFieldGenerator(x, y, cov_model, cond_sites = data.frame(x = 0, y = 0))
    save(gen, file = generator_path)
  }
  
  
  grid = expand.grid(gen$x, gen$y)
  names(grid) = c("x", "y")
  origin = which(grid$x == 0 & grid$y == 0)
  
  # us = qnorm(pexp(seq(2,12,length.out = 10)))
  us = seq(10,300, length.out = 10)
  # us = seq(2,6,0.5)
  
  N_fields = 500
  ep = extent_profile(chi2RF(gen, thresh = 10), 9, max_dist = 0.2)
  eps = lapply(us, function(u) return(matrix(rep(NA, length(ep$x)*N_fields), nrow = N_fields)))
  for (i in 1:N_fields){
    if (i %% 10 == 0) print(paste(i, "/", N_fields))
    for (j in 1:length(us)){
      u = us[j]
      # RF = studentRF(gen, thresh = u)
      # RF = generateRF(gen, cond_val = gauss_above(u))
      # RF = chi2RF(gen, thresh = u)
      RF = mixtureRF(gen, thresh = u, alpha = 2)
      
      ep = extent_profile(RF, u, max_dist = 0.2)
      eps[[j]][i,] = ep$y
    }
  }
  x = ep$x
  save(eps, x, us, file = "matrix_of_extent_profiles_mix_2.RData")
}

# ------------------------------------------
# Plotting results
# ------------------------------------------

library(stats)  # For gamma and normal distribution functions

# Function to compute E[L_k]
compute_C_star <- function(k, u, alpha, lambda) {
  omega = function(j){
    return(pi^(j/2) / gamma(1  + j/2))
  }
  
  special_choose = function(n, k){
    return(choose(n, k) * omega(n) / omega(j) / omega(n - j))
  }
  
  # Function for the lower incomplete gamma function
  lower_gamma <- function(a, x) {
    if (x <= 0) return (0)
    integrand <- function(t) {
      t^(a - 1) * exp(-t)
    }
    result <- integrate(integrand, lower = 0, upper = x)
    return(result$value)
  }
  
  # Function to compute E[rho_j(h(u, Lambda))]
  E_rho <- function(j, u, alpha) {
    if (j == 0) {
      return(u^(-2*alpha) * 2 ^ (alpha - 1) * (pi^(-1/2)) * lower_gamma(alpha + 1/2, u^2/2) + 1 - pnorm(u))
    } else if (j == 1) {
      return(u^(-2*alpha) * 2 ^ (alpha - 1) * (pi^(-1)) * alpha * lower_gamma(alpha, u^2/2))
    } else if (j == 2) {
      return(u^(-2*alpha) * 2 ^ (alpha - 1) * (pi^(-3/2)) * alpha * lower_gamma(alpha + 1/2, u^2/2))
    } else {
      stop("Invalid j")
    }
  }
  
  j = 2 - k
  return(special_choose(2, j) * lambda^(j/2) * E_rho(j, u, alpha))
}

create_plot_information = function(type){
  if (type == "gauss"){
    k = Inf
  } else {
    k = 3
  }
  ds = ds_upper = ds_lower = rep(0, length(us))
  for (i in 1:length(us)){
    u = us[i]
    x = results[,i]
    x = x[which(!is.na(x))]
    if (is.infinite(min(x))){
      ds[i] = ds_upper[i] = ds_lower[i] = 0
    } else {
      ds[i] = density_at_zero(x)
      N_bootstraps = 200
      bootstraps = rep(0, N_bootstraps)
      for (j in 1:N_bootstraps){
        bootstrap = sample(1:length(x), replace = TRUE)
        bootstraps[j] = density_at_zero(x[bootstrap])
      }
      alpha = 0.95
      ds_upper[i] = quantile(bootstraps, 0.5 + alpha / 2)
      ds_lower[i] = quantile(bootstraps, 0.5 - alpha / 2)
    }
  }
  data <- data.frame(
    u = us,
    density = ds,
    lower = ds_lower,
    upper = ds_upper
  )
  
  nu = 2.5
  lambda = nu/(nu-1)
  
  if (type == "chi2"){
    C1_star = function(u){
      u = (u-k)/sqrt(2*k) # unnormalized
      return(sqrt(lambda*pi)/2^((k+1)/2)/gamma(k/2)*(k+u*sqrt(2*k))^((k-1)/2)*exp(-(k + u*sqrt(2*k))/2))
    }
    C2_star = function(u){
      u = (u-k)/sqrt(2*k) # unnormalized
      return(1-pchisq(k + u*sqrt(2*k), df = k))
    }
    slope_star = function(u){
      return(2*C1_star(u)/C2_star(u))
    }
    x = seq(0,6.5,0.001)
    slope_data <- data.frame(
      x = x,
      y = slope_star(x)
    )
  } else if (type == "gauss" || type == "student"){
    C1_star = function(u, k = Inf){
      if (is.infinite(k)){
        return(1/4*sqrt(lambda)*exp(-u^2/2))
      }
      u = u * sqrt((k-2)/k) # unnormalized
      return(1/4*sqrt(lambda)*(1+u^2/(k-2))^((1-k)/2))
    }
    C2_star = function(u, k = Inf){
      if (is.infinite(k)){
        return(1-pnorm(u))
      }
      # return(1-pt(u*sqrt(k/(k-2)), df = k)) # normalized
      return(1-pt(u, df = k)) # unnormalized
    }
    slope_star = function(u, k = Inf){
      return(2*C1_star(u, k)/C2_star(u, k))
    }
    x = seq(-2.5,2.7,0.001)
    slope_data <- data.frame(
      x = x,
      y = slope_star(x, k)
    )
  } else if (type == "mixture"){
    slope_star = function(u){
      alpha = 1
      return(2*compute_C_star(1,u,alpha,lambda) / compute_C_star(2,u,alpha,lambda))
    }
    x = seq(-1.5,6.5,0.001)
    slope_data <- data.frame(
      x = x,
      y = sapply(x, slope_star)
    )
  }
  return(list(
    theoretical = slope_data,
    empirical = data
  ))
}


second_mix_x = seq(1,100,0.01)
nu = 2.5
lambda = nu/(nu-1)
wrapper = function(u){
  return(compute_C_star(2, u, alpha = 1, lambda = lambda))
}
C2 = sapply(second_mix_x, wrapper)
wrapper = function(u){
  return(compute_C_star(1, u, alpha = 1, lambda = lambda))
}
C1 = sapply(second_mix_x, wrapper)
second_mix_y = -2*C1/C2/pi

second_gauss_x = qnorm(pexp(seq(0.5,15,0.01)))
C1 = 1/4*sqrt(lambda)*exp(-second_gauss_x^2/2)
C2 = 1-pnorm(second_gauss_x)
second_gauss_y = -2*C1/C2/pi

second_stud_x = seq(1,18,0.01)
k = 3
u = second_stud_x * sqrt((k-2)/k) # unnormalized
C1 = 1/4*sqrt(lambda)*(1+u^2/(k-2))^((1-k)/2)
C2 = 1-pt(second_stud_x, df = k)
second_stud_y = -2*C1/C2/pi

second_chi_x = seq(1,30,0.01)
k = 3
u = (second_chi_x-k)/sqrt(2*k) # unnormalized
C1 = sqrt(lambda*pi)/2^((k+1)/2)/gamma(k/2)*(k+u*sqrt(2*k))^((k-1)/2)*exp(-(k + u*sqrt(2*k))/2)
C2 = 1-pchisq(second_chi_x, df = k)
second_chi_y = -2*C1/C2/pi

create_second_plot_information = function(){
  # us = 1:10
  slopes = slopes_upper = slopes_lower = rep(0, length(us))
  for (i in 1:length(us)){
    u = us[i]
    print(paste(
      "processing u =", u 
    ))
    ep = data.frame(
      x = x,
      y = apply(eps[[i]], MARGIN = 2, mean)
    )
    # plot(ep)
    tdf = tail_dependence_fun(ep)
    # lines(tdf)
    slopes[i] = tdf$slope_at_0
    N_bootstraps = 200
    bootstraps = rep(0, N_bootstraps)
    for (j in 1:N_bootstraps){
      bootstrap = sample(1:nrow(eps[[i]]), replace = TRUE)
      ep = data.frame(
        x = x,
        y = apply(eps[[i]][bootstrap,], MARGIN = 2, mean)
      )
      # plot(ep)
      tdf = tail_dependence_fun(ep)
      # lines(tdf)
      bootstraps[j] = tdf$slope_at_0
    }
    alpha = 0.95
    slopes_upper[i] = quantile(bootstraps, 0.5 + alpha / 2)
    slopes_lower[i] = quantile(bootstraps, 0.5 - alpha / 2)
  }
  data <- data.frame(
    u = us,
    slopes = slopes,
    lower = slopes_lower,
    upper = slopes_upper
  )
  return(data)
}

library(ggplot2)
set.seed(55)

# Function to compute the CDF of Z = X * Y
pmixture <- function(z, alpha = 1) {
  inner = function(z){
    # Integrand for the first term (Gaussian CDF)
    integrand <- function(x) {
      pnorm(z/x)*alpha/x^(alpha+1)
    }
    integral = integrate(integrand, 1, Inf)
    return(integral$value)
  }
  return(sapply(z, inner))
}

if (plot_first){
  load("matrix_of_extremal_ranges_student.RData")
  stud = create_plot_information("student")
  stud$empirical$u = stud$empirical$u * sqrt(3)
  if(normalize){
    stud$theoretical$x = qexp(pt(stud$theoretical$x, df = 3))
    stud$empirical$u = qexp(pt(stud$empirical$u, df = 3))
  }
  load("matrix_of_extremal_ranges.RData")
  gauss = create_plot_information("gauss")
  if(normalize){
    gauss$theoretical$x = qexp(pnorm(gauss$theoretical$x))
    gauss$empirical$u = qexp(pnorm(gauss$empirical$u))
  }
  load("matrix_of_extremal_ranges_chi2.RData")
  chi2 = create_plot_information("chi2")
  chi2$empirical$u = 3 + chi2$empirical$u * sqrt(6)
  if(normalize){
    chi2$theoretical$x = qexp(pchisq(chi2$theoretical$x, df = 3))
    chi2$empirical$u = qexp(pchisq(chi2$empirical$u, df = 3))
  }
  load("matrix_of_extremal_ranges_mixture.RData")
  us_temp = us
  results_temp = results
  load("matrix_of_extremal_ranges_mixture_2to6.RData")
  us = c(us_temp, us)
  results = cbind(results_temp, results)
  mix = create_plot_information("mixture")
  if(normalize){
    mix$theoretical$x = qexp(pmixture(mix$theoretical$x, alpha = 2))
    mix$empirical$u = qexp(pmixture(mix$empirical$u, alpha = 2))
  }
  
  size1 = 1.2
  size2 = 0.8
  type1 = "dotdash"
  type2 = "solid"
  
  if (normalize) {x_axis_name = ""} else x_axis_name = "Threshold"
  
  p = ggplot() +
    geom_point(data = gauss$empirical, aes(x = u, y = density, color = "Gaussian"), shape = 16, size = 2) +  # Gaussian points
    geom_errorbar(data = gauss$empirical, aes(x = u, ymin = lower, ymax = upper, color = "Gaussian"), 
                  width = 0.3, linewidth = size1, linetype = type2) +  # Gaussian error bars (thickest)
    geom_line(data = gauss$theoretical, aes(x = x, y = y, color = "Gaussian"), linewidth = size1, linetype = type2) +  # Gaussian theoretical line (thickest)
    
    geom_point(data = stud$empirical, aes(x = u, y = density, color = "Student"), shape = 18, size = 2) +  # Student points
    geom_errorbar(data = stud$empirical, aes(x = u, ymin = lower, ymax = upper, color = "Student"), 
                  width = 0.3, linewidth = size2, linetype = type1) +  # Student error bars (medium-thick)
    geom_line(data = stud$theoretical, aes(x = x, y = y, color = "Student"), linewidth = size2, linetype = type1) +  # Student theoretical line (medium-thick)
    
    geom_point(data = chi2$empirical, aes(x = u, y = density, color = "Chi-squared"), shape = 17, size = 2) +  # Chi-squared points
    geom_errorbar(data = chi2$empirical, aes(x = u, ymin = lower, ymax = upper, color = "Chi-squared"), 
                  width = 0.3, linewidth = size1, linetype = type1) +  # Chi-squared error bars (medium)
    geom_line(data = chi2$theoretical, aes(x = x, y = y, color = "Chi-squared"), linewidth = size1, linetype = type1) +  # Chi-squared theoretical line (medium)
    
    geom_point(data = mix$empirical, aes(x = u, y = density, color = "Mixture"), shape = 15, size = 2) +  # Mixture points
    geom_errorbar(data = mix$empirical, aes(x = u, ymin = lower, ymax = upper, color = "Mixture"), 
                  width = 0.3, linewidth = size2, linetype = type2) +  # Mixture error bars (thinnest)
    geom_line(data = mix$theoretical, aes(x = x, y = y, color = "Mixture"), linewidth = size2, linetype = type2) +  # Mixture theoretical line (thinnest)
    
    # ylim(0, 6) +  # Set y-axis limits
    labs(x = x_axis_name, y = "Slope of CDF", color = "Field Type", linewidth = "Field Type", linetype = "Field type") +  # Axis labels and legend titles
    scale_color_manual(values = c("Gaussian" = "red", "Student" = "blue", 
                                  "Chi-squared" = "orange", "Mixture" = "violet")) +  # Custom colors
    theme_minimal()  # Minimal theme
  
  print(p)
  
} else {
  x = seq(6,7, 0.1)
  y = qexp(pnorm(x))
  res = lm(log(y)~log(x))
  gauss_to_exp = function(u){
    inner = function(u){
      if (u < 7) return(qexp(pnorm(u)))
      ans = res$coefficients[1] + res$coefficients[2] * log(u)
      return(as.numeric(exp(ans)))
    }
    return(sapply(u, inner))
  }
  
  load("matrix_of_extent_profiles_gauss_2.RData")
  gauss = create_second_plot_information()
  gauss$u = gauss_to_exp(gauss$u)
  second_gauss_x = gauss_to_exp(second_gauss_x)
  load("matrix_of_extent_profiles_student.RData")
  us = 1:10
  stud = create_second_plot_information()
  stud$u = stud$u * sqrt(3)
  stud$u = qexp(pt(stud$u, df = 3))
  second_stud_x = qexp(pt(second_stud_x, df = 3))
  load("matrix_of_extent_profiles_chi2.RData")
  chi2 = create_second_plot_information()
  chi2$u = 3 + chi2$u * sqrt(6)
  chi2$u = qexp(pchisq(chi2$u, df = 3))
  second_chi_x = qexp(pchisq(second_chi_x, df = 3))
  load("matrix_of_extent_profiles_mix_2.RData")
  mix = create_second_plot_information()
  mix$u = qexp(pmixture(mix$u, alpha = 2))
  second_mix_x = qexp(pmixture(second_mix_x, alpha = 2))
  second_mix_x = c(second_mix_x, 12.5)
  second_mix_y = c(second_mix_y, second_mix_y[length(second_mix_y)])
  
  
  size1 = 1.2
  size2 = 0.8
  type1 = "dotdash"
  type2 = "solid"
  
  # Creating the plot
  p = ggplot() +
    geom_point(data = gauss, aes(x = u, y = slopes, color = "Gaussian"), shape = 16, size = 2) +
    geom_errorbar(data = gauss, aes(x = u, ymin = lower, ymax = upper, color = "Gaussian"), 
                  width = 0.3, linewidth = size1, linetype = type2) +
    geom_point(data = stud, aes(x = u, y = slopes, color = "Student"), shape = 18, size = 2) +
    geom_errorbar(data = stud, aes(x = u, ymin = lower, ymax = upper, color = "Student"), 
                  width = 0.3, linewidth = size2, linetype = type1) +
    geom_point(data = chi2, aes(x = u, y = slopes, color = "Chi-squared"), shape = 17, size = 2) +
    geom_errorbar(data = chi2, aes(x = u, ymin = lower, ymax = upper, color = "Chi-squared"), 
                  width = 0.3, linewidth = size1, linetype = type1) +
    geom_point(data = mix, aes(x = u, y = slopes, color = "Mixture"), shape = 15, size = 2) +
    geom_errorbar(data = mix, aes(x = u, ymin = lower, ymax = upper, color = "Mixture"), 
                  width = 0.3, linewidth = size2, linetype = type2) +
    
    # Adding curves to the plot
    geom_line(aes(x = second_gauss_x, y = second_gauss_y, color = "Gaussian"), linewidth = size1, linetype = type2) +
    geom_line(aes(x = second_stud_x, y = second_stud_y, color = "Student"), linewidth = size2, linetype = type1) +
    geom_line(aes(x = second_chi_x, y = second_chi_y, color = "Chi-squared"), linewidth = size1, linetype = type1) +
    geom_line(aes(x = second_mix_x, y = second_mix_y, color = "Mixture"), linewidth = size2, linetype = type2) +
    
    labs(x = "", y = "", color = "Field Type", linewidth = "Field Type", linetype = "Field Type") +
    scale_color_manual(values = c("Gaussian" = "red", "Student" = "blue", 
                                  "Chi-squared" = "orange", "Mixture" = "violet")) +
    # xlim(c(0,12)) +
    theme_minimal()
  
  print(p)
}
