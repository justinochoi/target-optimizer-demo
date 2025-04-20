library(tidyverse) 
library(mgcv) 
library(arrow)
library(mvtnorm) 

# four-seam fastballs only 
data <- read.csv('/Users/justinchoi/statcast_2024.csv') |>  
  filter(pitch_type == 'FF')

# add indicator variable for swing and count, re-define zones  
data <- data |> 
  drop_na(plate_x, plate_z, balls, strikes, description) |> 
  mutate(
    plate_x = plate_x * 12, 
    plate_z = plate_z * 12, 
    is_swing = case_when(
       description %in% c('foul','foul_tip','hit_into_play','swinging_strike') ~ 1, 
      .default = 0), 
    count = paste(balls, strikes, sep="-"), 
    count_type = case_when(
      count == '0-0' ~ 'first_pitch', 
      count == '3-2' ~ 'full_count', 
      count %in% c('0-1','1-0') ~ 'early', 
      count %in% c('0-2','1-2','2-2') ~ 'two_strikes', 
      .default = 'behind'
    ), 
    zone = case_when(
      plate_x < -3.32 & plate_z > 35.48 ~ 1, 
      between(plate_x, -3.32, 3.32) & plate_z > 35.48 ~ 2, 
      plate_x > 3.32 & plate_z > 35.48 ~ 3, 
      plate_x < -3.32 & between(plate_z, 26.89, 35.48) ~ 4, 
      plate_x > 3.32 & between(plate_z, 26.89, 35.48) ~ 6, 
      plate_x < -3.32 & plate_z < 26.89 ~ 7, 
      between(plate_x, -3.32, 3.32) & plate_z < 26.89 ~ 8, 
      plate_x > 3.32 & plate_z < 26.89 ~ 9, 
      .default = 5 
    )
  ) 

# convert to factor, re-level 
data$zone <- as.factor(data$zone) 
data$count_type <- as.factor(data$count_type) 
data$zone <- relevel(data$zone, ref = 5)
data$count_type <- relevel(data$count_type, ref = 'first_pitch')
soto <- filter(data, batter == 665742)

# build swing model 
soto_swing_mod <- gam(is_swing ~ zone + count_type, family = binomial, 
                       data = soto) 
summary(soto_swing_mod)

# table of mean run value by zone, count, swing/no swing for entire league 
mean_rv_table <- data |>  
  group_by(zone, count_type, is_swing) |>  
  summarize(
    mean_rv = mean(delta_run_exp, na.rm = T) 
  ) 

# add soto's predicted swing rate 
mean_rv_table$swing_rate = predict(soto_swing_mod, 
                                   select(mean_rv_table, zone, count_type), 
                                   type = 'response') 

# table of expected run value by zone and count for soto  
soto_rv_table <- mean_rv_table |> 
  mutate(
    pred_rv = case_when(
      is_swing == 1 ~ mean_rv * swing_rate,  
      .default = mean_rv * (1 - swing_rate)
    )
  ) |> 
  group_by(zone, count_type) |>  
  summarize(
    exp_rv = round(sum(pred_rv), 5) * 100 
  )

# find the lower and upper bounds of each zone 
# technically some of these limits should be +/-Inf but it's fine
zone_list <- c(1:9)
zone_bounds_finder <- function(data) {
  result = tibble() 
  for (i in zone_list) {
    zone_data <- filter(data, zone == i) %>% select(plate_x, plate_z) 
    bounds = data.frame(
      zone = i, 
      min_x = min(zone_data$plate_x), 
      max_x = max(zone_data$plate_x), 
      min_z = min(zone_data$plate_z), 
      max_z = max(zone_data$plate_z)
    )
    result = bind_rows(result, bounds)
  }
  return(result)
}
zone_bounds <- zone_bounds_finder(data) 
gilbert <- filter(data, pitcher == 669302) 

# calculate area under bivariate normal pdf 
riemann_approx <- function(mu, sigma, xlim, zlim, n) {
  
  x_seq <- seq(xlim[1], xlim[2], length.out = n)
  z_seq <- seq(zlim[1], zlim[2], length.out = n)
  dx <- diff(x_seq)[1]
  dz <- diff(z_seq)[1]
  
  grid <- expand.grid(plate_x = x_seq, plate_z = z_seq) 
  grid <- grid |> mutate(
    value = mapply(function(x, y) dmvnorm(c(x,y), mu, sigma, log = F), 
                   plate_x, plate_z) 
  )
  riemann_sum <- sum(grid$value) * dx * dz 
  return(riemann_sum)
}

# returns step_size^2 x 12 dataframe 
target_optimizer <- function(step_size) {
  
  # gather necessary data 
  pitcher_data <- filter(gilbert, count_type == 'full_count') |>  
    select(plate_x, plate_z)
  cov_matrix <- cov(pitcher_data) 
  batter_data <- filter(soto_rv_table, count_type == 'full_count') 
  
  # create list of targets
  min_x <- -10; max_x <- 10; min_z <- 18; max_z <- 44
  x_seq <- seq(min_x, max_x, length.out = step_size)
  z_seq <- seq(min_z, max_z, length.out = step_size) 
  targets <- expand.grid(x_seq, z_seq) 
  
  # for each available zone, calculate RV in each target 
  results <- tibble(targets[1], targets[2]) 
  for (j in zone_list) {
    rv <- filter(batter_data, zone == j)$exp_rv
    bounds <- filter(zone_bounds, zone == j) 
    min_x <- bounds$min_x; min_z <- bounds$min_z
    max_x <- bounds$max_x; max_z <- bounds$max_z 
    x_lim <- c(min_x, max_x); z_lim <- c(min_z, max_z)
    probs <- vector(length = step_size^2) 
    for (i in 1:nrow(results)) { 
      x <- targets[i,1]; z <- targets[i,2]
      probs[i] <- riemann_approx(c(x,z), cov_matrix, x_lim, z_lim, step_size)
    } 
    zone_rvs <- probs * rv  
    results <- bind_cols(results, zone_rvs) 
  } 
  
  # rename columns that are automatically altered by R (grrr)
  results <- results |> 
    rename(
      plate_x = Var1, 
      plate_z = Var2, 
      zone_1 = ...3, 
      zone_2 = ...4, 
      zone_3 = ...5, 
      zone_4 = ...6, 
      zone_5 = ...7, 
      zone_6 = ...8, 
      zone_7 = ...9, 
      zone_8 = ...10, 
      zone_9 = ...11, 
    ) |> 
    mutate(rv_sum = rowSums(across(zone_1:zone_9))) 
  
  return(results)
}

test <- target_optimizer(10) # works as expected! 

# find target corresponding to mimumum RV 
min_index <- which(test$rv_sum == min(test$rv_sum))
best_target <- tibble(
  plate_x = test$plate_x[min_index], 
  plate_z = test$plate_z[min_index] 
)

# this will be probably part of the optimizer function in the final product 
gilbert_cov <- cov(
  filter(gilbert, count_type == 'full_count') %>% select(plate_x, plate_z)
)

# simulate pitches from Gilbert based on target + covariance matrix 
set.seed(56)
gilbert_sim <- rmvnorm(n=500, mean = c(best_target$plate_x, best_target$plate_z), 
                       sigma = gilbert_cov) %>% 
  as_tibble() %>% rename(plate_x = V1, plate_z = V2) 

# visualize! 
ggplot(gilbert_sim, aes(x=plate_x, y=plate_z)) + 
  stat_density_2d(
    aes(x=plate_x, y=plate_z, fill=after_stat(density)), 
    geom = "raster", contour = F
  ) + 
  geom_rect(aes(xmin = -10, xmax = 10, ymin = 18, ymax = 44), 
            fill=NA, color='black') + 
  geom_point(data = best_target, shape = 4, size = 5, color='white') + 
  scale_fill_gradient2(low='yellow', mid='orange', high='red', midpoint = 0.0005) + 
  coord_fixed(ratio = 1, xlim = c(-20, 20), ylim = c(0, 50)) + 
  theme_bw() + 
  labs(
    x = 'Horizontal Location (in.)', 
    y = 'Vertical Location (in.)', 
    title = 'Where Should Logan Gilbert Aim vs. Juan Soto?', 
    subtitle = 'Full count, four-seam fastball', 
    fill = 'Simulated Density'
  ) + 
  theme(plot.title = element_text(face = 'bold'))


