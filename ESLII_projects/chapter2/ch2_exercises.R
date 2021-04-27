#### Chapter 2: Overview of Supervised Learning Methods ####
library(ggplot2) # For plotting
library(MASS) # For the bivariate normal function

# The basic linear equation
equation_2.1 = function(x, beta, beta_null)
{
  y = beta_null + x * beta_null
  
  return(y)
}

# Testing the function
# Setting beta as a random number between 0 and 1
beta <- runif(1)

# Setting the intercept/bias to a random number between 0 and 1
beta_null <- runif(1)

# Randomly sampling 100 values between 0 and 1
x <- runif(100)

# Calling the function and plotting results
y <- equation_2.1(x, beta, beta_null)
qplot(x, y) +
  lims(x = c(0,1), y = c(0,1))

# Clearly this results in a line/linear relationship between x and y

# Defining a function to calculate beta
equation_2.6 <- function(x_mat, y)
{
  return(beta)
}

##### Sim 2.3.3 #####
simulation_2.3.3 <- function(n)
{
  I <- cbind(c(1,0), c(0,1))
  blue_means <- mvrnorm(n=10, mu=c(1,0), Sigma = I)
  orange_means <- mvrnorm(n=10, mu=c(0,1), Sigma = I)
  
  # Simulate blue classes
  blue_table <- data.frame(matrix(nrow = ceiling(n/2), ncol = 3))
  names(blue_table) <- c("y", "x1", "x2")
  blue_table$y <- 0
  for (i in 1:nrow(blue_table))
  {
    mk_rand = sample(1:10, 1)
    blue_table[i,2:3] <- mvrnorm(1, blue_means[mk_rand,], I/5)
  }
  
  # Simulate orange_classes
  orange_table <- data.frame(matrix(nrow = floor(n/2), ncol = 3))
  names(orange_table) <- c("y", "x1", "x2")
  orange_table$y <- 1
  for (i in 1:nrow(orange_table))
  {
    mk_rand = sample(1:10, 1)
    orange_table[i,2:3] <- mvrnorm(1, orange_means[mk_rand,], I/5)
  }
  
  sim_results <- rbind(blue_table, orange_table)
  return(sim_results)
}

example_simulation <- simulation_2.3.3(200)
qplot(example_simulation$x1, example_simulation$x2, color = as.factor(example_simulation$y)) +
  scale_color_manual(values = c('lightskyblue3', 'orange1')) +
  labs(x = "X1", y = "X2", color = "Response")

# Testing a multivariate normal with identity matrix vs just two independent normals
mv_X1 <- mvrnorm(n=5000, mu=c(1,0), Sigma = I)[,1]
n_X1 <- rnorm(5000, 1, 1)
results <- data.frame(c(mv_X1, n_X1))
names(results) <- "mean"
results$dist_type <- NA
results$dist_type[1:5000] <- "multivariate_norm"
results$dist_type[5001:10000] <- "norm"
ggplot(results, aes(x = mean)) +
  geom_histogram(data = subset(results, dist_type == "multivariate_norm"), fill = "red", alpha = 0.2) +
  geom_histogram(data = subset(results, dist_type == "norm"), fill = "blue", alpha = 0.2)

# Implementing least squares equation
equation_2.6 <- function(y_vect, x_mat)
{
  x_null <- rep(1, nrow(x_mat))
  x_mat_full <- as.matrix(cbind(x_null, x_mat))
  betas <- (t(x_mat_full) %*% x_mat_full)^-1 %*% t(x_mat_full) %*% y_vect
  beta_null <- (t(x_null) %*% x_null)^-1 %*% t(x_null) %*% y_vect
  return(betas)
}

predict_eq_2.6 <- function(x_mat, betas)
{
  x_mat_full <- as.matrix(cbind(rep(1, nrow(x_mat)), x_mat))
  
  preds <- x_mat_full %*% betas
  
  return(preds)
}

calc_line_fig2.1 <- function(betas, x_mat)
{
  x1_max <- max(x_mat[,1])
  x1_min <- min(x_mat[,1])
  
  x2_min <- (midpt - betas[1] - x1_min * betas[2]) / betas[3]
  x2_max <- (midpt - betas[1] - x1_max * betas[2]) / betas[3]
  
  boundary_pts <- data.frame(c(x1_min, x1_max), c(x2_min, x2_max))
  names(boundary_pts) <- c("x1", "x2")
  
  return(boundary_pts)
}
midpt <- min(preds) + ((max(preds) - min(preds)) / 2)

test_data <- simulation_2.3.3(200)
lin_model <- lm(y ~ x1 + x2, data = test_data)
lin_preds <- predict(lin_model, test_data[,2:3])
betas <- equation_2.6(test_data$y, as.matrix(test_data[,2:3]))
preds <- predict_eq_2.6(test_data[,2:3], betas)
preds_shift <- (preds + 0 - min(preds))
preds_adj <- preds / max(preds)

qplot(test_data$x1, test_data$x2, color = test_data$y)
qplot(test_data$x1, test_data$x2, color = lin_preds)
test <- x_mat_full %*% betas
