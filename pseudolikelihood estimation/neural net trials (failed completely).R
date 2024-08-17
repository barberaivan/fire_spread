# Compare GAM vs Neural Network to emulate overlap function

library(tidyverse)
library(mgcv)
library(neuralnet)
library(cito)

target_dir <- file.path("files", "overlaps")

fire_name <- "1999_27j_N"
wlist <- readRDS(file.path(target_dir,
                           paste(fire_name, "-simulations_list2.rds", sep = "")))


df <- as.data.frame(
  cbind(
    wlist$like_sim$par_values,
    y = wlist$like_sim$overlap
  )
)
row.names(df) <- NULL
head(df)


maxx <- apply(df, 2, max)
minn <- apply(df, 2, min)
dfs <- as.data.frame(scale(df, center = minn, scale = maxx - minn))

nn.fit <- dnn(y ~ ., data = dfs, hidden = c(20, 20), lr = tune(0.0001, 0.1),
              loss = "mse") # tarda una banda!

nn.fit <- dnn(y ~ ., data = dfs, hidden = c(20, 20),
              loss = "gaussian")
# TIRA ERROR
plot(nn.fit)



# Con cito no anduvo, vamos con neuralnet ---------------------------------

# Split the data into training and testing set
index <- sample(1:nrow(dfs), round(0.75 * nrow(dfs)))
train_ <- dfs[index,]
test_ <- dfs[-index,]


# Build Neural Network
nn <- neuralnet(y ~ .,
                data = train_, hidden = c(20, 20),
                linear.output = TRUE)
# se toma un ratazo


# Predict on test data
pr.nn <- compute(nn, test_)

# Compute mean squared error
pr.nn_ <- pr.nn$net.result * (max(df$y) - min(df$y)) + min(df$y)
test.r <- (test_$y) * (max(df$y) - min(df$y)) + min(df$y)
MSE.nn <- sum((test.r - pr.nn_) ^ 2) / nrow(test_)

# Plot the neural network
plot(nn)