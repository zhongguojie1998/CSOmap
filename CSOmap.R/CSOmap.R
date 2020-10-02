CSOmap <- function(DataSetName) {
  library(plotly)
  # load data
  TPMpath <- paste0("./data/", DataSetName, "/TPM.txt")
  LRpath <- paste0("./data/", DataSetName, "/LR_pairs.txt")
  Labelpath <- paste0("./data/", DataSetName, "/label.txt")
  TPM <- read.table(TPMpath, header = TRUE, sep = "\t")
  LR <- read.table(LRpath, header = FALSE, sep = "\t")
  labelData <- read.table(Labelpath, header = TRUE, sep = "\t")
  # create output path
  dir.create(paste0("./results/", DataSetName))
  # set variables properly
  genenames <- TPM$X
  TPM <- TPM[,2:dim(TPM)[2]]
  TPM[is.na(TPM)] = 0
  cellnames <- colnames(TPM)
  labels <- labelData$labels[match(cellnames, labelData$cells)]
  labels[is.na(labels)] = "unlabeled"
  standards <- unique(labels)
  labelIx <- match(labels, standards)
  cellCounts <- table(labelIx)
  # find ligands and receptors TPM
  ligandsIndex <- match(c(LR$V1, LR$V2), genenames)
  receptorIndex <- match(c(LR$V2, LR$V1), genenames)
  ligandsTPM <- as.matrix(TPM[ligandsIndex[!is.na(ligandsIndex)&!is.na(receptorIndex)],])
  receptorTPM <- as.matrix(TPM[receptorIndex[!is.na(ligandsIndex)&!is.na(receptorIndex)],])
  LRscores <- LR$V3[!is.na(ligandsIndex)&!is.na(receptorIndex)]
  affinityMat <- t(ligandsTPM) %*% diag(LRscores) %*% receptorTPM
  # discret affinityMat
  for (i in 1:dim(affinityMat)) {
    affinityArray <- affinityMat[i,]
    affinityArray[i] = 0
    affinityArraySorted <- sort(affinityArray, decreasing = TRUE)
    affinityArray[affinityArray <= affinityArraySorted[50]] = 0
    affinityMat[i,] = affinityArray
  }
  # optimization
  coords <- optimization(affinityMat)
  coordsPath <- paste0("./results/", DataSetName, "/coordinates.txt")
  row.names(coords) <- cellnames
  colnames(coords) <- c('x', 'y', 'z')
  # write coords
  write.table(coords, coordsPath, quote = FALSE, sep = "\t")
  # do statistical tests
  dist <- as.matrix(dist(coords))
  counts <- matrix(0, nrow = length(standards), ncol = length(standards))
  colnames(counts) <- standards
  rownames(counts) <- standards
  k = 3 # default top k connections
  topKs <- c()
  diag(dist) <- Inf
  for (i in 1:dim(dist)[1]) {
    distSorted <- sort(dist[i,])
    topKs <- c(topKs, distSorted[3])
  }
  topK <- median(topKs)
  for (i in 1:dim(dist)[1]) {
    connects <- which(dist[i,]<=topK)
    for (j in connects) {
      counts[labelIx[i], labelIx[j]] = counts[labelIx[i], labelIx[j]] + 1
    }
  }
  diag(counts) <- diag(counts) / 2
  # write counts
  countsPath <- paste0("./results/", DataSetName, "/counts.txt")
  write.table(counts, countsPath, quote = TRUE, sep = "\t")
  # calculate stats
  countsN <- (sum(counts) + sum(diag(counts)))/2
  p_value <- copy(counts)
  K = countsN
  for (i in 1:dim(counts)[1]) {
    for (j in 1:dim(counts)[2]) {
      if (i == j) {
        M <- as.numeric(cellCounts[i]) * (as.numeric(cellCounts[i])-1) / 2
      } else {
        M <- as.numeric(cellCounts[i]) * (as.numeric(cellCounts[j]))
      }
      N <- sum(cellCounts) * (sum(cellCounts) - 1) / 2 - M
      p_value[i, j] <- phyper(counts[i, j], M, N, K, lower.tail = FALSE)
    }
  }
  q_value <- matrix(p.adjust(p_value), nrow = length(standards), ncol = length(standards))
  colnames(q_value) <- standards
  rownames(q_value) <- standards
  # write p and q value
  pvaluePath <- paste0("./results/", DataSetName, "/pvalue.txt")
  qvaluePath <- paste0("./results/", DataSetName, "/qvalue.txt")
  write.table(p_value, pvaluePath, quote = TRUE, sep = "\t")
  write.table(q_value, qvaluePath, quote = TRUE, sep = "\t")
  result = list()
  result$coords = coords
  result$counts = counts
  result$pvalue = p_value
  result$qvalue = q_value
  # # plot 3d and heatmap
  # plot_ly(x=coords[,1], y=coords[,2], z=coords[,3], type="scatter3d", mode="markers", color=labels)
  # heatmap(p_value, scale = "none", Colv = NA, Rowv = NA)
  # heatmap(q_value, scale = "none", Colv = NA, Rowv = NA)
  result
}

optimization <- function (affinityMat, initial_config = NULL, k = 3,  
                          max_iter = 1000, min_cost = 0, 
                          condition = "tight") 
{ # this function is inspired from tsne algorithm. We use similar gradient descent
  # method to optimize our target function specified in our paper.
  # condition can be loose or tight, we suggest using "loose" condition 
  # for dataset with over 10000 cells
  n = dim(affinityMat)[1]
  momentum = 0.5 # initial momentum
  final_momentum = 0.8 # final momentum
  mom_switch_iter = 250 # value to which momentum is changed
  epsilon = 1000 # initial learning rate
  min_gain = 0.01 # minimum gain for delta-bar-delta
  eps = 2^(-52)
  epoch = 100
  if (!is.null(initial_config) && is.matrix(initial_config)) {
    if (nrow(initial_config) != n | ncol(initial_config) != 
        k) {
      stop("initial_config argument does not match necessary configuration for X")
    }
    ydata = initial_config
  }
  else {
    ydata = matrix(rnorm(k * n), n)
  }
  P = 0.5 * (affinityMat + t(affinityMat))
  P[P < eps] <- eps
  P = P/sum(P)
  grads = matrix(0, nrow(ydata), ncol(ydata))
  incs = matrix(0, nrow(ydata), ncol(ydata))
  gains = matrix(1, nrow(ydata), ncol(ydata))
  for (iter in 1:max_iter) {
    
    sum_ydata = apply(ydata^2, 1, sum)
    d = sum_ydata + sweep(-2 * ydata %*% t(ydata), 2, -t(sum_ydata))
    num = 1/(1 + d)
    diag(num) = 0
    Q = num/sum(num)
    if (any(is.nan(num))) 
      message("NaN in grad. descent")
    Q[Q < eps] = eps
    P_Q = P - Q
    P_Q[P_Q > 0 & d <= 0.01] = -0.01;
    stiffnesses = 4 * (P - Q) * num
    for (i in 1:n) {
      grads[i, ] = apply(sweep(-ydata, 2, -ydata[i, ]) * 
                           stiffnesses[, i], 2, sum)
    }
    gains = ((gains + 0.2) * abs(sign(grads) != sign(incs)) + 
               gains * 0.8 * abs(sign(grads) == sign(incs)))
    gains[gains < min_gain] = min_gain
    incs = momentum * incs - epsilon * (gains * grads)
    ydata = ydata + incs
    ydata = sweep(ydata, 2, apply(ydata, 2, mean))
    if (iter == mom_switch_iter) 
      momentum = final_momentum
    if (iter%%epoch == 0) {
      cost = sum(apply(P * log((P + eps)/(Q + eps)), 1, 
                       sum))
      message("Iteration #", iter, " loss function cost is: ", 
              cost)
      if (cost < min_cost) 
        break
    }
    range = max(abs(ydata))
    if (condition == "tight") {
      if (range > 50 && iter%%10 == 0) {
        ydata = ydata * 50/range
      }
    } else {
      if (range > 50 && iter%%max_iter == 0) {
        ydata = ydata * 50/range
      }
    }
  }
  ydata
}
