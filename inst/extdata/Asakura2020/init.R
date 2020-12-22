# SOURCE: https://static-content.springer.com/esm/art%3A10.1038%2Fs42003-020-0863-y/MediaObjects/42003_2020_863_MOESM9_ESM.txt
# CITATION: https://www.nature.com/articles/s42003-020-0863-y#Sec16 

#
# miRNA_auto_FISHER
# (C) 2017-2019 DYNACOM., Co., Ltd.
#
# load libraries ####
library(hash)
library(MASS)
library(pROC)

# INPUT / OUTPUT / SETTINGS #####
fisher_pivot_num <- 20
fisher_turn_num <- 3
cancer_type <- "SC"
model_file <- "SC_20190523-160213-model.csv"
test_file <- "SC_20190523-160213-test.csv"
output_dir <- "output"
include_miRNA <- c()
include_miR_logical <- "and"
exclude_miRNA <- c()

# PROGRAM #####

print_log <- function(msg) {
  cat(paste0(msg, "\n"))
}

get_process_name <- function() {
  start_time <- format(Sys.time(), "%Y/%m/%d %H:%M:%OS")
  pname <- gsub("/", "", start_time)
  pname <- gsub(" ", "-", pname)
  pname <- gsub(":", "", pname)
  pname <- paste(pname, sep = "")
  pname
}

get_cols_str <- function(cols) {
  cols <- sort(cols)
  cols_str <- ""
  for (i in seq_along(cols)) {
    if (i > 1) {
      cols_str <- sprintf("%s:", cols_str)
    }
    cols_str <- sprintf("%s%d", cols_str, cols[i])
  }
  return(cols_str)
}

ck_cols_hash <- function(cols, cols_hash) {
  ret <- 1
  cols_str <- get_cols_str(cols)
  if (has.key(cols_str, cols_hash)) {
    ret <- 0
  } else {
    cols_hash[[cols_str]] <- 1
  }
  return(ret)
}

ck_include_miRNA <- function(cnames) {
  ret <- 0
  if (length(include_miRNA) == 0) {
    ret <- 1
  } else {
    i <- 1
    while (i <= length(cnames)) {
      if (length(grep(cnames[i], include_miRNA)) > 0) {
        ret <- ret + 1
      }
      i <- i + 1
    }
    if (length(cnames) >= length(include_miRNA)) {
      if (include_miR_logical == "and" && ret != length(include_miRNA)) {
        ret <- 0
      }
    }
  }
  return(ret)
}

get_fisher_poly_func <- function(a, cnames, b) {
  str <- ""
  for (i in seq_along(cnames)) {
    if (i > 1) {
      str <- paste(str, "+", sep = "")
    }
    str <- paste(str, "(", format(a[i], digit = 6), ")*", cnames[i], sep = "")
  }
  if (b >= 0) {
    str <- paste(str, "-", sep = "")
  }
  else {
    str <- paste(str, "+", sep = "")
  }
  str <- paste(str, format(abs(b), digit = 6), sep = "")
  return(str)
}

calc_spec <- function(idx, cols, learn, test) {
  cols <- cols[1:fisher_turn_num]
  cols <- cols[cols > 0]
  cnames <- colnames(learn)[cols]
  miR_num <- length(cnames)
  learn2 <- data.frame(learn[, c(cols, ncol(learn))])
  colnames(learn2)[ncol(learn2)] <- "Y"
  z <- lda(Y ~ ., learn2, prior = c(0.5, 0.5))
  a <- z$scaling
  b <- apply(z$means %*% z$scaling, 2, mean)

  rev <- 1
  case_s <- learn[, ncol(learn)]
  case_s_test <- test[, ncol(test)]
  pos <- grep(1, case_s)
  pos_test <- grep(1, case_s_test)
  model_predict <- c()
  test_predict <- c()
  model_predict <- as.vector(predict(z)$x)
  if (mean(model_predict[pos]) < 0) {
    a <- -1 * a
    b <- -1 * b
    model_predict <- -1 * model_predict
    rev <- -1
  }

  test2 <- as.data.frame(test[, cols])
  colnames(test2) <- colnames(learn2[-ncol(learn2)])
  test_predict <- as.vector(predict(z, test2)$x)
  test_predict <- rev * test_predict

  poly_func <- get_fisher_poly_func(a, cnames, b)
  Y <- learn[, ncol(learn)]
  mroc <- roc(Y, model_predict,
    plot = F,
    main = main, direction = "<", levels = c(0L, 1L)
  )
  auc <- as.numeric(split(mroc$auc, "curve: ")[1])
  thre_opt <- coords(mroc, "best", ret = c("threshold"), transpose = TRUE)[1]
  coords <- coords(mroc, thre_opt, "threshold",
    ret = c(
      "sensitivity", "specificity",
      "accuracy", "ppv", "npv"
    ),
    transpose = TRUE
  )
  model_sens <- format(coords[1], digits = 4)
  model_spec <- format(coords[2], digits = 4)
  model_accu <- format(coords[3], digits = 4)
  model_ppv <- format(coords[4], digits = 4)
  model_npv <- format(coords[5], digits = 4)
  model_auc <- format(auc, digits = 4)

  Y <- test[, ncol(test)]
  mroc <- roc(Y, test_predict,
    plot = F,
    main = main, direction = "<", levels = c(0L, 1L)
  )
  auc <- as.numeric(split(mroc$auc, "curve: ")[1])
  coords <- coords(mroc, thre_opt, "threshold",
    ret = c(
      "sensitivity", "specificity",
      "accuracy", "ppv", "npv"
    ),
    transpose = TRUE
  )
  test_sens <- format(coords[1], digits = 4)
  test_spec <- format(coords[2], digits = 4)
  test_accu <- format(coords[3], digits = 4)
  test_ppv <- format(coords[4], digits = 4)
  test_npv <- format(coords[5], digits = 4)
  test_auc <- format(auc, digits = 4)

  return(list(
    miR_num = miR_num,
    z = z,
    poly_func = poly_func,
    model_predict = model_predict,
    test_predict = test_predict,
    thre_opt = format(thre_opt, digits = 4),
    model_sens = model_sens,
    model_spec = model_spec,
    model_accu = model_accu,
    model_ppv = model_ppv,
    model_npv = model_npv,
    model_auc = model_auc,
    test_sens = test_sens,
    test_spec = test_spec,
    test_accu = test_accu,
    test_ppv = test_ppv,
    test_npv = test_npv,
    test_auc = test_auc
  ))
}

dot.plot <- function(vp, rp, case, x1, ctrl, x2, main0, thre_opt) {
  vp1 <- vp[rp == 1]
  vp2 <- vp[rp == 0]
  ovp1 <- order(vp1)
  ovp2 <- order(vp2)
  accu <- 0.02
  stp <- 0.01
  n1 <- length(vp1)
  n2 <- length(vp2)
  vp1 <- round(vp1 / accu) * accu
  vp2 <- round(vp2 / accu) * accu
  offset <- 0.5
  freq1 <- table(vp1)
  freq2 <- table(vp2)
  ylim <- c(-10.0, 10.0)
  xlim <- c(0.19, 0.8)
  str_x <- 0.3
  str_y1 <- 0.7
  str_y2 <- 0.3
  x1 <- c()
  y1 <- c()
  x2 <- c()
  y2 <- c()
  i <- 1
  for (j in seq(along = freq1)) {
    sp <- offset - stp * freq1[j] / 2
    for (k in 1:freq1[j]) {
      x1[i] <- sp + (k - 1) * stp - 0.2
      y1[i] <- as.numeric(names(freq1)[j])
      i <- i + 1
    }
  }

  i <- 1
  for (j in seq(along = freq2)) {
    sp <- offset - stp * freq2[j] / 2
    for (k in 1:freq2[j]) {
      x2[i] <- sp + (k - 1) * stp - 0.04
      y2[i] <- as.numeric(names(freq2)[j])
      i <- i + 1
    }
  }

  plot(x1, y1,
    xaxt = "n", xlim = xlim, ylim = ylim, xlab = "", ylab = "Predict Value",
    main = main0, col = "red"
  )
  points(x2, y2, col = rgb(0, 0.6, 0))

  abline(h = thre_opt, col = "blue")
  abline(h = 0.0)
}

output_fisher_file <- function(model, test, turn_n, all.res, process_name) {
  fisher_dat <- data.frame(
    id = 1, miR_num = 1, poly_func = "",
    model_sens = "", model_spec = "", model_accu = "",
    model_ppv = "", model_npv = "", model_auc = "",
    test_sens = "", test_spec = "", test_accu = "",
    test_ppv = "", test_npv = "", test_auc = "", thre = "",
    stringsAsFactors = FALSE
  )[integer(0), ]

  fname0 <- file.path(
    output_dir,
    paste0(cancer_type, "_", process_name, "-fisher")
  )
  csv_fname <- paste0(fname0, turn_n, ".csv")
  pdf_fname <- paste0(fname0, turn_n, ".pdf")
  pdf(file = pdf_fname, family = "Japan1GothicBBB")
  par(mfrow = c(2, 2), mar = c(4.5, 4.5, 4, 1))

  for (i in seq_along(all.res)) {
    fisher_dat[i, "id"] <- i
    fisher_dat[i, "miR_num"] <- all.res[[i]]$miR_num
    fisher_dat[i, "poly_func"] <- all.res[[i]]$poly_func
    fisher_dat[i, "model_sens"] <- all.res[[i]]$model_sens
    fisher_dat[i, "model_spec"] <- all.res[[i]]$model_spec
    fisher_dat[i, "model_accu"] <- all.res[[i]]$model_accu
    fisher_dat[i, "model_ppv"] <- all.res[[i]]$model_ppv
    fisher_dat[i, "model_npv"] <- all.res[[i]]$model_npv
    fisher_dat[i, "model_auc"] <- all.res[[i]]$model_auc
    fisher_dat[i, "test_sens"] <- all.res[[i]]$test_sens
    fisher_dat[i, "test_spec"] <- all.res[[i]]$test_spec
    fisher_dat[i, "test_accu"] <- all.res[[i]]$test_accu
    fisher_dat[i, "test_ppv"] <- all.res[[i]]$test_ppv
    fisher_dat[i, "test_npv"] <- all.res[[i]]$test_npv
    fisher_dat[i, "test_auc"] <- all.res[[i]]$test_auc
    fisher_dat[i, "thre"] <- all.res[[i]]$thre_opt

    z <- all.res[[i]]$z
    main <- paste0("model=", i, ")")
    Y <- model[, ncol(model)]
    predict <- unlist(all.res[[i]]$model_predict)
    main0 <- paste0("DOT PLOT (model=", i, " miR_num=", turn_n, ")")
    dot.plot(predict, Y, "CASE", 0.0, "CTRL", 0.0, main0, all.res[[i]]$thre_opt)
    main <- paste0(cancer_type, " Fisher ROC (model=", i, ")")
    msg <- c()
    msg[1] <- paste("Thre:", all.res[[i]]$thre_opt)
    msg[2] <- paste("Sens: ", all.res[[i]]$model_sens)
    msg[3] <- paste("Spec: ", all.res[[i]]$model_spec)
    msg[4] <- paste("Accu: ", all.res[[i]]$model_accu)
    msg[5] <- paste("ppv: ", all.res[[i]]$model_ppv)
    msg[6] <- paste("npv: ", all.res[[i]]$model_npv)
    legend("bottomright", legend = msg, cex = 0.8)
    mroc <- roc(Y, predict,
      plot = T, main = main,
      print.auc = TRUE, direction = "<", levels = c(0L, 1L)
    )

    Y <- test[, ncol(test)]
    predict <- unlist(all.res[[i]]$test_predict)
    main0 <- paste0("DOT PLOT (test=", i, " miR_num=", turn_n, ")")
    dot.plot(predict, Y, "CASE", 0.0, "CTRL", 0.0, main0, all.res[[i]]$thre_opt)
    main <- paste0(cancer_type, " Fisher ROC (test=", i, ")")
    msg <- c()
    msg[1] <- paste("Thre:", all.res[[i]]$thre_opt)
    msg[2] <- paste("Sens: ", all.res[[i]]$test_sens)
    msg[3] <- paste("Spec: ", all.res[[i]]$test_spec)
    msg[4] <- paste("Accu: ", all.res[[i]]$test_accu)
    msg[5] <- paste("ppv: ", all.res[[i]]$test_ppv)
    msg[6] <- paste("npv: ", all.res[[i]]$test_npv)
    legend("bottomright", legend = msg, cex = 0.8)
    mroc <- roc(Y, predict,
      plot = T, main = main,
      print.auc = TRUE, direction = "<", levels = c(0L, 1L)
    )
  }
  dev.off()
  write.csv(fisher_dat, csv_fname, quote = FALSE, row.names = FALSE)
}

output_turn_fisher <- function(model, test, turn_n, fisher_miR, process_name) {
  all.res <- list()
  for (i in seq_len(nrow(fisher_miR))) {
    if (is.na(fisher_miR[i, 1])) break
    sp <- fisher_miR[i, 1:fisher_turn_num]
    all.res[[i]] <- calc_spec(i, sp, model, test)
  }

  output_fisher_file(model, test, turn_n, all.res, process_name)
}

next_miR <- function(learn, pivot_n, cols, cols_hash) {
  miR_aics <- data.frame(stringsAsFactors = FALSE)
  learn2 <- data.frame()
  for (i in seq_len(ncol(learn) - 1)) {
    ng <- 0
    ck <- c(cols, i)
    if (ck_include_miRNA(colnames(learn)[ck])) {
      if (length(cols) > 0) {
        if (length(unique(ck)) > length(cols) &&
          ck_cols_hash(ck, cols_hash)) {
          learn2 <- data.frame(learn[, c(ck, ncol(learn))])
        } else {
          ng <- 1
        }
      }
      else {
        learn2 <- data.frame(learn[, c(ck, ncol(learn))])
      }
    } else {
      ng <- 1
    }
    if (ng == 0) {
      colnames(learn2)[ncol(learn2)] <- "Y"
      miR_aics[i, 1] <- i
      z <- lda(Y ~ ., learn2, CV = TRUE)
      z.tab <- table(learn2$Y, z$class)
      miR_aics[i, 2] <- sum(z.tab[row(z.tab) == col(z.tab)]) / sum(z.tab)
    }
    if (i %% 200 == 0) {
      if (length(cols) > 0) {
        cols_str <- get_cols_str(cols)
        print_log(sprintf("cols=%s %d miRNA processed.", cols_str, i))
      } else {
        print_log(sprintf("%d miRNA processed.", i))
      }
    }
  }

  return(miR_aics)
}

main <- function() {
  process_name <- get_process_name()
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  logfile <- file.path(output_dir, paste0(process_name, ".log"))
  cat(paste0("miRNA_auto_FISHER start: ", format(Sys.time()), "\n\n"),
      file = logfile, append = TRUE)

  all.params <- list(
    fisher_pivot_num = fisher_pivot_num,
    fisher_turn_num = fisher_turn_num,
    cancer_type = cancer_type,
    model_file = model_file,
    test_file = test_file,
    output_dir = output_dir,
    include_miRNA = include_miRNA,
    include_miR_logical = include_miR_logical,
    exclude_miRNA = exclude_miRNA
  )
  capture.output(all.params, file = logfile, append = TRUE)

  print_log(paste0("Fisher Process start (", process_name, ")"))

  model <- read.csv(model_file,
    header = TRUE,
    row.names = 1, check.names = F, stringsAsFactors = FALSE
  )
  test <- read.csv(test_file,
    header = TRUE,
    row.names = 1, check.names = F, stringsAsFactors = FALSE
  )

  if (!identical(colnames(model), colnames(test))) {
    stop("MODEL/TEST column name mismatched")
  }

  for (i in seq_len(ncol(model))) {
    miR2 <- unlist(strsplit(colnames(model)[i], ","))
    if (length(miR2) > 1) {
      colnames(model)[i] <- miR2[1]
      colnames(test)[i] <- miR2[1]
    }
  }

  for (miR in exclude_miRNA) {
    col <- charmatch(miR, colnames(model), 0)
    if (col > 0) {
      model <- model[-col]
      test <- test[-col]
    }
  }

  fisher_miR <- matrix(0, nrow = fisher_pivot_num, ncol = fisher_turn_num + 1)

  print_log(paste0("Fisher process turn : ", 1))

  cols <- c()
  cols_hash <- hash()
  miRs <- next_miR(model, fisher_pivot_num, cols, cols_hash)
  sortlist <- order(miRs$V2, decreasing = T)
  miRs <- miRs[sortlist, ]

  for (i in seq_len(nrow(fisher_miR))) {
    fisher_miR[i, 1] <- miRs$V1[i]
    fisher_miR[i, ncol(fisher_miR)] <- miRs$V2[i]
  }
  output_turn_fisher(model, test, 1, fisher_miR, process_name)

  for (i in seq(2, fisher_turn_num)) {
    print_log(paste0("Fisher process turn : ", i))
    k <- 0
    wk_fisher_miR <- matrix(0,
      nrow = fisher_pivot_num * ncol(model),
      ncol = ncol(fisher_miR)
    )
    for (j in seq_len(fisher_pivot_num)) {
      if (is.na(fisher_miR[j, 1])) break
      cols <- fisher_miR[j, 1:(i - 1)]
      miRs <- next_miR(model, fisher_pivot_num, cols, cols_hash)
      d <- k * ncol(model)
      for (l in seq_len(nrow(miRs))) {
        wk_fisher_miR[d + l, 1:(i - 1)] <- cols
        wk_fisher_miR[d + l, i] <- miRs$V1[l]
        wk_fisher_miR[d + l, ncol(fisher_miR)] <- miRs$V2[l]
      }
      k <- k + 1
    }
    sortlist <- order(wk_fisher_miR[, ncol(fisher_miR)], decreasing = T)
    wk_fisher_miR <- wk_fisher_miR[sortlist, ]
    wk_fisher_miR <- wk_fisher_miR[wk_fisher_miR[, ncol(fisher_miR)] > 0, ]
    for (j in seq_len(min(nrow(fisher_miR), nrow(wk_fisher_miR)))) {
      fisher_miR[j, ] <- wk_fisher_miR[j, ]
    }
    output_turn_fisher(model, test, i, fisher_miR, process_name)
    clear(cols_hash)
  }

  print_log(paste0("Fisher Process end (", process_name, ")"))

  capture.output(sessionInfo(), file = logfile, append = TRUE)
  cat(paste0("\nmiRNA_auto_FISHER end: ", format(Sys.time()), "\n"),
      file = logfile, append = TRUE)
}

# MAIN #####
# main()