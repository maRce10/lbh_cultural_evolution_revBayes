## --------------------------------------------------------------------------------------------

## ----load-packages, message=FALSE, results='hide'--------------------------------------------
x <-
  c(
    "tidyr",
    "ggplot2",
    "viridis",
    "coda"
    )

lapply(x, function(y) {
  # check if installed, if not install
  if (!y %in% installed.packages()[, "Package"])
    install.packages(y)
  
  # load package
  try(require(y, character.only = T), silent = T)
})


## ----useful-vars, eval=TRUE, echo=TRUE-------------------------------------------------------
# directory of output files
out_fp <- "../output/most_recent_revbayes_models/"

# all the different factors and their levels
lek <- c("BR1", "CCE", "HC1", "SUR", "TR1")
aln <- c("all.equal", "optimal", "prank")
record <- c("old", "new", "no.fossils")
sampling <- c("all", "early")
clock <- c("global", "Uexp")

# colours and labels for figures below
aln_cols <- viridis(3,direction = -1)
aln_labs <- c("MAFFT-agnostic", "MAFFT-optimal", "PRANK-TN93")

clock_cols <- viridis(2,direction = -1, begin = 0, end = 0.4)
clock_labs <- c("fixed", "relaxed")

record_cols <- viridis(3,direction = 1, begin = 0.7, end = 1)
record_labs <- c("old gapped", "recent complete", "no historical")

sampling_cols <- viridis(3,direction = 1, begin = 0.5, end = 0.8)
sampling_labs <- c("all historical", "earliest record")


## ----Uexp-all-old-data-----------------------------------------------------------------------
file.ls <- list.files(path=out_fp,pattern=paste0(glob2rx("*all_*Uexp*"),"_[0-9]{6}.log"))

dat_com <- data.frame()
for (i in lek){
  if (i == "BR1" | i == "TR1" )
  for (j in aln) {
    tmp_file <- grep(pattern = paste(i,j,"new", sep = "_"), file.ls)
    tmp_data <- read.table(paste0(out_fp,file.ls[tmp_file]), header = T)
    tmp_data$lek <- i
    tmp_data$aln <- j
    br_rates <- grep("branch_rates.", colnames(tmp_data))
    tmp_data <- tmp_data[,-br_rates]
    dat_com <- rbind(dat_com, tmp_data)
  } else {
    for (j in aln) {
    tmp_file <- grep(pattern = paste(i,j,"old", sep = "_"), file.ls)
    tmp_data <- read.table(paste0(out_fp,file.ls[tmp_file]), header = T)
    tmp_data$lek <- i
    tmp_data$aln <- j
    br_rates <- grep("branch_rates.", colnames(tmp_data))
    tmp_data <- tmp_data[,-br_rates]
    dat_com <- rbind(dat_com, tmp_data)
    }
  }
}


## ----plot-div-aln-h--------------------------------------------------------------------------
# what are the parameters of interest?
div_pars <- c("speciation_rate", "extinction_rate", "age_extant", "origin_time") 

# subset the data
plot_dat <- dat_com[,c("Iteration","lek", "aln", div_pars)]

# convert to long
plot_dat <- pivot_longer(plot_dat, cols = all_of(div_pars), names_to = "par" )

# reorder parameters
plot_dat$par_s = factor(plot_dat$par, levels = div_pars)

# plot
p_div_h_aln <- ggplot(data = plot_dat, aes(x = value, fill = aln)) + 
  geom_histogram(position = "identity", bins = 50, colour = "white", alpha = 0.6, size =0.05) +
  facet_grid(rows = vars(lek), cols = vars(par_s), scales = "free") + 
  theme_light(base_size = 12) + 
  scale_fill_manual(values = aln_cols, labels = aln_labs) +
  labs(x = "Parameter value", y = "Posterior frequency", fill = "Alignment")
p_div_h_aln


## ----aln-div-sensitivity, eval=TRUE, echo=FALSE----------------------------------------------
# create a data frame with n rows
n <-  length(lek)*length(aln)*length(div_pars)*2
overlap_df <-
  data.frame(
    lek = character(length = n)
  )

# fill it in with contrasts between query and target alignments
m <- 1
for (i in lek) {
  for (j in aln) {
    for (k in div_pars) {
      tmp_data <-
        dat_com[which(dat_com$lek == i & dat_com$aln == j), k]
      for (h in aln) {
        if (j != h) {
      tmp_query <-
        dat_com[which(dat_com$lek == i & dat_com$aln == h), k]
      l <- HPDinterval(mcmc(tmp_data))[1]
      u <- HPDinterval(mcmc(tmp_data))[2]
      below_l <- length(which(tmp_query < l))
      above_u <- length(which(tmp_query > u))
      No_overlap <- below_l  + above_u
      overlap_df$lek[m] <- i
      overlap_df$query[m] <- h
      overlap_df$target[m] <- j
      overlap_df$par[m] <- k
      overlap_df$contrast[m] <- paste0(h, " vs ", j)
      overlap_df$non_overlap[m] <- No_overlap / length(tmp_data)
      m = m+1
        }
      }
    }
  }
}

# reorder parameters
overlap_df$par_s = factor(overlap_df$par, levels = div_pars)

# plot non overlaps
p_div_s_aln <-
  ggplot(data = overlap_df, aes(x = target, y = non_overlap, colour = query)) +
  geom_point(size = 2, alpha = 0.7) +
  facet_grid(rows = vars(lek), cols = vars(par_s)) +
  theme_light(base_size = 12) +
  scale_colour_manual(values = aln_cols, labels = aln_labs) +
  scale_x_discrete(labels = aln_labs) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  geom_hline(yintercept = 0.05,
             colour = "grey60",
             lty = 2) +
  labs(x = "Target alignment", y = "Non overlap with 95% HPD interval", colour = "Query alignment")

p_div_s_aln


## ----plot-div-h-clock-prank------------------------------------------------------------------
file.ls <- list.files(path=out_fp,pattern=paste0(glob2rx("*prank*all*_*"),"[0-9]{6}.log"))

dat_prank_clock <- data.frame()
for (i in lek){
  if (i == "BR1" | i == "TR1" )
  for (j in clock) {
    tmp_file <- grep(pattern = paste(i,"prank","new","all",j, sep = "_"), file.ls)
    tmp_data <- read.table(paste0(out_fp,file.ls[tmp_file]), header = T)
    tmp_data <- tmp_data[,c("Iteration",div_pars)]
    tmp_data$lek <- i
    tmp_data$clock <- j
      dat_prank_clock <- rbind(dat_prank_clock, tmp_data)
  } else {
    for (j in clock) {
    tmp_file <- grep(pattern = paste(i,"prank","old","all",j, sep = "_"), file.ls)
    tmp_data <- read.table(paste0(out_fp,file.ls[tmp_file]), header = T)
    tmp_data <- tmp_data[,c("Iteration", div_pars)]
    tmp_data$lek <- i
    tmp_data$clock <- j
    dat_prank_clock <- rbind(dat_prank_clock, tmp_data)
    }
  }
}

plot_dat <- dat_prank_clock[,c("Iteration","lek", "clock", div_pars)]
plot_dat <- pivot_longer(plot_dat, cols = all_of(div_pars), names_to = "par" )

plot_dat$par_s = factor(plot_dat$par, levels = div_pars)

p_div_h_clock_prank <- ggplot(data = plot_dat, aes(x = value, fill = clock)) + 
  geom_histogram(position = "identity", bins = 50, colour = "white", alpha = 0.6, size =0.05) +
  facet_grid(rows = vars(lek), cols = vars(par_s), scales = "free") + 
  theme_light(base_size = 12) + 
  scale_fill_manual(values = clock_cols, labels = clock_labs) +
  labs(x = "Parameter value", y = "Posterior frequency", fill = "Clock model")

p_div_h_clock_prank


## ----plot-div-h-clock-optimal----------------------------------------------------------------
file.ls <- list.files(path=out_fp,pattern=paste0(glob2rx("*optimal*all*_*"),"[0-9]{6}.log"))

dat_optimal_clock <- data.frame()
for (i in lek){
  if (i == "BR1" | i == "TR1" )
  for (j in clock) {
    tmp_file <- grep(pattern = paste(i,"optimal","new","all",j, sep = "_"), file.ls)
    tmp_data <- read.table(paste0(out_fp,file.ls[tmp_file]), header = T)
    tmp_data <- tmp_data[,c("Iteration",div_pars)]
    tmp_data$lek <- i
    tmp_data$clock <- j
      dat_optimal_clock <- rbind(dat_optimal_clock, tmp_data)
  } else {
    for (j in clock) {
    tmp_file <- grep(pattern = paste(i,"optimal","old","all",j, sep = "_"), file.ls)
    tmp_data <- read.table(paste0(out_fp,file.ls[tmp_file]), header = T)
    tmp_data <- tmp_data[,c("Iteration", div_pars)]
    tmp_data$lek <- i
    tmp_data$clock <- j
    dat_optimal_clock <- rbind(dat_optimal_clock, tmp_data)
    }
  }
}

plot_dat <- dat_optimal_clock[,c("Iteration","lek", "clock", div_pars)]
plot_dat <- pivot_longer(plot_dat, cols = all_of(div_pars), names_to = "par" )

plot_dat$par_s = factor(plot_dat$par, levels = div_pars)

p_div_h_clock_optimal <- ggplot(data = plot_dat, aes(x = value, fill = clock)) + 
  geom_histogram(position = "identity", bins = 50, colour = "white", alpha = 0.6, size =0.05) +
  facet_grid(rows = vars(lek), cols = vars(par_s), scales = "free") + 
  theme_light(base_size = 12) + 
  scale_fill_manual(values = clock_cols, labels = clock_labs) +
  labs(x = "Parameter value", y = "Posterior frequency", fill = "Clock model")

p_div_h_clock_optimal


## ----plot-div-h-clock-agnostic---------------------------------------------------------------
file.ls <- list.files(path=out_fp,pattern=paste0(glob2rx("*all.equal*all*_*"),"[0-9]{6}.log"))

dat_all.equal_clock <- data.frame()
for (i in lek){
  if (i == "BR1" | i == "TR1" )
  for (j in clock) {
    tmp_file <- grep(pattern = paste(i,"all.equal","new","all",j, sep = "_"), file.ls)
    tmp_data <- read.table(paste0(out_fp,file.ls[tmp_file]), header = T)
    tmp_data <- tmp_data[,c("Iteration",div_pars)]
    tmp_data$lek <- i
    tmp_data$clock <- j
      dat_all.equal_clock <- rbind(dat_all.equal_clock, tmp_data)
  } else {
    for (j in clock) {
    tmp_file <- grep(pattern = paste(i,"all.equal","old","all",j, sep = "_"), file.ls)
    tmp_data <- read.table(paste0(out_fp,file.ls[tmp_file]), header = T)
    tmp_data <- tmp_data[,c("Iteration", div_pars)]
    tmp_data$lek <- i
    tmp_data$clock <- j
    dat_all.equal_clock <- rbind(dat_all.equal_clock, tmp_data)
    }
  }
}

plot_dat <- dat_all.equal_clock[,c("Iteration","lek", "clock", div_pars)]
plot_dat <- pivot_longer(plot_dat, cols = all_of(div_pars), names_to = "par" )

plot_dat$par_s = factor(plot_dat$par, levels = div_pars)

p_div_h_clock_all.equal <- ggplot(data = plot_dat, aes(x = value, fill = clock)) + 
  geom_histogram(position = "identity", bins = 50, colour = "white", alpha = 0.6, size =0.05) +
  facet_grid(rows = vars(lek), cols = vars(par_s), scales = "free") + 
  theme_light(base_size = 12) + 
  scale_fill_manual(values = clock_cols, labels = clock_labs) +
  labs(x = "Parameter value", y = "Posterior frequency", fill = "Clock model")
p_div_h_clock_all.equal


## ----plot-div-h-record-prank-----------------------------------------------------------------
file.ls <- list.files(path=out_fp,pattern=paste0(glob2rx("*prank*all*Uexp_*"),"[0-9]{6}.log","|",glob2rx("*prank*full.process*Uexp_*"),"[0-9]{6}.log"))

dat_prank_record <- data.frame()
for (i in lek) {
  if (i == "SUR") {
    for (j in record) {
      tmp_file <-
        grep(pattern = paste(i, "prank", j , sep = "_"), file.ls)
      tmp_data <-
        read.table(paste0(out_fp, file.ls[tmp_file]), header = T)
      if (j == "no.fossils") {
        tmp_data$age_extant <- NA
      }
      tmp_data <- tmp_data[, c("Iteration", div_pars)]
      tmp_data$lek <- i
      tmp_data$record <- j
      dat_prank_record <- rbind(dat_prank_record, tmp_data)
    }
  } else {
    if (i == "BR1" | i == "TR1") {
      for (j in record[2:3]) {
        tmp_file <-
          grep(pattern = paste(i, "prank", j , sep = "_"), file.ls)
        tmp_data <-
          read.table(paste0(out_fp, file.ls[tmp_file]), header = T)
        if (j == "no.fossils") {
          tmp_data$age_extant <- NA
        }
        tmp_data <- tmp_data[, c("Iteration", div_pars)]
        tmp_data$lek <- i
        tmp_data$record <- j
        dat_prank_record <- rbind(dat_prank_record, tmp_data)
      }
    } else {
      for (j in record[1:2]) {
        tmp_file <-
          grep(pattern = paste(i, "prank", j , sep = "_"), file.ls)
        tmp_data <-
          read.table(paste0(out_fp, file.ls[tmp_file]), header = T)
        tmp_data <- tmp_data[, c("Iteration", div_pars)]
        tmp_data$lek <- i
        tmp_data$record <- j
        dat_prank_record <- rbind(dat_prank_record, tmp_data)
      }
    }
  }
}

plot_dat <- dat_prank_record[,c("Iteration","lek", "record", div_pars)]
plot_dat <- pivot_longer(plot_dat, cols = all_of(div_pars), names_to = "par" )

plot_dat$par_s = factor(plot_dat$par, levels = div_pars)
plot_dat$record = factor(plot_dat$record, levels = c("old", "new", "no.fossils"))

p_div_h_record_prank <- ggplot(data = plot_dat, aes(x = value, fill = record)) + 
  geom_histogram(position = "identity", bins = 50, colour = "white", alpha = 0.6, size =0.05) +
  facet_grid(rows = vars(lek), cols = vars(par_s), scales = "free") + 
  theme_light(base_size = 12) + 
  scale_fill_manual(values = record_cols, labels = record_labs) +
  labs(x = "Parameter value", y = "Posterior frequency", fill = "Fossil record")

p_div_h_record_prank


## ----plot-div-h-record-optimal---------------------------------------------------------------
file.ls <- list.files(path=out_fp,pattern=paste0(glob2rx("*optimal*all*Uexp_*"),"[0-9]{6}.log","|",glob2rx("*optimal*full.process*Uexp_*"),"[0-9]{6}.log"))

dat_optimal_record <- data.frame()
for (i in lek) {
  if (i == "SUR") {
    for (j in record) {
      tmp_file <-
        grep(pattern = paste(i, "optimal", j , sep = "_"), file.ls)
      tmp_data <-
        read.table(paste0(out_fp, file.ls[tmp_file]), header = T)
      if (j == "no.fossils") {
        tmp_data$age_extant <- NA
      }
      tmp_data <- tmp_data[, c("Iteration", div_pars)]
      tmp_data$lek <- i
      tmp_data$record <- j
      dat_optimal_record <- rbind(dat_optimal_record, tmp_data)
    }
  } else {
    if (i == "BR1" | i == "TR1") {
      for (j in record[2:3]) {
        tmp_file <-
          grep(pattern = paste(i, "optimal", j , sep = "_"), file.ls)
        tmp_data <-
          read.table(paste0(out_fp, file.ls[tmp_file]), header = T)
        if (j == "no.fossils") {
          tmp_data$age_extant <- NA
        }
        tmp_data <- tmp_data[, c("Iteration", div_pars)]
        tmp_data$lek <- i
        tmp_data$record <- j
        dat_optimal_record <- rbind(dat_optimal_record, tmp_data)
      }
    } else {
      for (j in record[1:2]) {
        tmp_file <-
          grep(pattern = paste(i, "optimal", j , sep = "_"), file.ls)
        tmp_data <-
          read.table(paste0(out_fp, file.ls[tmp_file]), header = T)
        tmp_data <- tmp_data[, c("Iteration", div_pars)]
        tmp_data$lek <- i
        tmp_data$record <- j
        dat_optimal_record <- rbind(dat_optimal_record, tmp_data)
      }
    }
  }
}

plot_dat <- dat_optimal_record[,c("Iteration","lek", "record", div_pars)]
plot_dat <- pivot_longer(plot_dat, cols = all_of(div_pars), names_to = "par" )

plot_dat$par_s = factor(plot_dat$par, levels = div_pars)
plot_dat$record = factor(plot_dat$record, levels = c("old", "new", "no.fossils"))

p_div_h_record_optimal <- ggplot(data = plot_dat, aes(x = value, fill = record)) + 
  geom_histogram(position = "identity", bins = 50, colour = "white", alpha = 0.6, size =0.05) +
  facet_grid(rows = vars(lek), cols = vars(par_s), scales = "free") + 
  theme_light(base_size = 12) + 
  scale_fill_manual(values = record_cols, labels = record_labs) +
  labs(x = "Parameter value", y = "Posterior frequency", fill = "Fossil record")

p_div_h_record_optimal


## ----plot-div-h-record-agnostic--------------------------------------------------------------
file.ls <- list.files(path=out_fp,pattern=paste0(glob2rx("*all.equal*all*Uexp_*"),"[0-9]{6}.log","|",glob2rx("*all.equal*full.process*Uexp_*"),"[0-9]{6}.log"))

dat_all.equal_record <- data.frame()
for (i in lek) {
  if (i == "SUR") {
    for (j in record) {
      tmp_file <-
        grep(pattern = paste(i, "all.equal", j , sep = "_"), file.ls)
      tmp_data <-
        read.table(paste0(out_fp, file.ls[tmp_file]), header = T)
      if (j == "no.fossils") {
        tmp_data$age_extant <- NA
      }
      tmp_data <- tmp_data[, c("Iteration", div_pars)]
      tmp_data$lek <- i
      tmp_data$record <- j
      dat_all.equal_record <- rbind(dat_all.equal_record, tmp_data)
    }
  } else {
    if (i == "BR1" | i == "TR1") {
      for (j in record[2:3]) {
        tmp_file <-
          grep(pattern = paste(i, "all.equal", j , sep = "_"), file.ls)
        tmp_data <-
          read.table(paste0(out_fp, file.ls[tmp_file]), header = T)
        if (j == "no.fossils") {
          tmp_data$age_extant <- NA
        }
        tmp_data <- tmp_data[, c("Iteration", div_pars)]
        tmp_data$lek <- i
        tmp_data$record <- j
        dat_all.equal_record <- rbind(dat_all.equal_record, tmp_data)
      }
    } else {
      for (j in record[1:2]) {
        tmp_file <-
          grep(pattern = paste(i, "all.equal", j , sep = "_"), file.ls)
        tmp_data <-
          read.table(paste0(out_fp, file.ls[tmp_file]), header = T)
        tmp_data <- tmp_data[, c("Iteration", div_pars)]
        tmp_data$lek <- i
        tmp_data$record <- j
        dat_all.equal_record <- rbind(dat_all.equal_record, tmp_data)
      }
    }
  }
}

plot_dat <- dat_all.equal_record[,c("Iteration","lek", "record", div_pars)]
plot_dat <- pivot_longer(plot_dat, cols = all_of(div_pars), names_to = "par" )

plot_dat$par_s = factor(plot_dat$par, levels = div_pars)
plot_dat$record = factor(plot_dat$record, levels = c("old", "new", "no.fossils"))

p_div_h_record_all.equal <- ggplot(data = plot_dat, aes(x = value, fill = record)) + 
  geom_histogram(position = "identity", bins = 50, colour = "white", alpha = 0.6, size =0.05) +
  facet_grid(rows = vars(lek), cols = vars(par_s), scales = "free") + 
  theme_light(base_size = 12) + 
  scale_fill_manual(values = record_cols, labels = record_labs) +
  labs(x = "Parameter value", y = "Posterior frequency", fill = "Fossil record")

p_div_h_record_all.equal


## ----plot-div-h-sampling-prank---------------------------------------------------------------
file.ls <- list.files(path=out_fp,pattern=paste0(glob2rx("*prank*Uexp_*"),"[0-9]{6}.log"))

dat_prank_sampling <- data.frame()
for (i in lek) {
  if (i == "BR1" | i == "TR1") {
    for (j in sampling) {
      tmp_file <-
        grep(pattern = paste(i, "prank", "new" , j , sep = "_"), file.ls)
      tmp_data <-
        read.table(paste0(out_fp, file.ls[tmp_file]), header = T)
      tmp_data <- tmp_data[, c("Iteration", div_pars)]
      tmp_data$lek <- i
      tmp_data$sampling <- j
      dat_prank_sampling <- rbind(dat_prank_sampling, tmp_data)
    }
  } else {
    for (j in sampling) {
      tmp_file <-
        grep(pattern = paste(i, "prank", "old", j , sep = "_"), file.ls)
      tmp_data <-
        read.table(paste0(out_fp, file.ls[tmp_file]), header = T)
      tmp_data <- tmp_data[, c("Iteration", div_pars)]
      tmp_data$lek <- i
      tmp_data$sampling <- j
      dat_prank_sampling <- rbind(dat_prank_sampling, tmp_data)
    }
  }
}

plot_dat <- dat_prank_sampling[,c("Iteration","lek", "sampling", div_pars)]
plot_dat <- pivot_longer(plot_dat, cols = all_of(div_pars), names_to = "par" )

plot_dat$par_s = factor(plot_dat$par, levels = div_pars)

p_div_h_sampling_prank <- ggplot(data = plot_dat, aes(x = value, fill = sampling)) + 
  geom_histogram(position = "identity", bins = 50, colour = "white", alpha = 0.6, size =0.05) +
  facet_grid(rows = vars(lek), cols = vars(par_s), scales = "free") + 
  theme_light(base_size = 12) + 
  scale_fill_manual(values = sampling_cols, labels = sampling_labs) +
  labs(x = "Parameter value", y = "Posterior frequency", fill = "Fossil sampling")

p_div_h_sampling_prank


## ----plot-div-h-sampling-optimal-------------------------------------------------------------
file.ls <- list.files(path=out_fp,pattern=paste0(glob2rx("*optimal*Uexp_*"),"[0-9]{6}.log"))

dat_optimal_sampling <- data.frame()
for (i in lek) {
  if (i == "BR1" | i == "TR1") {
    for (j in sampling) {
      tmp_file <-
        grep(pattern = paste(i, "optimal", "new" , j , sep = "_"), file.ls)
      tmp_data <-
        read.table(paste0(out_fp, file.ls[tmp_file]), header = T)
      tmp_data <- tmp_data[, c("Iteration", div_pars)]
      tmp_data$lek <- i
      tmp_data$sampling <- j
      dat_optimal_sampling <- rbind(dat_optimal_sampling, tmp_data)
    }
  } else {
    for (j in sampling) {
      tmp_file <-
        grep(pattern = paste(i, "optimal", "old", j , sep = "_"), file.ls)
      tmp_data <-
        read.table(paste0(out_fp, file.ls[tmp_file]), header = T)
      tmp_data <- tmp_data[, c("Iteration", div_pars)]
      tmp_data$lek <- i
      tmp_data$sampling <- j
      dat_optimal_sampling <- rbind(dat_optimal_sampling, tmp_data)
    }
  }
}

plot_dat <- dat_optimal_sampling[,c("Iteration","lek", "sampling", div_pars)]
plot_dat <- pivot_longer(plot_dat, cols = all_of(div_pars), names_to = "par" )

plot_dat$par_s = factor(plot_dat$par, levels = div_pars)

p_div_h_sampling_optimal <- ggplot(data = plot_dat, aes(x = value, fill = sampling)) + 
  geom_histogram(position = "identity", bins = 50, colour = "white", alpha = 0.6, size =0.05) +
  facet_grid(rows = vars(lek), cols = vars(par_s), scales = "free") + 
  theme_light(base_size = 12) + 
  scale_fill_manual(values = sampling_cols, labels = sampling_labs) +
  labs(x = "Parameter value", y = "Posterior frequency", fill = "Fossil sampling")
p_div_h_sampling_optimal


## ----plot-div-h-sampling-agnostic------------------------------------------------------------
file.ls <- list.files(path=out_fp,pattern=paste0(glob2rx("*all.equal*Uexp_*"),"[0-9]{6}.log"))

dat_all.equal_sampling <- data.frame()
for (i in lek) {
  if (i == "BR1" | i == "TR1") {
    for (j in sampling) {
      tmp_file <-
        grep(pattern = paste(i, "all.equal", "new" , j , sep = "_"), file.ls)
      tmp_data <-
        read.table(paste0(out_fp, file.ls[tmp_file]), header = T)
      tmp_data <- tmp_data[, c("Iteration", div_pars)]
      tmp_data$lek <- i
      tmp_data$sampling <- j
      dat_all.equal_sampling <- rbind(dat_all.equal_sampling, tmp_data)
    }
  } else {
    for (j in sampling) {
      tmp_file <-
        grep(pattern = paste(i, "all.equal", "old", j , sep = "_"), file.ls)
      tmp_data <-
        read.table(paste0(out_fp, file.ls[tmp_file]), header = T)
      tmp_data <- tmp_data[, c("Iteration", div_pars)]
      tmp_data$lek <- i
      tmp_data$sampling <- j
      dat_all.equal_sampling <- rbind(dat_all.equal_sampling, tmp_data)
    }
  }
}

plot_dat <- dat_all.equal_sampling[,c("Iteration","lek", "sampling", div_pars)]
plot_dat <- pivot_longer(plot_dat, cols = all_of(div_pars), names_to = "par" )

plot_dat$par_s = factor(plot_dat$par, levels = div_pars)

p_div_h_sampling_all.equal <- ggplot(data = plot_dat, aes(x = value, fill = sampling)) + 
  geom_histogram(position = "identity", bins = 50, colour = "white", alpha = 0.6, size =0.05) +
  facet_grid(rows = vars(lek), cols = vars(par_s), scales = "free") + 
  theme_light(base_size = 12) + 
  scale_fill_manual(values = sampling_cols, labels = sampling_labs) +
  labs(x = "Parameter value", y = "Posterior frequency", fill = "Fossil sampling")
p_div_h_sampling_all.equal


## ----plot-evol-h-freq-aln--------------------------------------------------------------------
freq_pars <- c("sf.3.", "sf.2.", "sf.5.", "sf.4.", "sf.6.", "sf.1.") 

plot_dat <- dat_com[,c("Iteration","lek", "aln", freq_pars)]
plot_dat <- pivot_longer(plot_dat, cols = all_of(freq_pars), names_to = "par" )

plot_dat$par_s = factor(plot_dat$par, levels = freq_pars)
sf.labs <- c("medium trill", "fast trill", "slow trill", "pure tone", "upward tone", "downward tone")
names(sf.labs) <- freq_pars

p_freq_h_aln <- ggplot(data = plot_dat, aes(x = value, fill = aln)) + 
  geom_histogram(position = "identity", bins = 50, colour = "white", alpha = 0.6, size =0.05) +
  facet_grid(rows = vars(lek), cols = vars(par_s), scales = "free", labeller = labeller(par_s = sf.labs)) + 
  theme_light(base_size = 12) + 
  scale_fill_manual(values = aln_cols)

p_freq_h_aln


## ----plot-evol-h-evolw-aln-------------------------------------------------------------------
evol_pars <- c("er.6.", "er.11.", "er.8.", "er.14.", "er.3.", "er.5.") 

plot_dat <- dat_com[,c("Iteration","lek", "aln", evol_pars)]
plot_dat <- pivot_longer(plot_dat, cols = all_of(evol_pars), names_to = "par" )

plot_dat$par_s = factor(plot_dat$par, levels = evol_pars)
er.labs <- c("m <-> f", "m <-> s", "f <-> s", "p <-> u", "p <-> d", "u <-> d")
names(er.labs) <- evol_pars

p_evolw_h_aln <- ggplot(data = plot_dat, aes(x = value, fill = aln)) + 
  geom_histogram(position = "identity", bins = 50, colour = "white", alpha = 0.6, size =0.05) +
  facet_grid(rows = vars(lek), cols = vars(par_s), scales = "free", labeller = labeller(par_s = er.labs)) + 
  theme_light(base_size = 12) + 
  scale_fill_manual(values = aln_cols)

p_evolw_h_aln


## ----plot-evol-h-evolb-aln-------------------------------------------------------------------
evol_pars <- c("er.1.", "er.2.", "er.4.", "er.7.", "er.9.", "er.10.", "er.12.", "er.14.", "er.15.") 

plot_dat <- dat_com[,c("Iteration","lek", "aln", evol_pars)]
plot_dat <- pivot_longer(plot_dat, cols = all_of(evol_pars), names_to = "par" )

plot_dat$par_s = factor(plot_dat$par, levels = evol_pars)
er.labs <- c("d <-> f", "d <-> m", "d <-> s", "f <-> p", "f <-> u", "m <-> p", "m <-> u", "p <-> s", "s <-> u")
names(er.labs) <- evol_pars

p_evolb_h_aln <- ggplot(data = plot_dat, aes(x = value, fill = aln)) + 
  geom_histogram(position = "identity", bins = 50, colour = "white", alpha = 0.6, size =0.05) +
  facet_grid(rows = vars(lek), cols = vars(par_s), scales = "free", labeller = labeller(par_s = er.labs)) + 
  theme_light(base_size = 12) + 
  scale_fill_manual(values = aln_cols)

p_evolb_h_aln


## ----evol-aln-sensitivity--------------------------------------------------------------------
evol_pars <- colnames(dat_com)[grep(pattern = "[f,r]\\.", colnames(dat_com))]

n <-  length(lek)*length(aln)*length(evol_pars)*2
overlap_df <-
  data.frame(
    lek = character(length = n),
    query = character(length = n),
    target = character(length = n),
    par = character(length = n),
    contrast = character(length = n),
    non_overlap = numeric(length = n)
  )

m <- 1
for (i in lek) {
  for (j in aln) {
    for (k in evol_pars) {
      tmp_data <-
        dat_com[which(dat_com$lek == i & dat_com$aln == j), k]
      for (h in aln) {
        if (j != h) {
      tmp_query <-
        dat_com[which(dat_com$lek == i & dat_com$aln == h), k]
      l <- HPDinterval(mcmc(tmp_data))[1]
      u <- HPDinterval(mcmc(tmp_data))[2]
      below_l <- length(which(tmp_query < l))
      above_u <- length(which(tmp_query > u))
      No_overlap <- below_l  + above_u
      overlap_df$lek[m] <- i
      overlap_df$query[m] <- h
       overlap_df$target[m] <- j
      overlap_df$par[m] <- k
      overlap_df$contrast[m] <- paste0(h, " vs ", j)
      overlap_df$non_overlap[m] <- No_overlap / length(tmp_data)
      m = m+1
        }
      }
    }
  }
}

overlap_wide <- pivot_wider(overlap_df, names_from = par, values_from = non_overlap)

p_evol_s_aln <-
  ggplot(data = overlap_df, aes(x = target, y = non_overlap, colour = query)) +
  geom_point(size = 2) +
  facet_grid(rows = vars(lek), cols = vars(par)) +
  theme_light(base_size = 12) +
  scale_colour_manual(values = aln_cols)

p_evol_s_aln

