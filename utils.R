library(epitools)

my_riskratio <- function (x, y = NULL, conf.level = 0.95, rev = c("neither", 
                                                                  "rows", "columns", "both"), correction = FALSE, verbose = FALSE) 
  # with continuity correction:
  # https://stats.stackexchange.com/questions/21298/confidence-interval-around-the-ratio-of-two-proportions
  # http://www.zen103156.zen.co.uk/rr.pdf
  # https://stats.stackexchange.com/questions/3112/calculation-of-relative-risk-confidence-interval
{
  if (is.matrix(x) && !is.null(y)) {
    stop("y argument should be NULL")
  }
  if (is.null(y)) {
    x <- epitable(x, rev = rev)
  }
  else {
    x <- epitable(x, y, rev = rev)
  }
  tmx <- table.margins(x)
  p.exposed <- sweep(tmx, 2, tmx["Total", ], "/")
  p.outcome <- sweep(tmx, 1, tmx[, "Total"], "/")
  Z <- qnorm(0.5 * (1 + conf.level))
  nr <- nrow(x)
  small <- matrix(NA, nr, 3)
  small[1, 1] <- 1
  for (i in 2:nr) {
    a0 <- x[1, 2]
    b0 <- x[1, 1]
    a1 <- x[i, 2]
    b1 <- x[i, 1]
    n1 <- a1 + b1
    n0 <- a0 + b0
    m0 <- b0 + b1
    m1 <- a0 + a1

    est <- (a1/n1)/((a0)/(n0))
    logRR <- log(est)
    
    # Delta method
    # choose delta = 0.1 for only a slight correction, avoid divide by 0
    # SElogRR <- sqrt(1/(a1+0.1) - 1/(n1) + 1/(a0+0.1) - 1/(n0))
    # ci <- exp(logRR + c(-1, 1) * Z * SElogRR)
    
    # Poisson model method
    # https://books.google.com/books?id=oDlBZgyx54kC&pg=PP19&lpg=PP19
    
    alpha <- 1 - conf.level
    rr_low <- ifelse(a1 == 0, 0, n0/n1 * a1/(a0+1) / qf(1-alpha/2, 2*(a0+1), 2*a1))
    rr_high <- n0/n1 * (a1+1)/a0 * qf(1-alpha/2, 2*(a1+1), 2*a0)
    ci <- c(rr_low, rr_high)
    
    small[i, ] <- c(est, ci)
  }
  pv <- tab2by2.test(x, correction = correction)
  colnames(small) <- c("estimate", "lower", "upper")
  rownames(small) <- rownames(x)
  cn2 <- paste("risk ratio with", paste(100 * conf.level, "%", 
                                        sep = ""), "C.I.")
  names(dimnames(small)) <- c(names(dimnames(x))[1], cn2)
  rr <- list(x = x, data = tmx, p.exposed = p.exposed, p.outcome = p.outcome, 
             measure = small, conf.level = conf.level, p.value = pv$p.value, 
             correction = pv$correction)
  rrs <- list(data = tmx, measure = small, p.value = pv$p.value, 
              correction = pv$correction)
  attr(rr, "method") <- "small sample-adjusted UMLE & normal approx (Wald) CI"
  attr(rrs, "method") <- "small sample-adjusted UMLE & normal approx (Wald) CI"
  if (verbose == FALSE) {
    rrs
  }
  else rr
}

my_fish <- function(df) {
  res <- df %>%
    unlist() %>%
    matrix(ncol = 2, byrow = TRUE) %>%
    my_riskratio(correction = TRUE)
  
  res$measure[2, 1:3] %>% 
    c(p_value = res$p.value[2,2]) %>% 
    as.matrix() %>% t() %>% data.frame()
}


heat_theme_top <- function(){
  theme(panel.grid = element_blank(),
        strip.text.y = element_text(angle = 0),
        axis.text.x = element_blank(),
        legend.text = element_text(color = "white"),
        legend.title = element_text(color = 'white'),
        plot.margin = unit(c(1.5,0.2,1,0.5), "lines"),
        panel.spacing = unit(2, "points"))
}


heat_theme_bottom <- function(){
  theme(panel.grid = element_blank(),
        legend.title = element_text(),
        strip.text.y = element_text(angle = 0),
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        plot.margin = unit(c(1.5,0.2,0.5,0.5), "lines"),
        panel.spacing = unit(0, "points"))
}
