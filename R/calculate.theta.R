#' Calculates species habitat specialization using co-occurrence based metrics
#' @author David Zeleny (partly based on code written by Jason Fridley, Fridley et al. 2007)

#' @export
calculate.theta <- function (input.matrix, species.data = NULL, reps = 10, psample = 5, method = "multiplicative", beals.file = NULL, parallel = F, no.cores = 2, remove.out = F, verbal = F, juicer = F, tcltk = F) 
{
  METHODS <- c('additive', 'multiplicative', 'pairwise.jaccard', 'multi.sorensen', 'beals')
  method.n <- pmatch(method, METHODS)
  if (is.na (method.n)) stop ('invalid method')
  if (method.n == -1) stop ('ambiguous method')
  method <- METHODS[method.n]
  if (!verbal) win.pb <- NULL
  if ( is.na (reps) || reps < 2) {
    if (verbal) tkmessageBox (type = "ok", message = "Number of random subsamples must be integer >= 2") else stop ("Number of random subsamples must be integer >= 2")
    end.end <- F
  } else
    if (is.na (psample) || psample < 2) 
    {
      if (verbal) tkmessageBox (type = "ok", message = "Minimal frequency of species must be integer >= 2") else stop ("Minimal frequency of species must be integer >= 2")
      end.end <- F
    } else
    {
      if (!is.matrix (input.matrix)) input.matrix <- as.matrix (input.matrix)
      if (max (input.matrix) > 1) input.matrix <- ifelse (input.matrix > 0, 1, 0)
      Nplots <- dim (input.matrix)[1]
      plots.per.spp <- colSums (input.matrix)
      select.spp <- plots.per.spp[plots.per.spp >= psample]
      Nspp <- length (select.spp)
      
      if (method == "beals")
      {
        if (is.null (beals.file))
        {
          beals.matrix <- beals.2 (input.matrix, include = T, verbal = verbal)
          if (verbal) win.pb <- winProgressBar (title = "Beals smoothing", label = 'Start', min = 1, max = ncol (input.matrix), initial = 0, width = 300) 
          for (co in seq (1, ncol (input.matrix)))
          {
            if (verbal) setWinProgressBar (win.pb, co + 1, label = paste ("Prepare beals smoothed table: species ", co))
            beals.temp <- beals.matrix[,co][as.logical (input.matrix[,co])]
            stats.temp <- fivenum (beals.temp)
            iqr <- diff (stats.temp [c(2,4)])
            beals.thresh <- min (beals.temp[!(beals.temp < stats.temp[2] - 1.5 * iqr)])
            input.matrix[,co] <- as.numeric (beals.matrix[,co] >= beals.thresh)
          }
          if (verbal) close (win.pb)
          if (tcltk) write.table (input.matrix, file = 'beals-data.txt', sep = '\t', col.names = TRUE)
        } else  
        {
          if (verbal) win.pb <- winProgressBar (title = "Beals smoothing", label = 'Start', min = 1, max = ncol (input.matrix)+1, initial = 0, width = 300) 
          if (verbal) setWinProgressBar (win.pb, 1, label = "Reading Beals smoothed table")
          beals.matrix <- as.matrix (read.delim (file = beals.file, row.names = 1, check.names = F))
          if (! all (dim (beals.matrix) == dim (input.matrix))) {tkmessageBox (type = "ok", message = paste ("Selected Beals matrix has different size than species matrix! \nYou need to calculate new beals smoothed species pool data using current species data. Close the JUICE-R application and run it again from JUICE, and calculate the Multiplicative beta on species pool analysis without selecting the Beals smoothing table.")); stop ('Beals matrix has different size than species matrix!')}
          input.matrix <- beals.matrix
          if (verbal) close (win.pb)
        }
      }
      
      if (!parallel) 
      {
        if (verbal) win.pb <- winProgressBar (title = "Calculation progress", label = paste ("Species no. ", 1), min = 1, max = Nspp, initial = 0, width = 300) 
        temp.res <- lapply (1:Nspp, FUN = function (sp) calculate.theta.0 (input.matrix = input.matrix, sp = sp, select.spp = select.spp, remove.out = remove.out, psample = psample, reps = reps, method = method, parallel = parallel, win.pb = win.pb, verbal = verbal, juicer = juicer))
        if (verbal) close (win.pb)
        theta.out <- as.data.frame (matrix (unlist (temp.res), ncol = 9, byrow = T, dimnames = list (NULL, c('sci.name', 'local.avgS', 'occur.freq', 'meanco', 'meanco.sd', 'meanco.u', 'meanco.l', 'GS', 'GS.sd'))))
        if (!is.null (species.data)) theta.out <- cbind (sci.name = theta.out[,1], species.data[as.character (theta.out[,'sci.name']),1:2], theta.out[,-1])
        if (juicer) write.table (theta.out, file = 'theta_out.txt', sep = '\t', qmethod = 'double', col.names = T, row.names = F) else return (theta.out)
      }
      
      if (parallel)
      {
        require (parallel)
        workers <- makeCluster (no.cores)
        if (verbal) if (file.exists ('GS-progress.txt')) file.remove ('GS-progress.txt')
        clusterExport (workers, c('calculate.theta.0', 'input.matrix', 'select.spp', 'remove.out', 'psample', 'reps', 'method', 'parallel'), envir = environment ())
        temp.res <- clusterApplyLB (workers, 1:Nspp, fun = function (sp) calculate.theta.0 (input.matrix = input.matrix, sp = sp, select.spp = select.spp , remove.out = remove.out, psample = psample, reps = reps, method = method, parallel = parallel, win.pb = NULL, verbal = verbal, juicer = juicer))
        theta.out <- as.data.frame (matrix (unlist (temp.res), ncol = 9, byrow = T, dimnames = list (NULL, c('sci.name', 'local.avgS', 'occur.freq', 'meanco', 'meanco.sd', 'meanco.u', 'meanco.l', 'GS', 'GS.sd'))))
        if (!is.null (species.data)) theta.out <- cbind (sci.name = theta.out[,1], species.data[as.character (theta.out[,'sci.name']),1:2], theta.out[,-1])
        if (juicer) write.table (theta.out, file = 'theta_out.txt', sep = '\t', qmethod = 'double', col.names = T, row.names = F) else return (theta.out)
        stopCluster (workers)
      }
      
      if (juicer) write.table (file = "theta_import.species.data.via.clipboard.txt", theta.out[,c('full.sci.name', 'layer', 'GS')], quote = F, row.names = F, col.names = F, sep = '\t')
      if (juicer) write.table (file = "clipboard", theta.out[,c('full.sci.name', 'layer', 'GS')], quote = F, row.names = F, col.names = F, sep = '\t')
      
      
      if (verbal) tkmessageBox (type = "ok", message = paste ("Species theta values have been copied into clipboard - you can import them directly into JUICE (Edit > Paste Clipboard to BLACK species names).\n\nResult files were saved into", getwd (), "\n\nYou can also use the file theta_import.species.data.via.clipboard.txt to import the species theta values to JUICE (Edit > Paste Clipboard to BLACK species names)."))
      end.end <- T
    }
  cancel <- tclVar (1)
}

#' @name calculate.theta
#' @export
#' 
calculate.theta.0 <- function (input.matrix, species.data, sp, select.spp, remove.out, psample, reps, method, parallel, win.pb, verbal, juicer)
{
  if (verbal) if (parallel) write (paste (sp, '\n'), file = 'GS-progress.txt', append = T) else setWinProgressBar (win.pb, sp, label = paste ("Species no. ", sp))
  temp.matrix <- input.matrix[input.matrix [,colnames (input.matrix) == names (select.spp[sp])]>0,]
  temp.matrix <- temp.matrix[,colSums (temp.matrix) > 0]
  
  if (remove.out)
  {
    veg.dist <- as.matrix (dist (temp.matrix))
    diag (veg.dist) <- NA
    distances <- rowMeans (veg.dist, na.rm = T)
    outliers <- distances > (mean (distances) + 2*sd (distances))
    temp.matrix <- temp.matrix[!outliers,]
    temp.matrix <- temp.matrix[,colSums (temp.matrix) > 0]
  }
  
  if (!nrow (temp.matrix) < psample)  # not sure why there is this condition
  {
    rn.temp.matrix <- matrix (rownames (temp.matrix), ncol = reps, nrow = dim (temp.matrix)[1], byrow = F)
    sample.temp.matrix <- apply (rn.temp.matrix, 2, FUN = function (x) sample (x, psample))
    
    mc.mat <- array(0,dim=c(psample,dim (temp.matrix)[2],reps))  
    for(i in 1:reps) mc.mat[,,i] <- temp.matrix[sample.temp.matrix[,i],]
    total.rich <- colSums (apply (mc.mat, c(2,3), sum) > 0)
    mean.alpha <- colMeans (apply (mc.mat, c(1,3), sum))
    
    if (method == "multiplicative" | method == "beals") Wbeta.vec <- total.rich/mean.alpha
    if (method == "additive") Wbeta.vec <- total.rich-mean.alpha 
    if (method == "pairwise.jaccard") Wbeta.vec <- unlist (lapply (1:reps, FUN = function (i) mean (betapart::beta.pair (mc.mat[,,i], index = 'jaccard')$beta.jac)))
    if (method == "multi.sorensen") Wbeta.vec <- unlist (lapply (1:reps, FUN = function (i) betapart::beta.multi (mc.mat[,,i], index = 'sorensen')$beta.SOR))
    
    GS <- mean(Wbeta.vec)      #mean Whittaker beta value for all reps (= THETA: G-S metric)
    GS.sd <- sd(Wbeta.vec)			#s.d. of above
    meanco <- mean(total.rich)			#mean # cooccurrences in "psample" plots
    meanco.sd <- sd(total.rich)		#s.d. of above
    
    sci.name <- labels (select.spp[sp])	#scientific name
    local.avgS <- mean(mean.alpha)				#approximate mean local richness
    occur.freq <- as.vector (select.spp[sp])							#total number of plots
    
    meanco.u <- qnorm(.975,mean=meanco,sd=meanco.sd)			#97.5% confidence limit
    meanco.l <- qnorm(.025,mean=meanco,sd=meanco.sd)			#2.5% confidence limit
    result <- c(sci.name,local.avgS,occur.freq,meanco,meanco.sd,meanco.u,meanco.l,GS,GS.sd)
    return (result)
  }
}

#' @name calculate.theta
#' @export
beals.2 <- function (x, include = TRUE, verbal = FALSE) # method of beals from vegan, for only p/a data and with progress bar
{
  if (verbal) win.pb2 <- winProgressBar (title = 'Beals', label = 'start', min = 1, max = nrow (x)+2, initial = 0, width = 300)
  x <- as.matrix(x)
  x [x > 0] <- 1
  refX <- x
  incSp <- include
  refX <- as.matrix(refX)
  if (verbal) setWinProgressBar (win.pb2, 1, label = 'Crossprod')
  M <- crossprod(refX, refX)
  C <- diag(M)
  if (verbal) setWinProgressBar (win.pb2, 1, label = 'First sweep')
  M <- sweep(M, 2, replace(C, C == 0, 1), "/")
  if (!incSp) for (i in 1:ncol(refX)) M[i, i] <- 0
  S <- rowSums(x)
  b <- x
  for (i in 1:nrow(x)) {
    if (verbal) setWinProgressBar (win.pb2, i+1, label = i)
    b[i, ] <- rowSums(sweep(M, 2, x[i, ], "*"))
  }                       
  SM <- rep(S, ncol(x))
  if (!incSp) SM <- SM - x
  b <- b/replace(SM, SM == 0, 1)
  if (verbal) close (win.pb2)
  b
}

