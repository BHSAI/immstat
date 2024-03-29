"._a2r_hclu"       <- NULL # to receive an hclust object when 
                           # A2Rplot.hclust is called

"._a2r_counter"       <- NA # a counter used in A2Rplot.hclust
"._a2r_height_cut"    <- NA

"._a2r_envir"         <- NA
"._a2r_group"         <- NA


#===============================================================================
"A2Rplot" <- function(x,...){
  UseMethod("A2Rplot")
}
#===============================================================================
"A2Rplot.default" <- function(x,...){
  plot(x,...)
}
#===============================================================================
"A2Rplot.hclust" <- function(
  x ,             # an hclust object to draw
  k        = 2,   # the number of groups
  col.up   = "black",
  col.down = rainbow(k),
  lty.up   = 2,
  lty.down = 1,
  lwd.up   = 1,
  lwd.down = 2,
  type     = c("rectangle","triangle"),
  knot.pos = c("mean","bary","left","right","random"),
  criteria,
  fact.sup,
  show.labels=TRUE,
  only.tree=FALSE,
  main     = paste("Colored Dendrogram (",k," groups)"),
  boxes    = TRUE,
  members,
  ...
){

  if(missing(members)) members <- NULL
  opar <- par(no.readonly=TRUE)
  knot.pos <- match.arg(knot.pos)
  type     <- match.arg(type)
  # tests
  if(k<2) 
    stop("k must be at least 2")  
    
  ._a2r_counter    <<- 0
  ._a2r_hclu       <<- x

  ._a2r_envir      <<- environment()
  nn <- length(x$order) - 1

  ._a2r_height_cut <<- mean(x$height[nn-k+1:2])
  ._a2r_group      <<- 0
  
  n.indiv   <- length(x$order)
  groups.o  <- cutree.order(x, k=k)[x$order]
  
  bottom <- if(is.null(members)) 0 else x$height[nn] * -.2 
  
  if(only.tree){
    if(is.null(members)) plot(0,type="n",xlim=c(0.5,n.indiv+.5), ylim=c(bottom,x$height[nn]), xaxs="i", axes=FALSE, xlab="",ylab="") 
    else                 plot(0,type="n",xlim=c(0.5,sum(members)+.5), ylim=c(bottom,x$height[nn]), xaxs="i", axes=FALSE, xlab="",ylab="")
    #call to the ** recursive function ** .rec.hclust
    .rec.hclust(nn, col=col.up, lty=lty.up, lwd=lwd.up)
    
    if(boxes){
      axis(2)
      box()
    }
    return(NULL)
  }
  
  # prepare the layout
  matlayout <- matrix(c(2,4,6,1,3,5), nc=2, nr=3)
  widths    <- c(1,9)
  heights   <- c(8,1,1)
  if(!show.labels){
      matlayout <- matrix(c(2,4,1,3), nc=2, nr=2)
      widths    <- c(1,9)
      heights   <- c(9,1)
  }
  if(!missing(fact.sup) ) {
    heights   <- c(8,1,1)
  }
  if(missing(criteria) & missing(fact.sup)){
    matlayout <- matrix(c(2,4,1,3), nc=2, nr=2)
      widths    <- c(1,9)
      heights   <- c(9,1)
    
  }
  heights[1]=heights[1]-6
  layout(matlayout, width=widths, height=heights)
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ The tree (1)
  par(mar=c(0,0,3,4))
  if(is.null(members)) plot(0,type="n",xlim=c(0.5,n.indiv+.5), ylim=c(bottom,x$height[nn]), xaxs="i", axes=FALSE, xlab="",ylab="") 
  else plot(0,type="n",xlim=c(0.5,sum(members)+.5), ylim=c(bottom,x$height[nn]), xaxs="i", axes=FALSE, xlab="",ylab="") 
  #call to the ** recursive function ** .rec.hclust
  .rec.hclust(nn, col=col.up, lty=lty.up, lwd=lwd.up)
  title(main)
  if(boxes){
    box()
    axis(4)
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Criteria (2)
  if(!missing(criteria)){
    par(mar=c(0,0,3,0))
    plot(0,
         type="n",
         xlim=range(criteria), 
         ylim=c(0,x$height[nn]), 
         axes=FALSE, 
         xlab="",
         ylab="")
    par(las=2)
    n.crit <- length(criteria)
    heights.cut <- ( tail(x$height,n.crit) + 
                     tail(x$height,n.crit+1)[-(n.crit+1)] ) / 2
    heights.cut <- rev(heights.cut)
                   
    points(criteria   , heights.cut   , pch=21, bg="red", type="o")
    points(criteria[k-1], heights.cut[k-1], pch=21, cex=2, bg="blue", xpd=NA)
    if(boxes){
      axis(3)
      box()
    }
  }
  else{
    par(mar=c(0,0,3,0))
    plot(0,
         type="n",
         xlim=c(0,1), 
         ylim=c(0,1), 
         axes=FALSE, 
         xlab="",
         ylab="")
  }

  if(show.labels){
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Name of the observations (3)
    par(mar=c(0,0,0,4))
    par(srt=90)
    #obs.labels <- toupper(substr(x$labels[x$order],1,6))
	obs.labels <- toupper(substr(x$labels[x$order],1,22))
	output_list = x$labels[x$order]
	write.table(output_list, file= "dendro_xname.csv",  append=FALSE, sep =';', row.names=FALSE)
    if(is.null(members)) {
      plot(0,type="n",xlim=c(0.5,n.indiv+.5), ylim=c(0,1), xaxs="i", axes=FALSE, xlab="",ylab="") 
      text(0.5:(n.indiv-0.5)            , 0, obs.labels, pos=4, col=col.down[groups.o])
    }
    else{
      plot(0,type="n",xlim=c(0.5,sum(members)+.5), ylim=c(0,1), xaxs="i", axes=FALSE, xlab="",ylab="") 
      xo <-   members[x$order]
      text(cumsum(xo)-xo/2, 0, obs.labels, pos=4, col=col.down[groups.o])
    }
    par(srt=0)
    if(boxes){
      box()
    }
  
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Labels (4)
    par(mar=c(0,0,0,0))
    plot(0,type="n",xlim=c(0,1), ylim=c(0,1), xaxs="i", axes=FALSE, xlab="",ylab="") 
    #text(.5,.5,"Labels")
    if(boxes){
      box()
    }
      
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Quali (5,6)
  if(!missing(fact.sup)){
    quali  <- as.factor(fact.sup)[x$order]
    quanti <- as.numeric(quali)

    par(mar=c(1,0,0,4))
    n.levels <- length(levels(quali))
    plot(0,type="n",
         xlim=c(0.5,n.indiv+.5), 
         ylim=c(0,n.levels), 
         xaxs="i", yaxs="i",axes=FALSE, xlab="",ylab="") 
        
    rect(xleft    = (1:n.indiv)-.5,
         xright   = (1:n.indiv)+.5,
         ybottom  = quanti-1, 
         ytop     = quanti,
         col      = col.down[groups.o])
    par(las=1)
    axis(4, (1:n.levels)-.5,levels(quali), tick=FALSE)
      
    if(boxes){
      box()
    }
    
    
    par(mar=c(1,0,0,0))
    plot(0,type="n",xlim=c(0,1), ylim=c(0,1), xaxs="i", axes=FALSE, xlab="",ylab="") 
    text(.5,.5,deparse(substitute(fact.sup)))
    if(boxes){
      box()
    }
  }
  
  
  par(opar) # reset parameter
}

#===============================================================================

".rec.hclust" <- function(
  index, # index of the current tree to draw
  lwd = 1,
  lty = 1,
  col = "black"){

  members <- get('members', envir= ._a2r_envir) 
  bottom  <- get('bottom',  envir= ._a2r_envir) 
  if(index<0){ # it is a leaf
    if(is.null(members)){
       ._a2r_counter <<- ._a2r_counter + 1
       return(list( x = ._a2r_counter,
                    n = 1))       
    }
    else{
      cc <- ._a2r_counter
      mm <- members[-index]
      polygon(x  = c(cc, cc+mm/2, cc+mm),
              y  = c(bottom, 0, bottom),
              col= col, 
              border = col, 
              lwd=lwd)
      ._a2r_counter <<- ._a2r_counter + mm
      return(list(x = cc+mm/2,
                  n = mm))
    }
  }
  
  h.m   <- ._a2r_hclu$height[index]

  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~ do left
  index.l  <- ._a2r_hclu$merge[index,1]
  
  h.l <- if(index.l<0) 0 else ._a2r_hclu$height[index.l]
  if(h.l<._a2r_height_cut & h.m > ._a2r_height_cut){
      ._a2r_group <<- ._a2r_group + 1
      col.l <- get("col.down",envir=._a2r_envir)[._a2r_group]
      lwd.l <- get("lwd.down",envir=._a2r_envir)
      lty.l <- get("lty.down",envir=._a2r_envir)
  }
  else{
      col.l <- col
      lwd.l <- lwd
      lty.l <- lty
  }
  out.l   <- .rec.hclust(index.l, col=col.l, lty=lty.l, lwd=lwd.l)
  x.l     <- out.l$x
  n.l     <- out.l$n
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~ do right
  index.r  <- ._a2r_hclu$merge[index,2]
  h.r <- if(index.r<0) 0 else ._a2r_hclu$height[index.r]
  if(h.r<._a2r_height_cut & h.m > ._a2r_height_cut){
      ._a2r_group <<- ._a2r_group + 1
      col.r <- get("col.down",envir=._a2r_envir)[._a2r_group]
      lwd.r <- get("lwd.down",envir=._a2r_envir)
      lty.r <- get("lty.down",envir=._a2r_envir)
  }
  else{
      col.r <- col
      lwd.r <- lwd
      lty.r <- lty
  }
  out.r   <- .rec.hclust(index.r, col=col.r, lty=lty.r, lwd=lwd.r)
  x.r     <- out.r$x
  n.r     <- out.r$n
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~ draw what you have to draw
  
  type <- get("type",envir=._a2r_envir)
  x.m  <- (x.r + x.l) / 2  
  n    <- n.r + n.l
  x.b  <- (n.r * x.r + n.l * x.l) / n

  
  knot.pos <- get("knot.pos",envir=._a2r_envir) 
  
  x <- switch(knot.pos,
          mean = x.m,
          left = x.l,
          right= x.r,
          random = x.l + runif(1)*(x.r-x.l),
          bary   = x.b)

          
          
  if(type=="rectangle"){
    segments(x0  = c(x.l, x.l, x.r),
             x1  = c(x.l, x.r, x.r),
             y0  = c(h.l, h.m, h.r),
             y1  = c(h.m, h.m, h.m),
             col = col,
             lty = lty,
             lwd = lwd)
  }
  if(type =="triangle"){
    segments(x0  = c(x.l, x.r),
             x1  = c(x  , x),
             y0  = c(h.l, h.r),
             y1  = c(h.m, h.m),
             col = col,
             lty = lty,
             lwd = lwd)
  }
          
          
  list(x=x,n=n)
}
#===============================================================================
"cutree.order" <- function(hclu, k=NULL, h=NULL){
  
  coupe <- cutree(hclu,k=k, h=h)

  coupe.or <- coupe[hclu$order]
  coupe.out<- rep(NA,length(coupe))
  j <-  1 #
  k <-  coupe.or[1]
  for(i in 1:length(coupe)){
    if(coupe.or[i]==k) next
    else{
      coupe.out[which(coupe==k)] <- j
      j <- j + 1
      k <- coupe.or[i]
    }
  }
  coupe.out[is.na(coupe.out)] <- j
  names(coupe.out) <- names(coupe)
  coupe.out
}
#===============================================================================
