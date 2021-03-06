#' Functions for listing all possible partitions or just counting the number of them
#' 
#' The function \code{nkpartitions} lists all possible partitions of n objects in to k clusters.
#' 
# #' @usage nkpartitions(n, k, exact = TRUE, print = FALSE)
# #' nkpar(n, k)
#'
#' @param n Number of units/objects.
#' @param k Number of clusters/groups.
#' @param exact Search for partitions with exactly \code{k} or at most \code{k} clusters.
#' @param print Print results as they are found.
#'
#' @return The matrix or number of possible partitions.
#' 
#' @author Chris Andrews
#'
#' @examples
#' n <- 8 # If larger, the number of partitions increases dramatically,
#' # as does if we increase the number of clusters
#' net <- matrix(NA, ncol = n, nrow = n)
#' clu <- rep(1:2, times = c(3, 5))
#' tclu <- table(clu)
#' net[clu == 1, clu == 1] <- rnorm(n = tclu[1] * tclu[1], mean = 0, sd = 1)
#' net[clu == 1, clu == 2] <- rnorm(n = tclu[1] * tclu[2], mean = 4, sd = 1)
#' net[clu == 2, clu == 1] <- rnorm(n = tclu[2] * tclu[1], mean = 0, sd = 1)
#' net[clu == 2, clu == 2] <- rnorm(n = tclu[2] * tclu[2], mean = 0, sd = 1)
#' # Computation of criterion function with the correct partition
#' nkpar(n = n, k = length(tclu)) # Computing the number of partitions
#' all.par <- nkpartitions(n = n, k = length(tclu)) # Forming the partitions
#' all.par <- lapply(apply(all.par, 1, list), function(x) x[[1]])
#' # to make a list out of the matrix
#' res <- critFunC(M = net, clu = clu, approaches = "val",
#'                 blocks = c("nul", "com"), preSpecM = 4)
#' plot(res) # We get the original partition
#' 
#' @keywords cluster
#' 
#' @export

"nkpartitions" <-
function(n, k, exact=TRUE, print=FALSE) {
# n objects
# k subgroups
# exactly k or at most k?
# print results as they are found?
# Author: Chris Andrews <candrews@buffalo.edu>
   if (n != floor(n) | n<=0) stop("n must be positive integer")
   if (k != floor(k) | k<=0) stop("k must be positive integer")
   if (print) {
     printnkp <- function(a) {
       for (j in seq(max(a))) cat("{", seq(along=a)[a==j], "} ")
       cat("\n")
     }
   }

# How many?
   Stirling2nd <- function(n, k) {
     sum((-1)^seq(0,k-1) * choose(k, seq(0,k-1)) * (k-seq(0,k-1))^n) / 
factorial(k)
   }

   rows <- Stirling2nd(n,k)
   if (!exact & k>1) {
     for (i in seq(k-1,1)) {
       rows <- rows + Stirling2nd(n,i)
     }
   }

   if (print) cat("rows =",rows,"\n")

# Allocate space
   theparts <- matrix(NA, nrow=rows, ncol=n)

# begin counting
   howmany <- 0

# all in one group
   a <- rep(1,n)

# does this count?
   if (!exact | (k==1)) {
# increase count, store, and print
     howmany <- howmany + 1
     theparts[howmany,] <- a
     if (print) printnkp(a)
   }

# search for others
   repeat {

# start at high end
     last <- n
     repeat {

# increment it if possible
       if ((a[last] <= max(a[1:(last-1)])) & (a[last] < k)) {
         a[last] <- a[last]+1

# does this count?
         if (!exact | max(a)==k) {
# increase count, store, and print
           howmany <- howmany + 1
           theparts[howmany,] <- a
           if (print) printnkp(a)
         }
# start again at high end.
         break
       }

# otherwise set to 1 and move to a different object
       a[last] <- 1
       if (last>2) {
         last <- last-1
         next
       }

# report the partitions
       return(theparts)
     }
   }
}

