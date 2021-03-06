\name{print}
\encoding{latin1}
\alias{print}
\alias{print.rcc}
\alias{print.pls}
\alias{print.spls}
\alias{print.summary}
\alias{print.pca}
\alias{print.spca}

\title{Print Methods for CCA, (s)PLS and Summary objects}

\description{
Produce \code{print} methods for class \code{"rcc"}, 
\code{"pls"}, \code{"spls"} and \code{"summary"}.
}

\usage{
\method{print}{rcc}(x, \ldots)

\method{print}{pls}(x, \ldots)

\method{print}{spls}(x, \ldots)

\method{print}{pca}(x, \ldots)

\method{print}{spca}(x, \ldots)

\method{print}{summary}(x, \ldots)
}

\arguments{
  \item{x}{object of class inheriting from \code{"rcc"}, \code{"pls"}, 
    \code{"spls"}, \code{"pca"}, \code{"spca"} or \code{"summary"}.}
  \item{\ldots}{not used currently.}
}

\details{
\code{print} method for \code{"rcc"}, \code{"pls"} or \code{"spls"} class, 
returns a description of the \code{x} object including: the function used, the regularization 
parameters (if \code{x} of class \code{"rcc"}), 
the (s)PLS algorithm used (if \code{x} of class \code{"pls"} or \code{"spls"}),
the samples size,
the number of variables selected on each of the sPLS components (if \code{x} of class \code{"spls"})
and the available components of the object.

\code{print} method for \code{"summary"} class,
gives the (s)PLS algorithm used (if \code{x} of class \code{"pls"} or \code{"spls"}),
the number of variates considered,
the canonical correlations (if \code{x} of class \code{"rcc"}),
the number of variables selected on each of the sPLS components (if \code{x} of class \code{"spls"})
and the available components for Communalities Analysis,
Redundancy Analysis and Variable Importance in the Projection (VIP).
}

\author{S�bastien D�jean, Ignacio Gonz�lez and Kim-Anh L� Cao.}

\seealso{\code{\link{rcc}}, \code{\link{pls}}, \code{\link{spls}}, \code{\link{vip}}.}

\examples{
## print for objects of class 'rcc'
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
print(nutri.res)

## print for objects of class 'summary'
more <- summary(nutri.res, cutoff = 0.65)
print(more)

## print for objects of class 'pls'
data(linnerud)
X <- linnerud$exercise
Y <- linnerud$physiological
linn.pls <- pls(X, Y)
print(linn.pls)

## print for objects of class 'spls'
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50), 
                      keepY = c(10, 10, 10))
print(toxicity.spls)
}

\keyword{regression}
\keyword{multivariate}
