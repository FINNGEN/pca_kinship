\name{aberrant.plot}
\alias{aberrant.plot}
\title{
Scatter plot of summary statistics
}
\description{
Draws a scatter plot of the summary statistics showing outliers.
}
\usage{
aberrant.plot(x)
}

\arguments{
  \item{x}{
output from function aberrant
}
}

\details{
"Normal" individuals are coloured with a gradation from black to grey, with darker colours denoting higher density of individuals. Outliers are coloured with a gradation from orange to red, with darker colours denoting higher posterior probability of being an outlier. 99\% confidence ellipse of the inferred inlier distribution is shown as a dashed grey line.
}

\references{
Celine Bellenguez, Amy Strange, Colin Freeman, Wellcome Trust Case Control Consortium 2, Chris CA Spencer. A robust clustering algorithm for identifying problematic samples in genome-wide association studies. Bioinformatics.
}

\author{
Celine Bellenguez and Chris CA Spencer
}

\examples{
x<-rmvt(1000, sigma=matrix(c(1,0.5,0.5,1),2,2), df=3)
res<-aberrant(x, lambda=30, alpha=1, beta = 20)
aberrant.plot(res)
}