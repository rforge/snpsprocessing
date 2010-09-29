%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  fitCalMaTe.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{fitCalMaTe.matrix}
\alias{fitCalMaTe.matrix}
\alias{fitCalMaTe.matrix}


 \title{Calibrates SNP loci according to the CalMaTe method}

 \description{
  Calibrates SNP loci according to the CalMaTe method.
 }

 \usage{\method{fitCalMaTe}{matrix}(T, references, fB1=1/3, fB2=2/3, maxIter=50, ...)}

 \arguments{
  \item{T}{A 2xI \code{\link[base]{numeric}} \code{\link[base]{matrix}} of allele specific copy numbers (ASCNs),
     where 2 is the number alleles and I is the number of samples.}
  \item{references}{A \code{\link[base]{logical}} or \code{\link[base]{numeric}} \code{\link[base]{vector}} specifying which
     samples should be used as the reference set.}
  \item{fB1,fB2}{Thresholds for calling genotypes AA, AB, BB from the
     allele B fractions.}
  \item{maxIter}{The maximum number of iterations without converging
     before the algorithm quits.}
  \item{...}{Not used.}
 }

 \value{
   Returns a 2xI \code{\link[base]{numeric}} \code{\link[base]{matrix}} of calibrated ASCNs.
 }


\keyword{methods}
\keyword{internal}