\name{RFoptions}
\alias{RFoptions}
\alias{Shuffle}
\alias{Shuffle256}
\alias{TwoBit}
\alias{Packed}
\alias{Packed256}
\alias{Multiply}
\alias{Multiply256}
\alias{ThreeBit}
\alias{Hamming2}
\alias{Hamming3}
\alias{AutoCoding}
\alias{NoSNPcoding}
\alias{NoSNPcodingAVX}
\alias{NoSNPcodingR}
\alias{FirstGenuineMethod}
\alias{FirstMoBPSmethod}
\alias{LastGenuineMethod}
\alias{LastMoBPSmethod}
\title{Setting control arguments}
\alias{SNPCODING_NAMES}
\description{
 \command{\link{RFoptions}} sets and returns control arguments for
 diverse packages (\pkg{miraculix}, \pkg{RandomFields}).

 \command{\link{RFoptions}} should not be used within parallelizing R
 commands such as \command{\link[parallel]{mclapply}}
}
%\usage{
% RFoptions(..., no.readonly = TRUE)
%}
%\arguments{
% \item{...}{arguments in \code{tag = value} form, or a list of tagged values.}
% \item{no.readonly}{If \command{\link{RFoptions}} is called without
%   argument then all arguments are returned in a list. If
%   \code{no.readonly=TRUE} then only rewritable arguments are returned.
%   Currenlty all arguments are rewritable. So the list is empty.
% }
%}
\details{
%  Note that this function can only be used after
%  \code{library(RandomFieldsUtils)} has been made available explicitely.

  The specific parameters for \pkg{miraculix} are the following.
  See \code{\link[RandomFieldsUtils]{RFoptions}} in
  \pkg{RandomFieldsUtils}
  for further options.
  
  \describe{
    \item{any2bit}{ logical.
      If \code{TRUE} then always the most time efficient code
      is used among
      \itemize{
	\item{\code{TwoBit} (no SIMD needed)}
	\item{\code{Packed} (SSE2 needed)}
	\item{\code{Shuffle} (SSSE3 needed)}
	\item{\code{Shuffle256} (AVX2 needed)}
      }
      whatever is available.
      
      Default : \code{FALSE}. This value might change to \code{TRUE} in
      future. 
    }
      
    \item{centered}{logical or numerical.
      If \code{TRUE} the \eqn{P} matrix is substracted before the
      crossproduct of the the SNP matrix is calculated,
      see \code{\link{relationshipMatrix}} for the \eqn{P} matrix.

      If numeric, then the length of this vector must equal the number
      of SNPs per individual. Then this vector is substracted for each
      individual. 
      Furthermore, \code{normalized} is FALSE.
      As the size of \code{centered} can be large, this vector is never
      returned by \code{RFoption()}; instead \code{NA} is returned.

      Note that \code{centered} also sets the value of \code{normalized}.
      
      Default : \code{TRUE}
    }
   

    \item{\code{cores}}{
      Number of cores for multicore algorithms.
       
%      Default : 
    }

   \item{digits}{OBSOLETE. scalar.
      If \code{digits} is negative no rounding is performed.
      Else the matrix \eqn{P} when calculating the relationsship matrix
      \eqn{(M-P)^T(M-P)} is rounded to the given number of absolute (not
      significant) digits.
      
      Default : \code{3.0}.
    }
     
 
    \item{normalized}{logical.
      If \code{TRUE} the relationship matrix is normalized by
      \eqn{\sigma^2}, see \code{\link{relationshipMatrix}}.
      

      Its value is set to the value of \code{centered} whenever
      the value of \code{centered} is changed. So \code{normalized}
      must be set always after \code{centered},
      e.g. \code{RFoptions(centered=TRUE, normalized=FALSE)}, but not
      \code{RFoptions(normalized=FALSE, centered=TRUE)}.
      
      Default : \code{TRUE}
    }

  
     \item{snpcoding}{integer.
      Possible values are
      \describe{
	\item{\code{Shuffle}}{two bit mini hash table based on SSSE3}
	\item{\code{Shuffle256}}{two bit mini hash table based on AVX2}
	\item{\code{Packed}}{4-bit integer arithmetic based on SSE2}
	\item{\code{Packed256}}{4-bit integer arithmetic based on AVX2}
	\item{\code{Multiply}}{16-bit integer arithmetic based on SSE2}
	\item{\code{Multiply256}}{16-bit integer arithmetic based on AVX2}
	\item{\code{Twobit}}{two bit hash table}
	\item{\code{Threebit}}{three bit hash table}
	\item{\code{Hamming2}}{method used in PLINK}
 	\item{\code{Hamming3}}{method used in PLINK}
	\item{\code{AutoCoding}}{method is chosen by the programme itself}
	\item{\code{NoSNPcoding}}{no coding, i.e. 32 bit integer}
	\item{\code{NoSNPcodingR}}{No coding: 32 bit integer, R
	  code. Only for testing purposes.}
	\item{\code{NoSNPcodingAVX}}{No coding: 
	  AVX implementation if available (double precision or integer).}
      }

      In \code{for} loops that run through all available methods
      the constants \code{FirstGenuineMethod} and
      \code{LastGenuineMethod}
      might be useful.

      In case of the package \pkg{MoPBS} or if interest is in the 2 bit methods
      only, use the constants \code{FirstMoBPSmethod} and
      \code{LastMoBPSmethod}.

      In case the names of the method is needed, use 
      \code{SNPCODING_NAMES[snp_coding + 1]}.
      
      
      Default : \code{Shuffle}
    }

    
    \item{returnsigma}{
      logical. Whether \eqn{\sigma^2} shall be also returned when
      the relationship matrix is calculated.
    }
  }  
}


 
\value{
 \code{NULL} if any argument is given, and the full list of
 arguments, otherwise.
}

%\references{}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de/de/publications/software}}

\seealso{
\code{\link[RandomFieldsUtils]{RFoptions}},
}


\examples{% library(adoption)
RFoptions()$genetics
}

\keyword{models}

