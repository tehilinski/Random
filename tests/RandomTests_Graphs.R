#-------------------------------------------------------------------------------
# file: RandomTests_Graphs.R
# Display graphs of results produced by Test_Random.cpp
# Usage:
#   1. Generate random number files:
#      ./Test_Random
#   2. Plot random number graphs:
#      R -q
#      source("RandomTests_Graphs.R")
#
# Author: Thomas E. Hilinski <https:#github.com/tehilinski>
#
# Copyright 2020 Thomas E. Hilinski. All rights reserved.
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#    http:#www.apache.org/licenses/LICENSE-2.0
#    and in the accompanying file LICENSE.md.
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#-------------------------------------------------------------------------------

IUniform <- 1; INormal <- 2; IPoisson <- 3;
fileName <- c( "Random_Uniform.csv", "Random_Normal.csv", "Random_Poisson.csv" )
title <- c( "Random Uniform Distribution", "Random Normal Distribution", "Random Poisson Distribution" )

# For Poisson distribution
mu <- 2.0    # matches value in Test_Random.cpp


ReadCSV <- function( fileName )
{
    message( "Reading file: ", fileName )
    return( read.csv( fileName ) )
}

GraphTitle <- function ( title, count )
{
    return( paste( title, "\ncount =", count ) )
}

FitCurve.Normal <- function ( x, col )
{
    xfit <- seq( min(x), max(x), length=40 )
    yfit <- dnorm( xfit, mean=mean(x), sd=sd(x) )
    lines( xfit, yfit, col=col )
}

FitCurve.Poisson <- function ( x, freqs, mu, col )
{
    xfit <- sort( unique(x) )
    yfit <- dpois( xfit, lambda = mu )
    yscale = max( freqs ) / max( yfit )
    yfit <- yfit * yscale
    lines( xfit, yfit, col=col )
    print( xfit )
    print( yfit )
}

# read data
uniform.df <- ReadCSV( fileName[IUniform] )
normal.df <- ReadCSV( fileName[INormal] )
poisson.df <- ReadCSV( fileName[IPoisson] )

# configure plots
par( mfrow=c(3, 2) )
xlab.lagPlot <- "N - 1"
ylab.lagPlot <- "N"
xlab.hist <- "N"
color.contrast <- "red"

# plots: uniform

hist( uniform.df$N, breaks=40, freq=F,
      main=GraphTitle( title[IUniform], nrow(uniform.df) ), xlab=xlab.hist )
abline( 1, 0, col=color.contrast )

plot( uniform.df$Nprev, uniform.df$N, pch='.', col=color.contrast,
      main=GraphTitle( title[IUniform], nrow(uniform.df) ), xlab=xlab.lagPlot, ylab=ylab.lagPlot )

# plots: normal

hn <- hist( normal.df$N, breaks=40, freq=F,
      main=GraphTitle( title[INormal], nrow(normal.df) ), xlab=xlab.hist )
FitCurve.Normal( normal.df$N, color.contrast )

plot( normal.df$Nprev, normal.df$N, pch='.', col=color.contrast,
      main=GraphTitle( title[INormal], nrow(normal.df) ), xlab=xlab.lagPlot, ylab=ylab.lagPlot )

# plots: poisson

hist.poisson <- hist( poisson.df$N, breaks=40, freq=F,
		      main=GraphTitle( title[IPoisson], nrow(poisson.df) ), xlab=xlab.hist )
FitCurve.Poisson( poisson.df$N, hist.poisson$density, mu, color.contrast )

plot( poisson.df$Nprev, poisson.df$N, pch='+', col=color.contrast,
      main=GraphTitle( title[IPoisson], nrow(poisson.df) ), xlab=xlab.lagPlot, ylab=ylab.lagPlot )
