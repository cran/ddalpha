/*
  File:             stdafx.h
  Created by:       Oleksii Pokotylo
  First published:  28.02.2013
  Last revised:     28.02.2013
  
  Defines the Includes needed.
*/

#pragma once

#include <time.h>
#include <algorithm>
#include <math.h>
#include <float.h>
#include <vector>
#include <set>
#include <stdlib.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/random.hpp>

using namespace std;

#include "DataStructures.h"
#include "AlphaProcedure.h"
#include "TukeyDepth.h"
#include "ZonoidDepth.h"
#include "Knn.h"
#include "Polynomial.h"
#include "ProjectionDepth.h"
#include "Common.h"

#ifndef VStudio
#include <Rcpp.h> 
using namespace Rcpp;

#include <Rcpp/stats/random/runif.h>
static Rcpp::stats::UnifGenerator rndm;
#define ran(x) x*rndm()

#else
#define ran(x) rand()%x

#endif

int random(int x);
