/*
  File:             stdafx.h
  Created by:       Oleksii Pokotylo
  First published:  28.02.2013
  Last revised:     28.02.2013
  
  Defines the Includes needed.
*/

#pragma once

#define BOOST_UBLAS_NO_STD_CERR

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
#include <boost/random/linear_congruential.hpp>
#include <boost/random.hpp>

using namespace std;

#include "DataStructures.h"
#include "Common.h"
#include "AlphaProcedure.h"
#include "TukeyDepth.h"
#include "HD.h"
#include "ZonoidDepth.h"
#include "SimplicialDepth.h"
#include "OjaDepth.h"
#include "Knn.h"
#include "Polynomial.h"
#include "ProjectionDepth.h"


static boost::random::rand48 rEngine;
static boost::random::normal_distribution<double> normDist;

#define ran(x) rEngine()%x
#define setseed(x) rEngine.seed(x)

//#ifndef _MSC_VER
////#include <Rcpp.h> 
////using namespace Rcpp;
//#include <Rmath.h> 
//
////#include <Rcpp/stats/random/runif.h>
////static Rcpp::stats::UnifGenerator rndm;
//#define ran(x) ((int)runif(0,x))
//
//
//#define setseed(x) set_seed(x,x)
//
//#include <Rcpp.h> 
//using namespace Rcpp;
//
//#else
//#define ran(x) rand()%x
//
//#define setseed(x) srand(x)
//
//#endif

int random(int x);

