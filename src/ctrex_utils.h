// header file (.h) for ctrex_utils
// File stores C++ fucntions for direct use in R.

#ifndef ctrex_utils_H
#define ctrex_utils_H

// ----------------------------------------------------------------------------
#include <string.h>
#include "Rcpp.h"
// ----------------------------------------------------------------------------


//' @title
//' Hello ctrex from your C++ friends :).
//'
//' @name
//' hello_ctrex_from_Cpp
//'
//' @description
//' This function returns a greeting message.
//'
//' @return
//' A string with "Hello, ctrex!".
//'
//' @export
// // ------------------------------------------------------------------------
// [[Rcpp::export]]
std::string hello_ctrex_from_Cpp();




// ----------------------------------------------------------------------------
#endif
