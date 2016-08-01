/**
 * @file RcppBySum.cpp
 *
 * This file is part of Cyclops
 *
 * Copyright 2016 Observational Health Data Sciences and Informatics
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef __RcppBySum_cpp__
#define __RcppBySum_cpp__


#include <Rcpp.h>
#include "BySum.h"

using namespace Rcpp;

// [[Rcpp::export(".bySum")]]
DataFrame bySum(List ffValues, List ffBins) {

    using namespace ohdsi::cyclops;

    try {
        std::map<double,double> map = BySum::bySum(ffValues, ffBins);
        std::vector<double> bins;
        std::vector<double> sums;
        for(std::map<double,double>::iterator iter = map.begin(); iter != map.end(); ++iter){
            bins.push_back(iter->first);
            sums.push_back(iter->second);
        }
        return DataFrame::create(_["bins"] = bins, _["sums"] = sums);
    } catch (std::exception &e) {
        forward_exception_to_r(e);
    } catch (...) {
        ::Rf_error("c++ exception (unknown reason)");
    }
    return DataFrame::create();
}

#endif // __RcppBySum_cpp__
