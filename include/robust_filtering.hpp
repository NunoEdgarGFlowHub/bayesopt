/**  \file robust_filtering.hpp \brief Robust filtering for outliers */
/*
-------------------------------------------------------------------------
   This file is part of BayesOpt, an efficient C++ library for 
   Bayesian optimization.

   Copyright (C) 2011-2015 Ruben Martinez-Cantin <rmcantin@unizar.es>
 
   BayesOpt is free software: you can redistribute it and/or modify it 
   under the terms of the GNU Affero General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   BayesOpt is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Affero General Public License for more details.

   You should have received a copy of the GNU Affero General Public License
   along with BayesOpt.  If not, see <http://www.gnu.org/licenses/>.
------------------------------------------------------------------------
*/


#ifndef  _ROBUST_FILTERING_HPP_
#define  _ROBUST_FILTERING_HPP_

#include <boost/scoped_ptr.hpp>
#include "bayesopt/parameters.hpp"
#include "specialtypes.hpp"
#include "dataset.hpp"
#include "posteriormodel.hpp"


/**
 * Namespace of the library interface
 */
namespace bayesopt {


  //Forward declaration
  class PosteriorModel;
  class Dataset;

  /** \addtogroup BayesOpt
   *  \brief Filtering for outliers
   */
  /*@{*/

  class RobustFiltering
  {
  public:
    /** 
     * Constructor
     * @param params set of parameters (see parameters.hpp)
     */
    RobustFiltering(size_t dim, Parameters parameters, randEngine& eng);

    /** 
     * Default destructor
     */
    virtual ~RobustFiltering();

    const Dataset* filterPoints();

    void setSamples(const matrixd &x, const vectord &y)
    {
      mRobustModel->setSamples(x,y);
    };
    
    void addSample(const vectord &x, double y)
    {
      mRobustModel->addSample(x,y);
    };

  private:
    boost::scoped_ptr<PosteriorModel> mRobustModel;
    boost::scoped_ptr<Dataset> mFilteredData;

    double up_margin, low_margin;

    RobustFiltering();
  };

  /**@}*/



} //namespace bayesopt


#endif
