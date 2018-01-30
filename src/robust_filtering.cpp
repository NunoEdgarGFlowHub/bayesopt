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

#include "robust_filtering.hpp"

namespace bayesopt
{
  RobustFiltering::RobustFiltering(size_t dim, Parameters parameters, randEngine& eng)
  {
    Parameters par2 = parameters;
    up_margin = (100.0 - parameters.up_margin) / 100.0;
    low_margin = parameters.low_margin / 100.0;
    par2.surr_name = "sStudentTProcessNIG";
    par2.noise = 1e-3;
    mRobustModel.reset(PosteriorModel::create(dim,par2,eng));
  }

  RobustFiltering::~RobustFiltering(){}
  
  const Dataset* RobustFiltering::filterPoints()
  {
    mFilteredData.reset(new Dataset());
    vecOfvec XX = mRobustModel->getData()->mX;
    vectord YY = mRobustModel->getData()->mY;

    size_t n_points = mRobustModel->getData()->getNSamples();

    mRobustModel->updateHyperParameters();
    mRobustModel->fitSurrogateModel();

    for(size_t i = 0; i < n_points; ++i)
      {
        ProbabilityDistribution* pd = mRobustModel->getPrediction(XX[i]);
        double f_up = pd->quantile(up_margin);
        double f_low = pd->quantile(low_margin);
        if ((YY[i] < f_up) && (YY[i] > f_low))
          {
            mFilteredData->mX.push_back(XX[i]);
            utils::append(mFilteredData->mY, YY[i]);
            FILE_LOG(logINFO) << "Keeped value:" << YY[i] << " with thresholds:" << f_low << ", " << f_up << " and mean:" << pd->getMean();
          }
        else
          {
            FILE_LOG(logINFO) << "REMOVED value:" << YY[i] << " with thresholds:" << f_low << ", " << f_up;
          }
      }
    if (mFilteredData->getNSamples() <= n_points * 0.5)
      {
        mFilteredData->mX = XX;
        mFilteredData->mY = YY;
        FILE_LOG(logINFO) << "TOO MANY POINTS REMOVED.";
      }
    
    return mFilteredData.get();
  }
  
}
