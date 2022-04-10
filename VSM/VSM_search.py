#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
VSM - Version 1.0

   Volcanic and Seismic source Modelling

   Author:  Elisa Trasatti, elisa.trasatti@ingv.it
   Istituto Nazionale di Geofisica e Vulcanologia - Rome (Italy)

   Last update:  March 2021

   License:  E. Trasatti, covered by GNU-GPL License

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

      * Redistributions of source code must retain the above copyright 
        notice, this list of conditions and the following disclaimer.
      * Redistributions in binary form must reproduce the above copyright 
 	    notice, this list of conditions and the following disclaimer in 
 	    the documentation and/or other materials provided with the distribution.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPY RIGHT OWNER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
   
*******************************************************************************

   Neighborhood Algorithm direct-search optimization

   Modified from PIPy and GitHub 
   https://github.com/keithfma/neighborhood
   
   MIT License

   Copyright (c) 2018 Keith Ma

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:
    
   The above copyright notice and this permission notice shall be included in all
   copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.

*******************************************************************************

   BIBLIOGRAPHY

   Sambridge, M. (1999). Geophysical inversion with a Neighbourhood Algorithm -
       Searching a parameter space. Geophys. J. Int. 138, 479-494.
"""

from copy import copy, deepcopy
from random import uniform
import numpy as np
import pandas as pd


class Searcher():
    
    def __init__(self, objective, limits, num_samp, num_resamp, names=[], maximize=False, verbose=True):
        """
        Neighborhood algorithm direct-search optimization
        
        Arguments:
            objective: callable accepting a 1D numpy array of parameters, and
                returning a scalar misfit value
            limits: list of tuples defining range for each objective
                function parameter, as [(min_val, max_val), ...]
            num_samp: int, number of random samples taken at each iteration.
            num_resamp: int, number of best Voronoi polygons sampled at
                each iteration.
            names: list of strings, names for objective function parameters,
                used for plotting, and totally optional
            maximize: boolean, set True to maximize the objective function,
                or false to minimize it.
            verbose: set True to print verbose progress messages
        
        """
        
        # store and validate input args
        self._objective = objective
        self._limits = deepcopy(limits)
        self._num_samp = num_samp
        self._num_resamp = num_resamp
        self._names = copy(names)
        #self._names = ['x{}'.format(ii) for ii in range(len(limits))]
        self._maximize = maximize
        self._verbose = verbose
        self._validate_args()

        # init constants and working vars
        self._num_dim = len(limits)
        self._param_min = np.array([x[0] for x in limits])
        self._param_max = np.array([x[1] for x in limits])
        self._param_rng = np.array([x[1]-x[0] for x in limits])
        self._sample = []
        self._queue = [] 
        self._iter = 0
        
    @property
    def sample(self):
        return deepcopy(self._sample)

    @property
    def sample_dataframe(self):
        samps = self.sample
        for samp in samps:
            for name, val in zip(self._names, samp['param']):
                samp[name] = val
            del samp['param']
        return pd.DataFrame(samps)

    def update(self, num_iter=10):
        """
        Execute search algorithm for specified number of iterations
        
        Arguments:
            num_iter: int, number of iterations to run
        """
        
        for ii in range(num_iter):
            
            # generate new sample (populates queue)
            if not self._sample:
                self._random_sample()
            else:
                self._neighborhood_sample()
                        
            # execute forward model for all samples in queue
            while self._queue:
                param = self._queue.pop()
                result = self._objective(param)
                self._sample.append({
                    'param': param,
                    'func': result,
                    'iter': self._iter
                    })
             
            # prepare for next iteration
            self._sample.sort(key=lambda x: x['func'], reverse=self._maximize)
            self._iter += 1
            if self._verbose:
                print(self)

    def _random_sample(self):
        """Generate uniform random sample for initial iteration"""
        for ii in range(self._num_samp):
            pt  = np.random.rand(self._num_dim)*self._param_rng + self._param_min
            self._queue.append(pt)

    def _neighborhood_sample(self):
        """Generate random samples in best Voronoi polygons"""
        
        vv = np.array([x['param'] for x in self._sample])
        vv = (vv - self._param_min)/self._param_rng # normalize
        
        for ii in range(self._num_samp):
            
            # get starting point and all other points as arrays
            kk = ii % self._num_resamp  # index of start point            
            vk = vv[kk,:]
            vj = np.delete(vv, kk, 0)
            xx = vk.copy()
            
            # get initial distance to ith-axis (where i == 0)
            d2ki = 0.0
            d2ji = np.sum(np.square(vj[:,1:] - xx[1:]), axis=1)
            
            # random step within voronoi polygon in each dimension
            for ii in range(self._num_dim):
                
                # find limits of voronoi polygon
                xji = 0.5*(vk[ii] + vj[:,ii] + (d2ki - d2ji)/(vk[ii] - vj[:,ii]))
                try:
                    low = max(0.0, np.max(xji[xji <= xx[ii]]))
                except ValueError: # no points <= current point
                    low = 0.0
                try:
                    high = min(1.0, np.min(xji[xji >= xx[ii]]))
                except ValueError: # no points >= current point
                    high = 1.0

                # random move within voronoi polygon
                xx[ii] = uniform(low, high)
                
                # update distance to next axis
                if ii < (self._num_dim - 1):
                    d2ki += (np.square(vk[ii  ] - xx[ii  ]) - 
                             np.square(vk[ii+1] - xx[ii+1]))
                    d2ji += (np.square(vj[:,ii  ] - xx[ii  ]) - 
                             np.square(vj[:,ii+1] - xx[ii+1]))
                    
            # update queue
            xx = xx*self._param_rng + self._param_min # un-normalize
            self._queue.append(xx)
    
    def _validate_args(self):
        """Check argument types, throw informative exceptions"""
        # # objective
        if not callable(self._objective):
            raise TypeError('"objective" must be a callable')
        # # limits 
        if not isinstance(self._limits, list):
            raise TypeError('"limits" must be a list')
        for lim in self._limits:
            if not isinstance(lim, tuple) or len(lim) != 2:
                raise TypeError('"limits" elements must be length-2 tuples')
            if lim[1] < lim[0]:
                raise ValueError('"limits" elements must be increasing')
        # # num_samp
        if int(self._num_samp) != self._num_samp:
            raise TypeError('"num_samp" must be an integer')
        if self._num_samp < 1:
            raise ValueError('"num_samp" must be positive')
        # # num_resamp
        if int(self._num_resamp) != self._num_resamp: 
            raise TypeError('"num_resamp" must be an integer')
        if self._num_resamp < 1:
            raise ValueError('"num_resamp" must be positive')    
        if self._num_resamp > self._num_samp:
            raise ValueError('"num_resamp must be <= "num_samp"')
        # # maximize
        if not isinstance(self._maximize, bool):
            raise TypeError('maximize must be boolean: True or False')
        # # verbose
        if not isinstance(self._verbose, bool):
            raise TypeError('verbose must be boolean: True or False')

    def __repr__(self):
        try:
            out = '{}(iteration={}, samples={}, best={:.6e})'.format(
                self.__class__.__name__,
                self._iter,
                len(self._sample),
                self._sample[0]['func'])
        except IndexError:
            out = '{}(iteration=0, samples=0, best=None)'.format(self.__class__.__name__)
        return out
