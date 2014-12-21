;::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
; GEMA - Extensible and modular software for Monte Carlo simulations          :
;                                                                             :
; Copyright (c) 2014 Centro Nacional de Investigaciones Cardiovasculares      :
;                       (http://www.cnic.es)                                  :
; Copyright (c) 2014 by CIBERES (http://www.ciberes.org)                      :
; All rights reserved. This program and the accompanying materials            :
; are made available under the terms of the Eclipse Public License v1.        :
; which accompanies this distribution, and is available at                    :
; http://www.eclipse.org/legal/epl-v10.html                                   :
;                                                                             :
; Contributors:                                                               :
;    Jose M. Perez-Sanchez (CIBERES) - Initial implementation                 :
;::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


;:::::::::::::::::::::::::::::::::::::::::::
; Namespace definition and module import   :
;:::::::::::::::::::::::::::::::::::::::::::

(ns gema.models.cylinder
  (:require [gema.utils :as utils]
            [clojure.contrib.generic.math-functions :as gmath])
  (:import [cern.jet.random.engine MersenneTwister]
           [cern.jet.random Uniform Normal]))

;:::::::::::::::::::::::::
; Function definitions   :
;:::::::::::::::::::::::::

(defn check-model-pars
  "Returns true if pars are compatible with the model, false otherwise"
  [simul-info pars]
  (and 
    (if (contains? pars :r) (> (pars :r) 0.0) true)
    (if (contains? pars :l) (> (pars :l) 0.0) true)))

(defn build-model-info
  "Returns map with model info from input parameters"
  [simul-info pars]
  (let [good-pars (check-model-pars simul-info pars)
        r (if (contains? pars :r) (pars :r) 1.0)
        l (if (contains? pars :l) (pars :l) 1.0)
        r2 (* r r) 
        z-min (* l -0.5)
        z-max (* l 0.5)]
    {:dim 3 :r r :r2 r2 :l l :z-min z-min :z-max z-max}))

(defn create-pdf
  "Returns map containing the Probability Density Functions for the model"
  [simul-info model-info seed]
  (let [r (model-info :r)
        z-min (model-info :z-min)
        z-max (model-info :z-max)
        step-len (simul-info :step-len)
        mt19937 (MersenneTwister. ^int seed)
        pdf-uniform-xy (Uniform. ^double (* -1.0 r) 
                                 ^double r 
                                 ^cern.jet.random.engine.RandomEngine mt19937)
        pdf-uniform-z (Uniform. ^double z-min 
                                ^double z-max 
                                ^cern.jet.random.engine.RandomEngine mt19937)
        pdf-normal (Normal. 0.0 step-len mt19937)] 
    {:gen mt19937 :init-pdf-xy pdf-uniform-xy :init-pdf-z pdf-uniform-z 
     :step-pdf pdf-normal}))

(defn get-displacement
  "Returns a vector with a candidate displacement"
  [simul-info model-info model-pdfs from-point]
  (utils/get-rand-3 (model-pdfs :step-pdf)))

(defn check-displacement
  "Returns the state the particle would have if displacement specified by the 
  argument `disp` would be performed. State can be any non-nil value. It 
  should return nil if the displacement is incompatible with the model. This 
  functions acts thus as a condition for selecting valid displacements, since 
  displacements for which it returns `nil` never happen. For each particle, 
  the simulation engine keeps track of the state returned by this function and 
  passes the value returned for the displacement corresponding to the last 
  step, as argument `state`, each time this function is called."
  [simul-info model-info from-point ^doubles with-disp ^doubles new-pos]
  (let [r2 (model-info :r2)
        z-min (model-info :z-min)
        z-max (model-info :z-max)
        x-next (aget ^doubles new-pos 0)
        y-next (aget ^doubles new-pos 1)
        z-next (aget ^doubles new-pos 2)
        r2-next (+ (* x-next x-next) (* y-next y-next))]
    (if (and (< r2-next r2) 
             (> z-next z-min) 
             (< z-next z-max))
      (:state from-point)
      nil)))

(defn get-initpos
  "Returns a vector with a candidate for initial position"
  [simul-info model-info model-pdfs]
  (let [results (double-array 3)
        init-xy (utils/get-rand-2 (model-pdfs :init-pdf-xy))
        init-z (utils/get-rand (model-pdfs :init-pdf-z))
        _ (aset-double results 0 (aget ^doubles init-xy 0))
        _ (aset-double results 1 (aget ^doubles init-xy 1))
        _ (aset-double results 2 init-z)
        ]
    results))

(defn check-initpos
  "Returns the initial state the particle would have if it gets assigned the 
  initial position specified by the argument `pos`. The return value should 
  fulfill the same rules as in the `check-displacement` function."
  [simul-info model-info ^doubles pos]
  (check-displacement simul-info 
                      model-info 
                      {:step 0 :state 0 :pos (double-array 3 0.0)} 
                      pos 
                      pos))

;::::::::::::::::
; End of file   :
;::::::::::::::::
