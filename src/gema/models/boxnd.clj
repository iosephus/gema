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

(ns gema.models.boxnd
  (:require [gema.utils :as utils])
  (:import [cern.jet.random.engine MersenneTwister]
           [cern.jet.random Uniform Normal]))

;:::::::::::::::::::::::::
; Function definitions   :
;:::::::::::::::::::::::::
(defn check-model-pars
  "Returns true if pars are compatible with the model, false otherwise"
  [simul-info pars]
  (let [contains-dim? (contains? pars :dim)
        contains-size? (contains? pars :size)]
  (and 
    (if contains-dim? 
      (pos? (pars :dim)) 
      true)
     (if contains-size?
      (pos? (count (pars :size)))
      true)
    (if (and contains-dim? contains-size?) 
      (= (count (pars :size)) (pars :dim))
      true))))

(defn build-model-info
  "Returns map with model info from input parameters"
  [simul-info pars] 
  (let 
    [contains-dim? (contains? pars :dim)
     contains-size? (contains? pars :size)
     dim (if contains-dim? 
           (int (pars :dim)) 
           (if contains-size? (count (pars :size)) 2)) 
     size (if contains-size?
          (double-array (pars :size)) 
          (double-array dim 1.0))] 
    {:dim dim :size size :min (double-array (mapv #(* -0.5 %) size)) 
     :max (double-array (mapv #(* 0.5 %) size))}))

(defn create-pdf
  "Returns map containing the Probability Density Functions for the model"
  [simul-info model-info seed]
  (let [x-min (double -0.5)
        x-max (double 0.5)
        step-len (simul-info :step-len)
        mt19937 (MersenneTwister. ^int seed)
        pdf-uniform (Uniform. x-min 
                              x-max 
                              ^cern.jet.random.engine.RandomEngine mt19937)
        pdf-normal (Normal. 0.0 step-len mt19937)] 
    {:gen mt19937 :init-pdf pdf-uniform :step-pdf pdf-normal}))

(defn get-displacement
  "Returns a vector with a candidate displacement"
  [simul-info model-info model-pdfs from-point]
  (utils/get-rand-n (model-pdfs :step-pdf) (model-info :dim)))

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
  (let [check-min (map > new-pos (model-info :min))
        check-max (map < new-pos (model-info :max))
        check-both (map (fn [x y] (and x y)) check-min check-max)]
    (if (every? true? check-both)
      (:state from-point)
      nil)))

(defn get-initpos
  "Returns a vector with a candidate for initial position"
  [simul-info model-info model-pdfs]
  (utils/mult-vec-j (model-info :size) 
                    (utils/get-rand-n (model-pdfs :init-pdf) 
                                      (model-info :dim))))

(defn check-initpos
  "Returns the initial state the particle would have if it gets assigned the 
  initial position specified by the argument `pos`. The return value should 
  fulfill the same rules as in the `check-displacement` function."
  [simul-info model-info ^doubles pos]
  0)

;::::::::::::::::
; End of file   :
;::::::::::::::::
