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

(ns gema.models.pointnd)

;:::::::::::::::::::::::::
; Function definitions   :
;:::::::::::::::::::::::::

(defn check-model-pars
  "Returns true if pars are compatible with the model, false otherwise"
  [simul-info pars]
  (let [contains-dim? (contains? pars :dim)
        contains-pos? (contains? pars :pos)]
  (and 
    (if contains-dim? 
      (pos? (pars :dim)) 
      true)
     (if contains-pos?
      (pos? (count (pars :pos)))
      true)
    (if (and contains-dim? contains-pos?) 
      (= (count (pars :pos)) (pars :dim))
      true))))

(defn build-model-info
  "Returns map with model info from input parameters"
  [simul-info pars] 
  (let 
    [contains-dim? (contains? pars :dim)
     contains-pos? (contains? pars :pos)
     dim (if contains-dim? 
           (int (pars :dim)) 
           (if contains-pos? (count (pars :pos)) 3)) 
     pos (if contains-pos?
          (double-array (pars :pos)) 
          (double-array dim 0.0))] 
    {:dim dim :pos pos}))

(defn create-pdf
  "Returns map containing the Probability Density Functions for the model"
  [simul-info model-info seed]
  {:pos (model-info :pos) :disp (double-array (model-info :dim) 0.0)})

(defn get-displacement
  "Returns a vector with a candidate displacement"
  [simul-info model-info model-pdfs from-point]
  (model-pdfs :disp))

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
  (:state from-point))

(defn get-initpos
  "Returns a vector with a candidate for initial position"
  [simul-info model-info model-pdfs]
  (model-pdfs :pos))

(defn check-initpos
  "Returns the initial state the particle would have if it gets assigned the 
  initial position specified by the argument `pos`. The return value should 
  fulfill the same rules as in the `check-displacement` function."
  [simul-info model-info ^doubles pos]
  0)


;::::::::::::::::
; End of file   :
;::::::::::::::::
