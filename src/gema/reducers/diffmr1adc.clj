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

(ns gema.reducers.diffmr1adc
  (:require 
     [incanter.stats :as stats]
     [clojure.contrib.generic.math-functions :as gmath]
     [gema.utils :as utils]))


;:::::::::::::::::::::::::
; Function definitions   :
;:::::::::::::::::::::::::

(defn- compute-magnitude
  "Return the NMR signal magnitude for a set of phase values.
  Multiply phase values by scale before computing magnitude."
  [phase scale]
  (let [scaled-phase (map (partial * scale) phase)
        s-real (apply + (map gmath/cos scaled-phase))
        s-imag (apply + (map gmath/sin scaled-phase))]
    (gmath/sqrt (utils/sum-of-sq [s-real s-imag]))))

(defn- get-gxg-linb
  "Return a collection of n values for the product of the gamma and the 
  gradient so that b-values are linearly distributed and max b is 
  optimal."
  [phase n]
  (let [function (partial compute-magnitude phase)
        target (* (gmath/exp -1.0) (function 0.0))
        epsilon 1.0e-3
        init-start 0.0
        init-end 1.0
        gxg-max (utils/bsearch-invert function 
                                      target 
                                      epsilon 
                                      init-start 
                                      init-end)
        gxg-max2 (* gxg-max gxg-max)
        delta-gxg2 (/ gxg-max2 (- n 1))]
    (mapv gmath/sqrt (range 0.0 (+ gxg-max2 (* 0.5 delta-gxg2)) delta-gxg2))))

(defn- compute-value-1shape
  "Return a map with the values for one shape."
  [phase nb sqrtb-coeff]
  (let [phase1 (vec phase)
        gxg-linb (get-gxg-linb phase1 nb)
        b-gxg-one (/ 1.0 (* sqrtb-coeff sqrtb-coeff))
        b-values (mapv (fn [x] (* (* x x) b-gxg-one) ) gxg-linb)
        magnitudes (mapv (partial compute-magnitude phase1) gxg-linb)
        log-m (mapv gmath/log magnitudes)
        linear-model (stats/linear-model log-m b-values)
        adc-m (* -1.0 (second (linear-model :coefs)))
        mean-phase2 (/ (utils/sum-of-sq phase) (count phase1))
        adc-true (/ mean-phase2 (* 2.0 b-gxg-one))]
    {:true-adc adc-true :magnitude-adc adc-m :b-values b-values 
     :magnitudes magnitudes}))

(defn check-input
  "Returns true if pars are compatible with the reducer, false otherwise"
  [simul-info model-info pars]
  (if (contains? pars :nb) (> (pars :nb) 1) true))

(defn compute-pars
  "Returns map with reducer info from input parameters"
  [simul-info model-info pars] 
  (let [nb  (if (contains? pars :nb) (pars :nb) 10)]
    {:nb nb}))

(defn compute-value
  "Returns reducer value"
  [simul-info model-info reducer-info data]
  (let
    [nb (reducer-info :nb)
     sqrtb-coeff {:nmr-phase-delta 1.0 
                  :nmr-phase-square 3.46410161513775458705 
                  :nmr-phase-sin 5.13019932064745638218}
     total-steps (apply max (simul-info :n-steps))
     diffmr1-values (mapv :diffmr1 (filter (fn [x] (= (:step x) total-steps)) 
                                           (flatten data)))
     shapes (keys (first diffmr1-values))
     phase-map (zipmap shapes (map (fn [sh] (mapv sh diffmr1-values)) shapes))]
    (apply conj {} (map (fn [x] (let [shape (first x) 
                                      phase (second x)] 
                                  [shape 
                                   (compute-value-1shape 
                                     phase 
                                     nb 
                                     (sqrtb-coeff shape))])) 
                        phase-map))))

;::::::::::::::::
; End of file   :
;::::::::::::::::
