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

(ns gema.trackers.diffmr1
  (:require 
     [gema.utils :as utils]
     [clojure.contrib.generic.math-functions :as gmath]))

(defrecord tracker-value [^double nmr-phase-delta 
                          ^double nmr-phase-square 
                          ^double nmr-phase-sin])
(defrecord tracker-info [^doubles grad-dir 
                         ^double twopi-over-steps 
                         ^int total-steps 
                         ^double step-dur 
                         ^doubles t-sin 
                         ^clojure.lang.PersistentVector save-spec])

;:::::::::::::::::::::::::
; Function definitions   :
;:::::::::::::::::::::::::

(defn check-input
  "Returns true if pars are compatible with the tracker, false otherwise"
  [simul-info model-info pars]
  (if (contains? pars :grad-dir) 
    (= (count (pars :grad-dir)) (model-info :dim)) 
    true))

(defn compute-pars
  "Returns map with tracker info from input parameters"
  [simul-info model-info pars] 
  (let 
    [dim (model-info :dim)
     grad-dir (if (contains? pars :grad-dir) 
                (pars :grad-dir) 
                (vec (repeat dim 1.0)))
     grad-dir-m (gmath/sqrt (utils/sum-of-sq grad-dir))
     norm-grad (if (> grad-dir-m 0) (mapv #(/ % grad-dir-m) grad-dir) grad-dir)]
    (->tracker-info (double-array norm-grad)
                    (double 
                      (/ (* 2.0 (. Math PI)) (last (simul-info :n-steps))))       
                    (last (simul-info :n-steps))
                    (simul-info :step-dur)
                     (let [total-steps (last (simul-info :n-steps)) 
                           step-values (range 0 (inc total-steps)) 
                           twopi-over-steps (/ (* 2.0 (. Math PI)) total-steps)] 
                       
                       )
                       [:nmr-phase-delta :float64 :nmr-phase-square :float64 
                        :nmr-phase-sin :float64])
    ))

(defn compute-init-value
  "Returns tracker initial value, using model-info and init position"
  [simul-info model-info tracker-info init-point]
  (->tracker-value 0.0 0.0 0.0))

(defn compute-value
  "Returns tracker value for a given step using model-info, previous particle 
  state and next pos"
  [simul-info model-info tracker-info point prev-point 
   ^tracker-value prev-tracker-value]
  (let
    [pos (:pos point)
     step (:step point)
     prev-pos (:pos prev-point)
     prev-nmr-phase-delta (:nmr-phase-delta prev-tracker-value)
     prev-nmr-phase-square (:nmr-phase-square prev-tracker-value)
     prev-nmr-phase-sin (:nmr-phase-sin prev-tracker-value)
     dt (:step-dur tracker-info)
     total-steps (:total-steps tracker-info)
     sin-grad (. Math sin (* (:twopi-over-steps tracker-info) 
                             (- (double step) 0.5)))
     ;sin-grad (aget ^doubles (:t-sin tracker-info) step)
     grad-dir (:grad-dir tracker-info)
     middle-pos (utils/scale-vec-j (utils/sum-vec-j pos prev-pos) 0.5)
     delta-phase (* (utils/dot-prod-j middle-pos grad-dir) dt)
     nmr-phase-delta (if (== step 1)
                      (+ prev-nmr-phase-delta (/ delta-phase dt))
                      (if (== step total-steps)
                        (- prev-nmr-phase-delta (/ delta-phase dt))
                        prev-nmr-phase-delta))
     nmr-phase-square (if (<= (* 2 step) total-steps)
                       (+ prev-nmr-phase-square delta-phase) 
                       (- prev-nmr-phase-square delta-phase))
     nmr-phase-sin (+ prev-nmr-phase-sin (* sin-grad delta-phase))
     ;_ (println (format "Step in diffmr1: %d of %d" step total-steps))
     ]
    (->tracker-value nmr-phase-delta nmr-phase-square nmr-phase-sin)))

;::::::::::::::::
; End of file   :
;::::::::::::::::
