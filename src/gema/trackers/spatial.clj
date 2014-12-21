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

(ns gema.trackers.spatial
  (:require
     [gema.utils :as utils] 
     [clojure.contrib.generic.math-functions :as gmath]))

;:::::::::::::::::::::::::
; Datatype definitions   :
;:::::::::::::::::::::::::

(defrecord tracker-value [^double path-len 
                          ^doubles coord-len 
                          ^doubles coord-sum])

;:::::::::::::::::::::::::
; Function definitions   :
;:::::::::::::::::::::::::

(defn check-input
  "Returns true if pars are compatible with the tracker, false otherwise"
  [simul-info model-info pars]
  true)

(defn compute-pars
  "Returns map with tracker info from input parameters"
  [simul-info model-info pars]
  {:save-spec [:path-len :float64 
               :coord-len [:float64] 
               :coord-sum [:float64]]})

(defn compute-init-value
  "Returns tracker initial value, using model-info and init position"
  [simul-info model-info tracker-info init-point]
  (let [dim (count (:pos init-point))]
    (->tracker-value 0.0 (double-array dim 0.0) (double-array dim 0.0))))

(defn compute-value
  "Returns tracker value for a given step using model-info, previous particle 
  state and next pos"
  [simul-info model-info tracker-info point prev-point 
   ^tracker-value prev-tracker-value]
  (let
    [pos (:pos point)
     prev-pos (:pos prev-point)
     disp (utils/diff-vec-j pos prev-pos)
     prev-path-len (:path-len prev-tracker-value)
     prev-coord-len (:coord-len prev-tracker-value)
     prev-coord-sum (:coord-sum prev-tracker-value)
     delta-path-len (gmath/sqrt (utils/dot-prod-j disp disp))]
    (->tracker-value 
      (+ prev-path-len delta-path-len) 
      (utils/sum-vec-j prev-coord-len (utils/abs-vec-j disp)) 
      (utils/sum-vec-j prev-coord-sum disp))))

;::::::::::::::::
; End of file   :
;::::::::::::::::
