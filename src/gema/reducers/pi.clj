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

(ns gema.reducers.pi
  (:require 
     [gema.utils :as utils]
     [clojure.contrib.generic.math-functions :as gmath]))


;:::::::::::::::::::::::::
; Function definitions   :
;:::::::::::::::::::::::::

(defn check-input
  "Returns true if pars are compatible with the reducer, false otherwise"
  [simul-info model-info pars]
  true)

(defn compute-pars
  "Returns map with reducer info from input parameters"
  [simul-info model-info pars] 
  {})

(defn compute-value
  "Returns reducer value."
  [simul-info model-info reducer-info data]
  (let
    [r2 0.25
     positions (mapv :pos (filter (fn [x] (= (:step x) 0)) (flatten data)))
     part-r2 (map utils/sum-of-sq positions)
     num-total (count part-r2)
     num-in-circle (count (filter #(<= % r2) part-r2))
     circle-area (/ (float num-in-circle) (float num-total))
     pi-value (/ circle-area r2)]
  {:pi-value pi-value}))

;::::::::::::::::
; End of file   :
;::::::::::::::::
