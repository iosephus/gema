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

(ns gema.models.acinarduct
  (:require [gema.utils :as utils]
            [clojure.contrib.generic.math-functions :as gmath])
  (:import [cern.jet.random.engine MersenneTwister]
           [cern.jet.random Uniform Normal]))

;:::::::::::::::::::::::::
; Function definitions   :
;:::::::::::::::::::::::::

;; Auxiliary functions

(defn- cart2cylrosq
  "Returns the cylindrical coordinates with radius squared corresponding to the 
  cartesian coordinates specified in argument"
  [^doubles pos]
  (let 
    [x (aget pos 0)
     y (aget pos 1)
     z (aget pos 2)
     ro-sq (+ (* x x) (* y y))
     theta (. Math atan2 y x)
     results (double-array 3 ro-sq)
     _ (aset results 1 theta)
     _ (aset results 2 z)]
    results))

(defn- cell-info 
  "Returns an integer unique for each cell. Zero for inner cylinder"
  [model-info ^doubles pos-cylrosq] 
  (let
    [r-in-sq (:r-in-sq model-info)
     r-out-sq (:r-out-sq model-info)
     disc-len (:disc-len model-info)
     wedge-angle (:wedge-angle model-info)
     n-wedges (:n-wedges model-info )
     z-min (:z-min model-info )
     z-max (:z-max model-info )
     aperture-start (:aperture-start model-info )
     aperture-end (:aperture-end model-info )
     ro-sq (aget pos-cylrosq 0)
     theta (aget pos-cylrosq 1)
     z (aget pos-cylrosq 2)
     wedge-units-theta (/ theta wedge-angle)
     disc-units-z (/ (- z z-min) disc-len)
     wedge (. Math floor wedge-units-theta)
     disc (. Math floor disc-units-z)
     theta-reduced (- wedge-units-theta wedge)
     z-reduced (- disc-units-z disc)
     cell-number (int (if (and (< ro-sq r-out-sq)
                               (> z z-min)
                               (< z z-max)) 
                        (if (< ro-sq r-in-sq) 
                            0 
                            (inc (+ (* disc n-wedges) wedge)))
                        -1
                        ))
     within-aperture-wedge (and (> z-reduced aperture-start)
                                (< z-reduced aperture-end)
                                (> theta-reduced aperture-start)
                                (< theta-reduced aperture-end))]
    [cell-number within-aperture-wedge]))

;; Model functions

(defn check-model-pars
  "Returns true if pars are compatible with the model, false otherwise"
  [simul-info pars]
  (and 
    (let [default-r-in 0.5
          min-r-out (if (contains? pars :r-in) (pars :r-in) default-r-in)] 
      (if (contains? pars :r-out) (> (pars :r-out) min-r-out) true))
    (let [default-r-out 1.5
          max-r-in (if (contains? pars :r-out) (pars :r-out) default-r-out)] 
      (if (contains? pars :r-in) 
        (and (> (pars :r-in) 0.0) (< (pars :r-in) max-r-in)) 
        true))
    (if (contains? pars :l) (> (pars :l) 0.0) true)
    (if (contains? pars :n-discs) (> (pars :n-discs) 0) true)
    (if (contains? pars :n-wedges) (> (pars :n-wedges) 1) true)
    (if (contains? pars :aperture) 
      (and (>= (pars :aperture) 0.0) (<= (pars :aperture) 1.0)) 
      true))
    (if (contains? pars :allowed-cells) 
      (every? (fn [x] (>= x 0)) (pars :allowed-cells)) 
      true))

(defn build-model-info
  "Returns map with model info from input parameters"
  [simul-info pars]
  (let [pi (. Math PI)
        r-in (if (contains? pars :r-in) (pars :r-in) 0.5)
        r-out (if (contains? pars :r-out) (pars :r-out) 1.5)
        delta-r (- r-out r-in)
        chord-meas-r (* 0.5 (+ r-out r-in))
        l (if (contains? pars :l) (pars :l) 10.0)
        n-discs (if (contains? pars :n-discs) 
                  (pars :n-wedges) 
                  (int (. Math round ^double (/ l delta-r))))
        n-wedges (if (contains? pars :n-wedges) 
                   (pars :n-wedges)
                   (int (. Math round (/ pi (. Math asin 
                                               (/ delta-r 
                                                  (* 2.0 chord-meas-r)))))))
        aperture (if (contains? pars :aperture) (pars :aperture) 0.5)
        r-in-sq (* r-in r-in) 
        r-out-sq (* r-out r-out) 
        z-min 0.0
        z-max l
        aperture-start (* -0.5 (dec (. Math sqrt aperture)))
        aperture-end (* 0.5 (inc (. Math sqrt aperture)))
        disc-len (/ l n-discs)
        wedge-angle (/ (* 2.0 pi) n-wedges)
        allowed-cells (range 0 (inc (* n-discs n-wedges)))]
    {:dim (int 3) 
     :r-in (double r-in) 
     :r-in-sq (double r-in-sq)
     :r-out (double r-out)
     :r-out-sq (double r-out-sq)
     :l (double l)
     :z-min (double z-min) 
     :z-max (double z-max)
     :n-discs (int (max n-discs 1))
     :n-wedges (int (max n-wedges 4))
     :disc-len (double disc-len) 
     :wedge-angle (double wedge-angle)
     :aperture (double aperture)
     :aperture-start (double aperture-start)
     :aperture-end (double aperture-end)
     ;:allowed-cells (int-array (map int allowed-cells))
     :allowed-cells (vec (map int allowed-cells))
     }))

(defn create-pdf
  "Returns map containing the Probability Density Functions for the model"
  [simul-info model-info seed]
  (let [r-out (:r-out model-info )
        z-min (:z-min model-info )
        z-max (:z-max model-info )
        step-len (simul-info :step-len)
        mt19937 (MersenneTwister. ^Integer seed)
        pdf-uniform-xy (Uniform. ^double (* -1.0 r-out) 
                                 ^double r-out 
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

;; This function has to be implemented
(defn check-displacement
  "Returns the state the particle would have if displacement specified by the 
  argument `disp` would be performed. State can be any non-nil value. It should 
  return nil if the displacement is incompatible with the model. This functions 
  acts thus as a condition for selecting valid displacements, since 
  displacements for which it returns `nil` never happen. For each particle, the 
  simulation engine keeps track of the state returned by this function and 
  passes the value returned for the displacement corresponding to the last 
  step, as argument `state`, each time this function is called."
  [simul-info model-info from-point ^doubles with-disp ^doubles new-pos]
  (let [pos-cylrosq (cart2cylrosq (:pos from-point))
        pos-cylrosq-next (cart2cylrosq new-pos)
        [cell-number within-aperture] (cell-info model-info pos-cylrosq)
        [cell-number-next within-aperture-next] (cell-info 
                                                  model-info 
                                                  pos-cylrosq-next)]
    (if (or 
          (== cell-number cell-number-next) ; Movement within same cell 
          (and ; Movement btw cell 0 (inner cyl.) and any other (outer shell) 
               (some #(== cell-number-next %) (:allowed-cells model-info)) 
               (or (zero? cell-number) (zero? cell-number-next)) 
               (and within-aperture within-aperture-next)))
      (:state from-point)
      nil)))

(defn get-initpos
  "Returns a vector with a candidate for initial position"
  [simul-info model-info model-pdfs]
  (let [init-xy (utils/get-rand-2 (model-pdfs :init-pdf-xy))
        init-z (utils/get-rand (model-pdfs :init-pdf-z))
        results (double-array 3 (aget ^doubles init-xy 0))
        _ (aset-double results 1 (aget ^doubles init-xy 1))
        _ (aset-double results 2 init-z)
        ]
    results))

(defn check-initpos
  "Returns the initial state the particle would have if it gets assigned the 
  initial position specified by the argument `pos`. The return value should 
  fulfill the same rules as in the `check-displacement` function."
  [simul-info model-info ^doubles pos]
  (let [pos-cylrosq (cart2cylrosq pos) 
        [cell-number _] (cell-info model-info pos-cylrosq)]
    (if (some #(== cell-number %) (:allowed-cells model-info))
      0
      nil)))

;::::::::::::::::
; End of file   :
;::::::::::::::::
