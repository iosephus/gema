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


(ns gema.core-test
  (:require [clojure.test :refer :all]
            [gema.core :refer :all]
            [gema.utils :as utils]
            [gema.io :as io]
            [clojure.string :as string]
            [clojure.contrib.generic.math-functions :as gmath]
            )
  (:import [cern.jet.random.engine MersenneTwister]
           [cern.jet.random Uniform Normal]))

;:::::::::::::::::::::::::
; Function definitions   :
;:::::::::::::::::::::::::

(def null-array (double-array 3 0.0))
(def point-class (class (->point 0 0 null-array)))
(def mt19937 (MersenneTwister. (int 27527)))
(def pdf-uniform (Uniform. (double -0.5) 
                           (double 0.5) 
                           ^cern.jet.random.engine.RandomEngine mt19937))
(def pdf-normal (Normal. 0.0 1.0e-4 mt19937))

(defn particle-step-data?
  "Return true if `data` looks like a step particle data"
  [data]
  (instance? point-class data))

(defn particle-data?
  "Return true if `data` looks like data for a single particle"
  [data]
  (and
    (vector? data)
    (every? particle-step-data? data)))

(defn particle-dataset?
  "Return true if `data` looks like a particle dataset"
  [data]
  (and
    (vector? data) 
    (every? particle-data? data)))

(defn particle-has-steps?
  "Return true if the step values stored in `data` for one particle are 
  `steps`"
  [data steps]
  (= (map :step data) steps)
  )

(defn dataset-has-steps?
  "Return true if every particle in given particle dataset `data` contains 
  given `steps`"
  [data steps]
  (every? (fn [x] particle-has-steps? x steps) data))

(defn make-trajectory
  "Returns a trajectory"
  [simul-info model-info model-pdfs step-get-fn step-condition-fn 
   ^gema.core.point init-point n]
  (loop [counter n from-point init-point acc (transient [init-point])]
    (if (== counter 0)
      (persistent! acc)
      (let [new-point (get-new-point simul-info 
                                     model-info 
                                     model-pdfs 
                                     step-get-fn 
                                     step-condition-fn 
                                     from-point)]
        (recur (dec counter) new-point (conj! acc new-point))))))


(defn get-dummy-trajectory
  "Return a trajectory of size `size` with all positions set to 
  `[0.0 0.0 0.0]` and all states set to zero"
  [size]
  (let [state 0
        step-get-fn (constantly (double-array 3 0.0))
        step-condition-fn (constantly state)
        init-point (->point 0 state (double-array 3 0.0))]
    (make-trajectory {} {} {} step-get-fn step-condition-fn init-point size)))

(defn get-dummy-additive-trajectory
  "Return a trajectory with size `size` from initial position `[0.0 0.0 0.0]`, 
  constant displacements `disp` and `state` set to zero."
  [size, disp]
  (let [state 0
        step-get-fn (constantly disp)
        step-condition-fn (constantly state)
        init-point (->point 0 state (double-array 3 0.0))]
    (make-trajectory {} {} {} step-get-fn step-condition-fn init-point size)))

(defn trajectory-data?
  "Return true if argument `data` looks like data for a trajectory"
  [data] 
  (let [data-steps (map :step data)
        good-pattern (range (inc (apply max data-steps)))]
    (and 
      (every? particle-step-data? data)
      (= data-steps good-pattern))))

(defn rand-step-set
  "Return a set of randoms. step-max is the maximum allowed step value and 
  add-zero/add-max force the presence of 0/step-max in the set when set to true"
  [step-max add-zero add-max]
  (let [size (inc (rand-int step-max))
        base-set (take size (repeatedly (partial rand-int step-max)))
        base-set-1 (if add-zero (conj base-set 0) base-set)
        base-set-2 (if add-max (conj base-set-1 step-max) base-set-1)]
    (sort (distinct base-set-2)))
)

(defn gen-step-sets
  "Return a collection of step sets, some of them fixed and some random. 
  n-rand is the number of random sets and step-max is the maximum step value 
  for all sets."
  [n-rand step-max]
  (let [fixed-sets [[0] [step-max] [0 step-max]]
        rand-sets (take n-rand (repeatedly (partial rand-step-set 
                                                    step-max 
                                                    (> (rand) 0.5) 
                                                    (> (rand) 0.5))))]
    (apply conj fixed-sets rand-sets)))

(defn write-and-back
  "Return a boolean value indicating if particle data written and read-back"
  [run-results]
  (let [filename (string/join "" 
                              ["test-" 
                               (str (. System currentTimeMillis)) 
                               (str (rand-int 100)) ".bin"])
        pos-data (map :pos (flatten (run-results :particle-data)))
        _ (io/write-particle-data 
                                   (run-results :simul-info) 
                                   (run-results :model-info) 
                                   (run-results :particle-data) 
                                   filename)
        dim ((run-results :model-info) :dim)
        read-pos-data (map :pos (io/read-positions dim filename))
        ]
    (every? identity (map (fn [x y] (. java.util.Arrays equals 
                                       (:pos x) 
                                       (:pos y))) 
                          pos-data read-pos-data))))


;::::::::::::::::::::::
; Tests definitions   :
;::::::::::::::::::::::

(deftest ^:utils utils-test
  (testing "Utils"
    (testing "data format"
      (is (= Double 
             (class (utils/get-rand pdf-uniform))) 
          "Data format get-rand")
      (is (= (class null-array) 
             (class (utils/get-rand-1 pdf-uniform))) 
          "Data format get-rand-1")
      (is (= (class null-array) 
             (class (utils/get-rand-2 pdf-uniform))) 
          "Data format get-rand-2")
      (is (= (class null-array) 
             (class (utils/get-rand-3 pdf-uniform))) 
          "Data format get-rand-3")
      (is (= (class null-array) 
             (class (utils/get-rand-n pdf-uniform 10))) 
          "Data format get-rand-n"))

    (testing "array functions"
      (is (= [0.0 0.0 0.0] (vec (utils/abs-vec-j 
                                  (double-array [0.0 0.0 0.0])))) 
          "abs-vec-j with null")
      (is (= [1.0 1.0 1.0] (vec (utils/abs-vec-j 
                                  (double-array [1.0 1.0 1.0])))) 
          "abs-vec-j with unity")
      (is (= [1.0 1.0 1.0] (vec (utils/abs-vec-j 
                                  (double-array [-1.0 -1.0 -1.0])))) 
          "abs-vec-j with minus unity")
      (is (let [v (utils/get-rand-3 pdf-uniform)] 
            (= (vec (utils/abs-vec-j v)) (mapv gmath/abs v))) 
          "abs-vec-j with random")

      (is (let [v (utils/get-rand-3 pdf-uniform)] 
            (= (vec (utils/scale-vec-j v 1.0)) (vec v))) 
          "scale-vec-j with random v and one")
      (is (let [v (utils/get-rand-3 pdf-uniform)] 
            (= (vec (utils/scale-vec-j v 0.0)) [0.0 0.0 0.0])) 
          "scale-vec-j with random v and zero")
      (is (let [v (utils/get-rand-3 pdf-uniform) 
                s (utils/get-rand pdf-uniform)] 
            (= (vec (utils/scale-vec-j v s)) (mapv #(* s %) v))) 
          "scale-vec-j with random v and s")

      (is (let [v1 (double-array 3 0.0)
                v2 (utils/get-rand-3 pdf-uniform)
                ] 
            (= (utils/dot-prod-j v1 v2) 0.0)) 
          "dot-prod-j with random v1 and null")      

      (is (let [v1 (double-array 3 1.0)
                v2 (double-array 3 1.0)
                ] 
            (= (utils/dot-prod-j v1 v2) 3.0)) "dot-prod-j with unity and unity")      

      (is (let [v1 (utils/get-rand-3 pdf-uniform) 
                x1 (aget ^doubles v1 0)
                y1 (aget ^doubles v1 1)
                z1 (aget ^doubles v1 2)
                v2 (double-array 3 1.0)
                ] 
            (= (utils/dot-prod-j v1 v2) (+ x1 y1 z1))) 
          "dot-prod-j with random v1 and unity")      

      (is (let [v1 (utils/get-rand-3 pdf-uniform) 
                v2 (double-array 3 0.0)
                ] 
            (= (vec (utils/sum-vec-j v1 v2)) (vec v1))) 
          "sum-vec-j with random v1 and null")      

      (is (let [v1 (utils/get-rand-n pdf-uniform 10) 
                v2 (utils/get-rand-n pdf-uniform 10)
                ] 
            (= (vec (utils/sum-vec-j v1 v2)) (mapv + v1 v2))) 
          "sum-vec-j with random v1 and v2 length 10")      

      (is (let [v1 (utils/get-rand-3 pdf-uniform) 
                v2 (double-array 3 0.0)
                ] 
            (= (vec (utils/diff-vec-j v1 v2)) (vec v1))) 
          "diff-vec-j with random v1 and null")      

      (is (let [v1 (utils/get-rand-n pdf-uniform 10) 
                v2 (utils/get-rand-n pdf-uniform 10)
                ] 
            (= (vec (utils/diff-vec-j v1 v2)) (mapv - v1 v2))) 
          "diff-vec-j with random v1 and v2 length 10")      

      (is (let [v1 (utils/get-rand-3 pdf-uniform) 
                v2 (double-array 3 0.0)
                ] 
            (= (vec (utils/mult-vec-j v1 v2)) [0.0 0.0 0.0])) 
          "mult-vec-j with random v1 and null")      

      (is (let [v1 (utils/get-rand-3 pdf-uniform) 
                v2 (double-array 3 1.0)
                ] 
            (= (vec (utils/mult-vec-j v1 v2)) (vec v1))) 
          "mult-vec-j with random v1 and unity")      

      (is (let [v1 (utils/get-rand-n pdf-uniform 10) 
                v2 (utils/get-rand-n pdf-uniform 10)
                ] 
            (= (vec (utils/mult-vec-j v1 v2)) (mapv * v1 v2))) 
          "mult-vec-j with random v1 and v2 length 10")      

      (is (let [v1 (utils/get-rand-3 pdf-uniform) 
                v2 (utils/get-rand-3 pdf-uniform)
                x1 (aget ^doubles v1 0)
                y1 (aget ^doubles v1 1)
                z1 (aget ^doubles v1 2)
                x2 (aget ^doubles v2 0)
                y2 (aget ^doubles v2 1)
                z2 (aget ^doubles v2 2)
                ] 
            (= (utils/dot-prod-j v1 v2) (+ (* x1 x2) (* y1 y2) (* z1 z2)))) 
          "dot-prod-j with random v1 and v2")      

      (is (let [v1 (utils/get-rand-n pdf-uniform 10) 
                v2 (utils/get-rand-n pdf-uniform 10)
                ] 
            (= (utils/dot-prod-j v1 v2) (reduce + (map * v1 v2)))) 
          "dot-prod-j with random v1 and v2 length 10")      
      )
  )
)

(deftest ^:engine engine-test
  (testing "Engine"
    (testing "trajectory generation"
      (is (trajectory-data? (get-dummy-trajectory 100)) "Data format")
      (is (every? (fn [x] (= (vec x) [0.0 0.0 0.0])) 
                  (map :pos (get-dummy-trajectory 100))) 
          "Zero displacement trajectory")
      (is (= [1000.0 500.0 100.0] 
             (vec (last (map :pos 
                             (get-dummy-additive-trajectory 
                               100 (double-array [10.0 5.0 1.0])))))) 
          "Additive displacement trajectory"))
    (testing "step extraction"
      (is (let [step-sets (gen-step-sets 5 1000)
                state 0
                init-point (->point 0 state (double-array 3 0.0))
                step-get-fn (constantly (double-array 3 0.0))
                step-condition-fn (constantly state)
                ] 
            (every? identity 
                    (map (fn [x] 
                           (particle-has-steps? 
                             (walk-one-particle {} 
                                                {} 
                                                {} 
                                                step-get-fn 
                                                step-condition-fn 
                                                init-point 
                                                x
                                                [] 
                                                [] 
                                                [] 
                                                []) 
                             x)) step-sets))) 
          "At the single particle level")
      )
  )
)

(deftest ^:execution execution-test
  (let [base-opts {:verbosity (int 0) :particles 100 :steps [0 1000]  
                   :workers 6 :model "free3d" :trackers [] 
                   :random-seeds [159792212 1468112666 -329264389 -1455025928 
                                  -1932585782 190291931 1735603013 60503471 
                                  -1921729617 -2066404548 2097571644 
                                  -949912970 873618374 -1672046354 -1974365447 
                                  2144244798]}
        model-names (utils/find-ns-tails-by-prefix "gema.models.") 
        ]
    (testing "All models and tracker combinations"
      (are [model trackers] 
           (particle-dataset? ((run (assoc base-opts 
                                           :model model 
                                           :trackers trackers) 
                                    "gema.models." "gema.trackers." 
                                    "gema.reducers.") 
                               :particle-data))
           "free3d" []
           "free3d" ["spatial"]
           "free3d" ["diffmr1"]
           "free3d" ["spatial" "diffmr1"]
           "segment" []
           "segment" ["spatial"]
           "segment" ["diffmr1"]
           "segment" ["spatial" "diffmr1"]
           "boxnd" []
           "boxnd" ["spatial"]
           "boxnd" ["diffmr1"]
           "boxnd" ["spatial" "diffmr1"]
           "circle" []
           "circle" ["spatial"]
           "circle" ["diffmr1"]
           "circle" ["spatial" "diffmr1"]
           "sphere" []
           "sphere" ["spatial"]
           "sphere" ["diffmr1"]
           "sphere" ["spatial" "diffmr1"]
           "cylinder" []
           "cylinder" ["spatial"]
           "cylinder" ["diffmr1"]
           "cylinder" ["spatial" "diffmr1"]
           "acinarduct" []
           "acinarduct" ["spatial"]
           "acinarduct" ["diffmr1"]
           "acinarduct" ["spatial" "diffmr1"]
           
           )
      )))

(deftest ^:storage storage-test
  (let [base-opts {:verbosity (int 0) :particles 100 :steps [0 1000]  
                   :workers 6 :model "free3d" :trackers [] 
                   :random-seeds [159792212 1468112666 -329264389 -1455025928 
                                  -1932585782 190291931 1735603013 60503471 
                                  -1921729617 -2066404548 2097571644 
                                  -949912970 873618374 -1672046354 -1974365447 
                                  2144244798]}
        run-results01 (run base-opts 
                           "gema.models." 
                           "gema.trackers." 
                           "gema.reducers.")
        run-results02 (run (assoc base-opts :steps [0]) 
                           "gema.models." 
                           "gema.trackers." 
                           "gema.reducers.")
        run-results03 (run (assoc base-opts :steps [99]) 
                           "gema.models." "gema.trackers." 
                           "gema.reducers.")
        run-results04 (run (assoc base-opts :steps [0 105] :auto-steps true) 
                           "gema.models." 
                           "gema.trackers." 
                           "gema.reducers.")
        run-results05 (run (assoc base-opts :steps [0 101 10] :auto-steps true) 
                           "gema.models." 
                           "gema.trackers." 
                           "gema.reducers.")
        run-results05 (run (assoc base-opts :steps [0 101 20] :auto-steps true) 
                           "gema.models." 
                           "gema.trackers." 
                           "gema.reducers.")
        run-results06 (run (assoc base-opts 
                                  :steps 
                                  (rand-step-set 100 
                                                 (> (rand) 0.5) 
                                                 (> (rand) 0.5))) 
                           "gema.models." 
                           "gema.trackers." 
                           "gema.reducers.")
        run-results07 (run (assoc base-opts 
                                  :steps 
                                  (rand-step-set 100 
                                                 (> (rand) 0.5) 
                                                 (> (rand) 0.5))) 
                           "gema.models." 
                           "gema.trackers." 
                           "gema.reducers.")
        run-results08 (run (assoc base-opts 
                                  :steps 
                                  (rand-step-set 100 
                                                 (> (rand) 0.5) 
                                                 (> (rand) 0.5))) 
                           "gema.models." 
                           "gema.trackers." 
                           "gema.reducers.")
        run-results09 (run (assoc base-opts 
                                  :steps 
                                  (rand-step-set 100 
                                                 (> (rand) 0.5) 
                                                 (> (rand) 0.5))) 
                           "gema.models." 
                           "gema.trackers." 
                           "gema.reducers.")
        run-results10 (run (assoc base-opts 
                                  :steps 
                                  (rand-step-set 100 
                                                 (> (rand) 0.5) 
                                                 (> (rand) 0.5))) 
                           "gema.models." 
                           "gema.trackers." 
                           "gema.reducers.")]
    (testing "Storage" 
      (is (write-and-back run-results01))
      (is (write-and-back run-results02))
      (is (write-and-back run-results03))
      (is (write-and-back run-results04))
      (is (write-and-back run-results05))
      (is (write-and-back run-results06))
      (is (write-and-back run-results07))
      (is (write-and-back run-results08))
      (is (write-and-back run-results09))
      (is (write-and-back run-results10))
      )))
