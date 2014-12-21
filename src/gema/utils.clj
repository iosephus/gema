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

(ns gema.utils
  (:require [clojure.string :as string]
            [clojure.pprint :as pprint]
            [clojure.contrib.generic.math-functions :as gmath]))


;:::::::::::::::::::::::::
; Function definitions   :
;:::::::::::::::::::::::::

(defn get-rand
  "Returns single random value (double) from distribution"
  [^cern.jet.random.AbstractDistribution distribution]
  (. distribution nextDouble))

(defn get-rand-1
  "Returns a double array with single random value from distribution"
  [^cern.jet.random.AbstractDistribution distribution]
  (double-array 1 (. distribution nextDouble)))

(defn get-rand-2
  "Returns a double array with 2 random numbers from distribution"
  [^cern.jet.random.AbstractDistribution distribution]
  (let [results (double-array 2 (. distribution nextDouble))
        _ (aset results 1 (. distribution nextDouble))]
    results))

(defn get-rand-3
  "Returns a double array with 3 random numbers (double) from distribution"
  [^cern.jet.random.AbstractDistribution distribution]
  (let [results (double-array 3 (. distribution nextDouble))
        _ (aset results 1 (. distribution nextDouble))
        _ (aset results 2 (. distribution nextDouble))]
    results))

(defn get-rand-n
  "Returns a double array with n random numbers from distribution"
  [^cern.jet.random.AbstractDistribution distribution n]
  (loop [idx (dec n) r (double-array n)] 
    (if (== idx -1) 
      r
      (let [_ (aset-double r idx (. distribution nextDouble))] 
        (recur (dec idx) r)))))

(defn sum-of-sq
  "Returns the sum of squares of vector components"
  [v]
  (reduce + (mapv * v v)))

(defn dot-prod 
  "Returns the dot product of two vectors"
  [v1 v2]
  (reduce + (map * v1 v2)))

(defn test-vec4int-j
  "Returns true if integer is in vector"
  [^ints v item]
  (let [target (int item)] 
    (loop [idx (int (dec (alength v)))] 
      (if (= target (aget v idx)) 
        true 
        (if (= idx 0) 
          false 
          (recur (dec idx)))))))

(defn test-vec4double-j
  "Returns true if double is in vector"
  [^doubles v ^double item]
  (loop [idx (dec (alength v))]
    (if (== item (aget v idx))
      true
      (if (== idx 0)
        false
        (recur (dec idx))))))

(defn scale-vec-j
  "Returns the product of a double array by a scalar"
  [^doubles v1 ^double scale]
  (amap v1 idx ret (* (aget v1 idx) scale)))

(defn abs-vec-j
  "Returns a double array with the abs values of input double array"
  [^doubles v1]
  (amap v1 idx ret (. Math abs (aget v1 idx))))

(defn sum-vec-j
  "Returns the sum of two arrays of doubles"
  [^doubles v1 ^doubles v2]
  (amap v1 idx ret (+ (aget v1 idx) (aget v2 idx))))

(defn mult-vec-j
  "Returns the per-element product of two arrays of doubles"
  [^doubles v1 ^doubles v2]
  (amap v1 idx ret (* (aget v1 idx) (aget v2 idx))))

(defn diff-vec-j
  "Returns the difference of two arrays of doubles"
  [^doubles v1 ^doubles v2]
  (amap v1 idx ret (- (aget v1 idx) (aget v2 idx))))

(defn dot-prod-j
  "Returns the dot product of two arrays of doubles"
  [^doubles v1 ^doubles v2]
  (areduce v1 i ret 0.0
           (+ ret (* (aget v1 i)
                     (aget v2 i)))))

(defn find-ns-by-prefix
  "Returns a list of available namespaces whose names starts with `prefix`"
  [prefix]
  (filter (fn [x] (let 
                    [found-ns (name (ns-name x)) 
                     cmplen (min (count found-ns) (count prefix))] 
                    (= (subs found-ns 0 cmplen) prefix))) 
          (all-ns)))

(defn find-ns-tails-by-prefix
  "Returns a list of available namespaces whose names starts with `prefix`.
  The prefix is removed from the names."
  [prefix]
  (let [ns-list (find-ns-by-prefix prefix)]
    (map (fn [x] (subs ((comp name ns-name) x) (count prefix))) ns-list)))

(defn split-words
  "Returns a vector of strings with the space separated words present in the 
  input string"
  [input-str]
  (vec (filter seq (string/split (string/trim input-str) #"\s+"))))

(defn split-words-norep
  "Returns a vector of strings with the space separated words present in the 
  input string, include each distinct word only once"
  [input-str]
  (vec (distinct (split-words input-str))))

(defn double-vec-from-str
  "Returns a vector with doubles for each of the words in input"
  [input-str]
  (mapv (comp double read-string) (split-words input-str)))

(defn int-vec-from-str
  "Returns a vector with integers for each of the words in input"
  [input-str]
  (mapv (comp int read-string) (split-words input-str)))

(defn num-vec-from-str
  "Returns a vector containing the evaluation of each of the words in input"
  [input-str]
  (mapv (comp num read-string) (split-words input-str)))

(defn num-vec-from-str
  "Returns a vector containing the evaluation of each of the words in input"
  [input-str]
  (mapv #(num (read-string %)) (split-words input-str)))

(defn map-from-str
  "Returns a map from a string"
  [input-str]
  (let [words-or-vecs (re-seq #"\[[0-9\.\s]*?\]|[\S]+" input-str)]
    (when (odd? (count words-or-vecs)) 
      (throw (Exception. "Number of total elements should be even")))
    (let [pairs (partition 2 words-or-vecs)
          my-keys (mapv (comp keyword first) pairs)
          _ (when (not= (count (distinct my-keys)) (count my-keys)) 
              (throw (Exception. "Parameter names should be unique")))
          my-vals (mapv (comp #(let [value (read-string %)] 
                                 (if (or (number? value) (vector? value)) 
                                   value %)) 
                              second) 
                        pairs)
          ]
      (zipmap my-keys my-vals))))

(defn get-pars-by-prefix 
  "Returns a map with the same values and keys changed from 'prefix.a' to 'a'"
  [input-map prefix]
  (apply merge {} (for 
                    [[k v] input-map 
                     :when (= prefix (first (string/split (name k) #"\." 2)))] 
                    [(keyword (second (string/split (name k) #"\." 2))) v])))

(defn- find-bsearch-interval
  "Return an interval containing an x such that (= (function x) target) is true.
  The supplied function should be a monotonous function of x."
  [function target init-start init-end]
  (let [increasing? (< (function init-start) (function init-end))
        decrease-start? (if increasing?
                          (fn [x] (> (function x) target)) 
                          (fn [x] (< (function x) target)))
        increase-end? (complement decrease-start?)]
    (loop [start init-start
           end init-end] 
      (if (and 
            (not (decrease-start? start)) 
            (not (increase-end? end)))
        [start end] 
        (let [new-start (if (decrease-start? start) 
                          (- start (- end start)) 
                          start)
              new-end (if (increase-end? end) 
                        (+ end (- end start)) 
                        end)
              ] (recur new-start new-end))))))

(defn bsearch-invert
  "Returns the argument of a one-dimensional decreasing function for the 
  supplied value"
  [function target epsilon init-start init-end]
  (let [increasing? (< (function init-start) (function init-end))
        increase-start? (if increasing?
                          (fn [x] (< (function x) target)) 
                          (fn [x] (> (function x) target)))
        decrease-end? (complement increase-start?) 
        [loop-init-start loop-init-end] (find-bsearch-interval function 
                                                               target 
                                                               0.0 
                                                               1.0)]
    (loop [start loop-init-start 
           end loop-init-end]
      (let [x (* 0.5 (+ start end))] 
        (if (< (gmath/abs (- target (function x))) epsilon) 
          x
          (recur
            (if (increase-start? x) x start)
            (if (decrease-end? x) x end)))))))

(defn print-info
  "Send info to log or standard output if option quiet is not activated"
  [message message-level verbosity-level]
  (when (>= verbosity-level message-level) (println message))
  )

(defn pretty-str
  "Returns a nice string representation of the object"
  [object]
  (with-out-str (pprint/write object))
  )

(defn format-arrays
  "Return a map with any array converted to vector in the original map"
  [map]
  )

;::::::::::::::::
; End of file   :
;::::::::::::::::
