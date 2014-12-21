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

(ns gema.io
  (:require [clojure.string :as string]
            [clojure.contrib.io :as clio]
            [clojure.java.io :as jio]
            [overtone.byte-spec :as bytespec]))


;:::::::::::::::::::::::::
; Function definitions   :
;:::::::::::::::::::::::::

(defn read-bytes
  "Read n bytes from device, for example /dev/random"
  [device n]
  (with-open [stream (jio/input-stream device)] 
    (loop 
      [buffer (byte-array n)
       read-len 0] 
      (if (< read-len n) 
        (recur buffer 
               (+ read-len (.read stream buffer read-len (- n read-len)))) 
        buffer))))

(defn read-binary-int32
  "Reads n data types in binary mode from device"
  [device n]
  (let 
    [bytes-to-read (* 4 n)
     byte-data (map byte-array (partition 4 (read-bytes device bytes-to-read)))
     bspec (bytespec/make-spec "bspec" [:value :int32])]
    (map (comp :value (partial bytespec/spec-read-bytes bspec)) byte-data)))

(defn write-binary
  "Write a collection to a device in binary mode with byte conversion 
  according to data-type"
  [device coll data-type]
  (let [bspec (bytespec/make-spec "bspec" [:values [data-type]])] 
    (with-open [stream (jio/output-stream device)] 
      (.write stream ^bytes (bytespec/spec-write-bytes bspec {:values coll})))))

(defn write-particle-data
  "Write particle positions to a binary file as a stream of 64bit floats"
  [simul-info model-info particles filename]
  ;; Save basic data for particles
  (let [pos-spec (bytespec/make-spec "pos-spec"  [:pos [:float64]])] 
    (with-open [stream (jio/output-stream filename)] 
      (doseq [part (flatten particles)] 
        (.write stream ^bytes (bytespec/spec-write-bytes pos-spec part))))))

(defn read-positions
  "Return a collections of maps, each map will have a position vector under 
  the `:pos` key."
  [dim filename]
  (let [pos-spec (bytespec/make-spec "pos-spec"  [:value :float64])
        raw-data (map byte-array 
                      (partition 8 
                                 (clio/to-byte-array (clio/as-file filename))))
        f64-data (loop [input raw-data 
                        acc (transient [])] 
                   (if (first input) 
                     (recur (rest input) 
                            (conj! acc 
                                   ((bytespec/spec-read-bytes pos-spec 
                                                              (first input)) 
                                    :value)))
                     (persistent! acc)))]
    (vec (doall (for [pos (partition dim f64-data)] 
                  {:pos (double-array pos)})))))

(defn read-particle-data
  "Return a particle dataset similar to what you get from a simulation, 
  without state values."
  [simul-info model-info particles filename]
  ;; Save basic data for particles
  (loop [pos-data (read-positions (model-info :dim) filename)
         n-steps (simul-info :n-steps)
         acc []]
    (if (first pos-data)
      (let [n-steps-1 (if (first n-steps) n-steps (simul-info n-steps))] 
        (recur (rest pos-data) 
               (rest n-steps-1) 
               (conj acc (assoc (first pos-data) :step (first n-steps-1)))))
      (vec (for [p  (partition (count (simul-info n-steps)) acc)] ((vec p)))))))

(defn guess-save-spec
  ""
  [tracker-map]
  {})

(defn write-tracker-data
  "" 
  [simul-info model-info particles tracker-map filename]
  (let
    [t-name (tracker-map :name)
     save-spec (if (contains? (tracker-map :pars) :save-spec) 
                 (:save-spec (tracker-map :pars)) 
                 (guess-save-spec tracker-map)) 
     t-spec (bytespec/make-spec "t-spec" save-spec)
     t-ospec (bytespec/make-spec "t-ospec" [(keyword t-name) t-spec])]
    (with-open [stream (jio/output-stream filename)] 
      (doseq 
        [part (flatten particles)] 
        (.write stream ^bytes (bytespec/spec-write-bytes t-ospec part))))))

;::::::::::::::::
; End of file   :
;::::::::::::::::
