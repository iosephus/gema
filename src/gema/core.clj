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

(ns gema.core
  (:gen-class)
  (:require [clojure.string :as string]
            [clojure.tools.cli :refer [parse-opts]]
            [clojure.contrib.generic.math-functions :as gmath]
            [gema.io :as io]
            [gema.utils :as utils]
            [gema.models.pointnd]
            [gema.models.free3d]
            [gema.models.segment]
            [gema.models.boxnd]
            [gema.models.circle]
            [gema.models.sphere]
            [gema.models.cylinder]
            [gema.models.acinarduct]
            [gema.trackers.spatial]
            [gema.trackers.diffmr1]
            [gema.reducers.pi]
            [gema.reducers.diffmr1adc]))

;:::::::::::::::::::::::::
; Datatype definitions   :
;:::::::::::::::::::::::::

(defrecord point [^int step ^int state ^doubles pos])

;:::::::::::::::::::::::::
; Macro definitions   :
;:::::::::::::::::::::::::

(defmacro get-version []
  (System/getProperty "gema.version"))

;:::::::::::::::::::::::::
; Function definitions   :
;:::::::::::::::::::::::::


(defn extract-map-values
  "Return the existing tracker values in a particle map, given the list of 
  tracker keys and the map."
  [tracker-keys p]
  (loop [tk tracker-keys
         acc []]
    (if (seq tk)
      (recur (rest tk) (conj acc (p (first tk))))
      acc)))

(defn compute-init-tracker-vals
  "Returns new tracker values"
  [simul-info model-info tracker-infos tracker-init-fns init-point]
  (loop [tinfs tracker-infos
         tfuns tracker-init-fns 
         acc []]
    (let [pars (first tinfs)] 
      (if pars 
        (recur (rest tinfs) 
               (rest tfuns) 
               (conj acc ((first tfuns) simul-info model-info pars init-point)))
        acc))))

(defn compute-tracker-vals
  "Returns new tracker values"
  [simul-info model-info tracker-infos tracker-fns point prev-point 
   prev-tracker-values]
  (loop [tinfs tracker-infos
         tfuns tracker-fns 
         tvals prev-tracker-values 
         acc []]
    (let [pars (first tinfs)] 
      (if pars
        (recur (rest tinfs) (rest tfuns) (rest tvals) 
               (conj acc ((first tfuns) simul-info 
                                        model-info 
                                        pars 
                                        point 
                                        prev-point 
                                        (first tvals))))
        acc))))

(defn get-init-point
  "Returns a valid displacement"
  [simul-info init-model-info init-model-pdfs init-get-fn init-condition-fn]
  (loop [new-pos (init-get-fn simul-info init-model-info init-model-pdfs)]
    (let [new-state (init-condition-fn simul-info init-model-info new-pos)] 
      (if (nil? new-state)
        (recur (init-get-fn simul-info init-model-info init-model-pdfs)) 
        (->point 0 new-state new-pos)))))

(defn get-new-point
  "Returns a valid displacement"
  [simul-info model-info model-pdfs step-get-fn step-condition-fn 
   ^point from-point]
  (loop [disp (step-get-fn simul-info model-info model-pdfs from-point)]
    (let [new-pos (utils/sum-vec-j (:pos from-point) disp)
          new-state (step-condition-fn simul-info 
                                       model-info 
                                       from-point 
                                       disp 
                                       new-pos)] 
      (if (nil? new-state)
        (recur (step-get-fn simul-info model-info model-pdfs from-point)) 
        (->point (inc (:step from-point)) new-state new-pos)))))

(defn walk-one-particle
  "Walk one particle from its current position to the end of the trajectory 
  and retrieve values according to step specification"
  [simul-info model-info model-pdfs step-get-fn step-condition-fn 
   ^point init-point n-steps tracker-keys tracker-infos tracker-init-fns 
   tracker-fns]
  (let [total-steps (last n-steps)
        init-tracker-values (compute-init-tracker-vals simul-info model-info 
                                                       tracker-infos 
                                                       tracker-init-fns 
                                                       init-point)
        n-steps-1 (filter (fn [x] (> x (:step init-point))) n-steps)
        store-init? (== (:step init-point) (first n-steps))
        ] 
    (loop 
      [counter (int total-steps)
       point (get-new-point simul-info 
                            model-info 
                            model-pdfs 
                            step-get-fn 
                            step-condition-fn 
                            init-point)
       prev-point init-point
       prev-tracker-values init-tracker-values
       checkpoints-todo n-steps-1
       acc (if store-init?
             [(if (seq prev-tracker-values) 
               (apply assoc 
                      init-point 
                      (interleave tracker-keys init-tracker-values)) 
                init-point)]
             [])]
      (if (pos? counter)
        (let [pos (:pos point)
              state (:state point)
              step (:step point)
              prev-pos (:pos prev-point)
              prev-state (:state prev-point)
              tracker-values (compute-tracker-vals simul-info 
                                                   model-info 
                                                   tracker-infos 
                                                   tracker-fns 
                                                   point 
                                                   prev-point 
                                                   prev-tracker-values)
              reached-checkpoint? (== step (first checkpoints-todo))
              next-checkpoints-todo (if reached-checkpoint? 
                                      (rest checkpoints-todo) 
                                      checkpoints-todo)
              next-acc (if reached-checkpoint?
                       (let [p (if (seq tracker-values) 
                                 (apply assoc 
                                        point 
                                        (interleave tracker-keys 
                                                    tracker-values)) 
                                 point)] 
                         (conj acc p)) 
                         acc)]
        (recur (dec counter) 
               (get-new-point simul-info 
                              model-info 
                              model-pdfs 
                              step-get-fn 
                              step-condition-fn 
                              point) 
               point 
               tracker-values 
               next-checkpoints-todo 
               next-acc))
      acc))))


(defn walk-particles
  "While there is work to do, create new particles, move them n-steps, store 
  them in local transient vector"
  [simul-info init-model-info model-info init-model-pdf model-pdfs todo 
   n-steps init-get-fn init-condition step-get-fn step-condition trackers-maps]
  (let [tracker-keys (keys trackers-maps)
        trackers-maps-values (extract-map-values tracker-keys trackers-maps)
        tracker-infos (if tracker-keys (map :pars trackers-maps-values) nil)
        tracker-init-fns (if tracker-keys 
                           (map :create-fn trackers-maps-values) 
                           nil)
        tracker-compute-fns (if tracker-keys 
                              (map :compute-fn trackers-maps-values) 
                              nil)
        ]
    (loop [start-work (ref false) worker-results (transient [])] 
      (dosync (when (pos? @todo) 
                (alter todo dec)
                (ref-set start-work true)))
      (if @start-work
        (let [init-point (get-init-point simul-info 
                                         init-model-info 
                                         init-model-pdf 
                                         init-get-fn 
                                         init-condition)
              result (walk-one-particle simul-info 
                                        model-info 
                                        model-pdfs 
                                        step-get-fn 
                                        step-condition 
                                        init-point 
                                        n-steps 
                                        tracker-keys 
                                        tracker-infos 
                                        tracker-init-fns 
                                        tracker-compute-fns)]
          (recur (ref false) (conj! worker-results result))) 
        (persistent! worker-results)))))

(defn build-tracker-map
  "Returns a map for a given tracker that contains the tracker pars and the 
  functions for creating and computing tracker values"
  [simul-info model-info input-pars track-ns-prefix tracker-name]
  (let
    [full-tracker-name (string/join "" [track-ns-prefix tracker-name])
     tracker-ns (symbol full-tracker-name)
     check-fn (ns-resolve tracker-ns (symbol "check-input"))
     pars-fn (ns-resolve tracker-ns (symbol "compute-pars"))
     create-fn (ns-resolve tracker-ns (symbol "compute-init-value"))
     compute-fn (ns-resolve tracker-ns (symbol "compute-value"))
     check (check-fn simul-info model-info input-pars)
     pars (pars-fn simul-info model-info input-pars)]
    {:name tracker-name
     :check check
     :pars pars 
     :create-fn create-fn
     :compute-fn compute-fn}))

(defn build-reducer-map
  "Returns a map for a given reducer that contains the reducer pars and the 
  functions for creating and computing reducer values"
  [simul-info model-info input-pars reducer-ns-prefix reducer-name]
  (let
    [full-reducer-name (string/join "" [reducer-ns-prefix reducer-name])
     reducer-ns (symbol full-reducer-name)
     check-fn (ns-resolve reducer-ns (symbol "check-input"))
     pars-fn (ns-resolve reducer-ns (symbol "compute-pars"))
     compute-fn (ns-resolve reducer-ns (symbol "compute-value"))
     check (check-fn simul-info model-info input-pars)
     pars (pars-fn simul-info model-info input-pars)]
    {:name reducer-name
     :check check
     :pars pars 
     :compute-fn compute-fn}))

(defn compute
  "Launch the calculations concurrently"
  [todo simul-info model-info init-model-info get-initpos-fn 
   init-check-initpos-fn check-initpos-fn get-displacement-fn 
   check-displacement-fn pdf-maps tracker-maps is-init-different]
  (let
     [final-check-initpos-fn (if is-init-different 
                               (fn [initpos] (and 
                                               (init-check-initpos-fn initpos) 
                                               (check-initpos-fn initpos))) 
                               init-check-initpos-fn) 
      start-time (. System currentTimeMillis)
      futures-computing (doall (map 
                                 (fn [pdf-map] 
                                   (future (walk-particles
                                             simul-info
                                             init-model-info
                                             model-info
                                             (pdf-map :init-model)
                                             (pdf-map :model)
                                             todo 
                                             (simul-info :n-steps)
                                             get-initpos-fn
                                             final-check-initpos-fn
                                             get-displacement-fn
                                             check-displacement-fn 
                                             tracker-maps))) 
                                    pdf-maps))]
    ;; Add worker results to a global vector while dereferencing the futures
    (loop [source futures-computing acc []] 
      (if (first source) 
        (recur (rest source) (into acc (deref (first source))))
        {:data acc :start-time start-time 
         :duration (- (. System currentTimeMillis) start-time)}))))

(defn run-reducers
  "Return results from reducers after running them concurrently."
  [simul-info model-info reducer-maps particle-data]
  (let
     [start-time (. System currentTimeMillis)
      futures-reducers (doall (map 
                                 (fn [rmap] 
                                   [(first rmap) 
                                    (future (((second rmap) :compute-fn) 
                                             simul-info model-info 
                                             ((second rmap) :pars) 
                                             particle-data))]) 
                                    reducer-maps))]
    ;; Add worker results to a global vector while dereferencing the futures
    (loop [source futures-reducers acc {}] 
      (if (first source) 
        (recur (rest source) (assoc acc 
                                    (ffirst source) 
                                    (deref (second (first source)))))
        {:results acc :start-time start-time 
         :duration (- (. System currentTimeMillis) start-time)}))))


(defn run
  "Perform initialization and performs most of program tasks for a given set 
  of command line options. The main functions just parses command line 
  options, gets some values from the environment and calls this function."
  [cli-opts model-ns-prefix track-ns-prefix reducer-ns-prefix]

  (let 
    [;; Get simulation parameters from arguments
     verb-level (cli-opts :verbosity)
     n-particles (cli-opts :particles) ; Number of particles
     auto-steps? (cli-opts :auto-steps)
     n-steps (if auto-steps? 
               (doall (apply range (cli-opts :steps))) 
               (sort (distinct (cli-opts :steps)))) ; Numbers of steps

     total-steps (apply max n-steps)

     _ (utils/print-info (format "Work to do: %d particles in %d total steps" 
                                 n-particles total-steps) 
                         1 
                         verb-level)

     _ (utils/print-info (format "Values will be stored for %d steps between 
                                 %d and %d" 
                                 (count n-steps) 
                                 (apply min n-steps) 
                                 (apply max n-steps)) 
                         1 
                         verb-level)

     ;; Determine number of workers
     ;; If passed use what user requests, otherwise (+ n-cpu 2)
     ;; Never let n-workers to be greater than number of particles
     n-cpu (.availableProcessors (Runtime/getRuntime))
     n-workers (min 
                 n-particles 
                 (if (contains? cli-opts :workers) 
                   (max 1 (cli-opts :workers)) 
                   (+ 2 n-cpu)))
     _ (utils/print-info (format "Will spawn %d workers (found %d CPUs)" 
                                 n-workers 
                                 n-cpu) 
                         1 
                         verb-level)

     ;; Compute simulation derived parameters Tdiff and Dfree set to one
     ; Time delta associated to each step
     step-dur (if (pos? total-steps) (/ 1.0 total-steps) 0.0)
     ; Displacement associated to each step for each direction
     step-len (gmath/sqrt (* 2.0 step-dur)) 

     ;; Buils simul info map for passing to workers, models and trackers
     simul-info {:n-particles n-particles :n-steps n-steps 
                 :auto-steps auto-steps? :step-len step-len :step-dur step-dur}

     ;; Initialize models, trackers and PRNGs
     ;;::::::::::::::::::::::::::::::::::::::
     ;;
     ;; Find functions for selected model based on namespaces
     selected-model-ns (symbol 
                         (string/join [model-ns-prefix, (cli-opts :model)]))
   
     check-model-pars-fn (ns-resolve selected-model-ns 
                                     (symbol "check-model-pars"))
     build-model-info-fn (ns-resolve selected-model-ns 
                                     (symbol "build-model-info"))
     create-pdf-fn (ns-resolve selected-model-ns (symbol "create-pdf"))
     get-initpos-fn (ns-resolve selected-model-ns (symbol "get-initpos"))
     check-initpos-fn (ns-resolve selected-model-ns (symbol "check-initpos"))
     get-displacement-fn (ns-resolve selected-model-ns 
                                     (symbol "get-displacement"))
     check-displacement-fn (ns-resolve selected-model-ns 
                                       (symbol "check-displacement"))


     ;; Find functions for init model based on namespaces
     init-model-ns (if (contains? cli-opts :init-model)
                          (symbol (string/join 
                                    [model-ns-prefix (cli-opts :init-model)]))
                          selected-model-ns)
   
      init-check-model-pars-fn (ns-resolve init-model-ns 
                                           (symbol "check-model-pars"))
      init-build-model-info-fn (ns-resolve init-model-ns 
                                           (symbol "build-model-info"))
      init-create-pdf-fn (ns-resolve init-model-ns (symbol "create-pdf"))
   
      init-get-initpos-fn (ns-resolve init-model-ns (symbol "get-initpos"))
      init-check-initpos-fn (ns-resolve init-model-ns (symbol "check-initpos"))

      ;; Create boolean variable to check whether init conditions are different,
      ;; either because init model or init pars are different.
      is-init-different (or 
                          (not= init-model-ns selected-model-ns) 
                          (not (empty? (cli-opts :init-model-pars))))

      model-info (try 
                   (build-model-info-fn simul-info (cli-opts :model-pars)) 
                   (catch Exception e (utils/print-info 
                                        (format "Error building model -> %s" 
                                                (.getMessage e)) 
                                        1 
                                        verb-level) 
                     (System/exit 1)))
  
      init-model-info (if is-init-different 
                        (try 
                          (init-build-model-info-fn simul-info 
                                                    (cli-opts :init-model-pars))
                          (catch Exception e 
                            (utils/print-info 
                              (format "Error building init model -> %s" 
                                      (.getMessage e))) 
                            (System/exit 1))) 
                        model-info)


      _ (utils/print-info (format "Model in use: %s with pars %s" 
                                  (cli-opts :model) 
                                  (utils/pretty-str model-info)) 
                          1 
                          verb-level)
      _ (utils/print-info (if (not= model-info init-model-info) 
                  (format "Init model in use: %s with pars %s" 
                          (if (cli-opts :init-model) 
                            (cli-opts :init-model) 
                            (cli-opts :model)) 
                          (utils/pretty-str init-model-info))
                  "Init model in use is same as simulation model.") 
                          1 
                          verb-level)

     ;; Get list of requested trackers from arguments
     selected-trackers (cli-opts :trackers) 
     ;; Build tracker map for each tracker after separating input pars by 
     ;; tracker prefix
     tracker-input-pars
       (zipmap 
         selected-trackers
         (map (fn [t-name] (utils/get-pars-by-prefix 
                             (cli-opts :trackers-pars) 
                             t-name)) 
              selected-trackers))
   
     tracker-maps 
       (zipmap 
         (map keyword selected-trackers) 
         (map (fn [t-name] (build-tracker-map simul-info 
                                              model-info 
                                              (tracker-input-pars t-name) 
                                              track-ns-prefix 
                                              t-name)) 
              selected-trackers))

     tracker-formatter (fn [x] (utils/pretty-str (dissoc x :save-spec))) 

     _ (utils/print-info "Trackers in use:" 1 verb-level)
     _ (doseq [t-map (vals tracker-maps)] (utils/print-info 
                                            (format "%s with pars %s" 
                                                    (t-map :name) 
                                                    (tracker-formatter 
                                                      (t-map :pars))) 
                                            1 
                                            verb-level))
     _ (when-not (seq (vals tracker-maps)) 
         (utils/print-info "None." 1 verb-level))

     ;; Get list of requested reducers from arguments
     selected-reducers (cli-opts :reducers) 
     ;; Build reducer map for each reducer after separating input pars by 
     ;; reducer prefix
     reducers-input-pars
       (zipmap 
         selected-reducers
         (map (fn [t-name] 
                (utils/get-pars-by-prefix (cli-opts :reducers-pars) t-name)) 
              selected-reducers))
   
     reducer-maps 
       (zipmap 
         (map keyword selected-reducers) 
         (map (fn [t-name] (build-reducer-map simul-info model-info 
                                              (reducers-input-pars t-name) 
                                              reducer-ns-prefix 
                                              t-name)) 
              selected-reducers))

     reducer-formatter (fn [x] (utils/pretty-str (dissoc x :save-spec))) 

     _ (utils/print-info "Reducers in use:" 1 verb-level)
     _ (doseq [t-map (vals reducer-maps)] (utils/print-info 
                                            (format "%s with pars %s" 
                                                    (t-map :name) 
                                                    (reducer-formatter 
                                                      (t-map :pars))) 
                                            1 
                                            verb-level))
     _ (when-not (seq (vals reducer-maps)) (utils/print-info "None." 
                                                             1 
                                                             verb-level))
     

     ;; RNG Initialization
     seeds-device (cli-opts :seed-device)
     _ (utils/print-info (format "Initializing %d RNGs" n-workers) 1 verb-level)
     _ (utils/print-info (format "Getting random seeds from %s..." 
                                 (if (seq (cli-opts :random-seeds)) 
                                   "command line arguments" 
                                   seeds-device)) 
                         1 
                         verb-level)

     n-seeds (if is-init-different 
               (* 2 n-workers) 
               n-workers)

     rand-seeds  (if (seq (cli-opts :random-seeds)) 
                   (partition n-workers (cli-opts :random-seeds))
                   (let [seeds (io/read-binary-int32 seeds-device n-seeds)] 
                     (partition n-workers seeds)))

     _ (when (not= seeds-device "last-seeds.bin") 
        (io/write-binary "last-seeds.bin" (flatten rand-seeds) :int32) 
        (utils/print-info (format 
                            "Random seeds written to \"%s\"" "last-seeds.bin") 
                          1 
                          verb-level))

     ;; Create random number generators and Probability Density Functions
     pdfs (doall (map (fn [seed] (create-pdf-fn simul-info model-info seed)) 
                      (first rand-seeds)))
     init-pdfs (if is-init-different
                      (doall (map (fn [seed] (init-create-pdf-fn 
                                               simul-info 
                                               init-model-info 
                                               seed)) 
                                  (second rand-seeds)))
                      pdfs)
   
     pdf-maps (doall (map (fn [init-pdf pdf] 
                            {:init-model init-pdf :model pdf}) 
                          init-pdfs pdfs))
   
     ;; Create ref for tracking progress of computation
     todo (ref n-particles)
     ;; Launch the workers and wait for them to finish
     _ (utils/print-info "Starting computation..." 1 verb-level)

     ;; Perform concurrent computation
     ;;::::::::::::::::::::::::::::::::

     computing-results (compute todo 
                                simul-info 
                                model-info 
                                init-model-info 
                                init-get-initpos-fn 
                                init-check-initpos-fn 
                                check-initpos-fn 
                                get-displacement-fn 
                                check-displacement-fn 
                                pdf-maps 
                                tracker-maps 
                                is-init-different)

     particle-data (computing-results :data)
     computing-time-seconds (/ (double (computing-results :duration)) 1000)

     _ (utils/print-info (format "Computation done! Took %f seconds" 
                                 computing-time-seconds) 
                         1 
                         verb-level)

     n-particles-final (count particle-data)

     ;; Print info for user
     _ (utils/print-info (format "Number of particles computed: %d" 
                                 n-particles-final) 
                         1 
                         verb-level)
     _ (utils/print-info "Example particle:" 2 verb-level) 
     _ (utils/print-info (utils/pretty-str (last particle-data)) 2 verb-level)

     ;; Launch the workers and wait for them to finish
     _ (utils/print-info "Running reducers..." 1 verb-level)

     reducer-results (run-reducers simul-info 
                                   model-info 
                                   reducer-maps 
                                   particle-data)
     reducer-data (reducer-results :results)
     reducing-time-seconds (/ (double (reducer-results :duration)) 1000)

     _ (utils/print-info (format "Reducing done! Took %f seconds" 
                                 reducing-time-seconds) 
                         1 
                         verb-level)
     _ (utils/print-info "Reducers results:" 1 verb-level)
     _ (utils/print-info reducer-data 1 verb-level)
     ]


    ;; Builds model info map for passing to workers and trackers
    ;(when-not (try 
    ;            (check-model-pars-fn simul-info (cli-opts :model-pars)) 
    ;            (catch Exception e (utils/print-info (format "Error while checking model pars -> %s" (.getMessage e))) (System/exit 1)))
    ;  (utils/print-info "Supplied model parameters do not pass model input check! Quiting!")
    ;  (System/exit 1))
  
    ;(when is-init-different 
    ;  (when-not (try 
    ;              (init-check-model-pars-fn simul-info (cli-opts :init-model-pars)) 
    ;              (catch Exception e (utils/print-info (format "Error while checking init model pars -> %s" (.getMessage e))) (System/exit 1)))
    ;    (utils/print-info "Supplied init model parameters do not pass model input check! Quiting!")
    ;    (System/exit 1)))


     {:simul-info simul-info 
      :model-info model-info 
      :particle-data particle-data 
      :computing-time-seconds computing-time-seconds 
      :reducer-data reducer-data 
      :reducing-time-seconds reducing-time-seconds 
      :tracker-maps tracker-maps}))

(defn save-data
  "Save positions and tracker data to binary files"
  [simul-info model-info particle-data tracker-maps prefix]
  (let [futures-possave [(future (io/write-particle-data 
                                   simul-info 
                                   model-info 
                                   particle-data 
                                   (string/join "" [prefix "positions.bin"])))]
        futures-datasave (if (seq tracker-maps) 
                           (apply conj futures-possave 
                                  (doall 
                                    (map (fn [tmap] 
                                           (let 
                                             [tname (tmap :name) 
                                              filename (string/join 
                                                         "" 
                                                         [prefix tname ".bin"])] 
                                               (future (io/write-tracker-data 
                                                         simul-info 
                                                         model-info 
                                                         particle-data 
                                                         tmap 
                                                         filename)))) 
                                           (vals tracker-maps)))) 
                           futures-possave)]
    (doseq [f futures-datasave] (deref f))))


;::::::::::::::::::::::
; Program execution   :
;::::::::::::::::::::::

(defn -main
  "This is the main function of the program. It will be executed when a jar is 
  executed or invoked with lein run."
  [& args]

  (let
    [prog-info {:name "GEMA"
                :description "Thinking about a name..."
                :version (get-version)}

     jre-info {:name (System/getProperty "java.runtime.name")
               :vendor (System/getProperty "java.vendor")
               :version (System/getProperty "java.runtime.version")}

     jvm-info {:name (System/getProperty "java.vm.name")
               :version (System/getProperty "java.vm.version")
               :info (System/getProperty "java.vm.info")}

     os-info {:name (System/getProperty "os.name")
              :version (System/getProperty "os.version")
              :arch (System/getProperty "os.arch")}

     version-str (string/join "" 
                              [(println-str (format 
                                              "This is the Clojure 
                                              implementation of %s version %s" 
                                              (prog-info :name) 
                                              (prog-info :version)))
                               (println-str "(For help use the -h or --help 
                                            command line option)")
                               (println-str (format "JRE: %s %s on %s %s (%s)" 
                                                    (jre-info :name) 
                                                    (jre-info :version) 
                                                    (os-info :name) 
                                                    (os-info :version) 
                                                    (os-info :arch)))
                               (println-str (format "JVM: %s (build %s %s)" 
                                                    (jvm-info :name) 
                                                    (jvm-info :version) 
                                                    (jvm-info :info)))])

     ;; Find models, trackers and reducers
     ;;:::::::::::::::::::::::::::::::::::
     model-ns-prefix "gema.models."
     model-names (utils/find-ns-tails-by-prefix model-ns-prefix) 

     track-ns-prefix "gema.trackers."
     track-names (utils/find-ns-tails-by-prefix track-ns-prefix) 

     reducer-ns-prefix "gema.reducers."
     reducer-names (utils/find-ns-tails-by-prefix reducer-ns-prefix) 

     module-info-str (string/join "" [(println-str (format 
                                                     "Available models: %s" 
                                                     (string/join 
                                                       ", " 
                                                       model-names))) 
                                      (println-str (format 
                                                     "Available trackers: %s" 
                                                     (string/join 
                                                       ", " 
                                                       track-names)))
                                      (println-str (format 
                                                     "Available reducers: %s" 
                                                     (string/join 
                                                       ", " 
                                                       reducer-names)))])

     default-seed-device (if (= (os-info :name) "Linux") 
                           "/dev/random" 
                           "seeds.bin")

     cli-options-spec [["-h" "--help" "Show help and version info" :flag true 
                        :default false] 
                  ["-v" 
                   "--verbosity NUMBER" 
                   "Verbosity level (0->Quiet, 1->Normal, 2->Verbose)" 
                   :default 1 
                   :parse-fn (comp int read-string)] 
                  ["-q" 
                   "--quiet" 
                   "Suppress most of output" 
                   :flag true 
                   :default false] 
                  ["-p" 
                   "--particles NUMBER" 
                   "Number of particles" 
                   :default 1000 
                   :parse-fn (comp int read-string)] 
                  ["-s" 
                   "--steps NUMBERS" 
                   "Step values (example \"0 100 1000\")" 
                   :default [0 1000] 
                   :default-desc "\"0 1000\"" 
                   :parse-fn utils/int-vec-from-str] 
                  ["-a" 
                   "--auto-steps" 
                   "Compute steps automatically using 
                   --steps as [start stop every]" 
                   :flag true 
                   :default false] 
                  ["-m" 
                   "--model NAME" 
                   "Model (See available models below)" 
                   :default "circle"]
                  ["-M" 
                   "--model-pars PARLIST" 
                   "Parameters for selected model (example \"r 1.0 l 4.5\")" 
                   :default {} 
                   :default-desc "\"\"" 
                   :parse-fn utils/map-from-str] 
                  ["-i" 
                   "--init-model NAME" 
                   "Positions initialization model 
                   (Any available model, defaults to --model)"]
                  ["-I" 
                   "--init-model-pars PARLIST" 
                   "Parameters for positions initialization model" 
                   :default {} 
                   :default-desc "\"\"" 
                   :parse-fn utils/map-from-str]
                  ["-t" 
                   "--trackers NAMELIST" 
                   "Value trackers (See available trackers below)" 
                   :default ["spatial"] 
                   :default-desc "\"spatial\"" 
                   :parse-fn utils/split-words-norep]
                  ["-T" 
                   "--trackers-pars PARLIST" 
                   "Parameters for selected trackers 
                   (example \"trackername.mypar 0.0\")" 
                   :default {} 
                   :default-desc "\"\"" 
                   :parse-fn utils/map-from-str] 
                  ["-r" 
                   "--reducers NAMES" 
                   "Data reducers (See available reducers below)" 
                   :default [] 
                   :default-desc "\"\"" 
                   :parse-fn utils/split-words-norep]
                  ["-R" 
                   "--reducers-pars PARLIST" 
                   "Parameters for selected reducers 
                   (example \"reducername.mypar 0.0\")" 
                   :default {} 
                   :default-desc "\"\"" 
                   :parse-fn utils/map-from-str] 
                  ["-w" 
                   "--workers NUMBER" 
                   "Number of concurrent workers" 
                   :parse-fn (comp int read-string)]
                  ["-S" 
                   "--random-seeds NUMBERS" 
                   "Random seeds for PRNGs" 
                   :default [] 
                   :default-desc "\"\"" 
                   :parse-fn utils/int-vec-from-str]
                  ["-d" 
                   "--seed-device DEVICE" 
                   "Binary file containing seeds or device to read them" 
                   :default default-seed-device]
                  ]

     ;; Get command line arguments
     parse-opts-results (parse-opts args cli-options-spec)

     cli-opts (:options parse-opts-results)
     banner (:summary parse-opts-results)
     usage "Usage with jar file: java -jar <jarfile> [options]\nUsage with 
           Leiningen: lein run -m gema.core [options]"

     cli-help (string/join 
                "" 
                [usage "\n\n" "Options:" "\n" banner "\n\n" module-info-str])

    ;; Print help and exit when help switch is passed
    _ (when (cli-opts :help) 
        (println version-str) 
        (println cli-help) 
        (System/exit 0))

     verb-level (cli-opts :verbosity)

    ;; Print program and Java version info
    _ (utils/print-info version-str 1 verb-level)

    run-results (run cli-opts model-ns-prefix track-ns-prefix reducer-ns-prefix)
    
    _ (utils/print-info "Saving particle data..." 1 verb-level)

    _ (save-data (run-results :simul-info) 
                 (run-results :model-info) 
                 (run-results :particle-data) 
                 (run-results :tracker-maps) "gema-")

    _ (utils/print-info "Done. Quiting!" 1 verb-level)
    _ (shutdown-agents)] 
    0)
    
  ;; End of main function
  )


;::::::::::::::::
; End of file   :
;::::::::::::::::
