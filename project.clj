;::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
; GEMA - Extensible and modular software for Monte Carlo simulations          :
;                                                                             :
; Copyright (c) 2014 by Centro Nacional de Investigaciones Cardiovasculares   :
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


(defproject gema "0.4.2"
  :description "Geometric Model of the Acinus (GEMA) written in Clojure"
  :url ""
  :license {:name "Closed source"}
  :dependencies [[org.clojure/clojure "1.5.1"]
                 [org.clojure/clojure-contrib "1.2.0"]
                 [org.clojure/tools.cli "0.3.1"]
                 [colt "1.2.0"]
                 [overtone/byte-spec "0.3.1"]
                 [incanter "1.5.4"]]
  :main gema.core
  :global-vars {*warn-on-reflection* true}
  :profiles {:uberjar {:aot :all}}
  :jvm-opts [ "-Xms6G" "-Xmx6G"]
  :test-selectors {:default :all
                     :dummy :dummy
                     :utils :utils
                     :engine :engine
                     :execution :execution
                     :storage :storage
                     :all (constantly true)} 
  )
