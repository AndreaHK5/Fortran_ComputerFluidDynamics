validazione su S/L vortex

corrvel/corrvort sono le analisi di correlazione fatte con i parametri in
"inputcr"

plot 'corrvel' u 1:5 w l
plot 'corrvort' u 1:3 w l


slspot contiene l'analisi di spettro nella posizione dettata da inputcr1, che
conferma la frequenza di S/L rolling

plot 'slstudy' u 1:($6*0.0039) w l
