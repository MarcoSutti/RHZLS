# Validazione corrente

Dal 12 ottobre 2026 il pacchetto usa lo snapshot ufficiale Manopt
`6c369671d8834aa28f7e222c7d89c5f216e6953d`, datato 2026-03-03,
con `manopt_version` pari a 8.0.0. I 584 file upstream sono inclusi
senza modifiche; la provenienza e gli hash sono in `PROVENANCE.txt`
e `UPSTREAM_MANIFEST.json`. Il confronto Armijo lascia
`options.linesearch` non impostato: `steepestdescent` seleziona il proprio
default ufficiale. Si veda [MANOPT_BASELINE.md](MANOPT_BASELINE.md).

La verifica con MATLAB R2025a del 12 ottobre 2026 ha superato tutti i
**30 test**, senza fallimenti. Il log e `matlab_official_manopt_tests.log`.
I controlli includono il callback Armijo scelto dal solver, i suoi
contatori originali e il rifiuto di assegnare automaticamente ng=0 a
callback diversi da quello ufficiale.

L'esecuzione completa della migrazione a Manopt ufficiale, precedente alla
riduzione del limite Thomson, ha terminato tutte le **58 ottimizzazioni**,
salvato i tre MAT e rigenerato i nove PDF. Il log e
`matlab_official_manopt_experiments.log`, che termina con
`OFFICIAL_MANOPT_ALL_EXPERIMENTS_PASS`. I controlli confermano che tutte
le storie contengono ogni iterazione fino al maxiter previsto e che
ogni baseline usa l'Armijo ufficiale, con ng=0 nella ricerca lineare.
I totali, inclusi i controlli RHZLS `either` non disegnati, sono
nf=300326 e ng=32008. Questi contatori non includono costo e gradiente
calcolati dal solver esterno.

La verifica precedente della suite dei tre esempi candidati ha completato tutte le 18 esecuzioni
con Manopt ufficiale (`matlab_official_candidate_slopes.log`): 1000
iterazioni per Grassmann generalizzata e Thomson, 60 per Max-Cut.
Il limite di Thomson e stato successivamente ridotto a 50 per entrambi i
metodi. La nuova verifica dei 18 run, con Grassmann a 1000 e Max-Cut a 60,
e documentata in `matlab_thomson_50.log`; i MAT correnti e il PDF di Thomson
usano il nuovo limite.
Lo stesso log termina con `THOMSON_50_PASS` e verifica anche i 58 run
salvati, compresi quelli non rieseguiti: ogni storia arriva al proprio
maxiter e ogni baseline conserva il callback Armijo ufficiale.
Questa modifica non ha richiesto una nuova esecuzione delle altre suite;
il precedente run dei candidati a 100 iterazioni per Thomson resta nel
log storico `matlab_thomson_100.log`.
I controlli verificano le derivate del costo lungo le retrazioni con
differenze finite e l'equivalenza dei pairing metrici e ambientali.
Per Thomson e Max-Cut verificano anche tangenza delle velocita e norme
unitarie delle colonne. I pairing metrici riutilizzano lo stesso egrad
e la stessa velocita esatta, senza valutazioni aggiuntive del gradiente.

La suite `test_rhzls` copre ricerca lineare, derivate delle
retrazioni, dominio fixed-rank, solver e confronto con Manopt Armijo.
I tre nuovi test verificano le chiamate effettive a costo e gradiente,
l'origine gia disponibile e i tentativi di retrazione non validi.
I log `matlab_plateau_*` e le verifiche sotto documentano la dipendenza
precedente; non attestano una nuova validazione dello snapshot ufficiale.
La verifica del 2026-10-09 ha superato tutti i 29 test e completato le 58
esecuzioni con nove PDF e nove sole finestre finali nella modalita interattiva.
Il controllo dei contatori ha confermato che iterazioni, costi, norme del
gradiente, passi e distanze dalla soluzione sono identici alle 58 storie
precedenti. Il log e `matlab_evaluation_histories_validation.log`.
Successivamente sono stati aggiornati i limiti per famiglia: 200 per
Rayleigh positivo, 850 per Brockett positivo, 60 per Procrustes e Max-Cut,
250 per l'approssimazione di rango 2. Rayleigh base mantiene 250.
La nuova esecuzione completa e verificata in `matlab_plateau_experiments.log`.
La verifica dopo la pulizia delle cartelle del 2026-10-08 ha confermato
26 test superati, nessun fallimento e la rigenerazione di tutti i nove PDF
dai MAT correnti. Il log e `matlab_cleanup_validation.log`; `archive`
rimane escluso dal path MATLAB.
L'esempio di matrix completion, con e senza cost offset, e stato rimosso.
La geometria delle matrici di rango fisso resta disponibile per
l'approssimazione di rango 2 della suite di accuratezza.

Gli esempi spettrali predefiniti sono Brockett/Stiefel e Rayleigh/sfera.
La suite `run_accuracy_examples` aggiunge quattro famiglie confrontando
Armijo, RHZLS mista e RHZLS con Wolfe approssimate su tre seed prefissati.
Risultati e criteri sono in [ACCURACY_EXAMPLES.md](ACCURACY_EXAMPLES.md).

I nove PDF hanno due pannelli. Quello superiore mostra tre quantita
normalizzate: errore relativo del costo,
gradiente rispetto alla norma iniziale e distanza relativa dalla soluzione.
Quello inferiore mostra `nf(k)`, il numero di chiamate al costo durante
ogni ricerca lineare. `ng(k)` e salvato separatamente nei MAT. I contatori
escludono le valutazioni esterne del solver e le metriche diagnostiche.
Le figure mostrano tutte le iterazioni registrate, compresi i plateau.
La tolleranza sul gradiente e disabilitata (`tolgradnorm=0`) in tutti gli esperimenti.
I dati originali rimangono nei file MAT.

Tre ulteriori esempi derivati dal codice ufficiale Manopt sono descritti
in [MANOPT_CANDIDATES.md](MANOPT_CANDIDATES.md).

Per rigenerare i PDF disponibili, senza ripetere l'ottimizzazione:

```matlab
startup_rhzls
set(groot, 'defaultFigureVisible', 'off')
run('docs/plot_validation.m')
```

Per ripetere tutte le 58 esecuzioni e generare i nove PDF con un comando:

```matlab
startup_rhzls
results = run_all_experiments;
```

`run_all_experiments(false)` esporta senza finestre. Il runner ripristina
cartella corrente e visibilita delle figure; il log della verifica completa
di migrazione e `matlab_official_manopt_experiments.log`; il successivo
aggiornamento dei candidati e in `matlab_thomson_50.log`.

Per usare una copia esterna della release ufficiale Manopt 8.0:

```matlab
cd('C:/path/to/RHZLS')
clear functions
startup_rhzls('C:/path/to/Manopt_8.0')
run_all_experiments
```

I percorsi indicano le cartelle locali dei due pacchetti estratti. Lo
startup accetta anche la sottocartella interna `manopt` e mantiene la
selezione nelle chiamate successive senza argomenti della stessa sessione.
`clear functions` va eseguito prima della selezione, perche la azzera.
Il controllo facoltativo `report = test_rhzls; assert(report.failed == 0);`
puo essere eseguito dopo startup e prima del runner completo.
Il runner ricrea i problemi e i risultati, senza usare i callback dei MAT
preesistenti che contengono percorsi della macchina di origine.

La verifica separata della release numerata Manopt 8.0 del 12 ottobre 2026
ha superato tutti i 30 test e completato tutte le 58 ottimizzazioni con
i nove PDF in MATLAB R2025a. E stata eseguita in una copia di RHZLS
spostata in un'altra cartella, escludendo completamente `external` e
selezionando la release scaricata da
[Manopt_8.0.zip](https://www.manopt.org/downloads/Manopt_8.0.zip).
Il log `matlab_manopt_8_external_validation.log` termina con
`MANOPT_8_EXTERNAL_PORTABILITY_PASS`: ogni storia contiene tutte le
iterazioni fino a maxiter e ogni baseline usa il callback ufficiale
`@linesearch`, con ng=0. I totali della verifica esterna sono nf=223351
e ng=26283. I dati correnti distribuiti con lo snapshot GitHub non sono
stati sostituiti dalla verifica della release esterna.

Il precedente controllo `validate_plateau_iterations` aveva verificato tutti i 58 run:
nessuna iterazione mancante fino a maxiter. Il log e
`matlab_plateau_iterations.log`. In quella verifica i nove PDF avevano pagina 764 x 477 pt,
riempimenti trasparenti alpha=0.4 e legende interamente dentro la pagina.

La verifica attuale dei nove PDF conferma una dimensione comune di pagina
1086 x 612 pt, trasparenza alpha=0.4 e testo interamente dentro la pagina.
Il controllo visivo conferma le legende fuori dai pannelli e la loro
leggibilita. La selezione da root, percorso relativo e sottocartella
`manopt` e verificata in `matlab_official_manopt_paths.log`; i tre MAT
registrano versione e commit dello snapshot usato.
