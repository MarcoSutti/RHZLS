# Esempi di accuratezza

La suite comprende quattro famiglie con tre inizializzazioni prefissate
(seed 1, 2, 3), per 36 esecuzioni. I risultati correnti sono salvati in
`accuracy_examples.mat`. Dal 12 ottobre 2026 il confronto usa lo snapshot
ufficiale Manopt descritto in [MANOPT_BASELINE.md](MANOPT_BASELINE.md),
lasciando che `steepestdescent` selezioni direttamente il proprio Armijo
di default. I controlli della migrazione sono in
[VALIDAZIONE.md](VALIDAZIONE.md).
Il solver esterno e `steepestdescent` Manopt, con lo stesso problema, punto
iniziale, retrazione e `tolgradnorm=0` per i tre metodi.
Dal 9 ottobre 2026, `maxiter` vale 200 per Rayleigh positivo, 850 per
Brockett positivo, 60 per Procrustes e 250 per l'approssimazione di rango 2.

## La distinzione necessaria nel line search

La versione precedente usa `acceptance='either'`: termina se passa Wolfe
ordinaria oppure Wolfe approssimata. Vicino al minimo, il costo puo venire
arrotondato allo stesso valore a punti diversi e la diminuzione ordinaria
puo risultare soddisfatta numericamente anche per passi troppo grandi.
Il nuovo test scalare riproduce questo caso con un costo quadratico di base
`1e20`: la modalita approssimata rifiuta il passo che supera il minimo.

E stata aggiunta l'opzione esplicita `acceptance='approximate'`, che richiede
entrambi i limiti sulla derivata, oltre al limite sul costo. La modalita
mista resta il default e non viene presentata come equivalente a questa.
Hager e Zhang descrivono le modalita V1, V2 e V3 (Wolfe, Wolfe approssimate,
oppure passaggio permanente dall'una all'altra), nella sezione 3, pagina 121,
di [Algorithm 851: CG_DESCENT](https://people.clas.ufl.edu/hager/files/cg_code.pdf).
Questo progetto implementa ora la selezione esplicita tra la modalita mista
gia presente e la richiesta delle sole Wolfe approssimate; non implementa
la politica adattiva V3 o l'intero solver a gradiente coniugato CG_DESCENT.

## Problemi, senza traslazioni artificiali del costo

1. **Rayleigh positivo sulla sfera S^11**: `f(x)=x'*D*x`, con
   `D=diag(linspace(1,5,12))`, minimo 1 e minimizzatori +/-e1.
2. **Brockett positivo sulla Stiefel St(12,2)**: `trace(X'*D*X*N)`,
   con la stessa D e `N=diag([2,1])`. I primi due autovettori minimizzano
   il costo, a meno dei segni delle colonne.
3. **Procrustes ortogonale su SO(4)**: `f(Q)=0.5*norm(Q-A,'fro')^2`,
   con `A=diag(1:4)`. Il minimo e 7, raggiunto in Q=I. Si usa la fabbrica
   Stiefel quadrata, inizializzata con determinante positivo; la retrazione
   polare mantiene questa componente delle matrici ortogonali.
4. **Migliore approssimazione di rango 2**, sulla varieta di matrici 12x12
   di rango fisso: `f(X)=0.5*norm(X-A,'fro')^2`, con
   `A=diag(linspace(2,0.2,12))`. La soluzione e la SVD troncata; il minimo
   positivo e la meta della somma dei quadrati dei valori singolari scartati.

La tabella seguente riporta le norme finali delle 36 esecuzioni con
lo snapshot ufficiale Manopt. I valori sono confermati nei risultati
correnti e nel log `matlab_official_manopt_experiments.log`.

| Esempio | Gradiente Armijo, intervallo sui 3 seed | Gradiente RHZLS approx., intervallo sui 3 seed | Iterazioni RHZLS approx. |
|---|---:|---:|---|
| Rayleigh positivo | 2.981e-8 -- 1.025e-7 | 5.070e-18 -- 1.987e-15 | 200 per tutti i seed |
| Brockett positivo | 4.216e-8 -- 4.230e-7 | 2.311e-17 -- 1.233e-15 | 850 per tutti i seed |
| Procrustes ortogonale | 3.026e-8 -- 8.843e-8 | 5.636e-32 -- 5.163e-27 | 60 per tutti i seed |
| Approssimazione rango 2 | 2.441e-8 -- 3.623e-8 | 4.660e-15 -- 3.823e-13 | 250 per tutti i seed |

Tutti i metodi completano il budget previsto: la tolleranza sul gradiente
non causa arresti. Le norme piccole sono residui calcolati, non una misura
del numero di cifre accurate: nei problemi strutturati alcune coordinate
possono continuare a diminuire sotto eps. Non vengono imposti plateau
artificiali o aggiunte iterazioni alle storie.

Per una figura illustrativa si consigliano Procrustes e l'approssimazione
a rango fisso: hanno un minimo positivo naturale e una soluzione nota.

## Riproduzione

```matlab
cd('C:/path/to/RHZLS')
startup_rhzls
set(groot, 'defaultFigureVisible', 'on')
results = run_accuracy_examples;
```

Per rigenerare solo le nuove figure dai risultati salvati:

```matlab
plot_accuracy_examples
```

Per ciascuna famiglia, un solo `accuracy_*.pdf` mostra il seed 1 e le tre
metriche per Armijo e RHZLS approssimata (le curve della modalita `either`
sono escluse dai grafici):
errore relativo del costo, gradiente diviso per la norma iniziale comune,
distanza relativa dalla soluzione. Non si generano figure `_control` separate.
Gli altri due seed restano disponibili nei risultati e nel log completo.
Gli errori del costo vengono divisi per `abs(f_star)`; se il minimo e zero,
si usa `abs(f_0-f_star)`. I dati originali nei MAT restano assoluti.
I valori arrotondati a zero vengono
visualizzati a eps. Le distanze e i gradienti sono valutati indipendentemente
dalla sottrazione dei costi. Le curve terminano all'arresto effettivo.

La finestra visualizzata comprende tutte le iterazioni: 200 per la sfera,
850 per Stiefel, 60 per Procrustes e 250 per il rango fisso.

Per selezionare questa modalita in un esempio esistente:

```matlab
opts.solver.linesearch = struct('acceptance','approximate');
run_rhzls('sphere', opts);
```

La suite di regressione comprende il controllo del passo eccessivo
accettabile solo per arrotondamento del costo. L'esito corrente dei test
con Manopt ufficiale e riportato in [VALIDAZIONE.md](VALIDAZIONE.md).
