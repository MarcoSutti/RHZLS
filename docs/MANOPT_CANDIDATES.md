# Altri esempi dal codice di Manopt

Lo script `explore_manopt_examples.m` adatta tre problemi degli esempi
ufficiali di [Manopt](https://github.com/NicolasBoumal/manopt/tree/master/examples)
per studiare il confronto tra Armijo e RHZLS con Wolfe approssimate.
Sono esperimenti separati dalla suite già presente. Le funzioni costo hanno
minimi naturalmente diversi da zero: non viene aggiunto alcun cost offset.
La differenza di accuratezza va misurata, senza presumere che debba apparire
in ogni problema o per ogni punto iniziale.

## Protocollo riproducibile

Per ogni problema si usano i seed `1:3` con `rng(seed, 'twister')`.
Per ciascun seed, entrambi i metodi ricevono lo stesso problema, lo stesso
punto iniziale e la stessa retrazione. Entrambi chiamano `steepestdescent`
di Manopt, con `maxiter = 1000` per Grassmann generalizzata, 50 per Thomson
e 60 per Max-Cut su C5, `tolgradnorm = 0`,
`minstepsize = 0` e le stesse altre opzioni.

La sola differenza è la ricerca lineare selezionata:

- Manopt ufficiale: nessun override di `options.linesearch`, quindi
  `steepestdescent` seleziona il proprio Armijo `@linesearch` di default;
- `rhzls.manopt_linesearch`, con `acceptance = 'approximate'`: RHZLS con
  accettazione tramite Wolfe approssimate.

La provenienza del baseline Armijo è descritta in `MANOPT_BASELINE.md`.
Gli esempi ufficiali chiamano anche altri solutori: qui ne adattiamo le
funzioni costo e le geometrie, usando steepest descent per isolare il
contributo della ricerca lineare. Si eseguono 18 ottimizzazioni in totale.

Lo script verifica anche la velocità della retrazione e la derivata del
costo lungo la curva con differenze finite centrate (`t = 0.15`,
`h = 1e-6`, errore relativo richiesto inferiore a `1e-7`), e controlla
che la velocità a `t = 0` coincida con la direzione tangente.

Le quantità per i grafici sono sempre normalizzate:

\[
\frac{|f_k-f_\star|}{|f_\star|},\qquad
\frac{\|\operatorname{grad} f(x_k)\|}{\|\operatorname{grad} f(x_0)\|},
\qquad e_{\mathrm{solution}}(x_k).
\]

Le norme del gradiente e gli scarti assoluti nel costo restano disponibili
nei risultati grezzi. L'errore della soluzione è già relativo e tiene
conto delle simmetrie di ciascun problema.

## 1. Autospazio dominante su Grassmann generalizzata

Fonte: [generalized_eigenvalue_computation.m](https://github.com/NicolasBoumal/manopt/blob/master/examples/generalized_eigenvalue_computation.m).

Si minimizza

\[
f(X)=-\operatorname{tr}(X^T A X),\qquad X^T B X=I_2,
\]

su `grassmanngeneralizedfactory(12, 2, B)`. Nell'istanza controllata,
`B = diag(linspace(1,2,12))`, `L = chol(B,'lower')`,
`lambda = linspace(1,5,12)` e `A = L*diag(lambda)*L'`.
Le coordinate sbiancate sono `Q = L'*X`: il sottospazio ottimo è quello
generato dagli ultimi due vettori canonici. Pertanto

\[
f_\star=-(\lambda_{11}+\lambda_{12})=-\frac{106}{11}.
\]

Il punto iniziale deriva dal QR di una matrice gaussiana `12 x 2`.
L'errore della soluzione è

\[
e_{\mathrm{solution}}(X)=\frac{\|QQ^T-P_\star\|_F}{\sqrt{2}},
\]

dove `Pstar` è il proiettore sulle ultime due coordinate. Questo evita
di penalizzare una rotazione della base di uno stesso sottospazio.

Si usa la retrazione polare nelle coordinate sbiancate:

\[
R_X(tD)=L^{-T}\operatorname{polar}(L^T(X+tD)),\qquad
\frac{d}{dt}R_X(tD)=L^{-T}\frac{d}{dt}
\operatorname{polar}(L^T(X+tD)).
\]

La derivata polare viene calcolata da `rhzls.retractions.stiefel_curve`.
Lo script assegna esplicitamente questa stessa retrazione a entrambi
i metodi; non mescola retrazioni diverse nel confronto.

## 2. Problema di Thomson: tetraedro su un prodotto di sfere

Fonte: [thomson_problem.m](https://github.com/NicolasBoumal/manopt/blob/master/examples/thomson_problem.m).

Quattro punti su `S^2` sono rappresentati dalle colonne di
`X` su `obliquefactory(3,4)`. Si conserva il costo normalizzato ufficiale:

\[
f(X)=\frac{1}{16}\sum_{i<j}\frac{1}{1-x_i^T x_j}.
\]

Il minimo globale è un tetraedro regolare, con prodotti scalari
fuori diagonale pari a `-1/3`, e

\[
f_\star=\frac{6}{16}\frac{3}{4}=\frac{9}{32}.
\]

Infatti `norm(sum(X,2))^2 >= 0` implica che il prodotto scalare medio
delle coppie è almeno `-1/3`. La convessità di `1/(1-s)` fornisce il
limite inferiore precedente, raggiunto dal tetraedro.

L'esempio ufficiale propone inizialmente un gradiente a differenze
finite, inadatto a questo confronto di accuratezza. Qui il gradiente
ambiente è analitico:

\[
\nabla f(X)=\frac{XW}{16},\qquad
W_{ij}=(1-x_i^T x_j)^{-2}\ (i\ne j),\qquad W_{ii}=0.
\]

Il punto iniziale è casuale sulla varietà. Per l'errore della soluzione
si usa il Gram di riferimento `Gstar = (4/3)*eye(4)-ones(4)/3`:

\[
e_{\mathrm{solution}}(X)=
\frac{\|X^T X-G_\star\|_F}{\|G_\star\|_F}.
\]

Questo errore è invariante alle rotazioni e, per il tetraedro regolare,
alla rinumerazione dei vertici. I punti coincidenti danno costo infinito;
lo script gestisce esplicitamente tale caso.

## 3. Max-Cut rilassato sul ciclo C5

Fonti: [maxcut.m](https://github.com/NicolasBoumal/manopt/blob/master/examples/maxcut.m)
e [elliptope_SDP.m](https://github.com/NicolasBoumal/manopt/blob/master/examples/elliptope_SDP.m).

Si considera il Laplaciano `L` del ciclo non pesato con cinque vertici,
e cinque vettori unitari nel piano, memorizzati nelle colonne di
`X` su `obliquefactory(2,5)`. Il costo della rilassazione è

\[
f(X)=-\frac14\operatorname{tr}(XLX^T),\qquad
\nabla f(X)=-\frac12 XL.
\]

È la trasposta della rappresentazione con righe unitarie usata nel
Max-Cut ufficiale. Qui il rango resta fisso a due: non si applicano
gli aumenti di rango o il rounding del programma completo.

Il riferimento ha colonne

\[
x_{\star,j}=\begin{bmatrix}
\cos(4\pi j/5)\\ \sin(4\pi j/5)
\end{bmatrix},\qquad j=0,\ldots,4,
\]

e costo minimo della rilassazione

\[
f_\star=-\frac52(1+\cos(\pi/5))\simeq-4.522542486.
\]

Il limite spettrale segue da `lambda_max(L) = 2+2*cos(pi/5)`
e `trace(X*X') = 5`; la configurazione a pentagramma lo raggiunge.
Il punto iniziale è casuale. L'errore della soluzione è lo scarto
relativo del Gram rispetto a `Xstar'*Xstar`, così le trasformazioni
ortogonali nel piano non alterano l'indicatore.

## Retrazione e velocità per i due problemi oblique

Per ogni colonna, con `z_i = x_i+t*d_i` e `s_i = norm(z_i)`,

\[
y_i=\frac{z_i}{s_i},\qquad
\dot y_i=\frac{d_i-y_i(y_i^T d_i)}{s_i}.
\]

La funzione `oblique_curve` in `+rhzls/+retractions`, richiamata anche
dal wrapper locale in `explore_manopt_examples.m`,
valuta queste formule. La normalizzazione è indipendente per ciascuna
colonna: non si normalizza l'intera matrice con una singola norma.

## Candidati da valutare successivamente

- [positive_definite_intrinsic_mean.m](https://github.com/NicolasBoumal/manopt/blob/master/examples/positive_definite_intrinsic_mean.m):
  media intrinseca di matrici SPD. Con dati commutanti e ben condizionati
  esiste un riferimento esatto tramite la media dei logaritmi. La
  retrazione locale è `X+t*D+0.5*t^2*D*(X\D)`, con velocità
  `D+t*D*(X\D)`. Il codice ufficiale avverte che il calcolo del gradiente
  può limitare l'accuratezza con molto rumore: questa limitazione va
  distinta da quella della ricerca lineare.
- [doubly_stochastic_denoising.m](https://github.com/NicolasBoumal/manopt/blob/master/examples/doubly_stochastic_denoising.m):
  proiezione in norma di Frobenius su matrici bistocastiche positive.
  Richiederebbe la derivata dello scaling Sinkhorn. La factory Manopt
  inclusa nel progetto usa PCG con tolleranza relativa `1e-6` per il
  gradiente e la proiezione, più uno scaling iterativo e un clipping
  finale. Questi errori numerici possono prevalere su quelli della
  ricerca lineare a precisioni molto elevate. Non è incluso nel pilot.

## Esecuzione e risultati

Da MATLAB:

```matlab
cd('C:/path/to/RHZLS')
startup_rhzls;
results = explore_manopt_examples();
plot_manopt_examples(results);
```

I risultati correnti sono salvati in `manopt_candidate_results.mat`.
Dal 12 ottobre 2026 si usa lo snapshot ufficiale descritto in
`MANOPT_BASELINE.md`; i controlli della migrazione sono riportati in
`VALIDAZIONE.md`. Le storie complete comprendono anche eventuali casi
che non convergono entro il limite e gli errori espliciti del solutore.
La validazione con Manopt ufficiale e il limite Thomson aggiornato, documentata in
`matlab_thomson_50.log`, ha eseguito con successo tutte
le 18 ottimizzazioni e i
controlli delle derivate. Entrambi i metodi completano 1000 iterazioni per Grassmann,
50 per Thomson e 60 per Max-Cut. La tolleranza sul gradiente e disabilitata.
Il precedente esperimento a 100 iterazioni e documentato nel log storico
`matlab_thomson_100.log`.

La tabella seguente riporta i risultati correnti con Manopt ufficiale.

| Problema | Seed | Norma finale Armijo | Iterazioni RHZLS | Norma finale RHZLS | Errore relativo soluzione RHZLS |
|---|---:|---:|---:|---:|---:|
| Grassmann generalizzata | 1 | `3.686e-8` | 1000 | `1.374e-92` | `1.342e-92` |
| Grassmann generalizzata | 2 | `2.422e-7` | 1000 | `4.779e-15` | `3.945e-16` |
| Grassmann generalizzata | 3 | `1.391e-8` | 1000 | `1.633e-15` | `3.729e-16` |
| Thomson, tetraedro | 1 | `5.408e-9` | 50 | `1.719e-17` | `1.799e-16` |
| Thomson, tetraedro | 2 | `3.991e-9` | 50 | `1.252e-17` | `1.862e-16` |
| Thomson, tetraedro | 3 | `2.667e-11` | 50 | `8.845e-18` | `5.888e-17` |
| Max-Cut, ciclo C5 | 1 | `3.191e-8` | 60 | `4.839e-16` | `2.243e-16` |
| Max-Cut, ciclo C5 | 2 | `2.657e-8` | 60 | `6.131e-16` | `3.621e-16` |
| Max-Cut, ciclo C5 | 3 | `1.757e-8` | 60 | `1.205e-15` | `4.844e-16` |

Per Grassmann, la derivata lungo la stessa retrazione viene valutata come
`M.inner(y, M.egrad2rgrad(y, egrad), velocity)`. La velocita esatta e
orizzontale: questo e matematicamente equivalente al pairing ambientale,
ma riduce la cancellazione della componente normale vicino al minimo.
La proiezione riusa egrad, senza una seconda chiamata al gradiente.
Le differenze finite e il confronto tra i due pairing sono verificati
per tutti i seed. Coordinate nulle esattamente rappresentabili possono
produrre residui sotto eps; questo non significa accuratezza generale
oltre la precisione macchina.

Anche per Thomson e Max-Cut si usa
`M.inner(y, M.egrad2rgrad(y, egrad), velocity)`. Sui prodotti di sfere
la velocita esatta e tangente a ciascuna colonna: eliminare la componente
normale del gradiente prima del pairing preserva la derivata e riduce
la cancellazione vicino a un punto stazionario. Egrad e riutilizzato,
quindi non si aggiungono chiamate al gradiente. Lo script verifica il
pairing contro quello ambientale e contro differenze finite del costo,
oltre a tangenza della velocita e norme unitarie delle colonne.

Il limite finale di Armijo non è sempre `1e-8`: per Thomson, seed 3,
raggiunge `2.667e-11`. La differenza di accuratezza osservata in questi
casi è quindi documentata dai risultati di tutti e tre i seed, non
assunta come universale. Vicino al minimo, lo scarto del costo può
arrotondarsi a zero anche quando il gradiente e l'errore della soluzione
sono ancora distinguibili: il confronto usa tutti e tre gli indicatori.

`plot_manopt_examples.m` produce una sola figura per problema, usando il
seed 1 e le storie già salvate. Ogni figura contiene le sei curve dei
tre indicatori per Armijo e RHZLS, con colori, marker e dimensioni comuni
al resto del progetto. I file PDF sono `candidate_generalized_eigen.pdf`,
`candidate_thomson_tetrahedron.pdf` e `candidate_maxcut_cycle5.pdf`.
Le figure mostrano tutte le iterazioni registrate, senza tagliare i plateau.

Per rigenerare e visualizzare soltanto questi PDF, senza rieseguire le
ottimizzazioni:

```matlab
cd('C:/path/to/RHZLS')
startup_rhzls;
set(groot, 'defaultFigureVisible', 'on')
plot_manopt_examples;
```

Le implementazioni delle retrazioni e delle derivate sono raccolte in
`+rhzls/+retractions`: `generalized_grassmann_curve.m` e `oblique_curve.m`
richiamano, dove opportuno, la derivata polare in `stiefel_curve.m`.
Per riutilizzare un `problem` caricato da un MAT in una nuova ottimizzazione,
riassegnare `problem.M.retr = problem.curve`: MATLAB puo ricostruire due
handle anonimi numericamente equivalenti ma non identici al caricamento.
