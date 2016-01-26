1) INTRODUCCIO
==============

COM FUNCIONA:

- Per compilar tots els programes:

$ sh compileAll --all

- Per compilar nomes el programa principal:

$ sh compileAll



QUE FA EL MAIN:

Donat un poliedre simplicial com a entrada, retorna les classes d'homologia.
Per a cada grau, per a cada clase, retorna una matriu (en format dispers) d'una fila
amb els coefficients de la classe, en base els simplexs de la dimensio corresponent.

Exemple:
Classes del tor:
$ ./s1 | ./times | ./main


2) FORMATS
==========


FORMAT DELS POLIEDRES SIMPLICIALS:

 - Dimensio + 1 del poliedra (ie, nombre de vertex dels simplexs)
 - Nombre de simplexs
 - Per columnes, la llista de vertexs de cada simplex.

 Exemple: (Poliedre corresponent a una S^1)

2 3
 0 1 2
 1 2 0



FORMAT DE LES MATRIUS (DISPERSES)

  Estan en format CRS (https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_row_.28CSR.2C_CRS_or_Yale_format.29)

  Concretament:

  - Nombre de files N
  - Nombre de columnes M
  - Llista R de N+1 numeros indicant els indexs de cada fila
  - Llista C de R[N+1] numeros indicant la columna de cada valor no nul
  - Llista V de R[N+1] numeros indicant valor de cada valor no nul

Exemple:

3 4
0 2 4 5
1 3 0 3 2
-1 5 -9 7 4

Correspon a la matriu

 0 -1  0  5
-9  0  0  7
 0  0  4  0



FORMAT DELS COMPLEXOS SIMPLICIALS

  - Un numero N indicant la dimensio + 1 del complex
  - Per a cada grau, desde 0 fins a N-1, el poliedre simplicial dels simplex
  - Per a cada grau, desde 1 fins a N-1, la matriu (dispersa) de l'operador vora
  - Un nombre que indica nombre de simplexs de dimensio maxima M.
  - M valors (1's o -1's) indicant l'orientancio dels simplex maximals

Exemple: (Complex simplicial corresponent al poliedre de una S^1)

2
1 3
 0 1 2
2 3
 0 0 1
 1 2 2
3 3
 0 2 4 6
 0 1 0 2 1 2
 -1 -1 1 -1 1 1
3
 1 -1 1


3) PROGRAMES
============


FUNCIONS AMB POLIEDRES

$ ./s1 [N]
  retorna la triangulacio de una S^1 triangulada en N troc,os. Per defecte, n = 3.

$ ./simplex N
  Retorna un simplex de dimensio N.

$ ./boundary
  Espera com a entrada un poliedre. Retorna la seva vora (orientada).

  Exemple:

  Triangulacio de l'esfera:
  $ ./simplex 3 | ./boundary

$ ./cone
  Espera com a entrada un poliedre. Retorna el seu con.

$ ./suspension
  Espera com a entrada un poliedre. Retorna la seva suspensio.

  Exemple:

  Triangulacio de l'esfera (octaedre):
  $ ./s1 4 | ./supension

$ ./times [-d]
  Espera com a entrada dos poliedres, i retorna el poliedre producte cartesia.
  Si te l'opcio -d, nomes espera un poliedre i retorna el producte per ell
  mateix.

  Exemples:

  Tor de 18 simplexs:
  $ ./s1 | ./times -d

  Tor de 40 simplexs:
  $ (./s1 4 ; ./s1 5) | ./times

$ ./sum [-d]
  Espera com a entrada dos poliedres, i retorna el poliedre de la suma connexa,
  fet per l'ultim simplex del primer poliedre i el primer simplex del segon.

  Si te l'opcio -d, nomes espera un poliedre i retorna la suma amb una copia d'ell
  mateix.

  Nota: Fer una suma connexa vol dir identificar dos simplexs i eliminarlos. Per
  tant, fer una suma d'un poliedre P amb un sol simplex equival a eliminar
  l'ultim simplex de P. 

  Exemples:

  Superficie de genere 2:
  $ ./s1 | ./times -d | ./sum -d

  Tor de 18 + 3 = 21 simplexs:
  $ (./s1 | ./times -d ; ./simplex 3 | ./boundary ) | ./sum

  Tor amb un forat:
  $ (./s1 | ./times -d ; ./simplex 2) | ./sum