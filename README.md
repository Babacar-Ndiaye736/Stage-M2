# Enumeration

## code  pour les extensions galoisiennes

C'est dans le fichier [`Extension_galoisienne.py`](Extension_galoisienne.py)

Ce code énumère tous les polynômes d’un degré donné, dont les coefficients appartiennent à l’ensemble {0, 1, -1}. Pour chacun de ces polynômes, il vérifie s’il définit une extension galoisienne. Un polynôme satisfait cette condition s’il est irréductible modulo un certain premier p qui ne divise pas le discriminant du polynôme, ou s’il se factorise en des facteurs de même degré modulo p.

#### Lancer le code sur allmight

Pour lancer le code sur allmight on procède comme suit:

```bash

ssh allmight
cd stage-m2/sage/
sage -python Extension_galoisienne -d (ici on met le degré choisi ex:12)

```

#### Commentaire des resultats de sortie


 Les fichiers de sortie de ce code contiennent l’ensemble des polynômes , à coefficients dans {0, 1, -1}, de degré la longueur de la liste moins un, qui sont irréductibles et définissent des extensions galoisiennes. Les résultats sont présentés sous forme de listes de coefficients correspondant aux polynômes considérés.
 En conséquence, les seuls polynômes trouvés sont les polynômes cyclotomiques.
 On a pu montré que les polynômes phi_n tels que n different de q^a, 2q^a, 2, 4 avec q premier impair et a un entier naturel non nul ne sont pas utiles pour TNFS.

## code  pour les extensions cycliques


C'est dans le fichier [`Extension_cyclique.py`](Extension_cyclique.py)

Après avoir défini le test permettant de vérifier si un polynôme f à coefficient dans  {0, 1, -1} définit une extension galoisienne, le test visant à déterminer si f définit une extension cyclique consiste à vérifier d’abord que f est galoisien (c'est à dire f passe le teste de extension_galoisienne), puis à s’assurer qu’il existe un nombre premier p(parmi ceux générés jusqu’à une borne B) tel que soit f irréductible modulo p.

#### Lancer le code sur allmight

Pour lancer le code sur allmight on procède comme suit:

```bash

ssh allmight
cd stage-m2/sage/
sage -python Extension_cyclique -d (ici on met le degré choisi ex:12)

```

#### Commentaire des resultats de sortie


 Les fichiers de sortie de ce code contiennent l’ensemble des polynômes , à coefficients dans {0, 1, -1}, de degré la longueur de la liste moins un, qui sont irréductibles et définissent des extensions cycliques. Les résultats sont présentés sous forme de listes de coefficients correspondant aux polynômes considérés.
 En conséquence, les seuls polynômes trouvés sont les polynômes cyclotomiques de la forme phi_n tels que n où n=q^k, 2q^k, 2, 4 avec q premier impair et k un entier naturel non nul. Il s'avère également que ce sont les seuls polynomes utiles pour notre algorithme TNFS.

 ## code  pour les extensions ni cycliques ni galoisiennes 

C'est dans le fichier [`Extension_with_automorphismes.py`](Extension_with_automorphismes.py)


Ce code énumère tous les polynômes d’un degré donné, dont les coefficients appartiennent à l’ensemble {0, 1, -1}. Pour chacun de ces polynômes, je determine l'ensemble K des automorphismes générés. Un polynôme est candidat s’il est irréductible modulo un certain premier p qui ne divise pas le discriminant du polynôme et que le cardinal de l'ensemble K est strictement supérieur à 1 (c'est-à-dire K contient au moins un automorphisme non trivial). Le code détermine l'ensemble des  polynômes définissant des extensions ni cycliques ni galoisiennes.

#### Lancer le code sur allmight

Pour lancer le code sur allmight on procède comme suit:

```bash

ssh allmight
cd stage-m2/sage/
sage -python Extension_with_automorphismes -d (ici on met le degré choisi ex:12)

```

#### Commentaire des resultats de sortie

Les resultats obtenus contiennent plein de  polynômes palindromes définissant un seul automorphisme d'ordre 2 et d'autres polynomes définissant un seul automorphisme d'ordre 2. 
