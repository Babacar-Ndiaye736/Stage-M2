p = 1019
p = next_prime(p)
while not ((p-1)//2).is_prime():
    p = next_prime(p)
p
Fp = GF(p)
Fp(2)^((p-1)/2)
(p-1).factor()
Fp(2)^((p-1)//2)
# 2 est generateur
for i in range(p-1):
    t = ZZ(Fp(2)^i)
    tf = t.factor()
    if len(tf) > 1 and tf[-1][0] <= 5:
        i, t, tf
M = Matrix(Integers(p-1), 4, 3, [[0, 2, 1], [1,2,1], [2,2,1], [0,1,3]])
M
M.right_kernel()
Y = vector(Integers(p-1), [31,32,33, 421])
M.solve_right(Y)
Fp(2)^646
Fp(2)^1111

# Normalement on trouve p=1187 premier
#Maintenant étant donné h=314, on veut trouver x mod p-1 tel que h=g^x modulo p

h=314
for i in range(p-1):
    t = ZZ(h*Fp(2)^i)
    tf = t.factor()
    if len(tf) > 1 and tf[-1][0] <= 5:
        i, t, tf


# Resultat x=248