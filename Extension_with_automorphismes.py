""" author: Babacar Ndiaye email: babacar.ndiaye@inria.fr

Ce petit code énumère des polynomes unitaires de degré fixé donné, à
coefficients bornés en valeur absolue par une limite fixée donnée.

Exemple : degré 6, coefficients bornés par 2.

Source : projet https://gitlab.inria.fr/tnfs-alpha/alpha
fichier sage/tnfs/simul/polyselect_utils.py
https://gitlab.inria.fr/tnfs-alpha/alpha/-/blob/master/sage/tnfs/simul/polyselect_utils.py
"""
from sage.all_cmdline import *
from polyselect_utils import *
#from enumerate_poly import *
from sage.misc.functional import log # import_statements(log)
from sage.rings.integer import Integer
from sage.rings.fast_arith import prime_range
from sage.rings.finite_rings.finite_field_constructor import FiniteField
from sage.rings.number_field.number_field import NumberField

QQx = QQ['x']; (x,) = QQx._first_ngens(1)


# alternative : 
# from enumerate_poly import number_poly
# fonction deplacee de enumerate_poly_Babacar ici
def number_poly(deg, max_coeff, monic=True, irreducible=True):
    """
    return the number of candidate polynomials (not necessarily irreducible)
    of degree deg, whose coefficients are bounded by by max_coeff,
    of leading coefficient positive, or equals to 1 if monic

    The formula is:
    if monic: the leading coefficient is 1, there are deg other coefficients
    the constant coefficient is non-zero (otherwise it is not irreducible)
    -> (2*max_coeff) possible values for the constant coeff
    the coefficients from indice 1 to deg-2 can take any value
    -> (2*max_coeff+1)**(deg-2) possible values
    the coefficient of indice deg-1 is non-negative, because of the change of variables
    x -> -x
    if deg even: f(x) with f_deg > 0, f_{deg-1} < 0 -> f(-x) has f_deg > 0, f_{deg-1} > 0
    if deg odd: -f(-x) with f_deg > 0, f_{deg-1} < 0 -> -f(-x) has f_deg > 0, f_{deg-1} > 0

    if non-monic:
    in addition, the leading coefficient is strictly positive, because
    of the change of variables x->1/x or x-> -1/x (reverse polynomial)
    -> max_coeff possible values for the leading term
    """
    if monic and irreducible:
        return (2*max_coeff)*(2*max_coeff+1)**(deg-2)*(max_coeff+1)
    else:
        return (2*max_coeff)*(2*max_coeff+1)**(deg-2)*(max_coeff+1)*max_coeff



# Fonction pour calculer l'ordre d'un automorphisme
def ordre_automorphisme(sigma, max_iter=100):
    """ il faut mettre max_iter = degre du corps de nombres"""
    f = sigma
    for n in range(1, max_iter + 1):
        if f.is_identity():
            return n
        f = f * sigma
    raise ValueError("Ordre trop grand (>{})".format(max_iter))


def get_list_irr_poly_automorphisme(deg, max_coeff=1, monic=True, output_file="", start_counter=None, stop_counter=None, verbose=False, only_attained_max_coeff=False):
    """ get the list of (monic) irreducible univariate polynomials of degree ``deg`` and coefficients bounded by ``max_coeff``
    Moreover, only keep those which do have an automorphism
    (use K.automorphisms(), slow but exact)
    :param             deg: degree
    :param       max_coeff: bound on the coefficients, inclusive
    :param           monic: monic poly (True/False)
    :param     output_file: filename for output, no extension (will generate .py and .gp)
    :param   start_counter: for parallel running of this function
    :param    stop_counter: for parallel running of this function
    :param         verbose: for verbose output
    :param only_attained_max_coeff: for parallel running of this function
    :returns: None

    This variant does not check irreducibility over ZZ, only irreducibility
    modulo many primes p --> if it stays irreducible for at least one p, then it should be irreducible over ZZ.
    """
    #ZZy = PolynomialRing(ZZ, names=('y',)) # needed to check irreducibility
    #(y,) = ZZy._first_ngens(1)
    nb_polys_cyclic=0
    tab_h = []
    non_detected_duplicates = 0
    number_irr = 0 # number of irreducible polynomials h
    # initialise the counters, consider the input values if the script is run in parallel
    counter = 0
    max_counter = number_poly(deg, max_coeff, monic=monic)
    if start_counter != None and start_counter >= 0 and start_counter < max_counter:
        counter = start_counter
    if stop_counter != None and stop_counter >= 0 and stop_counter < max_counter:
        max_counter = stop_counter
    min_counter = counter
    
    if len(output_file) == 0:
        print("Error please provide a filename to write the data on-the-fly")
        return
    out_file = open(output_file, 'w+')
    if not out_file:
        print("# error opening file"+str(output_file))
        return
    print("# results written in file "+str(output_file))
    out_file.write("auto_tab_h_{}_{} = [ \\\n".format(deg,max_coeff))

    # estimate the number of polynomials being enumerated
    total_no_polys = number_poly(deg, max_coeff, monic)
    print("# the estimated number of polynomials to be enumerated is {} approx 2^{:.2f}".format(total_no_polys, float(log(total_no_polys, 2))))
    while counter < max_counter:
        hc, counter = get_coeffs_from_counter(counter, deg, max_coeff, monic=monic)
        if hc is None:
                continue

        if only_attained_max_coeff and max(abs(ai) for ai in hc) != max_coeff:
                continue

        #for hc in tab_h_10_1 :
        h = QQx(hc)
        if not h.is_irreducible():
            continue

        # on cherche une extension avec au moins un automorphisme non trivial de Galois mais pas forcément galoisienne
        # par contre il faut que le polynome soit irréductible modulo au moins un nombre premier
        number_irr += 1
        K = NumberField(h, names=('a',)); (a,) = K._first_ngens(1)
        AK = K.automorphisms()
        with_automorphisms=False   
        if len(AK)>1:
            # il y a un automorphisme non trivial, cherchons si h reste irreductible mod p pour un p:
            L = prime_range(1000)
            for p in L:
                Fp = FiniteField(p)
                Fpz= Fp['z']
                (z,) = Fpz._first_ngens(1)
                poly = Fpz(hc)
                if poly.is_irreducible():
                    with_automorphisms=True
                    break

        if with_automorphisms:
            hc_string = pretty_print_coeffs_from_coeffs(hc)
            h_string = pretty_print_poly_from_coeffs(hc)
            # ordre_sigma = ordre_automorphisme(AK[1])
            # regarder quels sont les ordres
            out_file.write(" {}, ".format(hc_string))
            out_file.write("# Polynôme trouvé : {}, ".format(h_string))
            out_file.write("Automorphismes d'ordre : {}\n".format([ordre_automorphisme(sigma) for sigma in AK]))
            nb_polys_cyclic+=1 #compteurs de polynomes cycliques

        if verbose:
            print("    {}, # {}".format(hc_string, h))

    print("# There were {} entries for counter from {} to {}".format(number_irr, min_counter, max_counter))

    out_file.write("]\n\n")
    out_file.flush()
    out_file.close()
    #print(f"Nombre de polynomes cycliques trouvés {nb_polys_cyclic}")
    return None

if __name__ == "__main__":
    # one example with degree 6 and coefficients 0, 1, -1
    #get_list_irr_poly_galois(deg=4, max_coeff=1, monic=True, output_file="test_list_polys_monic_cyclic_deg4.py", start_counter=None, stop_counter=None, verbose=False, only_attained_max_coeff=False)
    # one example with degree 6 and coefficients -10, ..., 10 --> 42785820 polynomials
    args = sys.argv

    # args[0]
    if len(args) >= 3 and args[1] == "-d":
        d = Integer(args[2])
    else:
        d=6
    print("d={}".format(d))
    get_list_irr_poly_automorphisme(deg=d, max_coeff=1, monic=True, output_file="List_polys_monic_with_automorphism_deg"+str(d)+".py", start_counter=None, stop_counter=None, verbose=False, only_attained_max_coeff=False)
