""" author: Babacar Ndiaye email: babacar.ndiaye@inria.fr

Ce petit code énumère des polynomes unitaires de degré fixé donné, à
coefficients bornés en valeur absolue par une limite fixée donnée.

Exemple : degré 6, coefficients bornés par 2.

Source : projet https://gitlab.inria.fr/tnfs-alpha/alpha
fichier sage/tnfs/simul/polyselect_utils.py
https://gitlab.inria.fr/tnfs-alpha/alpha/-/blob/master/sage/tnfs/simul/polyselect_utils.py
"""
import sys
from sage.misc.functional import log # import_statements(log)
from sage.rings.integer import Integer
from sage.rings.fast_arith import prime_range
from sage.rings.finite_rings.finite_field_constructor import FiniteField
from sage.rings.number_field.number_field import NumberField

from polyselect_utils import *

QQx = QQ['x']; (x,) = QQx._first_ngens(1)

def extension_galoisienne_poly(hc, bound_p) :
    """ cette fonction teste si une extension est galoisienne
    :param           hc: polynomial 
    :param           bound_p: upper bound of prime number 
    """

    L = prime_range(bound_p) #ici on definit une borne sup à choisir 

    extension_galoisienne = True  # On suppose vrai au début
    
    
    for p in L:
        
        Fp = FiniteField(p)
        Fpz= Fp['z']
        (z,) = Fpz._first_ngens(1)

        poly = Fpz(hc)

        
        if poly.is_irreducible():
            continue  # OK, on passe au suivant
        else:
            factors = poly.factor()
            degrees = [f[0].degree() for f in factors]
            if not all(deg == degrees[0] for deg in degrees):
                extension_galoisienne = False
                break
                    
    return extension_galoisienne


def get_list_irr_poly_galois(deg, max_coeff=1, bound_p=1000, monic=True, output_file="", onthefly=False, start_counter=None, stop_counter=None, verbose=False, only_attained_max_coeff=False):
    """ get the list of (monic) irreducible univariate polynomials of degree ``deg`` and coefficients bounded by ``max_coeff``
    Moreover, only keep those which potentially have an automorphism
    because of identified splitting patterns modulo many primes p
    :param             deg: degree
    :param       max_coeff: bound on the coefficients, inclusive
    :param           bound_p: upper bound of prime number 
    :param           monic: monic poly (True/False)
    :param     output_file: filename for output, no extension (will generate .py and .gp)
    :param        onthefly: do not store in a table, only write them on-the-fly to a file (useful for deg >= 12)
    :param   start_counter: for parallel running of this function
    :param    stop_counter: for parallel running of this function
    :param         verbose: for verbose output
    :param only_attained_max_coeff: for parallel running of this function
    :returns: list of polynomials h (list of coefficients), None if onthefly

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
    
    if onthefly and len(output_file) == 0:
        print("Error please provide a filename to write the data on-the-fly")
        return
    elif onthefly:
        out_file = open(output_file, 'w+')
        if not out_file:
            print("# error opening file"+str(output_file))
            return
        print("# results written in file "+str(output_file))
        out_file.write("cyclic_tab_h_{}_{} = [ \\\n".format(deg,max_coeff))


    # estimate the number of polynomials being enumerated
    total_no_polys = number_poly(deg, max_coeff, monic)
    print("# the estimated number of polynomials to be enumerated is {} approx 2^{:.2f}".format(total_no_polys, float(log(total_no_polys, 2))))
    while counter < max_counter:
        # get the next candidate polynomial
        hc, counter = get_coeffs_from_counter(counter, deg, max_coeff, monic=monic)
        if hc is None: # the polynomial is a duplicate, or it is not irreducible, detected with a quick test
            continue
        
        mc = max([abs(ai) for ai in hc])
        # only if run in parallel or re-run:
        if only_attained_max_coeff and mc != max_coeff:
            # do not consider the polynomials with all abs(coeffs) < max_coeff
            continue

        # Lignes de codes testant si une extension est galoisiènne

        number_irr += 1 # increment counter 
        h=QQx(hc)
        
        if h.is_irreducible() :

            Test_poly=extension_galoisienne_poly(h, bound_p)

       
            hc_string = pretty_print_coeffs_from_coeffs(hc)
            h_string = pretty_print_poly_from_coeffs(hc)
            
            if Test_poly==True:
                if not onthefly:
                    tab_h.append(hc)
                elif len(output_file) > 0:
                    out_file.write("  {},".format(hc_string))
                    out_file.write(f" # Polynôme trouvé : {h_string}\n")
                    nb_polys_cyclic+=1 #compteurs de polynomes cycliques




        if verbose:
            print("    {}, # {}".format(hc_string, h))

    print("# There were {} entries for counter from {} to {}".format(number_irr, min_counter, max_counter))

    if not onthefly and len(output_file) > 0:
        write_tab_h(tab_h, deg, max_coeff, output_file, with_zeta=False)
    elif onthefly and len(output_file) > 0:
        out_file.write("]\n\n")
        out_file.flush()
        out_file.close()
        print(f"Nombre de polynomes cycliques trouvés {nb_polys_cyclic}")
        return None
    return tab_h


if __name__ == "__main__":
   
    
    args = sys.argv

    # args[0]
    if len(args) >= 3 and args[1] == "-d":
        d = Integer(args[2])
    else:
        d=4
    print("d={}".format(d))


    get_list_irr_poly_galois(deg=d, max_coeff=1, monic=True, output_file="test_list_polys_extension_galoisienne_deg"+str(d)+".py", onthefly=True, start_counter=None, stop_counter=None, verbose=False, only_attained_max_coeff=False)
