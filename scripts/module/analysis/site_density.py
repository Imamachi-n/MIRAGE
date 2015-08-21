#!/usr/bin/env python
'''
site_density.py:
Calculate target site density in 3'UTR of targetRNA.

<Reference>
[1] http://biochem218.stanford.edu/Projects%202013/Middleton.pdf
[2] Tay Y, Kats L, Salmena L, Weiss D, Tan SM, Ala U, Karreth F, Poliseno L, Provero P, Di Cunto F, Lieberman J, Rigoutsos I, Pandolfi PP.
    Coding-independent regulation of the tumor suppressor PTEN by competing endogenous mRNAs. Cell. 2011 Oct 14;147(2):344-57.

'''

def site_density_simple(MRE_number, targetRNA_length):
    site_density = MRE_number / targetRNA_length
    return site_density

def site_density_MuTaMe(pairs_key, mirna_length):
    MRE_site = []
    if len(pairs_key) == 1:
        site_density = 0
        return site_density
    #print(pairs_key)
    for line in pairs_key:
        data = line.split('||')
        st = int(data[2])
        ed = int(data[2]) + int(mirna_length) - 1
        MRE_site.append(st)
        MRE_site.append(ed)

    MRE_site.sort()
    d_mostleft = MRE_site[0]
    d_mostright = MRE_site[-1]
    MRE_distance_most = (d_mostright - d_mostleft + 1)**2
    MRE_distance_each = 0

    MRE_site_for = list(range(len(MRE_site))[::2])
    MRE_site_for.pop()
    for x in MRE_site_for:
        MRE_distance_each += (MRE_site[x+2] - MRE_site[x+1] + 1)**2

    site_density = MRE_distance_most/MRE_distance_each

    return site_density
