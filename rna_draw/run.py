# import rna_draw as rd
from draw import rna_draw

# rd.rna_draw(
#     ss="(((.(((..((((........))))..((((.......))))......(((((.......)))))))).)))....",
#     seq="GCGGAUUUAGCUCAGCCGGGAGAGCGCCAGACUGGACACCUGGAGGUCCUGUGCCCGAUCCACAGAAUUCGCACCA",
#     render_type="res_type")

fig = rna_draw(ss=".(((.....))).")
fig.savefig("Yes.png")