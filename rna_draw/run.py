# import rna_draw as rd
from draw import rna_draw

fig =rna_draw(
    ss="(((.(((..((((........))))..((((.......))))......(((((.......)))))))).)))....",
    seq="GCGGAUUUAGCUCAGCCGGGAGAGCGCCAGACUGGACACCUGGAGGUCCUGUGCCCGAUCCACAGAAUUCGCACCA",
    render_type="res_type")
fig.savefig("Test.png")

# fig = rna_draw(ss=".(((.....))).")
# fig.savefig("Yes.png")

# rd.rna_draw(ss=".(((.....))).")