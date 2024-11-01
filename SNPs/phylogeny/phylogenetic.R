library("treeio")
library("ggtree")

nwk = "/mnt/hdd/dminh/Cauris/SNPs/phylogeny/fasttree/fasttree_phylogeny.nwk"
tree = read.tree(nwk)
ggtree(tree, color="firebrick", size=2, linetype="dotted")

## Not working, need to check the tree format again