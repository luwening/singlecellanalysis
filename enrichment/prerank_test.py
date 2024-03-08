import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt
import sys
import os

name=sys.argv[1]
gene_set=sys.argv[2]
#gene_set='KEGG_2021_Human'
#rnk = pd.read_csv(name, header=None, index_col=0, sep="\t")
pre_res = gp.prerank(rnk=name, # or rnk = rnk,
                     gene_sets=gene_set,
                     threads=4,
                     min_size=5,
                     max_size=1000,
                     permutation_num=1000, # reduce number to speed up testing
                     outdir=None, # don't write to disk
                     seed=6,
                     verbose=True, # see what's going on behind the scenes
                    )

## easy way
from gseapy import gseaplot
terms = pre_res.res2d.Term
for i in range(100):
    # to save your figure, make sure that ofname is not None
    gseaplot(rank_metric=pre_res.ranking,
            term=terms[i],
            **pre_res.results[terms[i]])

    # save figure
    #remove / from string
    gseaplot(rank_metric=pre_res.ranking, term=terms[i], ofname='%s_%s_%s_%s.png' % (os.path.basename(name),gene_set,i,'_'.join(terms[i].replace('/', '').split())), **pre_res.results[terms[i]])


from gseapy import dotplot
# to save your figure, make sure that ``ofname`` is not None
#remove GO

pre_res.res2d['Term']=[s.split('(GO:')[0] for s in pre_res.res2d['Term']]
ax = dotplot(pre_res.res2d,
             column="FDR q-val",
             title=gene_set,
             cmap=plt.cm.viridis,
             size=3, # adjust dot size
             figsize=(4,5), cutoff=0.25, show_ring=False,ofname='%s_%s_dotplot.png' % (os.path.basename(name),gene_set))

'''
from gseapy import dotplot
# to save your figure, make sure that ``ofname`` is not None
ax = dotplot(pre_res.res2d,
             column="FDR q-val",
             title='KEGG_2016',
             cmap=plt.cm.viridis,
             size=3, # adjust dot size
             figsize=(4,5), cutoff=0.25, show_ring=False)

from gseapy import enrichment_map
# return two dataframe
nodes, edges = enrichment_map(pre_res.res2d)

import networkx as nx
# build graph
G = nx.from_pandas_edgelist(edges,
                            source='src_idx',
                            target='targ_idx',
                            edge_attr=['jaccard_coef', 'overlap_coef', 'overlap_genes'])

fig, ax = plt.subplots(figsize=(8, 8))

# init node cooridnates
pos=nx.layout.spiral_layout(G)
#node_size = nx.get_node_attributes()
# draw node
nx.draw_networkx_nodes(G,
                       pos=pos,
                       cmap=plt.cm.RdYlBu,
                       node_color=list(nodes.NES),
                       node_size=list(nodes.Hits_ratio *1000))
# draw node label
nx.draw_networkx_labels(G,
                        pos=pos,
                        labels=nodes.Term.to_dict())
# draw edge
edge_weight = nx.get_edge_attributes(G, 'jaccard_coef').values()
nx.draw_networkx_edges(G,
                       pos=pos,
                       width=list(map(lambda x: x*10, edge_weight)),
                       edge_color='#CDDBD4')
plt.savefig('enrichment_map.png', dpi=300, bbox_inches='tight')

'''