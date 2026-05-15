"""
phase2_gene_profile_clustering.py

Cluster all 129 network genes by a feature profile oriented so that high =
toolkit-like (peripheral, narrow pathway membership, tissue-specific, skin
top tissue, LoF-tolerant). Overlay canonical pigmentation-gene membership
(Baxter 2018 union with GWAS catalog) as ANNOTATION, not filter.

Question: when genes are clustered on their data profile, do canonical
'pigmentation' genes form one tight cluster, or do they distribute across
several? And which non-pigmentation-annotated genes share the toolkit
profile?

Outputs (in output/):
  - gene_profile_matrix.csv             oriented + percentile-ranked features
  - gene_profile_clusters.csv           cluster assignment per gene
  - figure_gene_profile_clustering.{png,pdf}
"""
import os
import gzip
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.cluster.hierarchy import linkage, fcluster, leaves_list
from scipy.spatial.distance import pdist

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR    = os.path.join(PROJECT_DIR, 'data')
OUT_DIR     = os.path.join(PROJECT_DIR, 'output')

SKIN_TISSUES = ['Skin - Sun Exposed (Lower leg)',
                'Skin - Not Sun Exposed (Suprapubic)']
FIBROBLAST_TISSUE = 'Cells - Cultured fibroblasts'

# ── Load per-gene data ─────────────────────────────────────────────────────
df = pd.read_csv(os.path.join(OUT_DIR, 'peripheral_pigmentation_pergene.csv'))
df['gene'] = df['gene'].str.upper()

# ── Canonical "pigmentation gene" reference (Baxter 2018 union GWAS) ───────
baxter = pd.read_excel(
    os.path.join(DATA_DIR, 'baxter2018_650_pigmentation_genes_tableS7.xlsx'),
    sheet_name='650 Pigmentation Genes',
)
baxter_set = set(baxter['Human gene symbol'].dropna().str.upper())
print(f"Baxter 2018 pigmentation set: {len(baxter_set)} genes")

gwas = pd.read_csv(os.path.join(DATA_DIR, 'gwas_pigmentation_associations.csv'))
gwas_genes = set()
for col in gwas.columns:
    if 'gene' in col.lower():
        gwas_genes.update(
            gwas[col].dropna().astype(str).str.upper().str.split(r'[,;\s]+').sum())
gwas_genes.discard('')
print(f"GWAS pigmentation gene mentions: {len(gwas_genes)}")

canonical = baxter_set | gwas_genes
df['canonical_pigmentation'] = df['gene'].isin(canonical)
print(f"Network genes flagged canonical: {df['canonical_pigmentation'].sum()}/{len(df)}")

# ── Feature matrix (oriented so toolkit-like = high) ───────────────────────
# Definitions:
#   centrality_inv = 1 - rank(centrality)           low centrality → high score
#   kegg_inv       = 1 - rank(kegg_pathway_count)   low KEGG cnt   → high score
#   tau_pct        = rank(tau)                      high tau       → high score
#   skin_inv       = 1 - rank(skin_rank)            rank 1         → high score
#   skin_z_pct     = rank(skin_z)                   high skin_z    → high score
#   loeuf_pct      = rank(LOEUF)                    high LOEUF     → high score
#   pli_inv        = 1 - rank(pLI) (not available — use LOEUF only)
#
# Each percentile-rank is in [0,1]; mean across columns = composite toolkit score.

def pct(s):
    return s.rank(pct=True, na_option='keep')

feat = pd.DataFrame({
    'gene':            df['gene'],
    'centrality_inv': 1 - pct(df['betweenness_centrality']),
    'kegg_inv':       1 - pct(df['kegg_pathway_count']),
    'tau_pct':        pct(df['tau']),
    'skin_top_inv':   1 - pct(df['skin_rank']),
    'skin_z_pct':     pct(df['skin_z']),
    'loeuf_pct':      pct(df['LOEUF']),
})
feat['canonical'] = df['canonical_pigmentation'].values
feat['composite'] = feat[['centrality_inv', 'kegg_inv', 'tau_pct',
                          'skin_top_inv', 'skin_z_pct', 'loeuf_pct']].mean(axis=1)

# Drop genes missing the toolkit features
profile_cols = ['centrality_inv', 'kegg_inv', 'tau_pct',
                'skin_top_inv', 'skin_z_pct', 'loeuf_pct']
feat_complete = feat.dropna(subset=profile_cols).reset_index(drop=True)
print(f"  {len(feat_complete)} genes with complete profiles")

# ── Hierarchical clustering on rows ────────────────────────────────────────
mat = feat_complete[profile_cols].values
link = linkage(pdist(mat, metric='euclidean'), method='ward')
order = leaves_list(link)
# K-cut for cluster labels
K = 4
clusters = fcluster(link, t=K, criterion='maxclust')
feat_complete['cluster'] = clusters
print(f"  cluster sizes: {dict(pd.Series(clusters).value_counts())}")

# Relabel clusters by mean composite score so cluster 1 = most toolkit-like
clust_means = feat_complete.groupby('cluster')['composite'].mean().sort_values(
    ascending=False)
relabel = {old: new for new, old in enumerate(clust_means.index, start=1)}
feat_complete['cluster'] = feat_complete['cluster'].map(relabel)

# Per-cluster summaries
print("\nPer-cluster mean profile + canonical-pigmentation count:")
for c in sorted(feat_complete['cluster'].unique()):
    sub = feat_complete[feat_complete['cluster'] == c]
    canon = sub['canonical'].sum()
    means = sub[profile_cols].mean()
    print(f"  cluster {c}: n={len(sub):3d}  canon={canon:3d}  "
          f"composite={sub['composite'].mean():.2f}")
    for col in profile_cols:
        print(f"     {col:18s} {means[col]:.2f}")

feat_complete.to_csv(os.path.join(OUT_DIR, 'gene_profile_matrix.csv'),
                     index=False)
feat_complete[['gene', 'cluster', 'composite', 'canonical']].to_csv(
    os.path.join(OUT_DIR, 'gene_profile_clusters.csv'), index=False)

# ── Cluster-1 ("toolkit-profile") members ──────────────────────────────────
toolkit_cluster = feat_complete[feat_complete['cluster'] == 1].sort_values(
    'composite', ascending=False)
print(f"\n=== Cluster 1 (toolkit-profile, n={len(toolkit_cluster)}) ===")
print(toolkit_cluster[['gene', 'canonical', 'composite']
                      + profile_cols].to_string(index=False))

# Non-canonical genes in toolkit cluster
nonpig = toolkit_cluster[~toolkit_cluster['canonical']]
print(f"\n  non-canonical-pigmentation genes in toolkit cluster: "
      f"{len(nonpig)}")
print(nonpig[['gene', 'composite']].to_string(index=False))

# ── Where do canonical pigmentation genes land? ────────────────────────────
print("\n=== Canonical pigmentation-network genes by cluster ===")
canon_sub = feat_complete[feat_complete['canonical']].sort_values(
    ['cluster', 'composite'], ascending=[True, False])
print(canon_sub[['gene', 'cluster', 'composite'] + profile_cols].to_string(
      index=False))

# ── Figure: clustered heatmap with canonical annotation ────────────────────
print("\nDrawing figure...")
# Order rows by (cluster, -composite) so the heatmap reads top-to-bottom from
# most-toolkit-like to most-hub-like.
ordered = feat_complete.sort_values(
    ['cluster', 'composite'], ascending=[True, False]).reset_index(drop=True)
heat = ordered[profile_cols].values
gene_labels = ordered['gene'].values
canon_strip = ordered['canonical'].values.astype(int)
cluster_strip = ordered['cluster'].values
composite_bar = ordered['composite'].values

n_genes = len(gene_labels)
fig = plt.figure(figsize=(12, max(20, n_genes * 0.20)))
gs = fig.add_gridspec(1, 5,
                      width_ratios=[0.10, 0.08, 1.0, 0.40, 0.08],
                      wspace=0.04)

# Y-positions where the cluster changes — used to draw separators
cluster_breaks = []
for i in range(1, n_genes):
    if cluster_strip[i] != cluster_strip[i - 1]:
        cluster_breaks.append(i - 0.5)

# (1) cluster colour strip
ax_clu = fig.add_subplot(gs[0, 0])
cluster_colors = {1: '#D94040', 2: '#E89060', 3: '#7AB6E8', 4: '#2A4878'}
clu_img = np.array([[cluster_colors[c] for c in cluster_strip]]).T
ax_clu.imshow([[0]], aspect='auto', alpha=0)  # placeholder
ax_clu.clear()
for i, c in enumerate(cluster_strip):
    ax_clu.add_patch(Rectangle((0, i - 0.5), 1, 1,
                               color=cluster_colors[c], lw=0))
ax_clu.set_xlim(0, 1); ax_clu.set_ylim(-0.5, n_genes - 0.5)
ax_clu.invert_yaxis()
ax_clu.set_xticks([]); ax_clu.set_yticks([])
ax_clu.set_title('cluster', fontsize=9, pad=4)
for sp in ax_clu.spines.values():
    sp.set_visible(False)

# (2) canonical-pigmentation strip
ax_can = fig.add_subplot(gs[0, 1])
for i, c in enumerate(canon_strip):
    if c:
        ax_can.add_patch(Rectangle((0, i - 0.5), 1, 1,
                                   color='#222222', lw=0))
ax_can.set_xlim(0, 1); ax_can.set_ylim(-0.5, n_genes - 0.5)
ax_can.invert_yaxis()
ax_can.set_xticks([]); ax_can.set_yticks([])
ax_can.set_title('canonical\npigment.', fontsize=8, pad=4)
for sp in ax_can.spines.values():
    sp.set_visible(False)

# (3) main heatmap
ax_h = fig.add_subplot(gs[0, 2])
im = ax_h.imshow(heat, aspect='auto', cmap='magma_r', vmin=0, vmax=1)
ax_h.set_xticks(range(len(profile_cols)))
short_names = {
    'centrality_inv': 'centrality\n(low→high)',
    'kegg_inv':       'KEGG cnt\n(low→high)',
    'tau_pct':        'tau\n(high→high)',
    'skin_top_inv':   'skin rank\n(1→high)',
    'skin_z_pct':     'skin z\n(high→high)',
    'loeuf_pct':      'LOEUF\n(high→high)',
}
ax_h.set_xticklabels([short_names[c] for c in profile_cols],
                     fontsize=9, rotation=30, ha='right')
ax_h.set_yticks(range(n_genes))
ax_h.set_yticklabels(gene_labels, fontsize=7)
ax_h.tick_params(axis='y', length=0, pad=2)
# bold canonical gene labels
for i, c in enumerate(canon_strip):
    if c:
        ax_h.get_yticklabels()[i].set_fontweight('bold')
        ax_h.get_yticklabels()[i].set_color('#7C1A1A')
# cluster separator lines
for b in cluster_breaks:
    ax_h.axhline(b, color='white', lw=1.5)
ax_h.set_title('Per-gene profile (percentile rank, oriented so high = toolkit-like)',
               fontsize=11, fontweight='bold', loc='left')

# (4) composite-score bar
ax_bar = fig.add_subplot(gs[0, 3], sharey=ax_h)
y = np.arange(n_genes)
bar_colors = [cluster_colors[c] for c in cluster_strip]
ax_bar.barh(y, composite_bar, color=bar_colors,
            edgecolor='black', linewidth=0.2, height=0.85)
ax_bar.axvline(0.5, color='gray', ls='--', lw=0.6)
ax_bar.set_xlim(0, 1); ax_bar.invert_yaxis()
ax_bar.set_yticks([])
ax_bar.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
ax_bar.tick_params(axis='x', labelsize=8)
ax_bar.set_xlabel('composite\ntoolkit score', fontsize=9)
for sp in ('top', 'right'):
    ax_bar.spines[sp].set_visible(False)

# (5) colorbar
cax = fig.add_subplot(gs[0, 4])
cb = fig.colorbar(im, cax=cax, orientation='vertical')
cb.set_label('percentile rank', fontsize=9)
cb.ax.tick_params(labelsize=8)

# Legend for cluster colours
from matplotlib.patches import Patch
clu_handles = [Patch(color=cluster_colors[c],
                     label=f'cluster {c}  (n={int((cluster_strip==c).sum())})')
               for c in sorted(cluster_colors)]
clu_handles.append(Patch(color='#222222', label='canonical pigmentation'))
fig.legend(handles=clu_handles, loc='lower center', ncol=5,
           fontsize=8, frameon=False,
           bbox_to_anchor=(0.5, -0.005))

fig.suptitle('Gene-profile clustering across the melanogenesis network\n'
             '(rows = genes, ordered by Ward clustering on 6 oriented features)',
             fontsize=12, fontweight='bold', y=0.995)
fig.tight_layout(rect=[0, 0.01, 1, 0.98])

for ext in ('png', 'pdf'):
    fig.savefig(os.path.join(OUT_DIR, f'figure_gene_profile_clustering.{ext}'),
                dpi=180, bbox_inches='tight', facecolor='white')
plt.close(fig)
print("  saved figure_gene_profile_clustering.png/pdf")
print("\nDone.")
