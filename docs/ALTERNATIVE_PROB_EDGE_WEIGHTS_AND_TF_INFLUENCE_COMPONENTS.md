# Alternative ways of calculating "prob" edge weights and TF influence components

This document compares native ANANSE 0.5.1 with the custom STARE + ANANSE method, and records proposed comparisons for future testing. It covers two distinct decisions: **how a state-specific TF–gene edge is scored**, and **how the resulting network and TF differential-expression components are combined for ranking**. Changing one does not automatically correct limitations in the other.

The alternative edge formulas, expression transformations, eligibility policies and influence weights in this document remain proposals. State-specific TPM selection is now implemented for the existing multiplicative method; it is the current default. No alternative edge or influence formula has been implemented or run. The [main scoring explanation](TF_GENE_AFFINITIES_BINDING_PROBABILITIES_AND_RANKING.md) gives the STARE affinity equations and more extensive examples. The [historical audit](LABEL_IDENTITY_AUDIT_v01.md#historical-tpm-state-selection-l12) records the pooled-TPM mistake separately.

## 1. Define the inputs before comparing the formulas

For a TF $f$, target gene $g$ and cell state $s$, let:

| Symbol | Meaning | Variation across edges |
|---|---|---|
| $B^s_{f,g}$ | Scaled TF–gene binding support | TF-, gene- and state-specific |
| $T^s_f$ | Scaled TF RNA expression | Constant across the outgoing edges of that TF in a state |
| $E^s_g$ | Scaled target-gene RNA expression | Constant across the incoming edges of that gene in a state |
| $A^s_f$ | Native ANANSE scaled TF motif activity | Constant across the outgoing edges of that TF in a state |

All four are bounded between zero and one, but their transformations and underlying evidence differ. Here $A_f$ means **TF motif activity**, not the region-level ATAC activity $A_r$ in the STARE equations.

RNA abundance is not the same as differential expression. $T_f$ measures relative transcript availability within a state; the downstream `G_score` measures the focal TF's significant expression change between states. Neither directly measures protein abundance, nuclear localization, phosphorylation, cofactor availability or physical TF occupancy.

For every proposed edge formula, use only the RNA samples belonging to the network's own `cell_state`. Average replicates within that state, select TPM columns by exact sample identity, and record those columns. Native ANANSE also averages supplied expression columns/files: the caller must provide the appropriate state-specific inputs. It does not automatically separate a mixed source-and-target table. See the archived [expression loading code](../historical/software/ananse/network.py) and the [official separate-state example commands](https://github.com/vanheeringen-lab/ANANSE/blob/master/docs/examples.md).

### Scaling is part of the method

| Input | Native ANANSE 0.5.1 | Custom STARE conversion |
|---|---|---|
| Binding | Predicted region binding is aggregated to genes using TSS-distance weights; the resulting TF–gene scores are ranked with `method="min"` and min–max scaled across the combined edge set | Replicate STARE affinity matrices are aligned and averaged; raw affinities are min–max scaled **separately for each TF** across its genes |
| TF expression | Expression is log-transformed, then ranked and min–max scaled across TFs | Raw mean TPM is min–max scaled across TFs |
| Target expression | Expression is log-transformed, then ranked and min–max scaled across targets | Raw mean TPM is min–max scaled across targets |
| TF activity | Motif-based activity is ranked and min–max scaled across TFs | No separately estimated native-ANANSE activity component |

Native expression uses average ranks for ties; binding and activity use minimum ranks for ties in the inspected implementation. The monotonic log transform ordinarily preserves expression order; rank transformation, rather than the log transform alone, removes dependence on large abundance gaps. Constant inputs scale to zero.

Per-TF STARE scaling preserves relative targeting but removes uniform amplitude changes in a TF's affinity profile and comparability of absolute support between TFs. Native global binding ranks retain cross-TF ordering but are also relative to the edge universe. A one is therefore not the same biological quantity in the two methods. Switching the formula alone does not reproduce native ANANSE.

## 2. Native ANANSE averages four scaled components

With all four components available, native ANANSE computes:

$$
P_{native}=(B_{native}+T_{native}+E_{native}+A_{native})/4.
$$

The exact aggregation in the archived [network.py](../historical/software/ananse/network.py) is:

```python
columns = [
    "tf_expression",
    "target_expression",
    "weighted_binding",
    "activity",
]
columns = [col for col in columns if col in result]
result["prob"] = result[columns].mean(1).astype(np.float32)
```

Each available component has the same additive coefficient. With four components, the nominal coefficient allocation is 25 points each; with three, it is 33.3 points each. The implementation averages **available columns** and does not always divide by four.

| Component | Four-component coefficient budget | Consequence |
|---|---:|---|
| Weighted binding | 25 | Binding contributes evidence but is not a compulsory multiplicative gate |
| TF expression | 25 | Higher relative expression adds the same amount to every outgoing edge of the TF |
| Target expression | 25 | Higher relative expression adds the same amount to every incoming edge of the gene |
| TF motif activity | 25 | Region-level motif/activity evidence adds TF-wide support independently of a particular target |

The combined RNA-expression coefficients occupy 50 of 100 nominal points. These are mathematical coefficients, not measured biological or causal contributions. Real effects depend on variation, correlations, scaling and filtering. An available edge can have positive `prob` when its binding component is zero because other terms supply support.

Native binding/network construction integrates enhancer and expression evidence, but **the final mean is not a regression of RNA expression on binding**. The published method uses rank aggregation for TF–gene interaction scoring. [Xu et al., 2021](https://academic.oup.com/nar/article/49/14/7966/6318498).

## 3. What the native regression step actually estimates

The regression discussed in ANANSE occurs upstream in `ananse binding`, in [peakpredictor.py](../historical/software/ananse/peakpredictor.py), `predict_factor_activity()`. It estimates motif-associated TF activity from regional chromatin signals. It is distinct from the upstream region-binding prediction models and from the final network mean.

A conceptual regression model is:

$$
y_r\approx\alpha+\sum_m X_{r,m}\beta_m+\epsilon_r,
$$

where $X_{r,m}$ is motif evidence in region $r$, and $y_r$ is a standardized ATAC, H3K27ac or CAGE signal. Bayesian ridge regression regularizes coefficients; related motif predictors can share correlated evidence. The coefficients provide motif-activity information, not proof of direct TF binding or unique identification among TFs with similar motifs.

The key archived call is:

```python
signal = pd.DataFrame({col: scale(signal)}, index=df.index)
activity = activity.join(
    moap(
        f.name,
        genome=self.genome,
        method="bayesianridge",
        pfmfile=self.pfmfile,
        ncpus=self.ncore,
        random_state=state,
    ),
    how="outer",
    rsuffix=f"_{i}",
)
```

For each available signal column, the implementation performs three runs using up to 50,000 regions, ranks motif activities across runs/signals, averages those ranks and scales them. It assigns a TF the maximum activity among its associated motifs. The network step subsequently rank-scales that TF activity component again.

This activity term can supply chromatin-based evidence that is different from the TF's own transcript level. **It is not a fitted correction for highly expressed TFs, and it does not remove the direct `G_score` contribution in influence ranking.** The final network coefficients are still equal when four components are present. Introducing a genuine activity term into the STARE hybrid would require a defined, independently reproducible estimator; duplicating a mean STARE affinity or ATAC count as an extra component would add correlated evidence without recreating this estimator. The [upstream ANANSE source](https://github.com/vanheeringen-lab/ANANSE/blob/master/ananse/peakpredictor.py) provides the public implementation; the archived copy defines the version discussed here.

## 4. The hybrid makes expression modify binding support

The custom formula is:

$$
P_{mult}=B(0.5+0.5T)(0.5+0.5E).
$$

The existing archived calculation is:

```python
gt = float(tgt_expr_s.get(gene, 0.0))
boost_target = 0.5 + 0.5 * gt
for tf, b in bind.loc[gene].to_dict().items():
    b = float(b)
    if b <= min_binding:
        continue
    te = float(tf_expr_s.get(tf, 0.0))
    prob = b * (0.5 + 0.5 * te) * boost_target
```

### The previously decided rationale, with precise limits

Build each network directly from the state's STARE genes×TF affinity matrices, accepting aligned replicate matrices. Average replicate affinities, scale the TF profiles to zero–one, and apply light expression modulation. STARE already supplies motif-, accessibility- and Hi-C-informed gene connectivity; expression should improve plausibility while retaining that connectivity as the primary signal. Write the two columns `tf_target` and `prob`, joining TF and gene by the U+2014 em dash, for ANANSE influence.

This is a reasonable **binding-first heuristic** for preserving STARE's regulatory evidence. Both modifiers lie between 0.5 and one, so:

$$
0.25B\le P_{mult}\le B.
$$

Expression therefore attenuates binding support less for relatively expressed TFs and targets; it does not boost the score above $B$. Each expression term can change a fixed-binding edge by twofold, and the two together by fourfold. Nonzero binding remains necessary, while low expression is a soft penalty rather than a strict exclusion.

The original rationale used phrases such as “likelihood,” “expression prior,” and “mimics ANANSE's regression step.” These are biological intuitions, not literal statistical properties: this formula is not a calibrated probability, a Bayesian posterior or native ANANSE regression. ANANSE's interface accepts it, but equivalent column names do not make the inferred networks or rankings equivalent.

“Never zeroes a biologically meaningful edge” is also too strong. An expression floor does not zero a positive retained binding edge; however, per-TF min–max scaling assigns zero to the lowest affinity and to constant profiles, and `min_binding` excludes edges. Missing expression and the minimum positive measured TPM can also become scaled zero. That value should not automatically be interpreted as confirmed absence of the transcript.

### Why the method can favor highly expressed TFs

A high $T_f$ preserves more support across **all** outgoing edges. A large positive significant TF logFC subsequently supplies a high `G_score`. These are distinct abundance and change measurements, but they reuse RNA information. With state-specific TPM, relative TF induction can itself strengthen many edges, so it can affect both ranking branches. Pooling TPM did not eliminate this preference; it merely removed state-specific expression information.

High expression is not inherently a defect in a candidate. A TF with strong network and target-DE support may be biologically important. The concern is whether a TF ranks highly mainly through its own RNA signal, or whether network breadth and reused RNA information inflate its apparent support. The final component columns are needed to distinguish these possibilities. Native ANANSE also directly rewards TF expression and TF DE; its regression was not designed to solve this specific concern.

## 5. Coefficients, realized shares and local sensitivities are different

For a differentiable positive score, define proportional sensitivity:

$$
e_i=\frac{\partial\log P}{\partial\log x_i},\qquad
points_i=100\frac{e_i}{\sum_j e_j}.
$$

These points describe an infinitesimal proportional change to **already-scaled inputs**, holding the other inputs fixed. They do not attribute dataset variance, raw TPM changes after rescaling, changes across thresholds or causality. At zero, use the continuous limit; a finite increase from zero can still matter.

For native four-component addition:

$$
e_B=\frac{B}{B+T+E+A},\quad e_T=\frac{T}{B+T+E+A},\quad
e_E=\frac{E}{B+T+E+A},\quad e_A=\frac{A}{B+T+E+A}.
$$

The normalized sensitivity shares equal the realized component shares. Equal 25% coefficients do not imply equal realized shares.

| Native illustrative $(B,T,E,A)$ | $P_{native}$ | Binding points | TF-expression points | Target-expression points | Activity points |
|---|---:|---:|---:|---:|---:|
| $(0.5,0.5,0.5,0.5)$ | 0.500 | 25.0 | 25.0 | 25.0 | 25.0 |
| $(1,0.1,0.1,0.1)$ | 0.325 | 76.9 | 7.7 | 7.7 | 7.7 |
| $(0.1,1,1,0.1)$ | 0.550 | 4.5 | 45.5 | 45.5 | 4.5 |
| $(0.1,0.1,0.1,1)$ | 0.325 | 7.7 | 7.7 | 7.7 | 76.9 |
| $(0,1,1,1)$ | 0.750 | 0.0 | 33.3 | 33.3 | 33.3 |

For the multiplicative hybrid:

$$
e_B=1,\qquad e_T=\frac{T}{1+T},\qquad e_E=\frac{E}{1+E}.
$$

| Hybrid illustrative $(B,T,E)$ | $P_{mult}$ | Binding points | TF-expression points | Target-expression points |
|---|---:|---:|---:|---:|
| $(0.5,0,0)$ | 0.12500 | 100.0 | 0.0 | 0.0 |
| $(0.5,0.1,0.1)$ | 0.15125 | 84.6 | 7.7 | 7.7 |
| $(0.5,0.5,0.5)$ | 0.28125 | 60.0 | 20.0 | 20.0 |
| $(0.5,1,1)$ | 0.50000 | 50.0 | 25.0 | 25.0 |
| $(0.5,1,0)$ | 0.25000 | 66.7 | 33.3 | 0.0 |
| $(0.5,0.1,1)$ | 0.27500 | 62.9 | 5.7 | 31.4 |

Binding has the greatest local proportional sensitivity in this hybrid, but there is no universal fixed binding/RNA contribution budget. The same formula can still yield RNA-associated dominance in final rankings through cohort variation, outgoing-edge breadth and the independent `G_score` branch.

## 6. Three-component addition is a useful comparator, not a native replacement

Without an independently estimated activity term, an ANANSE-inspired alternative is:

$$
P_{additive}=(B+T+E)/3.
$$

```python
# Proposed comparator; not installed in the workflow.
prob = (binding_scaled + tf_expr_scaled + target_expr_scaled) / 3.0
```

Its coefficient budget is 33.3 binding, 33.3 TF expression and 33.3 target expression: **two thirds of the nominal additive weight is RNA availability**. This does not make it a preferred remedy for excessive RNA-driven prioritization. Its local sensitivities are $B/(B+T+E)$, $T/(B+T+E)$ and $E/(B+T+E)$; normalized points are the corresponding component shares.

| Illustrative $(B,T,E)$ | $P_{mult}$ | $P_{additive}$ | Additive binding / TF / target points | Interpretation |
|---|---:|---:|---|---|
| $(0.5,0.5,0.5)$ | 0.28125 | 0.50000 | 33.3 / 33.3 / 33.3 | Balanced inputs; the numerical ranges already differ |
| $(1,0,0)$ | 0.25000 | 0.33333 | 100 / 0 / 0 | Strong binding with expression floors |
| $(0.1,1,1)$ | 0.10000 | 0.70000 | 4.8 / 47.6 / 47.6 | Addition allows expression to compensate for weak binding |
| $(0,1,1)$ | Excluded | 0.66667 before a binding mask | 0 / 50 / 50 | Addition can supply support without binding evidence |
| $(0.8,0.2,0.2)$ | 0.28800 | 0.40000 | 66.7 / 16.7 / 16.7 | Addition gives more of its realized score to binding here |

For an initial controlled comparison, retain the same positive-binding candidate-edge mask in both variants. This prevents expression-only edges from changing the hypothesis being tested. The masked additive variant is then an explicitly **STARE-supported comparator**, not literal native ANANSE. Native global binding ranks, RNA rank scaling and motif activity remain different.

Holding $T$ and $E$ fixed across states gives additive $\Delta P=\Delta B/3$, whereas the multiplicative method attenuates $\Delta B$ by expression-dependent factors. With state-specific expression, additive changes sum as $(\Delta B+\Delta T+\Delta E)/3$; multiplicative changes interact. Both can allow an expression fall to counteract a binding increase. Neither supplies an activating/repressing sign for an edge.

## 7. Other plausible edge formulas and what each tests

These are proposals, not validated settings. Choice of transformation and common gene/TF universe must be recorded alongside the formula.

| Proposal | Formula | What it tests / limitation |
|---|---|---|
| Binding-only ablation | $P=B$ | How much expression modulation changes the network; downstream target DE and `G_score` still use RNA |
| Adjustable soft modulation | $P=B[c_T+(1-c_T)T][c_E+(1-c_E)E]$ | Larger floors weaken expression modulation; $c_T=c_E=0.5$ recovers the custom formula, while setting $c_T=1$ removes TF availability from edges |
| Target-expression-only modulation | $P=B[0.5+0.5E]$ | Avoids TF abundance multiplying all its edges; it also loses a TF-availability plausibility signal |
| Weighted addition | $P=w_BB+w_TT+w_EE$, $w_i\ge0$, $\sum w_i=1$ | Can give binding more weight; expression can still supply positive support at $B=0$ unless masked |
| Weighted geometric combination | $P=B^{w_B}T^{w_T}E^{w_E}$, positive weights summing to one | Requires joint support and has explicit log-sensitivity weights; zeros can erase real edges and exponents can increase weak positive values, so it does not preserve the custom bounds |
| Saturating expression availability | Replace expression modifiers' input by $\psi(x)=x/(x+K)$ using within-state raw TPM $x$ | Treats expression as availability that saturates instead of unlimited relative abundance; $K$ needs biological/data justification and prespecified sensitivity checks |
| Rank-scaled expression | Replace raw-TPM min–max scaling by rank–min–max scaling | Reduces dependence on expression outlier magnitudes; retains relative-abundance preferences and cohort/tie dependence |

For adjustable soft modulation, proportional TF-expression sensitivity is:

$$
e_T=\frac{(1-c_T)T}{c_T+(1-c_T)T}.
$$

For example, with $c_T=0.75$, even $T=1$ gives $e_T=0.25$, compared with 0.5 at the original floor. These floors are explicit modelling choices; 0.75 is illustrative, not an empirically optimal value.

A saturating transform should be used directly as an availability input if absolute availability is the goal. Automatically min–max scaling it again can restore dependence on the cohort's extrema. Likewise, ranking each state separately can hide a uniform expression shift. A shared transformation across states may preserve such changes but is a distinct method that needs its own definition.

### Two edge methods to test first

The recommended first comparison is **state-specific multiplicative scoring versus state-specific masked three-component addition**, holding STARE affinities, expression transformations, candidate-edge universe and influence settings fixed. This is a proposed comparison, not a decision to change the default.

1. Multiplication preserves the project's binding-first rationale and tests the intended method once TPM selection is corrected.
2. Addition tests a substantially different assumption: independent expression evidence can compensate for weak binding. It provides an ANANSE-inspired reference without pretending that native motif activity has been estimated.

Both should initially use identical $B,T,E$ so differences are attributable to the combination rule. A later comparison with native-style expression ranks must be labelled separately. If the narrower question is specifically how to suppress TF-abundance dominance, compare the multiplicative baseline with weaker or removed TF modulation; equal addition is not the recommended fix for that concern. The historical pooled-TPM output remains an audit baseline, not a third competing biological method.

## 8. Why edge changes affect influence beyond a simple mean

ANANSE influence takes the strongest selected edges from each network first, then constructs the positive target-minus-source graph:

$$
\Delta P_{f,g}=P^{target}_{f,g}-P^{source}_{f,g},\qquad \Delta P_{f,g}>0.
$$

In the study setting, selection keeps 100,000 edges per state. Absent source edges after selection count as zero. A formula change can alter top-edge membership, connectivity, the apparent gains and the strongest paths. Matching the number of retained edges does not guarantee matching their identities.

For a reachable gene $g$, the significant DE contribution is:

$$
D_g=|\log_2FC_g|\,\mathbf1(padj_g<0.05).
$$

The [archived path scorer](../historical/software/ananse/nx.py) combines positive differential weights along a selected path, with a length penalty. For a path of $L$ edges:

$$
W_{f,g}=\frac{\prod_{e\in path(f,g)}\Delta P_e}{L},\qquad
target\_score_f=\sum_g D_gW_{f,g}.
$$

The inspected influence implementation constructs an outward radius-two ego subgraph and searches within it; the path search does not separately impose a two-edge limit on paths within that subgraph. It does not sum every possible path to a gene. See [influence.py](../historical/software/ananse/influence.py).

Because paths multiply weights, even a common rescaling can affect longer paths more strongly than direct edges. Subsequent rank scaling cannot undo changes in the network's topology or path support. A large `target_score` can reflect many moderate targets or a smaller set of strongly changing, strongly connected targets. It is a sum, not a degree-normalized measure or a causal attribution.

## 9. The focal TF's RNA contributes through `G_score` and eligibility

The second component is:

$$
G\_score_f=|\log_2FC_f|\,\mathbf1(padj_f<0.05).
$$

A logFC of 4 with `padj=0.001` gives four; a logFC of 0.2 with `padj=0.2` gives zero. Greater statistical significance below the cutoff does not further increase the component at a fixed logFC.

There is also a separate eligibility gate in native ANANSE 0.5.1:

```python
tfs = [node for node in self.grn.nodes() if self.grn.out_degree(node) > 0]
de_tfs = [tf for tf in tfs if tf in self.expression_change]
de_tfs = set(tf for tf in de_tfs if self.expression_change[tf].realfc > 0)
```

Thus a small positive, non-significant TF change can enter the ranking with zero `G_score`. A TF with exactly zero or negative directional logFC cannot enter, however strong its network evidence. Increasing the weight on `target_score` alone does not bypass this rule.

This matters biologically: accessibility, chromatin contacts, cofactors, protein modification and localization can alter regulatory activity without inducing the transcript. The current network data can suggest altered regulatory opportunity; it does not directly measure these protein mechanisms. A future broadened candidate policy could retain TFs with usable expression/network evidence without requiring positive TF logFC, and report transcript direction separately. Detectability thresholds, motif ambiguity and uncertainty would need explicit rules. That is a separate departure from native eligibility, not a consequence of changing `prob` or its ranking coefficients.

Downregulated TFs with gaining network support could then be candidates for increased activity despite reduced RNA, but this does not establish which perturbation would work. Similarly, loss of regulatory support is not scored in the positive-gain graph for that direction; reverse comparisons are necessary, and they are not simply a sign flip of the final ranking.

## 10. Raw influence and final ranks combine different scales

Native raw influence is:

$$
influence\_score\_raw=target\_score+G\_score.
$$

Equal coefficients here do not imply equal numerical contributions: a sum over many genes can greatly exceed a single TF's fold-change score. The final ranking first dense-ranks and min–max scales each component across eligible TFs:

$$
R_T=\operatorname{minmax}(\operatorname{denseRank}(target\_score)),\qquad
R_G=\operatorname{minmax}(\operatorname{denseRank}(G\_score)).
$$

Dense ties share a rank; the next distinct value gets the next integer. With $K>1$ distinct values, rank $r$ becomes $(r-1)/(K-1)$; a constant component becomes zero. Native influence then computes:

$$
I=\operatorname{minmax}(\operatorname{denseRank}(R_T+R_G)).
$$

This equalizes component ranges, rather than preserving their raw magnitude imbalance. A much larger target sum does not help once its ordering is unchanged. Many zeros in `G_score` occupy one lowest dense-rank level; a small significant score can receive a substantial scaled rank if there are few distinct levels. Scores depend on the candidate cohort and tie pattern, and the final transform preserves order and ties while removing magnitude gaps.

The nominal 50/50 coefficient interpretation applies to the scaled-component sum. The target branch already contains motif, ATAC, Hi-C, state-specific expression modifiers and target DE; it is not a pure binding score. TF expression can appear indirectly in that branch and directly through `G_score`. These sources of RNA evidence are not independent.

## 11. A network-focused alternative to equal influence weights

Consider:

$$
U_\lambda=\lambda R_T+(1-\lambda)R_G,\qquad 0\le\lambda\le1,
$$

followed, if desired for ANANSE-style presentation, by:

$$
I_\lambda=\operatorname{minmax}(\operatorname{denseRank}(U_\lambda)).
$$

At $\lambda=0.5$, halving the native sum leaves its ranking and final rank-scaled values unchanged. At $\lambda=0.75$, the nominal budget is 75 points network/target-DE and 25 points focal TF-DE. This is a reasonable prespecified sensitivity scenario for prioritizing regulators whose networks help explain a transition, **not a validated optimal biological weighting**. Preserve raw components, $R_T,R_G$, $U_\lambda$, the final rank and the chosen $\lambda$ in any future outputs.

```python
# Proposed alternative; not installed.
u = lambda_target * target_score_scaled + (1 - lambda_target) * G_score_scaled
influence_alternative = minmax_scale(rankdata(u, method="dense"))
```

| Illustrative scaled components | $U_{0.5}$ | $U_{0.75}$ | What changes |
|---|---:|---:|---|
| Strong network, no direct TF-DE: $(R_T,R_G)=(1,0)$ | 0.500 | 0.750 | Network evidence receives more priority |
| Weak network, strong TF-DE: $(0.2,1)$ | 0.600 | 0.400 | Large TF induction no longer compensates as strongly |
| Both strong: $(0.8,0.9)$ | 0.850 | 0.825 | Remains well supported |
| Both midrange: $(0.5,0.5)$ | 0.500 | 0.500 | Equal components are unaffected |
| Both weak: $(0.1,0.1)$ | 0.100 | 0.100 | Still weak relative component support |
| TF-DE alone: $(0,1)$ | 0.500 | 0.250 | A potential expression marker receives less priority |

These are pre-final-transform sums, not final influence values. Final values require ranking all candidates together. At $\lambda=1$, only the target component discriminates, but its network-size bias and RNA dependence remain, and the native positive-logFC gate still applies unless independently revised.

Whether a network-focused ranking is more biologically useful depends on the experiment. For finding drivers rather than transcript markers, stronger emphasis on target/network evidence is defensible. Equal weighting remains the native comparator. Network evidence is inferred and noisy, and a broad predicted regulon can inflate `target_score`; simply increasing its coefficient is not a guarantee of improved driver prediction.

## 12. A controlled future evaluation

The first comparison should isolate edge construction: state-specific multiplicative and masked additive networks, with identical inputs, scaling and candidate support. Then compare $\lambda=0.5$ with $\lambda=0.75$ on each network. This two-by-two design separates the effects of edge scoring and final-component weighting. Keep native eligibility fixed for that comparison; any broader candidate policy belongs in a separately labelled analysis.

Useful outputs include component ranks, top-20 overlap, changes in direct/indirect targets, edge overlap and weight distributions, TF abundance versus rank, TF logFC versus rank, outdegree versus `target_score`, and stability across biological replicate resampling. Component plots can reveal expression-supported candidates; they cannot by themselves prove a technical bias or causality.

Network-size adjustment should be evaluated against a suitable null, for example networks or target sets matched for degree, expression and regulatory-region properties. Naively dividing by target number changes the biological question and can penalize true broad regulators. Where possible, perturbation-response benchmarks and matched motif-enrichment/TF-DE baselines should evaluate predictive utility. Do not select settings solely because familiar TFs rise in the ranking.

The working recommendation is to retain the binding-first method, use the implemented state-specific TPM selection, and use addition and network-weighted influence as separately labelled comparisons. If high-expression dominance remains, investigate expression transformation/modulation and candidate eligibility rather than assuming native regression or three-component addition solves it. Motif enrichment and TF DE remain useful complementary baselines; the custom hybrid's superiority must be demonstrated rather than inherited from native ANANSE benchmarks.

## References and implementation anchors

* [Main STARE affinity, expression-modulation and ranking explanation](TF_GENE_AFFINITIES_BINDING_PROBABILITIES_AND_RANKING.md).
* [Historical TPM-selection exception and offending snippet](LABEL_IDENTITY_AUDIT_v01.md#historical-tpm-state-selection-l12).
* [Archived custom converter](../historical/scripts/07_network_conversion/01_STARE_TF_Gene_Affinities_to_ananse_network.py), especially `load_tpm()` and `make_network()`.
* [Archived ANANSE network construction](../historical/software/ananse/network.py), especially `aggregate_binding()`, `_save_temp_expression()` and the final component mean.
* [Archived ANANSE motif-activity regression](../historical/software/ananse/peakpredictor.py), `predict_factor_activity()`.
* [Archived influence components and eligibility](../historical/software/ananse/influence.py), and [path weighting](../historical/software/ananse/nx.py).
* Xu et al. (2021), *ANANSE: an enhancer network-based computational approach for predicting key transcription factors in cell fate determination*. [Nucleic Acids Research](https://doi.org/10.1093/nar/gkab598).
* Hecker et al. (2023), *The adapted Activity-By-Contact model for enhancer–gene assignment and its application to single-cell data*. [Bioinformatics](https://doi.org/10.1093/bioinformatics/btad062).
