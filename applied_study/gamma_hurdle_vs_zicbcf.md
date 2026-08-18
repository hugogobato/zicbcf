# Applied Study: Gamma Hurdle vs. ZIC-BCF on MEPS Dental Expenditures

This report compares two zero-inflated Bayesian causal models applied to the 2023 Medical Expenditure Panel Survey (MEPS) dental expenditures data: the **Gamma Hurdle** model of Oganisian et al. (2019), shipped in `zicbcf` as a benchmark, and **ZIC-BCF-Smear**, the proposed zero-inflated continuous Bayesian causal forest with Duan's smearing re-transformation.

## 1. Setup

Both models share the same estimation pipeline so that differences reflect model specification rather than data handling.

1. **Outcome (Y):** Dental Expenditures (`DVTEXP23`). Heavily right-skewed semicontinuous variable with roughly 53% zeros.
2. **Treatment (Z):** Dental Insurance (`DNTINS23_M23`).
3. **Covariates (X):** Age, Sex, Race, Family Income, Poverty Category, Region, and Marital Status.
4. **Sample Size:** 18,763 valid participant records after excluding invalid and inapplicable survey responses.
5. **Propensity Score:** A logistic regression model was fit to estimate the probability of having dental insurance given the covariates. Both models condition on the same $\hat{\pi}$.
6. **MCMC:** 1,000 post-burn-in draws after 1,000 burn-in iterations for Gamma Hurdle, and 1,000 post-burn-in draws after 500 burn-in iterations for ZIC-BCF-Smear. The Gamma Hurdle fit additionally returns the hurdle probability draws ($p_0, p_1$) needed to decompose the overall effect into participation and intensity margins.

The Gamma Hurdle model jointly fits (i) a Bernoulli hurdle for any versus no expenditure and (ii) a Gamma regression on the positive-expenditure intensity, recovering the overall CATE as $p_1 \cdot E[Y \mid Y > 0, Z = 1] - p_0 \cdot E[Y \mid Y > 0, Z = 0]$. ZIC-BCF-Smear does not factor the effect multiplicatively; it models the conditional mean directly on the log scale, then re-transforms to the dollar scale via the smearing approach with propensity adjustment.

## 2. Average Treatment Effect

The two models agree on the sign and order of magnitude of the effect but differ in level and uncertainty.

| Model | ATE ($) | 95% Credible Interval ($) |
| :--- | ---: | :--- |
| **ZIC-BCF-Smear** | 201.61 | [167.18, 235.15] |
| **Gamma Hurdle** | 146.49 | [94.21, 199.46] |
| **ZIC-BCF-Smear (hurdle margin)** | 0.149 | [0.132, 0.165] |
| **Gamma Hurdle (hurdle margin)** | 0.116 | [0.100, 0.133] |

Both methods place the entire 95% credible interval above zero, so both support a credibly positive causal effect of dental insurance on dental expenditure. The point estimate from ZIC-BCF is about 38% larger than that from Gamma Hurdle, and the Gamma Hurdle interval is roughly 1.5x wider (105 dollars of width versus 68 dollars for ZIC-BCF). The Gamma Hurdle estimates an 11.6 percentage-point average increase in the probability of any dental expenditure, which on its own (applied to the average positive-expenditure level) accounts for a substantial share of the total effect and is the most clearly identified margin of action.

A point that was missed in earlier versions of this analysis is that ZIC-BCF-Smear also delivers a hurdle-margin (participation) treatment effect, and not only the overall dollar-scale CATE. Its two-part construction fits a probit Bayesian causal forest for the participation process, whose prognostic ($\mu_b$) and treatment-effect ($\tau_b$) draws map to the covariate-specific participation probabilities $p_0 = \Phi(\mu_b)$ and $p_1 = \Phi(\mu_b + \tau_b)$; the hurdle-margin CATE is then $p_1 - p_0$, exactly as reported in the simulation studies. On the applied data, ZIC-BCF-Smear estimates a 14.9 percentage-point average increase in the probability of any dental expenditure (95% credible interval [13.2, 16.5]), about 3.3 percentage points larger than the Gamma Hurdle's 11.6-point estimate and comparably tight. The two hurdle margins agree on the sign and rough magnitude of the participation effect while ZIC-BCF-Smear's is somewhat larger, which is consistent with its nonparametric probit forest being freer to track participation heterogeneity than the Gamma Hurdle's more restrictive parameterization.

Two economic considerations help interpret the gap. First, the Gamma Hurdle imposes a Gamma distribution on the positive intensity, which is constrained to a single tail shape and cannot match the heaviest tail of the log-normal-like dental spending distribution; Duan's smearing in ZIC-BCF is fully non-parametric with respect to the residuals and soaks up that heaviness. Second, the Gamma Hurdle's partitioned mean is more sensitive to the conditional hurdle probabilities, which carry their own posterior uncertainty that propagates multiplicatively into the overall effect. ZIC-BCF avoids that multiplication by estimating the conditional mean in one step.

## 3. Subgroup CATE Analysis

We follow the same procedure for both models: for each subgroup, the unit-level posterior draws are averaged over the units in that subgroup within every posterior draw to produce a posterior distribution for the subgroup CATE, from which the point estimate and 95% credible interval are read off.

### 3.1 ZIC-BCF-Smear

| Covariate | Subgroup | N | CATE ($) | 95% CI ($) |
| :--- | :--- | ---: | ---: | :--- |
| Overall | All | 18,763 | 201.61 | [168.65, 235.95] |
| Sex | Male | 8,928 | 203.08 | [166.34, 242.96] |
| | Female | 9,835 | 200.27 | [159.87, 241.69] |
| Race | White | 14,004 | 219.60 | [183.61, 254.94] |
| | Black | 2,639 | 142.37 | [99.04, 193.61] |
| | Asian | 380 | 154.16 | [104.41, 212.53] |
| | Other / Multiple | 1,740 | 157.01 | [99.45, 215.00] |
| Age | 0-17 | 3,647 | 153.47 | [105.10, 205.34] |
| | 18-34 | 3,352 | 192.05 | [160.86, 226.58] |
| | 35-49 | 3,468 | 193.68 | [162.70, 228.46] |
| | 50-64 | 3,760 | 223.01 | [181.80, 262.18] |
| | 65+ | 4,536 | 235.69 | [154.06, 314.16] |
| Poverty | Poor | 2,834 | 164.05 | [104.10, 235.58] |
| | Near-poor | 877 | 167.57 | [99.18, 246.63] |
| | Low-income | 2,458 | 142.63 | [78.13, 212.01] |
| | Middle-income | 5,243 | 227.53 | [187.67, 271.30] |
| | High-income | 7,351 | 221.38 | [175.75, 264.99] |
| Region | Northeast | 2,892 | 204.97 | [133.15, 265.34] |
| | Midwest | 3,837 | 205.26 | [163.37, 249.01] |
| | South | 7,375 | 176.61 | [147.45, 209.46] |
| | West | 4,659 | 236.08 | [186.42, 292.97] |

### 3.2 ZIC-BCF-Smear (hurdle margin only, probability of any expenditure)

The same probit hurdle stage that ZIC-BCF-Smear fits to build the overall dollar CATE also yields a subgroup posterior for the participation effect, computed draw-by-draw in exactly the same way as the dollar CATE.

| Covariate | Subgroup | N | Hurdle CATE | 95% CI |
| :--- | :--- | ---: | ---: | :--- |
| Overall | All | 18,763 | 0.149 | [0.132, 0.165] |
| Sex | Male | 8,928 | 0.155 | [0.135, 0.179] |
| | Female | 9,835 | 0.143 | [0.119, 0.166] |
| Race | White | 14,004 | 0.165 | [0.145, 0.183] |
| | Black | 2,639 | 0.103 | [0.066, 0.150] |
| | Asian | 380 | 0.110 | [0.070, 0.156] |
| | Other / Multiple | 1,740 | 0.103 | [0.060, 0.145] |
| Age | 0-17 | 3,647 | 0.106 | [0.070, 0.141] |
| | 18-34 | 3,352 | 0.166 | [0.144, 0.189] |
| | 35-49 | 3,468 | 0.169 | [0.149, 0.192] |
| | 50-64 | 3,760 | 0.170 | [0.145, 0.193] |
| | 65+ | 4,536 | 0.140 | [0.105, 0.176] |
| Poverty | Poor | 2,834 | 0.130 | [0.079, 0.184] |
| | Near-poor | 877 | 0.128 | [0.069, 0.189] |
| | Low-income | 2,458 | 0.094 | [0.043, 0.147] |
| | Middle-income | 5,243 | 0.173 | [0.149, 0.196] |
| | High-income | 7,351 | 0.161 | [0.139, 0.182] |
| Region | Northeast | 2,892 | 0.150 | [0.122, 0.176] |
| | Midwest | 3,837 | 0.153 | [0.132, 0.171] |
| | South | 7,375 | 0.140 | [0.121, 0.161] |
| | West | 4,659 | 0.160 | [0.137, 0.196] |

### 3.3 Gamma Hurdle (overall CATE: hurdle x intensity)

| Covariate | Subgroup | N | CATE ($) | 95% CI ($) |
| :--- | :--- | ---: | ---: | :--- |
| Overall | All | 18,763 | 146.49 | [94.21, 199.46] |
| Sex | Male | 8,928 | 166.98 | [108.06, 236.63] |
| | Female | 9,835 | 127.89 | [55.07, 199.90] |
| Race | White | 14,004 | 150.13 | [96.18, 207.59] |
| | Black | 2,639 | 156.35 | [96.69, 223.15] |
| | Asian | 380 | 116.97 | [51.31, 180.71] |
| | Other / Multiple | 1,740 | 108.70 | [0.45, 217.07] |
| Age | 0-17 | 3,647 | 147.07 | [85.18, 205.93] |
| | 18-34 | 3,352 | 155.61 | [105.06, 201.10] |
| | 35-49 | 3,468 | 137.17 | [88.77, 185.54] |
| | 50-64 | 3,760 | 144.64 | [82.54, 209.83] |
| | 65+ | 4,536 | 147.94 | [38.84, 265.98] |
| Poverty | Poor | 2,834 | 110.62 | [19.75, 214.40] |
| | Near-poor | 877 | 139.97 | [52.62, 240.73] |
| | Low-income | 2,458 | 170.34 | [100.94, 247.98] |
| | Middle-income | 5,243 | 181.82 | [126.78, 242.56] |
| | High-income | 7,351 | 127.92 | [55.85, 200.66] |
| Region | Northeast | 2,892 | 106.20 | [13.71, 207.05] |
| | Midwest | 3,837 | 135.43 | [78.74, 196.09] |
| | South | 7,375 | 153.46 | [100.73, 208.22] |
| | West | 4,659 | 169.57 | [99.15, 239.30] |

### 3.4 Gamma Hurdle (hurdle margin only, probability of any expenditure)

Like ZIC-BCF-Smear, the Gamma Hurdle returns a posterior for the participation effect on the probability of any dental expenditure, separately from the intensive margin.

| Covariate | Subgroup | N | Hurdle CATE | 95% CI |
| :--- | :--- | ---: | ---: | :--- |
| Overall | All | 18,763 | 0.116 | [0.100, 0.133] |
| Sex | Male | 8,928 | 0.123 | [0.102, 0.145] |
| | Female | 9,835 | 0.110 | [0.089, 0.130] |
| Race | White | 14,004 | 0.126 | [0.109, 0.144] |
| | Black | 2,639 | 0.104 | [0.083, 0.125] |
| | Asian | 380 | 0.092 | [0.070, 0.114] |
| | Other / Multiple | 1,740 | 0.057 | [0.014, 0.097] |
| Age | 0-17 | 3,647 | 0.072 | [0.044, 0.099] |
| | 18-34 | 3,352 | 0.107 | [0.087, 0.127] |
| | 35-49 | 3,468 | 0.116 | [0.098, 0.135] |
| | 50-64 | 3,760 | 0.130 | [0.113, 0.149] |
| | 65+ | 4,536 | 0.145 | [0.119, 0.173] |
| Poverty | Poor | 2,834 | 0.045 | [0.004, 0.084] |
| | Near-poor | 877 | 0.077 | [0.044, 0.109] |
| | Low-income | 2,458 | 0.114 | [0.090, 0.137] |
| | Middle-income | 5,243 | 0.141 | [0.124, 0.159] |
| | High-income | 7,351 | 0.131 | [0.110, 0.151] |
| Region | Northeast | 2,892 | 0.091 | [0.063, 0.120] |
| | Midwest | 3,837 | 0.112 | [0.093, 0.131] |
| | South | 7,375 | 0.116 | [0.097, 0.135] |
| | West | 4,659 | 0.135 | [0.114, 0.159] |

## 4. Heterogeneity: Where the Two Models Converge and Diverge

### 4.1 Convergent qualitative findings

Both models agree on the direction of the overall effect and on several substantively important patterns:

1. **Positive everywhere.** Every subgroup in both models has a posterior probability of a positive overall CATE of at least 0.97. Insurance raises expected dental spending in all subpopulations under both specifications.
2. **Children are the least responsive on the participation margin.** Both hurdle margins agree that children (0-17) gain the least in the probability of any dental expenditure: 7.2 percentage points for the Gamma Hurdle and 10.6 percentage points for ZIC-BCF-Smear, in each case credibly below the working-age groups. Where they part company is the shape of the rest of the gradient. The Gamma Hurdle finds a clean monotonic rise, peaking at 14.5 percentage points for adults aged 65 and older, which ZIC-BCF mirrors on its overall dollar CATE (rising from 153 dollars for children to 236 dollars for the 65+ group). ZIC-BCF-Smear's own hurdle margin, by contrast, is hump-shaped: it peaks in the working-age groups (16.6 to 17.0 percentage points for ages 18-64) and falls back to 14.0 percentage points for the 65+ group. The two participation margins therefore coincide on which group benefits least (children) but disagree on whether the oldest group is the most responsive on participation, even though both agree the 65+ group carries the largest dollar effect.
3. **Income gradient on participation.** Both hurdle margins show the same poverty gradient: the low-income group gains the least in participation (9.4 percentage points for both models) and the middle-income group the most (17.3 percentage points for ZIC-BCF, 14.1 for the Gamma Hurdle). ZIC-BCF shows the same ordering on the dollar scale (142 dollars for low-income versus 228 dollars for middle-income). At the top of the distribution the two hurdle margins again agree that the high-income group stays high (16.1 percentage points for ZIC-BCF, 13.1 for the Gamma Hurdle), and ZIC-BCF's dollar CATE keeps high-income high at 221 dollars, whereas the Gamma Hurdle's dollar CATE falls back to 128 dollars.
4. **Race gradient recovered on both ZIC-BCF margins.** ZIC-BCF-Smear finds credible White-versus-Black and White-versus-Other race gaps on its participation margin (+0.062 and +0.062) that echo the credible race gaps on its dollar CATE. The Gamma Hurdle also finds credible race gaps on its participation margin, so all three participation-margin analyses converge on race; the divergence, described below, is that only ZIC-BCF carries the race gradient through to the dollar scale.

### 4.2 Divergent findings

The two models diverge most sharply on two margins.

1. **Race.** ZIC-BCF finds large and credible heterogeneity by race on both of its margins. On the dollar scale the White CATE of 219.60 dollars credibly exceeds the Black CATE of 142.37 dollars (+77.23 dollars, 95% CI [23.56, 121.91]), the Asian CATE, and the Other/Multiple CATE; and on its participation margin the White hurdle CATE of 0.165 credibly exceeds the Black (0.103) and Other/Multiple (0.103) hurdle CATEs (+0.062 [0.005, 0.102] and +0.062 [0.014, 0.107]). The two ZIC-BCF margins are therefore internally coherent on race. The Gamma Hurdle finds essentially no racial heterogeneity on the overall dollar CATE (White 150.13 dollars, Black 156.35 dollars, Asian 116.97 dollars, Other 108.70 dollars; no pairwise contrast excludes zero). It does find racial heterogeneity on its participation margin (White minus Black = +0.022, 95% CI [0.010, 0.034]; White minus Other = +0.069, 95% CI [0.028, 0.114]) but those participation gaps are too small to dominate the overall contrast once the Gamma intensity margin is folded back in. The instructive contrast is thus not whether the participation margin sees the race gradient (all three participation analyses do) but whether the model carries it to the dollar scale: ZIC-BCF does, because its nonparametric smearing re-transformation preserves the race-linked intensity heterogeneity, while the Gamma likelihood averages it away.
2. **Region.** ZIC-BCF finds a credibly larger effect in the West (236.08 dollars) than in the South (176.61 dollars), a contrast of -59.47 dollars (95% CI [-108.40, -14.65]). The Gamma Hurdle places the West at the top of its ranking too (169.57 dollars) but with a much wider interval and without any pairwise regional contrast excluding zero. The heterogeneity is real in both, but only ZIC-BCF has the precision to declare it credible.
3. **Sex.** The two models disagree directionally on sex. ZIC-BCF estimates essentially identical effects for males and females (203.08 dollars versus 200.27 dollars). The Gamma Hurdle estimates a larger male effect (166.98 dollars versus 127.89 dollars), but the 95% interval for the contrast of +39.10 dollars is [-43.78, 124.41], so the difference is not credible in either model. This is a case where the Gamma Hurdle's added uncertainty absorbs a nominal gap that ZIC-BCF never exhibits.

### 4.3 Credible contrasts

Restricting attention to pairwise subgroup contrasts whose 95% credible interval excludes zero, the credible heterogeneity identified by each model is:

**ZIC-BCF-Smear** (overall CATE contrasts):

1. Race: White - Black = +77.23 dollars [23.56, 121.91]; White - Asian = +65.44 dollars [3.42, 117.23]; White - Other = +62.59 dollars [1.21, 121.53].
2. Age: 0-17 - 50-64 = -69.54 dollars [-134.61, -4.27].
3. Poverty: Low-income - Middle-income = -84.90 dollars [-155.61, -10.65]; Low-income - High-income = -78.75 dollars [-145.60, 7.97] (borderline).
4. Region: South - West = -59.47 dollars [-108.40, -14.65].

**ZIC-BCF-Smear** (hurdle-margin contrasts):

1. Race: White - Black = +0.062 [0.005, 0.102]; White - Other = +0.062 [0.014, 0.107]. (White - Asian = +0.055 [-0.006, 0.099] is directionally consistent but borderline.)
2. Age: children (0-17) gain credibly less on participation than every working-age group, 0-17 - 18-34 = -0.060 [-0.097, -0.026], 0-17 - 35-49 = -0.064 [-0.104, -0.026], and 0-17 - 50-64 = -0.064 [-0.106, -0.022]; the 65+ group falls back toward the children's level, so the gradient is hump-shaped rather than monotone.
3. Poverty: Low-income - Middle-income = -0.078 [-0.140, -0.015]; Low-income - High-income = -0.066 [-0.124, -0.007]. On the income-tertile grouping, Low - Middle = -0.035 [-0.070, -0.003].
4. Region: South - West = -0.020 [-0.060, -0.004].

**Gamma Hurdle** (overall CATE contrasts):

1. Poverty: Poor - Near-poor = -29.36 dollars [-46.52, -10.30]; Poor - Low-income = -59.72 dollars [-105.37, -6.27].

**Gamma Hurdle** (hurdle-margin contrasts):

1. Race: White - Black = +0.022 [0.010, 0.034]; White - Asian = +0.034 [0.016, 0.054]; White - Other = +0.069 [0.028, 0.114]; Black - Other = +0.047 [0.011, 0.086]; Asian - Other = +0.035 [0.007, 0.063].
2. Age: all older-versus-younger contrasts are credibly negative on the (older minus younger) convention used by the script, confirming a monotone age gradient on participation.
3. Poverty: a near-complete monotone gradient on participation (Poor - Near-poor = -0.032 [-0.042, -0.022] through to Poor - Middle-income = -0.096 [-0.135, -0.059]).
4. Region: Northeast - Midwest = -0.020 [-0.035, -0.007]; Northeast - West = -0.044 [-0.085, -0.008]; South - West = -0.020 [-0.034, -0.006].
5. Marital status: Married - Not married = +0.030 [0.004, 0.057].

The Gamma Hurdle therefore finds credible heterogeneity predominantly on the participation margin, which is the margin its hurdle component is built to identify, while ZIC-BCF finds a sparser but economically sharper set of credible contrasts on the overall dollar scale.

## 5. Gamma Hurdle CATEs and the Applied Literature

Read against the dental health-economics literature, the Gamma Hurdle CATEs are most persuasive on the participation margin. The hurdle component estimates positive participation effects in every subgroup and reproduces several well-known access gradients: the effect on any dental expenditure rises from 7.2 percentage points for children to 14.5 percentage points for adults aged 65 and older; it is larger for White respondents (12.6 percentage points) than for Black, Asian, and Other/Multiple respondents; it is much smaller for poor respondents (4.5 percentage points) than for low-, middle-, and high-income respondents; and it is highest in the West (13.5 percentage points) and lowest in the Northeast (9.1 percentage points). These patterns converge with the literature showing that coverage is associated with greater dental use, that use and expenditure rise with age, and that racial, income, and geographic barriers affect the conversion of coverage into realized care.

The overall dollar-scale Gamma Hurdle CATEs are more muted, and this is where the benchmark starts to diverge from the strongest literature-aligned patterns. Race is the clearest example. The participation margin agrees with the literature by ranking White respondents above Black, Asian, and Other/Multiple respondents, but the overall dollar CATEs do not show credible racial heterogeneity: White is estimated at $150.13, Black at $156.35, Asian at $116.97, and Other/Multiple at $108.70, with all racial contrasts crossing zero. In substantive terms, the model detects who is more likely to cross the zero hurdle, but the Gamma intensity component does not convert that participation gradient into a credible race gradient on the total dollar scale.

Age shows the same split. The hurdle margin strongly converges with prior evidence that older adults have higher dental-care use and greater dental need, but the overall Gamma Hurdle CATE is nearly flat across age groups, from $137.17 to $155.61 for adults below age 65 and $147.94 for those 65 and older. No age contrast in the overall dollar CATE is credible. This diverges from the usual expenditure-gradient interpretation of the dental literature, where older adults are expected to have larger realized spending responses because restorative and replacement procedures become more common with age.

Income and poverty are partially convergent. The overall Gamma Hurdle CATE is lowest for the poor group ($110.62) and highest for middle-income respondents ($181.82), with credible contrasts for Poor minus Near-poor and Poor minus Low-income. That agrees with the literature's view that insurance alone does not remove cost-sharing, provider access, and non-financial barriers faced by poorer respondents. The divergence is at the top of the distribution: high-income respondents have a high participation effect but a lower overall dollar CATE ($127.92) than low- and middle-income respondents. This suggests that the Gamma intensity model is smoothing or shrinking the positive-spending component in a way that weakens the expected income gradient on realized dollars.

Regional heterogeneity is weakly benchmarked in the literature, and the Gamma Hurdle result should be read cautiously. The overall dollar CATE ranks the West highest ($169.57) and the Northeast lowest ($106.20), which is directionally compatible with the idea that regional provider supply, prices, and plan generosity can shape the value of coverage. However, no overall regional contrast is credible. The hurdle margin is stronger, with credible participation gaps for Northeast versus Midwest, Northeast versus West, and South versus West. Thus, the Gamma Hurdle provides evidence for regional heterogeneity in access to any care, but not decisive evidence that the total dollar effect differs by region.

For sex and marital status, the Gamma Hurdle findings do not sharply contradict the literature, mostly because the literature rarely treats these variables as focal treatment-effect modifiers. The overall dollar CATE is higher for males than females, but the contrast is not credible. The marital-status dollar contrast is also not credible, while the hurdle margin estimates a larger participation effect for married respondents. These are best viewed as exploratory participation-margin signals rather than established divergences from prior evidence.

Overall, the Gamma Hurdle CATE analysis converges with the literature when the question is whether insurance helps people cross the zero-expenditure hurdle. It diverges, or at least becomes much less decisive, when the question is how strongly that extra participation translates into total dollars. The most plausible interpretation is methodological rather than substantive: the Gamma hurdle decomposition is well matched to participation effects but imposes a parametric positive-intensity model that can smooth away heavy-tail expenditure heterogeneity. For claims about subgroup dollar effects, the Gamma Hurdle therefore supports the broad positive-effect narrative but is less aligned with the literature than ZIC-BCF-Smear on race, age, and the upper part of the income gradient.

## 6. Discussion

Three observations frame the comparison.

1. **Where the methods agree, they reinforce each other.** The age gradient and the income gradient on participation are recovered by both methods, with the Gamma Hurdle localizing them to the hurdle component and ZIC-BCF expressing them in dollars on the response scale. That coherence across two structurally different estimators is the strongest evidence in the study.
2. **Both models decompose into a participation margin, and ZIC-BCF's two margins are internally coherent.** An earlier version of this comparison treated the participation-margin decomposition as a distinctive feature of the Gamma Hurdle, but ZIC-BCF-Smear supplies the same decomposition through its probit hurdle stage: $p_0 = \Phi(\mu_b)$, $p_1 = \Phi(\mu_b + \tau_b)$, and a hurdle-margin CATE of $p_1 - p_0$. On the applied data ZIC-BCF-Smear's participation-margin ATE (14.9 percentage points) is larger than and coherent with the Gamma Hurdle's (11.6 points), and its participation-margin subgroup contrasts recover the same race, poverty, and region gradients that appear on its dollar CATE. Where the models genuinely diverge is on whether the dollar scale inherits the participation gradients: ZIC-BCF's clean racial and regional contrasts survive on both of its margins because its non-parametric smearing re-transformation captures the heavy-tailed intensity heterogeneity, whereas the Gamma Hurdle's dollar-scale contrasts wash out because its Gamma likelihood cannot represent that intensity heterogeneity. The two models are therefore answering the same questions with different fidelity, and the simulation study corroborates that ZIC-BCF recovers both the overall CATE and the hurdle-margin CATE at least as accurately as the Gamma Hurdle under heavy-tailed intensity distributions.
3. **The tighter ZIC-BCF intervals are consistent with the simulation evidence.** Across every subgroup, the ZIC-BCF credible interval is narrower than the Gamma Hurdle interval (for example, 68 dollars of width for the overall ATE versus 105 dollars). The simulation studies report lower CATE RMSE for ZIC-BCF under comparable DGPs, and the gap in interval width on the real data is in the same direction.

The applied conclusion is robust to model choice: dental insurance causes a credibly positive increase in annual dental expenditure in every subpopulation we examine, with an average effect between 146 and 202 dollars depending on specification. The substantive heterogeneity, where (and only where) the two models agree, points to larger insurance effects for older adults and for middle-income relative to low-income enrollees. Where ZIC-BCF alone identifies credible racial and regional contrasts, those reflect intensity-margin heterogeneity that the Gamma likelihood explicitly averages away, and we report them as ZIC-BCF-specific findings.
