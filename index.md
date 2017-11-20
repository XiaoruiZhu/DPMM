---
title       : Dirichlet Process Mixture Models
subtitle    : Models and Inferences
author      : Xiaorui Zhu
job         : Business Analytics
date        : 11/12/2017
framework   : io2012        # {io2012, html5slides, shower, dzslides, ...}
highlighter : highlight.js  # {highlight.js, prettify, highlight}
hitheme     : tomorrow      # 
widgets     : [mathjax, bootstrap, quiz, shiny, interactive]
ext_widgets : {rCharts: [libraries/nvd3, libraries/highcharts]}
mode        : selfcontained  # {selfcontained, standalone, draft}
knit        : slidify::knit2slides
logo        : logo.png
fig_caption : true
---

---
## Contents

1. Stick-Breaking and Chinese Restaurant Process
2. Dirichlet Process Mixture Models
3. Conjugate Prior
3. Gibbs Sampling Algorithms
4. Simulation results

<!-- --- &radio -->

<!-- ## Who has higher creativity? -->

<!-- Who has higher creativity? -->

<!-- 1. Man -->
<!-- 2. Woman -->
<!-- 3. Engineer -->
<!-- 4. Artist -->

<!-- *** .hint -->
<!-- Creativity Diversity -->

<!-- *** .explanation -->

--- 
## Dirichlet Process

$DP$ is a random measure defined as: $$\mu = \sum^{\infty}_{n=1} p_n \delta_{\phi_n}, $$ where:

- $(p_n)_{n\in N}$ are random weights by stick-breaking construction with parameter $\theta$
- and $(\phi_n)_{n\in N} \overset{iid}{\sim} G_0$ is "the base measure". 

Therefore, $\mu \sim DP(\theta, G_0)$ has following repressentation: 

$$\begin{array}
  {rl}
  \mu & = \;  \displaystyle \sum^{\infty}_{i=1} \Big[ V_i \prod^{i-1}_{j=1}(1-V_j) \Big] \delta_{\phi_i} \\
  V_i & \overset{iid}{\sim} \; Beta(1, \theta) \\
  \phi_i & \overset{iid}{\sim} \; G_0
  \end{array}$$


---
## Stick-Breaking construction

Let $(V_n)_{n\in N}$ be i.i.d. $\text{Beta}(1,\theta)$ random variables.

That is, $P(V_1\in dx) = \theta (1 − x)^{\theta−1} \textbf{1}_{{x\in (0,1)} }dx.$

Consider:

$$\begin{array}
  {rl}
  P_1 & := \;  V_1 \\
  P_2 & := (1-V_1)V_2 \\
  P_3 & := (1-V_1)(1-V_2)V_3 \\
      & \vdots \\
  P_{n+1} & := \displaystyle V_n \prod^{n-1}_{j=1}(1-V_j)
  \end{array}$$


---
## Stick-Breaking construction

<!-- <center>![SB](figure/SB.png)</center> -->
<center><img width=800px height=700px src="figure/SB.png"></img></center>

---
## Chinese Restaurant Process (CRP)

- Imagine a Chinese restaurant that has unlimited number of tables.
- First customer sits at the first table.
- Customer $n$ sits at: 
  - Table $k$ with probability $n_k/(\alpha_0+n−1)$, where $n_k$ is the number of customers
at table $k$.
  - A new table $k + 1$ with probability $\alpha_0/(\alpha_0+n−1)$

In this metaphor, customers are analogies of integers and tables of clusters. This process can also be summarized as follows:
$$p(c_n=k|c_{1:(n-1)}) = \{ \begin{array}
  {l}
  \frac{n_k}{\alpha_0+n−1 },  \; \text{if occupied table;} \\
  \frac{\alpha_0}{\alpha_0+n−1}, \; \text{if new table} 
  \end{array} $$

---
## Animation of CRP

<!-- <center>![CRP](figure/example_1.gif)</center> -->
<center><img width=500px height=500px src="figure/example_1.gif"></img></center>


---
## Simulation of Asymptotics

Asymptotics of $K_n$: Number of clusters

- **Theorem:** $\displaystyle \text{lim}_{n\rightarrow\infty}K_n/\text{log}n = \theta$ almost surely.







