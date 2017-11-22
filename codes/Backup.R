---
  ## Simulation of Asymptotics
  
  Asymptotics of $K_n$: Number of clusters

- **Theorem:** $\displaystyle \text{lim}_{n\rightarrow\infty}K_n/\text{log}n = \theta$ almost surely.

```{r echo = F, results = 'asis', fig.width=6, fig.height=5, fig.align='center'}
load(file = "data/Them10_5.Rdata")
plot(Them10_5$ratio~Them10_5$n, type="l", xlab="Sample Size: n", ylab = "Kn/log(n)", main=expression("Simulation of "*theta[0]*" is 2"))
```

<!-- --- -->
  <!-- ## Stick-Breaking construction -->
  
  <!-- Let $(V_n)_{n\in N}$ be i.i.d. $\text{Beta}(1,\theta)$ random variables. -->
  
  <!-- That is, $P(V_1\in dx) = \theta (1 − x)^{\theta−1} \textbf{1}_{{x\in (0,1)} }dx.$ -->
    
    <!-- Consider: -->
    
    <!-- $$\begin{array} -->
    <!--   {rl} -->
    <!--   P_1 & := \;  V_1 \\ -->
    <!--   P_2 & := (1-V_1)V_2 \\ -->
    <!--   P_3 & := (1-V_1)(1-V_2)V_3 \\ -->
    <!--       & \vdots \\ -->
    <!--   P_{n+1} & := \displaystyle V_n \prod^{n-1}_{j=1}(1-V_j) -->
    <!--   \end{array}$$ -->
    
    
    <!-- --- -->
    <!-- ## Stick-Breaking construction -->
    
    <!-- <!-- <center>![SB](figure/SB.GIF)</center> --> -->
    <!-- <center><img width=800px height=700px src="figure/SB.GIF"></img></center> -->
      
      <!-- --- -->
      <!-- ## Chinese Restaurant Process (CRP) -->
      
      <!-- - Imagine a Chinese restaurant that has unlimited number of tables -->
      <!-- - First customer sits at the first table -->
      <!-- - Customer $n$ sits at:  -->
      <!--   - Table $k$ with probability $n_k/(\alpha_0+n−1)$ -->
      <!--   - A new table $k + 1$ with probability $\alpha_0/(\alpha_0+n−1)$ -->
      
      <!-- In this metaphor, customers are analogies of integers and tables of clusters. This process can also be summarized as follows: -->
      <!-- $$p(c_n=k|c_{1:(n-1)}) = \{ \begin{array} -->
          <!--   {l} -->
          <!--   \frac{n_k}{\alpha_0+n−1 },  \; \text{if occupied table;} \\ -->
          <!--   \frac{\alpha_0}{\alpha_0+n−1}, \; \text{if new table}  -->
          <!--   \end{array} $$ -->
          
              
              <!-- --- -->
              <!-- ## Simulation of Asymptotics -->
              
              <!-- - **Theorem:** Asymptotic distribution of $K_n$: $$\frac{K_n-\mathbb{E}K_n}{\sqrt{\text{Var}(K_n)}} \Rightarrow \mathcal{N}(0,1)$$ -->
              
              <!-- <center>![Theorem 10.6](figure/Them10.6.png) -->
              
              