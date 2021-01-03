## Best approximation

### Approximation

Approximating a set of data or a function in $[a, b]$ consists in finding a suitable function $f$ that represents them with enough accuracy.

We can use Taylor polynomials to approximate complex functions, but they require many computations and have unpredictable behaviors on the side of the domain.

**Definition.** If $X$ is a Banach space and $M \subseteq X$ is a subset, we say that $p \in M$ is the **<font color=#00ADEF>best approximation</font>** of a function $f \in M$ when
$$
||f - p|| = E(f) := \inf_{q \in M} ||f - q||
\nonumber
$$
**Theorem** (Existence theorem) **.** *If $M$ is a finite-dimensional subspace of $X$, so if $\exists \{v_n\}$ s.t. $M = \text{span}\{v_n\}$, then $\exists p$ the best approximation of $f$ in $M$.*

**Definition.** A vector space $X$ is **<font color=#00ADEF>strictly convex</font>** if $f \neq g$ with $||f||=||g||=1$ and $0 < \theta < 1$, then $||\theta f + (1 - \theta)g|| < 1$. Geometrically, this means that if any $x, y$ on the unit sphere $\partial B$ are joined by a segment that touches $\partial B$ only in $x$ and $y$.

**Theorem** (Uniqueness theorem) **.** *If $X$ is **strictly convex**, then the best approximation $p$ is **unique**.*

**Proof.** We should prove the existence and the in Hilbert spaces uniqueness. Let's prove the uniqueness: let $p_1 \neq p_2$, we have that
$$
\begin{align*}
E(f) &=  ||f - p_1|| = ||f - p_2||\\
&\le \left|\left|f - \frac12(p_1 + p_2) \right|\right| =  \left|\left|\frac12(f - p_1) + \frac12(f - p_2) \right|\right| \\
&< \frac12 ||f - p_1|| + \frac12 ||f - p_2|| =E(f)
\end{align*}
$$
(the first inequality holds because $\frac12 (p_1 + p_2)$ is any function while $p_1$ and $p_2$ are the best approximations) that is impossible, so it must be $p_1 = p_2$.
<img src="http://latex.codecogs.com/gif.latex?\square" border="0" align="right"/>

### Best approximation in Hilbert Spaces

Recall that an Hilbert space is a Banach space plus a scalar product $(\cdot, \cdot)$ with norm defined by $||u||^2:=(u, u)$ and that $p$ is B.A. of $f$ w.r.t. a chosen norm if we have that $||f - p|| \le ||f - q||$ $\forall q \in \mathbb P^n$.

**Theorem** (Best approximation theorem)**.** *Let $H$ be a Hilbert space. Given a function $f \in H$, then $p$ is best approximation of $f$ in $H$ (from $V \subset H$) if and only if*
$$
(f, q) = (p, q) \; \forall q \in V.
\nonumber
$$
**Proof.** ($\Rightarrow$) If $p$ is B.A. of $f$ in $H$, then $||f-p||^2 = E(f)^2 := \inf_{q \in V} ||f-q||^2$ and so we have that for any perturbation $p + tq$ of $p$, it holds
$$
||f - p||^2 \le ||f - p + tq||^2 \quad \forall t>0, \forall q \in V
\nonumber
$$
Consider that $||a+b||^2 - ||a - b||^2 = 4(a, b)$, then we have that, if $a = f- p +\frac t2q$ and $b = \frac t2 q$,
$$
\begin{align*}
0 &\le \left|\left|f - p + tq\right|\right|^2 - \left|\left|f - p \right|\right|^2 = \left|\left|f - p + \frac t2q + \frac t2 q \right|\right|^2 - \left|\left|f - p + \frac t2q - \frac t2q \right|\right|^2 \\
&= 4 \left(f - p + \frac t2q, \frac t2q \right) = 4 \left( f - p, \frac t2q \right) + 4 \left(\frac t2q, \frac t2q \right) = 2t \left(f - p, q \right) + t^2 ||q||^2
\end{align*}
$$
so we have that $(f-p, q) \ge - \frac t2 ||q||^2$.

By adding instead of $tq$ the term $-tq$, with the same reasoning we have that $(f- p, q) \le \frac t2 ||q||^2$, so we have that $\forall t>0, \forall q \in V$
$$
- \frac t2 ||q||^2 \le (f-p, q) \le \frac t2 ||q||^2
\nonumber
$$
which implies that $(f - p, q) = 0$ since a $t$ can be chosen to bound it on both sides. Then we have that
$$
(f - p, q) = 0 \Longleftrightarrow (f, q) - (p, q) = 0 \Longleftrightarrow (f,q) = (p, q)
\nonumber
$$
for each $q \in V$, so we have our thesis.

($\Leftarrow$) If $(f, q) = (p, q) \; \forall q \in V \; \Longleftrightarrow \; (f - p, q) = 0 \; \forall q \in V$, then we have that
$$
||f-q||^2 = ||f-p+p-q||^2 = ||f-p||^2 + ||p - q||^2 + 2(f-p, p-q)
\nonumber
$$
where $2(f-p, p - q) = 0$ since $p - q \in V$ (both $p, q \in V$).

So we have that $||f-q||^2 = ||f-p||^2 + ||p - q||^2 \; \forall q \in V$, which implies that $||f- p||^2 \le ||f-q||^2$ and so $||f-p|| \le ||f - q|| \; \forall q \in V$, and this means that $p$ is B.A. of $f$ in $H$.
<img src="http://latex.codecogs.com/gif.latex?\square" border="0" align="right"/>



<!--In particular, if we consider $H = L^2(\mathbb P^n)$ where $\mathbb P = \text{span}\{v_i\}$ with $i=0, \ldots, n$, we have that, since $(f - p, q) = 0 \; \forall q \Longleftrightarrow (f - p, v_i) \; \forall i=0, \ldots, n$, then $(f, v_i) = (p, v_i) \Longrightarrow (f, v_i) = (\sum_{j=0}^n p_jv_j, v_i)$.-->

<!--Computing integrals is easier than performing interpolation, and it yields better results.-->

Now consider $L^2(0,1)$, that is a Hilbert space where the scalar product between vectors is $(a, b):= \int_0^1 ab \, ds$ for $a, b \in L^2(0, 1)$ and the norm is $||a|| := \sqrt{\int^1_0 |a|^2 ds}$, and take the space $V = \text{span}\{v_i\}_{i=1}^n$.

**The best approximation in $L^2$.** We want the best approximation (in Hilbert Spaces) of the function $f$, on the space $V = \mathrm{span}\{v_i\}$. Then we have seen that $p\in V$ is best approximation of $f$ if and only if:
$$
(f, v) = (p, v), \quad \forall v \in V.
\nonumber
$$
In particular, for every basis functions $v_i \in V$ we have $(p, v_i) = (f, v_i)$, for $i = 1,\ldots,n$. Since $p \in V$, we have that it can be expressed as a linear combination of the basis functions $v_i$ and is uniquely defined by the coefficients $p_j$: $p = \sum_{j=1}^n p_j v_j$, so we have that
$$
(p, v_i) = (\sum_j p_j v_j, v_i) = \sum_j p_j (v_j, v_i).
\nonumber
$$
Collecting this informations together we get:
$$
\sum_{j=1}^n p_j (v_j, v_i) = (f, v_i), \; \forall v_i \in V, i=1,\ldots,n \quad \Longleftrightarrow \quad Mp = F
\nonumber
$$
where $M$ and $F$ are matrices such that $M_{ij} := (v_j, v_i) = \int_0^1 v_jv_i \, ds$ and $F_i := (f, v_i) = \int_0^1 fv_i \, ds$.

If we set $V^n := \text{span}\{x^i\}_{i=0}^{n-1}$ we will then have that $v_i = x^i$, then
$$
M_{ij} := \int_0^1x^jx^i dx = \left. \frac{x^{j+i+1}}{j+i+1} \right|_0^1 = \frac1{j+i+1}.
\nonumber
$$
that is called the $n \times n$ **<font color=#00ADEF>Hilbert matrix</font>** $H$, which is invertible but it is very ill conditioned, so it is difficult to invert, because of collinear lines. Its condition number is
$$
K(H) = ||H|| \,||H^{-1}|| \sim O\left(\left(1+{\sqrt {2}}\right)^{4n}/{\sqrt {n}}\right).
\nonumber
$$
When $n$ increases $K$ explodes, which is very bad. We would like to have $H=I$ the identity, which means $M_{ij} = \delta_{ij}$, so we use the **Legendre basis function**, which are orthonormal basis, to make it orthonormal (perpendicular, a.k.a. diagonal) w.r.t. $L^2$.

We want $v_i \in \mathbb P^n = V^n = \text{span} \{x^i\}$ s.t. $M_{ij} = (v_i, v_j) = \delta_{ij}$. To build it, we use the **Gram-Schmidt process**:
$$
\begin{cases}
v_0 = 1 & f \text{ s.t. } \int_0^1 f = 1 \\
v_{i+1} = \frac{k_{i+1}}{||k_{i+1}||}  & \text{where } k_{i+1} = [x \,v_i - \sum_i (x \, v_i, v_i)v_i]
\end{cases}
\nonumber
$$
or, alternatively

$$
\begin{cases}
\begin{align*}
p_0(x) &= 1 \\
p_k(x) &= x^k - \sum_{j=0}^{k-1} \frac{(x^k,p_j(x))}{(p_j(x),p_j(x))}\\
&= x\,p_{k-1}(x) - \sum_{j=0}^{k-1} \frac{(x p_{k-1}(x),p_j(x))}{(p_j(x),p_j(x))}
\end{align*}
\end{cases}
$$
The set of additive basis having unity as first element. This ensures orthogonality between basis functions.

As the degree $i$ increases, we have that $v_{i+1} = \frac{k_{i+1}}{||k_{i+1}||} \to \infty$ since $x^{i+1} \to \infty$ and thus $k^{i+1} \to 0$. We can avoid instability by using $v_{i+1} = \frac{k_{i+1}}{k_{i+1}(0)}$ instead.

The points created with Gram-Schmidt represent the **Legendre basis**. The make the best approximation $p$ easy to compute, since $M$ becomes easy to invert and we have a diagonal matrix formed by orthogonal basis $p = M^{-1}F$.
