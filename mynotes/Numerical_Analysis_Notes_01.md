# Numerical Analysis notes



## Foundations of matrix analysis

Numerical Analysis is the *“art of approximating”*. Quoting Wikipedia:

> An approximation is an inexact representation of something that is still close enough to be useful. Although approximation is most often applied to numbers, it is also frequently applied to such things as mathematical functions, shapes, and physical laws.

Approximations should be used when incomplete information prevents use of exact representations. Many problems in physics are either too complex to solve analytically, or impossible to solve using the available analytical tools. Thus, even when the exact representation is known, an approximation may yield a sufficiently accurate solution while reducing the complexity of the problem significantly.

The type of approximation used depends on the available information, the degree of accuracy required, the sensitivity of the problem to this data, and the savings (usually in time and effort) that can be achieved by approximation.

In this course we focus on three main aspects:

- Methodologies (or approximation algorithms)
- Analysis (estimate errors and convergence properties)
- Implementation (through $\texttt{Python}$ and `numpy` notebooks)

Approximation is a matter of 

- Representation (floating point values VS real numbers, finite dimensional spaces VS infinite dimensional ones, etc.)
- Measure of the Error (how do we know that we did a good job in approximating?)

In general, we will end up working with $\mathbb R^n$.

<!--[...pages 1-3 Gabri notes...]-->

### Matrices

Denoting by $A_{ij}$ the matrix of order $n - 1$ obtained from $A$ by eliminating the $i$-th row and the $j$-th column, we call the **<font color=#00ADEF>complementary minor</font>** associated with the entry $a_{ij}$ the determinant of the matrix $A_{ij}$. We call the **<font color=#00ADEF>$k$-th principal (dominating) minor</font>** of $A$, $d_k$, the determinant of the principal submatrix of order $k$, $A_k = A(1 : k, 1 : k)$.

----

The computational cost is O (ops) and can be constant, linear, polynomial, exponential, factorial, etc.

*Total error*: $f(x) - \hat f(\hat x) = \underbrace{\hat f(\hat x) - f(\hat x)}_{\text{computational error }e_c = e_t + e_r} + \underbrace{f(\hat x) - f(x)}_{\text{propagated data error (independent from $f$)}}$

**Example:** finite difference approximation $f'(x) = \lim_{h \to 0} \frac{f(x - h) - f(x)}{h}$

- *truncation error* (obtained through Taylor) $\sim \frac12 |f''(x)|h + O(h^2)$
- *rounding error* $\sim \frac{2\varepsilon}h$, with $\varepsilon =$ machine precision

The optimal $h$ is thus $h = 2 \sqrt{\frac{\varepsilon}{|f''(x)|}}$.

### Norms of vectors, matrices and functions

Given a vector space $V$ over the field of real ${\mathbb{R}^{}}$ or complex numbers ${\mathbb{C}^{}}$ ($V$ might be infinite dimensional), a **<font color=#00ADEF>semi-norm</font>** on $V$ is a function ${|\cdot|}: V\rightarrow
{\mathbb{C}^{}}$ satisfying:

1.  ${|c f|}={|c|}{|f|}$, for all $c\in {\mathbb{C}^{}}$   [homogeneity];

2.  ${|f+g|}\le {|f|}+{|g|}$   [triangle inequality].

As it can be easily seen (1)–(2) imply that the norm is always non-negative:
$$
0 = 0 \cdot{|f|} = {|0\cdot f|} = {|(1-1)f|} = {|f-f|} \le {|f|}+{|(-1)f|} = 2{\|f|}.
\nonumber
$$
The semi-norm becomes a **<font color=#00ADEF>norm</font>** if in addition to (1)–(2) we have also that for all $f\in V$
$$
{|f|}=0 \text{ if and only if } f=0.
\nonumber
$$
A vector space is said to be **<font color=#00ADEF>complete</font>** if every Cauchy sequence in that space converges to one of the space's elements.

A complete vector space with a norm is called a **<font color=#00ADEF>Banach space</font>**.

A scalar product is a function $(\cdot,\cdot):V\times V\mapsto {\mathbb{C}^{}}$ which is:

1. **symmetric**: $(f,g)=\overline{(g,f)}$;

2. **linear**:

   - $(\alpha f,g)=\alpha (f,g)$ for all $\alpha\in {\mathbb{C}^{}}$;
   - $(f+g,h)=(f,h)+(g,h)$;

   equivalent to

   - $(\alpha_1 f + \alpha_2 g,h) = \alpha_1 (f,h) + \alpha_2 (g,h)$

3. **positive definite**: $(f,f)\ge 0$ and $(f,f)=0$ if and only if $f=0$.

The norm is then defined as $\|f\|^2=(f,f)$. That this is a norm (i.e. satisfies the triangle inequality) is proved by using the **<font color=#00ADEF>Cauchy-Schwarz</font>** inequality
$$
|(f,g)| \le \sqrt{(f,f)(g,g)}.
\nonumber
$$
A Banach space with a scalar product and a norm induced by it is called a **<font color=#00ADEF>Hilbert space</font>**.

**Example:** Some norms in Banach spaces
$$
\begin{aligned}[c]
\|x\|_p &= \left(\sum_{i=1}^n {|x_i|}^p \right)^{\frac 1p}, \quad 1\le p  < \infty, &\|x\|_1 = \sum_{i=1}^n |x_i| \\
\|x\|_2 &= \sqrt{x_1^2 + x_2^2} \quad \text{[Euclidean norm]} & \|x\|_\infty = \sup_{1\le i\le n} {|x_i|}.
\end{aligned}
$$
**Definition:** Two norms $||\cdot||_p, ||\cdot||_q$ on a vector space $V$ are **<font color=#00ADEF>equivalent</font>** if there exist two positive constants $c_1$ and $c_2$ such that:
$$
c_1 ||x||_q \le ||x||_p \le c_2 ||x||_q \quad \forall x \in V
\nonumber
$$
In a finite-dimensional normed vector space (the dimension is given by the number of vectors in the basis) all norms are equivalent.