# Finding the Conjectural Bound of the a-function

See https://arxiv.org/abs/2212.10781 for the notation (particular $J$-parameter systems in Definition 3.5).

This code produces the conjectural bound of Lusztig's a-function using the formula in https://arxiv.org/abs/2212.10781 (Conjecture 6.80 combined with Conjecture 6.81). 

The formula is as follows: for the representation $(\pi_{J,\mathsf{v}},M_{J,\mathsf{v}},B_{J,\mathsf{v}})$ (see Section 4 of the paper for the construction of this representation) define

$$\mathbf{a}_{J,\mathsf{v}} = L(\mathsf{w}_0) - \frac{1}{2} \mathrm{deg}{\prod_{\alpha\in \Phi}}'\frac{1-\mathsf{v}^{\alpha^\vee}}{1-\mathsf{q}_{\alpha}^{-2}\mathsf{v}^{\alpha^\vee}}$$

where ${\prod}'$ indicates that we remove factors that are zero from the numerator and the denominator. This gives a conjectural formula for Lusztig's a-function for the elements that are recognised by this representation. 

To find the conjectural bound input conjecturalbound(type,J,sign) where 
- type is the type of reduced, irreducible root system as a string (Eg, "A6", "F4")
- J is a subset of the set of generator indices (eg J = {1,2,3})
- sign is a tuple that represents the chosen $J$-parameter system. Let $a$ be the weight of the short roots and $b$ be the weight of the long roots, that is $\mathsf{q}_\alpha = \mathsf^a$ if $\alpha$ short and $\mathsf{q}_\alpha = \mathsf^b$ if $\alpha$ long. Then $\mathsf{v}_\alpha \in \{\mathsf{q}_\alpha,-\mathsf{q}_\alpha^{-1}\}$. Then sign := $[i,j]$ where $\mathsf{v}_\alpha = i\mathsf{q}^{i}$ for short roots $\alpha$ and $\mathsf{v}_\alpha = j\mathsf{q}^{j}$ for long roots $\alpha$. For example, for type $F4$ if we take $\mathsf{v}_{\alpha_1} = \mathsf{v}_{\alpha_2} = \mathsf{q}^{a}$ and $\mathsf{v}_{\alpha_3} = \mathsf{v}_{\alpha_4} = -\mathsf{q}^{-b}$ then sign = $[1,-1]$. For the simply laced case only the longroot input will matter but sign should still be inputed as a tuple (input [-1,-1] to get the bound).


For a classification of which of these representations are bounded (and thus when this formula makes sense) see Proposition 5.11 and Theorem 5.9 of the paper. 

The output of the conjecturalbound function is a list of ranges and then their corresponding bounds. The ranges are of the value of $a/b$. The ranges go up to 100 however not all of these bounds make sense (check with the boundedness classification from the paper). 

For example, let type = "F4", $J = {1,2,3}$ and sign = $[1,-1]$. Then the output is 

[
    [ 0, 2/3 ],
    [ 2/3, 1 ],
    [ 1, 6/5 ],
    [ 6/5, 2 ],
    [ 2, 3 ],
    [ 3, 4 ],
    [ 4, 6 ],
    [ 6, 100 ]
]
[
    2*b + 2*a,
    5*a,
    -b + 6*a,
    -7*b + 11*a,
    b + 7*a,
    4*b + 6*a,
    12*b + 4*a,
    12*b + 4*a
]

By Proposition 5.11, the representation is bounded when $a/b\leq 2$ therefore we have the following conjectural bounds:

$$\mathbf{a}_{J,\mathsf{v}} = \begin{cases} 2a+2b &\text{ if } 0\leq \frac{a}{b}\leq \frac{2}{3}\\
5a &\text{ if } \frac{2}{3}\leq \frac{a}{b}\leq 1 \\
6a-b &\text{ if } 1\leq \frac{a}{b}\leq \frac{6}{5}\\
11a-7b &\text{ if } \frac{6}{5}\leq \frac{a}{b}\leq 2 \end{cases} $$

