# Finding the Conjectural Bound of the a-function

See https://arxiv.org/abs/2212.10781 for the notation (particular $J$-parameter systems and in Definition 3.5).

This code produces the conjectural bound of Lusztig's a-function using the formula in https://arxiv.org/abs/2212.10781 (Conjecture 6.80 combined with Conjecture 6.81). 

The formula is as follows: for representation $(\pi_{J,\mathsf{v}},M_{J,\mathsf{v}},B_{J,\mathsf{v}})$ (see Section 4 of the paper for the construction of this representation) define

$$ a_{J,\mathsf{v}} = L(\mathsf{w}_0) - \frac{1}{2} \mathrm{deg}{\prod_{\alpha\in \Phi}}'\frac{1-\mathsf{v}^{\alpha^\vee}}{1-\mathsf{q}_{\alpha}^{-2}\mathsf{v}^{\alpha^\vee}}$$

where ${\prod}'$ indicates that we remove factors that are zero from the numerator and the denominator. This gives a conjectural formula for Lusztig's a-function for the elements that are recognised by this representation. 


