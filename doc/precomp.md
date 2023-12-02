# Precomputation & multiple blinds per group element

To further save client-side commitment time **beyond the paper**, it is possible to store precomputation tables of group elements w<sub>1</sub>, ..., w<sub>d</sub>. The problem with optimization is space constraint: the precomputation table of a single element takes 120KB. To reach a balance between time and space, a solution is to use multiple blinds per group element. 

Here are the details. Client i: 
- generates secrets r<sub>i</sub> = (r<sub>i1</sub>, ..., r<sub>ie</sub>), e is **the number of blinds per group element**
- commits u<sub>i</sub> = (u<sub>i1</sub>, ..., u<sub>1d</sub>) with y<sub>i</sub> = (y<sub>i1</sub>, ..., y<sub>1d</sub>), where y<sub>ij</sub> = g<sup>u<sub>ij</sub></sup> w<sub>p</sub><sup>r<sub>iq</sub></sup>, where j - 1 = (p - 1) e + (q - 1), 1 ≤ q ≤ e
- shares r<sub>i</sub> with Shamir

We only need ceiling(d / e) many group elements to make commitments at the cost of sharing e elements via Shamir. 
The paper becomes a special case with e = 1. The proof generation and verification parts are modified accordingly. For instance, when d = 1M, choosing e = 100 makes it possible to store precomputation tables in memory, while the extra cost of Shamir sharing and proof generation/verification is small. Experiments show that the overall speedup is small (<1.5x), so this optimization is not mentioned in the paper. 