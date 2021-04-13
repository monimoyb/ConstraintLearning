# Learning to Satisfy Unknown Constraints in Iterative MPC

These set of codes replicate the results of the paper https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=9303765. Two constraint learning strategies are presented for iterative tasks solving robust MPC problems with a linear system. The first strategy constructs an inner approximation of the true (unknown) constraints using CVX hull of the recorded states that satisfy the true constraints. The second strategy is based on a kernelized SVM classifier, which constructs constraint estimates based on a user-specified tolerable probability of constraint violations. A set of example estimates are shown below. 

<img width="1171" alt="con" src="https://user-images.githubusercontent.com/12418616/114485677-2fed3400-9bc1-11eb-8239-f9bcecab5e6f.png">

After using the SVM constructed constraint estimates for control design, the actual failure (i.e., violation of the true constraints) is significantly lower than the user-specified acceptable limits, as shown in the table below. This trend is obtained from 100 Monte Carlo simulations. A safety vs performance trade-off is also seen in this data. The more the allowable probability of constraint violations, the lower is the average cost of successful iterations. Thus, the CVX hull estimate is the safest, but least cost effective.

![tab](https://user-images.githubusercontent.com/12418616/114486202-3def8480-9bc2-11eb-994c-e0e50f3b356c.png)
