## Doubly Fair Dynamic Pricing implementation

This group project is an attempt to implement a reaseach paper named : " **Doubly Fair Dynamic Pricing** " in python published.

In the `**code**` folder, you will find the actual code implementation of the **algorithms** crafted in the paper.
There is also a `**presentataion**` in the codebase to understand the **crux** of the research paper and the algorithms in it.

### Abstract

We study the problem of online dynamic pricing with two types of fairness constraints: a "procedural fairness" which requires the proposed prices to be equal in expectation among different groups, and a "substantive fairness" which requires the accepted prices to be equal in expectation among different groups. A policy that is simultaneously procedural and substantive fair is referred to as "doubly fair". We show that a doubly fair policy must be random to have higher revenue than the best trivial policy that assigns the same price to different groups. In a two-group setting, we propose an online learning algorithm for the 2-group pricing problems that achieves 
{{<raw>}}
\( \tilde{O}(\sqrt{T}) \)
{{</raw>}}
regret, zero procedural unfairness and
{{<raw>}}
\( \tilde{O}(\sqrt{T}) \)
{{</raw>}}
substantive unfairness over rounds of learning. We also prove two lower bounds showing that these results on regret and unfairness are both information-theoretically optimal up to iterated logarithmic factors. To the best of our knowledge, this is the first dynamic pricing algorithm that learns to price while satisfying two fairness constraints at the same time.
