Absolutely. Let us start with the **abstract**, but unpack it as if we are rebuilding the paper from first principles.

## 1) What problem is the paper actually about?

The paper studies **online dynamic pricing** with **two groups of customers** and **two fairness constraints**. In plain language:

* You repeatedly offer prices over time.
* You do **not** know in advance how much customers are willing to pay.
* You observe only a binary response: buy or not buy.
* You want to maximize revenue.
* At the same time, you want the pricing process to be fair across groups.  

This is not just “set one good price.” It is a **sequential learning problem**: the seller must learn demand from feedback while also respecting fairness.

---

## 2) The classical dynamic pricing setup behind the abstract

The standard online pricing model is:

1. A customer arrives.
2. The seller proposes a price (v_t) without knowing the customer’s private valuation (y_t).
3. The customer buys if (v_t \le y_t).
4. The seller earns revenue (r_t = v_t \cdot \mathbf{1}(v_t \le y_t)). 

The important point is that the seller sees only **censored feedback**: buy or no-buy. The seller does **not** observe (y_t). So the true demand curve must be learned indirectly from repeated interaction.

That is why this problem belongs to **online learning / bandits / sequential decision-making**, not just static optimization.

---

## 3) What is new in this paper: fairness in pricing

The paper argues that pricing can become unfair when different groups are treated differently, especially when those groups differ by sensitive attributes such as gender, race, or age. The authors motivate two distinct notions of fairness:

* **Procedural fairness**: the **proposed prices** should be equal in expectation across groups.
* **Substantive fairness**: the **accepted prices** should be equal in expectation across groups.

The distinction matters.

### Procedural fairness

This is about the **price offered** by the seller.

If Group 1 is offered a systematically higher price than Group 2, the process itself is unfair, even if some buyers refuse and the final transaction outcomes look similar.

Mathematically, the paper defines procedural unfairness as the absolute difference between the expected proposed prices of the two groups. 

### Substantive fairness

This is about the **price actually paid by buyers who accept**.

Even if you offer both groups similar average prices, the buyers who end up buying might still face different average accepted prices. This is a deeper, outcome-based fairness notion. The paper defines substantive unfairness using the difference in expected accepted prices. 

So the paper is not merely asking for equal treatment at the offer stage. It wants equality in the realized outcome as well.

---

## 4) Why the word “doubly fair”?

“Doubly fair” means the policy is fair in **both** senses at once:

* the average proposed price is the same across groups,
* the average accepted price is also the same across groups.

This is a strong requirement.

A useful way to think about it:

* Procedural fairness controls **what the seller does**.
* Substantive fairness controls **what the market outcome becomes**.

The paper’s central claim is that satisfying both simultaneously is hard, and in general the optimal fair policy must be **randomized**.

---

## 5) Why randomization is necessary

This is one of the most important ideas in the abstract.

The paper says that a doubly fair policy must be random if it wants to outperform the trivial policy that offers the same fixed price to every group.

### Why is that?

If you force deterministic prices and also require fairness across groups, then the simplest way to satisfy both fairness conditions is often to offer the **same fixed price** to everyone. But that can be revenue-suboptimal.

Randomization lets you mix prices in a way that preserves fairness **in expectation** while still extracting more revenue.

This is a subtle but important mathematical idea:

* A single deterministic price may lie in the intersection of fairness constraints, but that intersection may be too restrictive.
* A **distribution over prices** can satisfy fairness constraints on average while achieving a higher expected revenue.

This is exactly what the example in the paper demonstrates.

---

## 6) The paper’s toy example, unpacked

The paper gives a concrete example with:

* two groups,
* three candidate prices: (0.625), (0.7), and (1),
* different acceptance rates across groups.

The key result in that example is:

* If you insist on deterministic pricing and both fairness constraints, the only feasible behavior is essentially to propose the same price to both groups.
* The best such fixed price is (1), yielding expected revenue (0.5).

But the paper then constructs a **random policy**:

* Group 1 gets (0.625) with probability (20/29), and (1) with probability (9/29).
* Group 2 gets (0.7) with probability (25/29), and (1) with probability (4/29).

Under this randomized policy:

* the expected proposed price is the same across groups,
* the expected accepted price is the same across groups,
* and the expected revenue is (74/145), which is greater than (0.5).

### Why this matters mathematically

This example proves a conceptual point:

> Fairness does not force you into the lowest-revenue trivial solution, provided you allow randomization and measure fairness in expectation.

That is one of the paper’s biggest contributions.

---

## 7) What does “learn online” mean here?

The seller does **not know** the demand distributions (D_1) and (D_2) in advance. That is why the problem becomes a learning problem rather than a simple optimization problem. The seller must estimate acceptance behavior over time from binary feedback.

This is crucial:

* If you knew the full demand curves, you could try to solve the fair optimization problem directly.
* But in the paper, the demand curves are unknown.
* Therefore the algorithm must balance:

  * **exploration**: try prices to learn demand,
  * **exploitation**: use the learned information to earn revenue,
  * **fairness constraints**: never drift too far from the fair objective.

So the abstract is really describing a **constrained online learning** problem.

---

## 8) What is the goal of the algorithm?

The paper’s target is an optimal pricing policy (\pi^*) under the fairness constraints. The policy itself is defined as a pair of distributions, one for each group, and it maximizes expected revenue subject to both fairness conditions. 

The abstract says the algorithm achieves:

* ( \tilde O(\sqrt{T}) ) regret,
* zero procedural unfairness,
* ( \tilde O(\sqrt{T}) ) substantive unfairness.

Let us unpack each of these.

---

## 9) What is regret here?

Regret is the performance gap between the algorithm’s revenue and the revenue of the best fair policy. The paper defines cumulative regret as the sum over time of the difference between the optimal fair revenue and the revenue achieved by the policy used at each round.

So if regret is small, the algorithm is learning a policy nearly as good as the best possible fair policy.

### Why ( \sqrt{T} ) is important

In online learning, ( \sqrt{T} )-type regret is the classical sign of a good no-regret algorithm. It means average regret per round goes to zero:

[
\frac{\sqrt{T}}{T} = \frac{1}{\sqrt{T}} \to 0.
]

So over a long horizon, the algorithm becomes asymptotically optimal.

The paper emphasizes that its regret rate is optimal up to iterated logarithmic factors.

---

## 10) What does “zero procedural unfairness” mean operationally?

This means the algorithm is designed so that the **expected proposed prices** are exactly equal across the groups. The paper can guarantee this exactly, not just approximately.

This is easier than substantive fairness because procedural fairness depends directly on the policy the seller chooses. It is a property of the decision rule itself, not of the hidden demand environment.

---

## 11) Why is substantive unfairness harder?

Substantive fairness depends on **accepted prices**, and accepted prices depend on the unknown customer valuation distributions (D_1) and (D_2). That makes it intrinsically harder because:

* you do not observe valuations directly,
* you only observe acceptance or rejection,
* the fairness metric contains a ratio of expectations,
* so the fairness constraint is non-linear and non-convex.

This is one of the core technical difficulties in the paper.

That is why the paper does not promise perfect substantive fairness at every round. Instead, it proves an ( \tilde O(\sqrt{T}) ) cumulative substantive unfairness guarantee.

---

## 12) What does “information-theoretically optimal” mean?

The paper says the regret and unfairness rates are optimal up to iterated logarithmic factors.

This means:

* no algorithm can fundamentally beat the stated scaling in (T),
* the paper proves matching lower bounds,
* so the algorithm is not just good empirically; it is near-best possible in the worst case.

The abstract is signaling a strong theoretical result, not just a heuristic one.

---

## 13) Why are there lower bounds?

The paper proves two kinds of lower bounds:

1. a regret lower bound,
2. a substantive unfairness lower bound for algorithms that are regret-optimal.

This is important because it shows the fairness constraint creates a genuine trade-off.

In particular, the paper argues that if you want near-optimal regret, you may still be forced to accept some substantive unfairness during learning. That is not a flaw in the algorithm; it is a structural property of the problem.

---

## 14) A first-principles intuition for the abstract

Here is the whole abstract in one clean mental model:

You are selling a product over time to two groups. You do not know how much each group values the product. You want to learn the best prices from repeated feedback. But you also want the process to be fair in two ways:

* the groups should not be offered systematically different prices,
* and the groups should not end up paying systematically different accepted prices.

If you insist on fairness deterministically, you may lose too much revenue. If you allow randomized pricing, you can do better.

So the paper builds an online learning algorithm that learns these pricing distributions over time, with provable guarantees that:

* revenue is nearly optimal,
* procedural fairness is exact,
* substantive fairness becomes small over time,
* and all of this is essentially best possible.

---

## 15) What the abstract is *really* claiming, between the lines

Reading carefully, the abstract is making four deep claims:

### (a) Fair pricing is not just “equal prices”

The paper is not about naive equality. It distinguishes between fairness in the offer stage and fairness in the outcome stage.

### (b) Randomization is not an implementation trick; it is mathematically necessary

The optimal fair policy can lie outside the set of deterministic prices. Random mixtures are essential.

### (c) Fairness and learning interact in a nontrivial way

You are not just optimizing a known constrained problem. You are learning the constraints themselves from feedback.

### (d) The paper is about a near-optimal theoretical frontier

The regret and fairness bounds are matched by lower bounds, so the algorithm is close to best possible.

---

## 16) What to understand before moving to the next section

Before we move past the abstract, you should have these objects very clear:

* a **policy** is a distribution over prices for each group,
* **revenue** is the expected accepted price times the price,
* **procedural fairness** is equality of proposed prices in expectation,
* **substantive fairness** is equality of accepted prices in expectation,
* **regret** is the revenue gap to the best fair policy,
* the optimal fair policy is generally **randomized**.

Once these are solid, the rest of the paper becomes much easier.

The next natural step is to go through the **Introduction and Problem Setup**, because that is where the paper turns these ideas into precise mathematical definitions.

Great. The **Introduction** is where the paper turns the abstract into a full research problem. It does four things in order:

1. it recalls the classic online pricing model,
2. it explains why fairness becomes unavoidable in pricing,
3. it shows why deterministic fairness can be too restrictive, and
4. it previews the key technical idea: **randomized pricing with online learning**.

I will unpack each part carefully.

## 1) The classical online pricing model

The paper begins with the standard dynamic pricing setup. At each round (t):

* a customer arrives with a private valuation (y_t),
* the seller chooses a price (v_t) without seeing (y_t),
* the customer buys if (v_t \le y_t),
* the seller earns revenue (r_t = v_t \cdot \mathbf{1}(v_t \le y_t)). 

This is a **censored-feedback** problem. The seller does not observe the valuation itself, only whether the customer accepted the price. That means the seller must learn the demand distribution indirectly from yes/no responses. This is why the problem belongs to **online learning** and not just ordinary pricing optimization.

### First-principles interpretation

If you think of valuation as “maximum willingness to pay,” then pricing is a repeated experiment:

* too high a price gives rejections and little learning,
* too low a price gives acceptance but less revenue,
* intermediate prices help reveal the shape of demand.

So the seller is solving a two-way tradeoff:

* **explore** to learn demand,
* **exploit** to earn revenue.

That is the classic dynamic pricing problem the paper starts from.

---

## 2) Why fairness enters the pricing problem

The paper then says: dynamic pricing becomes socially sensitive once customers are split into groups such as gender, race, age, or any other observable segment. If the seller can price differently across groups, the system can drift into price discrimination, which raises fairness concerns and reputation risk.

This is not a superficial concern. In pricing, different groups may have different learned behavior or different willingness to pay. If the algorithm simply maximizes revenue without fairness constraints, it may end up charging one group more than another in a systematic way. That is where the paper’s fairness definitions begin to matter.

### Why pricing fairness is tricky

A pricing algorithm is not like a classifier where you output a label once. Here the system **interacts over time**, and the choices it makes today affect:

* what it learns tomorrow,
* how customers perceive the seller,
* and which outcomes become statistically visible.

So fairness has to be defined in a way that is compatible with repeated decision-making and learning from feedback.

---

## 3) The two fairness notions introduced in the introduction

The paper introduces fairness using the language of **procedural** and **substantive** unconscionability, adapted into pricing fairness.

### Procedural fairness

This is about the **prices offered**.

If two groups are offered systematically different average prices, the process itself is unfair. The paper later formalizes this as the difference in expected proposed prices across groups.

### Substantive fairness

This is about the **prices actually accepted**.

Even if the seller offers similar prices, the set of customers who accept in each group may differ, and the average accepted price can still differ. The paper defines this as the difference in expected accepted prices across groups.

### Why the paper separates them

This separation is important because it distinguishes between:

* fairness of the **decision rule**,
* fairness of the **realized outcome**.

That is a much stronger and more nuanced statement than just “do not discriminate.” It says the algorithm must be fair both in what it proposes and in what the market ends up paying on average.

---

## 4) The paper’s core warning: fairness and revenue can conflict

The introduction makes a key point: if you try to satisfy both fairness constraints with deterministic pricing, you may be forced into a trivial solution where both groups receive the same price. The paper explicitly says that, under deterministic pricing, the only way to guarantee both fairness constraints in general is to set the same price across groups.

That sounds fair, but it can be very expensive in revenue.

### Intuition

Suppose:

* Group 1 has one demand curve,
* Group 2 has another demand curve,
* the revenue-maximizing price differs across groups.

If you force the same deterministic price everywhere, you may lose revenue because one group’s optimal price is not the other group’s optimal price. So fairness can collapse the feasible set too much.

This is the reason the paper does **not** stop at deterministic fairness.

---

## 5) Example 1: why randomization matters

The introduction’s example is the paper’s key intuition engine. It shows two groups, three candidate prices, and group-specific acceptance rates:

* Group 1: acceptance rates (3/5, 1/2, 1/2)
* Group 2: acceptance rates (4/5, 4/5, 1/2)

for prices (0.625, 0.7, 1).

The paper shows that if you require both fairness constraints deterministically, the best feasible price is the same for both groups, and the best such fixed price is (1), giving expected revenue (0.5). 

But then the paper does something deeper: it constructs a **randomized policy**.

* For Group 1: offer (0.625) with probability (20/29), and (1) with probability (9/29).
* For Group 2: offer (0.7) with probability (25/29), and (1) with probability (4/29). 

Under this policy:

* the expected proposed price is equal across groups,
* the expected accepted price is equal across groups,
* and the expected revenue becomes (74/145), which is greater than (0.5). 

### What this example is proving

It proves that:

1. fairness does **not** force the lowest-revenue trivial strategy,
2. but you must allow **randomization**,
3. and fairness must be interpreted **in expectation**.

This is one of the paper’s central conceptual contributions.

---

## 6) Why random pricing is mathematically natural here

The introduction is not saying randomness is a hack. It is saying randomness is the correct mathematical object once fairness is defined in expectation.

A pricing policy is no longer “choose one price.” It becomes:

* a **distribution over prices** for Group 1,
* a **distribution over prices** for Group 2. 

That matters because expected values behave linearly, so the seller can tune the distribution to satisfy constraints that a single price cannot satisfy.

### Intuition from convex geometry

Deterministic prices are points. Randomized policies are mixtures of points, so they live in a convexified action space. Many fairness constraints become easier to satisfy in the convex hull than at the vertices. That is the deeper reason randomization unlocks better fair outcomes. The paper’s example is a concrete demonstration of this phenomenon.

---

## 7) The paper’s formal viewpoint introduced in the introduction

The introduction then states the policy more formally:

* a policy is a pair (\pi = (\pi_1, \pi_2)),
* where (\pi_1) and (\pi_2) are distributions over prices for Group 1 and Group 2,
* the objective is to maximize weighted expected revenue,
* subject to equality of expected proposed prices and equality of expected accepted prices. 

This is important because it shifts the problem from “choose a price” to “choose a randomized rule.”

### Why this is harder than ordinary pricing

Even if the true demand distributions (D_1) and (D_2) were known, the fairness-constrained optimization is not trivial, because the substantive fairness constraint involves conditional expectations of the form

[
\mathbb{E}[v \mid \text{accepted}],
]

which is a ratio of two expectations and is therefore nonlinear and generally non-convex. The paper explicitly points this out. 

So the introduction is already warning you: this is not a simple linear program over known demand curves. It is a hard constrained optimization problem even before learning enters.

---

## 8) Why learning makes it harder

The seller does not know (D_1) and (D_2) in advance. The only data comes from acceptance/rejection feedback. So the algorithm has to **learn** the demand structure while also maintaining fairness.

That creates three simultaneous difficulties:

* the seller must estimate demand,
* the seller must optimize revenue,
* the seller must keep fairness under control.

And because the fairness constraint depends on unknown acceptance probabilities, it cannot be enforced perfectly at every single time step without additional assumptions. That is why the paper moves to an online-learning framing with cumulative guarantees.

---

## 9) The introduction’s preview of the algorithmic idea

The introduction already hints at the structure of the algorithm:

* use an initial phase to estimate a lower bound on acceptance probabilities,
* divide time into epochs,
* in each epoch, select good-and-exploratory policies,
* estimate acceptance probabilities,
* solve an empirical fair optimization problem,
* eliminate policies that are too unfair or too suboptimal.

This is the core design principle of FPA, the Fairly Pricing Algorithm.

### Why epochs?

Epochs are used because the algorithm needs periods long enough to gather stable estimates. If policies changed every round, the seller would keep resetting the learning process and the fairness signals would be harder to estimate. The paper even comments later that customers should “feel” fairness through sustained policy windows rather than constant switching.

---

## 10) How the introduction positions the contribution

The introduction also frames the contribution relative to earlier pricing-fairness work:

* some prior work models fairness as a soft penalty or utility tradeoff,
* some uses hard fairness constraints,
* some fairness notions are impossible to satisfy simultaneously with deterministic prices,
* randomization is one way around that. 

The paper’s claim is more specific:

* it studies **two-group online pricing**,
* with **two fairness constraints at once**,
* in a **learning** setting,
* and proves near-optimal regret and unfairness guarantees.

That is the novelty being prepared in the introduction.

---

## 11) The assumptions that appear early in the paper

The introduction also foreshadows a few technical assumptions that will matter later:

* the number of candidate prices (d) is finite,
* the highest price has a positive lower bound on acceptance probability,
* the fair optimization is carried out over a fixed price set.

These assumptions are not just formalities. They are needed so the fairness ratio is well-defined and the learning problem is tractable.

---

## 12) What you should take away from the introduction

The introduction is telling you the following story:

### First

Dynamic pricing is a learning problem under partial feedback.

### Second

When customers belong to different groups, pricing can become unfair in two distinct ways: in what is offered and in what is ultimately paid.

### Third

If you require both fairness constraints deterministically, you may lose too much revenue.

### Fourth

Randomized pricing policies can satisfy both fairness constraints in expectation and still outperform the best fixed-price fair policy.

### Fifth

The resulting learning problem is hard because fairness depends on unknown acceptance probabilities, so the algorithm must learn and optimize simultaneously.

### Sixth

The paper’s main contribution is an epoch-based policy-elimination algorithm with near-optimal theoretical guarantees.

---

## 13) Mental model to carry forward

A very good way to think about the introduction is this:

* **Classic dynamic pricing** asks: “What price should I set?”
* **This paper** asks: “What randomized pricing rule should I learn so that revenue is high and fairness is preserved in two different senses?”

That is the conceptual leap.

Once this is clear, the later sections become much easier:

* **Problem Setup** will define the exact math,
* **Algorithm** will build the epoch machinery,
* **Analysis** will prove regret and unfairness guarantees,
* **Discussion** will explain limitations and extensions.

If you want, I will continue next with the **Problem Setup section**, and derive every definition carefully from first principles so the notation becomes completely transparent.


## Problem Setup: what the paper is actually formalizing

This section turns the story into a mathematical decision problem. The key move is simple:

* every customer belongs to one of two groups,
* the seller chooses prices from a finite set,
* the seller does not know the customer’s valuation,
* the seller only observes buy / no-buy,
* and the seller wants both revenue and fairness. 

So the whole paper is building a framework for **learning a randomized pricing rule** under two fairness constraints.

---

## 1) The online pricing process, step by step

At round (t):

1. A customer arrives.
2. Their group (G_e) is observed, where (e \in {1,2}).
3. The seller picks a price from a finite set
   [
   V={v_1,v_2,\dots,v_d}, \quad 0 < v_1 < \cdots < v_d \le 1.
   ]
4. The customer has a private valuation (y_t^e), which the seller does not see.
5. The customer buys if the offered price is at most their valuation:
   [
   \mathbf{1}(v_t^e \le y_t^e).
   ]
6. If they buy, the seller earns revenue equal to the price; otherwise revenue is (0). 

### Why this is important

The seller never sees the valuation directly. They only see whether the customer accepted the price. So the seller is learning from **censored feedback**: accept/reject, not the underlying willingness to pay.

That is what makes this a dynamic learning problem rather than a static pricing problem.

---

## 2) What a “policy” means here

The paper does **not** mean “one price.” It means a **distribution of prices**.

For each group (e), the seller uses a probability vector
[
\pi^e \in \Delta^d,
]
where (\Delta^d) is the probability simplex: all length-(d) vectors with nonnegative entries summing to 1. So (\pi^e(i)) is the probability of offering price (v_i) to group (e). The full policy is

[
\pi=(\pi^1,\pi^2).
]

So a policy is: “for Group 1, use this price distribution; for Group 2, use that price distribution.” 

### Why distributions instead of fixed prices?

Because the paper later shows that deterministic prices are often too restrictive under fairness. A mixture over prices can satisfy fairness while earning more revenue.

---

## 3) Acceptance probabilities: what the unknown demand looks like

For each group (e) and each price (v_i), define

[
F_e(i)=\Pr(y_t^e \ge v_i).
]

This is the probability that a customer from group (e) accepts price (v_i). Since prices are ordered from low to high, these acceptance probabilities are nonincreasing in (i): higher prices are accepted less often. 

The paper packages these into a diagonal matrix:

[
F_e = \operatorname{diag}(F_e(1),F_e(2),\dots,F_e(d)).
]

That is just a compact way to store the acceptance probabilities.

### Intuition

Think of (F_e(i)) as the demand curve, but only sampled at the finite candidate prices.

---

## 4) Expected proposed price: the average price you offer

If group (e) receives prices according to (\pi^e), then the expected proposed price is

[
\mathbb{E}[v_t^e] = \sum_{i=1}^d v_i \pi^e(i) = v^\top \pi^e.
]

This is just the weighted average of prices under the chosen price distribution. 

### Why this matters

This quantity is what procedural fairness will compare across groups. It measures the **offer stage**, not the actual transaction stage.

---

## 5) Expected reward: how revenue is computed

Now the crucial derivation.

If the seller offers price (v_i), then the customer buys with probability (F_e(i)). The expected revenue from that price is therefore

[
v_i \cdot F_e(i).
]

If the seller randomizes over prices using (\pi^e), then by the law of total expectation:

[
\mathbb{E}[r_t^e]
= \sum_{i=1}^d \pi^e(i), v_i, F_e(i).
]

Using matrix notation, this becomes

[
\mathbb{E}[r_t^e] = v^\top F_e \pi^e.
]



### Step-by-step intuition

* pick a price (v_i) with probability (\pi^e(i)),
* that price is accepted with probability (F_e(i)),
* revenue is (v_i) if accepted, else (0),
* so expected revenue is “price × acceptance probability,” averaged over the price distribution.

That is the basic revenue formula used throughout the paper.

---

## 6) Expected acceptance rate

The acceptance rate is the probability that the customer buys, regardless of what price was offered. So for group (e),

[
\Pr(\text{accept} \mid e)
= \sum_{i=1}^d \pi^e(i)F_e(i)
= \mathbf{1}^\top F_e \pi^e.
]

This is the denominator that later appears in accepted-price calculations. 

### Intuition

This is just: “if I sample a price from my policy, how often do buyers accept it?”

---

## 7) Expected accepted price: the average price among buyers

This is the subtle one, and it is where the paper’s substantive fairness comes from.

The expected accepted price for group (e) is:

[
\mathbb{E}[v \mid \text{accepted}, e]
=====================================

\frac{\mathbb{E}[v \cdot \mathbf{1}(\text{accepted}) \mid e]}
{\Pr(\text{accepted} \mid e)}.
]

The numerator is the same as expected revenue:

[
\mathbb{E}[v \cdot \mathbf{1}(\text{accepted}) \mid e]
= v^\top F_e \pi^e.
]

The denominator is the acceptance rate:

[
\Pr(\text{accepted} \mid e)=\mathbf{1}^\top F_e \pi^e.
]

So the expected accepted price is

[
\frac{v^\top F_e \pi^e}{\mathbf{1}^\top F_e \pi^e}.
]



### Why this is important

This quantity answers: among the customers who actually buy, what average price did they pay?

That is not the same as the average offered price. A group could be offered lower prices on average but still end up paying more among those who accept, depending on the demand profile.

This is exactly why the paper introduces a second fairness notion.

---

## 8) Overall expected revenue across both groups

Let (q) be the probability that an arriving customer belongs to Group 1, so Group 2 has probability (1-q). Then the total expected revenue of a policy (\pi) is the weighted average of the group-specific expected revenues:

[
R(\pi;F_1,F_2)
==============

q,v^\top F_1 \pi^1
+
(1-q),v^\top F_2 \pi^2.
]



### Intuition

If 30% of your customers are Group 1 and 70% are Group 2, then your total revenue should reflect that mix. The policy is evaluated on the population as a whole, not group-by-group in isolation.

---

## 9) Procedural unfairness: fairness of the offer process

The paper defines procedural unfairness as

[
U(\pi)=|v^\top \pi^1 - v^\top \pi^2|.
]

This is the absolute difference in the **expected proposed prices** between the two groups. 

### Meaning

If Group 1 is offered a higher average price than Group 2, then the process is procedurally unfair.

### Why this is easy to compute

This depends only on the policy (\pi), not on the unknown demand (F_e). So it is tractable.

---

## 10) Substantive unfairness: fairness of the outcome

The paper defines substantive unfairness as

[
S(\pi;F_1,F_2)
==============

\left|
\frac{v^\top F_1 \pi^1}{\mathbf{1}^\top F_1 \pi^1}
--------------------------------------------------

\frac{v^\top F_2 \pi^2}{\mathbf{1}^\top F_2 \pi^2}
\right|.
]

This is the absolute difference in the **expected accepted prices** between the two groups. 

### Meaning

Even if both groups are offered fair average prices, the buyers who accept in each group may still end up paying different average amounts. This metric captures that second-layer fairness.

### Why it is harder

Unlike procedural fairness, substantive fairness depends on the unknown acceptance probabilities (F_1) and (F_2). That makes it difficult to enforce directly.

---

## 11) The optimization goal

The ideal fair policy is:

[
\pi^*
=====

\arg\max_{\pi=(\pi^1,\pi^2)\in\Pi}
R(\pi;F_1,F_2)
\quad
\text{s.t. }
U(\pi)=0,;
S(\pi;F_1,F_2)=0.
]



### What this means in words

Among all policies that are perfectly fair in both senses, choose the one with the highest expected revenue.

### Why this is hard

The constraint on substantive fairness contains ratios, so the feasible set is non-convex and non-linear. This is one reason the paper later needs a careful algorithm rather than a direct closed-form solution.

---

## 12) Why the optimal policy can still be feasible

The paper notes that feasibility is not empty: a trivial policy exists where both groups are always offered the same fixed price. That automatically makes procedural unfairness zero, and with the right price it can also satisfy substantive fairness. 

### But there is a catch

That feasible solution may be far from revenue-optimal. So the existence of a feasible point does not mean the problem is easy or economically attractive.

---

## 13) Regret: how learning performance is measured

The cumulative regret of an algorithm (A) is

[
\mathrm{Reg}_T(A)
=================

\sum_{t=1}^T
\left(
R(\pi^*;F_1,F_2)-R(\pi_t;F_1,F_2)
\right).
]

Here (\pi_t) is the policy used at time (t). 

### Intuition

At each round, compare the revenue of the policy you used with the revenue of the best fair policy. Add up the gaps over time.

If this total grows sublinearly, then average regret per round goes to zero, which means the algorithm is learning well.

### A subtle point

The paper allows the per-round regret expression to be negative if the current policy violates fairness and earns more than the fair optimum. So regret is measured against the **fair benchmark**, not against any arbitrary unconstrained policy. That keeps the objective aligned with the fairness-constrained problem.

---

## 14) Cumulative substantive unfairness

The paper also defines the cumulative substantive unfairness as

[
S_T(A)=\sum_{t=1}^T S(\pi_t;F_1,F_2).
]



### Meaning

This measures how much total fairness violation accumulates over time.

Even if each round is only slightly unfair, small violations can add up. The paper wants this quantity to stay controlled.

---

## 15) The algorithmic target

The paper’s goal is not to guarantee perfect substantive fairness at every round, because that would be unrealistic when the demand curves are unknown. Instead, it aims for:

* optimal-order regret,
* zero procedural unfairness,
* minimal cumulative substantive unfairness. 

This is the key online-learning formulation.

---

## 16) Technical assumptions, explained carefully

### Assumption 1: acceptance probability is bounded away from zero at the highest price

The paper assumes there is a constant (F_{\min}>0) such that

[
F_e(d)\ge F_{\min}, \quad e=1,2.
]



#### Why this is needed

If the highest price is never accepted, then some acceptance rates can become zero, and the accepted-price ratio can break down. The assumption ensures the denominator in the substantive fairness term stays well-defined and stable.

### Assumption 2: the number of prices is controlled

The paper treats (d) as independent of (T), and assumes

[
d=O(T^{1/3}).
]



#### Why this appears

This is mainly a technical condition used in the regret analysis. It keeps the price-grid complexity from growing too quickly relative to the horizon.

---

## 17) A small example to make the formulas concrete

Suppose there are only two prices:

[
V={0.6, 1.0}.
]

Suppose Group 1 uses policy

[
\pi^1=(0.75, 0.25),
]

so it offers 0.6 with probability 0.75 and 1.0 with probability 0.25.

If the acceptance probabilities are

[
F_1(1)=0.9,\quad F_1(2)=0.5,
]

then:

### Expected proposed price

[
0.6(0.75)+1(0.25)=0.7.
]

### Expected revenue

[
0.6(0.75)(0.9)+1(0.25)(0.5)
=0.405+0.125=0.53.
]

### Acceptance rate

[
0.75(0.9)+0.25(0.5)=0.675+0.125=0.8.
]

### Expected accepted price

[
\frac{0.53}{0.8}=0.6625.
]

So:

* the seller proposes average price (0.7),
* buyers who accept pay average price (0.6625),
* and total expected revenue is (0.53).

This example shows why proposed price and accepted price are different objects.

---

## 18) The big picture of the problem setup

All of the notation in this section serves one purpose:

* define how revenue is computed,
* define how fairness is measured,
* define what “best” means under fairness,
* define what the learning algorithm is trying to approximate over time. 

Once these definitions are clear, the later algorithm becomes much easier to understand. The next section will build the actual procedure: how the paper learns acceptance probabilities, forms candidate policies, and eliminates bad ones epoch by epoch.


## Algorithm section: how FPA works

The algorithm is the paper’s main engineering idea. It is built to do three things at once:

* **learn** the acceptance probabilities,
* **keep procedural fairness exactly zero**,
* **keep substantive unfairness as small as possible while still learning well**.

The paper does this with an **epoch-based policy elimination scheme**. Think of it as: learn a little, estimate the market, throw away bad policies, and repeat with a better candidate set.

---

## 1) What the algorithm receives as input

FPA takes:

* the time horizon (T),
* the finite price set (V={v_1,\dots,v_d}),
* an error probability (\epsilon),
* a constant (L),
* and the group proportion (q).

### What these mean

* (T): how long the selling process lasts.
* (V): the only prices you are allowed to use.
* (\epsilon): how confident you want the guarantees to be.
* (L): a technical constant controlling how much reward tolerance is needed when fairness is relaxed slightly.
* (q): how often Group 1 appears in the population.

So the algorithm is not guessing the customer mix from scratch; it is told the population share (q), but it still does **not** know the demand curves (F_1, F_2).

---

## 2) The “before epochs” phase

Before the main epoch loop starts, the algorithm repeatedly offers the **highest price** (v_d) for about

[
\tau_0 = O(\log T)
]

rounds. The purpose is to estimate a **lower bound** on the minimum acceptance probability, denoted (\hat F_{\min}).

### Why the highest price?

Because if even the highest price is accepted sometimes, then all lower prices should be at least as acceptable. That gives a safe baseline. The paper uses this to prevent division by values close to zero later in the substantive-fairness ratio.

### What is being estimated?

For each group (e), the algorithm counts:

* (M_{0,e}): number of times that group appeared,
* (N_{0,e}): number of times the highest price was accepted.

Then it computes a conservative estimate like

[
\hat F_{\min}=\min\left(\frac{N_{0,1}}{2M_{0,1}},\frac{N_{0,2}}{2M_{0,2}}\right).
]

The factor (1/2) is a safety margin. It is not trying to estimate the exact acceptance rate; it is trying to get a **lower bound that is unlikely to be too optimistic**.

### Intuition

If you later divide by acceptance probabilities when computing accepted prices, a too-small denominator would make the estimates unstable. So the algorithm first builds a floor under those probabilities.

---

## 3) The main loop: epochs

After the warm-up, the algorithm splits time into epochs (k=1,2,\dots). Each epoch has a length (\tau_k) that grows like a doubling schedule. The paper writes this as a doubling-epoch structure.

### Why epochs?

Because in online learning, you want time blocks long enough to estimate quantities reliably. If you changed the policy every round, the estimates would be too noisy. If epochs get longer over time, the algorithm can exploit more as it learns more.

### Intuition

Early on, the model is uncertain, so it should explore more. Later, once it has learned enough, it should exploit better policies more heavily. The epoch structure is exactly how the paper balances those two goals.

---

## 4) Candidate policy set (\Pi_k)

At the start of epoch (k), the algorithm keeps a set (\Pi_k) of candidate policies. Initially, (\Pi_1) is the set of all procedurally fair policies:

[
\Pi_1=\Pi={\pi=(\pi^1,\pi^2):U(\pi)=0}.
]

This means the policy space is already restricted so that the two groups have equal expected proposed prices.

### Important point

The algorithm does **not** search over all possible price rules forever. It gradually eliminates policies that look too unfair or too suboptimal. That is why it is called a **policy-elimination** method.

---

## 5) Good-and-exploratory policies

This is one of the cleverest parts of the paper.

A policy may be good for revenue but bad for learning, because it might almost never try some prices. The paper therefore builds a smaller set (A_k) of **good-and-exploratory policies**.

### What does “good-and-exploratory” mean?

For each group (e) and each price (v_i), the algorithm looks inside the remaining candidate set (\Pi_k) and finds a policy that maximizes the probability of proposing (v_i) to that group. If that probability is large enough, the policy is kept in (A_k). If not, that price index may be dropped from further exploration.

### Why is this needed?

Suppose a policy almost never proposes price (v_i). Then the algorithm cannot estimate the acceptance probability of (v_i) well. Since the whole revenue and fairness calculation depends on acceptance probabilities, those prices need enough exploration exposure.

### Threshold idea

The paper uses a threshold of about (1/\sqrt{T}). If a price is selected with probability below that level even by the best remaining policy, the algorithm treats it as not worth exploring further.

### Intuition

This is like saying: “Do not waste time on actions that are so rare they will not meaningfully affect learning or the final policy.”

---

## 6) Estimating acceptance probabilities inside an epoch

Once (A_k) is built, the algorithm runs each policy in (A_k) for a batch of rounds, roughly (\tau_k/|A_k|) rounds per policy. During this phase it records:

* (M_{k,e}(i)): how many times price (v_i) was proposed to group (e),
* (N_{k,e}(i)): how many of those proposals were accepted.

Then it estimates the acceptance probability by

[
\bar F_{k,e}(i)=\max\left(\frac{N_{k,e}(i)}{M_{k,e}(i)},\hat F_{\min}\right)
]

for the prices still being tracked. If a price is too rarely observed, the algorithm falls back to (\hat F_{\min}).

### Why take a max with (\hat F_{\min})?

Because raw empirical ratios can be unstable when sample sizes are small. If (M_{k,e}(i)) is small, the estimate might be artificially tiny or noisy. Flooring it at (\hat F_{\min}) keeps the optimization stable.

### Simple example

Suppose for Group 1 and price (v_2):

* the price was proposed 20 times,
* accepted 12 times.

Then the raw estimate is (12/20=0.6). If (\hat F_{\min}=0.2), the algorithm stores (0.6). If instead the sample ratio were (0.05), it would store (0.2), not (0.05). That is a conservative safety step.

---

## 7) Building the empirical revenue and fairness model

Once (\bar F_{k,1}) and (\bar F_{k,2}) are estimated, the algorithm constructs diagonal matrices

[
\hat F_{k,e}=\mathrm{diag}(\bar F_{k,e}(1),\dots,\bar F_{k,e}(d)).
]

These are then plugged into the revenue and fairness formulas as if they were the true demand curves. That gives the algorithm a **learned surrogate world** inside epoch (k).

### Why this matters

The algorithm cannot optimize with respect to the real (F_1,F_2) because they are unknown. So it optimizes the estimated versions, then controls the error through concentration and conservative elimination.

---

## 8) Algorithm 2: empirical optimal policy

Now comes the core optimization step.

The paper wants the best policy under the estimated model, subject to fairness. But the substantive-fairness constraint has a ratio, so the problem is not directly convex. The trick is to introduce a scalar variable (w), which represents the common accepted-price level. The paper then searches over possible values of (w) in steps of

[
\epsilon = \delta_{k,s}/2.
]

For each fixed (w), the problem becomes a linear program.

### Why does fixing (w) help?

The substantive fairness condition is based on accepted prices. If you say “the accepted price level for Group 1 should be (w), and Group 2 should be close to (w),” then the ratio constraint becomes equivalent to linear inequalities once (w) is fixed.

That is the central mathematical simplification.

### What Algorithm 2 does

For each (w_\ell):

* solve the LP maximizing estimated revenue,
* subject to:

  * procedural fairness on the accepted-price level for Group 1,
  * approximate substantive fairness for Group 2.

Then it keeps the best solution over all tested (w_\ell).

### Intuition

This is like scanning a one-dimensional axis of “possible fair accepted-price levels” and, at each level, finding the best policy consistent with that level.

---

## 9) Why the search over (w) is only approximate

The paper does not search over all real values of (w). It uses a grid with step (\delta_{k,s}/2). This is enough because the objective and constraints are Lipschitz under the paper’s assumption that acceptance probabilities are bounded away from zero. So discretization error stays controlled.

### First-principles intuition

If a function changes smoothly, you do not need to test every point on a continuum. Checking a fine enough grid is enough to get near-optimal performance. The paper’s Lipschitz argument is what justifies this.

---

## 10) Algorithm 3: policy elimination

After finding the empirical optimal policy, the algorithm prunes the candidate set. A policy is removed if it is:

* too unfair, or
* too far below the empirical optimum in revenue.

The elimination rule is:

[
S(\pi,\hat F_{k,1},\hat F_{k,2})>\delta_{k,s}
]

or

[
R(\pi,\hat F_{k,1},\hat F_{k,2})
<
R(\hat \pi_{k,*},\hat F_{k,1},\hat F_{k,2})
-\delta_{k,r}-L\delta_{k,s}.
]

### Why the second subtraction has two pieces

* (\delta_{k,r}) handles estimation error in revenue.
* (L\delta_{k,s}) handles the fact that allowing a small fairness slack can change the optimal revenue by a bounded amount. 

### Why this is conservative

The point is not to remove every policy that looks slightly bad. The point is to remove only those that are confidently bad, so the true optimal fair policy is still likely to remain in the candidate set.

### Intuition

This is a statistical version of “keep only policies that are plausible winners.”

---

## 11) Why the algorithm is efficient

The paper emphasizes that FPA is oracle-efficient. The reason is:

* the number of epochs is only logarithmic,
* the search over (w) is a one-dimensional grid,
* each fixed-(w) problem is a linear program.

So even though the original optimization is non-convex, the algorithm converts it into a sequence of manageable LPs.

### Big picture

This is a standard but powerful research strategy:
turn a hard constrained online problem into a sequence of easier offline convex subproblems.

---

## 12) How the whole algorithm fits together

You can think of the full loop like this:

1. **Warm up**: estimate a safe lower bound on acceptance.
2. **Form candidate policies**: keep policies that are still plausible and fair.
3. **Select exploratory policies**: ensure enough coverage of prices.
4. **Estimate acceptance probabilities** from observed buys.
5. **Solve an estimated fair optimization problem**.
6. **Eliminate bad policies**.
7. **Repeat with a tighter candidate set and longer epoch.**

That is the entire algorithmic logic.

---

## 13) A small intuition example

Imagine only three prices:

[
V={0.5,0.7,1.0}
]

and two groups.

Suppose after an epoch, the estimated acceptance probabilities look like this:

* Group 1: (F_1=(0.9,0.6,0.3))
* Group 2: (F_2=(0.95,0.5,0.2))

Now a policy that puts too much mass on price (1.0) might still be fair in proposed-price average, but its accepted-price average may diverge across groups because Group 2 rejects (1.0) more often. So the empirical optimal step will search for a randomized mixture that keeps accepted prices aligned while maximizing revenue.

Then elimination removes policies that cannot compete under those learned estimates.

This is the paper’s logic in miniature.

---

## 14) What to remember for implementation

If you later code this correctly, the key objects are:

* price set (V),
* policy distributions (\pi^1,\pi^2),
* acceptance estimates (\hat F_{k,e}),
* candidate policy set (\Pi_k),
* exploratory set (A_k),
* epoch lengths (\tau_k),
* tolerance levels (\delta_{k,r},\delta_{k,s}),
* and the (w)-grid search used in Algorithm 2 and 3.

Everything else is bookkeeping around these core pieces.

The next natural step is the **regret and unfairness analysis**, because that is where the paper proves why this whole construction is theoretically justified.



## Regret and Unfairness Analysis: what the paper proves and why it matters

This section is where the paper stops describing the algorithm and starts proving that it is actually good. The main question is:

**Does FPA really learn a nearly optimal fair pricing rule, and is that rate the best possible?**

The answer the paper gives is:

* **Upper bound:** FPA achieves low regret and low substantive unfairness.
* **Lower bound:** no algorithm can do fundamentally better in the worst case, up to small logarithmic factors.

So this section is not just “analysis.” It is the paper’s proof that its algorithm sits near the theoretical limit.

---

## 1) The main upper bound result

Theorem 6 says that FPA guarantees

[
O!\left(\sqrt{T}, d^{3/2}\log d \log T/\epsilon \right)
]

regret, with **zero procedural unfairness** and the same order of **substantive unfairness**, with probability at least (1-\epsilon). 

### What this means in plain language

Over a long horizon (T):

* total revenue loss relative to the best fair policy grows only like (\sqrt{T})-type growth,
* the average loss per round therefore goes to zero,
* the algorithm never violates procedural fairness,
* and substantive unfairness also becomes small at the same asymptotic rate.

So the algorithm is both **learning** and **staying fair**.

---

## 2) Why (\sqrt{T}) is the “right” kind of result

A total regret of order (\sqrt{T}) means average regret is

[
\frac{\sqrt{T}}{T}=\frac{1}{\sqrt{T}},
]

which vanishes as (T) grows. This is the standard no-regret regime in online learning.

So the theorem says the algorithm does not merely work in a fixed finite run; it becomes more accurate as data accumulates. That is the correct notion of asymptotic learning. 

---

## 3) How the proof works at a high level

The paper says the proof is by **induction over epochs**. That is the key structural idea. 

Here is the logic in simple steps.

### Step 1: Start with the true optimal fair policy inside the candidate set

At epoch (k=1), the candidate set (\Pi_1) contains the optimal fair policy (\pi^*). This is the induction base case. 

### Step 2: Estimate demand from data

Inside epoch (k), the algorithm estimates the acceptance probabilities (F_e(i)) using the observed accept/reject feedback. Because it only has finite samples, these estimates have error. The paper controls that error using concentration inequalities. 

### Step 3: Convert estimation error into revenue and fairness error

Once the estimated acceptance probabilities are good enough, the paper shows that the estimated revenue and estimated substantive unfairness are close to their true values. This step is where the math matters most: the uncertainty in acceptance probabilities propagates into uncertainty in both objective and constraint.

### Step 4: Keep only policies that are confidently good

The algorithm then eliminates policies that appear too unfair or too low-revenue under the estimated model. The elimination is conservative, so the true optimal fair policy is not thrown out.

### Step 5: Preserve the induction hypothesis

The paper shows that if (\pi^*) is in (\Pi_k), then after elimination it still remains in (\Pi_{k+1}). That lets the argument repeat for the next epoch. 

### Step 6: Sum across epochs

Finally, because epoch lengths grow in a doubling fashion, the total regret and total unfairness across all epochs add up to the stated bound.

---

## 4) Why the epoch structure gives the right scaling

The doubling-epoch design is not cosmetic. It is what makes the proof and the learning curve work together.

Early epochs are short, so the algorithm explores more and learns fast. Later epochs are longer, so once the policy set has narrowed, the algorithm can exploit a better candidate for a longer time. This is why the paper can both learn and keep fairness under control.

### Intuition with a toy picture

Suppose the seller is unsure between two pricing mixtures:

* Policy A: slightly more revenue, but maybe unfair.
* Policy B: slightly safer and more fair.

If you only test them for a handful of rounds, the data is too noisy. But if you gradually make epochs longer, you can estimate which one is truly better and remove the worse one. That is the logic behind the doubling schedule.

---

## 5) Why substantive unfairness is not zero

The paper is careful: procedural unfairness is exactly zero, but substantive unfairness is only bounded, not eliminated. Why?

Because substantive fairness depends on the unknown acceptance probabilities (F_1,F_2). You do not observe those directly. You can only estimate them from binary feedback, so some residual mismatch is unavoidable during learning.

This is why the paper’s guarantee is asymptotic: the unfairness becomes small over time, but it need not be exactly zero at every round.

---

## 6) The remark after Theorem 6: what it really says

Remark 7 says the algorithm guarantees (O(\sqrt{T}\log\log T)) regret and unfairness simultaneously, and that this matches the generic (O(\sqrt{1/T})) estimation error when averaged over time. It also says these upper bounds are tight up to (O(\log\log T)) factors. 

### What that means intuitively

The algorithm is limited by how fast it can estimate unknown acceptance probabilities from noisy binary feedback. Since estimation error typically scales like (1/\sqrt{n}), and (n) grows with (T), a (\sqrt{T})-type cumulative bound is exactly what you would expect in a hard online learning problem.

So the theorem is saying: the algorithm is not leaving obvious performance on the table.

---

## 7) Regret lower bound: why no algorithm can beat (\Omega(\sqrt{dT}))

Theorem 8 states that any algorithm must suffer at least

[
\Omega(\sqrt{dT})
]

regret under the stated two-group fair pricing problem, assuming (d \le T^{1/3}). 

### Why is this true?

The paper reduces the fair two-group problem to the ordinary online pricing problem by taking the special case where the two groups are identical:

[
F_1(i)=F_2(i), \quad \forall i.
]

Then fairness is not really adding structure, because any policy with (\pi^1=\pi^2) is both procedurally and substantively fair. The problem collapses to standard online pricing with finitely many prices, for which bandit-style lower bounds already exist.

### Intuition

If the groups are identical, the algorithm is essentially just learning the best price from accept/reject feedback. That is a classic exploration-exploitation problem, and it cannot be solved faster than the known bandit lower bound.

So the fair problem cannot be easier than the ordinary pricing problem.

---

## 8) What the (d) dependence is telling you

The (\Omega(\sqrt{dT})) bound says the difficulty grows with:

* (T): because you need time to learn,
* (d): because more candidate prices mean more actions to distinguish. 

### Intuition

If there are only two prices, learning is easier. If there are many prices, the seller must determine which prices are good and which are not, and each extra action adds uncertainty. So the lower bound scales with both horizon and action-space size.

---

## 9) The fairness lower bound: why optimal regret may force unfairness

This is the paper’s deepest result.

Theorem 9 says that any algorithm with cumulative regret (C_x\sqrt{T}) and zero procedural unfairness must suffer at least (C_u\sqrt{T}) substantive unfairness. 

### What is the message?

Even if you are very good at maximizing revenue, you may still be forced to incur some substantive unfairness while learning.

That is not because the algorithm is bad. It is because the information needed to make the fair decision is itself hard to extract quickly.

---

## 10) Why the paper needs two nearby problem settings

To prove Theorem 9, the paper constructs two nearly identical environments:

* Setting A: the original Example 1.
* Setting B: the same example, except some acceptance probabilities are perturbed by a small amount (\zeta).

The idea is that the optimal fair policy differs between these two settings, but the difference is hard to detect from samples.

### Why does this matter?

If two environments are very close statistically, then any algorithm needs many rounds to tell them apart. During that time, it must either:

* act cautiously and give up revenue, or
* commit to one policy and risk unfairness in the wrong environment.

That is the standard indistinguishability argument.

---

## 11) A simple intuition for the indistinguishability argument

Imagine two doors that look almost identical.

* Behind Door A, the fair action is to mix prices one way.
* Behind Door B, the fair action is to mix prices slightly differently.

But the observed accept/reject feedback is almost the same under both doors. So to know which door you are in, you need lots of data.

If you do not wait long enough, you may choose the wrong mixture and become substantively unfair. If you wait too long, you lose revenue. That is the tradeoff the theorem formalizes.

---

## 12) Why the theorem uses “optimal algorithms” in the statement

The paper’s lower bound is not just saying “some bad algorithms are unfair.” It is saying that even algorithms with regret as low as (O(\sqrt{T})) cannot avoid substantive unfairness altogether. 

That is important because it shows the unfairness is not a design flaw; it is an information-theoretic necessity under the online-learning constraints.

---

## 13) The tradeoff the lower bound reveals

The paper explicitly says:

* if you try not to distinguish the two settings, you lose revenue or fairness,
* if you try to distinguish them, you need time and data,
* and because the regret budget is tight, the unavoidable consequence is (\Omega(\sqrt{T})) substantive unfairness.

### In plain words

You cannot learn the right fair randomization instantly. The learning process itself creates some unfairness.

That is a subtle and very realistic point.

---

## 14) Why this is not the same as ordinary regret lower bounds

The paper stresses that Theorem 9 is different from standard regret lower bounds because it also requires the algorithm to be **optimal** or near-optimal in regret. 

So the lower bound is not merely saying “learning is hard.”

It is saying:

> If you insist on learning fast enough to get optimal regret, then some substantive unfairness is unavoidable.

This is a stronger and more interesting statement.

---

## 15) How to read the upper and lower bounds together

Taken together, Theorems 6, 8, and 9 say:

* FPA achieves the correct rate.
* You cannot do substantially better in regret.
* You cannot simultaneously have optimal regret and zero substantive unfairness in general.

That is why the paper claims near-optimality in both revenue and fairness.

---

## 16) The role of Assumption 1 in the analysis

The paper later explains that the lower bound (F_{\min}>0) is important because it ensures the substantive unfairness function is well-behaved and Lipschitz. Without that, the accepted-price ratio could become unstable or even undefined if some group almost never accepts.

### Intuition

If acceptance probability is near zero, then “average accepted price” becomes a fragile quantity. Small estimation errors can create huge swings. Bounding acceptance away from zero prevents that blow-up.

---

## 17) What the discussion adds to the analysis

The discussion section says that the algorithm cannot guarantee perfect any-time substantive fairness, but it can guarantee asymptotic fairness. It also notes a practical perspective: longer epoch batches help customers experience fairness more consistently, instead of seeing a rapidly changing policy every round.

This is a nice bridge between theory and practice:
the same epoch structure that helps the proof also makes the fairness feel more stable in a real system.

---

## 18) A compact intuition summary

Here is the whole analysis in one chain:

1. The algorithm learns acceptance probabilities from binary feedback.
2. Those estimates get better over time.
3. Better estimates allow better fair policy selection.
4. The candidate set is pruned conservatively so the true optimum survives.
5. Doubling epochs ensure enough data for each refinement step.
6. This yields (\tilde O(\sqrt{T}))-type regret and unfairness.
7. Standard pricing lower bounds show you cannot beat the regret order.
8. A two-environment indistinguishability argument shows that optimal regret forces some substantive unfairness.

---

## 19) What you should remember for implementation

For coding this later, the analysis tells you exactly what must be preserved:

* the candidate set must shrink conservatively,
* the optimal fair policy must remain feasible,
* acceptance probabilities must be estimated carefully,
* the epoch schedule must grow,
* the fairness-revenue tradeoff must be checked by oracle solving,
* and the evaluation must measure both regret and substantive unfairness.

If any of those pieces is broken, the theory no longer applies.

The next section worth studying is the **Discussion**, because it explains the paper’s assumptions, limitations, and possible extensions, which are very important if you want to implement this correctly on your own.

## Discussion section: what the paper is really saying beyond the theorem

This section is not about proving the main algorithm again. It is about the **limits, assumptions, extensions, and interpretation** of the result. In other words, it answers:

* What did we assume to make the theory work?
* What does the algorithm mean from a human/customer perspective?
* What happens if we relax the fairness target?
* Can the method extend to harder settings?

This section is especially important if you want to implement the paper correctly, because it tells you which parts are essential and which parts are technical conveniences.

---

## 1) The technical assumption (F_{\min}>0): why the paper needs it

The paper assumes that for both groups, the acceptance probability of the highest price is bounded below:

[
F_e(d)\ge F_{\min}>0,\quad e=1,2.
]

The discussion explains that this assumption is stronger than what one would ideally want, but it is useful for two reasons. First, if some proposed price is never accepted, then the expected accepted price can become undefined or unstable. Second, a positive lower bound lets the substantive unfairness function behave nicely, in particular making it Lipschitz.

### Why does zero acceptance break things?

Substantive unfairness uses

[
\frac{v^\top F_e\pi^e}{\mathbf{1}^\top F_e\pi^e}.
]

The denominator is the acceptance rate. If that is zero, then the “average accepted price” is not defined. Even if it is very small, the ratio becomes numerically unstable.

### Simple example

Suppose Group 1 never accepts the highest price, but Group 2 does. Then for a policy that puts a lot of mass on that price, Group 1’s accepted-price average may be undefined, while Group 2’s is well-defined. The paper asks a natural question: is such a policy even meaningfully fair? Its answer is that the theory becomes messy unless we exclude this case.

### What the paper says about this assumption

The authors explicitly say that one could try to design an algorithm that works for (F_e(i)>0) in full generality, but that remains open. So this is not a philosophical claim; it is a technical limitation of the current proof framework.

### Practical interpretation

When implementing the paper, do not treat (F_{\min}) as a cosmetic constant. It is one of the pillars that keeps the fairness ratio well-defined and the optimization stable.

---

## 2) “Feelings of fairness”: why batch length matters

This is a very interesting and very practical part of the discussion.

The paper notes that FPA does not change the policy at every single round. Instead, it runs each exploratory policy for a **continuous batch** of rounds of length roughly

[
\frac{\tau_k}{|A_k|} = \Omega(\sqrt{T},2^k).
]

The authors argue that this is long enough for customers to **experience** the fairness. Customers can compare their average offered/accepted prices with others in the other group over a meaningful time window.

### Why does this matter?

Fairness in practice is not only a mathematical property. It is also a **perceived consistency** property.

If an algorithm changes policy every round, then customers may feel they are being treated by a random, unstable system, even if the long-run averages are fair.

### Example

Suppose:

* on Monday, Group 1 gets one kind of policy,
* on Tuesday, another,
* on Wednesday, yet another.

Even if the long-run averages are fair, customers may never notice the balance because each individual’s comparison window is too short.

The paper is suggesting that fairness should not be so “fast-changing” that it becomes invisible to the people it is supposed to protect.

### Why FPA is better in that sense

Because FPA uses batches, customers see a stable pricing regime long enough to compare. That is not just a cosmetic implementation detail. It is part of the paper’s fairness story.

---

## 3) Relaxing substantive fairness with a tolerance (\delta)

The paper then asks: do we really need perfect substantive fairness at every point? Maybe a small tolerance is more realistic.

So it introduces

[
\pi_{\delta,*}
==============

\arg\max_{\pi\in\Pi} R(\pi;F_1,F_2)
\quad \text{s.t.}\quad
U(\pi)=0,;
S(\pi;F_1,F_2)\le \delta.
]

This means:

* procedural fairness must still be exact,
* substantive fairness only needs to be within tolerance (\delta).

### Why is this useful?

In real systems, “exact equality” is often too strict. A small acceptable gap may be operationally enough, especially if the fairness quantity is an average over many customers.

### What changes mathematically?

The unfairness measure becomes

[
\max{0, S(\pi;F_1,F_2)-\delta}.
]

So if the accepted-price gap is below (\delta), it is treated as zero violation.

### The paper’s key observation

The authors show a tradeoff:

[
R(\pi^*) \le R(\pi_{\delta,*}),
\quad
R(\pi^*) \ge R(\pi_{\delta,*}) - L\delta.
]

This means the strictly fair optimum (\pi^*) and the relaxed fair optimum (\pi_{\delta,*}) are close in revenue, up to a controlled linear loss in (\delta).

### Intuition

Allowing a tiny bit of fairness slack should not drastically improve revenue, otherwise the fairness constraint would be very brittle. The lemma formalizes that the improvement is bounded.

### Why this matters for implementation

This relaxation is a natural direction if you want a practical version of the algorithm. In a real system, you may prefer a policy that is almost fair every round instead of exactly fair only asymptotically.

---

## 4) What the paper says about the rate when (\delta) changes

The discussion then explores how regret and unfairness might behave under different (\delta) regimes.

The paper states:

* If (\delta=0), both regret and unfairness rates are (\Theta(\sqrt{T})).
* If (\delta\ge 1), the problem becomes essentially unconstrained, so optimal regret is still (\Theta(\sqrt{T})) and unfairness can be (0).
* If (\delta=O(1/\sqrt{T})), the authors believe one may still get (O(\sqrt{T})) regret and unfairness.
* For (\delta > 1/\sqrt{T}), they conjecture the unfairness could scale like
  [
  \Theta!\left(\frac{1}{\sqrt{T},\delta}\right).
  ]

### How to interpret this

As you loosen the fairness constraint, the system gets more freedom, so revenue can improve. But the exact tradeoff curve is not fully characterized.

### Practical reading

This is the paper saying: “We know the strict case well. We have partial intuition for the relaxed case, but the full theory is open.”

For implementation, this means you should be careful not to overclaim what happens when you loosen the substantive fairness threshold.

---

## 5) Trade-off between procedural and substantive fairness

The paper makes a conjecture that there is probably **not much to gain** by intentionally trading substantive fairness for procedural fairness. The reason is conceptual:

* the substantive unfairness lower bound comes from the difficulty of distinguishing two very similar environments,
* making prices intentionally unfair does not necessarily make those environments easier to distinguish.

### Why is that plausible?

Suppose two worlds differ only slightly in acceptance probabilities. You might think that by charging different prices more aggressively, you could learn the difference faster. The authors argue that this probably does not help enough, because the information bottleneck is still the same: you only see binary accept/reject feedback.

So the paper suggests that procedural unfairness is not a useful “resource” to spend in order to eliminate substantive unfairness.

### Intuition

The fairness violation is not the thing slowing learning down. The problem is the hidden demand ambiguity. So violating one fairness notion does not magically resolve the learning problem.

---

## 6) Continuous pricing space: why the paper restricts to a finite set

The paper then asks: what if prices were allowed to vary continuously in ([0,1]) instead of being chosen from a finite grid (V={v_1,\dots,v_d})?

Its answer is: the problem becomes much harder. The optimal policy may be a pair of continuous distributions, and even if customer valuations only lie on a finite set, the fairness-constrained optimum may still not live on that same set. 

### Why is this hard?

In the finite case, the policy is just a probability vector over (d) prices. In the continuous case, the policy is a distribution over infinitely many prices. That is a much bigger space.

The paper says methods like continuous discretization might help, but they would likely create exponential time complexity. 

### Example

If you allow every price in ([0,1]), then instead of choosing from 3 or 10 points, you are choosing from a continuum. Even if the answer “looks” simple, the optimization landscape is much more complicated because the fairness constraint interacts with the demand curve in a nonlinear way.

### Why this matters for your implementation

This is the reason the paper’s algorithm is designed for a finite candidate price set. If you later want a more realistic implementation, you will usually discretize the price space first.

---

## 7) Extending from two groups to many groups

The paper is explicit that it studies two groups only to keep the presentation clean. But it suggests a multi-group extension.

A simple extension would define multi-group unfairness as the sum of pairwise unfairness across all (O(G^2)) pairs of groups. Then the algorithm would need to be adjusted, for instance by lengthening epochs by a factor of about (G/2). The paper says this would worsen the upper regret bound by a factor on the order of (G^3), while still leaving the dependence on (T) near-optimal up to iterated logs. 

### Why is this nontrivial?

Because with many groups:

* there are more fairness constraints,
* more pairwise comparisons,
* more places where the accepted-price average can differ,
* and more data needed to estimate all group-specific acceptance probabilities well.

### A simple example

If there are three groups instead of two, you now need fairness not just between Group 1 and Group 2, but also between Group 1 and Group 3, and Group 2 and Group 3. The number of fairness relations grows quickly.

### Why the paper calls this controversial

Pairwise fairness in multi-group settings is not always the right notion. In some real systems, pairwise equality may overconstrain the policy or fail to reflect the true equity objective. The paper acknowledges that this extension is not fully settled and would need further discussion. 

---

## 8) What the discussion reveals about the paper’s philosophy

The discussion is not just a list of caveats. It reveals the paper’s main philosophy:

1. **Fairness should be measurable and operational.**
2. **Fairness should be stable enough to matter to actual users.**
3. **Perfect fairness may be too strong, but asymptotic fairness is meaningful.**
4. **The theory is cleanest for finite prices and two groups, but the ideas may extend.**
5. **Some tradeoffs are real, but some expected tradeoffs may not actually exist.**

---

## 9) A concrete way to think about the main discussion ideas

### If (F_{\min}) is too small

Then fairness ratios become unstable. This is like dividing by a nearly zero acceptance rate.

### If epochs are too short

Customers see a changing policy and fairness becomes hard to perceive.

### If (\delta) is larger

You allow a small fairness gap and may get a slightly better revenue-fairness tradeoff.

### If the price space is continuous

You need a much more complex optimization scheme.

### If there are many groups

The fairness constraint scales badly unless you redesign the metric carefully.

---

## 10) Why this section matters for implementation

If you want to implement the paper “professionally,” this section tells you what not to miss:

* include a safeguard for very low acceptance rates,
* use batch/epoch behavior instead of per-round random switching,
* decide whether you want exact substantive fairness or a relaxed (\delta)-fair version,
* restrict to finite price grids first,
* test two-group behavior before attempting multi-group generalization.

---

## 11) Bottom line

The discussion section says that the paper’s algorithm is powerful but still bounded by real structural issues:

* the fairness ratio needs nonzero acceptance,
* fairness is more meaningful when experienced over a time window,
* the strict fair optimum can be relaxed at controlled revenue cost,
* procedural and substantive fairness may not be a meaningful tradeoff pair,
* continuous prices and many groups make the problem much harder.

That is the real “take-home” from the discussion: the paper is not just proving a theorem; it is carving out the boundary of what is feasible in fair dynamic pricing.

The natural next step is the **Conclusion**, where the authors compress all of this into the final research claims and the open problem direction.


Below is the algorithmic core of the paper, rewritten as a **working mental model** for implementation.

The paper is built around one central idea: **learn fair pricing by repeatedly estimating demand, solving a constrained optimization problem on the estimate, and then eliminating policies that are now provably implausible**. The whole method is an epoch-based refinement loop. 

---

# 1) The big picture first

The algorithm has to do three things at once:

1. **Learn** unknown acceptance probabilities from buy / no-buy feedback.
2. **Keep procedural fairness exactly satisfied**.
3. **Keep substantive fairness controlled while still earning near-optimal revenue**. 

The reason the paper does not just solve one static optimization problem is that the true demand curves (F_1, F_2) are unknown. So the algorithm needs to alternate between:

* **exploration**: gather information,
* **estimation**: infer acceptance rates,
* **optimization**: find the best fair policy under the estimate,
* **elimination**: throw away policies that are now clearly bad. 

If you keep that loop in mind, the whole paper becomes much easier to implement.

---

# 2) Before Epochs

## What it does

Before the main learning loop starts, the algorithm offers the **highest price** (v_d) for a short warm-up period of length

[
\tau_0 = O(\log T).
]

It does this to estimate a **safe lower bound** on the smallest acceptance probability, called (\hat F_{\min}). 

## Why the highest price?

Because the highest price is the hardest one to accept. If even that price is accepted sometimes, then lower prices should be at least as acceptable. So this gives a conservative baseline. 

## What the counters mean

For each group (e\in{1,2}):

* (M_{0,e}) = how many times group (e) appeared,
* (N_{0,e}) = how many times the highest price was accepted. 

Then the paper uses something like

[
\hat F_{\min} = \min\left{\frac{N_{0,1}}{2M_{0,1}}, \frac{N_{0,2}}{2M_{0,2}}\right}.
]

The factor (1/2) is a safety margin. The idea is not to estimate the true minimum acceptance rate precisely, but to get a **conservative floor** that is unlikely to be too optimistic. 

## Why this matters in implementation

Later, substantive fairness uses a ratio:

[
\frac{\text{expected accepted revenue}}{\text{expected acceptance rate}}.
]

If the denominator becomes too small, the ratio becomes unstable. So this warm-up is a safety step, not a cosmetic one. 

---

# 3) Doubling Epochs

## What it does

After warm-up, time is split into epochs (k=1,2,\dots), where each epoch is longer than the previous one. The paper uses a doubling-style schedule. 

## Why epochs are needed

If you update the policy every single round, your estimates stay noisy and you never give a candidate policy enough time to reveal how it performs. Epochs solve this by giving the algorithm blocks of time where it can:

* collect enough samples,
* build reliable estimates,
* then refine the candidate policy set. 

## Why epoch length grows

Early on, you know very little, so you need exploration. Later, once many bad policies have been removed, you want to exploit the good ones longer. The increasing epoch length is how the algorithm gradually shifts from exploration to exploitation. 

## Implementation view

Think of the algorithm as repeatedly doing:

1. estimate demand,
2. optimize on the estimate,
3. keep only policies that survive,
4. make the next epoch longer so the estimates improve.

That is the entire logic of the doubling-epoch design. 

---

# 4) Good-and-Exploratory Policies

## Why this step exists

A good policy for revenue may ignore some prices almost completely. But if a price is never tried, you cannot estimate its acceptance probability well. So the paper constructs a smaller set of policies that are both:

* **good** for revenue,
* **exploratory** enough to learn the demand curve. 

## What the paper does

For each group (e) and each price index (i), it looks inside the current candidate policy set (\Pi_k) and finds a policy that maximizes the probability of proposing price (v_i) in group (e). If that maximum probability is large enough, the policy is kept in the exploratory set (A_k). If not, that price index may be removed from future consideration. 

## Why the threshold is around (1/\sqrt{T})

The paper uses a threshold of this scale because prices proposed less frequently than that are too rare to matter much asymptotically. In other words, if a price is almost never used by any plausible optimal policy, then spending learning budget on it is not worth much. 

## Intuition

This is like saying: do not waste samples on actions that are so rare they will not influence the final answer in a meaningful way.

## Example

Suppose price (v_3) is only proposed with probability (0.001) by every remaining policy. Even if you learned its acceptance probability perfectly, it would barely change revenue. So the algorithm is justified in dropping it from active exploration. 

---

# 5) Estimate Acceptance Probabilities

## What is being estimated

For every epoch (k), group (e), and price (v_i), the algorithm estimates

[
F_e(i) = \Pr(y^e \ge v_i).
]

This is the acceptance probability of price (v_i) for group (e). 

## How the estimation is done

During the epoch, the algorithm runs each exploratory policy (\pi \in A_k) for a batch of rounds and records:

* (M_{k,e}(i)): number of times price (v_i) was proposed to group (e),
* (N_{k,e}(i)): number of those proposals that were accepted. 

Then the empirical estimate is

[
\tilde F_{k,e}(i) = \frac{N_{k,e}(i)}{M_{k,e}(i)}
]

and the paper floors it at (\hat F_{\min}):

[
\bar F_{k,e}(i) = \max\left(\tilde F_{k,e}(i), \hat F_{\min}\right).
]

This prevents tiny or unstable estimates from breaking the optimization. 

## Why the floor matters

Suppose in a small batch, price (v_i) is accepted only once out of 100 by accident, so the empirical rate is (0.01). If the true rate is actually larger, the estimate is just noisy. Flooring avoids overreacting to such noise. 

## Implementation intuition

This step is just frequentist estimation with a conservative safety guard. It is simple, but it is one of the most important parts of the whole algorithm.

---

# 6) Policy Elimination

## What it does

At the end of each epoch, the algorithm removes policies that are now clearly bad. There are two reasons to remove a policy:

1. it looks too unfair,
2. it looks too low-revenue compared to the best current empirical policy. 

## The fairness elimination rule

If a policy has estimated substantive unfairness larger than the allowed tolerance (\delta_{k,s}), it is removed. 

## The revenue elimination rule

A policy is also removed if its estimated revenue is too far below the estimated best policy:

[
R(\pi,\hat F_{k,1},\hat F_{k,2})
<
R(\hat\pi_{k,*},\hat F_{k,1},\hat F_{k,2})
-\delta_{k,r}-L\delta_{k,s}.
]

The two subtractions mean:

* (\delta_{k,r}): allow for estimation error in revenue,
* (L\delta_{k,s}): allow for the fact that relaxing fairness a little can increase optimal revenue a bit. 

## Why this is conservative

The elimination is designed so the true optimal fair policy is not thrown away just because the current estimates are imperfect. So the algorithm removes only those policies that are confidently bad. 

## Example

Suppose the best empirical policy has estimated revenue (0.80), and a candidate policy has estimated revenue (0.74). If the tolerance terms sum to (0.03), then the cutoff is (0.77). Since (0.74<0.77), the policy is eliminated. If the tolerance were larger, it would survive. This is exactly how the conservative pruning works.

---

# 7) Computational Cost

The paper stresses that the algorithm is efficient enough to run because the hard optimization is converted into a small sequence of linear programs. 

## Why the original problem is hard

The substantive fairness constraint involves a ratio:

[
\frac{v^\top F_e \pi^e}{\mathbf{1}^\top F_e \pi^e}.
]

That ratio makes the direct optimization non-convex. So you cannot just hand it to a simple convex solver. 

## The trick

The paper introduces a scalar (w), which represents the common accepted-price level. Once (w) is fixed, the ratio constraint becomes linear:

[
v^\top F_e \pi^e = w,\mathbf{1}^\top F_e \pi^e.
]

That is a linear constraint in (\pi^e). So for each fixed (w), the problem becomes a linear program. 

## Why a grid search over (w) works

The paper searches over (w) in steps of (\epsilon = \delta_{k,s}/2). Because the functions are Lipschitz under the paper’s assumptions, the discretization error is controlled. 

## Implementation takeaway

For every epoch, you do not solve one huge nonlinear fairness problem. You solve many small LPs over a one-dimensional grid of (w)-values.

---

# 8) Algorithm 2: Empirical Optimal Oracle

This oracle finds the best policy under the current empirical model.

## What it does

For each candidate (w_\ell):

1. fix the target accepted-price level,
2. solve a linear program maximizing estimated revenue,
3. subject to procedural fairness and approximate substantive fairness,
4. keep the best solution across all tested (w_\ell). 

## Why fixing (w) helps

The fairness ratio becomes linear once you say:

[
\text{accepted-price level} = w.
]

Then the accepted-price condition becomes a linear equality or inequality. That is the core mathematical simplification. 

## Intuition

Imagine the accepted-price level is a dial. For each dial setting (w), ask: what is the most profitable fair policy? Then search over the dial settings.

## Example

Suppose (w=0.7). Then the optimizer looks for a price mixture that makes both groups’ expected accepted prices equal to (0.7), while also making expected proposed prices equal. Among all such policies, it picks the one with highest estimated revenue.

That is exactly what Algorithm 2 is doing. 

---

# 9) Algorithm 3: Policy Elimination Oracle

This oracle takes the empirical optimum from Algorithm 2 and removes policies that are now clearly inferior.

## What it does

For each (w_\ell), it solves a linear program over the candidate set and collects policies that satisfy:

* substantive fairness up to (\delta_{k,s}),
* revenue at least near the empirical best, up to (\delta_{k,r}+L\delta_{k,s}). 

The union of all such feasible policies becomes the next candidate set (\Pi_{k+1}). 

## Why this matters

This is the pruning step that makes the algorithm progressively sharper. The candidate set gets smaller, but the optimal policy is preserved with high probability. That is the whole point of the elimination scheme. 

---

# 10) Algorithm 1: Fairly Pricing Algorithm (FPA)

This is the full algorithm that combines everything.

## The pipeline

1. **Before epochs**: estimate (\hat F_{\min}).
2. **Initialize**: start with all procedurally fair policies.
3. **For each epoch**:

   * select good-and-exploratory policies,
   * estimate acceptance probabilities,
   * compute empirical optimal policy,
   * eliminate bad policies,
   * move to the next epoch. 

## Why it works conceptually

FPA is basically “learn a rough market model, optimize on the rough model, prune bad rules, repeat.” The doubling epochs make the estimates more reliable over time, and the conservative elimination keeps the true optimum alive. 

## What to implement

In code, this should be a loop over epochs with clear modules:

* warm-up estimation,
* exploratory policy selection,
* empirical estimation,
* LP-based optimization,
* elimination,
* update candidate set.

That is the clean implementation structure.

---

# 11) Regret Upper Bound: why the algorithm is good

Theorem 6 says FPA achieves

[
\tilde O(\sqrt{T})
]

regret and the same order of substantive unfairness, while keeping procedural unfairness exactly zero. 

## Why (\sqrt{T}) is the right scale

In online learning, (\sqrt{T}) cumulative regret means average regret per round goes to zero. So the algorithm becomes asymptotically optimal. 

## Proof idea in plain language

The proof is by induction over epochs:

1. assume the optimal fair policy is still in the candidate set,
2. show that the empirical estimates are accurate enough,
3. show the elimination step does not remove the true optimum,
4. bound the regret and unfairness incurred in that epoch,
5. sum over all epochs. 

## Why this induction works

Because the algorithm is conservative. It only removes policies when the estimates are strong enough to justify it. So the true fair optimum stays inside the feasible set as long as the concentration bounds hold. 

## Intuition

The algorithm is not trying to be brave; it is trying to be safe. Safety is what makes the proof go through.

---

# 12) Regret Lower Bound

Theorem 8 says any algorithm must suffer at least

[
\Omega(\sqrt{dT})
]

regret under the setting of the paper. 

## Why this is true

The proof reduces the fair two-group problem to ordinary online pricing by considering the special case where the two groups are identical:

[
F_1(i)=F_2(i)\quad \forall i.
]

Then fairness no longer helps you, and the problem becomes a standard finite-action pricing problem. That problem is known to have a bandit-style lower bound. 

## Intuition

If the groups are the same, then no clever fairness mechanism can make learning easier. So the fair problem cannot beat the classical pricing difficulty.

## Why (d) appears

More candidate prices mean more possible actions to distinguish, so learning becomes harder as the action space grows. That is why the lower bound scales with (\sqrt{dT}). 

---

# 13) Unfairness Lower Bound with Optimal Revenue

Theorem 9 says that any algorithm with near-optimal regret and zero procedural unfairness must still suffer (\Omega(\sqrt{T})) substantive unfairness. 

## What this means

Even if you are very good at maximizing revenue, you cannot always keep substantive unfairness tiny during learning.

## Proof idea

The paper constructs two very similar environments:

* one is the example from the introduction,
* the other is a slightly perturbed version. 

These environments are hard to distinguish from limited binary feedback. So any algorithm must either:

* commit too early and risk fairness errors, or
* wait too long and lose revenue. 

## Intuition

This is an information problem, not just an optimization problem.

If two worlds look almost the same from the data you observe, then learning which one you are in takes time. During that time, some substantive unfairness is unavoidable if you still want to keep regret low.

## Why this theorem matters

It tells you the fairness cost is not just because of a bad algorithm. It is a fundamental limitation of the online setting. 

---

# 14) What you need to preserve in implementation

If you want to implement this paper correctly, the following are non-negotiable:

* policy is a **pair of distributions**, not one fixed price,
* procedural fairness is enforced exactly,
* substantive fairness is enforced through the accepted-price ratio,
* acceptance probabilities must be estimated conservatively,
* policy elimination must be conservative,
* the optimization must be turned into LPs via the (w)-search,
* epochs must grow over time,
* candidate sets must shrink without removing the true optimum too early. 

---

# 15) Practical implementation roadmap

A clean implementation order would be:

1. implement the pricing environment simulator,
2. implement policy representation ((\pi^1,\pi^2)),
3. implement revenue, procedural unfairness, substantive unfairness,
4. implement warm-up and (\hat F_{\min}),
5. implement exploratory policy selection,
6. implement acceptance estimation,
7. implement Algorithm 2 LP search over (w),
8. implement Algorithm 3 elimination,
9. wrap them in the epoch loop,
10. test on the toy example from the paper before running larger experiments.

That is the safest path.

If you want, the next useful step is to turn this into a **code design blueprint**: data structures, function signatures, and a clean module-by-module implementation plan.
