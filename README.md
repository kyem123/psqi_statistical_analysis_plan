# PSQI_Statistical_Analysis_Plan
Authors: K.M. and R.V.R.
Date: 28/05/2026

The purpose of this script is to address the secondary outcomes of subjective sleep quality data, as measured by the Pittsburgh Sleep Quality Index (PSQI).  

The study is a double-blinded placebo-controlled crossover RCT, in which participants were allocated to one of two groups. Each group received a placebo or active treatment for six weeks, before undergoing a four-week washout period and receiving the other treatment for a further six weeks. Wristwatch actigraphy measures were collected for the seven days preceding each trial arm and for the last seven days of each arm. PSQI data were also collected pre and post trial arms.

Secondary outcomes include subjective sleep quality as assessed using the PSQI. Specifically, this includes the Global PSQI Score (0-21) and seven individual component scores (0-3).

The primary outcomes were objective actigraphy-derived sleep parameters, specifically sleep period time (SPT), sleep duration, sleep efficiency and wake after sleep onset (WASO). These measures capture complementary aspects of sleep quantity and quality and are therefore analysed concurrently. Primary outcomes are analysed in a separate script.

**Statistical Analysis**
The primary objective of the analysis is to assess whether pre-to-post changes in global and component PSQI scores differ between treatment conditions, expressed as a difference-in-differences (DID) effect. 

The dataset contains repeated observations for each participant across:
-	Treatment (within-subject): A vs B
-	Time (within-subject): pre vs post
-	Treatment Period (within-subject): period 1 vs period 2
-	Treatment Order (between subject): A-first vs B-first

The global PSQI score is analysed as a continuous outcome, while individual component scores are analysed as ordinal outcomes. Component 6 (sleep medication use) is recorded as a binary outcome due to a bimodal distribution and is analysed separately using a logistic mixed model. 

For the global PSQI score, linear mixed effects models with a random intercept (1|id) are fitted. Fixed effects structures, including treatment x time with optional adjustments for order and period are compared using AICc, with the best-fitting model selected. The models will then be tested for non-normality, heteroscedasticity and collinearity. If the model assumptions are violated, transformations are applied sequentially (square-root then log), followed by a gamma generalised linear mixed model. Assumptions are also evaluated manually through visual inspection, as formal statistical tests can be overly sensitive with small sample sizes and sleep data. There is an option to override the automatic model selection based on visual diagnostics.

For ordinal component scores, cumulative link mixed models (CLMMs) are fitted, with fixed-effects structures selected via AICc. For component 6, a binary logistic mixed model (glmer) is fitted, once again with model selection based on AICc.
