#abstract writing exercise

## Original abstract
Wind pollinated conifer species rely on synchronized pollen shed and cone receptivity for mating. Fecundity, gene flow, and range limits of many widespread North American conifers are likely to be negatively affected as the climate changes and their pollination phenology shifts. In this study, we model the relationship between temperature and pollination phenology in *Pinus contorta* var. *latifolia* (lodgepole pine) and calculate shifts in pollination phenology and length of phenological period under different warming scenarios. Using pollination phenology data from trees grown at 8 sites over 15 years and sourced from across British Columbia, we built a hierarchical Bayesian model predicting the timing and length of pollen shed and cone receptivity and extended it across entire range of lodgepole pine in the present and under three climate change scenarios. We find that there is no significant genetic effect on pollination phenology and only heat sum is required for prediction. This leads to strong spatial patterns in pollination phenology that may create and maintain local adaptation. As temperatures warm, the pollination phenology period advances and shortens, desynchronizing populations and restricting gene flow. Pollen shed advances more rapidly than receptivity, especially in southern populations, which may lead some populations to become pollen limited and increases outcrossing. Pollination phenology is driven by temperature, and the phenological shifts expected under climate change are likely to restrict gene flow in the future.

## Notes on abstract from Graydon
    Wind pollinated conifer species rely on synchronized pollen shed and cone receptivity for mating.


I think this sentence should introduce -- or be followed by a short
sentence that introduces -- the term "pollination phenology", since it
is used in the following sentences as though it is defined. Since
you're starting from a kinda "first principles" framing here, it'd be
good to reiterate the terms.

    Fecundity, gene flow, and range limits of many widespread North American conifers are likely to be negatively affected as the climate changes and their pollination phenology shifts.


Two issues: clause ordering should reflect causal ordering ("As the
climate changes ... fecundity etc. are likely to be effected.") and
also this sentence ends with a qualifier "as the climate changes and
their pollination phenology shifts" that actually presents as an
assumption the specific question you're investigating: whether or not
climate change will cause the pollination phenology to shift.
Suggested rephrasing that brings that question to the forefront:
"Changes in temperature due to climate change may cause changes to
pollination phenology of many widespread North American conifers,
which could subsequently affect their fecundity, gene flow and range
limits."

    In this study, we model the relationship between temperature and pollination phenology in *Pinus contorta* var. *latifolia* (lodgepole pine) and calculate shifts in pollination phenology and length of phenological period under different warming scenarios. Using pollination phenology data from trees grown at 8 sites over 15 years and sourced from across British Columbia, we built


"built" is past tense, previous sentence was present tense.

    a hierarchical Bayesian model predicting the timing and length of pollen shed and cone receptivity and extended it across entire range of lodgepole pine in the present and under three climate change scenarios. We find that there is no significant genetic effect on pollination phenology and only heat sum is required for prediction.


FIrst sentence is overlong. End at "receptivity" and make a new
sentence for the extended part and possibly another for the 3 climate
change scenarios. Further: I don't know from this sentence and the
following one (which offers a conclusion drawn from the model) where
the modelling of observed-data stops and the predicting of
non-observed data starts. I _think_ you're saying that you built a
model from 8 sites / 15 years, and then concluded from that model that
heat sum was the only requirement for prediction, and then applied a
model-with-only-heat-sum to the entire range of the species in the
present, and then further applied it to 3 climate change scenarios in
the future. If that's the order, then this should be rewritten as
such:

"... a hierarchical Bayesian model predicting the timing and length of
pollen shed and cone receptivity. We find that there is no significant
genetic effect on pollination phenology and only heat sum is required
for prediction. We then use the model to predict pollination phenology
across the entire range of lodgepole pine in the present. Finally, we
extend this prediction over time into three possible climate change
scenarios."

Also curious: you did not mention genetic effects until the findings
sentence. Maybe describe the bayesian model as "a hierarchical
bayesian model including both temperature and genetic factors" (and
any other factors you included-and-disproved).

    This leads to strong spatial patterns in pollination phenology that may create and maintain local adaptation.


"This leads to" => "Our model predictions show".

    As temperatures warm, the pollination phenology period advances and shortens, desynchronizing populations and restricting gene flow. Pollen shed advances more rapidly than receptivity, especially in southern populations, which may lead some populations to become pollen limited and increases outcrossing.


Tense/sense mixup: "and increases outcrossing" vs. "and to increase
outcrossing" (I think? reads odd)

    Pollination phenology is driven by temperature, and the phenological shifts expected under climate change are likely to restrict gene flow in the future.


"We conclude that pollination phenology is driven by ..."

Mostly reads good! Just a little tinkering.

#Abstract attempt 2

The northern hemisphere is home to vast coniferous forests, and the ranges of many species in these forests encompass remarkably diverse climates.
Despite high levels of gene flow due to wind-dispersed pollen, adaptation to local climates is common. Phenological traits like budburst and budset have been studied extensively due to their influence on growth and survival. However, pollination phenology is less studied in conifers, despite its influence on fecundity, adaptation, and range limits, and its potential role in the creation and maintainance of local adaptation.

*Pinus contorta* spp. *latifolia* (lodgepole pine) is a widespread, wind-pollinated conifer species with local adaptation across its range. In this study, we model the relationship between temperature and pollination phenology in lodgepole pine, then calculate shifts in pollination phenology and length of phenological period under different warming scenarios. 

Using pollination phenology data from trees grown at 8 sites over 15 years and sourced from across British Columbia with climate data from local weather stations, we construct a hierarchical Bayesian model predicting the timing and length of pollen shed and cone receptivity period and accounting for population source. We then extend the predictions across the entire range of lodgepole pine in the present and under three climate change scenarios.

While fall phenological traits like budset in lodgepole pine show strong signals of local adaptation, bud burst in the spring does not. In keeping with these results, we find that there is no significant genetic effect on pollination phenology and only heat sum is required for prediction. 
Applying our model across the species range shows strong spatial patterns in pollination phenology that may create and maintain local adaptation. A
s temperatures warm, the pollination phenology period advances and shortens, desynchronizing populations and restricting gene flow. Pollen shed advances more rapidly than receptivity, especially in southern populations, which may lead some populations to become pollen limited and increase outcrossing.

In conclusion, lodgepole pine pollination phenology is driven by temperature and, as climate changes, the timing, length, and synchronization of pollen shed and cone receptivity is likely to shift with implications for evolution and conservation.