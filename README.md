# Social-context-matters---Characterizing-adolescent-cooperation-strategies

This repository contains code and data for the paper:
**"Social-context-matters---Characterizing-adolescent-cooperation-strategies"**, by Wenda Liu, Isabella Mark, Christoph W. Korn, Gabriela Rosenblau, published in *Journal of Adolescense*
[DOI: 10.1002/jad.70025](https://doi.org/10.1002/jad.70025)

---

## Abstract

Introduction

Cooperation is a crucial prosocial skill refined during adolescence, a period marked by increased interactions with peers and artificial intelligence agents. While adolescents prioritize peer relationships, it remains unclear whether this translates into increased trust and reciprocity during social exchange. This study explores whether adolescents differentiate their cooperation strategies when interacting with human peers versus adaptive or fixed computer partners.

Methods

Adolescents (N = 67, 36 female, mean age: 12.3 ± 3.0 years) were recruited and invited into the laboratory to play a multiround trust game, first with a peer. In subsequent games, they were told that they continued to play with a peer (i.e., social condition) but they were playing with an adaptive algorithm. In the remaining games, they were correctly informed that they were playing with a computer (i.e., nonsocial condition).

Results

Adolescents sustained cooperation more when interacting with adaptive versus fixed computer partners. In the social condition, older adolescents and those with higher IQs and social skills were more prosocial than in the non-social condition. We fitted direct reciprocity strategies, reinforcement learning (RL) models, and a combination thereof to participants' choices. Direct reciprocity best captured adolescents' cooperation behavior. In the social condition, however, adolescents reinitiated cooperation more often after previously defecting. The winning model in the social condition was a more generous direct reciprocity strategy than the standard model. It used an RL forgiveness term that prescribed an evolving tendency to reinitiate cooperation.

Conclusion

Overall, we show that both partner adaptivity and social context play an important role in adolescents' cooperation decisions.

---

## Repository Structure
```bash
├── Data/                # processed datasets
├── Scripts/             # Scripts: data prep, analysis, modeling
    ├── CooperationBehavior/             # Scripts: data cleaning, analysis
    ├── MOdelling_analysis/              # Scripts: modeling
└── README.md
```

## Acknowledgments

This study was supported by the Bridge to Independence Award, awarded by the Simons Foundation for Autism Research. Christoph W. Korn was supported by the German Research Foundation (DFG); specifically by an Emmy Noether Research Group [392443797]. We thank Archana Venkataraman for her help with designing the adaptive computer algorithms and Kevin Pelphrey and Brent Vander Wyk for providing the research infrastructure to conduct this study at Yale. Gabriela Rosenblau was supported by the Hilibrand Postdoctoral Fellowship during the initial data collection at the Yale Child Study Center.

## Contact
[wenda_liu@gwu.edu](wenda_liu@gwu.edu) [grosenblau@gwu.edu](grosenblau@gwu.edu) 
