*EJ Replicate Results*

Here I put the replicate results.

**Black 2021**
This exercise compare the first and second child born in a three plus children family between 1980 and 2015 and aged between 25 to 65. All variables, except for education (which is observed cross-sectionally in 2021), are constructed as panel data covering the period from 2000 to 2020.

The list of disabilities is the one in Black et al. (2021). 

Note that in this exercise, disabilities are identified based on outpatient records from 2000 to 2020, which might pose some limitations on the estimation if using similar treatment designs in the paper. See more details below.

In the final constructed panel, there are approximately 610,000 families. Of these, around 270,000 have a third child identified with a "disability" according to the applied definition. 

This prevalence are unusually high, and I think that the disability criteria is looser than the definitions commonly used in Taiwan. For additional context, approximately 60,000 third children are diagnosed with a disability by age 5, and 10,000 by age 10.

***Model***

There are four treatments (interactions) and two regression designs, so eight models in total for each outcome. 

The first treatmen t(tb_disab) is whether the third child is disabled. The other three are defined as in the EJ paper.

The first regression design only includes family fe (and year fe if any); the second design is the standard form in the paper which controls for gender and includes family fe plus birth year and birth month fe (and year fe if any). 

Interactions: 
- int: tb_disab x second
- int_by_5: tb_disab_by_5 x second
- int_in_5_10: tb_disab_in_5_10 x second
- int_by_10: tb_disab_by_10 x second

***Outcomes***

Labor
- work: whether she worked in that year
- income: whether she earned labor income in that year
- ft_work: whether she worked full-time every month in that year
- ft_income: whether she earned labor income from a full-time job every month in that year

Edu
- highest_edu: her highest education level in year term in 2021

Marr
- marr: whether she was married in that year

Fert
- fertility: number of new born child in that year

Health
- total_dot: a measure of health expenditure in that year
- outpatient_count: the number of outpatient visits in that year
- chronic: whether she was diagnosed as any chronic disease in that year
- major: whether she was diagnosed as any major disease in that year
- md:  whether she was diagnosed as any mental illness in that year


Additionally, two alternative definitions of disability are available:

**Catastrophic Illness Registry**

This dataset includes individuals diagnosed with severe diseases or disabilities that require long-term or life-long treatment. It enables the identification of disabilities from 1995 onwards.

**Disability Certification Registry**

This dataset includes individuals holding an officially recognized disability certificate as defined under Taiwan’s legal framework. 
Disability information is available from 1985 onward, though I'd say the reliability and completeness of the records improve significantly after 1990.