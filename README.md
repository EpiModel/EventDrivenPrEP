# Modeling the Potential Impact of Scaling Up Event-Driven PrEP Among Gay, Bisexual, and Other Men who have Sex with Men


This repository contains the analysis scripts and tables for the following study:

> Chandra, C., Maloney, K.M., Le Guillou, A., Hoover, K.W., Jenness, S.M. (2025). Modeling the Potential Impact of Scaling Up Event-Driven PrEP Among Gay, Bisexual, and Other Men who have Sex with Men. Under review.

<img src="https://github.com/EpiModel/Mean-Degree-Analysis/blob/master/Figures/Figure1.png](https://github.com/EpiModel/EventDrivenPrEP/blob/main/figures/Fig1.png">

## Abstract

### Background

Event-driven preexposure prophylaxis/PrEP (EDP), or PrEP taken around the time of sex, could be a preferred option over daily PrEP among gay, bisexual, and other men who have sex with men (MSM). However, recommendations for EDP versus daily dosing require consideration of optimal strategies for HIV prevention based on heterogenous behavioral and biological PrEP indications.

### Methods

We used network-based mathematical HIV transmission modeling to understand the potential population-level of impact of implementing EDP alongside daily PrEP using the 2-1-1 dosing schedule: 2 pills on the day of sex, 1 pill a day after, and 1 pill 2 days after. We simulated HIV transmission among 100,000 MSM over 10 years. Compared to a baseline scenario with only daily PrEP and no EDP, we tested counterfactual scenarios that increased EDP use by changing initiation rates, expanding EDP to MSM not indicated for daily PrEP, and targeting EDP among MSM that discontinued PrEP. 

### Results

Increasing EDP initiation such that overall PrEP use was nearly double that at baseline (28.5% vs. 16.1%) yielded a percent of infections averted (PIA) of 18.3%. Expanding eligibility for EDP to any sexually active MSM resulted in a PIA of 15.8%. Targeting EDP for MSM discontinuing PrEP did not yield population-level improvements. Increasing initiation of EDP had a greater impact on PIA compared to improving EDP adherence.

### Conclusions

Using EDP to improve overall PrEP use resulted in reductions in HIV incidence, while targeted EDP initiation for MSM with low daily PrEP adherence did not improve outcomes over 10 years.


## Data

We used data from ARTnet - an anonymous cross-sectional online survey of HIV-related risk behaviors, testing, and use of prevention services among MSM in the United States. MSM were recruited from the American Mens’ Internet Survey (AMIS) Survey, so the dataset also includes variables from AMIS.

Additional documentation on ARTnet and information to accessing the data can be found [here](https://github.com/EpiModel/ARTnetData). Code to install the “ARTnetData” package can be found below, but it may require a [Github Personal Access Token](https://help.github.com/en/articles/creating-a-personal-access-token-for-the-command-line) since it is a private repository.

```r
install.packages("remotes")
remotes::install_github("EpiModel/ARTnetData")
```

## Code Organization

These models are written and executed in the R statistical software language. To run these files, it is necessary to first install our epidemic modeling software, [EpiModel](https://github.com/EpiModel/EpiModel/), and our extension package specifically for modeling HIV transmission dynamics, [EpiModelHIV](https://github.com/EpiModel/EpiModelHIV). The branch of the EpiModelHIV repository associated with this research project is EventDrivenPrEP.

The installation is accomplished with the renv package in R. First install `renv` (if you do not already have it installed) and run:

```r
renv::init()
```
in your project directory. Select the option to restore the package set from the `renv.lock` file.
