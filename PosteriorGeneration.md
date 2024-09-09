# Introduction

This is an introduction to the components and the generation of the
posterior distributions for the foodweb landscapes of Little Moose.
Below you will find snippits of code and examples of the data frames
that are used during this process. Here we will look at the data frames
that are associated with the `step1_LML_source_L.R` file. However, the
graphs will display the results of the posterior generation for both the
early and late periods of this experiment.

------------------------------------------------------------------------

## Define starting conditions for the model snippits.

For all periods of this experiment, we will be looking at 20 posterior
draws `n.posts`, 10 species `n_species`, 32 sites around the shore of
Little Moose Lake `n_sites`, and 4 years `n_years`. The early period
includes the years 2001-2004 and the late period includes the years
2019-2022.

    # Model conditions
    n.posts = 20
    n_species = 10
    n_sites = 32
    n_years = year_max - year_min -1 # seems missing 2023 data?

We begin by setting names for the data frame that holds the posterior
distributions. These names are composed of useful information that we
can later filter by. For example, `N_1_1_1` is the posterior
distribution for the estimated true abundance of species 1 at site 1
during the first year of that period.

    ##       N_1_1_1     N_2_1_1    N_3_1_1    N_4_1_1   N_5_1_1
    ## 1 0.035337345 0.108678941 0.05596628 0.11475847 0.3600164
    ## 2 0.008856450 0.046882053 0.08252995 0.03635013 0.2221019
    ## 3 0.038505263 0.042515693 0.15140276 0.21532147 0.1323097
    ## 4 0.041942243 0.006066467 0.08543462 0.20881001 0.3668654
    ## 5 0.005216185 0.079404302 0.01210736 0.01455679 0.1950053
    ## 6 0.072742179 0.025864737 0.20561266 0.17720656 0.5387817

But the entire posterior distribution is a large data frame. At this
point, we are not subsetting the posterior distribution and the
abundances are estimated at each site, at each year, for each species.
This results in a large number of columns.

    dim(chain_dat)

    ## [1] 4000 1030

## Subsetting the compontents of the landscape generation.

We begin with sigma. Sigma represents the shape of the isotope ellipse
for each species.It is composed of a covariance matrix.

Below we see an example of the sigma data from the posterior
distribution.

    ## # A tibble: 6 × 6
    ##    post species Sigma_1_1 Sigma_2_1 Sigma_1_2 Sigma_2_2
    ##   <int> <chr>       <dbl>     <dbl>     <dbl>     <dbl>
    ## 1     1 1           0.478   -0.417    -0.417      1.88 
    ## 2     1 2           1.10    -0.236    -0.236      0.206
    ## 3     1 3           0.147   -0.0368   -0.0368     1.45 
    ## 4     1 4           0.312    0.168     0.168      0.287
    ## 5     1 5           1.26    -0.462    -0.462      0.592
    ## 6     1 6           2.14     1.10      1.10       1.36

For `species == 1` and `post == 1`, see the below matrix as an example
of what the covariance matrix looks like.

    ##            [,1]       [,2]
    ## [1,]  0.4783136 -0.4167384
    ## [2,] -0.4167384  1.8750666

The mass data is generated through the use of length and weight
relationships for species in Little Moose Lake. The one exception here
is brown bullhead. For brown bullhead, we use weight-length
relationships that are available on FishBase from previously published
research.

    #mass weight graph

See below the structure of the mass data.

    ## # A tibble: 6 × 3
    ##    post species mass.avg
    ##   <int> <chr>      <dbl>
    ## 1     1 1        124.   
    ## 2     1 2        -15.0  
    ## 3     1 3        -65.4  
    ## 4     1 4        854.   
    ## 5     1 5         -0.846
    ## 6     1 6         11.0

Comparing the masses of fishes bewteen the two periods. Data from the
posterior for both the late and the early period is loaded in here.

![](PosteriorGeneration_files/figure-markdown_strict/unnamed-chunk-10-1.png)

See below for the snippit of the abundance data.

    ## `summarise()` has grouped output by 'post', 'species'. You can override using
    ## the `.groups` argument.

    ## # A tibble: 6 × 4
    ## # Groups:   post, species [2]
    ##    post species  year tot_abund
    ##   <dbl> <chr>   <dbl>     <dbl>
    ## 1     1 1           1      5.26
    ## 2     1 1           2      5.26
    ## 3     1 1           3      5.29
    ## 4     1 2           1    192.  
    ## 5     1 2           2    192.  
    ## 6     1 2           3    191.

See below for a comparison between the estimated total abundance for
both periods.

    ## Joining with `by = join_by(group)`

![](PosteriorGeneration_files/figure-markdown_strict/unnamed-chunk-12-1.png)

    ## # A tibble: 6 × 4
    ##    post species   mu_1   mu_2
    ##   <int> <chr>    <dbl>  <dbl>
    ## 1     1 1        2.28  -3.25 
    ## 2     1 2        0.636 -0.153
    ## 3     1 3       -0.459 -0.613
    ## 4     1 4       -0.573  0.605
    ## 5     1 5        0.313 -0.997
    ## 6     1 6        0.752  0.716

Finally, all the previous steps are combined. When applicable
(i.e. where abundances are estimated across multiple years) the
estimates are summarized by mean.

See below for the structure of the final posterior data frame is.

    ## # A tibble: 6 × 10
    ## # Groups:   post [1]
    ##    post species Sigma_1_1 Sigma_2_1 Sigma_1_2 Sigma_2_2 mass.avg tot_abund
    ##   <dbl> <chr>       <dbl>     <dbl>     <dbl>     <dbl>    <dbl>     <dbl>
    ## 1     1 1           0.478   -0.417    -0.417      1.88   124.         5.27
    ## 2     1 10          1.77     0.0294    0.0294     0.664  217.        92.3 
    ## 3     1 2           1.10    -0.236    -0.236      0.206  -15.0      191.  
    ## 4     1 3           0.147   -0.0368   -0.0368     1.45   -65.4      401.  
    ## 5     1 4           0.312    0.168     0.168      0.287  854.        17.6 
    ## 6     1 5           1.26    -0.462    -0.462      0.592   -0.846     26.4 
    ## # ℹ 2 more variables: mu_1 <dbl>, mu_2 <dbl>
