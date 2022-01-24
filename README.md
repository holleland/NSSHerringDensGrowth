Predicting density dependent somatic growth in Norwegian spring spawning
herring
================



*Contributors: Erling Kåre Stenevik<sup>1\*</sup>, Sondre
Hølleland<sup>1†</sup>, Katja Enberg<sup>2</sup>, Åge
Høines<sup>1</sup>, Are Salthaug<sup>1</sup>, Aril Slotte<sup>1</sup>,
Sindre Vatnehol<sup>1</sup>, Sondre Aanes<sup>3</sup>*

<sup>1</sup> Institute of Marine Research, Norway.<br> <sup>2</sup>
University of Bergen, Norway.<br> <sup>3</sup> Norwegian Computing
Center, Norway.<br> <sup>\*</sup> Corresponding author;
[erling.stenevik@hi.no](emailto:erling.stenevik@hi.no)<br> <sup>†</sup>
Responsible for the code. Correspondance related to this to:
[sondre.hoelleland@hi.no](emailto:sondre.hoelleland@hi.no)

[Paper reference goes here](https://academic.oup.com/icesjms)

In Stenevik et.al (2020a), we investigate how population density
influences individual growth of Norwegian spring spawning herring. A
parametrization of von Bertalanffy growth function (VBGF) is used to
model the growth in terms of length of the fish.

Since not all the data is available, rerunning everything is not
possible, but the main model and results can be reproduced by the user.
See details on how to get access to the necessary data below.

## Abstract

Density dependent growth, which might influence the effects of fisheries
on a population are often ignored when management strategies are
evaluated, mainly due to a lack of appropriate models readily available
to be implemented. To improve on this, we investigated if somatic growth
in Norwegian spring spawning herring (Clupea harengus) depend on cohort
density using a formulation of the von Bertalanffy growth function on
cohorts from 1921 to 2014 and found a significant negative correlation
between estimated asymptotic length and density. This clearly indicates
density dependent effects on growth, and we propose a model which can be
used to predict size-at-age of Norwegian spring spawning herring as
function of herring density (the abundance of two successive cohorts) in
future estimation of reference points (FMSY) and short-term predictions
of catch advice.

## Data

The main data is available at [Stenevik et al
(2022b)](https://doi.org/10.21335/NMDC-496562593). We do not have
permission to publish the temperature data, and these are therefore not
publicly available.

In our *R/1\_data.R* script, the individual herring data is downloaded
by running the following code:

``` r
if(!("HerringData.csv" %in% list.files(path = "inputdata/") )) {
  download.file(url = "https://ftp.nmdc.no/nmdc/IMR/Herring/HerringData.csv", 
                destfile = "inputdata/HerringData.csv")
}
```

We are not at liberty to publish the XSAM series here on github, but the
user can download it from the paper supplementary material (Table S3).
If you save it as *inputdata/N.txt*, the *R/1\_data.R* script will run
as intended without adjustments to the code.

## Authors’ github accounts

**Sondre Hølleland** - [holleland](https://github.com/holleland)

**Sindre Vatnehol** -
[sindrevatnehol](https://github.com/sindrevatnehol)

## R version

The code has been run on the following R version.

``` r
R.version
```

    ##                _                           
    ## platform       x86_64-w64-mingw32          
    ## arch           x86_64                      
    ## os             mingw32                     
    ## system         x86_64, mingw32             
    ## status                                     
    ## major          4                           
    ## minor          1.1                         
    ## year           2021                        
    ## month          08                          
    ## day            10                          
    ## svn rev        80725                       
    ## language       R                           
    ## version.string R version 4.1.1 (2021-08-10)
    ## nickname       Kick Things

## References

Erling Kåre Stenevik (HI), Sondre Hølleland (HI), Katja Enberg (UiB),
Åge Høines (HI), Are Salthaug (HI), Aril Slotte (HI), Sindre Vathehol
(HI), Sondre Aanes (NR) (2022a) Predicting density dependent somatic
growth in Norwegian spring spawning herring. Under review in *ICES
Journal of Marine Science*.

Erling Kåre Stenevik (HI), Sondre Hølleland (HI), Katja Enberg (UiB),
Åge Høines (HI), Are Salthaug (HI), Aril Slotte (HI), Sindre Vathehol
(HI), Sondre Aanes (NR) (2022b) Individual samples of Norwegian Spring
Spawning herring 1935-2019 <https://doi.org/10.21335/NMDC-496562593>

## License

The code for this project is licensed under the [MIT
licence](LICENCE.md).

[<img src="https://www.hi.no/en/hi/resources/layout/HI-logo-farger-engelsk.svg/original"
alt="Institute of Marine Research" width="200"/>](https://www.hi.no/en)
