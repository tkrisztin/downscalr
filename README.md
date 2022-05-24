# downscalr
An R package for downscaling of land-use and land-use change projections. You find the same information (and example [vignettes](https://tkrisztin.github.io/downscalr/articles/downscalr_tutorial.html)) on the [pkgdown](https://pkgdown.r-lib.org/) website [here](https://tkrisztin.github.io/downscalr/).

This package was developed in the context of the SPIPA Argentina project, in joint cooperation with colleagues from the [National Agricultural Technology Institute](https://www.argentina.gob.ar/inta). The main purpose of this was to port IIASA's [DownScale](https://github.com/iiasa/DownScale) model, for easier comparison with INTA's Argentina specific implementation of [Dinamica EGO](https://csr.ufmg.br/dinamica/) for the purpose of downscaling [FABLE Calculator](https://www.abstract-landscapes.com/fable-calculator) output. This enables a more in-depth representation of high-resolution land-use change projections.

Originally downscalr was developed to provide high-resolution projections of the [GLOBIOM](https://iiasa.ac.at/web/home/research/GLOBIOM/GLOBIOM.html) model. However, within the course of SPIPA Argentina we have expanded it to be fully compatible with FABLE Calculator output. Moreover, the code was ported to R and documented thoroughly, to enable ease of use by INTA, our Argentinian project partner.

# Installation:

You need R to run the scripts. In R the following commands install the `devtools` and `downscalr` packages. Devtools is required to install packages from Github. 
For now, if you want to install the package with its vignette -- this provides a quick guide on how to use `downscalr` --, `knitr` and `rmarkdown` are required as well. 

      # just install the downscalr package
      ## install.packages("devtools")
      devtools::install_github("tkrisztin/downscalr", ref="HEAD", repos = "http://cran.us.r-project.org")
      
      # install the downscalr package with vignette 
      # (note: this may take a couple minutes since some computing is required to render the vignette). 
      ## install.packages(c("devtools", "knitr", "rmarkdown"))
      devtools::install_github("tkrisztin/downscalr", ref="HEAD", build_vignettes = TRUE, repos = "http://cran.us.r-project.org")
      
# Running downscalr

If you have have followed the install instructions with vignette above, you can view a step by step tutorial via the following command in R: 

      browseVignettes("downscalr")

This should open a tab in your browser, click on "HTML" to access the tutorial. For additional guidance see the DownScale package documentation here: https://bit.ly/3fiLG3u
      
# Additional information
Supported by:

<p float="left">
  <img src="https://user-images.githubusercontent.com/9318979/152320807-db56e9d0-b0bb-4a5a-8f5a-51d1ab679de6.png" height="100" />
  <img src="https://user-images.githubusercontent.com/9318979/152320782-f64da019-96dd-4e3c-beb2-97989d01602d.png" height="100" /> 
</p>

*Downscalr was developed with the financial support of the European Union’s Partnership Instrument and the German Federal Ministry for the Environment, Nature Conservation, and Nuclear Safety (BMU) in the context of the International Climate Initiative (IKI). Its contents are the sole responsibility of International Institute for Applied Systems Analysis (IIASA) and do not necessarily reflect the views of the funders.*

