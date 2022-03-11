# downscalr
A package for high-resolution simulation of land-use change projections, with support for downscaling. 

This package was developed in the context of the SPIPA Argentina project, in joint cooperation with colleagues from the [National Agricultural Technology Institute](https://www.argentina.gob.ar/inta). The main purpose of this package is to combine the IIASA's [DownScale](https://github.com/iiasa/DownScale) model, with INTA's Argentina specific implementation of the [Dinamica EGO](https://csr.ufmg.br/dinamica/) for the purpose of downscaling [FABLE Calculator](https://www.abstract-landscapes.com/fable-calculator) output. This will enable a more in-depth representation of high-resolution land-use change projections.

Originally DownScale was developed to provide high-resolution projections of the [GLOBIOM](https://iiasa.ac.at/web/home/research/GLOBIOM/GLOBIOM.html) model. However, within the course of this project we have expanded it to be fully compatible with FABLE Calculator data. Moreover the code was ported to R and documented thoroughly, to enable ease of use by INTA, our Argentinian project partner.

# Installation:

You need R to run the scripts. In R the following commands install the devtools and downscalr packages. Devtools is required to install packages from Github.

      # install.packages("devtools")
      devtools::install_github("tkrisztin/downscalr", ref="HEAD", build_vignettes = TRUE)
      
# Running downscalr

For now you can view a step by step tutorial on how to run "downscalr' -- if you have have followed the install instructions above -- with the following command in R: 

      browseVignettes("downscalr")

This should open a tab in your browser, click on "HTML" to access the tutorial. For additional guidance see the DownScale package documentation here: https://bit.ly/3fiLG3u
      
# Additional information
Supported by:

<p float="left">
  <img src="https://user-images.githubusercontent.com/9318979/152320807-db56e9d0-b0bb-4a5a-8f5a-51d1ab679de6.png" height="100" />
  <img src="https://user-images.githubusercontent.com/9318979/152320782-f64da019-96dd-4e3c-beb2-97989d01602d.png" height="100" /> 
</p>

*Downscalr was developed with the financial support of the European Union’s Partnership Instrument and the German Federal Ministry for the Environment, Nature Conservation, and Nuclear Safety (BMU) in the context of the International Climate Initiative (IKI). Its contents are the sole responsibility of International Institute for Applied Systems Analysis (IIASA) and do not necessarily reflect the views of the funders.*

