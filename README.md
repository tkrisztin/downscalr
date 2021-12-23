# downscalr
A package for high-resolution simulation of land-use change projections, with support for downscaling. Based on IIASA's [DownScale](https://github.com/iiasa/DownScale) package, which was developed to provide high-resolution projections of the [GLOBIOM](https://iiasa.ac.at/web/home/research/GLOBIOM/GLOBIOM.html) model.

# Installation:

You need R to run the scripts. In R the following commands install the devtools and downscalr packages. Devtools is required to install packages from Github.

      # install.packages("devtools")
      devtools::install_github("tkrisztin/downscalr",ref="HEAD")
      
# Running downscalr

First a prior model needs to be run. This can be done using the mnlogit() function. For further help type in the R-console:

      help(downscalr)

For additional guidance see the DownScale package documentation here: https://bit.ly/3fiLG3u
      
# Additional informatio 

Downscalr was developed  with the financial support of the European Union’s Partnership Instrument and the German Federal Ministry for the Environment, Nature Conservation, and Nuclear Safety (BMU) in the context of the International Climate Initiative (IKI), particularly the SPIPA Argentina project.

