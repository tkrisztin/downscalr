# downscalr
A package for high-resolution simulation of land-use change projections, with support for downscaling. 

This package was developed in the context of the SPIPA Argentina project, in joint cooperation with colleagues from the [National Agricultural Technology Institute](https://www.argentina.gob.ar/inta). The main purpose of this package is to combine the IIASA's [DownScale](https://github.com/iiasa/DownScale) model, with INTA's Argentina specific implementation of the [Dinamica EGO](https://csr.ufmg.br/dinamica/) for the purpose of downscaling [FABLE Calculator](https://www.abstract-landscapes.com/fable-calculator) output. This will enable a more in-depth representation of high-resolution land-use change projections.

Originally DownScale was developed to provide high-resolution projections of the [GLOBIOM](https://iiasa.ac.at/web/home/research/GLOBIOM/GLOBIOM.html) model. However, within the course of this project we have expanded it to be fully compatible with FABLE Calculator data. Moreover the code was ported to R and documented thoroughly, to enable ease of use by Argentinian project partner.

# Installation:

You need R to run the scripts. In R the following commands install the devtools and downscalr packages. Devtools is required to install packages from Github.

      # install.packages("devtools")
      devtools::install_github("tkrisztin/downscalr",ref="HEAD")
      
# Running downscalr

First a prior model needs to be run. This can be done using the mnlogit() function. For further help type in the R-console:

      help(downscalr)

For additional guidance see the DownScale package documentation here: https://bit.ly/3fiLG3u
      
# Additional informatio 

*Downscalr was developed  with the financial support of the European Union’s Partnership Instrument and the German Federal Ministry for the Environment, Nature Conservation, and Nuclear Safety (BMU) in the context of the International Climate Initiative (IKI), particularly the SPIPA Argentina project.*

![SPIPA+wordmark](https://user-images.githubusercontent.com/9318979/150418047-76d8477c-20be-459a-8e3d-bee6c6bb5032.jpeg)

*The programme Strategic Partnerships for the Implementation of the Paris Agreement (SPIPA) is jointly commissioned by the European Commission as a Foreign Policy Instrument Action and the German Federal Ministry for the Environment, Nature Conservation and Nuclear Safety in the context of the International Climate Initiative (IKI). SPIPA is implemented by the Deutsche Gesellschaft für Internationale Zusammenarbeit (GIZ) GmbH.*

