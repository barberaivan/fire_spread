# Fire spread modelling

## Overview

The aim of this project is to fit a spatially-explicit fire spread model able to represent general regional-scale patterns of fire activity in northwestern Patagonia. In short, we aim to improve the model fitted by [Morales et al. (2015)](https://ri.conicet.gov.ar/bitstream/handle/11336/38304/CONICET_Digital_Nro.d2f95f9f-ea7f-49ea-8ac4-593883434965_A.pdf?sequence=2) by including more data and increasing the model's complexity.

Fire spread is a well-studied process and there are many models able to produce quite reallistic fire polygons and predict fire behaviour metrics. In general fire spread models can be grouped in two main and difuse categories: empirical, simple and coarse-resolution vs. mechanistic, complex, and fine-resolution. Mechanistic ones (e.g., FARSITE) are generally better at predicting the advance of particular fires based on fine-scale meteorology, fuels condition and topography. These model types are usually transferable to many ecosystems if the data to parameterize them is available, and their parameters, which are many, are usually based on laboratory tests and fire physics knowledge. On the other side, mechanistic models are usually simpler and make more general fire spread predictions, running at coarser temporal scales and requiring more simple data. Empirical models are frequently fitted from data, as the lower number of parameters makes it possible.

We aim to include the fire spread model in a dynamic vegetation model that simulates the trayectory of northwestern patagonian forests and is used to better understand forest dynamics in the context of global change. This landscape model (ALLADYNS (Kitzberger et al., in preparation)) runs at annual time steps, and fires are simulated in an atemporal scheme within each fire season. In this context, we need a model able to represent how regional patterns of fire activity are going to respond to climate change in the 21st century. As parameterizing mechanistic models requires data that we do not have and would imply a higher computation time, we prefer to fit a simpler one.

## Model structure

We defined the fire spread model representing the landscape with a lattice, as a cellular automata. Each bunable cell (i.e., not water, rock, or bare ground) can be not burned, burning of burned in every time step, and the burning state only lasts one step. The fire starts at an ignition point (cell) and spreads towards its 8-pixel neighbourhoood. The spread from fire i  to cell j ($y_{ij} \in \{0, 1\}$) follows a Bernoulli distribution with probability $p_{ij}$, defined as a linear function of the landscape covariates ($\boldsymbol{\text{X}}_{ij}$) and with an upper limit $u \in [0, 1]$:

\begin{aligned}
y_{ij} &\sim \text{Bernoulli}(p_{ij})\\
p_{ij} &= \frac{u}{1 + \exp(-\eta_{ij})} \\
\eta_{ij} &= \alpha \ + \boldsymbol{\text{X}}_{ij} \boldsymbol{\beta }\\
\end{aligned}

We fix $u = 0.5$ to avoid fires spreading through the whole landscape, and the estimation is aimed to find an approximate posterior distribution for $\boldsymbol{\beta}$. $\boldsymbol{\text{X}}_{ij}$ is a row vector from the design matrix, and its elements are:  

\begin{aligned}
&\text{subalpine}_j \in \{0, 1\} \\
&\text{wet}_j \in \{0, 1\} \\
&\text{dry}_j \in \{0, 1\} \\
&\text{FWI}_j \in (-\infty, \infty) \\
&\text{northing}_j \in [-1, 1] \ \ \text{called aspect in the code}\\
&\text{wind}_{ij} \in [-1, 1] \\
&\text{elevation}_{j} \in  (-\infty, \infty) \\
&\text{slope}_{ij} \in [-1, 1] \\
\end{aligned}

subalpine, wet and dry variables are dummies indicating the forest type. If all of them are zero, then the jth cell is shrubland, the reference category. FWI is the pixel-level anomaly in the summer-average daily values of the Fire Weather Index, intended to represent climatic interannual variability. Northing is computed as $\sin(\text{aspect}_j)$, taking 1 if the cell looks to the north and -1 if the cell looks to the south. wind takes 1 if the wind direction is the same as the direction from cell i to j, and -1 if it is exactly the opposite (**INSERT EQUATION HERE**). Elevation is the standardized elevation of pixel j, and slope is computed from the elevation of i and j and considering the distance between cells, taking 1 if the target cell is exactly above burning cell (impossible condition) and -1 if it is exactly below (**INSERT EQUATION HERE**). The first four elements in $\boldsymbol{\beta}$ lay in $(-\infty, \infty)$, but the remaining ones are constrained to reasonable regions of the parameter space to decrease the computation time and avoid non-sense estimations:

\begin{aligned}
&\beta_{\text{FWI}} \in (0, \infty) \\
&\beta_{\text{northing}} \in (0, \infty) \\
&\beta_{\text{wind}} \in (0, \infty) \\
&\beta_{\text{elevation}} \in  (-\infty, 0) \\
&\beta_{\text{slope}} \in (0, \infty) \\
\end{aligned}

In principle, one could choose the parameters values by eye (yes, by eye) to make a reasonable simulator, but serious people prefer to estimate them from data. The problem here is that the data about how real fires have spread pixel by pixel is not available (and, if possible, it would be really hard to get). So we cannot evaluate the likelihood of each fire transition. This takes us to likelihood-free methods. In our case, we know the likelihood function (derivable from eq. 1), but we do not have the data to compute it. An appealing solution is to use Approximate Bayesian Computation. This method replaces the likelihood function in the Bayes' theorem by a similarity measure between summaries of observed and simulated data. In our case, the summary of the (unobserved) spread data would be the final fire polygon and/or the number of pixels burned by vegetation type. Using this information we can compute a similarity measure between observed and simulated fire polygons $\rho(F_{\text{obs}}, F_{\text{sim}})$. As it is hard to know in advance which similarity function more suitable to emulate the likelihood, we test a few. The simplest is the spatial overlap, taking 1 if fires are exactly the same and 0 if they are disjoint (it cannot happen because they share at least the ignition point):

spatial overlap 
\rho_1(F_{\text{obs}}, F_{\text{sim}}) = \frac{F_{\text{obs}} \cap F_{\text{sim}}}{F_{\text{obs}} \cup F_{\text{sim}}}

where the numerator is the number of pixels of the intersection between both fire polygons, and the denominator is the number of pixels of their union.

Then this metric is combined with others that evaluate in a non-spatial manner the similarity in burned area by vegetation type:




  
    normalized difference
    [EQ]  
    squared normalized difference
    [EQ]  
    expquad normalized difference
    [EQ]  

These three non-spatial metrics are averaged with the spatial overlap with weights 0.25 and 0.5 to conform the following similarity measures:

    chequear orden en cpp

\begin{aligned}
\rho_2(F_{\text{obs}}, F_{\text{sim}}) &= 0.75 \rho_1 + 0.25 nd \\
\rho_3(F_{\text{obs}}, F_{\text{sim}}) &= 0.75 \rho_1 + 0.25 snd \\
\rho_4(F_{\text{obs}}, F_{\text{sim}}) &= 0.75 \rho_1 + 0.25 end \\
\\
\rho_5(F_{\text{obs}}, F_{\text{sim}}) &= 0.5 \rho_1 + 0.5 nd \\
\rho_6(F_{\text{obs}}, F_{\text{sim}}) &= 0.5 \rho_1 + 0.5 snd \\
\rho_7(F_{\text{obs}}, F_{\text{sim}}) &= 0.5 \rho_1 + 0.5 end \\
\end{aligned}

[to be continued...]

As simulating fire spread is computationally expensive we write the simulation function in ```C++``` and use the ```Rcpp``` package to make it available in ```R```.

We use large files storing the landscape data, so we don't upload them to the repo. The ```data``` folder is available on request to ivanbarbera93@gmail.com, and should be placed inside the cloned repo to work (e.g.: ```/path_to_the_repo/fire_spread/data```).
