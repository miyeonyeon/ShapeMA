# ShapeMA

Shape mediation analysis designed to explore the causal relationships between genetic exposures and clinical outcomes, whether mediated or unmediated by shape-related factors while accounting for potential confounding variables. Within our framework, we apply the square-root velocity function to extract elastic shape representations, which reside within the linear Hilbert space of square-integrable functions. Subsequently, we introduce a two-layer shape regression model to characterize the relationships among neurocognitive outcomes, elastic shape mediators, genetic exposures, and clinical confounders.


## R Functions

To estimate coefficients of a shape-on-scalar regression model and a scalar-on-shape partial linear single-index model, we employ a multivariate varying coefficient model (Zhu et al., 2012) and a partially linear single-index model (Liang et al., 2010). The resampling methods (Loh et al., 2022) are applied to calculate causal estimands. Our example codes can be used for a simplified version of 2SRM, where the number of features is 1.

- `SMA_MVCMprep`: prepare multivariate varying coefficient model formats to read a dataset.
- `SMA_PLSIMread`: transform a scalar-on-shape partial linear single-index model.
- `SMA_PLSIMest`: get profile least squares estimator in a scalar-on-shape partial linear single-index model.
- `SMA_CEdupdata`: create duplicated data for one subject.
- `SMA_CEest`: estimate the direct and indirect effects using an unpenalized outcome model.
- `SMA_SAIEest`: estimate the spatial average indirect effects.

Converted MATLAB to R to apply estimation procedure for a multivariate varying coefficient model.

- `MVCM_read`: read raw data and generate, arc length, standardized design and response matrices, and related dimension parameters.
- `MVCM_lpks_wob`: read arc length, design and response matrices and generate optimal bandwidth for weighted least squares estimation.
- `MVCM_lpks_wb1`: read arc length, design and response matrices, and optimal bandwidth and generate the estimated coefficient functions, their first derivatives, and fitted responses using weighted least squares estimation.


<!--- ## 2-layer Shape Regression Model --->

<!--- \usepackage{amsmath} --->
<!--- $\psi_{i}(s_k) = x_i\alpha(s_k)+\boldsymbol{w}_i^\top\boldsymbol{\lambda}(s_k)+\eta_{i}(s_k)+\epsilon_{i}(s_k)$ and --->
<!--- $y_i = x_i\tau+\boldsymbol{w}_i^\top\boldsymbol{\kappa}+g\left(\int_0^1 \psi_i(s)\beta(s)ds\right)+\delta_i,$ --->
<!--- where $k=1,\cdots,n_\upsilon,~i=1,\cdots,n$.--->
