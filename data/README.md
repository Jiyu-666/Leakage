# NANOGrav 15-year sky positions

`nanograv15_67_positions.csv` contains the 67-pulsar GWB analysis array from
the public NANOGrav 15-year tutorial Feather metadata.  The angles are the
published `theta` (colatitude) and `phi` (right ascension), in radians.

`nanograv15_excluded_position.csv` contains J0614-3329.  The official NG15
analyses exclude it because its timing baseline is shorter than three years.
The target paper nevertheless explicitly says that its signal-generation step
uses all 68 wideband pulsars.  Therefore the reproduction keeps two named
branches: 68 pulsars is the literal-paper branch and 67 pulsars is the standard
NG15-analysis sensitivity branch.

Sources:

- https://github.com/nanograv/15yr_stochastic_analysis
- https://doi.org/10.3847/2041-8213/acdc91
- https://fermi.gsfc.nasa.gov/ssc/data/access/lat/3rd_PSR_catalog/3PC_HTML/J0614-3329.html
