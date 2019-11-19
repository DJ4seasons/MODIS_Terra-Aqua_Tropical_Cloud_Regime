Tropical Cloud Regime from MODIS Terra-Aqua
---

#### Tropics defined for deriving cloud regime: 15S-15N

#### Data used for the cloud regime
- 2-D Joint Histogram of cloud optical thickness (COT or TAU) and top height (CTP)
- nCTP=7, nCOT=6, so total 42 elements
- Both MODIS Terra and Aqua
- 1deg x 1deg resolution, 2002.12 - 2016.11

#### Method for deriving the tropical cloud regime (TCR)
- k-means clustering, Final k=10
- One of 10 regimes is composed of low cloud fraction with significantly large population (about 40%)
- This one (we set as TCR10) is used for nested clustering (i.e., one more round of k-means only for TCR10s)
- Final sub_k=4

#### Troical Cloud Regime (TCR) Centroids (= means of similar histograms assigned for each regime)  
- 10 Centroids
- "Data/Centroid.MODIS_T+A_b42_DTR_CR10.f64dat"
- [k=10, nCTP=7, nTAU=6], binary, float64 (=8 Byte)

- 4 Sub-Centroids of TCR10
- "Data/Centroid.MODIS_T+A_b42_DTR_CR10.subCR4.f64dat"
- [k=4, nCTP=7, nTAU=6], binary, float64 (=8 Byte)

#### Centroid Related Python Program
- "centroid_read_bin+write_nc.py3.py": Transform binary centroid file to NetCDF format
- "cent_display_42bins_obs.py3.py": Display centroids
- "cent_display_42bins_obs.subCR.py3.py": Display sub-centroids of TCR10

#### TCR Relative Frequency of Occurrence (RFO) Map
- Each grid cell for all time/all domain is assigned to one of TCRs (a.k.a. CR_map file)
- By averaging over time, we can examine the geographycal preferance of each TCR.
- This is called as RFO map.
- For this purpose, the MODIS data is updated as C6 -> C6.1, 15S-15N -> 25S-25N, [2002.12-2016.11] -> [2003.01-2017.12]

#### RFO related Python Program
- "CRmap_read_bin+write_nc.py3.py": Transform binary CR_map file to NetCDF format
- "RFO_Map_plot_obs.py3.py": Display RFO maps for 10 TCRs
- "RFO_Map_plot_obs.subCR.py3.py": Display RFO maps for 4 sub-TCRs of TCR10

#### TCR Aggregate
- Liberally identified meaning of TCRs: TCR1(convective core-dominant), TCR2(various stages of convections), TCR3(anvil-dominant)
- For TCR1 or TCR1, 2, and 3 as a group, how many they are aggregated is identified. 
- In the case of TCR123 aggregate, the existence of TCR1 is prerequisite.
- For a daily map, once a grid cell of TCR1 is given, recursively find neighbors if they satisfy the condition, TCR1 or TCR123.
- Neighbor is defined as one of four grid cells shares a side of center grid cell. 
- "find_TCR_aggregate.py3.py": Identify aggregates and print out characteristics to a text file
