# Mammal diversity on islands worldwide

---

Data from: Barreto Elisa, Rangel Thiago F., Pellissier Loc and Graham Catherine H. 2021Area, isolation and climate explain the diversity of mammals on islands worldwideProc. R. Soc. B.2882021187920211879 [http://doi.org/10.1098/rspb.2021.1879](http://doi.org/10.1098/rspb.2021.1879)

This data contains a presence and absence matrix of mammal species per island and a compilation of mammalian species richness, endemism and a series of physical and environmental characteristics of each island. A total of 5,592 islands (out of the ~17,000 islands larger than 1km worldwide; Weigelt et al., 2013) are included in the data.

## Versioning

June 2024 update: Appendix 1 with the presence-absence matrix was updated to remove 15 duplicated entries. The file "mammal_insularity_from_IUCN.csv" was added.

## Description of the Data

We derived aglobal database of mammal species on over 5,500 islands  worldwide by overlapping island shapefiles from theGlobal Administrative areas version 3.6 (GADM, 2018) and mammalian range maps from IUCN (IUCN, 2017). We considered as islands all the land masses smaller than Greenland (2,166,000 km) that are surrounded by salty water. We carefully and extensivelly inspected and manually corrected any alignment inconsistencies between both shapefiles (i.e., islands and mammal range) using QGIS 3.6 (Open Source Geospatial Foundation Project, 2019). We opted for a highly conservative approach of excluding any island with the slightest doubt about species attribution and ignoring islands where no mammal species occurs according to the IUCN data (i.e., our dataset only includes islands with at least one species). For example, regions with clusters of nearby islands - e.g., Patagonia and Scandinavia.We removed introduced species from the database by excluding (1) species polygons recorded as introduced by IUCN and (2) species listed as invasive for each particular island in the Database of Island Invasive Species Eradication (DIISE, 2015). We also removed fully aquatic and marine semi-aquatic species. We ensured that native species that were extinct due to human activity were included in the database by adding occurrence records from (Faurby & Svenning, 2016; Upham, 2017; Faurby et al., 2018). For each island in the database, we gathered a series of environmental and physical characteristics expected to influence biodiversity.

## File structure

Appendix 1: presence and absence matrix of species per island, with island ID in rows and mammal species in columns.

Appendix 2: mammalian biodiversity and a series of phisical and environmental characteristic of each island. ID column matches the row names in the presence and absence matrix (appendix 1). This data set includes: richness of native mammal species, number of single island mammal endemics (SIE), and proportion of SIE, island's mean annual temperature (in degrees Celsius), annual precipitation (in millimeters), standard deviation of mean annual temperature and precipitation within the island, standard deviation in elevation within the island, area (in km), surrounding landmass proportion (SLMP), island connectivity to the mainland during the last glacial maximum (GMMC), climate change velocity in temperature since the last glacial maximum (CCVT, in meters/year), and realm.

We derived temperature and precipitation data from CHELSA using monthly estimates across the years 1979 to 2013 (Karger et al., 2017) and elevation from the Global Digital Elevation Model GTOPO30 (USGS, 1996) and calculated mean and standard deviation per island using QGIS 3.6 (Open Source Geospatial Foundation Project, 2019).We obtained island area, SLMP, GMMC and CCVT from a public island characterization database (Weigelt et al., 2013) by matching the centroid coordinates to the island polygons. SLMP is a proxy of island isolation with great predictive power and was calculated as the log10-transformed sum of the proportion of surrounding landmass within buffer distances of 100, 1,000 and 10,000 km around each island perimeter (Weigelt & Kreft, 2013). GMMC is a binary descriptor of historical isolation that uses past and present global bathymetry data to infer if islands were connected to the continent during the last glacial maximum by assuming the estimated sea level decrease of 122m at 18,000 years ago (more details in Weigelt et al., 2013). We multiplied SLMP by 1 and coded GMMC as 0 being connected and 1 being disconnected to the mainland during the LGM, so both metrics represent isolation (i.e., higher SLMP and GMMC represent greater isolation). CCVT over the past 21,000 years was calculated by dividing the difference in mean annual temperature between past and present by the spatial change in present mean temperature (Loarie et al., 2009; Weigelt et al., 2013). CCVT is interpreted as the speed at which the organism would have to move to keep pace with historical temperature change, assuming no change in topography (Loarie et al., 2009; Weigelt et al., 2013). Islands were classified into the 12 global mammalian zoogeographical regions (Holt et al., 2013), hereafter realm. We removed 505 islands from the dataset because it was not possible to derive all environmental variables or to assign a realm with confidence, usually because they were small (< 1km) or located on a biogeographical boundary.

mammal_insularity_from_IUCN.csv: list of all the mammal species in the presence and absence matrix and a categorization indicating whether the species is only "continental", "continental_and_insular" or only "insular".

## Sharing/access Information

DIISE (2015) The Database of Island Invasive Species Eradications, developed by Island Conservation, Coastal Conservation Action Laboratory UCSC, IUCN SSC Invasive Species Specialist Group, University of Auckland and Landcare Research New Zealand.

Faurby, S., Davis, M., Pedersen, R., Schowanek, S.D., Antonelli, A. & Svenning, J.C. (2018) PHYLACINE 1.2: The Phylogenetic Atlas of Mammal Macroecology. Ecology, 99, 2626.

Faurby, S. & Svenning, J.-C. (2016) Resurrection of the Island Rule: Human-Driven Extinctions Have Obscured a Basic Evolutionary Pattern. The American Naturalist, 187, 812820.

GADM (2018) GADM database of Global Administrative Areas, version 3.6.

Hanna, E. & Cardillo, M. (2014) Island mammal extinctions are determined by interactive effects of life history, island biogeography and mesopredator suppression. Global Ecology and Biogeography, 23, 395404.

Hbert, K., Millien, V., Lessard, J. & Masters, J. (2021) Source pool diversity and proximity shape the compositional uniqueness of insular mammal assemblages worldwide. Journal of Biogeography, jbi.14156.

Holt, B.G., Lessard, J.-P., Borregaard, M.K., Fritz, S.A., Arajo, M.B., Dimitrov, D., Fabre, P.-H., Graham, C.H., Graves, G.R., Jnsson, K.A., Nogus-Bravo, D., Wang, Z., Whittaker, R.J., Fjelds, J. & Rahbek, C. (2013) An update of Wallaces zoogeographic regions of the World. Science, 339, 7478.

IUCN (2017) IUCN Red List of threatened species  mammal range polygons.

Karger, D.N., Conrad, O., Bhner, J., Kawohl, T., Kreft, H., Soria-Auza, R.W., Zimmermann, N.E., Linder, H.P. & Kessler, M. (2017) Climatologies at high resolution for the earths land surface areas. Scientific Data, 4, 120.

Lavery, T.H., Olds, A.D., Seddon, J.M. & Leung, L.K.-P. (2016) The mammals of northern Melanesia: Speciation, ecology, and biogeography. Mammal Review, 46, 6076.

Loarie, S.R., Duffy, P.B., Hamilton, H., Asner, G.P., Field, C.B. & Ackerly, D.D. (2009) The velocity of climate change. Nature, 462, 10521055.

Meiri, S., Cooper, N. & Purvis, A. (2008) The island rule: made to be broken? Proceedings of the Royal Society B: Biological Sciences, 275, 141148.

Millien-Parra, V. & Jaeger, J.-J. (1999) Island biogeography of the Japanese terrestrial mammal assemblages: an example of a relict fauna. Journal of Biogeography, 26, 959972.

Open Source Geospatial Foundation Project (2019) QGIS 3.6.

Upham, N.S. (2017) Past and present of insular Caribbean mammals: Understanding Holocene extinctions to inform modern biodiversity conservation. Journal of Mammalogy, 98, 913917.

USGS (1996) GTOPO 30 - Global Digital Elevation Model. Sioux Falls, SD: US Geological Survey.

Weigelt, P., Jetz, W. & Kreft, H. (2013) Bioclimatic and physical characterization of the worlds islands. Proceedings of the National Academy of Sciences, 110, 1530715312.

Weigelt, P. & Kreft, H. (2013) Quantifying island isolation - insights from global patterns of insular plant species richness. Ecography, 36, 417429.

## Queries

Drop me an email in case of questions: [elisabpereira@gmail.com](mailto:elisabpereira@gmail.com)
