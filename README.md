# Work on Thüringia existing Tiles

A high-resolution, consistent, and nationwide Airborne Laser Scanning (ALS) dataset is essential for accurately deriving forest metrics such as **canopy cover**, **top height**, **growing stock**,
and **vertical layering** across Germany, as emphasized by the **ForestPulse** project. ALS is uniquely suited to this task because it penetrates vegetation layers and captures multiple echoes from
vertical forest profiles, providing the detailed three-dimensional structural information necessary to analyze complex forest canopies and understory capabilities that are critical 
for comprehensive forest assessment and monitoring at scale.


However, the anticipated nationwide ALS dataset from the **Federal Agency for Cartography and Geodesy (BKG)** is expected to be publicly available 
**soon**. Given this delay, we have initiated a **state-level approach** by retrieving freely accessible ALS data from the 
**geoportal of Thueringen** (*Geoportal Thüringen*), which serves as a valuable pilot region for early methodological development and demonstration.

The data were accessed via the following ATOM feed provided by the Thüringen State authority for Land Management and Geoinformation:

[https://geoportal.thueringen.de/gdi-th/download-offene-geodaten/download-hoehendaten]

From this portal, we systematically downloaded all available LAZ tiles that spatially cover the **entire state of Thüringia**. The tiles were organized by acquisition periods,
and our first priority was to utilize the most **recent and high-quality data (2020–2025 folder)** to ensure up-to-date assessments. However, since this folder is not yet 
complete in coverage, we supplemented the missing spatial areas using the **older dataset (2014–2019 folder)**. 
This hybrid temporal approach ensures complete coverage of Thuringia while maintaining as much recency as possible.

Due to the large volume of data involved, covering an entire federal state at high resolution (about 30 points/m²), we used storage and computational resources 
provided by the **GWDG Göttingen** platform. This infrastructure enabled both the **storage of raw LAS files** and their subsequent **preprocessing workflows**, 
which included:

* **Harmonization** of data formats and metadata,
* **Denoisation** to remove erroneous returns (e.g., atmospheric particles, birds),
* **Classification** of point clouds into ground, vegetation, and other relevant classes, and
* **Reformatting** the data into **COPC (Cloud Optimized Point Cloud)** format for efficient downstream processing and visualization.

![Flowchart-2](https://github.com/user-attachments/assets/ab211933-0c7b-4a6b-8f26-ba8e01b64629)


The preprocessed dataset, now organized in COPC format, forms the basis for our calculation of the targeted forest attributes:

* **Canopy cover**, estimated through point cloud or raster-based metrics derived from first returns,
* **Top height**, calculated as the 95th percentile of vegetation height (H95),
* **Growing stock**, inferred using statistical models that relate ALS-derived structural metrics to field inventory data, and
* **Vertical layering**, characterized using vertical foliage distribution profiles and relative density metrics.

This state-level analysis serves as both a proof-of-concept and preparatory work for the upcoming nationwide implementation, once the BKG dataset becomes available. 
The methodology developed here can be readily scaled to other federal states as data availability improves.

**References**

* Goodbody, T. R. H., Coops, N. C., Senf, C., & Seidl, R. (2023). Airborne laser scanning to optimize the sampling efficiency of a forest management inventory in South-Eastern Germany. Ecological Indicators, 157, Artikel 111281. https://doi.org/10.1016/j.ecolind.2023.111281

* Margaret Penner, Joanne C White, Murray E Woods, Automated characterization of forest canopy vertical layering for predicting forest inventory attributes by layer using airborne LiDAR data, Forestry: An International Journal of Forest Research, Volume 97, Issue 1, January 2024, Pages 59–75, https://doi.org/10.1093/forestry/cpad033

