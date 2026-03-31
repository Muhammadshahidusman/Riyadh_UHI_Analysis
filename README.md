# Riyadh Urban Heat Island & Land Use/Land Cover Analysis

A comprehensive geospatial analysis project examining the relationship between Land Use/Land Cover (LULC) changes and Urban Heat Island (UHI) effects in Riyadh, Saudi Arabia.

## Project Overview

This project analyzes multi-temporal satellite imagery (1990-2025) to:
- Classify and track Land Use/Land Cover (LULC) changes
- Quantify Urban Heat Island (UHI) intensity variations
- Correlate LULC transformations with temperature anomalies
- Perform spatial autocorrelation analysis (Getis-Ord Gi* hot spot analysis)

## Author

**Muhammad Shahid Usman**
- GeoSpatial Data Analyst
- Riphah International University

## Dataset Structure

```
Datasets/
├── LUL_Classified/       # Classified LULC images (1990-2025)
│   ├── LULC_1990.tif
│   ├── LULC_1995.tif
│   ├── LULC_2001.tif
│   ├── LULC_2014.tif
│   ├── LULC_2019.tif
│   └── LULC_2025.tif
├── LST/                  # Land Surface Temperature
│   ├── lst_1990.tif
│   ├── lst_1995.tif
│   ├── lst_2001.tif
│   ├── lst_2014.tif
│   ├── lst_2019.tif
│   └── lst_2025.tif
├── UHI_Classified/       # Reclassified UHI images
│   └── Reclass_UHI_*.tif
└── Riyadh_Boundary/      # Study area boundary shapefile
    └── Riyadh City.shp
```

### LULC Classification Scheme

| Code | Class |
|------|-------|
| 0 | Water |
| 1 | Built-Up |
| 2 | Vegetation |
| 3 | Barren Land |

### UHI Classification Scheme

| Code | Class |
|------|-------|
| 1 | Extreme Cool Island |
| 2 | Strong Cool Island |
| 3 | Moderate Cool |
| 4 | Moderate Hot |
| 5 | Strong Hot |
| 6 | Extreme Hot |

## Analysis Methods

### Statistical Analysis (main01.py)

- Area calculation for each LULC class per year
- LULC transition matrix (1990-2025)
- Cross-tabulation analysis (LULC vs UHI)
- Pearson correlation analysis
- Land Surface Temperature (LST) integration
- Getis-Ord Gi* hot spot analysis for spatial clustering

### Google Earth Engine Scripts (GEE_Scripts/)

Remote sensing scripts for LULC and UHI classification:
- `LULC_1990.js` through `LULC_2025.js` - Land cover classification
- `UHI_1990.js` through `UHI_2025.js` - Urban Heat Island classification

## Requirements

- Python 3.8+
- rasterio
- geopandas
- numpy
- pandas
- matplotlib
- scipy
- libpysal
- esda
- openpyxl

## Usage

Run the main statistical analysis:

```bash
python main01.py
```

This generates:
- `Riyadh_Statistical_Analysis_5_6.xlsx` - Comprehensive statistical results
- `maps/` directory - Generated visualizations

## Key Findings

The analysis covers temporal changes over 35 years (1990-2025), examining:
- Urban expansion patterns
- Vegetation loss/gain dynamics
- Heat island intensity trends
- Spatial hot spot evolution

## License

This project is for academic research purposes.