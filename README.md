# README

João Gonçalves 2024-03-27

# ThermalRaster

![](man/figures/ThermalRaster_Logo_250.png)

## Package description

**ThermalRaster** package is designed to enhance the processing and analysis of FLIR thermal and RGB images. It offers a comprehensive suite of functions for extracting, processing, and visualizing **thermal data** embedded within the metadata of FLIR images. Key features also include retrieving full and cropped **RGB images**. Imagery is returned as **`SpatRaster`** objects compatible with the **`terra`** package, enabling to profit from this package features/toolkit.

Using the overlap between low-resolution thermal imagery with high-resolution RGB images, the package enables the creation of synthetic or predicted thermal images for either cropped or entire RGB images. This is achieved through the application of the Random Forest algorithm (via the **`ranger`** package) or Deep Learning methodologies (utilizing **`keras`**/**`tensorflow`**).

The package also enables the handling of JSON annotations/masks from Roboflow (<https://roboflow.com>), enabling the extraction of ROIs from the images for further analysis, making possible to assess, plot, analyze and model fine-scale thermal variation in micro-habitats (i.e., TReMs - Tree Related Micro-habitats). Roboflow advantages include improved ROI digitization using manual or the SAM (Facebook Segment Anything) algorithm. ROIs can also be generated with the terra package as **`SpatVector`** objects. EXIFtool ([https://exiftool.org/](https://roboflow.com)) is broadly used for retrieving metadata from FLIR imagery.

------------------------------------------------------------------------

``` r
library(ThermalRaster)
library(terra)

image_path <- system.file("extdata", "BEECH1_EINBRGB_IPUQ4327.JPG", package = "ThermalRaster")

exiftool_path = "C:/MyFiles/R-dev/Giant_Trees"
```

Plot the FLIR image thumbnail

``` r
thumb <- flir_thumbnail_to_rast(image_path)

plotRGB(thumb)
```

![](README_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
rgb_lr <- flir_rgb_to_rast(image_path, exiftool_path, crop=TRUE)
plotRGB(rgb_lr)
```

![](README_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->
