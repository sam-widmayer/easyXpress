# easyxpress

easyxpress is specialized for use with worms and is specific to use in the Andersen Lab and, therefore, is not available from CRAN. To install easyxpress you will need the [`devtools`](https://github.com/hadley/devtools) package. You can install `devtools` and `easyxpress` using the commands below:

```r
install.packages("devtools")
devtools::install_github("AndersenLab/easyxpress")
```

The functionality of the package can be broken down into three main goals:

+ Reading data and diagnositic images generated from CellProfiler pipelines alongside information about strains, conditions, controls, and experimen tal design.

+ Flagging and pruning anomalous data points

## Directory structure

Because so much information must be transferred alongside the plate data, the directory structure from which you are reading is critically important. Below is an example of a correct directory structure. The `data` directory contains`.csv` files from the CellProfiler pipeline copied directly from the default output folder for CellProfiler. Depending on the CellProfiler pipeline you run you may have multiple NonOverlappingWorms files. In the example below there are two such files. The `processed_images` directory contains `.png` files from the CellProfiler pipeline. There should be one `.png` file for each well included in your analysis. The `design.csv` file contains all the variables necessary to describe your experiment, i.e., drug names, drug concentrations, strain names, food types, ... The `.cpproj` file is not used by easyxpress, but is necessary to reproduce the original analysis, so easyxpress checks to make sure this file is present before it processes the data.

```
processed_images/projects/20190920_15hHB101foodpreps_RUN1
├── data
│   ├── dual_modelNonOverlappingWorms_control.csv
│   └── dual_modelNonOverlappingWorms_full.csv
└── processed_images
    └── 20190913_15hHB101foodpreps_01_2x_A01_w1_overlay.png
│   ├── 20190913_15hHB101foodpreps_01_2x_A02_w1_overlay.png
│   ├── 20190913_15hHB101foodpreps_01_2x_A03_w1_overlay.png
│   ├── ...    
├── design.csv
├── 20190920_15hHB101foodpreps_RUN1.cpproj
```

### Experiment directory

The experiment directory contains all of the files attached to a specific experiment conducted on a specific date. The naming convention for these folders should include the date in the format 4-digit year::2-digit month::2-digit day and experiment name separated by underscores. 

```
# Example directory name
# Date is January 1st, 2020
# Experiment name is "ExperimentName"

20200101_ExperimentName/
```

### File naming

The processed image files should be formatted with the experiment data, name of the experiment, the plate number, the magnification used for imaging, the experimental timepoint (if needed), the well name, and the wavelength used.name of the conditions template file, and name of the controls template file all separated by underscores. All processed image files must be saved as `.png` files. In the file named `20190913_15hHB101foodpreps_01_2x_A01_w1_overlay.png` the first section `20190913` is the experiment date, `15hHB101preps` is the name of the experiment, `01` is the plate number, `2x` is the magnification used for imaging, `A01` is the well name, and `w1` is the wavelength.

## Pipeline

The complete easyxpress package consists of only eight functions: `read_data`, `remove_contamination`, `flag_data`, `view_worm`, `sum_plate`, `bamf_prune`, and `outlier_prune`, and `regress`.

### `read_data()`

`read_data()` can take as an argument a path to a single data file, a directory with sorter data files, or an experiment directory with both setup and score subdirectories. If the function is given either a single file or directory of files, it will output a single data frame of all of the data that were read. If the function is given an experiment directory with both setup and score subdirectories, it will output a two element list with the first element being the score data and the second element being the setup data.

For further information use the command `?read_plate` to access the documentation.

### `read_data()`

### `remove_contamination()`

### `flag_data()`

### `view_well()`
If the `remote_server` parameter for this function is set to `TRUE` the `RCurl` package with SFTP protocol support is required. To ensure SFTP protocol support follow these instructions. 

1. In terminal install home brew if necessary (https://brew.sh/)
    * —  /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
2. In terminal install curl with SFTP support (https://github.com/marcelomazza/homebrew-curl-libssh2)
    * — brew install marcelomazza/homebrew-curl-libssh2/curl
3. In R, update PATH before installing RCurl. This only effects R session.
    * — Sys.setenv(PATH=paste('/usr/local/opt/curl/bin', Sys.getenv('PATH'), sep=":"))
4. In R, confirm that new PATH looks for curl in /usr/local/opt/curl/bin first.
    * — Sys.getenv("PATH")
5. In R, install RCurl from source
    * — install.packages("RCurl", type = "source")
6.  In R, load RCurl package and check for sftp protocol support
    * — library(RCurl)
    * — RCurl::curlVersion()

### `view_dose()`

### `sum_plate()`

### `bamf_prune()`

### `outlier_prune()`

### `regress()`

### Overview



### Example

```r
library(easyxpress)
library(dplyr)

# Define a vector of your experiement directories
dirs <- c("your_directory1",
          "your_directory2")

# Read in the data
raw <- read_data(dirs)

# Remove all data from the contaminated wells
raw_nocontam <- remove_contamination(raw)

# Flag suspect data points within wells
raw_nocontam_flagged <- flag_data(raw_nocontam)

# Review suspect data with image overlay
p1 <- view_worms(raw_nocontam_flagged, directory = "your_directory", plate = "your_plate", well = "your_well")
p1

...

# Save the final data frame
...
```# easyXpress