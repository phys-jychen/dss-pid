# Change Log

## 23 February 2024 (Version 1.1.0)

Creation of the project.

## 27 February 2024

Added event display. The staggered structure has been taken into consideration.

## 2 March 2024

Cluster information is saved and used in BDT.

## 10 April 2024

Introduced energy projection to three axes and three planes.

## 18 July 2024

Significant modifications of the definitions.

- The scales of FD are restricted to below the cell number in $x$ and $y$ directions.
- Some 2D variables have been changed to 3D, e.g. energy ratios, shower radius;
- Staggered structure has been taken into account, e.g. 3 × 3 cells in the same layer of the hit, but 2 × 2 cells in consecutive layers;
- The event axis is obtained from 3D linear fitting instead of simply defined as parallel to $z$ axis.

## 22 July 2024 (Version 1.2.0)

Added a function to dump the true positive (TP) and true negative (TN) values for plotting ROC curves.

## 28 July 2024 (Version 1.3.0)

- Modified the BDT part to accelerate signal–background BDT.
- Split the `EventNumber` into two parts, one for storing the 5 lowest digits, and the other for storing the highest digits. This is designed to overcome the shortcomings of TMVA package that spectators are also automatically converted to `float` type.

## 3 August 2024 (Version 2.0.0)

- Changed from signal–background separation to multi-class BDT, and added weights for each process. The weights come from [the publication in 2023](https://doi.org/10.1007/s11433-022-1983-8).
- Removed some variables that are strongly correlated with others.
- Modified the contents dumped to the CSV file.

## 31 August 2024

Combine the ROOT files of each class before sending into BDT to speed up and reduce memory consumption.

## 9 September 2024 (Version 2.1.0)

- Simplified the selection process, i.e. removed `iPID -s`.
- Added more variables from other sub-detectors for BDT.
- To-do: check the selection criteria of HCAL energy terms.