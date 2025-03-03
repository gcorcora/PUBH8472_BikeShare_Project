# PUBH8472_BikeShare_Project

Files for Project 1 by Ting-Chien Huang and Gretchen Corcoran:
  1. Initial processing file:
     File: bikeshare_data.csv.zip
     Code: init_data_processing_part1.R
  3. Secondary processing data files:
     File: from_station_summary_data.csv.zip (This is the one we analyzed)
     File: end_station_summary_data.csv (We didn't analyze this one)
     Code: init_data_processing_part1.R
  3. Main analysis:
     main_analysis.R - contains 3 models, in order within the file:
         1. First model: Analysis with outlier
         2. Second model: Intercept only model, used for kriging - this takes a long time to run, be careful!
         3. Third model: Analysis without outlier
  4. Supplementary files:
       Two maps made with QGIS to show the bike stations in their real locations
         qgis_map_fix_with_outlier.png
         qgis_map_no_outlier.png
