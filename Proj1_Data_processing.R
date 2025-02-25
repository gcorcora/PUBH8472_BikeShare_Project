#loading divvy station data
bikeshare_stations <- read.csv('/Users/gretchen/Desktop/divvy_data/Divvy_Stations.csv')
trips_q1_data <- read.csv('/Users/gretchen/Desktop/divvy_data/Divvy_Trips_2015-Q1.csv')
#trips from each station - collect records that are trips originating from stations
#Origin station ID and trip length, lat/long at station location
#the thing on the left is individual trip records
#join station info onto trips
#left join, where left component is all trips originating from stations and right piece is station location
#more records of trips than stations, replicate station locations to all trip records
#then do a count with a group_by
#count # of trips, group_by station id - then get origin trips by station
#then can get same for trips that end at stations

library(dplyr)

#merging from_station information
bikeshare_merged_data <- trips_q1_data %>% left_join(
  bikeshare_stations, by = c("from_station_id" = "id"))

#data_processing, renaming variables, removing redundant info
bikeshare_merged_data <- bikeshare_merged_data %>% select(-name)

bikeshare_merged_data <- bikeshare_merged_data %>% rename(from_station_latitude = latitude)
bikeshare_merged_data <- bikeshare_merged_data %>% rename(from_station_longitude = longitude)
bikeshare_merged_data <- bikeshare_merged_data %>% rename(from_station_dpcapacity = dpcapacity)
bikeshare_merged_data <- bikeshare_merged_data %>% rename(from_station_online_date = online_date)

#merging to_station information
bikeshare_merged_data <- bikeshare_merged_data %>% left_join(bikeshare_stations, by = c("to_station_id" = "id"))

#data processing - rename vars, removing redundant info
bikeshare_merged_data <- bikeshare_merged_data %>% select(-name)

bikeshare_merged_data <- bikeshare_merged_data %>% rename(to_station_latitude = latitude)
bikeshare_merged_data <- bikeshare_merged_data %>% rename(to_station_longitude = longitude)
bikeshare_merged_data <- bikeshare_merged_data %>% rename(to_station_dpcapacity = dpcapacity)
bikeshare_merged_data <- bikeshare_merged_data %>% rename(to_station_online_date = online_date)

#calculating distance with manhattan distance because street grids etc. make sense
manhattan_distance <- function(lat1, long1, lat2, long2){
  abs(lat2 - lat1) + abs(long2 - long1)
}

#creating trip distances (will later have to calculate distance between stations for analysis)
bikeshare_merged_data <- bikeshare_merged_data %>%
  mutate(trip_distance_manhattan = manhattan_distance(from_station_latitude, from_station_longitude, to_station_latitude, to_station_longitude))

#saving locally
write.csv(bikeshare_merged_data, "bikeshare_data.csv", row.names = FALSE)
