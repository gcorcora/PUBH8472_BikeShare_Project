library(dplyr)

data <- read.csv("~/Documents/Spatial Biostats/bikeshare_data.csv")

unique(data$from_station_id)

# create dataset with counts or averages for each starting station
summary_data <- data %>% 
  mutate(subscribed = ifelse(usertype == "Subscriber", 1, 0)) %>% 
  mutate(customer = ifelse(usertype == "Customer", 1, 0)) %>% 
  mutate(dependent = ifelse(usertype == "Dependent", 1, 0)) %>%
  mutate(male = ifelse(gender == "Male", 1, 0)) %>% 
  mutate(female = ifelse(gender == "Female", 1, 0)) %>% 
  mutate(age = 2015 - birthyear) %>% 
  group_by(from_station_id) %>% 
  summarise(trip_count = n(),
            capacity = first(from_station_dpcapacity),
            longitude = first(from_station_longitude),
            latitude = first(from_station_latitude),
            avg_duration = mean(tripduration),
            subscribers = sum(subscribed),
            customers = sum(customer),
            dependents = sum(dependent),
            male = sum(male),
            female = sum(female),
            avg_age = mean(age, na.rm = TRUE),
            avg_trip_distance = mean(trip_distance_manhattan))

write.csv(summary_data, "from_station_summary_data.csv", row.names = FALSE)

# create dataset with counts or averages for each ending station
end_summary_data <- data %>% 
  mutate(subscribed = ifelse(usertype == "Subscriber", 1, 0)) %>% 
  mutate(customer = ifelse(usertype == "Customer", 1, 0)) %>% 
  mutate(dependent = ifelse(usertype == "Dependent", 1, 0)) %>%
  mutate(male = ifelse(gender == "Male", 1, 0)) %>% 
  mutate(female = ifelse(gender == "Female", 1, 0)) %>% 
  mutate(age = 2015 - birthyear) %>% 
  group_by(to_station_id) %>% 
  summarise(trip_count = n(),
            capacity = first(to_station_dpcapacity),
            longitude = first(to_station_longitude),
            latitude = first(to_station_latitude),
            avg_duration = mean(tripduration),
            subscribers = sum(subscribed),
            customers = sum(customer),
            dependents = sum(dependent),
            male = sum(male),
            female = sum(female),
            avg_age = mean(age, na.rm = TRUE),
            avg_trip_distance = mean(trip_distance_manhattan))

write.csv(end_summary_data, "end_station_summary_data.csv", row.names = FALSE)

