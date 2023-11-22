# Packages and data import
library(dplyr)
library(tidyverse)
library(ggplot2)
library(readr)
library(lubridate)
library(sf)
library(rgdal)
library(cartography)
library(reshape2)
library(RColorBrewer)
library(R2WinBUGS)
library(ggspatial)
library(gganimate)
library(transformr)
library(gifski)
library(reshape2)
library(readxl)
library(geosphere)
library(factoextra)
library(conclust)

setwd("./data/")
Baboons <- read_excel("Baboons.xlsx")

# Data manipulation 
colnames(Baboons) <- c("Event", "StationID", "Nbp", "Start", "End", "Eventd", "Date", "Station", "Ind", "Male", "Female", "Adult", "SubAdult", 
                       "Juvenile", "UknA", "Activity", "Mouvement", "Comments", "Photo/Video", "Observation")
Baboons = Baboons %>%
  filter(!is.na(Baboons$Male))
  
events_per_station <- table(Baboons$`StationID`)
events_per_station <- as.data.frame(events_per_station)
names(events_per_station) <- c("Station ID", "Number of Events")

events_per_station$`Station ID` <- with(events_per_station, reorder(`Station ID`, -`Number of Events`))

# Number of events per station histogram
ggplot(events_per_station, aes(x = `Station ID`, y = `Number of Events`)) +
  geom_bar(stat = "identity", fill = "#8888ff") +
  theme_minimal() +
  labs(x = "Station ID", y = "Number of Events", title = "Histogram of the Number of Events per Station") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# Reducing the dataset to columns of interest
Baboons2 <- Baboons %>%
  select(Event, Station, Date, Ind) %>%
  mutate(
    Date = ymd(Date),
    Week = floor(as.numeric(difftime(Date, first(Date), units = "weeks"))),
    Year = year(Date),
    Month = month(Date),
    Day = day(Date),
    Nbevent = 1  
  ) %>%
  arrange(Date) 

# Per day
by_day <- Baboons2 %>%
  count(Station, Date, name = "Nbevent") %>%
  spread(key = Date, value = Nbevent, fill = 0)

# Per week
by_week <- Baboons2 %>%
  count(Station, Week, name = "Nbevent") %>%
  spread(key = Week, value = Nbevent, fill = 0)

# Per month
by_month <- Baboons2 %>%
  count(Station, Month, Year, name = "Nbevent") %>%
  unite("Month_Year", c("Month", "Year"), sep = "-") %>%
  spread(key = Month_Year, value = Nbevent, fill = 0) %>%
  mutate(Tot = rowSums(.[, -c(1:2)], na.rm = TRUE))  

# Per year
by_year <- Baboons2 %>%
  count(Station, Year, name = "Nbevent") %>%
  spread(key = Year, value = Nbevent, fill = 0) %>%
  mutate(Tot = rowSums(.[, -c(1)], na.rm = TRUE))

# Coordinates per station
Stations = st_read(dsn="QGIS_data", layer="Stations")
coords <- Stations %>%
  group_by(Stations) %>%
  summarise(Latitude = mean(Latitude, na.rm = TRUE),
            Longitude = mean(Longitude, na.rm = TRUE)) %>%
  rename(Station = Stations)

calculate_distance <- function(lat1, long1, lat2, long2) {
  distm(c(long1, lat1), c(long2, lat2), fun = distHaversine)
}

# Map
Terrain = st_read(dsn= "QGIS_data", layer = "Terrain")
Stations = Stations %>%
  filter(Nbevent != 0)
colnames(Stations) = c("ID","Station_ID","Date","Nbevent","Latitude","Longitude","geometry")

ggplot()+
  geom_sf(data = Terrain) +
  geom_sf(data=Stations) + 
  geom_sf_label(data = Stations, aes(label = Station_ID), 
                nudge_x = 0, nudge_y = 0.01, 
                check_overlap = TRUE, color = '#139339', size = 2) +
  labs(x="Longitude", y="Latitude", title = "Carte de la zone d'étude et des stations")+
  annotation_north_arrow(location = "br", which_north = "true", pad_x = unit(0.1, "in"), pad_y = unit(0.1, "in"), style = north_arrow_fancy_orienteering) +
  theme_minimal()

# For spatial distance matrix
spatial_distances <- data.frame(Event = integer(), Distance = numeric(), stringsAsFactors = FALSE)
for(i in 1:(nrow(Baboons) - 1)) {
  lat1 <- coords$Latitude[coords$Station == Baboons$Station[i]]
  long1 <- coords$Longitude[coords$Station == Baboons$Station[i]]
  lat2 <- coords$Latitude[coords$Station == Baboons$Station[i + 1]]
  long2 <- coords$Longitude[coords$Station == Baboons$Station[i + 1]]
  distance <- calculate_distance(lat1, long1, lat2, long2)
  spatial_distances <- rbind(spatial_distances, data.frame(Event = Baboons$Event[i], Distance = distance))
}

# For temporal distance matrix
extract_time <- function(datetime) {
  format(as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S"), "%H:%M:%S")
}

Baboons <- Baboons %>%
  mutate(
    StartTime = paste(Date, extract_time(Start)),
    EndTime = paste(Date, extract_time(End)),
    # Convert the combined strings to POSIXct datetime objects
    StartTime = as.POSIXct(StartTime, format = "%Y-%m-%d %H:%M:%S"),
    EndTime = as.POSIXct(EndTime, format = "%Y-%m-%d %H:%M:%S")
  )

temporal_distances <- data.frame(Event = integer(), TimeDiff = numeric(), stringsAsFactors = FALSE)
for(i in 1:(nrow(Baboons) - 1)) {
  time1 <- Baboons$EndTime[i]
  time2 <- Baboons$StartTime[i + 1]
  # Calculate the difference in minutes
  timediff <- as.numeric(difftime(time2, time1, units = "mins"))
  temporal_distances <- rbind(temporal_distances, data.frame(Event = Baboons$Event[i], TimeDiff = timediff))
}

event_pairs <- data.frame(
  Event1 = Baboons$Event[-nrow(Baboons)],      
  Event2 = Baboons$Event[-1],                     
  SpatialDistance = spatial_distances$Distance,   
  TemporalDistance = temporal_distances$TimeDiff 
)

# To determine if consecutive events are trigerred by the same troop or not
event_pairs <- event_pairs %>%
  mutate(
    SameTroop = ifelse(SpatialDistance / 1000 <= TemporalDistance / 60 * 1.26, "Possible", "Non possible")
  )

# Graph
ggplot(event_pairs, aes(x = TemporalDistance / 60, y = SpatialDistance / 1000)) +
  geom_point(aes(color = SameTroop)) +
  geom_abline(intercept = 0, slope = 1.26, color = "#000000", size = 1) +
  labs(x = "Distance temporelle (heures)", y = "Distance spatiale (km)", 
       title = "Distance spatiale par distance temporelle pour chaque paire d'événements successifs") +
  xlim(0, 100) +
  ylim(0,25)

# K-means Clustering
Baboons$Prop_adult = Baboons$Adult / Baboons$Ind
Baboons$Prop_sub = Baboons$SubAdult / Baboons$Ind
Baboons$Prop_juv = Baboons$Juvenile / Baboons$Ind

Event_ID = Baboons$Event

Clustering_data = Baboons %>%
  select(Prop_adult,Prop_sub,Prop_juv) %>%
  scale()

set.seed(123)
rownames(Clustering_data) = Event_ID

fviz_nbclust(Clustering_data, kmeans, method = "wss") +
  geom_vline(xintercept = 3, color = "#ff0000")

kmeans_result = kmeans(Clustering_data, centers = 3, nstart = 25)

fviz_cluster(list(data = Clustering_data, cluster = kmeans_result$cluster)) 

########## MOD BAYESIEN ###########
# NOTE: To run the Bayesian model in this script, ensure that JAGS (Just Another Gibbs Sampler) is installed on your system. 
# Be advised that running the Bayesian model might take a significant amount of time. 
# If you prefer to avoid the download and the wait, you can directly refer to the 'Bayesian_results.csv' file in Data folder for the model's results.

library(rjags)
library(coda)

# Data import
read_data <- function(yobs_file, id_file) {
  event_name <- basename(yobs_file)
  event_name <- sub("\\.csv$", "", event_name)
  yobs <- read.csv2(yobs_file, sep = ";", header = TRUE)
  ID <- read.csv2(id_file, sep = ";", header = TRUE)
  return(list(event_name = event_name,yobs = yobs, ID = ID))
}

path = "CMR"
files <- list.files(path = path, pattern = ".csv", full.names = TRUE)
yobs = files[!grepl("ID\\.csv$", files)]
ID = list.files(path = path, pattern = "ID.csv", full.names = TRUE)
pairs = lapply(seq_along(yobs), function(i) {
  list(yobs=yobs[i], ID = ID[i]) })

all_data <- lapply(pairs, function(pair) {
  read_data(pair$yobs, pair$ID)
})

# Age class transformation (Adult : 1, Subadult : 2, Juvenils : 3, Unknown : 1)
for (i in 1:length(all_data)) {
  all_data[[i]][[2]] = all_data[[i]][[2]][,-1] 
  all_data[[i]][["ID"]][["Class"]] = as.character(all_data[[i]][["ID"]][["Class"]]) 
  for (j in 1:nrow(all_data[[i]][["ID"]])){  
    if (all_data[[i]][["ID"]][["Class"]][[j]] == "A") 
      {all_data[[i]][["ID"]][["Class"]][[j]] = "1"}
    
    if (all_data[[i]][["ID"]][["Class"]][[j]] == "S") 
      {all_data[[i]][["ID"]][["Class"]][[j]] = "2"}
    
    if (all_data[[i]][["ID"]][["Class"]][[j]] == "J") 
      {all_data[[i]][["ID"]][["Class"]][[j]] = "3"}
    
    if (all_data[[i]][["ID"]][["Class"]][[j]] == "U") 
      {all_data[[i]][["ID"]][["Class"]][[j]] = "1"}
  }
  IDaug = data.frame(ID=1:150, 
                     Class=sample(1:3, 150, replace=TRUE, prob=c(0.6, 0.3, 0.1)), 
                     Sex=rep("U",150)) # Data augmentation
  colnames(IDaug) = colnames(all_data[[i]][["ID"]])
  all_data[[i]][["ID"]] = rbind(all_data[[i]][["ID"]], IDaug) 
  yAug = as.data.frame(matrix(0,nrow = 150, ncol = ncol(all_data[[i]][["yobs"]])))
  colnames(yAug) = colnames(all_data[[i]][["yobs"]])
  all_data[[i]][["yobs"]] = rbind(all_data[[i]][["yobs"]],yAug)
}

# Jags model
model_string = "
model {
  # Priors
  omega ~ dunif(0, 1) # Uniform prior for the probability an individual in the augmented population belongs to the true population N
  for (i in 1:3) { 
    p[i] ~ dunif(0, 1)     # Uniform prior for the probability of capture if not previously captured
    c[i] ~ dunif(0, 1)      # Uniform prior for the probability of capture if previously captured
  }   
  
  # Likelihood
  for (i in 1:M) {     # M is the number of rows in the data table (win.data)
    z[i] ~ dbern(omega) # Bernoulli distribution for the indicator of whether an individual represents a real individual or not
    
    # First occasion
    yaug[i, 1] ~ dbern(p.eff[i, 1])   # Bernoulli distribution, p.eff[i, 1] for the probability of success
    p.eff[i, 1] <- z[i] * p[class[i]]
    
    # All subsequent occasions
    for (j in 2:T) {
      yaug[i, j] ~ dbern(p.eff[i, j])
      p.eff[i, j] <- z[i] * ((1 - yaug[i, (j-1)]) * p[class[i]] + yaug[i, (j-1)] * c[class[i]])
    } # j
  } # i
  
  # Derived quantities
  N <- sum(z) 
  for (i in 1:M) {
    y[i] <- step(1 - class[i])
    z1[i] <- z[i] * y[i]
    z3[i] <- z[i] * step(class[i] - 3)
  }
  N1 <- sum(z1)
  N3 <- sum(z3)
  N2 <- N - N1 - N3
} # end model
"

all_results = list()
params <- c("N", "p", "c", "omega", "N1", "N2", "N3")

ni <- 1000 # Nb of iterations
nt <- 2 # Thinning interval (improve independance of samples)
nb <- 500 # Burn-in, discard the first 500 iterations
nc <- 3 # Nb of chains, 3 random initial values, allow to check the convergence of the MCMC algorithm

# Running the JAGS model
for (i in 1:length(all_data)) {
  data <- all_data[[i]] 
  yobs = data[["yobs"]]
  class = as.integer(data[["ID"]][["Class"]])
  M <- nrow(yobs)
  T <- ncol(yobs)
  jags_data <- list(yaug = yobs, class = class, M = M, T = T)

  inits <- function() {
    list(z = rep(1, M), p = runif(3, 0, 1)) 
  }

  jags.model <- jags.model(textConnection(model_string), data = jags_data, inits = inits, n.chains = nc)
  update(jags.model, n.iter = nb)
  samples <- coda.samples(jags.model, variable.names = params, n.iter = ni, thin = nt)
  results_summary <- summary(samples)
  

  # Display results
  event_name = all_data[[i]][["event_name"]]
  cat("Résultats pour l'événement:", event_name, "\n")
  all_results[[event_name]] = results_summary
  print(results_summary) 
}
