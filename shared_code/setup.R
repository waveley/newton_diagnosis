# #######
# 
# This file contains all our default knitr settings, libraries, and beeper. 
# Intended to be run before every program.
#
# ########

library(tidyverse)
library(progress) 
library(beepr)
library(gtsummary)
library(matrixcalc)


# set knitr defaults
knitr::opts_chunk$set(
  echo      = TRUE
  , fig.align = "center"
  , fig.width = 6
  , fig.asp   = .6
  , out.width = "90%"
)

# set theme defaults
theme_set(
  theme_bw() +
    theme(
      legend.position = "bottom"
      , plot.title    = element_text(hjust = 0.5)
      , plot.subtitle = element_text(hjust = 0.5)
      , plot.caption  = element_text(hjust = 0.0)
    )
)

# set color scale defaults
options(
  ggplot2.continuous.colour = "viridis",
  ggplot2.continuous.fill   = "viridis"
)
scale_colour_discrete = scale_colour_viridis_d
scale_fill_discrete   = scale_fill_viridis_d


########### Beeper #############

options(error = function(){    # Beep on error
  beepr::beep()
  Sys.sleep(1)
}
)

.Last <- function() {          # Beep on exiting session
  beepr::beep()
  Sys.sleep(1)
}