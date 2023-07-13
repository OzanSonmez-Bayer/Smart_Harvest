# Dynamic way of checking required libraries in R
list.of.packages <- c("tidyverse", 
                      "readxl", 
                      "kableExtra",
                      "DT",
                      "plotly",
                      "tidyr",
                      "shiny",
                      "kableExtra",
                      "flexdashboard",
                      "lme4",
                      "stringr",
                      "progress",
                      "parallel"
)


dynamic_require <- function(package){
  if(eval(parse(text=paste("require(",package,")")))){return(TRUE)}
  install.packages(package)
  return(eval(parse(text=paste("require(",package,")"))))
}
sapply(list.of.packages, dynamic_require)

# make sure rename function is called from dplyr
rename = dplyr::rename
filter = dplyr::filter
select = dplyr::select
