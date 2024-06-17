library(dplyr)

# generate suggestions:
suggests <- 
  CTVsuggest::CTVsuggest(taskview = "FunctionalData", n = 40)
# get dl stats:
popularity <- cranlogs::cran_downloads(suggests$Packages, 
                                                from = Sys.Date() - 365, 
                                                to = Sys.Date()) |> 
  group_by(package) |> summarize(popularity = sum(count))
# combine & sort
suggests <- 
  merge(suggests, popularity, by.x = "Packages", by.y = "package") |> 
  dplyr::arrange(desc(round(FunctionalData, 2)), desc(popularity))


# read descriptions:
sapply(paste0("https://cran.r-project.org/web/packages/",suggests$Packages), browseURL)


# after modifications to FunctionalData.md: 


ctv::read.ctv("FunctionalData.md", cran = TRUE)
ctv::check_ctv_packages("FunctionalData.md")
ctv::ctv2html("FunctionalData.md")
browseURL("FunctionalData.html")
