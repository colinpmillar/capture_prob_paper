
{

library(CLdata)

if (Sys.info()["user"] == "millaco") {
  setwd("~/work/SMFS-report")    
} else 
if (Sys.info()["user"] == "millarc") {
  setwd("C:/work/repos/papers/capture_prop_paper/")
} else 
if (Sys.info()["user"] == "Millarc") {
  setwd("C:/work/SMFS-report")
}


# subset data to remove low Runs and no area.
data(ef)
ef <- subset(ef, !is.na(Site_OBJECTID) & !is.na(Area) & Runs >= 2)

# restructure ef
ef <- do.call(rbind, 
        lapply(c("S0", "SP", "T0", "TP"), 
            function(x) {
              out <- ef[c("Site_OBJECTID", "Site.Name", "Dataset", "Width", "Date", "Runs", "Area", "Trust", paste0(x, "_R", 1:6))]
              names(out) <- gsub(x, "n", names(out))
              out $ Species <- if (substring(x, 1, 1) == "S") "Salmon" else "Trout"
              out $ LifeStage <- if (substring(x, 2, 2) == "0") "Fry" else "Parr"
              out
            }
        ))

# drop rows with all NAs
ef <- ef[apply(ef[paste0("n_R", 1:6)], 1, function(x) !all(is.na(x))), ]


# tag on HMA data
data(hma)
hma <- hma[!(hma $ HAName %in% c("Shetlands", "Orkneys")),]
hma $ hmidx <- 1:nrow(hma)
gis <- CLdata::gis
gis @ data <- cbind(gis @ data, sp::over(gis, hma))

# tag on gis data
ef <- cbind(ef, gis[ef $ Site_OBJECTID,])
# fix missing value
ef $ n_R3[is.na(ef $ n_R3) & ef $ Runs == 3] <- 0
# add on data summaries
ef <- cbind(ef, CLmodel::getData(ef, passnames = paste0("n_R", 1:6)))
# add in some dates
ef $ pDate <- as.POSIXlt(ef $ Date, tz = "GMT", format = "%d/%m/%Y")
# try using lubridate
ef $ year <- ef $ pDate $ year + 1900
ef $ doy <- ef $ pDate $ yday

ef $ Trust <- as.character(ef $ Trust)
ef $ Trust[ef $ Trust == "Youngson"] <- "Other"
ef $ Trust[ef $ Trust == "Walker"] <- "MSS"

trustnames <- read.csv("trust_names.csv", stringsAsFactors = FALSE)
rownames(trustnames) <- trustnames $ code
bool <- ef $ Trust %in% c("MSS", "Other", "SEPA")
ef $ Trust[!bool] <- trustnames[ef $ Trust[!bool],"Trust"]
ef $ Trust <- factor(ef $ Trust, levels = unique(ef $ Trust))

ef $ CATCH_ID <- factor(ef $ CATCH_ID)

ef $ fyear <- factor(ef $ year)
ef $ sinSlope <- sin(ef $ Slope_deg/180*pi)

landuse <- c("CTrees", "Urban", "NCTrees", "Mixed", "Marsh", "Other")
ef $ totlanduse <- rowSums(ef[landuse])

ef[landuse] <- ef[landuse] / ef $ totlanduse

ef $ Distance_s <- ef $ Distance_s / 1000

# trim data
ef $ keep <- 
        with(ef, doy > 150 & doy < 325 & 
                 year >= 1997 & year <= 2013 & 
                 Area < 5000 &
                 Water_W < 35 & 
                 Elevation_ < 500 &  
                 sinSlope < 0.4 & 
                 Upcatch_km < 600 # &
                 #!barrier
            )

save(ef, file = "rData/modelData.rData")
}

# some data summaries

length(unique(ef $ Site_OBJECTID))
length(unique(paste(ef $ Site_OBJECTID, ef$Date)))


length(unique(ef $ DESCRIPTIO))
length(unique(ef $ CATCH_ID))

table(ef $ year)

table(ef $ Runs)
table(ef $ Runs)/nrow(ef)*100