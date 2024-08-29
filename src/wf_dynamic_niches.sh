#######################################
#--     Dynamic Niches Scratch      --#
#--           1000 Cranes           --#
#--         Scott Yanco, Phd        --#
#--       scott.yanco@yale.edu      --#
#######################################

# Workflow script for cleaning and annotating GPS tracks


#-----------------------------------------------------------------#

#set working directory
wd=~/projects/dynamic_niches

#move to WD
cd $wd


#-- Clean Data --#

# Run clean script
Rscript $wd/src/workflow/clean_movement.r $wd/data/anno_move.db event event_clean


#-- Annotate Data with GEE --#
# Use GEE to annotate temp, evi, and water proximity

# conda activate annotate

src=$wd/src
db=$wd/data/anno_move.db

export MOSEYENV_SRC=~/projects/dynamic_niches/src/mosey #for mosey_anno_gee.sh


#-- Create study.csv file to act as ctf (ONLY RUN ONCE)
# 
# sql="select study_id, study_name, 1 as run
# from study
# where study_id in (select distinct study_id from event_clean)
# order by study_id"
# 
# sqlite3 -header -csv $db "$sql;" > ctfs/study.csv

#-- Set up variables

#GEE
geePtsP=users/syanco/cranes/tracks #folder holding the gee point datasets

#GCS
gcsBucket=1000cranes-bucket
gcsInP=ingest_gee #This holds the csvs that will be imported to gee
gcsOutP=annotated #This is the output folder for annotated csvs (excluding bucket)
gcsInURL=gs://${gcsBucket}/${gcsInP} #This is the url to the gee ingest folder
gcsOutURL=gs://${gcsBucket}/${gcsOutP} #This is the url to the output folder (includes bucket)

#LOACL
csvP=out/anno/individual-files #local folder that holds the csv files to be ingested into gee
annoP=out/anno/annotated #local folder that holds the annotated csv files
envP=ctfs/env.csv
sesid=full_wf

#---- Import studies into GEE

$MOSEYENV_SRC/gee_ingest.sh $sesid $geePtsP $gcsInURL $csvP --db $db


#---- Annotate 
$MOSEYENV_SRC/mosey_anno_gee.sh $geePtsP $gcsOutP 


# Import to db
sqlite3 $db "alter table event_clean add column dist2water REAL;"
sqlite3 $db "alter table event_clean add column EVI REAL;"
sqlite3 $db "alter table event_clean add column lst REAL;"


# $MOSEYENV_SRC/import_anno_hab.sh $gcsOutURL $annoP $db 
$MOSEYENV_SRC/import_anno.sh $gcsOutURL $annoP $db 


#---- Annotate Data with STOAT
# Use STOAT annotations to include ESA CCI
# Unpack the STOAT-annotated ESA-CCI data and then join up to the `event_clean` table in db, save back into db as `event_clean_new`

Rscript $wd/src/workflow/unpack_habitat.r
Rscript $wd/src/workflow/add_unpacked.r


#---- Data Prep

# This prepares the data for use in the PCA scripts.  Saves rsults to "niche_prep" table
Rscript $wd/src/workflow/prep_data.r $wd/data/anno_move.db event_clean_new


#-- Analysis --#
# (these scripts also produce main results visualizations)

# Time-ordered PCA (Niche position)
Rscript $wd/src/workflow/pca_time_position.r 

# Time-ordered PCA (Niche breadth)
Rscript $wd/src/workflow/pca_time_var.r 

# PCA Static (Niche position)
Rscript $wd/src/workflow/pca_all_sp_position.r 

# PCA Static (Niche breadth)
Rscript $wd/src/workflow/pca_all_sp_var.r 

# PCA Stats
Rscript $wd/src/workflow/prop_var.r 

#-- Maps and other plots --##
$wd/src/plots/map_fig.r # Makes track map for Figure 1 
$wd/src/plots/plot_scaled_vars.r # Makes example env tradeoff panel for Figure 1 
$wd/src/plots/raw_vars_plots.r # Makes Supplemental figures showing raw variable traces 

