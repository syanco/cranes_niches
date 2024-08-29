#######################################
#--     Dynamic Niches Scratch      --#
#--           1000 Cranes           --#
#--         Scott Yanco, Phd        --#
#--       scott.yanco@yale.edu      --#
#######################################

# Workflow script for processing individual, time-dynamic
# Grinellian niches.  


#-----------------------------------------------------------------#

#set working directory
wd=~/projects/dynamic_niches

#move to WD
cd $wd


##-- Join annotations --##

# Activate conda env
# conda activate db

# Run annotation join script - results returned to db
# NOTE: control file must be updated to select which variables are joined
# Rscript $wd/src/workflow/join_annos.r $wd/data/anno_move.db $wd/data/anno_move.db $wd/ctfs/anno_vars.csv


#-- Clean Data --#

# Activate conda env
# conda activate niche_mix

# Run clean script
Rscript $wd/src/workflow/clean_movement.r $wd/data/anno_move.db event event_clean


#-- Add Water Annotation --#

# conda activate annotate

src=$wd/src
db=$wd/data/anno_move.db

export MOSEYENV_SRC=~/projects/dynamic_niches/src/mosey #for mosey_anno_gee.sh


#-- Create study.csv file since I don't have one
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
# $MOSEYENV_SRC/mosey_anno_hee_test.sh $geePtsP $gcsOutP 

# $MOSEYENV_SRC/mosey_anno_gee_hab.sh $geePtsP $gcsOutP 

# Import to db
sqlite3 $db "alter table event_clean add column dist2water REAL;"
sqlite3 $db "alter table event_clean add column EVI REAL;"
sqlite3 $db "alter table event_clean add column lst REAL;"
# sqlite3 $db "alter table event_clean add column hab TEXT;"
# sqlite3 $db "alter table event_clean add column pfor_h REAL;"

# $MOSEYENV_SRC/import_anno_hab.sh $gcsOutURL $annoP $db 
$MOSEYENV_SRC/import_anno.sh $gcsOutURL $annoP $db 


#---- Unpack ESA CCI
# UNpack the STOAT-annotated ESA-CCI data and then join up to the `event_clean` table in db, save back into db as `event_clean_new`

# Warning:  bad code, takes ~5 hrs...
Rscript $wd/src/workflow/unpack_habitat.r
Rscript $wd/src/workflow/add_unpacked.r


#---- Data Imputation and Prep

Rscript $wd/src/workflow/prep_data.r $wd/data/anno_move.db event_clean_new


#-- Analysis --#
########          DEPRECATED              #######
#################################################
# #---- Q1: Individual Consistency ax Years ----#
# 
# # Compare season_matched niches across years (within individuals)
# 
# # Activate conda env
# # conda activate niche_mix
# 
# # Run script
# Rscript $wd/src/workflow/within_ind_across_years.r $wd/data/anno_move.db niche_prep $wd/out/wi_ind_ax_yrs.csv
# 
# #---- Q2: Individual Consistency ax Seasons ----#
# 
# # Compare season_matched niches across years (within individuals)
# 
# # Activate conda env
# # conda activate niche_mix
# 
# # Run clean script
# Rscript $wd/src/workflow/within_ind_across_seasons.r $wd/data/anno_move.db event_clean $wd/out/wi_ind_ax_seas.csv