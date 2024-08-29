# dynamic_niches

This code accompanies the paper "Migratory birds modulate niche tradeoffs in rhythm with seasons and life history" by Yanco et al. PNAS http://doi.org/10.1073/pnas.2316827121

The main workflow for this analysis can be run by following the steps oultined in the file `wf_dynamic_niches.sh`.  This produces all anlayses and then primary results plots.  The remaining ancilary plots can be found in the `./src/plots` directory.  The data used in this analysis are publically available (see below) but have not been assembled in the SQLite db format used in this analysis.  Please see https://benscarlson.github.io/mosey for instructions on how to generate a database.  With minimal manipulations, the code can also be made to run on raw CSV files.

Environmental annotations rely on the Mosey system - code for that can be found in `src/mosey` and more details can be found at https://benscarlson.github.io/mosey An associated (but nut strictly necessary) conda environment is included in `src/conda` to assist with packages for annotation.  Necessary packages can also be installed independently (and may be an easier choice, especially if you don't already have conda installed and initialized).

`src/funs` and `src/init` contian convenience code called by other scripts.

## Data

The data used in this study can be found at the following locations:

Yuriy Andryushchenko, Ivan Pokrovsky, Bernd Vorneweg, Wolfgang Fiedler, Martin Wikelski. 2024. Data from: Study "1000 Cranes - Ukraine" Movebank Data Repository. https://doi.org/10.5441/001/1.590

Elena I. Ilyashenko, Ivan Pokrovsky, Andrey E. Gavrilov, Syrymgul Zaripova, Wolfgang Fiedler, Martin Wikelski. 2024. Data from: Study "1000 Cranes. Southern  Kazakhstan." Movebank Data Repository. https://doi.org/10.5441/001/1.591

Elena I. Ilyashenko, Ivan Pokrovsky, Oleg A. Goroshko, Elena A. Mudrik, Dmitry Politov, Valentin Yu. Ilyashenko, Wolfgang Fiedler, Martin Wikelski. 2024. Data from: Study "1000 Cranes. Russia. Transbaikalia." Movebank Data Repository. https://doi.org/10.5441/001/1.592

Elena I. Ilyashenko, Ivan Pokrovsky, Valentin Yu. Ilyashenko, Mikhail Korepov, Wolfgang Fiedler, Martin Wikelski. 2024. Data from: Study "1000 Cranes. Russia. Common Crane" Movebank Data Repository. https://doi.org/10.5441/001/1.593

Elena I. Ilyashenko, Ivan Pokrovsky, Valentin Yu. Ilyashenko, Elena A. Mudrik, Mikhail Korepov, Mnatsekanov R., Dmitry Politov, Wolfgang Fiedler, Martin Wikelski. 2024. Data from: Study "1000 Cranes. Russia. Taman. Azov." Movebank Data Repository. https://doi.org/10.5441/001/1.594

Elena I. Ilyashenko, Ivan Pokrovsky, Valentin Yu. Ilyashenko, Elena A. Mudrik, Mikhail Korepov, Dmitry Politov, Elena Gugueva, Wolfgang Fiedler, Martin Wikelski. 2024. Data from: Study "1000 Cranes. Russia. Volga-Ural Interfluve." Movebank Data Repository. https://doi.org/10.5441/001/1.595

Elena I. Ilyashenko, Ivan Pokrovsky, Elena A. Mudrik, Dmitry Politov, Kirill Postelnykh, Valentin Yu. Ilyashenko, Wolfgang Fiedler, Martin Wikelski. 2024. Data from: Study "1000 Cranes. Russia. Altai." Movebank Data Repository. https://doi.org/10.5441/001/1.596

Mindaugas Dagys & Ramūnas Žydelis. 2024. Data from: Study "Common Crane Lithuania GPS, 2015-2016" Movebank Data Repository. https://doi.org/10.5441/001/1.604

Mindaugas Dagys & Ramūnas Žydelis. 2024. Data from: Study "Common Crane Lithuania GPS, 2016" Movebank Data Repository. https://doi.org/merged

Ramūnas Žydelis, Mark Desholm, Johan Månsson, Lovisa Nilsson, Henrik Skov. 2024. Data from: Study "GPS telemetry of Common Cranes, Sweden" Movebank Data Repository. https://doi.org/10.5441/001/1.597

Nyambayar Batbayar, Batbayar Galtbalt, Tseveenmyadag Natsagdorj, Tuvshintugs Sukhbaatar, Martin Wikelski. 2024. Data from: Study "1000 Cranes. Mongolia" Movebank Data Repository. https://doi.org/10.5441/001/1.598

Nyambayar Batbayar, Batbayar Galtbalt, Tseveenmyadag Natsagdorj, Tuvshintugs Sukhbaatar, Martin Wikelski. 2024. Data from: Study "LifeTrack Mongolia Demoiselle cranes" Movebank Data Repository. https://doi.org/10.5441/001/1.599

Nyambayar Batbayar, Batbayar Galtbalt, Tseveenmyadag Natsagdorj, Tuvshintugs Sukhbaatar, Martin Wikelski. 2024. Data from: Study "White-naped crane Mongolia WSCC" Movebank Data Repository. https://doi.org/10.5441/001/1.600

Sherub Sherub, Martin Wikelski. 2024. Data from: Study "Grus nigricollis - BHUTAN-MPIAB  GSM" Movebank Data Repository. https://doi.org/10.5441/001/1.601

## Questions?

If you have any questions about this code base or the paper itself, please email:  
Scott Yanco. 
syanco@umich.edu
