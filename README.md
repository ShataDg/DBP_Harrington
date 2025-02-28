# DBP_Mary_Harrington_Circadian

Circadian cycle analysis of bioluminescent intravital microscopy of mice brain (SNC region).

The Harrington Lab uses bioluminescent markers to study circadian rhythms in suprachiasmatic nuclei (SCN). These intravital images are challenging to segment because of the difficulty to align cells/regions between time-points since the marker/signal can remain off for hours. 

To analyze the movies you need to first run the SNC_Cell_Detection_tracking.py script inside FIJI, make sure you have Trackmate v.7.6.1 installed and runnning.

After, you should run the Cell_Tracks_processing notebook to organize the tracks and filter out the short ones, this notebook will export a csv file ready to use in CIRCADIA app.


https://docs.google.com/document/d/1Fl0RbnoJP__zSHtG2HTdJrppsOX7NIvlbibB7ktUJx4/edit
