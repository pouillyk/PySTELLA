# PySTELLA
A Python tool for SpecTral Emission Lines (variabiLity) Analysis

/!\ README under construction for a better shape /!\

This package allows you to perform a variability analysis of stellar emission lines (but can be adapted for absorption lines as well).
It was developed from ESPaDOnS observations, but it should work using any instrument.
It was designed for a terminal command launch, and text editor to fill the input files.
Please cite Pouilly et al. 2024 (submitted) if you use it for publication, and add the following to the acknoledgment:
The emission line variability was perfomed use the PySTELLA tool, available at https://github.com/pouillyk/PySTELLA

The following library are needed:
- numpy
- matplotlib
- PyAstronomy
- scipy

## Input files
First carrefully fill the `In/stellarParams.dat` and `In/infoSpec.dat` that contain all information allowing the code to run properly.

### stellarParams.dat
Enter the name of the star your are studying, without space (for a design purpose in output file headers), `y` or `n` if it is an SB2 system or not, and if yes please provide the luminosity ratio (`LR`=LA/LB)
Fill the parameters of star(s) (rotation period, vsini, inclination to the rotation axis, veiling, and same for secondary if SB2 + orbital period)
The number of observations.
Then enter the path of the observation file, the associated heliocentric julian date, the radial velocity, and if SB2 the radial velocity of the secondary. One line per observation.
Enter an origin of time (for phase folded plots and phase computation).
Enter the similar informations for the photospheric template (used for residual line computation). If SB2, you can provide different information on the two separated columns if you want to use different template for each component. A template with a lower (or equal) vsini than the studied star is mandatory.

### infoSpec.dat
For each star (Object of interest + Template 1+2), enter the information about the shape of the observation files. This means the number of header lines, the unit of the wavelength (if angstrom, it will be converted to nm), the column where it is located each parameter (0 being the first column). Only `wl`, `I` and `sigI` are needed for now, the others are there for a futur update including polarimetry.

## Produce a line profile file
### makeProfile.py
Now you can produce your first line profile file using `makeProfile.py`
In a terminal, run it using the following command:
```
python makeProfile.py wavelength outFile
```
The wavelength in nanometer and outFile without '' nor extension. Example:
```
python makeProfile.py 656.279 ha
```
This command will produce the `Out/Data/ha.out` file, containing the Ha lines of all observations, corrected from there radial velocity. Additional keywords are availables for different options:

 - `res=1` : computed the residual line, meaning it substracts the rotationally broadened and veiling-corrected photospheric template. It saves the Out/Data/outFile_res.out file (_res.out added automatically). With this mode, a figure will appear displaying the first line profile with the photospheric template as it will be substracted. This allows you to check that the various parameters and normalisation are correct. Once the figure close type `y` in the terminal to continue the procedure, `n` will quit to let you perform the modifications. In SB2 mode the template of both component are displayed, as well as the composition of the template (respecting the provided `LR`) that will be substracted.
 - With `res=1`, you can also add `inFile='outFile.out'`. If you already produced an outFile.out file, you can use it to compute the outFile_res.out file. It is a bit quicker to run and usefull if you modify the `outFile.out` file (as example, a re-normalisation, or a spike correction)
 - `corrVel=0` : produce an Out/Data/ouFile_unCorrVel.out file, containing the line profile not corrected from the radial velocity. WARNING: `corrVel=0` is not compatible with `res=1`.
 - `velMin=`/`velMax=` : to choose de velocity range on wich you want to compute the line profiles (default -500/+500 km/s)
 - `velStep=` : to change the velocity sampling (default 1.5 km/s)

In a terminal, enter the following command to get a reminder on its usage:
```
python makeProfile.py
```

### Continuum normalisation
If your line profile is not well normalised, ou can use the code `normLine.py`
```
python normLine.py outFile.out
```

For each profile, this opens a figure were you can double click on the continuum to select some continuum points (right click on it to delete the point). Then it fit a polynomial curve to these points to re-normalise the line profile. 
You can use the additional keyword `deg=` to adapt the polynome's degree.
The polynome is shown on the plot when you select more than 4 points.
At the end the normalised file is save in `Out/Data/norm_outFile.out`
To get a reminder on the usage:
```
python normLine.py
```
### Visualisation
Now you produced your line profiles file, you can visualise it using `plotProfilesSup.py` (superimposed visualisation) or `plotProfilesCol.py` (in columns visualisation)
```
python plotProfilesSup.py outFile.out
python plotProfilesCol.py outFile.out
```

For both you have the additional keywords `save=1` (to save the figure in `Out/Figure/outFileSup.pdf`, or `outFileCol.pdf`), and `showPlot=0` (to not show the plot, e.g. if you only want to save it)
For `plotProfilesCol.py`, you also have `colNB=` (the number of columns you want in the plot, default 2) and `profNB=` (the number of line profiles you want in each column, default 10).

To get a reminder on the usage:
```
python plotProfilesCol.py
```
## Analysis
### 2-Dimensional periodogram
You can also compute a 2D periodogram of the line, meaning and Lomb-Scargle periodogram of each velocity channel of the line using `periodo2d.py`.
```
python periodo2d.py outFile.out
```

You also have the `save` and `showPlot` keywords
It will show a plot with the frequency on y-axis, the line velocity on x-axis, and the periodogram power in color. A white-dotted line indicate the rotation period provided.
A lower panel shows the mean line profile and its variance.
A second figure (never saved, even with `save=1`) displays the corresonding False Alarm Probability map.

### Correlation matrices
Finally, you can computed the correlation matrix of lines with c`orrMatrix.py`. This consist of a computation of a linear correlation coefficient (here Pearson) between the velocity channels of two lines. This allows you to identify the various variability regions within the line and there correlation, and helps to identify the physical processes in place.
```
python corrMatrix.py outFileX.out outFileY.out
```
The `save` and `showPlot` keywords are still available.
This shows a figure with the velocity if `outFileX.out` on the x-axis, the velocity of `outFileY.out` on the y-axis, and the correlation coefficient in color.
The panels along the axis show the correponding mean line profile and it variance.
A second figure (never saved, even with `save=1`) displays only the correlation matrix, with contours highlighting the strongly correlated or anti-correlated regions.




