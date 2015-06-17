## Use EAZY and BPZ to yield PDF(z) plots

Setup seechange: (this is where the majority of images come from)
    ```
    https://github.com/scplbl/seechange
    ```
Download EAZY:
    ```
    wget http://www.astro.yale.edu/eazy/download/eazy-1.00.tar.gz
    ```
Download BPZ:
    ```
    wget http://www.stsci.edu/~dcoe/BPZ/bpz-1.99.3.tar.gz
    ```
Setup current working directory with BPZ and EAZY tools: ($BPZPATH = path to BPZ, $EAZYPATH = path to EAZY)
    ```
    ln -s $BPZPATH/bpz.py BPZSRC
    ln -s $EAZYPATH/filters/FILTER.RES.latest FILTER.RES.latest
    ln -s $EAZYPATH/filters/FILTER.RES.latest.info FILTER.RES.latest.info
    ln -s $EAZYPATH/templates templates
    ln -s $EAZYPATH/src/eazy EAZYSRC
    ```
    
With all of this setup check the BPZ and EAZY test cases
    $BPZPATH/test
    $EAZYPATH/doc/PythonDemo/
    
Now you need to ingest some hst images into the seechange_db
    Ask Kyle for instructions

First get images from various filters(seechange)
Apply aperature photometry to the objects detected by sep (Kyle Barbary)
Write the AB magnitudes for each filter to a catalog
    Catalog needs to include the object ID number (arbitrary), magnitude, and magnitude error in each filter 
    I have also included the x,y and RA,DEC.

1. Access seechange database
    from shell, run
    ```
    see_change_db --data_path data/
    ```
    this opens an ipython prompt
    
2. Initialize cluster
    ```
    import makeCat as MC
    cl = MC.ClusterData('cluster name',detectionFilter='F105W')
    ```
    Default detectionFilter is F105W
    Change if you want to use a different filter to detect objects
    Can change filter at any time using cl.SetDetectionFilter('filter')

  1. If the cluster hasn't been analyzed yet run (recommended)
    ```
    cl.FromScratch()
    ```
    This will run all of the necessary commands to yield P(z) plots
    Continue to step 3
  2. Check requirements to run properly, this will create all required but non existant directories
        ```
        cl.CheckRequiredDirectories()
        cl.CheckDetectionFilter()
        ```
  3. Write all required files to run BPZ and EAZY
        ```
        cl.WriteCatalog() # for BPZ and EAZY
        cl.WriteColumns()
        cl.ParamBPZ()
        cl.CheckBPZRequirements() # make sure everything is ready to run BPZ
        cl.ParamEAZY(ASCII=False)
        cl.TranslateEAZY()
        cl.ZeroPointEAZY()
        cl.CheckEAZYRequirements() # make sure everything is ready to run EAZY
        ```
  4. Run BPZ and EAZY
        ```
        cl.RunBPZ() # requires symlink setup
        cl.RunEAZY() # requires symlink setup
        ```
  5. Turn the BPZ and EAZY output into a more useable format
        Saves outputs as a pickled (binary or ASCII) dictionary
        {'z':zrange, 'obID':P(z), ...}
        ```
        cl.ProcessBPZOutput(ASCII=False)
        cl.ProcessEAZYOutput(ACII=False)
        # the best P(z) come from combining the BPZ and EAZY
        cl.CombinePZ(ASCII=False)
        ```
  6. Plot all P(z)'s
        ```
        cl.PlotPickledPZ('BPZ')
        cl.PlotPickledPZ('EAZY')
        cl.PlotPickledPZ('COMBINED')
        cl.PlotAllMethods()
        ```

3. With all the plots and catalogs made, it is useful to be able to get the P(z)
    from the closest object to a given RA and DEC.
    You specify an RA and DEC, get the P(z) from the object closest to that RA and DEC
    ```
    objInfo, sep = cl.MatchRADECtoOBJ(ra,dec) # ra and dec should be unitless in degree format
    # objInfo contains all information in the catalog for this object
    # sep is the separation in arcsec from the given ra and dec
    ```
    
4. To view the plots
    ```
    cd pzplots
    display ID#_clusterName_typeOfPlot.png
    ```
    
