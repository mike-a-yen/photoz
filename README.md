## Use EAZY and BPZ to yield PDF(z) plots

First get images from various filters(seechange). Apply aperature photometry to the objects detected by sep (Kyle Barbary)
Write the AB magnitudes for each filter to a catalog. Catalog needs to include the object ID number (arbitrary), magnitude, and magnitude error in each filter. I have also included the x,y and RA,DEC.

1. Setup seechange: (this is where the majority of images come from)
    ```
    https://github.com/scplbl/seechange
    ```

2. Download EAZY:
    ```
    wget http://www.astro.yale.edu/eazy/download/eazy-1.00.tar.gz
    ```

3. Download BPZ:
    ```
    wget http://www.stsci.edu/~dcoe/BPZ/bpz-1.99.3.tar.gz
    ```

4. Setup current working directory with BPZ and EAZY tools:
    ```
    # $BPZPATH = path to BPZ
    # $EAZYPATH = path to EAZY
    ln -s $BPZPATH/bpz.py BPZSRC
    ln -s $EAZYPATH/filters/FILTER.RES.latest FILTER.RES.latest
    ln -s $EAZYPATH/filters/FILTER.RES.latest.info FILTER.RES.latest.info
    ln -s $EAZYPATH/templates templates
    ln -s $EAZYPATH/src/eazy EAZYSRC
    ```

5. With all of this setup check the BPZ and EAZY test cases
    ```
    $BPZPATH/test
    $EAZYPATH/doc/PythonDemo/
    ```
6. Now you need to ingest some hst images into the seechange_db
    Ask Kyle for instructions

7. Access seechange database
    from shell, run
    ```
    see_change_db --data_path data/
    ```
    this opens an ipython prompt
    
8. Initialize cluster
    ```
    import makeCat as MC
    cl = MC.ClusterData('cluster name',detectionFilter='F105W')
    ```
    Default detectionFilter is F105W.
    
    Change if you want to use a different filter to detect objects.
    
    Can change filter at any time using
    ```
    cl.SetDetectionFilter('filter')
    ```

  1. If the cluster hasn't been analyzed yet run (recommended)
    ```
    cl.FromScratch()
    ```
    
    This will run all of the necessary commands to yield P(z) plots.
    
    **Continue to step 3**
    
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
  7. Plot SED fits from BPZ and EAZY (only have EAZY right now)
      ```
      cl.PlotEAZYSED()
      ```
9. With all the plots and catalogs made, it is useful to be able to get the P(z)
    from the closest object to a given RA and DEC.
    You specify an RA and DEC, get the P(z) from the object closest to that RA and DEC
    ```
    objInfo, separation = cl.MatchRADECtoOBJ(ra,dec,threshold=10.) # ra and dec should be unitless in degree format
    # objInfo contains all information in the catalog for this object
    # separation in arcsec from the given ra and dec
    ```
    
10. To view the plots
    ```
    cd pzplots
    display ID#_clusterName_typeOfPlot.png
    ```

* In addition to objects detected by sep (Kyle Barbary), you can specify your own RA, DEC, and aperature size, and perform photo z measurements on these new 'objects'. This might be useful if something in the image is not detected by sep, or you know you want a photo z measurement of a particular RA and DEC.
  * For the aperture size you need to specify a major (a) and minor (b) axis as well as an orientation angle (theta).
  * Arguments for ra, dec, a, b, theta can be number types, lists, or np.ndarrays. They must be of the same type and have the same length.
 
1. You can add your aperture when initializing ClusterData.
  ```
  cl = ClusterData(ra=...,dec=...,a=...,b=...,theta=...)
  cl.WriteCatalog() # write the aperture info to catalog
  ```
  
2. You can also add an aperture on the fly.
  ```
  cl.ManualApertureAdd(ra=...,dec=...,a=...,b=...,theta=...,apFile=None)
  # apFile is a path to a file with aperture params ra dec a b theta
  # the apertures in this file will be appended to the apertures you specify by hand
  # This will automatically run ManApFillObjInfo() and update cl.objInfo
  cl.WriteCatalog() # you still need to update the catalog explicitly
  ```
  
3. To remove your apertures
  ```
  cl.ManualApertureRemove(ra=...,dec=...)
  # Removes the aperture at specified RA and DEC and deletes aperture info from objInfo
  cl.ManualApertureRemove(purge=True)
  # Removes all apertures and deletes aperture info from objInfo
  cl.WriteCatalog() # you still need to update the catalog explicitly
   ```
4. To save your apertures in a text file for later use
  ```
  cl.ManualAperturesSave('path/to/aperture_file')
  ```
