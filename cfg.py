version_name = "OzFluxQC"
version_number = "V2.9.6b"
# V2.9.6b  - bug fix of ustar implementation at L5
#            - gfSOLO was picking the target data from the L4 (not filtered)
#              data structure not the L5 (filtered) data structure so
#              filtered data in L5 was being overwritten by gap filled,
#              unfiltered data
#            - changed source of target data in gfSOLO_runsolo,
#              gfSOLO_runseqsolo and gfSOLO_plot from dsa to dsb
# V2.9.6a  - implementation of ustar filtering at L5 or L6
#            - previous versions applied the ustar filter at L6
#              after gap filling at L5 which meant the NN used
#              for gap filling at L5 was being trained on Fc
#              observations from periods when ustar was below the
#              threshold.
# V2.9.5   - implementation of new respiration options
#            - removed NN related code from qcrp.py and placed this
#              in a stand-alone module qcrpNN.py.
#            - implemented Ian McHugh's code for Lloyd-Taylor
# V2.9.4   - major bug fix
#            - a bug was introduced in V2.8.7 on 15/04/2015 that caused
#              Fg corrected for heat storage in the layer above the
#              ground heat flux plates to be replaced with uncorrected
#              Fg during the L3 processing
#            - this release fixes the bug
# V2.9.3   - updates and bug fixes
#            - implemented batch processing for L1 to L6 including climatology,
#              CPD, conatenation
#            - completed implementation of plot_path in control files
#            - fixed bug that caused the units of NEE, NEP, GPP and ER
#              in the L6 output file to be gC/m2
#            - fixed bug on gfalternate_matchstartendtimes
# V2.9.2   - updates and bug fixes
#            - implemented summary output to Excel file and plots at L6
#            - implemented "ols_thru0", "rma" and "odr" fit types at L4
#            - fixed bug in qcrp.GetERFromFc that let gap filled Fc data
#              through when estimating ecosystem respiration (ER) from
#              u*-filtered, nocturnal Fc
# V2.9.1   - hopefully completed the major re-write of the gap filling
#            routines for L4 and L5
#            - much testing and tweaking of gfalternate_autocomplete
#              to get it to run, the logic is rather tortuous at present
#              and needs a re-working, there is a gfalternate_autocomplete_rewrite
#              routine, just needs completion
#            - updated QCCPD code with Ian McHugh's latest version and added
#              code to trap empty results data frame before plotting histograms
# V2.9.0   - major re-write of gap filling routines to simplify workflow
#            - will document later
# V2.8.7   - fixed several bugs in the gap filling routine and improved
#            the gap filling workflow, implemented ability to split a
#            netCDF file at specified dates
# V2.8.6   - added a new switch "UseL2Fluxes" to the L3 processing:
#            - if true, skip calculating fluxes from covariances and skip corrections
#            - if false (default), use covariances as normal
# V2.8.5   - miscellaneous changes arising from use of V2.8.4 at the
#            2014 OzFlux Data Workshop
# V2.8.4   - changes as follows:
#            - split gap filling into L4 (meteorological drivers) and
#              L5 (fluxes), partitioning is now L6
#            - associated changes to the template control files
#            - implemented gap filling from BIOS2
#            - implemented "Import" at L4 to allow importing MODIS data
#              into OzFluxQC data path
# V2.8.3   - implemented estimation of u* threshold by CPD (Barr et al)
# V2.8.2   - implemented;
#            - estimation of ecosystem respiration from nocturnal Fc
#              using SOLO and FFNET.
#            - several bug fixes
# V2.8.1   - refactor of gfACCESS_plotdetailed and associated code
# V2.8.0   - implemented;
#            - gap filling using ACCESS data (works for any alternate site file)
#            - menu at top of OzFluxQC GUI
# V2.7.2   - several enhancements, mainly to do with gap filling
#            - implemented "interpolated daily" method of
#              gap filling from climatology
#            - implemented gap filling at L3
# V2.7.1   - fixed bug in CorrectFcForStorage, Fc_storage_in typo
# V2.7.0   - major bug fixes as for V2.6.3 above
#          - minor bug fixes to clean up use of switches in ['Options'] section
#          - minor fixes to check that requested files exist
#          - implemented compare_ep.py to automate comparison of OzFluxQC and
#            EddyPro results
# V2.6.3 - clean up of code after comparing EddyPro and OzFluxQC
#          - fix bugs in mf.vapourpressure and mf.RHfromabsolutehumidity
#          - deprecated WPLcov after finding an error in Fc_WPLcov (units of wA term)
#          - rationalised calculation of dry and moist air densities and parial
#            density of water vapour
#          - implemented EddyPro method of calculating rho*Cp
#          - tidyied up use of densities in WPL correction and used flux
#            form (WPL80 Eqn 42a and 44) for Fe and Fc
#          - implemented EddyPro method of calculating Fh from Fhv
#          - rationalised use of densities and rho*Cp when calculating fluxes
# V2.6.2 - implement Lloyd-Taylor ER
# V2.6.1 - fix 2D coordinate rotation for momentum covariances
# V2.6.0 - fix ConvertCO2 units bug
# V2.5.2 - implement sofm/solo/seqsolo
# V2.5.1 - post-Cairns 2013 version
