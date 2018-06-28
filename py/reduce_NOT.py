

# Plotting
import matplotlib; matplotlib.use('TkAgg')
import matplotlib.pyplot as pl
import seaborn as sns; sns.set_style('ticks')

# Utils
import os
from astropy.io import fits
import astropy.units as u

import numpy as np
import pandas
import glob
import ccdproc
from reproject import reproject_exact, reproject_interp
import astroscrappy


class NOTdataset(object):

    def __init__(self, input_dir, output_dir, instrument):

        self.input_dir = input_dir
        self.output_dir = output_dir
        self.instrument = instrument

        try:
            os.mkdir(self.input_dir)
        except:
            pass

        self.construct_datasets()


    def construct_datasets(self):
        # Construct datasets based on file categories and image sizes
        self.filelist = glob.glob(self.input_dir + "*.fits")
        self.bias_files = [ii for ii in self.filelist if "bias" in fits.open(ii)[0].header["OBJECT"]]

        if self.instrument == "AlFOSC":
            self.flat_files = [ii for ii in self.filelist if "FLAT" in fits.open(ii)[0].header["OBJECT"]]
        elif self.instrument == "StanCam":
            self.flat_files = [ii for ii in self.filelist if "flat" in fits.open(ii)[0].header["OBJECT"]]

        self.science_files = [ii for ii in self.filelist if "SCIENCE" in fits.open(ii)[0].header["IMAGECAT"]]

        # Loop through science files and find unique filter names
        filters = []
        sizes, xbins, ybins = [], [], []
        for ii in self.science_files:
            if self.instrument == "AlFOSC":
                filt = fits.open(ii)[0].header["FAFLTNM"]
            elif self.instrument == "StanCam":
                filt = fits.open(ii)[0].header["STFLTNM"]

            size = fits.open(ii)[0].header["DETWIN1"]
            xbin = fits.open(ii)[0].header["DETXBIN"]
            ybin = fits.open(ii)[0].header["DETYBIN"]
            if filt not in filters:
                filters.append(filt)
                sizes.append(size)
                xbins.append(xbin)
                ybins.append(ybin)

        self.filters = filters
        self.sizes = sizes
        self.xbins = xbins
        self.ybins = ybins

        self.datasets = {}
        for ii, ll, jj, hh in list(zip(self.filters, self.sizes, self.xbins, self.ybins)):

            if self.instrument == "AlFOSC":
                science_files = [kk for kk in self.science_files if ii in fits.open(kk)[0].header["FAFLTNM"]]
                flat_files = [kk for kk in self.flat_files if ii in fits.open(kk)[0].header["FAFLTNM"] and ll == fits.open(kk)[0].header["DETWIN1"] and jj == fits.open(kk)[0].header["DETXBIN"] and hh == fits.open(kk)[0].header["DETYBIN"]]
                bias_files = [kk for kk in self.bias_files if ll == fits.open(kk)[0].header["DETWIN1"] and jj == fits.open(kk)[0].header["DETXBIN"] and hh == fits.open(kk)[0].header["DETYBIN"]]
            elif self.instrument == "StanCam":
                science_files = [kk for kk in self.science_files if ii in fits.open(kk)[0].header["STFLTNM"]]
                flat_files = [kk for kk in self.flat_files if ii in fits.open(kk)[0].header["STFLTNM"] and ll == fits.open(kk)[0].header["DETWIN1"] and jj == fits.open(kk)[0].header["DETXBIN"] and hh == fits.open(kk)[0].header["DETYBIN"]]
                bias_files = [kk for kk in self.bias_files if ll == fits.open(kk)[0].header["DETWIN1"] and jj == fits.open(kk)[0].header["DETXBIN"] and hh == fits.open(kk)[0].header["DETYBIN"]]


            self.datasets[ii] = [science_files, flat_files, bias_files]


    def make_reduction(self):


        for ii in self.filters:
            self.filename = ii.split()[0]

            dataset = self.datasets[ii]
            self.make_master_bias(dataset[2])
            self.make_master_flat(dataset[1])
            self.make_science(dataset[0])
            # final_image = self.combine_images()



    def make_master_bias(self, list_of_biasfiles):
        bias_list = [0]*len(list_of_biasfiles)
        for ii, kk in enumerate(list_of_biasfiles):
            fitsfile = fits.open(kk)
            bias = ccdproc.CCDData(data=fitsfile[1].data, meta=fitsfile[1].header, unit="adu")
            bias_list[ii] = ccdproc.subtract_overscan(bias, fits_section='[5:35, :]')

        self.master_bias = ccdproc.combine(bias_list, method="median")
        self.master_bias.write(self.output_dir+"/"+self.filename+"_masterbias.fits", overwrite=True)

        # print(np.percentile(self.master_bias, (1, 50, 99)))
        # pl.imshow(master_bias, vmax=maxlim, vmin=minlim, cmap="viridis")
        # pl.show()

    def make_master_flat(self, list_of_flatfiles):
        flat_list = [0]*len(list_of_flatfiles)
        for ii, kk in enumerate(list_of_flatfiles):
            fitsfile = fits.open(kk)
            flat = ccdproc.CCDData(data=fitsfile[1].data, meta=fitsfile[1].header, unit="adu")
            flat_scan = ccdproc.subtract_overscan(flat, fits_section='[5:35, :]')
            flat_bias = ccdproc.subtract_bias(flat_scan, self.master_bias)
            flat_bias.data /= np.median(flat_bias.data[100:-100, 100:-100])
            flat_list[ii] = flat_bias

        self.master_flat = ccdproc.combine(flat_list, method="average", sigma_clip=True)


        self.master_flat.write(self.output_dir+"/"+self.filename+"_masterflat.fits", overwrite=True)

        # minlim, maxlim = np.percentile(self.master_flat, (1, 99))
        # pl.imshow(self.master_flat, vmax=maxlim, vmin=minlim, cmap="viridis")
        # pl.show()


    def make_science(self, list_of_sciencefiles):


        science_list = [0]*len(list_of_sciencefiles)
        science_names = [0]*len(list_of_sciencefiles)
        for ii, kk in enumerate(list_of_sciencefiles):
            fitsfile = fits.open(kk)
            science = ccdproc.CCDData(data=fitsfile[1].data, meta=fitsfile[1].header, unit="adu")
            # science_cos = ccdproc.cosmicray_lacosmic(science, verbose = True)
            overscan_sub = ccdproc.subtract_overscan(science, fits_section='[5:35, :]')
            bias_sub = ccdproc.subtract_bias(overscan_sub, self.master_bias)
            flat_corr = ccdproc.flat_correct(bias_sub, self.master_flat)

            # science_list[ii] = flat_corr
            flat_corr.write(self.output_dir+"/"+self.filename+"_"+str(ii)+".fits", overwrite=True)
            science_names[ii] = self.output_dir+"/"+self.filename+"_"+str(ii)+".fits"



        # science_list = [0]*len(list_of_sciencefiles)
        coverages = [0]*len(list_of_sciencefiles)
        for ii, kk in enumerate(science_names):
            fitsfile = fits.open(kk)
            fitsfile_common, coverage = reproject_interp(fitsfile, fits.open(science_names[0])[0].header, hdu_in=0)
            science = ccdproc.CCDData(data=fitsfile_common, meta=fitsfile[0].header, unit="adu")
            science_list[ii] = science
            coverages[ii] = ccdproc.CCDData(data=coverage, meta=fitsfile[0].header, unit="adu")

        self.combined_science = ccdproc.combine(science_list, method="median")
        # self.combined_coverage = ccdproc.combine(coverages, method="sum")


        header0 = fits.open(list_of_sciencefiles[0])[0].header
        for key in header0:
            # print(str(key), header0[str(key)])
            if "COMMENT" in key:
                continue
            self.combined_science.header[str(key)] = header0[str(key)]

        self.combined_science.write(self.output_dir+"/"+self.filename+".fits", overwrite=True)
        # minlim, maxlim = np.percentile(self.master_science, (14, 86))
        # pl.imshow(self.master_science, vmax=maxlim, vmin=minlim, cmap="viridis")
        # pl.show()

def main():

    input_dir = "/Users/jselsing/Work/work_rawDATA/NOT/GRBs/GRB171010A/"
    outdir = input_dir + "reduced_data"

    instrument = "AlFOSC"


    datasets = NOTdataset(input_dir, outdir, instrument)

    datasets.make_reduction()




if __name__ == '__main__':
    main()