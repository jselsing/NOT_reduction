

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

    def __init__(self, input_dir, output_dir):

        self.input_dir = input_dir
        self.output_dir = output_dir

        try:
            os.mkdir(self.input_dir)
        except:
            pass

        self.construct_datasets()


    def construct_datasets(self):
        # Construct datasets based on file categories and image sizes
        self.filelist = glob.glob(self.input_dir + "*.fits")
        self.bias_files = [ii for ii in self.filelist if "bias" in fits.open(ii)[0].header["OBJECT"]]
        self.flat_files = [ii for ii in self.filelist if "FLAT" in fits.open(ii)[0].header["OBJECT"]]
        self.science_files = [ii for ii in self.filelist if "SCIENCE" in fits.open(ii)[0].header["IMAGECAT"]]

        # Loop through science files and find unique filter names
        filters = []
        sizes, xbins, ybins = [], [], []
        for ii in self.science_files:
            filt = fits.open(ii)[0].header["FAFLTNM"]
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
            science_files = [kk for kk in self.science_files if ii in fits.open(kk)[0].header["FAFLTNM"]]
            flat_files = [kk for kk in self.flat_files if ii in fits.open(kk)[0].header["FAFLTNM"] and ll == fits.open(kk)[0].header["DETWIN1"] and jj == fits.open(kk)[0].header["DETXBIN"] and hh == fits.open(kk)[0].header["DETYBIN"]]
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
            bias_list[ii] = ccdproc.CCDData(data=fitsfile[1].data, meta=fitsfile[1].header, unit="adu")

        self.master_bias = ccdproc.combine(bias_list)
        self.master_bias.write(self.output_dir+"/"+self.filename+"_masterbias.fits", clobber=True)

        # minlim, maxlim = np.percentile(master_bias, (1, 99))
        # pl.imshow(master_bias, vmax=maxlim, vmin=minlim, cmap="viridis")
        # pl.show()

    def make_master_flat(self, list_of_flatfiles):
        flat_list = [0]*len(list_of_flatfiles)
        for ii, kk in enumerate(list_of_flatfiles):
            fitsfile = fits.open(kk)
            flat = ccdproc.CCDData(data=fitsfile[1].data, meta=fitsfile[1].header, unit="adu")
            flat_list[ii] = ccdproc.subtract_bias(flat, self.master_bias)

        self.master_flat = ccdproc.combine(flat_list)
        self.master_flat.write(self.output_dir+"/"+self.filename+"_masterflat.fits", clobber=True)

        # minlim, maxlim = np.percentile(self.master_flat, (1, 99))
        # pl.imshow(self.master_flat, vmax=maxlim, vmin=minlim, cmap="viridis")
        # pl.show()


    def make_science(self, list_of_sciencefiles):
        science_list = [0]*len(list_of_sciencefiles)
        for ii, kk in enumerate(list_of_sciencefiles):
            fitsfile = fits.open(kk)
            fitsfile_common, _ = reproject_interp(fitsfile, fits.open(list_of_sciencefiles[0])[1].header, hdu_in=1)
            science = ccdproc.CCDData(data=fitsfile_common, meta=fitsfile[1].header, unit="adu")
            science_list[ii] = ccdproc.ccd_process(science, master_bias=self.master_bias, master_flat=self.master_flat, gain_corrected=False)



        self.combined_science = ccdproc.combine(science_list, method="median")

        self.master_science = ccdproc.cosmicray_lacosmic(self.combined_science, sigclip=5)

        self.master_science.write(self.output_dir+"/"+self.filename+".fits", clobber=True)

        # minlim, maxlim = np.percentile(self.master_science, (14, 86))
        # pl.imshow(self.master_science, vmax=maxlim, vmin=minlim, cmap="viridis")
        # pl.show()

def main():

    input_dir = "/Users/jselsing/Work/work_rawDATA/NOT/GRBs/GRB1710101A/"
    outdir = input_dir + "reduced_data"




    datasets = NOTdataset(input_dir, outdir)

    datasets.make_reduction()




if __name__ == '__main__':
    main()