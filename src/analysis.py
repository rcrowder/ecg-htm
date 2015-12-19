import sys
import datetime
import os
import math
import argparse
import pandas as pd
import numpy as np
from scipy.signal import butter, lfilter, freqz, find_peaks_cwt
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Set CSV data directory
DATA_DIR = "data/"


# Setup data and time parsing, including milliseconds
def date_parse(dates):
  return [pd.datetime.strptime(d, '%Y-%m-%d %H:%M:%S.%f') for d in dates]


def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a


def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y


def analyse_csv_file(file_name, grid_spec, column_num):
    print("DEBUG: Reading ", file_name)
  
    # https://www.sparkfun.com/products/12650 - Part #SEN-12650
    # The SparkFun heart monitor sensor is capturing at ~9600 baud rate

    # Read csv with no header row, using comma as delimiter and first column as index
    df = pd.read_csv(DATA_DIR+file_name, sep=',', header=None, index_col=0, parse_dates=[0], date_parser=date_parse)
  
    # Re-sample data to milliseconds
    df.resample('L', how='mean')
  
    # Grab some sample data
    samples = df[2500:7500][1]

    # Or all the samples
    #samples = df[:][1]

    samples -= samples.min()
    samples /= samples.max()
    samples -= samples.mean()

    # Filter requirements.
    order = 6
    fs = 30.0     # sample rate, Hz

    # Heart beats never exceed (220-age) beats per minute (60 seconds)
    cutoff_frequency = 220.0 / 60.0  # desired cutoff frequency of the filter, Hz

    # Get the filter coefficients so we can check its frequency response.
    b, a = butter_lowpass(cutoff_frequency, fs, order)

    # Plot the frequency response.
    w, h = freqz(b, a, worN=8000)
    ax0 = plt.subplot(grid_spec[0, column_num])
    ax0.plot(0.5*fs*w/np.pi, np.abs(h), 'b')
    ax0.plot(cutoff_frequency, 0.5*np.sqrt(2), 'ko')
    ax0.axvline(cutoff_frequency, color='k')
    ax0.set_xlim(0, 0.5*fs)
    ax0.set_title("Low-pass Filter Response")
    ax0.set_xlabel('Frequency [Hz]')
    ax0.grid()

    # Filter the data.
    filtered_samples = butter_lowpass_filter(samples, cutoff_frequency, fs, order)

    # Find peaks in sample data
    indexes = find_peaks_cwt(
        filtered_samples,
        np.arange(1, int(fs)),
        max_distances=np.arange(1, int(fs))*2)
    indexes = np.array(indexes) - 1

    ax1 = plt.subplot(grid_spec[1, column_num])
    ax1.set_title(file_name)
    ax1.plot(filtered_samples)
    ax1.plot(indexes, filtered_samples[indexes], '+', mfc=None, mec='r', mew=2, ms=8)
    ax1.grid()

    window_size = 64
    kaiser_filter = np.kaiser(window_size, 14)

    ax2 = plt.subplot(grid_spec[2, column_num])
    ax2.set_title("Specgram")
    ax2.specgram(filtered_samples, NFFT=window_size, noverlap=16, window=kaiser_filter)
  
    plt.draw()
    return


# Parse arguments
parser = argparse.ArgumentParser(description = "")
parser.add_argument('--source', required=True)
args = parser.parse_args()

plt.figure(1)
plt.grid()
# 3 graphs per column, and one column per CSV file analysed
gs = gridspec.GridSpec(3, 4)

analyse_csv_file("healthy_person1.csv", gs, 0)
analyse_csv_file("healthy_person2.csv", gs, 1)
analyse_csv_file("disease_person1.csv", gs, 2)
analyse_csv_file("disease_person2.csv", gs, 3)

plt.show()

raw_input("Finished. Press any key to exit.")
exit()


#------------------------------------------------------------------------------
class ECG:

    def trim(self, data, degree=100):
        print 'trimming'
        i, data2 = 0, []
        while i < len(data):
            data2.append(sum(data[i:i+degree])/degree)
            i += degree
        return data2

    def smooth(self, list, degree=15):
        mults = [1]
        s = []
        for i in range(degree):
            mults.append(mults[-1]+1)
        for i in range(degree):
            mults.append(mults[-1]-1)
        for i in range(len(list)-len(mults)):
            small = list[i:i+len(mults)]
            for j in range(len(small)):
                small[j] *= mults[j]
            val = sum(small) / sum(mults)
            s.append(val)
        return s

    def smoothWindow(self, list, degree=10):
        list2 = []
        for i in range(len(list)):
            list2.append(sum(list[i:i+degree])/float(degree))
        return list2

    def invertYs(self):
        print 'inverting'
        self.ys *= -1

    def takeDeriv(self, dist=5):
        print 'taking derivative'
        self.dys = []
        for i in range(dist,len(self.ys)):
            self.dys.append(self.ys[i]-self.ys[i-dist])
        self.dxs = self.xs[0:len(self.dys)]

    def genXs(self, length, hz):
        print 'generating Xs'
        step = 1.0 / hz
        xs = []
        for i in range(length):
            xs.append(step*i)
        return xs

    def loadFile(self, fname, startAt=None, length=None, hz=1000):
        print 'loading', fname
        self.ys = np.memmap(fname, dtype='h', mode='r') * -1
        print 'read %d points.'%len(self.ys)
        self.xs = self.genXs(len(self.ys),hz)
        if startAt and length:
            self.ys = self.ys[startAt:startAt+length]
            self.xs = self.xs[startAt:startAt+length]

    def findBeats(self):
        print 'finding beats'
        self.bx, self.by = [], []
        for i in range(100, len(self.ys)-100):
            if self.ys[i] < 15000:
                continue # SET THIS VISUALLY
            if self.ys[i] < self.ys[i+1] or self.ys[i] < self.ys[i-1]:
                continue
            if self.ys[i]-self.ys[i-100] > 5000 and self.ys[i]-self.ys[i+100] > 5000:
                self.bx.append(self.xs[i])
                self.by.append(self.ys[i])
        print "found %d beats" % (len(self.bx))

    def genRRIs(self, fromText=False):
        print 'generating RRIs'
        self.rris = []
        if fromText:
            mult = 1
        else:
            1000.0
        for i in range(1, len(self.bx)):
            rri = (self.bx[i]-self.bx[i-1]) * mult
            #if fromText==False and len(self.rris)>1:
                #if abs(rri-self.rris[-1])>rri/2.0: continue
            #print i, "%.03f\t%.03f\t%.2f"%(bx[i],rri,60.0/rri)
            self.rris.append(rri)

    def removeOutliers(self):
        beatT = []
        beatRRI = []
        beatBPM = []
        for i in range(1, len(self.rris)):
            if self.rris[i] < 0.5 or self.rris[i] > 1.1:
                continue #CHANGE THIS AS NEEDED
            if abs(self.rris[i]-self.rris[i-1]) > self.rris[i-1]/5:
                continue
            beatT.append(self.bx[i])
            beatRRI.append(self.rris[i])
        self.bx = beatT
        self.rris = beatRRI

    def graphTrace(self):
        plt.plot(self.xs,self.ys)
        #plt.plot(self.xs[100000:100000+4000],self.ys[100000:100000+4000])
        plt.title("Electrocardiograph")
        plt.xlabel("Time (seconds)")
        plt.ylabel("Potential (au)")

    def graphDeriv(self):
        plt.plot(self.dxs,self.dys)
        plt.xlabel("Time (seconds)")
        plt.ylabel("d/dt Potential (au)")

    def graphBeats(self):
        plt.plot(self.bx,self.by,'.')

    def graphRRIs(self):
        plt.plot(self.bx,self.rris,'.')
        plt.title("Beat Intervals")
        plt.xlabel("Beat Number")
        plt.ylabel("RRI (ms)")

    def graphHRs(self):
        #HR TREND
        hrs = (60.0 / np.array(self.rris)).tolist()
        bxs = (np.array(self.bx[0:len(hrs)])/60.0).tolist()
        plt.plot(bxs,hrs,'g',alpha=.2)
        hrs = self.smooth(hrs,10)
        bxs = bxs[10:len(hrs)+10]
        plt.plot(bxs,hrs,'b')
        plt.title("Heart Rate")
        plt.xlabel("Time (minutes)")
        plt.ylabel("HR (bpm)")

    def graphPoincare(self):
        #POINCARE PLOT
        plt.plot(self.rris[1:],self.rris[:-1],"b.",alpha=.5)
        plt.title("Poincare Plot")
        plt.ylabel("RRI[i] (sec)")
        plt.xlabel("RRI[i+1] (sec)")

    def graphFFT(self):
        #PSD ANALYSIS
        fft = np.fft.fft(np.array(self.rris)*1000.0)
        fftx = np.fft.fftfreq(len(self.rris), d=1)
        fftx, fft = fftx[1:len(fftx)/2], abs(fft[1:len(fft)/2])
        fft = self.smoothWindow(fft, 15)
        plt.plot(fftx[2:],fft[2:])
        plt.title("Raw Power Sprectrum")
        plt.ylabel("Power (ms^2)")
        plt.xlabel("Frequency (Hz)")

    def graphFFT2(self):
        #PSD ANALYSIS
        fft = np.fft.fft(np.array(self.rris)*1000.0)
        fftx = np.fft.fftfreq(len(self.rris), d=1)
        fftx, fft = fftx[1:len(fftx)/2], abs(fft[1:len(fft)/2])
        fft = self.smoothWindow(fft,15)
        for i in range(len(fft)):
            fft[i] = fft[i]*fftx[i]
        plt.plot(fftx[2:],fft[2:])
        plt.title("Power Sprectrum Density")
        plt.ylabel("Power (ms^2)/Hz")
        plt.xlabel("Frequency (Hz)")

    def graphHisto(self):
        plt.hist(self.rris,bins=20,ec='none')
        plt.title("RRI Deviation Histogram")
        plt.ylabel("Frequency (count)")
        plt.xlabel("RRI (ms)")
        #pdf, bins, patches = plt.hist(self.rris,bins=100,alpha=0)
        #plt.plot(bins[1:],pdf,'g.')
        #y=self.smooth(list(pdf[1:]),10)
        #x=bins[10:len(y)+10]
        #plt.plot(x,y)

    def saveBeats(self,fname):
        print "writing to",fname
        np.save(fname,[np.array(self.bx)])
        print "COMPLETE"

    def loadBeats(self,fname):
        print "loading data from",fname
        self.bx = np.load(fname)[0]
        print "loadded",len(self.bx),"beats"
        self.genRRIs(True)


def snd2txt(fname):
    ## SND TO TXT ##
    a = ECG()
    a.loadFile(fname)#,100000,4000)
    a.invertYs()
    plt.figure(figsize=(7,4), dpi=100)
    plt.grid(alpha=.2)
    a.graphTrace()
    a.findBeats()
    a.graphBeats()
    a.saveBeats(fname)
    plt.show()

def txt2graphs(fname):
    ## GRAPH TXT ##
    a = ECG()
    a.loadBeats(fname+'.npy')
    a.removeOutliers()
    plt.figure(figsize=(7,4), dpi=100)
    plt.grid(alpha=.2)
    a.graphHRs()
    plt.subplots_adjust(left=.1, bottom=.12, right=.96)
    plt.savefig("DIY_ECG_Heart_Rate_Over_Time.png")
    plt.figure(figsize=(7,4), dpi=100)
    plt.grid(alpha=.2)
    a.graphFFT()
    plt.subplots_adjust(left=.13, bottom=.12, right=.96)
    plt.savefig("DIY_ECG_Power_Spectrum_Raw.png")
    plt.figure(figsize=(7,4), dpi=100)
    plt.grid(alpha=.2)
    a.graphFFT2()
    plt.subplots_adjust(left=.13, bottom=.12, right=.96)
    plt.savefig("DIY_ECG_Power_Spectrum_Weighted.png")
    plt.figure(figsize=(7,4), dpi=100)
    plt.grid(alpha=.2)
    a.graphPoincare()
    plt.subplots_adjust(left=.1, bottom=.12, right=.96)
    plt.savefig("DIY_ECG_Poincare_Plot.png")
    plt.figure(figsize=(7,4), dpi=100)
    plt.grid(alpha=.2)
    a.graphRRIs()
    plt.subplots_adjust(left=.1, bottom=.12, right=.96)
    plt.savefig("DIY_ECG_RR_Beat_Interval.png")
    plt.figure(figsize=(7,4), dpi=100)
    plt.grid(alpha=.2)
    a.graphHisto()
    plt.subplots_adjust(left=.1, bottom=.12, right=.96)
    plt.savefig("DIY_ECG_RR_Deviation_Histogram.png")
    plt.show()

fname='publish_05_10min.snd'
#raw_input("\npress ENTER to analyze %s..."%(fname))
snd2txt(fname)
#raw_input("\npress ENTER to graph %s.npy..."%(fname))
txt2graphs(fname)
