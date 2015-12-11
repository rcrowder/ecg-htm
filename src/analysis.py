import sys
import datetime
import os
import pandas as pd
import numpy
import matplotlib.pyplot as plt
from collections import deque

# setting for fft conversion
SEGMENT_SIZE = 501

# setting CSV data directory
DATA_DIR = "data/"

# setting file
import argparse
parser = argparse.ArgumentParser(description = "")
parser.add_argument('--target', required=True)
args = parser.parse_args()

t = args.target
targetPath = DATA_DIR + t + ".csv"
outputPath = DATA_DIR + t + "_fft_converted.csv"

# setting input file
targetFile = open(targetPath, "r")
print("DEBUG: TARGET FILE PATH:", targetPath)

outputFile = open(outputPath, "w")
print("DEBUG: OUTPUT FILE PATH:", outputPath)

currentSegmentTime = deque(maxlen = SEGMENT_SIZE)
currentSegmentValue = deque(maxlen = SEGMENT_SIZE)

dateparse = lambda dates: [pd.datetime.strptime(d, '%Y-%m-%d %H:%M:%S.%f') for d in dates]

df = pd.read_csv(targetFile, parse_dates=[0], date_parser=dateparse, header=None)

# make initial segment
for i in range(SEGMENT_SIZE):
    date_time = df[0][i]
    value = df[1][i]
    currentSegmentTime.append(int(datetime.datetime(date_time).strftime('%s')))
    currentSegmentValue.append(int(df[1][i]))

plt.plot(currentSegmentTime, currentSegmentValue)

# write csv headers
'''
csvWriter.writerow(["timestamp", "raw_value", "wavelet_value"])
csvWriter.writerow(["datetime", "int", "float"])
csvWriter.writerow(["T", "", ""])
'''

# generate headers
header_row_1 = ["timestamp"]
for i in range(SEGMENT_SIZE / 2):
    header_row_1.append("f" + str(i))
header_row_2 = ["datetime"] + (["int"] * (SEGMENT_SIZE / 2))
header_row_3 = ["T"] + ([""] * (SEGMENT_SIZE / 2))

csvWriter.writerow(header_row_1)
csvWriter.writerow(header_row_2)
csvWriter.writerow(header_row_3)

# https://en.wikipedia.org/wiki/Kaiser_window
# http://docs.scipy.org/doc/numpy-1.10.0/reference/generated/numpy.kaiser.html
kaiser_filter = numpy.kaiser(SEGMENT_SIZE, 14)
#plt.plot(kaiser_filter)
#plt.title("Kaiser window")
#plt.ylabel("Amplitude")
#plt.xlabel("Sample")
#plt.show()

for row in csvReader:
    date, value = row[0], int(row[1])
    currentSegment.append(value)
    window = currentSegment * kaiser_filter
    # FFTValue = numpy.log(numpy.abs(numpy.fft.fft(window).real))
    # FFTValue = numpy.abs(numpy.fft.fft(window).real)
    # FFTValue = pywt.dwt(currentSegment ,"db1")[0]
    # FFTValue = numpy.log([numpy.sqrt(c.real ** 2 + c.imag ** 2) for c in numpy.fft.fft(window)][0:(FFT_SEGMENT_SIZE/2)])
    FFTValue = map(lambda x:float(x),
                   [numpy.sqrt(c.real ** 2 + c.imag ** 2) for c in numpy.fft.fft(window)][0:(SEGMENT_SIZE / 2)])

    # csvWriter.writerow(FFTValue)
    csvWriter.writerow([date] + FFTValue)
    # csvWriter.writerow((date, value, FFTValue))

