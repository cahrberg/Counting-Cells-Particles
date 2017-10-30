# Christian D. Ahrberg - BNTL Sogang University - 2017
#
# Programm opening a series of images, finding and marking circles in them. Afterwards brightness of circles is measured
# and plotted in a histogram. A threshold fluorescence is entered manually and the concentration of dPCR experiments determined
# Program file should be in the same folder with the image files. The images, showing the wells should be labeled numerically usiing integers
# starting at 1. The Code creates the following files in the same folder:
#
# Xmarked_circles.jpg       -   Image X with the circles found circles marked in it
# fittedXmarked_circles.jpg -   Image X with the fitted array of circles marked in it
# circles.txt               -   Tab delimited text file with coordinates of found circles
#                               x-coordinate -- y-coordinate -- radius -- #image number -1   -   All values in pixels
# cirlces-fitted.txt        -   Tab delimited text file with coordinates of found circles, formating as before
# intensity_table.txt       -   Tab delimited text file with measured brightness of found circles
#                               brighness (from 0 to 1) -- #image number -1
# intensity_table-fitted.txt-   Tab delimited text file with measured brightness of fitted circle arrays
#                               Formating as before
# temp-plot.html            -   Histogram as created by Plotly
#
###########################################################################################

# Importing required packages

import numpy as np
import cv2
import plotly.graph_objs as go
import plotly as py
import math


# ==============================================================================================================
# Defining Variables can be changed by user

# Variables defining well array
NPan        = 3                         # Number of pannels that are analysed (i.e. the number of image files)
ArraySize   = 30                        # Number of wells in one column of a pannel (assumption of square arrays is made here)
VolWell     = 1.257 * math.pow(10,-5)   # Volume of a single well in uL
xdist       = 24.5                      # Horizontal distance between wells in array in pixcles
ydist       = 24.5                      # Vertical distance between wells in array in pixcles

# Variables for circle detection
dp1         = 1     # Ratio of accumulator, one means same resolution as input image, bigger numbers mean image is reduced 
minDist1    = 20    # Minimum distance between two circle centers in um
param11     = 1     # Threshold passed to Canny edge detector
param21     = 7    # Accumulatoir threshold for circles, the lower the more false circles are recognised
minRadius1  = 5     # Minimum radius of found circles in pixcels
maxRadius1  = 13    # Maximum radius of found circles in pixcels

# ==============================================================================================================

# Function for esimating Concentration

def ConcCallculation(pHat,Npart,Vol):
    # Function callculating the concentration of outcome of dPCR
    #
    # Parameters:
    # pHat  =   estimated propability of positive partition
    # Npart =   total number of partitions
    # Vol   =   Volume of partitions in uL
    #
    # Returns:
    # C_est = callculated concentration in #particles/uL
    # C_low = lower confidence intervall of calculated concentration (95% z-distribution) in # particles/uL
    # C_upp = upper confidence intervall of calculated concentration (95% z-distribution) in # particles/uL
    #
    #######################################################

    # Defingin constants
    zc = 1.96   # 95% confidence intervall z-distribution

    # Callculation of confidence interval on pHat
    pHat_Dev = zc * math.sqrt((pHat * (1-pHat))/Npart)  # Deviation on expected result
    p_hat_low = pHat - pHat_Dev  # Lower bound of p_hat
    p_hat_upp = pHat + pHat_Dev  # Upper bound of p_hat

    # Callculating mean number of molecules per patition including 95%
    # confidence intervall
    lambda1 = -math.log(1-pHat)     # average number of molecules per division as per Poission distribution
    lambda_low = -math.log(1-p_hat_low)  # lower bound of average number of molecules per division
    lambda_upp = -math.log(1-p_hat_upp)  # upper bound of average number of molecules per division

    # Callculating concentrations in mol/uL from lambda values including
    # confidence intervalls
    C_est = lambda1 / Vol       # Esitmated concentration
    C_low = lambda_low / Vol    # Estimated lower bound of concentration
    C_upp = lambda_upp / Vol    # Estimated higher bound of concentration

    return C_est, C_low, C_upp

# ==============================================================================================================

# Function for Callculating flourescence intensity inside of cicle

def CircleIntensity(centerx,centery,radius,image):
    # Callculates the intensity indside of a circle
    #
    # Parameters:
    # centerx = x-coordinte of circle
    # centery = y-coordinate of circle
    # radius = Radius of circle
    # image = image for callculating intensity (brightness saved as 8 bit, i.e. from 0 to 255
    #
    # Returns:
    # Intensity = average intensity value of circle in the tested image

    # Definging required parameters
    npixels = 0.0     # Count for pixels to find average
    brightness = 0.0  # Count for brightness of pixel

    # Creating squrare around circle
    for x in range(int(round(centerx-radius-2)),int(round(centerx+radius+2))):        # Varying through x
            for y in range((int(centery-radius-2)),int(round(centery+radius+2))):       # Varying through y
                 if x <= image.shape[0]-1 and y <= image.shape[1]-1 and x>=0 and y>=0:  # Making sure coordinate in image
                     pixeldistance = ((centerx - x)**2 + (centery - y)**2)**0.5         # Pythagoras to find radius from iterated pixcel to center of circle
                     if pixeldistance <= radius:                                        # If Pixel is in circle add to intensity callculation
                         brightness = brightness + float(image[x,y])/255                # Updating total brightness
                         npixels = npixels + 1                                          # Updating total pixcel count
    if npixels == 0:
        npixels = 1 # Preventing error, division by zero

    Intensity = brightness / npixels                                                    # Callculating average intesnity of circle

    return Intensity

# ==============================================================================================================

# Function for plotting histograms
    # Plots hiostogram of data using the Plotly ploting toolbox, opens new window in browser with interactive histogram
    #
    # Parameters:
    # Filename = Text file containing circle intensities

def Histogram(Filename):
    # Importing data for histogram without fitting
    ff = open(Filename,'r')
    lines = ff.readlines()
    f = []
    # Dividing imported data into format for plotting
    for x in lines:
        f.append(x.split('\t')[0])
    ff.close()
        
    plotdata = map(float, f)


    # Making histogram
    data = [
        go.Histogram(
            x=plotdata, xbins=dict(start=0, size=0.005, end=1) 
        )
    ]

    layout = go.Layout(
        title='Histogram of flourescence intensity of wells',
        xaxis=dict(title='Flourescence intensity',
                   range=[0,1]),
        yaxis=dict(title='Well count'), 
        )

    fig = go.Figure(data=data, layout=layout)

    py.offline.plot(fig)

# ===============================================================================================================

# Main code file calling the previous defined functions

# Creating file for saving circles
file=open("circles.txt","w")
file.close()

# Creating file for saving intensities
intensity_table = open("intensity_table.txt","w") #Creating file for results
intensity_table.close()

# Images filetype
name2 = '.jpg'
path = '../circles/'
tag = 'marked_circles'

# Looop going through the number of images defined by NPan, finding wells, fitting array of wells, and measuring well intensities
for x in range (0,NPan):

    # Debugging
    print "Currently on image number %d" % (x+1)

    # Creating filename for saving images with marked circles
    y = x+1
    name1 = str(y)
    name = name1 + name2

    # Converting image for different functions
    img = cv2.imread(name,0) # Opening image
    img = cv2.medianBlur(img,5) # Smothening data from image
    # Converting image to greyscale, plus copies for different applications
    cimg = cv2.cvtColor(img,cv2.COLOR_GRAY2BGR)
    cimg2 = cv2.cvtColor(img,cv2.COLOR_GRAY2BGR)
    gimg = cv2.cvtColor(cimg,cv2.COLOR_BGR2GRAY)
    
    # Fitting circles using Hough transform from package cv2
    # function calls as follows: cv2.HoughCircles(image, method, dp, minDist, circles, param1, param2, minRadius, maxRadius)
    circles = cv2.HoughCircles(img,cv2.HOUGH_GRADIENT,dp1,minDist1,param1=param11,param2=param21,minRadius=minRadius1,maxRadius=maxRadius1)

    # List for saveing x and y values of grid
    xlist = [] 
    ylist = []
    circlesfitted = [] 

    # Saving circels to tab delimited text file
    file=open("circles.txt","a")
    for j in circles[0,:]:
        file.write(str(j[0]) +'\t' + str(j[1]) + '\t' + str(j[2]) + '\t' + str(x))
        file.write("\n")
        xlist.append(j[0])
        ylist.append(j[1])
    file.close()

    # Finding bottom left corner of array, done by finding the y-coordinate of the lowest circle and the x-coordinate of the most left circle
    xmin = min(xlist)
    ymin = min(ylist)

    # Fitting array
    # Directly writing to file
    file=open("circles-fitted.txt","a")
    for xcount in range(ArraySize): # Going through all columns
        for ycount in range(ArraySize): # Going through all rows
            file.write(str(xmin+24.5*(xcount)) +'\t' + str(ymin+24.5*(ycount)) + '\t' + str(6) + '\t' + str(x))
            file.write("\n")
            circlesfitted.append([int(xmin+xdist*(xcount)),int(ymin+ydist*(ycount)),6]) # For plotting fitted circels to image
    file.close()

    # Drawing fitted circles
    circles = np.uint16(np.around(circles)) # Preparing data for plotting
    for i in circles[0,:]:
        # draw the outer circle
        cv2.circle(cimg,(i[0],i[1]),i[2],(0,255,0),2)
        # draw the center of the circle
        cv2.circle(cimg,(i[0],i[1]),2,(0,0,255),3)


    # Creating path for saving
    target = name1 + tag + name2
    # Saving image with marked circles
    cv2.imwrite(target,cimg)

    # Drawing fitted circles, inverted colours to image with fitted circles for differentiation
    circles = np.uint16(np.around(circles)) # Preparing data for plotting
    for k in circlesfitted[:]:
        # draw the outer circle
        cv2.circle(cimg2,(k[0],k[1]),k[2],(0,0,255),2)
        # draw the center of the circle
        cv2.circle(cimg2,(k[0],k[1]),2,(0,255,0),3)

    # Creating path for saving
    target = 'fitted' + name1 + tag + name2
    # Saving image with marked circles
    cv2.imwrite(target,cimg2)                             

    # Measuring intesity found cirlces
    intensity_table = open("intensity_table.txt","a")
    for i in circles[0,:]:
        intensity = CircleIntensity(i[0],i[1],i[2],gimg)
        # Wirting to file
        intensity_table.write(str(intensity) + '\t' + str(x))
        intensity_table.write("\n")
    intensity_table.close()

    # measuring intensity fitted circles
    intensity_table = open("intensity_table-fitted.txt","a")
    for i in circlesfitted[:]: # Measuring intensity of all circles using previous function
        intensity = CircleIntensity(i[0],i[1],i[2],gimg)
        # Wirting to file
        intensity_table.write(str(intensity) + '\t' + str(x))
        intensity_table.write("\n")
    intensity_table.close()

# Plotting results and calculating number of particles

# Histogram without fitting
Histogram("intensity_table.txt")

# Histogram with fitting
Histogram("intensity_table-fitted.txt")

# Callculating Concentration

Threshold = input('Please enter fluorescence threshold for positive call: ')
BeadsOrBacteria = input('Are you counting beads (1) or bacteria (2), please enter according number: ')

# Defining variables
Npos = 0.0 # Counter for positive wells, float so probability can be calculated
Nneg = 0.0 # Counter for negative wells, float so probability can be calculated

# Importing Analysing positive and negative wells:
ff = open("intensity_table-fitted.txt",'r')
lines = ff.readlines()
f = []
for x in lines:
    f.append(x.split('\t')[0])
ff.close()
        
plotdata = map(float, f)

# Counting positive and negative partitions
for x in plotdata:
    if x < Threshold:       # If well is positive average brightness is low
        Npos = Npos + 1
    if x >= Threshold:      # Wells contianing now particle/cell have a higher intensity
        Nneg = Nneg + 1

if BeadsOrBacteria == 2: # Exchanging Npos and Nneg for Flourescent labeled bacteria
    Nsave = Npos
    Npos  = Nneg
    Nneg  = Nsave


# Callculation of concentrations according to Poission distribution
pHat = Npos / (Npos + Nneg) # Estimated propability of positive well
Npart = Npos + Nneg # Total number of detected wells
C_est, C_low, C_upp = ConcCallculation(pHat,Npart,VolWell) # Callculation of concentrations according to predifined function

print('Number of Poitive calls:                 {0:.2f}' .format(Npos))
print('Number of Negative calls:                {0:.2f}' .format(Nneg))
print('Propability of positive partition:       {0:.2f} \n' .format(pHat))
# Outputting results
print('Estimated concentration:                 {0:.3f} copies/uL' .format(C_est))
print('Lower bound of 95% Confidence intervall: {0:.3f} copies/uL' .format(C_low))
print('Upper bound of 95% Confidence intervall: {0:.3f} copies/uL' .format(C_upp))










